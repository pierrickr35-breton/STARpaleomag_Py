"""
Lecteur Python équivalent à la subroutine Fortran importtexte()
pour remplir la structure Pmag définie dans starmac_OSX.inc

Format supporté : fichiers .ren.txt "nouveau format" (newformat=.TRUE. en Fortran)
"""

import re
from dataclasses import dataclass, field
from typing import List, Optional
import csv


# ---------------------------------------------------------------------------
# Structures équivalentes au Fortran STRUCTURE /pmagdata/
# ---------------------------------------------------------------------------

@dataclass
class Measurement:
    """Equivalent d'une ligne de mesure (etape, cod1, cod2, x, y, z, q, ins, s).

    `etape` est maintenant un FLOAT representant directement la valeur
    physique du pas (mT pour A/F, degC pour D/S/T/K...) - demande
    explicite utilisateur ("on importing legacy files from Rennes, we
    can divide by 10 from oersted to mT and in the software convert all
    step integer to float"), qui REVIENT sur une decision anterieure
    ("modify the whole starmac" avait ete rejete au profit du seul champ
    `step_value` cible - voir git history). Avant ce changement, `etape`
    restait un ENTIER portant, pour un pas AF (A/F), l'equivalent
    Oersted historique (mT*10, 1 Oe = 0.1 mT) plutot que le mT reel -
    .prmag lui-meme a toujours stocke le vrai mT dans sa colonne "step"
    (voir convert_ren_to_r.py/convert_magic_to_r.py), seule la
    representation EN MEMOIRE gardait cette echelle heritee. `etape`
    remplace desormais aussi `step_value` (supprime, redondant)."""
    etape: float
    cod1: str
    cod2: str
    x: float
    y: float
    z: float
    q: int
    ins: str
    s: float
    xech: Optional[float] = None
    yech: Optional[float] = None
    zech: Optional[float] = None
    heuremes: Optional[str] = None
    # champ labo (uT) - absent des .ren (deduit du commentaire par
    # paleointensity.parse_com_field), rempli directement pour les
    # fichiers .prmag depuis treat_dc_lowfield (voir read_prmag_file) -
    # demande explicite utilisateur ("extract the field value from
    # treat_dc_field").
    treat_dc_field: Optional[float] = None
    # 'g'/'b' (good/bad) - colonne "quality" du .prmag, jusque-la jamais
    # lue (seule "error" -> `q` l'etait) - demande explicite utilisateur
    # ("je voudrais utiliser le critere de qualite b/g dans les donnees
    # pour ne pas prendre en compte cette etape... un des problemes en
    # paleointensite est eventuellement qu'une serie d'echantillons ne
    # soient pas mis correctement dans le four... on peut mettre la
    # qualification b"). "g" par defaut (mesure normale) si la colonne
    # est absente/vide - voir datatools.remove_bad_quality_steps.
    quality: str = "g"


@dataclass
class Pmag:
    """Equivalent de RECORD /pmagdata/ pmag(5000)"""
    id: str = ""
    cin: float = 0.0
    caz: float = 0.0
    dip: float = 0.0
    str_: float = 0.0          # 'str' est un mot réservé Python -> str_
    norme: str = ""            # 'v' ou 'm'
    vol: float = 0.0
    com: str = ""
    lat: float = 0.0
    rlong: float = 0.0
    altitude: float = 0.0
    year: int = 0
    month: int = 0
    day: int = 0
    hour: int = 0
    minute: int = 0
    azmag: float = 0.0
    azsun: float = 0.0
    outilorient: str = ""
    roche: str = ""
    # Champs MagIC decodes depuis `roche` (voir _decode_roche) : la ligne
    # "roche" est ecrite par extract_magic.py au format
    # `Site: "..." Sample: "..." Fm: "..." Age: "..." GC: "..." SMT: "..."
    # Li: "..." Loc: "..." Obs: "..."` lors d'un import MagIC - GC =
    # Geologic Classes, SMT = Geologic Type (ER_Site.class / ER_Site.type1
    # dans fichiers_magic.f), Fm = Formation, Li = Lithology, Loc = Location.
    # Vides si la ligne roche ne suit pas ce format (anciens fichiers .ren).
    magic_site: str = ""
    magic_sample: str = ""
    magic_fm: str = ""
    magic_age: str = ""
    magic_gc: str = ""
    magic_smt: str = ""
    magic_li: str = ""
    magic_loc: str = ""
    magic_obs: str = ""
    # Position stratigraphique du specimen dans une section/carotte -
    # champ MagIC reel (`samples.height`, "Stratigraphic Height", metres,
    # positif vers le haut - verifie dans pmagpy/data_model/data_model.json)
    # plutot qu'un nom invente : demande explicite utilisateur ("add an
    # additional variable: stratigraphic_position (or the Magic
    # equivalent)... this field will replace the need to load a file for
    # magnetostratigraphic studies" - voir List and depth.../
    # ouvrir_lismesdepth_dialog, qui exigeait jusqu'ici un fichier externe
    # specimen/depth separe). None si non renseigne (specimen sans
    # position connue, ou format d'origine - .ren/Utrecht - qui n'a pas ce
    # concept).
    stratigraphic_height: Optional[float] = None
    mesures: List[Measurement] = field(default_factory=list)

    @property
    def nbmes(self) -> int:
        return len(self.mesures)


_ROCHE_KEY_FIELDS = {
    "site": "magic_site",
    "sample": "magic_sample",
    "fm": "magic_fm",
    "age": "magic_age",
    "gc": "magic_gc",
    "smt": "magic_smt",
    "li": "magic_li",
    "loc": "magic_loc",
    "obs": "magic_obs",
}
_ROCHE_PAIR_RE = re.compile(r'(\w+):\s*"([^"]*)"')
# Ancien format (fichiers campagne Corbieres, ~2016) : valeurs SANS guillemets,
# separees par '|' et collees en fin de ligne 'L:' plutot que sur une ligne a
# part - ex. "Fm:Sediments|Age:Paleocene|Li:carbonated_silts|Loc:Lairiere_1_16"
# (cf. commentaire modele dans orient_paleomag.f : roche="Fm:?|Age:?|Li:?|Loc:?").
_ROCHE_PAIR_RE_PIPE = re.compile(r'(\w+):([^|]*)')


def decode_roche(p: Pmag) -> None:
    """Decode `p.roche` dans les champs magic_* de `p` - mute `p` en place.
    Deux conventions acceptees (essayees dans l'ordre) : `Cle: "valeur"`
    (guillemets, ecrite par extract_magic.py, voir _ROCHE_KEY_FIELDS) et
    l'ancienne `Cle:valeur|Cle:valeur` sans guillemets (fichiers Corbieres).
    Cle insensible a la casse (constate dans des fichiers .ren reels : "obs"
    en minuscule alors que extract_magic.py ecrit toujours "Obs" en
    majuscule). N'importe quel sous-ensemble/ordre des cles est accepte ;
    les cles absentes ou non reconnues sont ignorees, et sur un ancien
    fichier sans aucune de ces conventions (0 correspondance), tous les
    champs magic_* restent a leur valeur par defaut (chaine vide)."""
    text = p.roche or ""
    pairs = _ROCHE_PAIR_RE.findall(text)
    if not pairs:
        pairs = _ROCHE_PAIR_RE_PIPE.findall(text)
    for key, value in pairs:
        attr = _ROCHE_KEY_FIELDS.get(key.strip().lower())
        if attr is not None:
            setattr(p, attr, value.strip())


# ---------------------------------------------------------------------------
# Fonctions utilitaires de parsing
# ---------------------------------------------------------------------------

def _safe_float(txt: str, default: float = 0.0) -> float:
    try:
        return float(txt)
    except (ValueError, TypeError):
        return default


def _safe_int(txt: str, default: int = 0) -> int:
    try:
        return int(float(txt))
    except (ValueError, TypeError):
        return default


def parse_id_line(line: str) -> Pmag:
    """
    Parse une ligne du type :
    'Id:19DN0101A    in: 49.0 az:195.0 dip:  0.0 str: 12.0 v: 10.8 com: ...'
    (équivalent de la boucle 'do llk=1,len(chaine)' en Fortran)
    """
    p = Pmag()

    m = re.search(r'Id:\s*(\S+)', line)
    if m:
        p.id = m.group(1)[:12]

    m = re.search(r'\bin:\s*([-\d.]+)', line)
    if m:
        p.cin = _safe_float(m.group(1))

    m = re.search(r'\baz:\s*([-\d.]+)', line)
    if m:
        p.caz = _safe_float(m.group(1))

    m = re.search(r'\bdip:\s*([-\d.]+)', line)
    if m:
        p.dip = _safe_float(m.group(1))

    m = re.search(r'\bstr:\s*([-\d.]+)', line)
    if m:
        p.str_ = _safe_float(m.group(1))

    # 'v:' (volume, norme='v') ou 'm:' (masse, norme='m')
    m = re.search(r'(?<!a)\bv:\s*([-\d.]+)', line)
    if m:
        p.norme = 'v'
        p.vol = _safe_float(m.group(1))
    else:
        m = re.search(r'\bm:\s*([-\d.]+)', line)
        if m:
            p.norme = 'm'
            p.vol = _safe_float(m.group(1))

    m = re.search(r'com:\s*(.*)$', line)
    if m:
        p.com = m.group(1).strip()

    return p


def parse_L_line(p: Pmag, line: str) -> None:
    """
    Parse la ligne d'orientation/date, ex :
    'L: -55.01404 G: -67.75302  H:   0  T:2019  1 13  0  0   azm: 93.0 azs:  0.0  Or:A12_0_3_90'
    """
    m = re.search(r'L:\s*([-\d.]+)', line)
    if m:
        p.lat = _safe_float(m.group(1))

    m = re.search(r'G:\s*([-\d.]+)', line)
    if m:
        p.rlong = _safe_float(m.group(1))

    m = re.search(r'H:\s*([-\d.]+)', line)
    if m:
        p.altitude = _safe_float(m.group(1))

    # T: year month day hour minute
    m = re.search(r'T:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
    if m:
        p.year, p.month, p.day, p.hour, p.minute = (int(x) for x in m.groups())

    m = re.search(r'azm:\s*([-\d.]+)', line)
    if m:
        p.azmag = _safe_float(m.group(1))

    m = re.search(r'azs:\s*([-\d.]+)', line)
    if m:
        p.azsun = _safe_float(m.group(1))

    m = re.search(r'Or:\s*(\S+)', line)
    if m:
        p.outilorient = m.group(1)
        # ancien format (fichiers Corbieres) : metadonnees roche collees a
        # la suite de 'Or:...' sur cette meme ligne, plutot que sur une
        # ligne separee - voir decode_roche pour la convention sans
        # guillemets utilisee ici ("Fm:...|Age:...|Li:...|Loc:...").
        trailing = line[m.end():].strip()
        if trailing:
            p.roche = trailing


# cod1 pour lesquels un fichier .ren HERITE encode le pas AF en Oersted
# (entier brut, ex. "1300F+" = 1300 Oe) - convertis en mT reel (/10, 1 Oe
# = 0.1 mT) DES LA LECTURE par parse_measure_line, pour que Measurement.
# etape soit uniformement la valeur physique reelle partout dans
# l'application, plus seulement dans le fichier .prmag sur disque -
# demande explicite utilisateur ("on importing legacy files from Rennes,
# we can divide by 10 from oersted to mT... convert all step integer to
# float").
_LEGACY_AF_CODES = {"A", "F"}


def parse_measure_line(line: str) -> Optional[Measurement]:
    """
    Parse une ligne de mesure, ex :
    '    0N0  6.918E-07  2.046E-08 -1.014E-06   0 C1  213.0'
    ou (format long avec xech,yech,zech,heuremes) :
    '  1300F+  9.406E-08 -1.638E-08 -8.055E-08   0 C1    0.0  1.2E-05 3.4E-05 5.6E-05  2019/01/13 12:00:00'
    Equivalent aux formats Fortran 201 et 2011. Le pas AF (A/F) est
    stocke en Oersted dans ce format herite - converti en mT reel (/10)
    ici, voir _LEGACY_AF_CODES.
    """
    parts = line.split()
    if len(parts) < 7:
        return None

    step_field = parts[0]  # ex '0N0', '150D+', '1300F+'
    m = re.match(r'(\d+)([A-Za-z])(\S)', step_field)
    if not m:
        return None

    cod1 = m.group(2)
    cod2 = m.group(3)
    etape = int(m.group(1)) / 10.0 if cod1 in _LEGACY_AF_CODES else float(m.group(1))

    try:
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        q = _safe_int(parts[4])
        ins = parts[5]
        s = float(parts[6])
    except ValueError:
        return None

    meas = Measurement(etape, cod1, cod2, x, y, z, q, ins, s)

    # champs optionnels (format long, équivalent format 2011)
    if len(parts) >= 10:
        meas.xech = _safe_float(parts[7])
        meas.yech = _safe_float(parts[8])
        meas.zech = _safe_float(parts[9])
        if len(parts) >= 12:
            meas.heuremes = parts[10] + " " + parts[11]
        elif len(parts) == 11:
            meas.heuremes = parts[10]

    return meas


# ---------------------------------------------------------------------------
# Lecture du fichier complet
# ---------------------------------------------------------------------------

def _read_text_auto(filepath: str, encoding: Optional[str]) -> str:
    """Lit le contenu texte d'un fichier .ren en devinant l'encodage si
    non precise explicitement (`encoding=None`) : tente UTF-8 d'abord
    (convention des .ren modernes, ecrits par convert_legacy_ren.py/
    convert_ren_to_r.py, et de plus en plus frequente aussi pour un
    fichier edite aujourd'hui), puis latin-1 en repli (convention des
    vieux fichiers Rennes pre-Unicode) - demande explicite utilisateur
    ("when importing the ren files the ± is not well decoded") : un
    caractere UTF-8 valide (ex. ± = octets 0xC2 0xB1) relu a tort en
    latin-1 (l'ancien defaut INCONDITIONNEL de read_ren_file) se
    retrouvait mal-decode en deux caracteres "Â±", et ENCORE une fois si
    ce texte deja corrompu repassait par la meme lecture erronee plus
    loin dans le pipeline (chaine complete "Import Starmac legacy
    files..." : convert_legacy_ren.py ecrit TOUJOURS son .ren
    intermediaire en UTF-8, jamais relu comme tel jusqu'ici). Un octet
    latin-1 isole hors ASCII (ex. 0xE9 pour 'é') n'est quasiment jamais
    une sequence UTF-8 valide a lui seul, ce qui rend cette detection
    fiable en pratique - un `encoding` explicite (ancien comportement)
    reste bien sur prioritaire et jamais devine."""
    if encoding is not None:
        with open(filepath, "r", encoding=encoding) as f:
            return f.read()
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return f.read()
    except UnicodeDecodeError:
        with open(filepath, "r", encoding="latin-1") as f:
            return f.read()


def read_ren_file(filepath: str, encoding: Optional[str] = None) -> List[Pmag]:
    """
    Lit un fichier .ren.txt et retourne la liste des échantillons (Pmag),
    équivalent au tableau pmag(1:nb_ech) rempli par importtexte().

    Lecture des 3 lignes d'entete (Id:/L:/roche) guidee par le CONTENU de
    chaque ligne plutot que par sa seule position, comme le fait le Fortran
    d'origine (`if(chaine(1:1)=="L")` avant de consommer la ligne suivante
    comme entete coordonnees/roche - importtexte, fichiers_mod_magic.f) :
    la ligne 'L:' et/ou la ligne site/roche peuvent etre absentes pour un
    echantillon donne (vieux fichiers, edition manuelle) sans faire perdre
    la premiere mesure qui suivrait immediatement 'Id:' - demande explicite
    de l'utilisateur, sur le modele deja tolerant de `decode_roche` pour la
    ligne 3 (n'importe quel sous-ensemble/ordre des cles est accepte).

    `encoding` : None (par defaut) devine UTF-8/latin-1 - voir
    _read_text_auto ; passer une valeur explicite pour forcer un
    encodage precis (comportement inchange pour les appelants qui le
    faisaient deja)."""
    pmag_list: List[Pmag] = []
    current: Optional[Pmag] = None

    # attendu : "L" = ligne coordonnees, "roche" = ligne site/roche, "mesures"
    expect = "L"

    content = _read_text_auto(filepath, encoding)
    for raw_line in content.splitlines():
        line = raw_line.rstrip("\n")

        if line.strip() == "":
            continue

        if line.lstrip().startswith("Id:"):
            if current is not None:
                pmag_list.append(current)
            current = parse_id_line(line)
            expect = "L"
            continue

        if current is None:
            continue

        if expect in ("L", "roche") and parse_measure_line(line) is not None:
            # ligne 'L:' et/ou roche absente(s) pour cet echantillon :
            # ce qui suit 'Id:' est deja une mesure - ne pas la perdre.
            expect = "mesures"

        if expect == "L":
            parse_L_line(current, line)
            if current.roche:
                # roche deja decodee depuis la fin de cette meme ligne
                # (ancien format Corbieres) - pas de ligne separee a
                # attendre en plus.
                decode_roche(current)
                expect = "mesures"
            else:
                expect = "roche"
            continue

        if expect == "roche":
            current.roche = line.strip()
            decode_roche(current)
            expect = "mesures"
            continue

        meas = parse_measure_line(line)
        if meas is not None:
            current.mesures.append(meas)

    if current is not None:
        pmag_list.append(current)

    return pmag_list


def _prmag_nd(txt: Optional[str], default: float = 0.0) -> float:
    txt = (txt or "").strip()
    if not txt or txt == "n.d":
        return default
    return _safe_float(txt, default)


def _prmag_nd_opt(txt: Optional[str]) -> Optional[float]:
    """Comme _prmag_nd mais retourne None (pas 0.0) si absent - pour
    distinguer "champ non applicable a ce pas" de "champ reellement nul"."""
    txt = (txt or "").strip()
    if not txt or txt == "n.d":
        return None
    try:
        return float(txt)
    except ValueError:
        return None


def _prmag_kv_line(line: str) -> dict:
    """Parse une ligne 'cle: valeur\\tcle2: valeur2...' (bloc d'entete
    .prmag) en dict {cle: valeur}."""
    result = {}
    for chunk in line.split("\t"):
        if ":" not in chunk:
            continue
        k, _sep, v = chunk.partition(":")
        result[k.strip()] = v.strip()
    return result


def _prmag_text(v: Optional[str]) -> str:
    v = (v or "").strip()
    return "" if v == "n.d" else v


def read_prmag_file(filepath: str, encoding: str = "utf-8") -> List[Pmag]:
    """Lit un fichier .prmag (nouveau format, voir convert_ren_to_r.py /
    convert_magic_to_r.py) et retourne la meme structure `List[Pmag]` que
    `read_ren_file`, pour que TOUT le reste de l'application (Zijderveld,
    stereo, statistiques...) fonctionne sans changement, qu'un fichier
    .ren ou .prmag ait ete ouvert - demande explicite utilisateur ("the
    open file still look for .ren").

    Bloc par specimen : 4 lignes d'entete (specimen/sample/site/volume/
    mass/lat/lon/elevation/stratigraphic_height/comment ; azimuth/dip/
    date/magnetic_azimuth/solar_azimuth/orient_tool ; bed_dip_strike/
    bed_dip ; formation/age/geologic_classes/geologic_types/lithologies/
    location/obs/method_codes)
    puis une ligne d'entete de mesures (step/cod1/cod2/x/y/z/...) suivie
    des lignes de mesure - voir docstring de convert_ren_to_r.py pour le
    detail exact du format et des unites."""
    with open(filepath, "r", encoding=encoding) as f:
        lines = [raw.rstrip("\n") for raw in f]

    pmag_list: List[Pmag] = []
    i, n = 0, len(lines)

    def _skip_blank_and_comments():
        nonlocal i
        while i < n and (lines[i].strip() == "" or lines[i].lstrip().startswith("#")):
            i += 1

    _skip_blank_and_comments()
    while i < n:
        if i + 4 >= n:
            break  # bloc incomplet en fin de fichier
        line_a = _prmag_kv_line(lines[i]); i += 1
        line_b = _prmag_kv_line(lines[i]); i += 1
        line_c = _prmag_kv_line(lines[i]); i += 1
        line_d = _prmag_kv_line(lines[i]); i += 1
        header_cols = [c.strip() for c in lines[i].split("\t")]; i += 1
        col_idx = {name: idx for idx, name in enumerate(header_cols)}

        p = Pmag()
        p.id = _prmag_text(line_a.get("specimen"))
        p.magic_sample = _prmag_text(line_a.get("sample"))
        p.magic_site = _prmag_text(line_a.get("site"))
        p.com = _prmag_text(line_a.get("comment"))
        p.lat = _prmag_nd(line_a.get("lat"))
        p.rlong = _prmag_nd(line_a.get("lon"))
        p.altitude = _prmag_nd(line_a.get("elevation"))
        p.stratigraphic_height = _prmag_nd_opt(line_a.get("stratigraphic_height"))
        vol_raw, mass_raw = line_a.get("volume", "n.d"), line_a.get("mass", "n.d")
        if _prmag_text(vol_raw):
            p.norme, p.vol = "v", _prmag_nd(vol_raw)
        elif _prmag_text(mass_raw):
            p.norme, p.vol = "m", _prmag_nd(mass_raw)

        p.caz = _prmag_nd(line_b.get("azimuth"))
        p.cin = _prmag_nd(line_b.get("dip"))
        m = re.match(r"(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2})", line_b.get("date") or "")
        if m:
            p.year, p.month, p.day, p.hour, p.minute = (int(g) for g in m.groups())
        p.azmag = _prmag_nd(line_b.get("magnetic_azimuth"))
        p.azsun = _prmag_nd(line_b.get("solar_azimuth"))
        p.outilorient = _prmag_text(line_b.get("orient_tool"))

        p.str_ = _prmag_nd(line_c.get("bed_dip_strike"))
        p.dip = _prmag_nd(line_c.get("bed_dip"))

        p.magic_fm = _prmag_text(line_d.get("formation"))
        p.magic_age = _prmag_text(line_d.get("age"))
        p.magic_gc = _prmag_text(line_d.get("geologic_classes"))
        p.magic_smt = _prmag_text(line_d.get("geologic_types"))
        p.magic_li = _prmag_text(line_d.get("lithologies"))
        p.magic_loc = _prmag_text(line_d.get("location"))
        p.magic_obs = _prmag_text(line_d.get("obs"))

        while i < n and lines[i].strip() != "":
            parts = lines[i].split("\t")
            i += 1

            def col(name: str) -> str:
                idx = col_idx.get(name)
                return parts[idx].strip() if idx is not None and idx < len(parts) else ""

            cod1 = col("cod1") or "?"
            cod2 = col("cod2") or "0"
            # `step` du .prmag est DEJA la valeur physique reelle (mT/degC -
            # voir convert_ren_to_r.py/convert_magic_to_r.py, qui l'ecrivent
            # ainsi pour TOUS les cod1, AF inclus) : plus de mise a l'echelle
            # Oersted (*10) a la lecture, `etape` = cette valeur directement.
            etape = _prmag_nd(col("step"))

            meas = Measurement(
                etape=etape, cod1=cod1, cod2=cod2,
                x=_prmag_nd(col("x")), y=_prmag_nd(col("y")), z=_prmag_nd(col("z")),
                q=int(round(_prmag_nd(col("error")))),
                ins=_prmag_text(col("instrument")),
                s=_prmag_nd(col("s")),
                treat_dc_field=_prmag_nd_opt(col("treat_dc_lowfield")),
                quality=_prmag_text(col("quality")) or "g",
            )
            p.mesures.append(meas)

        pmag_list.append(p)
        _skip_blank_and_comments()

    return pmag_list


# ---------------------------------------------------------------------------
# Export CSV (échantillons + mesures) - pratique pour vérification / usage externe
# ---------------------------------------------------------------------------

def export_specimens_csv(pmag_list: List[Pmag], filepath: str) -> None:
    fields = ["id", "cin", "caz", "dip", "str_", "norme", "vol", "com",
              "lat", "rlong", "altitude", "year", "month", "day", "hour",
              "minute", "azmag", "azsun", "outilorient", "roche", "nbmes"]
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        for p in pmag_list:
            writer.writerow([getattr(p, f) if f != "nbmes" else p.nbmes for f in fields])


def export_measurements_csv(pmag_list: List[Pmag], filepath: str) -> None:
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["specimen_id", "etape", "cod1", "cod2", "x", "y", "z",
                          "q", "ins", "s", "xech", "yech", "zech", "heuremes"])
        for p in pmag_list:
            for m in p.mesures:
                writer.writerow([p.id, m.etape, m.cod1, m.cod2, m.x, m.y, m.z,
                                  m.q, m.ins, m.s, m.xech, m.yech, m.zech, m.heuremes])


# ---------------------------------------------------------------------------
# Programme principal
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    fichier = "testfile.ren.txt"

    donnees = read_ren_file(fichier)

    print(f"{len(donnees)} sample(s) read\n")

    for p in donnees:
        print(f"ID: {p.id!r}")
        print(f"  cin={p.cin}  caz={p.caz}  dip={p.dip}  str={p.str_}  "
              f"norme={p.norme!r}  vol={p.vol}")
        print(f"  lat={p.lat}  rlong={p.rlong}  altitude={p.altitude}")
        print(f"  date={p.year}/{p.month}/{p.day} {p.hour}:{p.minute}")
        print(f"  azmag={p.azmag}  azsun={p.azsun}  outil={p.outilorient!r}")
        print(f"  roche: {p.roche}")
        print(f"  nb mesures: {p.nbmes}")
        for m in p.mesures:
            print(f"    etape={m.etape:>5} cod1={m.cod1} cod2={m.cod2} "
                  f"x={m.x: .3e} y={m.y: .3e} z={m.z: .3e} q={m.q} "
                  f"ins={m.ins} s={m.s}")
        print()

    # Export optionnel en CSV
    export_specimens_csv(donnees, "specimens.csv")
    export_measurements_csv(donnees, "mesures.csv")
    print("Export CSV : specimens.csv, mesures.csv")
    