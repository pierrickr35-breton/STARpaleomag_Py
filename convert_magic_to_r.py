"""
Convertit un fichier de "contribution MagIC" (telechargement combine
locations+sites+samples+specimens+measurements en un seul fichier texte,
tel que magic_contribution_XXXXX.txt) vers le nouveau format .r - test
d'interoperabilite explicitement demande par l'utilisateur ("prepare a
new import to check the interoperability"), sur des donnees REELLES
d'autres laboratoires (contributeur @ltauxe, PAS les propres donnees
Roperch/Rennes) plutot que sur les .ren internes deja converti par
convert_ren_to_r.py.

Contrairement a convert_ren_to_r.py (qui PART d'un cod1/cod2 Rennes deja
existant et le convertit vers les champs MagIC), ce module fait le
chemin INVERSE : il n'existe aucun cod1/cod2 dans la source, il faut les
INVENTER a partir de method_codes/treat_*, pour que le fichier .r
resultant garde la meme lisibilite (colonnes step/cod1/cod2) que les
fichiers convertis depuis un .ren.

Convention de lettrage des pas IZZI sans champ applique (LT-T-Z, "no
field") - demande explicite utilisateur ("convert also the paleointensity
codes to SA, SB etc when there is no field for the IZZI") : cod1='S'
(reprend le code Rennes deja existant pour un pas thermique zero-field
IZZI - voir magic_export.py cod1=='S'), cod2=lettre sequentielle (A,B,C..)
par temperature de premiere apparition parmi les pas LT-T-Z/LT-T-I/
LT-PTRM-I d'un meme specimen (heuristique : PAS garanti identique a la
convention Rennes originale pour les pas P, ceux-ci reutilisent la lettre
de LEUR PROPRE temperature si deja vue, sinon une nouvelle - a verifier
sur le fichier produit).

Champs manquants dans les donnees source (frequent : ces contributions
n'ont souvent ni azimuth/dip/bed_dip_strike/bed_dip au niveau sample -
echantillons archeologiques sans orientation de terrain, ex. tessons de
poterie) : "n.d", memes conventions que convert_ren_to_r.py."""

import argparse
import math
import os
import re
from typing import Dict, List, Optional, Tuple

from extract_magic import (
    _OUT_OF_SCOPE_PROTOCOLS, magic_results_to_redo_lines, magic_site_means,
    magic_pint_results_to_redo_lines, magic_anisotropy_rows,
)
from testlect import read_prmag_file
from calcul import (
    FitResult, archivres, fit_from_redo_file, results_path_for, dp_dm_from_a95,
    ani_path_for, write_ani_tensors_from_magic_rows,
)

FORMAT_HEADER = (
    "#STARpaleomag_Py .prmag v1  angles=deg  fields in milliTesla (mT) for strong "
    "fields AF or IRM and in microTesla (uT) for low field paleointensity "
    "or ARM  temperatures in degC  date=ISO8601"
)
S_FIELD_HEADER = "#s = magnetic susceptibility in 1e-5 SI"

_MEAS_HEADER = (
    "step\tcod1\tcod2\tx\ty\tz\terror\tquality\tinstrument\ts\t"
    "treat_temp\ttreat_ac_field\ttreat_dc_strongfield\ttreat_dc_lowfield\t"
    "treat_dc_field_phi\ttreat_dc_field_theta\t"
    "method_codes\tinstrument_codes\ttreat_step_num"
)


def parse_magic_contribution(path: str) -> Dict[str, List[Dict[str, str]]]:
    """Coupe le fichier combine en tables sur les marqueurs
    "tab delimited\\t<nom>", lit chaque table (ligne d'en-tete + lignes de
    donnees) en liste de dict {colonne: valeur}.

    Chaque table (sauf la derniere du fichier) se termine par une ligne
    ">>>>>>>>>>" - marqueur de fin de table du format "contribution
    combinee" MagIC reel (verifie sur magic_contribution_19491.txt :
    present a l'identique apres locations/sites/samples/specimens/
    measurements/criteria, absent seulement apres la derniere table -
    "ages" - qui s'arrete a la fin du fichier). Bug reel corrige, signale
    par l'utilisateur ("à la fin de l'importation d'une contribution
    Magic, il y a une ligne supplémentaire") : ce marqueur n'etait pas
    filtre et se retrouvait traite comme une derniere ligne de DONNEES de
    la table measurements, produisant un "specimen" fantome (id vide, 1
    mesure bidon) dans le .prmag converti."""
    tables: Dict[str, List[Dict[str, str]]] = {}
    with open(path, "r", encoding="utf-8") as f:
        raw = f.read()
    blocks = re.split(r"^tab delimited\t(\S+)\s*$", raw, flags=re.MULTILINE)
    # blocks[0] est vide (avant le 1er marqueur) ; puis alternance nom,corps
    for i in range(1, len(blocks), 2):
        name = blocks[i]
        body = blocks[i + 1].strip("\n")
        lines = [l for l in body.split("\n") if l.strip() != ">>>>>>>>>>"]
        if not lines or not lines[0]:
            tables[name] = []
            continue
        header = lines[0].split("\t")
        rows = []
        for line in lines[1:]:
            if not line:
                continue
            parts = line.split("\t")
            row = {header[j]: (parts[j] if j < len(parts) else "") for j in range(len(header))}
            rows.append(row)
        tables[name] = rows
    return tables


def _index_by(rows: List[Dict[str, str]], key: str) -> Dict[str, Dict[str, str]]:
    return {r[key]: r for r in rows if r.get(key)}


def _f(row: Dict[str, str], key: str) -> Optional[float]:
    v = (row or {}).get(key, "")
    if v in ("", None):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _s(row: Dict[str, str], key: str) -> str:
    v = (row or {}).get(key, "")
    return v.strip() if v else "n.d"


def _temp_k(row: Dict[str, str]) -> Optional[float]:
    """Temperature (K) du palier - `treat_temp` (cible du TRAITEMENT) en
    priorite, repli sur `meas_temp` (temperature REELLE de la MESURE,
    "Room temperature is 293" selon le data model MagIC 3 - verifie en
    ligne, earthref.org/MagIC/data-models/3.0.json) pour les fichiers NON
    CONVENTIONNELS ou treat_temp est absent mais le contributeur a
    neanmoins enregistre la vraie temperature de palier dans meas_temp -
    demande explicite utilisateur ("use these meas_field_ac meas_temp for
    AF and thermal to import the step and code in unconventional
    files"). Sans effet sur un fichier bien forme (treat_temp deja
    present, meas_temp jamais consulte)."""
    v = _f(row, "treat_temp")
    return v if v is not None else _f(row, "meas_temp")


def _ac_field_t(row: Dict[str, str]) -> Optional[float]:
    """Champ AF (Tesla) du palier - `treat_ac_field` (cible du TRAITEMENT)
    en priorite, repli sur `meas_field_ac` (champ REELLEMENT applique
    pendant la MESURE) pour les fichiers non conventionnels - MEME
    logique que `_temp_k`. -1 y est une valeur SENTINELLE ("ambient
    field", pas une vraie mesure de palier - verifie sur le data model
    MagIC 3 en ligne, "No field equals 0 and ambient field equals -1") :
    jamais recuperee comme un palier reel."""
    v = _f(row, "treat_ac_field")
    if v is not None:
        return v
    v2 = _f(row, "meas_field_ac")
    if v2 is None or v2 < 0:
        return None
    return v2


def _nd(value, fmt: Optional[str] = None) -> str:
    if value is None:
        if fmt:
            width_digits = "".join(ch for ch in fmt.split(".")[0] if ch.isdigit())
            if width_digits:
                return "n.d".rjust(int(width_digits))
        return "n.d"
    return format(value, fmt) if fmt is not None else str(value)


def _dir2cart(dec: Optional[float], inc: Optional[float], mag: Optional[float]):
    """Inverse de `selection.polere` (dec,inc,mag -> x,y,z) - necessaire
    car ces contributions MagIC stockent souvent dir_dec/dir_inc/
    magn_moment plutot que magn_x/y/z bruts."""
    if dec is None or inc is None or mag is None:
        return None, None, None
    dr, ir = math.radians(dec), math.radians(inc)
    x = mag * math.cos(ir) * math.cos(dr)
    y = mag * math.cos(ir) * math.sin(dr)
    z = mag * math.sin(ir)
    return x, y, z


def _iso_from_magic_ts(ts: str) -> str:
    m = re.match(r"(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2})", ts or "")
    return f"{m.group(1)}T{m.group(2)}" if m else "n.d"


def scan_dip_sign(path_in: str) -> Dict[str, object]:
    """Verifie le fichier AVANT conversion (demande explicite utilisateur
    "can you check the file before the import?") : recense, pour chaque
    specimen EFFECTIVEMENT converti (present dans measurements.txt -
    meme population que `by_specimen` dans convert_magic_file), le dip
    BRUT (non encore negate) qui sera utilise - specimen-level en
    priorite, sample-level en repli (voir _orientation_source).

    Retourne {"n_total", "n_positive", "n_negative", "n_zero",
    "n_from_specimen", "majority_positive"}. "majority_positive" est le
    signal d'alerte - REGLE METIER donnee explicitement par l'utilisateur
    ("en pratique on ne devrait avoir que des valeurs negatives <= a 0
    pour le dip donc quand la majorite des echantillons ont un dip >= a
    0 c'est qu'il y a un probleme") : PAS une condition d'unanimite
    (une premiere version exigeait "tous positifs" ou "tous negatifs",
    qui manquait exactement le cas reel signale par l'utilisateur -
    magic_contribution_20476.txt, 848 positifs sur 886, donc PAS "tous"
    mais bien la VASTE MAJORITE - "still not right for the sign when the
    magic has all the wrong positive dip"). Un dip brut MagIC positif
    est TOUJOURS anormal en pratique (jamais juste "un signe parmi
    d'autres") ; le seuil est donc simplement n_positive > n_total/2,
    jamais une histoire d'uniformite. Un dip negatif ou nul est SAIN,
    jamais signale (contrairement a une premiere version qui traitait
    aussi "tous negatifs" comme suspect - a tort, c'est le cas normal
    attendu)."""
    tables = parse_magic_contribution(path_in)
    samples = _index_by(tables.get("samples", []), "sample")
    specimens = _index_by(tables.get("specimens", []), "specimen")
    converted_ids = {row.get("specimen", "") for row in tables.get("measurements", [])}
    converted_ids.discard("")

    n_total = n_positive = n_negative = n_zero = n_from_specimen = 0
    for specimen_id in converted_ids:
        spec_row = specimens.get(specimen_id)
        sample_row = samples.get(spec_row.get("sample", "")) if spec_row else None
        _az_raw, dip_raw, from_specimen = _orientation_source(spec_row, sample_row)
        if dip_raw is None:
            continue
        n_total += 1
        if from_specimen:
            n_from_specimen += 1
        if dip_raw > 0:
            n_positive += 1
        elif dip_raw < 0:
            n_negative += 1
        else:
            n_zero += 1

    return {
        "n_total": n_total,
        "n_positive": n_positive,
        "n_negative": n_negative,
        "n_zero": n_zero,
        "n_from_specimen": n_from_specimen,
        "majority_positive": n_total >= 2 and n_positive > n_total / 2.0,
    }


def _orientation_source(spec_row, sample_row) -> Tuple[Optional[float], Optional[float], bool]:
    """(azimuth_raw, dip_raw, from_specimen) - MagIC autorise azimuth/dip
    au niveau SPECIMEN (specimens.txt), pas seulement au niveau sample
    (samples.txt) - demande explicite utilisateur ("in some Magic
    contributions, the orientation of the core azimuth and dip are given
    at the specimen level, when this is the case, please take it into
    account"). Priorite au niveau specimen quand au moins l'un des deux
    champs y est renseigne (les deux valeurs viennent alors ENSEMBLE du
    meme niveau, jamais un azimuth specimen combine a un dip sample ou
    l'inverse) ; repli sur le niveau sample sinon - AUCUN des deux niveau
    n'est jamais fabrique."""
    spec_az = _f(spec_row, "azimuth") if spec_row else None
    spec_dip = _f(spec_row, "dip") if spec_row else None
    if spec_az is not None or spec_dip is not None:
        return spec_az, spec_dip, True
    return (_f(sample_row, "azimuth"), _f(sample_row, "dip"), False)


_DIP_MODES = {
    "negate": lambda d: -d,   # standard MagIC convention (defaut)
    "raw": lambda d: d,       # dip .prmag = dip MagIC tel quel, sans transformation
    "utrecht": lambda d: 90.0 - d,  # convention Utrecht/PMAG2 (voir convert_utrecht_to_r.py)
}


def _sample_header_block(
    specimen: str, spec_row, sample_row, site_row, loc_row, dip_mode: str = "negate",
) -> str:
    sample_name = _s(spec_row, "sample") if spec_row else (sample_row.get("sample", "n.d") if sample_row else "n.d")
    site_name = _s(sample_row, "site") if sample_row else (site_row.get("site", "n.d") if site_row else "n.d")
    lat = _f(sample_row, "lat") if sample_row and sample_row.get("lat") else _f(site_row, "lat")
    lon = _f(sample_row, "lon") if sample_row and sample_row.get("lon") else _f(site_row, "lon")
    elevation = _f(site_row, "elevation")
    # convention STARpaleomag_Py != convention MagIC (voir DIFF_WITH_MAGIC) :
    # azimuth = azimuth de X MagIC +90 ; dip = -(dip MagIC) par defaut - meme
    # transformation que extract_magic.py (az_trans = core_az_raw+90,
    # inc_trans = -core_dip_raw), appliquee ici a l'import plutot qu'a
    # l'export - demande explicite utilisateur ("correct azimuth and dip
    # to be consistent with the comment in the header"). `dip_mode` rend
    # cette transformation CONTROLABLE par l'utilisateur - demande
    # explicite utilisateur suite a la decouverte de contributions dont
    # la MAJORITE des dip sont positifs (voir scan_dip_sign/docstring
    # module), signe que cette contribution n'a peut-etre pas enregistre
    # le dip dans la convention MagIC brute. Trois choix (voir
    # ouvrir_convert_magic_to_r_dialog, qui ne suppose jamais silencieusement) :
    # "negate" (standard MagIC, defaut), "raw" (dip .prmag = dip MagIC tel
    # quel), "utrecht" (dip .prmag = 90 - dip MagIC, MEME convention que
    # convert_utrecht_to_r.py pour coreDip - une contribution peut avoir
    # ete uploadee par un logiciel suivant cette convention plutot que la
    # convention MagIC "X depuis l'horizontale").
    azimuth_raw, dip_raw, _from_specimen = _orientation_source(spec_row, sample_row)
    azimuth = (azimuth_raw + 90.0) % 360.0 if azimuth_raw is not None else None
    dip = _DIP_MODES.get(dip_mode, _DIP_MODES["negate"])(dip_raw) if dip_raw is not None else None
    bed_dip_direction = _f(sample_row, "bed_dip_direction")
    bed_dip_strike = (bed_dip_direction - 90.0) % 360.0 if bed_dip_direction is not None else None
    bed_dip = _f(sample_row, "bed_dip")
    # `samples.height` (MagIC data model, "Stratigraphic Height", metres,
    # positif vers le haut - verifie dans pmagpy/data_model/
    # data_model.json) - demande explicite utilisateur ("add an
    # additional variable: stratigraphic_position (or the Magic
    # equivalent)"), reprend directement le nom/l'unite MagIC plutot
    # qu'un nom invente.
    stratigraphic_height = _f(sample_row, "height") if sample_row else None

    formation = _s(site_row, "formation") if site_row else "n.d"
    age_low = _f(site_row, "age_low")
    age_high = _f(site_row, "age_high")
    age_unit = _s(site_row, "age_unit") if site_row else "n.d"
    age = f"{_nd(age_low)} - {_nd(age_high)} {age_unit}" if (age_low is not None or age_high is not None) else "n.d"
    geologic_classes = _s(site_row, "geologic_classes") if site_row else "n.d"
    geologic_types = _s(site_row, "geologic_types") if site_row else "n.d"
    lithologies = _s(site_row, "lithologies") if site_row else "n.d"
    location = _s(loc_row, "location") if loc_row else (_s(site_row, "location") if site_row else "n.d")

    # specimens.txt "volume" (m3) / "weight" (kg), convention MagIC data
    # model - convertis vers les unites .prmag (cm3/g, voir
    # testlect.read_prmag_file : "volume:"/"mass:" sont lues telles
    # quelles par _prmag_nd, sans conversion, donc DEJA en cm3/g ici,
    # comme ech.vol partout ailleurs dans le code - meme facteurs que
    # magic_export.build_specimens_rows au sens inverse (*1e-6/*1e-3).
    # Beaucoup de contributions ne renseignent ni l'un ni l'autre -
    # demande explicite utilisateur ("assume a volume of 10, with warning
    # in the obs field") : volume par defaut avec l'anomalie signalee dans
    # "obs" plutot que silencieusement absorbee comme 1.0 (ce qui
    # fausserait toute intensite/susceptibilite normalisee). Valeur
    # 10.8 cm3 (pas 10.0) - demande explicite utilisateur ("le volume
    # n'est toujours pas correct. Voir les sources en Fortran") : meme
    # defaut que le Fortran d'origine (`pmag(i).vol=10.8`, fichiers.f:1172
    # et fichiers_mod_magic.f:2474 - taille de carotte standard, pas un
    # chiffre rond arbitraire).
    volume_raw = _f(spec_row, "volume") if spec_row else None
    weight_raw = _f(spec_row, "weight") if spec_row else None
    obs = "n.d"
    if volume_raw is not None:
        volume_field, mass_field = f"{volume_raw * 1.0e6:.4f}", "n.d"
    elif weight_raw is not None:
        volume_field, mass_field = "n.d", f"{weight_raw * 1.0e3:.4f}"
    else:
        volume_field, mass_field = "10.8000", "n.d"
        obs = "no volume"

    line_a = (
        f"specimen: {specimen}\tsample: {sample_name}\tsite: {site_name}\t"
        f"volume: {volume_field}\tmass: {mass_field}\t"
        f"lat: {_nd(lat, '.5f')}\tlon: {_nd(lon, '.5f')}\televation: {_nd(elevation, '.1f')}\t"
        f"stratigraphic_height: {_nd(stratigraphic_height, '.2f')}\t"
        f"comment: n.d"
    )
    line_b = (
        f"azimuth: {_nd(azimuth, '.1f')}\tdip: {_nd(dip, '.1f')}\tdate: n.d\t"
        f"magnetic_azimuth: n.d\tsolar_azimuth: n.d\torient_tool: n.d"
    )
    line_c = f"bed_dip_strike: {_nd(bed_dip_strike, '.1f')}\tbed_dip: {_nd(bed_dip, '.1f')}"
    line_d = (
        f"formation: {formation}\tage: {age}\t"
        f"geologic_classes: {geologic_classes}\tgeologic_types: {geologic_types}\t"
        f"lithologies: {lithologies}\tlocation: {location}\tobs: {obs}\t"
        f"method_codes: {_s(spec_row, 'method_codes') if spec_row else 'n.d'}"
    )
    return "\n".join([line_a, line_b, line_c, line_d])


def _experiment_signature(method_codes: str) -> str:
    """Regroupe les mesures d'un specimen par experience plutot que par
    specimen entier - un meme specimen peut porter plusieurs experiences
    distinctes (ex. ATRM (LP-AN-TRM) PUIS paleointensite IZZI) qui ne
    doivent pas partager le meme pool de lettres A/B/C.../demande
    explicite utilisateur, exemple reel HP01-01 : l'ATRM a 600degC est
    faite avant l'IZZI, les deux doivent etre detectees et lettrees
    separement.

    Retourne "SKIP:<protocole>" pour un protocole hors perimetre (voir
    _OUT_OF_SCOPE_PROTOCOLS) - teste avant LP-AN-TRM (LP-AN-ARM/LP-AN-IRM
    partagent le prefixe "LP-AN-" mais ne sont PAS l'ATRM geree ici)."""
    for prefix in _OUT_OF_SCOPE_PROTOCOLS:
        if prefix in method_codes:
            return f"SKIP:{prefix}"
    return "AN-TRM" if "LP-AN-TRM" in method_codes else "PI"


def _experiment_groups(meas_rows: List[Dict[str, str]]):
    """Coupe `meas_rows` en groupes contigus de meme signature
    d'experience (voir _experiment_signature) - preserve l'ordre."""
    groups = []
    cur_sig, cur = None, []
    for row in meas_rows:
        sig = _experiment_signature(row.get("method_codes", ""))
        if cur and sig != cur_sig:
            groups.append((cur_sig, cur))
            cur = []
        cur_sig, cur = sig, cur + [row]
    if cur:
        groups.append((cur_sig, cur))
    return groups


def _atrm_axis_sign(phi: Optional[float], theta: Optional[float]):
    """(axe, signe) a partir de phi/theta pour un pas ATRM (LP-AN-TRM) -
    meme convention que `_measurement_treatment` (magic_export.py:541-547)
    pour les cod1 X/Y/Z natifs Rennes : theta=+/-90 -> Z, sinon phi=0/180
    -> X, phi=90/270 -> Y. Fonctionne aussi bien pour un pas principal que
    pour un controle pTRM ATRM (LT-PTRM-I), qui rapporte le meme phi/theta
    que l'axe qu'il revisite - PAS besoin de "ligne precedente" ici,
    contrairement au controle pTRM IZZI (voir _derive_cod)."""
    def _close(a, b):
        return a is not None and abs(a - b) < 1e-6
    if _close(theta, 90.0):
        return "Z", "+"
    if _close(theta, -90.0):
        return "Z", "-"
    if _close(phi, 0.0):
        return "X", "+"
    if _close(phi, 180.0):
        return "X", "-"
    if _close(phi, 90.0):
        return "Y", "+"
    if _close(phi, 270.0):
        return "Y", "-"
    return "?", "0"


def _assign_izzi_letters(meas_rows: List[Dict[str, str]]) -> Dict[int, str]:
    """id(row) -> lettre A,B,C... Deux strategies de regroupement selon
    le protocole :
    - thermique (LT-T-Z/LT-T-I/LT-PTRM-I) : par VALEUR de `treat_temp`
      (discrete, identique pour Z et I d'un meme palier).
    - micro-onde (LT-M-Z/LT-M-I, protocole LP-PI-M) : par PAIRE
      POSITIONNELLE (compteur, pas par valeur) - `treat_temp` y reste
      fige a 273 sur toutes les lignes (inutilisable), et
      `treat_mw_integral` differe legerement entre Z et I d'un meme
      palier (energie reellement delivree, pas une valeur nominale) - un
      appariement par valeur les separerait a tort en deux groupes.
    Les controles (LT-PTRM-I/LT-PTRM-MD/LT-PMRM-I) ne recoivent PAS de
    lettre ici - la leur vient de `prev_cod2` dans _derive_cod."""
    result: Dict[int, str] = {}
    temp_letters: Dict[str, str] = {}
    mw_pair_count = 0
    for row in meas_rows:
        codes = row.get("method_codes", "")
        if "LT-M-Z" in codes or "LT-M-I" in codes:
            result[id(row)] = chr(ord("A") + mw_pair_count // 2)
            mw_pair_count += 1
        elif ("LT-T-Z" in codes or "LT-T-I" in codes) and "LP-PI-" in codes:
            # letter uniquement pour la paleointensite (S/R IZZI) - un
            # simple demag thermique directionnel (LP-DIR-T, pas de champ
            # applique, pas de paire Z/I a distinguer) n'en a pas besoin -
            # demande explicite utilisateur.
            temp = row.get("treat_temp", "")
            if temp not in temp_letters:
                temp_letters[temp] = chr(ord("A") + len(temp_letters))
            result[id(row)] = temp_letters[temp]
    return result


def _derive_cod(row: Dict[str, str], letters: Dict[int, str], prev_cod2: str):
    """(cod1, cod2, step_value) a partir de method_codes/treat_* - pas de
    cod1/cod2 natif dans une source MagIC, invente ici pour garder la
    meme forme step/cod1/cod2 que les fichiers convertis depuis .ren.

    Pour un controle pTRM (P/M), `cod2` reprend celui de la mesure
    PRECEDENTE dans la sequence (`prev_cod2`), PAS un lookup par sa
    propre temperature (celle-ci est souvent une temperature deja passee,
    revisitee pour le controle) - precision explicite de l'utilisateur :
    la lettre indique a quel point de la sequence (quel groupe R/S en
    cours) le controle a ete fait, pas quelle temperature il verifie
    (deja portee par `step`)."""
    codes = row.get("method_codes", "")
    temp = _temp_k(row)
    letter = letters.get(id(row), "")

    if "LP-AN-TRM" in codes:
        # ATRM : baseline zero-field -> 'D' (code Rennes deja existant,
        # magic_export.py cod1=='D'), pas en champ -> axe/signe X/Y/Z
        # +/- derives du phi/theta DE CETTE LIGNE (main step ou controle
        # pTRM, meme regle - voir _atrm_axis_sign).
        if "LT-T-Z" in codes:
            return "D", "0", (temp - 273.0) if temp is not None else 0.0
        axis, sign = _atrm_axis_sign(_f(row, "treat_dc_field_phi"), _f(row, "treat_dc_field_theta"))
        return axis, sign, (temp - 273.0) if temp is not None else 0.0

    if "LT-NO" in codes:
        return "N", "0", 0.0
    if "LT-T-Z" in codes:
        if "LP-PI-" in codes:
            return "S", letter or "0", (temp - 273.0) if temp is not None else 0.0
        # demag thermique directionnel simple (ex. LP-DIR-T) : pas de
        # paleointensite, pas de paire Z/I -> 'D' sans lettre - demande
        # explicite utilisateur ("do not use S for thermal demag").
        return "D", "0", (temp - 273.0) if temp is not None else 0.0
    if "LT-T-I" in codes:
        if "LP-PI-II" in codes:
            # BUG CORRIGE (signale par l'utilisateur - "there is still a
            # problem in the calculation of the PTRM check") : LP-PI-II
            # ("Original Thellier-Thellier method", verifie sur le
            # vocabulaire en ligne) est le protocole "champ en Z puis Z-"
            # (deja documente/gere cote paleointensity_magic.py,
            # build_magic_dataframe - "there is no choice than (R+V)/2",
            # confirme utilisateur sur le specimen "02B") : DEUX mesures
            # EN CHAMP par etape, orientations opposees (theta=+90 puis
            # -90), PAS de pas zero-field distinct (aucun LT-T-Z dans la
            # source, verifie sur magic_contribution_19987.txt,
            # kr01_02B3). Avant ce fix, les DEUX mesures d'un meme etape
            # recevaient cod1='R' (le code ne distinguait jamais Z+ de
            # Z-), produisant deux lignes cod1/cod2/etape IDENTIQUES -
            # aucune n'etait alors reconnaissable comme 'V', donc le
            # mecanisme (R+V)/2 de build_magic_dataframe ne se declenchait
            # JAMAIS et `pmag.sortarai` ne trouvait aucun pas zero-field
            # (echec total : "not enough IZZI/Thellier steps found (1)").
            # Ici : theta=+90 (Z+, la reference) -> 'R', theta=-90 (Z-) ->
            # 'V' - meme convention que magic_export.py (sens export) et
            # build_magic_dataframe (sens analyse).
            theta = _f(row, "treat_dc_field_theta")
            cod1 = "V" if (theta is not None and theta < 0) else "R"
            return cod1, letter or "0", (temp - 273.0) if temp is not None else 0.0
        return "R", letter or "0", (temp - 273.0) if temp is not None else 0.0
    if "LT-PTRM-I" in codes or "LT-PTRM-Z" in codes:
        # LT-PTRM-I ("after zero field step, in field cooling", controle
        # classique) et LT-PTRM-Z ("after in-field step, zero field
        # cooling at a lower temperature", pTRM tail check a temperature
        # plus basse) sont TOUS LES DEUX un pas 'P' cote Rennes - le
        # format natif STARpaleomag_Py ne distingue pas les deux variantes par
        # cod1 (voir magic_export.py cod1=='P', qui les distingue plutot
        # via dc_field selon que le pas precedent est 'S' ou 'R' - demande
        # explicite utilisateur, "there is two ways of doing this ...").
        # Confondre les deux ici (avant ce fix, LT-PTRM-Z tombait dans le
        # cod1='?' non reconnu ci-dessous) romprait l'import de donnees
        # IZZI reelles telechargees depuis MagIC utilisant le controle en
        # champ nul.
        return "P", prev_cod2 or "0", (temp - 273.0) if temp is not None else 0.0
    if "LT-PTRM-MD" in codes:
        # pTRM-tail check (MD = multidomain tail check, protocole IZZI) -
        # aucun code Rennes existant pour ce pas, 'M' invente ici.
        return "M", prev_cod2 or "0", (temp - 273.0) if temp is not None else 0.0
    if "LT-M-Z" in codes:
        # Paleointensite micro-onde (LP-PI-M) : meme mecanique IZZI que
        # LT-T-Z/LT-T-I/LT-PTRM-I (zero-field/in-field/controle), reprend
        # donc les memes codes S/R/P - `treat_temp` reste fige a 273 pour
        # TOUTES les lignes micro-onde (pas de vrai palier thermique),
        # `step` est ici l'integrale micro-onde (treat_mw_integral),
        # seule valeur qui croit reellement d'un palier au suivant.
        mw = _f(row, "treat_mw_integral")
        return "S", letter or "0", mw if mw is not None else 0.0
    if "LT-M-I" in codes:
        mw = _f(row, "treat_mw_integral")
        return "R", letter or "0", mw if mw is not None else 0.0
    if "LT-PMRM-I" in codes:
        mw = _f(row, "treat_mw_integral")
        return "P", prev_cod2 or "0", mw if mw is not None else 0.0
    if "LT-AF-Z" in codes:
        af = _ac_field_t(row)
        return "F", "0", (af * 1.0e3 if af else 0.0)
    if "LT-AF-I" in codes:
        af = _ac_field_t(row)
        return "A", "0", (af * 1.0e3 if af else 0.0)
    if "LT-IRM" in codes:
        dc = _f(row, "treat_dc_field")
        return "I", "0", (dc * 1.0e3 if dc else 0.0)
    # BUG UTILISATEUR REEL dans certaines contributions MagIC (signale par
    # l'utilisateur - "in some Magic contribution, the users did not
    # write well the code, they use LP instead of LT") : LP-DIR-T/
    # LP-DIR-AF sont des codes de PROTOCOLE (niveau experience -
    # "Directional data: Step-wise thermal/AF demagnetization", verifie
    # sur le vocabulaire MagIC en ligne), PAS des codes de TRAITEMENT PAR
    # PALIER (LT-T-Z/LT-AF-Z) - certains contributeurs les posent quand
    # meme SEULS sur chaque ligne de mesure, sans jamais indiquer le bon
    # LT-. Repli, UNIQUEMENT si aucun LT- reconnu ci-dessus n'a matche :
    # deduire tout de meme le type de demagnetisation plutot que de
    # laisser tomber ces pas en cod1='?' (jamais exploitable ensuite,
    # Zijderveld/ajuslig les ignorent). LP-DIR-T/LP-DIR-AF designent PAR
    # DEFINITION une demagnetisation directionnelle simple, donc TOUJOURS
    # en champ nul (comme LT-T-Z/LT-AF-Z, jamais LT-T-I/LT-AF-I - aucune
    # ambiguite possible ici).
    if "LP-DIR-T" in codes and temp is not None:
        return "D", "0", temp - 273.0
    if "LP-DIR-AF" in codes:
        af = _ac_field_t(row)
        if af is not None:
            return "F", "0", af * 1.0e3
    # non reconnu : conserve quand meme la temperature si presente
    return "?", "0", (temp - 273.0) if temp is not None else 0.0


_DC_LOWFIELD_COD1 = {"R", "P", "X", "Y", "Z"}


def _measurement_rows(specimen: str, meas_rows: List[Dict[str, str]], problems: list) -> Optional[List[str]]:
    """Construit les lignes de mesure - les groupes d'experience HORS
    PERIMETRE (voir _OUT_OF_SCOPE_PROTOCOLS) sont ECARTES ici (pas de
    ligne ecrite pour eux) et journalises dans `problems` (tuples
    (reason, count, specimen)) plutot que de faire echouer la conversion
    - demande explicite utilisateur ("provide a list of the problems
    encountered, but proceed with the demagnetization data")."""
    rows = [_MEAS_HEADER]
    j = 0
    unrecognized_codes: Dict[str, int] = {}
    for sig, group in _experiment_groups(meas_rows):
        if sig.startswith("SKIP:"):
            prefix = sig[len("SKIP:"):]
            label = _OUT_OF_SCOPE_PROTOCOLS.get(prefix, prefix)
            problems.append((f"{label} ({prefix}) - not processed by STARpaleomag_Py", len(group), specimen))
            continue

        # letters/prev_cod2 scopes a l'experience (voir _experiment_groups)
        # - une nouvelle experience sur le meme specimen (ex. ATRM puis
        # IZZI) ne doit pas hereter du pool de lettres de la precedente.
        letters = _assign_izzi_letters(group)
        prev_cod2 = "0"
        for row in group:
            cod1, cod2, step_val = _derive_cod(row, letters, prev_cod2)
            prev_cod2 = cod2
            if cod1 == "?":
                key = (row.get("method_codes", "") or "(empty method_codes)").strip()
                unrecognized_codes[key] = unrecognized_codes.get(key, 0) + 1

            dec = _f(row, "dir_dec")
            inc = _f(row, "dir_inc")
            mag = _f(row, "magn_moment")
            x, y, z = _dir2cart(dec, inc, mag)
            if x is None:
                x, y, z = _f(row, "magn_x"), _f(row, "magn_y"), _f(row, "magn_z")

            # _ac_field_t/_temp_k (repli meas_field_ac/meas_temp pour les
            # fichiers non conventionnels, voir leur docstring) reutilises
            # ICI AUSSI - pas seulement dans _derive_cod - pour que les
            # colonnes treat_ac_field/treat_temp ECRITES dans le .r
            # restent coherentes avec le cod1/etape deja derive plus haut
            # (sinon un pas recupere via meas_temp afficherait un etape
            # correct mais une colonne treat_temp vide/a 0, incoherence
            # visible dans le fichier converti).
            af = _ac_field_t(row)
            dc = _f(row, "treat_dc_field")
            af_mT = af * 1.0e3 if (af and cod1 in ("A", "F")) else None
            dc_strong = dc * 1.0e3 if (dc and cod1 == "I") else None
            dc_low = dc * 1.0e6 if (dc and cod1 in _DC_LOWFIELD_COD1) else 0.0

            phi = _f(row, "treat_dc_field_phi")
            theta = _f(row, "treat_dc_field_theta")
            temp = _temp_k(row)
            temp_c = (temp - 273.0) if temp is not None else 0.0

            ins_full = _s(row, "instrument_codes")
            codes = row.get("method_codes", "n.d")
            s_val = None  # pas de susceptibilite dans ces contributions IZZI

            line = "\t".join([
                f"{step_val:6.1f}", cod1, cod2,
                _nd(x, "11.3E"), _nd(y, "11.3E"), _nd(z, "11.3E"),
                "n.d", _s(row, "quality"), "n.d", _nd(s_val, "7.1f"),
                _nd(temp_c, "6.1f"),
                _nd(af_mT, "7.2f"),
                _nd(dc_strong, "7.2f"),
                _nd(dc_low, "7.2f"),
                _nd(phi, "6.1f"), _nd(theta, "6.1f"),
                f"{codes:<40}", ins_full, str(j + 1),
            ])
            rows.append(line)
            j += 1

    for codes_str, count in unrecognized_codes.items():
        problems.append((f"unrecognized method_codes: {codes_str}", count, specimen))

    if len(rows) == 1:
        # tous les groupes de ce specimen etaient hors perimetre - rien
        # a ecrire (deja journalise ci-dessus, un par groupe).
        return None
    return rows


def convert_magic_file(
    path_in: str, path_out: str, dip_mode: str = "negate",
) -> Tuple[int, str, int, int, int, Optional[str], int]:
    tables = parse_magic_contribution(path_in)
    locations = _index_by(tables.get("locations", []), "location")
    sites = _index_by(tables.get("sites", []), "site")
    samples = _index_by(tables.get("samples", []), "sample")
    specimens = _index_by(tables.get("specimens", []), "specimen")

    by_specimen: Dict[str, List[Dict[str, str]]] = {}
    for row in tables.get("measurements", []):
        by_specimen.setdefault(row.get("specimen", ""), []).append(row)

    lines = [
        FORMAT_HEADER,
        f"#converted from {os.path.basename(path_in)} by convert_magic_to_r.py "
        "(MagIC contribution -> .r, interoperability test)",
        S_FIELD_HEADER,
        "",
    ]
    blocks = []
    problems: List[str] = []
    for specimen, meas_rows in by_specimen.items():
        spec_row = specimens.get(specimen)
        sample_row = samples.get(spec_row.get("sample", "")) if spec_row else None
        site_row = sites.get(sample_row.get("site", "")) if sample_row else None
        loc_row = locations.get(site_row.get("location", "")) if site_row else None

        meas_lines = _measurement_rows(specimen, meas_rows, problems)
        if meas_lines is None:
            continue  # tout le specimen etait hors perimetre (voir problems)

        header_block = _sample_header_block(
            specimen, spec_row, sample_row, site_row, loc_row, dip_mode=dip_mode)
        blocks.append(header_block + "\n" + "\n".join(meas_lines))

    with open(path_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        f.write("\n\n".join(blocks) + "\n")

    report = format_problems_report(path_in, problems)
    nb_results, nb_means = _convert_magic_results(path_in, path_out)
    nb_pint, redo_pint_path = _write_pint_redo_file(path_in, path_out)
    nb_aniso = _write_ani_tensors_from_magic(path_in, path_out)
    return len(blocks), report, nb_results, nb_means, nb_pint, redo_pint_path, nb_aniso


def _write_pint_redo_file(path_in: str, path_out: str) -> Tuple[int, Optional[str]]:
    """Ecrit un fichier "redo_pint_contribution_NUM.txt" (une ligne
    "specimen tmin tmax" par determination, voir
    magic_pint_results_to_redo_lines) a cote de `path_out`, directement
    utilisable par "View Paleoint Results..." (ouvrir_openfilepint_dialog)
    - demande explicite utilisateur ("during import of Magic paleoint
    data, can you write a redo_pint_contribution_num.txt file with the
    specimen number and step1 and step2 found in specimens.txt file. This
    redo file can be used in view Paleointensity data"). NUM est le
    numero de contribution extrait du nom du fichier source
    (magic_contribution_NUM.txt, cf. .gitignore) ; a defaut (fichier
    renomme), reprend le nom de base de `path_out`. Retourne (0, None)
    sans creer de fichier si aucune ligne exploitable (specimens.txt
    absent, ou aucune determination LP-PI-TRM de qualite 'g')."""
    lines = magic_pint_results_to_redo_lines(path_in, combined=True)
    if not lines:
        return 0, None
    match = re.search(r"magic_contribution_(\w+)", os.path.basename(path_in))
    num = match.group(1) if match else os.path.splitext(os.path.basename(path_out))[0]
    redo_path = os.path.join(os.path.dirname(path_out) or ".", f"redo_pint_contribution_{num}.txt")
    with open(redo_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return len(lines), redo_path


def _write_ani_tensors_from_magic(path_in: str, path_out: str) -> int:
    """Archive dans le .pmagani compagnon de `path_out` les tenseurs
    d'anisotropie DEJA calcules par la contribution MagIC source (voir
    extract_magic.magic_anisotropy_rows) - demande explicite utilisateur
    ("lors de l'importation des fichiers Magic, archiver dans le fichier
    pmagani les donnees de tenseurs d'anisotropie avec le code A0 si
    c'est un tenseur de TRM, N0 pour l'AMS, AA pour l'ARM"). Retourne le
    nombre de lignes ecrites (0 si specimens.txt n'a aucune ligne
    d'anisotropie exploitable - pas une erreur, la plupart des
    contributions n'en ont pas)."""
    rows = magic_anisotropy_rows(path_in, combined=True)
    if not rows:
        return 0
    return write_ani_tensors_from_magic_rows(ani_path_for(path_out), rows)


def _convert_magic_results(path_in: str, path_out: str) -> Tuple[int, int]:
    """Genere AUSSI le .pmagres compagnon a partir de la MEME contribution
    MagIC - demande explicite utilisateur ("is the pmagres file also
    generated during the magic import") : jusqu'ici seul convert_ren_to_r.py
    le faisait (pour un vieux fichier .r legacy compagnon), pas cet
    importeur MagIC. Deux sources independantes dans sites.txt/
    specimens.txt (les dialogues dedies app.py qui faisaient la meme chose
    separement - import_magic_results_direct_dialog/
    import_magic_site_means_dialog/import_magic_results_dialog - ont ete
    retires a la demande explicite de l'utilisateur, cette generation
    automatique les rendant redondants pour le cas d'usage normal) :

    1) specimens.txt deja interprete (dir_dec/dir_inc/method_codes DE-BFL/
       DE-BFP calcules par pmagpy/demag_gui) -> RECALCULE via le fit natif
       STARpaleomag_Py (calcul.fit_from_redo_file) sur les mesures qu'on vient
       d'ecrire dans path_out (relues via testlect.read_prmag_file, pas
       recopiees depuis MagIC) - memes raisons qu'avant :
       rester coherent avec ce que STARpaleomag_Py calculerait lui-meme.
    2) sites.txt (dir_dec/dir_inc/dir_alpha95/dir_k/vgp_lat/vgp_lon deja
       calcules) -> archive TEL QUEL (rien a recalculer, STARpaleomag_Py n'a pas
       les directions specimen individuelles utilisees pour ce calcul),
       `liste` cross-reference les specimens de la colonne "specimens"
       contre les resultats FRAICHEMENT archives a l'etape 1 ci-dessus
       (memes id, meme session de conversion).

    Retourne (0, 0) sans erreur si aucune des deux tables ne porte de
    resultat exploitable - un fichier de mesures pures reste un cas normal,
    pas une erreur."""
    path_results = results_path_for(path_out)
    if os.path.exists(path_results):
        # meme raison que convert_ren_to_r.py/convert_utrecht_to_r.py :
        # une re-conversion vers le MEME chemin de sortie dupliquerait les
        # resultats si on ne repartait pas d'un fichier vide (archivres
        # ajoute a un fichier existant, ne l'ecrase pas).
        os.remove(path_results)

    donnees = read_prmag_file(path_out)
    existing_ids: Optional[set] = None
    nb_results = 0

    lines = magic_results_to_redo_lines(path_in, combined=True)
    specimen_fits: Dict[str, FitResult] = {}
    if lines:
        for fit in fit_from_redo_file(donnees, lines):
            _, existing_ids = archivres(fit, path_results, existing_ids)
            specimen_fits[fit.id] = fit
            nb_results += 1

    nb_means = 0
    site_rows = magic_site_means(path_in, combined=True)
    for row in (site_rows or []):
        cs = [specimen_fits[s].c for s in row["specimens"] if s in specimen_fits]
        # Ovale de confiance du VGP (dp/dm) - demande explicite utilisateur
        # ("during the archive of the mean direction in pmagres, the VGP
        # is calculated. Is it possible to add the dp dm ellipse of
        # confidence derived from the a95, available in Stereo_Py") :
        # reprend dp/dm PUBLIES par la contribution MagIC elle-meme s'ils
        # existent (plus fiable, ex. issu d'un bootstrap - voir
        # magic_site_means), sinon les derive de alpha95/inc (verifie
        # contre un exemple reel publiant les deux, calcul.dp_dm_from_a95).
        if row["vgp_dp"] or row["vgp_dm"]:
            vgp_dp, vgp_dm = row["vgp_dp"], row["vgp_dm"]
        else:
            vgp_dp, vgp_dm = dp_dm_from_a95(row["alpha95"], row["inc"])
        mean_fit = FitResult(
            id=f"mean: {row['site']}", cat1="F", cat2="i",
            nb=row["n"], dec=row["dec"], inc=row["inc"], mad=row["alpha95"],
            tx=(row["k"], 0.0), par3_mean=row["orientation"],
            lat=row["lat"], rlong=row["lon"],
            par4=row["vgp_lat"], par5=row["vgp_lon"],
            vgp_dp=vgp_dp, vgp_dm=vgp_dm,
            # Etiquette de composante de magnetisation (dir_comp_name,
            # "A" par defaut) - voir calcul.FitResult.component et
            # extract_magic.magic_site_means.
            component=row["component"],
            liste="codes:" + ":".join(str(c) for c in cs),
        )
        _, existing_ids = archivres(mean_fit, path_results, existing_ids)
        nb_means += 1

    return nb_results, nb_means


def format_problems_report(path_in: str, problems: list) -> str:
    """Rapport de fin de conversion (demande explicite utilisateur :
    "provide a list of the problems encountered") - regroupe par raison
    (protocole hors perimetre / method_codes non reconnu) plutot que de
    lister chaque specimen individuellement (fichiers reels : jusqu'a
    ~900 specimens, une ligne par specimen serait illisible). Retourne le
    texte plutot que de l'imprimer directement, pour etre reutilisable
    tel quel par l'interface graphique (app.py)."""
    if not problems:
        return f"{os.path.basename(path_in)}: no problems encountered."
    by_reason: Dict[str, Dict[str, int]] = {}
    for reason, count, specimen in problems:
        by_reason.setdefault(reason, {})[specimen] = by_reason.get(reason, {}).get(specimen, 0) + count
    lines_out = [f"{os.path.basename(path_in)}: {len(problems)} problem group(s) found -"]
    for reason, per_specimen in by_reason.items():
        total_rows = sum(per_specimen.values())
        lines_out.append(f"  - {reason}: {total_rows} row(s) across {len(per_specimen)} specimen(s)")
    return "\n".join(lines_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a MagIC contribution to the .r format")
    parser.add_argument("magic_file")
    parser.add_argument("-o", "--output")
    parser.add_argument(
        "--dip-mode", choices=("negate", "raw", "utrecht"), default="negate",
        help="How to convert MagIC dip to .prmag dip: 'negate' (standard MagIC "
             "convention, default), 'raw' (dip .prmag = dip MagIC as-is), 'utrecht' "
             "(dip .prmag = 90 - dip MagIC, Utrecht/PMAG2 convention) - see scan_dip_sign")
    args = parser.parse_args()
    out = args.output or os.path.splitext(args.magic_file)[0] + ".prmag"
    dip_scan = scan_dip_sign(args.magic_file)
    if dip_scan["majority_positive"]:
        pct = 100.0 * dip_scan["n_positive"] / dip_scan["n_total"]
        print(
            f"WARNING: {dip_scan['n_positive']}/{dip_scan['n_total']} dip value(s) "
            f"({pct:.0f}%) in this file are positive ({dip_scan['n_from_specimen']} "
            "from specimen-level orientation) - in practice raw MagIC dip should only "
            "ever be negative or zero. Check --dip-mode raw/utrecht if this "
            "contribution did not use the standard MagIC convention.\n"
        )
    n, report, nb_results, nb_means, nb_pint, redo_pint_path, nb_aniso = convert_magic_file(
        args.magic_file, out, dip_mode=args.dip_mode)
    print(report)
    print(f"{n} specimen(s) converted -> {out}")
    print(f"{nb_results} result(s) and {nb_means} site mean(s) -> {results_path_for(out)}")
    if redo_pint_path:
        print(f"{nb_pint} paleointensity determination(s) -> {redo_pint_path}")
    if nb_aniso:
        print(f"{nb_aniso} anisotropy tensor(s) -> {ani_path_for(out)}")
