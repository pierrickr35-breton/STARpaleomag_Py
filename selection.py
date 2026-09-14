"""
Equivalent Python des subroutines Fortran `selmes` et `lismes` (select.f.txt),
plus les rotations `corfor`, `corpen`, `polere` (portees depuis calcul.f, le
seul endroit du depot ou leur code source a ete retrouve - starmac_OSX.inc ne
contient que les structures de donnees, pas ces subroutines).

Travaille sur les objets Pmag / Measurement produits par testlect.read_ren_file().
"""

import math
import re
import sys
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, TextIO, Tuple

from testlect import Measurement, Pmag
from orient_sample import compute_local_declination

# Voir app.HEADER_MARK - marque une ligne de titres de colonnes pour un
# affichage en gras cote console, meme valeur, definie localement pour
# eviter un import circulaire depuis app.py.
_HEADER_MARK = "\x01"


# ---------------------------------------------------------------------------
# Equivalent de la structure /Echantillons/ ech(3000) (starmac_OSX.inc)
# ---------------------------------------------------------------------------

@dataclass
class SelectedSample:
    """Un echantillon selectionne, avec uniquement ses mesures retenues."""
    id: str = ""
    cin: float = 0.0
    caz: float = 0.0
    dip: float = 0.0
    str_: float = 0.0
    vol: float = 0.0
    norme: str = ""
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
    magic_site: str = ""
    magic_sample: str = ""
    magic_fm: str = ""
    magic_age: str = ""
    magic_gc: str = ""
    magic_smt: str = ""
    magic_li: str = ""
    magic_loc: str = ""
    magic_obs: str = ""
    # Position stratigraphique (metres, positif vers le haut) - voir
    # testlect.Pmag.stratigraphic_height ; None si non renseignee.
    stratigraphic_height: Optional[float] = None
    mesures: List[Measurement] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Segmentation en experiences distinctes - PAS dans le Fortran (un
# echantillon y ne suit toujours qu'un seul protocole) : un meme
# echantillon peut enchainer plusieurs protocoles dans sa liste de
# mesures (ex. desaimantation AF, puis acquisition IRM, puis
# desaimantation thermique de cette IRM acquise - cas reel signale par
# l'utilisateur, ex. 19DN1607B). Chaque nouvelle experience redemarre son
# etape a une valeur STRICTEMENT plus basse que la derniere etape de la
# precedente - signal fiable confirme par l'utilisateur.
# ---------------------------------------------------------------------------

def split_experiments(mesures: List[Measurement]) -> List[List[Measurement]]:
    """Decoupe `mesures` en runs d'experiences distinctes (voir commentaire
    de section) - retourne toujours au moins 1 run (jamais vide) si
    `mesures` est non vide, dans l'ordre du fichier."""
    if not mesures:
        return []
    runs = [[mesures[0]]]
    for prev, cur in zip(mesures, mesures[1:]):
        if cur.etape < prev.etape:
            runs.append([])
        runs[-1].append(cur)
    return runs


def experiment_kind(run: List[Measurement]) -> str:
    """Type d'experience d'un run (voir split_experiments) : 'F' (AF),
    'D' (thermique de la NRM), 'D_IRM' (thermique d'une IRM deja acquise -
    cod2='I' sur les etapes D/S/T/K, confirme par l'utilisateur : "0DI",
    "130DI"... marquent explicitement la desaimantation thermique de
    l'IRM, pas de la NRM), 'I' (acquisition IRM), ou le code cod1 brut si
    non reconnu. Le premier point d'un run est presque toujours 'N' (NRM
    initiale), pas representatif - meme logique que
    xygraph._sample_demag_code (cherche le premier code reel)."""
    for m in run:
        if m.cod1 == "N":
            continue
        if m.cod1 in ("A", "F"):
            return "F"
        if m.cod1 in ("D", "S", "T", "K"):
            return "D_IRM" if m.cod2 == "I" else "D"
        if m.cod1 == "I":
            return "I"
        return m.cod1
    return "N"


def zijderveld_measurements(ech) -> List[Measurement]:
    """Memes mesures que zijderveld.draw_zijderveld affiche : concatene
    les runs (split_experiments) depuis le debut, arrete au premier run
    IRM/thermique-d'IRM (experiment_kind in ('I','D_IRM')) - partage entre
    zijderveld.py, auto_interpretation.py et interpretation_quality.py
    (ici plutot que dans l'un des deux pour eviter un import circulaire
    entre ces deux derniers, qui s'importent deja l'un l'autre)."""
    runs = split_experiments(ech.mesures)
    mesures: List[Measurement] = []
    for run in runs:
        if experiment_kind(run) in ("I", "D_IRM"):
            break
        mesures.extend(run)
    return mesures


# ---------------------------------------------------------------------------
# Rotations geometriques (portees depuis Appli_Paleomag/2G_Magneto/calcul.f)
# ---------------------------------------------------------------------------

def polere(x: float, y: float, z: float) -> tuple:
    """Cartesien (x,y,z) -> (magnitude, declinaison, inclinaison) en degres."""
    horiz = math.hypot(x, y)
    mag = math.hypot(horiz, z)
    if mag == 0.0:
        return 0.0, 0.0, 0.0
    dec = math.degrees(math.atan2(y, x))
    if dec < 0.0:
        dec += 360.0
    inc = math.degrees(math.atan2(z, horiz))
    return mag, dec, inc


def corfor(n: float, e: float, vv: float, g: float, h: float) -> tuple:
    """Correction de forage : repere echantillon -> repere in-situ, a partir
    du pendage (g) et de l'azimut (h) de l'axe de la carotte."""
    gr, hr = math.radians(g), math.radians(h)
    xx = n * math.cos(gr) * math.sin(hr) + e * math.cos(hr) + vv * math.sin(hr) * math.sin(gr)
    yy = -n * math.cos(gr) * math.cos(hr) + e * math.sin(hr) - vv * math.cos(hr) * math.sin(gr)
    zz = -n * math.sin(gr) + vv * math.cos(gr)
    return xx, yy, zz


def corpen(x: float, y: float, z: float, j: float, k: float) -> tuple:
    """Correction de pendage (tectonique) : rotation autour de la direction
    de la strike (k) de l'angle de pendage (j)."""
    jr, kr = math.radians(j), math.radians(k)
    cj, sj = math.cos(jr), math.sin(jr)
    ck, sk = math.cos(kr), math.sin(kr)
    xx = x * (ck ** 2 + sk ** 2 * cj) + y * (ck * sk * (1.0 - cj)) - z * sk * sj
    yy = x * ck * sk * (1.0 - cj) + y * (sk ** 2 + cj * ck ** 2) + z * sj * ck
    zz = x * sj * sk - y * sj * ck + z * cj
    return xx, yy, zz


def angle(p: float, a: float, b: float) -> tuple:
    """Distance angulaire (deg) entre deux directions sur la sphere, a partir
    de leur difference d'azimut `p` et de leurs deux colatitudes `a`, `b`
    (tous en degres). Retourne (del, ang) ; seul `del` (l'angle) est utilise
    par lismesdepth. Portee depuis calcul.f (SUBROUTINE angle)."""
    pr, ar, br = math.radians(p), math.radians(a), math.radians(b)
    cs = math.cos(ar) * math.cos(br) + math.sin(ar) * math.sin(br) * math.cos(pr)
    cs = min(cs, 1.0)
    y = math.sqrt(max(0.0, 1.0 - cs * cs))
    delr = math.atan2(y, cs)
    dell = math.sin(delr) * math.sin(br)
    if dell == 0.0:
        dell = 0.000001
    cs2 = (math.cos(ar) - math.cos(delr) * math.cos(br)) / dell
    cs2 = min(cs2, 1.0)
    y2 = math.sqrt(max(0.0, 1.0 - cs2 * cs2))
    ang = math.degrees(math.atan2(y2, cs2))
    return math.degrees(delr), ang


def apply_orientation(x: float, y: float, z: float, ech: "SelectedSample", orientation: int) -> tuple:
    """Applique la correction d'orientation demandee (1=echantillon,
    2=in-situ, 3=apres correction de pendage)."""
    if orientation == 1:
        return x, y, z
    if orientation == 2:
        return corfor(x, y, z, ech.cin, ech.caz)
    if orientation == 3:
        x2, y2, z2 = corfor(x, y, z, ech.cin, ech.caz)
        return corpen(x2, y2, z2, ech.dip, ech.str_)
    raise ValueError(f"invalid orientation: {orientation}")


def _s_value(ech: "SelectedSample", m: Measurement) -> float:
    if ech.norme == "m":
        return m.s * 1.0e-7 / ech.vol if ech.vol else 0.0
    return m.s * 1.0e-4 / ech.vol if ech.vol else 0.0


# cod1 -> unite d'etape affichable - meme convention que xygraph._AF_CODES/
# _THERMAL_CODES (etape = valeur physique reelle directement, mT pour A/F,
# degC pour D/S/T/K) - R/V/P (paleointensite Thellier/IZZI) INCLUS ici,
# egalement en degC - demande explicite utilisateur ("when the cod1 is R,
# V, P, the unit is also °C"), confirme apres une premiere version qui les
# excluait par prudence (unite alors non confirmee). PAS etendue aux
# autres codes (I, N...) qui n'ont pas d'unite de temperature/champ
# univoque.
_AF_CODES = {"A", "F"}
_THERMAL_CODES = {"D", "S", "T", "K", "R", "V", "P"}


def _fmt_step(m: "Measurement") -> str:
    """Affiche l'etape avec son unite quand elle est sans ambiguite (degC
    pour un pas thermique, mT pour un pas AF - etape est directement la
    valeur physique reelle, plus une echelle Oersted a diviser - voir
    testlect.Measurement) suivie du code - demande explicite utilisateur
    ("can we clean the list data... 210D+ -> 210 degC D+" puis "convert
    all step integer to float")."""
    if m.cod1 in _AF_CODES:
        return f"{m.etape:5.1f} mT {m.cod1}{m.cod2}"
    if m.cod1 in _THERMAL_CODES:
        return f"{m.etape:5.0f} °C {m.cod1}{m.cod2}"
    return f"{m.etape:5.0f}    {m.cod1}{m.cod2}"


def _selection_norme_labels(selected: List["SelectedSample"]) -> Tuple[str, str]:
    """(libelle colonne intensite, libelle colonne susceptibilite) -
    determines depuis la normalisation (volume/masse) REELLEMENT utilisee
    par la selection, plutot que d'afficher systematiquement les deux
    unites cote a cote avec un "--" pour celle qui ne s'applique pas -
    demande explicite utilisateur ("remove Am2/kg if normalization by
    volume... susceptibility: K (with its unit following the type of
    normalization)"). Si la selection melange volume ET masse (rare, mais
    possible), retombe sur le double libelle - la seule situation ou les
    deux unites peuvent legitimement apparaitre dans la meme liste."""
    normes = {ech.norme for ech in selected if ech.vol}
    if normes == {"v"}:
        return "A/m", "K (SI)"
    if normes == {"m"}:
        return "Am2/kg", "K (m3/kg)"
    return "A/m or Am2/kg", "K (SI or m3/kg)"


def normalized_intensity(ech: "SelectedSample", raw_magnitude: float) -> Tuple[float, str]:
    """(valeur, unite) - normalise `raw_magnitude` (Am2, deja la norme
    d'un vecteur brut) par le volume/masse de `ech`, comme partout
    ailleurs (facteur 1e6 pour un volume -> A/m, 1e3 pour une masse ->
    Am2/kg). Si `ech.vol` est manquant (0/None), retourne `raw_magnitude`
    INCHANGE avec l'unite "Am2" (moment total brut, PAS de facteur
    applique) - demande explicite utilisateur ("when the volume or the
    mass of the sample is not given, best to show the data in total
    moment as done in the list... the previous Fortran was systematically
    dividing by mass and volume and was not expecting no vol and no
    mass") : le Fortran d'origine appliquait quand meme le facteur 1e3/
    1e6 a un volume/masse absent (repli implicite sur 1.0), affichant une
    valeur sans rapport avec l'unite revendiquee plutot que le moment
    total honnete - meme principe que zijderveld._scale_factor."""
    if not ech.vol:
        return raw_magnitude, "Am2"
    if ech.norme == "m":
        return raw_magnitude * 1.0e3 / ech.vol, "Am2/kg"
    return raw_magnitude * 1.0e6 / ech.vol, "A/m"


def _parse_heuremes(s: Optional[str]) -> Optional[datetime]:
    if not s:
        return None
    m = re.match(r"(\d{4})\D(\d{1,2})\D(\d{1,2})\D(\d{1,2})\D(\d{1,2})\D(\d{1,2})", s)
    if not m:
        return None
    year, month, day, hour, minute, second = (int(x) for x in m.groups())
    try:
        return datetime(year, month, day, hour, minute, second)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# selmes : selection des echantillons / mesures
# ---------------------------------------------------------------------------

def _build_selected_sample(p: Pmag, matched_mesures: List[Measurement]) -> SelectedSample:
    return SelectedSample(
        id=p.id,
        cin=p.cin,
        caz=p.caz,
        dip=p.dip,
        str_=p.str_,
        vol=p.vol,
        norme=p.norme,
        com=p.com,
        lat=p.lat,
        rlong=p.rlong,
        altitude=p.altitude,
        year=p.year,
        month=p.month,
        day=p.day,
        hour=p.hour,
        minute=p.minute,
        azmag=p.azmag,
        azsun=p.azsun,
        outilorient=p.outilorient,
        roche=p.roche,
        magic_site=p.magic_site,
        magic_sample=p.magic_sample,
        magic_fm=p.magic_fm,
        magic_age=p.magic_age,
        magic_gc=p.magic_gc,
        magic_smt=p.magic_smt,
        magic_li=p.magic_li,
        magic_loc=p.magic_loc,
        magic_obs=p.magic_obs,
        stratigraphic_height=p.stratigraphic_height,
        mesures=matched_mesures,
    )


def select_samples(
    pmag_list: List[Pmag],
    pattern: str = "*",
    step_min: float = 0,
    step_max: float = 9999,
    demag1: str = "*",
    demag2: str = "*",
    verbose: bool = True,
) -> List[SelectedSample]:
    """Equivalent de la subroutine `selmes`.

    - pattern : numero d'echantillon. '*' = tous. Sinon comparaison sur les
      `nlen` premiers caracteres (nlen = longueur de pattern) ; si pattern
      contient un '?' avant son dernier caractere, un echantillon est aussi
      retenu si son dernier caractere (position nlen) correspond a celui de
      pattern (comportement d'origine, y compris son cote approximatif : la
      position du '?' n'est pas utilisee, seul le dernier caractere compte).
    - step_min / step_max : bornes (incluses) sur le numero d'etape.
    - demag1 / demag2 : code de demagnetisation (ex 'N','0' pour NRM ; 'T'/'A'
      pour thermique/AF). '*' = tous.
    """
    pattern = pattern.upper().strip()[:12] or "*"
    nlen = len(pattern)
    has_question = "?" in pattern[: nlen - 1]

    selected: List[SelectedSample] = []
    total_mes = 0

    for p in pmag_list:
        pid = (p.id or "").ljust(12).upper()

        sample_match = (
            pattern[0] == "*"
            or pid[:nlen] == pattern[:nlen]
            or (has_question and pid[nlen - 1] == pattern[nlen - 1])
        )
        if not sample_match:
            continue

        matched_mesures = [
            m
            for m in p.mesures
            if step_min <= m.etape <= step_max
            and (demag1 == "*" or m.cod1 == demag1)
            and (demag2 == "*" or m.cod2 == demag2)
        ]
        if not matched_mesures:
            # Porteur de metadonnees de SITE sans specimen reel (id vide,
            # 0 mesure - "specimen: n.d" dans le .prmag, voir
            # magic_export.build_sites_rows) : jamais un candidat de
            # FIT (rien a demagnetiser), mais doit rester selectionnable
            # en mode "tout" (pattern='*', ET aucun filtre demag actif)
            # pour que son lat/lon/bed_dip/formation/age atteignent
            # l'export MagIC - sans ce cas particulier, `if not
            # matched_mesures: continue` l'excluait TOUJOURS de
            # self.selection, donc de tout export, meme apres avoir ete
            # correctement complete dans le .prmag (demande explicite
            # utilisateur : "the export of results without data lack the
            # information from the site that was added to the .prmag
            # files"). Le filtre demag1/demag2 (et pas seulement pattern)
            # doit AUSSI etre '*' - un porteur sans mesure ne "matche"
            # jamais un code demag precis, il ne doit donc pas apparaitre
            # des qu'une selection FILTRE reellement sur les mesures (ex.
            # demag='I') - demande explicite utilisateur ("when there are
            # sites without samples, do not select data", suite a une
            # selection filtree sur le code 'I' qui remontait quand meme
            # ces porteurs de site).
            if not (p.id or "").strip() and pattern == "*" and demag1 == "*" and demag2 == "*":
                selected.append(_build_selected_sample(p, []))
            continue

        selected.append(_build_selected_sample(p, matched_mesures))
        total_mes += len(matched_mesures)

    if verbose:
        print(f" nb samples: {len(selected):5d}  nb measurements: {total_mes:5d}")

    return selected


def select_samples_by_site(
    pmag_list: List[Pmag],
    site: str = "*",
    step_min: float = 0,
    step_max: float = 9999,
    demag1: str = "*",
    demag2: str = "*",
    verbose: bool = True,
) -> List[SelectedSample]:
    """Pas dans le Fortran (aucun equivalent "selection par site" n'existe
    dans le menu d'origine) : selectionne tous les echantillons dont
    `magic_site` (decode depuis la ligne roche, voir testlect.decode_roche)
    correspond exactement a `site` (insensible a la casse, espaces ignores).
    '*' = tous les echantillons ayant un site MagIC decode (non vide).
    Meme logique de filtrage etape/demag que `select_samples`."""
    site = (site or "*").strip()
    want_all = site in ("", "*")
    site_upper = site.upper()

    selected: List[SelectedSample] = []
    total_mes = 0

    for p in pmag_list:
        psite = (p.magic_site or "").strip()
        if not psite:
            continue
        if not want_all and psite.upper() != site_upper:
            continue

        matched_mesures = [
            m
            for m in p.mesures
            if step_min <= m.etape <= step_max
            and (demag1 == "*" or m.cod1 == demag1)
            and (demag2 == "*" or m.cod2 == demag2)
        ]
        if not matched_mesures:
            # Meme cas particulier que select_samples (voir sa docstring,
            # y compris la restriction demag1/demag2=='*' - "when there
            # are sites without samples, do not select data") : un
            # porteur de metadonnees de site (id vide, 0 mesure) reste
            # selectionnable si ce site precis (ou tous, want_all) est
            # demande, mais seulement hors filtre demag.
            if not (p.id or "").strip() and demag1 == "*" and demag2 == "*":
                selected.append(_build_selected_sample(p, []))
            continue

        selected.append(_build_selected_sample(p, matched_mesures))
        total_mes += len(matched_mesures)

    if verbose:
        print(f" nb samples: {len(selected):5d}  nb measurements: {total_mes:5d}")

    return selected


def _read_int(prompt: str, default: int) -> int:
    """Equivalent des boucles Fortran 7777/7778 : redemande tant que la
    valeur entree n'est pas un entier valide ; vide -> default."""
    while True:
        chaine = input(prompt).strip()
        if not chaine:
            return default
        try:
            return int(chaine)
        except ValueError:
            continue


def select_samples_interactive(pmag_list: List[Pmag], entete: str = "") -> List[SelectedSample]:
    """Equivalent interactif de `selmes` : pose les memes questions que le
    Fortran (numero d'echantillon, step min, step max, code demag NRM/Th/AF)
    puis appelle select_samples()."""
    raw = input(f" Sample number:{entete}").strip()
    if entete:
        combined = entete + raw
    else:
        combined = raw
    pattern = combined.strip()[:12].upper() or "*"

    step_min = _read_int(" Step min value:", 0)
    step_max = _read_int(" Step max value:", 9999)

    chaine = input(" Code NRM Th ou AF (example N0):").strip()
    if not chaine:
        demag1, demag2 = "*", "*"
    elif len(chaine) == 1:
        demag1, demag2 = chaine[0], "*"
    else:
        demag1, demag2 = chaine[0], chaine[1]

    return select_samples(
        pmag_list,
        pattern=pattern,
        step_min=step_min,
        step_max=step_max,
        demag1=demag1,
        demag2=demag2,
    )


def init_selection() -> List[SelectedSample]:
    """Equivalent de `initmes` : reinitialise la selection (liste vide)."""
    print("\n-- Init the list - no sample selected- ")
    print("--  Please select new samples  -- ")
    return []


# ---------------------------------------------------------------------------
# effmes : retrait d'echantillons/mesures d'une selection existante
# ---------------------------------------------------------------------------

def delete_measurements(
    selected: List[SelectedSample],
    pattern: str = "*",
    step_min: float = 0,
    step_max: float = 9000,
    demag1: str = "*",
    demag2: str = "*",
    occurrence: str = "*",
    verbose: bool = True,
) -> List[SelectedSample]:
    """Equivalent de la subroutine `effmes` : retire de `selected` les mesures
    qui correspondent aux criteres (et l'echantillon entier s'il ne lui reste
    plus aucune mesure). Contrairement a selmes/select_samples, `pattern` ne
    supporte pas le joker '?' (comportement d'origine de effmes).

    `occurrence` : PAS dans le Fortran (impossible d'y distinguer un doublon
    - meme etape/code - dont l'un est parfois une mesure erronee) - '*' =
    supprime TOUTES les mesures correspondantes (comportement d'origine),
    '1'/'2'/... = ne supprime QUE la 1ere/2eme/... mesure rencontree parmi
    celles qui correspondent, dans l'ordre du fichier, pour chaque
    echantillon - permet de retirer un doublon precis sans toucher a
    l'autre. Demande explicite de l'utilisateur."""
    pattern = (pattern or "*").upper().strip()[:12] or "*"
    nlen = len(pattern)
    occurrence = (occurrence or "*").strip()

    remaining: List[SelectedSample] = []
    for ech in selected:
        pid = (ech.id or "").ljust(12).upper()
        if not (pattern[0] == "*" or pid[:nlen] == pattern[:nlen]):
            remaining.append(ech)
            continue

        def matches(m: Measurement) -> bool:
            return (
                step_min <= m.etape <= step_max
                and (demag1 == "*" or m.cod1 == demag1)
                and (demag2 == "*" or m.cod2 == demag2)
            )

        if occurrence == "*":
            kept = [m for m in ech.mesures if not matches(m)]
        else:
            try:
                target_n = int(occurrence)
            except ValueError:
                target_n = 1
            kept = []
            count = 0
            for m in ech.mesures:
                if matches(m):
                    count += 1
                    if count == target_n:
                        continue
                kept.append(m)

        if kept:
            ech.mesures = kept
            remaining.append(ech)

    if verbose:
        nb_mesures = sum(len(e.mesures) for e in remaining)
        print(f" nb samples: {len(remaining):5d}  nb measurements: {nb_mesures:5d}")

    return remaining


def delete_measurements_interactive(selected: List[SelectedSample]) -> List[SelectedSample]:
    """Equivalent interactif de `effmes`."""
    pattern = input(" Sample to delete:").strip().upper()[:12] or "*"
    step_min = _read_int(" min step value:", 0)
    step_max = _read_int(" max step value:", 9000)

    chaine = input(" code NRM Th or AF (example N0):").strip()
    if not chaine:
        demag1, demag2 = "*", "*"
    elif len(chaine) == 1:
        demag1, demag2 = chaine[0], "*"
    else:
        demag1, demag2 = chaine[0], chaine[1]

    return delete_measurements(
        selected, pattern=pattern, step_min=step_min, step_max=step_max,
        demag1=demag1, demag2=demag2,
    )


# ---------------------------------------------------------------------------
# lismes : listage des mesures selectionnees
# ---------------------------------------------------------------------------

# iorient : 1=coordonnees echantillon, 2=in-situ, 3=apres correction de
# pendage, 4=apres correction d'orientation de forage (carottes de mine)

# idiom conserve bilingue (F/E, port fidele du Fortran d'origine), mais le
# defaut est passe de "F" a "E" - demande explicite utilisateur ("check
# and change the french words by english in the text output") : app.py
# n'a JAMAIS passe `idiom` explicitement a list_measurements/list_xyz/
# list_measurements_vrm/list_measurements_depth/diff_measurements,
# reposant donc silencieusement sur le defaut - qui produisait du texte
# FRANCAIS dans l'appli reelle (plusieurs elements de menu actifs) tant
# qu'il valait "F".
_HEADERS_FR = {
    1: "donnees en coordonnes echantillon",
    2: "donnees en coordonnes in situ",
    3: "donnees apres correction de pendage",
    4: "donnees apres correction de forage",
}
_HEADERS_EN = {
    1: "data in sample coordinate",
    2: "data in sample in situ",
    3: "data after tilt correction",
    4: "data after mining drill core orientation",
}


def _write_report_header(out: TextIO, idiom: str, orientation: int, columns_fr: str, columns_en: str) -> None:
    # HEADER_MARK ("\x01") entoure la ligne de titres de colonnes pour que
    # STARpaleomagApp._afficher (app.py) l'affiche en gras - demande explicite
    # utilisateur ("throughout the software, is it possible to write the
    # header in bold"). Meme caractere que app.HEADER_MARK, redefini
    # localement (pas d'import depuis app.py - creerait une dependance
    # circulaire, ce module etant importe par app.py, jamais l'inverse).
    if idiom == "F":
        out.write("\n ---- liste mesures selectionnees -----\n")
        out.write(f" {_HEADERS_FR[orientation]}\n")
        out.write(f"{_HEADER_MARK}{columns_fr}{_HEADER_MARK}\n")
    else:
        out.write("\n ---- list selected measurements -----\n")
        out.write(f" {_HEADERS_EN[orientation]}\n")
        out.write(f"{_HEADER_MARK}{columns_en}{_HEADER_MARK}\n")


def list_measurements(
    selected: List[SelectedSample],
    orientation: int = 1,
    idiom: str = "E",
    out: TextIO = sys.stdout,
) -> None:
    """Equivalent de la subroutine `lismes` : imprime les mesures des
    echantillons selectionnes, avec correction d'orientation optionnelle.

    Colonnes nettoyees (demande explicite utilisateur, "can we clean the
    list data") par rapport a l'ancien format brut Fortran :
    - "Etape"/"Step" porte maintenant l'unite du pas (degC thermique, mT
      AF - voir _fmt_step) plutot qu'un code compact "210D+" a decoder.
    - UNE SEULE colonne d'intensite normalisee (A/m OU Am2/kg, jamais les
      deux avec un "--" pour celle qui ne s'applique pas) - voir
      _selection_norme_labels/normalized_intensity.
    - Susceptibilite relabelee "K", avec son unite suivant elle aussi la
      normalisation reellement utilisee (SI pour un volume, m3/kg pour une
      masse) plutot que le double libelle "(s: SI ou m3/kg)" fige.
    - "Mag" (nom de colonne historique, sans rapport avec la magnetisation -
      c'etait deja le code INSTRUMENT) renomme "Inst"."""
    mag_label, k_label = _selection_norme_labels(selected)
    _write_report_header(
        out, idiom, orientation,
        f"       Numero        Etape          Mtot Am2      {mag_label:<15s}Dec    Inc   q  Inst   {k_label}",
        f"       Specimen      Step          Mtot Am2      {mag_label:<15s}Dec    Inc   q  Inst   {k_label}",
    )

    ij = 1
    for ech in selected:
        for m in ech.mesures:
            xx, yy, zz = apply_orientation(m.x, m.y, m.z, ech, orientation)
            mtot, dec, inc = polere(xx, yy, zz)
            intensity, _unit = normalized_intensity(ech, mtot)
            s_conv = _s_value(ech, m)

            out.write(
                f"{ij:5d}: {ech.id:<12s}  {_fmt_step(m)}  "
                f"{mtot:10.3E}  {intensity:10.3E}  "
                f"{dec:6.1f} {inc:6.1f}{m.q:4d}  {m.ins:<4s}  {s_conv:10.3E}  "
                f"{m.heuremes or '':<19s}\n"
            )
            ij += 1

    out.write("\n")


def list_xyz(
    selected: List[SelectedSample],
    orientation: int = 1,
    idiom: str = "E",
    out: TextIO = sys.stdout,
) -> None:
    """Equivalent de `listeXYZ` : comme list_measurements, mais affiche les
    composantes X,Y,Z du vecteur oriente au lieu de Dec/Inc (pas de polere).

    Memes nettoyages que list_measurements (demande explicite utilisateur,
    "the same treatment applied there too") : etape avec son unite,
    susceptibilite relabelee "K" avec son unite suivant la normalisation
    reellement utilisee, "Mag" -> "Inst". X/Y/Z n'ont pas la duplication
    A/m + Am2/kg de list_measurements (une seule colonne normalisee par
    composante) - leur unite commune est indiquee une fois, dans la ligne
    de description au-dessus du tableau plutot que repetee 3 fois dans
    l'en-tete de colonnes."""
    mag_unit, k_label = _selection_norme_labels(selected)
    _write_report_header(
        out, idiom, orientation,
        f"       Numero      Etape          Mtot Am2          X          Y           Z ({mag_unit})   q  Inst   {k_label}",
        f"       Specimen      Step          Mtot Am2          X          Y           Z ({mag_unit})   q  Inst   {k_label}",
    )

    ij = 1
    for ech in selected:
        # Pas de facteur 1e3/1e6 sans volume/masse (0/None) - affiche les
        # composantes brutes (Am2) inchangees plutot qu'un faux "0.0" -
        # demande explicite utilisateur ("when the volume or the mass of
        # the sample is not given, best to show the data in total moment
        # as done in the list"), meme principe que normalized_intensity.
        factor = (1.0e3 if ech.norme == "m" else 1.0e6) if ech.vol else 1.0
        vol = ech.vol or 1.0
        for m in ech.mesures:
            xx, yy, zz = apply_orientation(m.x, m.y, m.z, ech, orientation)
            rxx = xx * factor / vol
            ryy = yy * factor / vol
            rzz = zz * factor / vol
            s_conv = _s_value(ech, m)

            out.write(
                f"{ij:5d}: {ech.id:<12s}  {_fmt_step(m)}  "
                f"{xx:10.3E}  {rxx:10.3E}  {ryy:10.3E}  {rzz:10.3E}{m.q:4d}  {m.ins:<4s}  {s_conv:10.3E}\n"
            )
            ij += 1


def list_measurements_vrm(
    selected: List[SelectedSample],
    orientation: int = 1,
    idiom: str = "E",
    out: TextIO = sys.stdout,
) -> None:
    """Equivalent de `lismesVRM` : comme list_measurements, avec une colonne
    d'heures ecoulees depuis la premiere mesure de l'echantillon (utile pour
    les etudes de viscosite / VRM). Memes nettoyages que list_measurements
    (demande explicite utilisateur, "the same treatment applied there
    too")."""
    mag_label, k_label = _selection_norme_labels(selected)
    _write_report_header(
        out, idiom, orientation,
        f"       Numero      Etape          heures     Mtot Am2      {mag_label:<15s}Dec    Inc   q  Inst   {k_label}",
        f"       Specimen      Step          hours      Mtot Am2      {mag_label:<15s}Dec    Inc   q  Inst   {k_label}",
    )

    ij = 1
    for ech in selected:
        if not ech.mesures:
            continue
        t0 = _parse_heuremes(ech.mesures[0].heuremes)
        for m in ech.mesures:
            tj = _parse_heuremes(m.heuremes)
            elapsed_h = (tj - t0).total_seconds() / 3600.0 if (t0 and tj) else 0.0

            xx, yy, zz = apply_orientation(m.x, m.y, m.z, ech, orientation)
            mtot, dec, inc = polere(xx, yy, zz)
            intensity, _unit = normalized_intensity(ech, mtot)
            s_conv = _s_value(ech, m)

            out.write(
                f"{ij:5d}: {ech.id:<12s}  {_fmt_step(m)}  {elapsed_h:9.4f}  "
                f"{mtot:10.3E}  {intensity:10.3E}  "
                f"{dec:6.1f} {inc:6.1f}{m.q:4d}  {m.ins:<4s}  {s_conv:10.3E}  "
                f"{m.heuremes or '':<19s}\n"
            )
            ij += 1


def list_measurements_depth(
    selected: List[SelectedSample],
    depths: Dict[str, float],
    expected_dec: float = 0.0,
    expected_inc: float = 0.0,
    orientation: int = 1,
    idiom: str = "E",
    out: TextIO = sys.stdout,
) -> None:
    """Equivalent de `lismesdepth` : comme list_measurements, avec la
    profondeur de l'echantillon (depuis `depths`, {sample_id: depth} ; -999.0
    si absent) et l'angle par rapport a la direction attendue (D,I) au site,
    calcule via `angle()` (equivalent de la subroutine `angle` de calcul.f).
    Memes nettoyages que list_measurements (demande explicite utilisateur,
    "the same treatment applied there too")."""
    mag_label, k_label = _selection_norme_labels(selected)
    _write_report_header(
        out, idiom, orientation,
        f"  Prof     Numero      Etape          Mtot Am2      {mag_label:<15s}Dec    Inc   q  Inst   {k_label}  Angle",
        f"  depth      Specimen      Step          Mtot Am2      {mag_label:<15s}Dec    Inc   q  Inst   {k_label}  Angle",
    )

    for ech in selected:
        depth = depths.get(ech.id, -999.0)
        for m in ech.mesures:
            xx, yy, zz = apply_orientation(m.x, m.y, m.z, ech, orientation)
            mtot, dec, inc = polere(xx, yy, zz)
            delta, _ang = angle(expected_dec - dec, 90.0 - inc, 90.0 - expected_inc)
            if dec > 270.0:
                dec -= 360.0

            intensity, _unit = normalized_intensity(ech, mtot)
            s_conv = _s_value(ech, m)

            out.write(
                f"{depth:8.2f}   {ech.id:<12s}  {_fmt_step(m)}  "
                f"{mtot:10.3E}  {intensity:10.3E}  "
                f"{dec:6.1f} {inc:6.1f}{m.q:4d}  {m.ins:<4s}  {s_conv:10.3E}  {delta:6.1f}\n"
            )


def diff_measurements(
    selected: List[SelectedSample],
    orientation: int = 1,
    idiom: str = "E",
    out: TextIO = sys.stdout,
) -> None:
    """Equivalent de `diffmes` : pour chaque echantillon, imprime la
    difference vectorielle entre mesures consecutives (etape j et j+1),
    apres correction d'orientation."""
    if idiom == "F":
        out.write("\n ---- liste mesures selectionnees -----\n")
    else:
        out.write("\n ---- list selected measurements -----\n")
    out.write(f" {(_HEADERS_FR if idiom == 'F' else _HEADERS_EN)[orientation]}\n")

    ij = 1
    for ech in selected:
        for j in range(len(ech.mesures) - 1):
            a, b = ech.mesures[j], ech.mesures[j + 1]
            dx, dy, dz = a.x - b.x, a.y - b.y, a.z - b.z
            xx, yy, zz = apply_orientation(dx, dy, dz, ech, orientation)
            mtot, dec, inc = polere(xx, yy, zz)
            # sans volume/masse (0/None) : moment total brut (Am2), pas un
            # faux "0.0" - demande explicite utilisateur ("when the volume
            # or the mass of the sample is not given, best to show the
            # data in total moment as done in the list").
            mtot, _unit = normalized_intensity(ech, mtot)

            out.write(
                f"{ij:4d}: {ech.id:<12s}  {_fmt_step(a)}  "
                f"{mtot:10.3E} {dec:6.1f} {inc:6.1f}{a.q:4d}  {a.ins:<2s}{a.s:7.1f}\n"
            )
            ij += 1


def subtract_measurement(selected: List[SelectedSample], row_number: int) -> List[SelectedSample]:
    """Equivalent de `soustra` : soustrait, dans l'unique echantillon
    selectionne, le vecteur (x,y,z) de la mesure `row_number` (1-indexe, tel
    qu'affiche par list_measurements) a toutes les mesures de cet
    echantillon, puis retire cette mesure de la selection. Modifie et
    retourne `selected`."""
    if len(selected) != 1:
        raise ValueError("soustra requires a selection limited to a single sample")

    ech = selected[0]
    idx = row_number - 1
    base = ech.mesures[idx]
    xsous, ysous, zsous = base.x, base.y, base.z

    for m in ech.mesures:
        m.x -= xsous
        m.y -= ysous
        m.z -= zsous

    del ech.mesures[idx]
    return selected


# (libelle, largeur, alignement) - demande explicite utilisateur
# ("est-ce possible d'aligner les colonnes avec la ligne de
# presentation") : la meme table de largeurs sert a construire l'en-tete
# ET chaque ligne de donnees (comme calcul._row/_cols_header), pour
# garantir l'alignement plutot que de dupliquer des f-strings a largeurs
# fixes independantes de l'en-tete (source du desalignement precedent,
# ex. "lithology" plus long que sa largeur decalait tout ce qui suivait).
_INFO_FIELDS = [
    ("specimen", 14, "<"), ("nb_meas", 9, ">"), ("Step[1]", 10, ">"), ("Step[n]", 10, ">"),
    ("Decli_IGRF", 12, ">"), ("Decli_local", 13, ">"), ("lithology", 18, "<"),
    ("core_dip", 10, ">"), ("azimuth", 10, ">"), ("bed_dip", 10, ">"), ("bed_strike", 12, ">"),
    ("vol/mass", 11, ">"), ("Lat", 11, ">"), ("Long", 11, ">"),
]


def _info_row(values: List[str]) -> str:
    # separateur "  " explicite entre colonnes (pas juste des largeurs
    # fixes accolees) : une largeur pile a la longueur de son libelle
    # (ex. "Decli_local", 11 caracteres) ne laissait sinon aucun espace
    # avant la colonne suivante des que sa PROPRE valeur faisait au moins
    # aussi large que le libelle.
    return "  ".join(
        (v.ljust(w) if align == "<" else v.rjust(w))
        for (_label, w, align), v in zip(_INFO_FIELDS, values)
    )


def sample_info(selected: List[SelectedSample]) -> str:
    """Equivalent de `infoech` (dataselect2.f:937-1064) : tableau
    recapitulatif des echantillons selectionnes. Colonnes/ordre choisis
    par l'utilisateur (demande explicite) : specimen, nb_meas, Step[1],
    Step[n], Decli_IGRF, Decli_local (voir orient_sample.py), lithology
    (magic_li - remplace le dump brut de la ligne "roche" d'une version
    precedente), core_dip (cin), azimuth (caz), bed_dip (dip), bed_strike
    (str_), vol/mass (vol, suffixe par la lettre `norme` - 'v' ou 'm' -
    pour lever l'ambiguite), Lat, Long.

    "n.d" (au lieu d'une valeur calculee) - demande explicite
    utilisateur ("when there is no solar data, n.d for solar data") :
    - Decli_IGRF/Decli_local/Lat/Long : specimen sans annee (`year==0`,
      meme condition que le Fortran - aucune position/date exploitable).
    - Decli_local SEULE (Decli_IGRF reste calculee - l'IGRF ne depend que
      de la position/date, pas d'une visee solaire) : aucune visee
      solaire faite (azsun/hour/minute tous nuls, meme convention que
      compute_local_declination pour la valeur Fortran 0.0 qu'elle
      retourne dans ce cas - reinterpretee ici comme "non determine"
      plutot que comme une declinaison locale reellement nulle)."""
    header_row = _info_row([label for label, _w, _a in _INFO_FIELDS])
    lines = [f"{_HEADER_MARK}{header_row}{_HEADER_MARK}"]
    for i, ech in enumerate(selected, start=1):
        if ech.mesures:
            first, last = ech.mesures[0], ech.mesures[-1]
            step1 = f"{first.etape}{first.cod1}{first.cod2}"
            stepn = f"{last.etape}{last.cod1}{last.cod2}"
        else:
            step1 = stepn = ""

        no_solar = ech.azsun == 0.0 and ech.hour == 0 and ech.minute == 0
        if ech.year:
            # meme convention que infoech (dataselect2.f:1029-1030) :
            # outil d'orientation commencant par 'A' ("aiguille verticale
            # sur platine tournante") = 1, tout le reste ("equerre
            # pivotante sur platine fixe") = 2.
            ioutil = 1 if ech.outilorient[:1] == "A" else 2
            decli_igrf, declin = compute_local_declination(
                ech.lat, ech.rlong, int(ech.year), int(ech.month), ech.day,
                ech.hour, ech.minute, ioutil, ech.azmag, ech.azsun,
            )
            decli_igrf_txt = f"{decli_igrf:.1f}"
            decli_local_txt = "n.d" if no_solar else f"{declin:.1f}"
            lat_txt = f"{ech.lat:.5f}"
            long_txt = f"{ech.rlong:.5f}"
        else:
            decli_igrf_txt = decli_local_txt = lat_txt = long_txt = "n.d"

        vol_mass = f"{ech.vol:.3f}{ech.norme or 'v'}"
        lines.append(_info_row([
            f"{i}: {ech.id}", str(len(ech.mesures)), step1, stepn,
            decli_igrf_txt, decli_local_txt, ech.magic_li,
            f"{ech.cin:.1f}", f"{ech.caz:.1f}", f"{ech.dip:.1f}", f"{ech.str_:.1f}",
            vol_mass, lat_txt, long_txt,
        ]))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from testlect import read_ren_file

    fichier = "testfile.ren.txt"
    donnees = read_ren_file(fichier)
    print(f"{len(donnees)} echantillon(s) lu(s)\n")

    selection = select_samples_interactive(donnees)
    list_measurements(selection, orientation=1)
