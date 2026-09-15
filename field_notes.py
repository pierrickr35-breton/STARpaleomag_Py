"""
Menu "Field notes to Pmag files" - port de `subroutine orientation`
(orient_paleomag.f:5-441, reference/Starmac_AWE_v22/orient_paleomag.f) -
demande explicite utilisateur ("ajouter le menu field notes to Pmag files.
From the field orientations, it prepares the prmag files (check the
Fortran file). It also prepares the .ged file for the AGICO instruments").

Format du fichier "simple" d'entree - REVU (2e demande explicite
utilisateur) pour etre plus rapide a remplir sur le terrain, voir
l'exemple reel fourni (simple.txt) :

- "!..." en debut de ligne = commentaire, ignoree.
- Une ligne commencant par "#" est soit une ligne d'ENTETE/legende (si le
  premier mot APRES le "#" est "site" ou "core" - ignoree), soit une ligne
  de site dont le "#" colle directement au code (ex. "#EX02...", vu tel
  quel dans l'exemple reel) - le "#" est alors simplement retire avant de
  lire la ligne normalement (tolerant, pas une erreur).
- Ligne de site (11 champs, PLUS de guillemets, PLUS de bloc "&year:" -
  chaque site est desormais AUTONOME) : "code_site lat_deg lat_min
  lon_deg lon_min annee utdif mois jour pendage_couche
  direction_ou_strike". Derniere colonne : direction de pendage OU strike
  selon la convention choisie (voir `strati` ci-dessous - PAS fige,
  demande explicite utilisateur : "we need the two options, some use a
  compass where they record bed dip direction"). Plus de champ altitude
  ni de champs geologie (Fm/Age/GC/SMT/Li/Loc/Obs) sur cette ligne - voir
  fichier complement ci-dessous.
- Une ligne par carotte (INCHANGE) : "code_carotte plongee_carotte
  az_boussole az_soleil heure minute" - az_soleil/heure/minute tous a 0
  si aucune visee solaire n'a ete faite (secours sur l'IGRF). Une ligne
  speciale dont le code commence par 'C' (ex. "CC 0 0 0 0 0") termine le
  site. Le commentaire "!..." en fin de ligne (ex. "!TR PMO") est
  toujours conserve (voir plus bas).
- `outil`/`ellipsoide`/`strati` (ASC vs autre outil, WGS84 vs autre datum,
  Strikedip vs Dipdip) ne sont PLUS dans le fichier du tout (ne variaient
  pas d'un site a l'autre dans les fichiers reels vus jusqu'ici) -
  demandes UNE FOIS a l'appelant (voir parametres `tool`/`ellipsoid`/
  `strati` de `parse_orientation_file`, prompts dans app.py), appliques
  uniformement a tous les sites du fichier plutot que repetes ligne apres
  ligne. `strati` : "Dipdip" (la colonne est une direction de pendage,
  boussole donnant l'azimut de plongement) -> strike=valeur-90 (regle
  main droite) ; "Strikedip" (la colonne EST deja le strike) -> aucune
  conversion - voir _parse_site_line_simple.

Fichier "complement" (OPTIONNEL, meme principe que le complement deja
utilise par `convert_legacy_auto` pour un vieux fichier .ren) - demande
explicite utilisateur : "as we now have the possibility to complete the
different fields later, we can have two files (a simple file written
following the field trip) and a complement file that can be written
later". Permet de fournir/completer APRES le terrain :

- la geologie d'un site : une ligne "code_site "Fm" "Age" "GC" "SMT"
  "Litho" "Location" ["Obs"]" (memes guillemets que l'ancien format .ori -
  du texte libre, une virgule dans "Location" ne doit pas couper le
  champ).
- des carottes OUBLIEES/ajoutees apres coup pour ce site : une ou
  plusieurs lignes au format carotte habituel juste apres la ligne de
  geologie du site concerne (avant le prochain code_site ou la fin du
  fichier) - voir parse_complement_file/apply_complement.

Le code_site en tete de chaque ligne de geologie est OBLIGATOIRE dans ce
port (necessaire pour savoir a quel site du fichier simple rattacher
l'information) - l'extrait fourni par l'utilisateur en montre le contenu
(texte geologie + une carotte supplementaire "03" manquante du fichier
simple pour EX02) sans code_site explicite devant, une fois complete par
l'utilisateur avec le code_site en tete de ligne.

Sortie :
- Un fichier .prmag (FORMAT MODERNE STARpaleomag_Py - `_sample_header_block`,
  directement rechargeable via read_prmag_file, PAS le texte brut "Id:/L:/
  Site:" du Fortran d'origine, obsolete/non relisible par cette
  application) - mesures VIDES (aucun instrument n'a encore mesure quoi
  que ce soit a ce stade, uniquement l'orientation de terrain).
- UN SEUL fichier .ged combinant TOUS les sites du fichier d'entree, pour
  les instruments AGICO (format SUFAR/Anisoft/Remasoft) - correspond au
  comportement REEL du Fortran (`unfichierged='Y'` code en dur - un mode
  par-site existe dans le source, complet et fonctionnel, mais
  inatteignable, gate derriere une saisie interactive entierement
  commentee). Revenu sur un choix initial de ce port (fichiers separes par
  site) a la demande explicite de l'utilisateur : "It is much more simple
  to handle a single .ged file but all samples need to share the same
  orientation scheme as in prmag. My own software (STARpaleomag_Py) handle the
  pmag or AMS data" - un seul schema d'orientation (12_0_3_90) s'applique
  deja uniformement a tous les specimens quel que soit leur site, rendant
  la separation par site inutile pour l'usage reel. Ligne d'entete
  (format 666) ecrite UNE SEULE FOIS a la fin, avec "????" comme code site
  (pas de site unique pour un fichier combine) et lat/lon du DERNIER site
  traite - meme convention que le Fortran, qui ne reinitialise jamais
  lat/rlong pour ce recapitulatif final (`write(21,666) iligne+1, "????",
  lat,rlong` au label 999, `iligne` n'etant JAMAIS remis a 1 en mode
  combine - seule la branche par-site desactivee le fait).
- Un rapport diagnostic (equivalent du fichier ".res" - azgeo/declinaison
  locale/IGRF par specimen, avec "err>3/5/10" si l'ecart depasse un seuil).

BUG Fortran CONFIRME (a signaler, pas silencieusement reproduit tel quel) :
`d_err` (ecart utilise pour "err>3/5/10") n'est JAMAIS assigne dans le
source (aucune ligne "d_err=...") - variable non initialisee. L'intention
est evidente (comparer la declinaison locale mesuree `declin` au modele
IGRF `decli_igrf`) : ce port calcule explicitement `d_err = declin -
decli_igrf`, le calcul manifestement voulu mais jamais ecrit dans le
Fortran.

AJOUT (pas dans le Fortran, qui ignore tout texte au-dela des 6 champs
numeriques attendus sur une ligne de carotte) : le commentaire "!..." en
fin de ligne de carotte (ex. "!TR PMO" - Testigo Roto/Posiblemente mal
orientado, per les abreviations documentees en tete du fichier reel) est
conserve dans le champ `comment` du specimen plutot que silencieusement
perdu - purement additif, ne change aucune valeur calculee."""

import math
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from testlect import Pmag
from convert_ren_to_r import _sample_header_block, _MEAS_HEADER
from orient_sample import igrf_declination, _sun_declination_and_eot, _spherical_azimuth
from calcul import ani_path_for, create_empty_pmagani_if_missing

_LETTERS = ("A", "B", "C")

# Correction fixe de datum (SAD69/PSAD56 vs WGS84) pour la zone Terre de
# Feu/Cap Horn - reprise TELLE QUELLE du Fortran (`if(wgs=='n') then if(lat
# dans cette boite) then lat=lat-0.0035 ; rlong=rlong-0.0020`), applicable
# uniquement quand l'ellipsoide du site n'est PAS "WGS84".
_DATUM_SHIFT_BOX = (-60.0, 0.0, -85.0, -60.0)  # lat_min, lat_max, lon_min, lon_max
_DATUM_SHIFT = (-0.0035, -0.0020)

# Volume de carotte fixe (cm3) - constante codee en dur dans le Fortran
# (`volume=10.8`), jamais demandee a l'utilisateur.
_DEFAULT_VOLUME = 10.8


@dataclass
class FieldSpecimenRow:
    car_sample: str
    dip_core: float
    azmag: float
    azsun: float
    hour: float
    minute: float
    comment: str = ""


@dataclass
class FieldSite:
    car_site: str
    lat: float
    rlong: float
    altitude: float
    month: int
    day: int
    plane_dip: float
    strike: float
    year: int
    utdif: float
    wgs: str
    outil: str
    # Geologie (Fm/Age/GC/SMT/Li/Loc/Obs) - VIDE par defaut dans le
    # fichier "simple" (voir docstring module), completee eventuellement
    # via apply_complement/parse_complement_file.
    fm: str = ""
    age: str = ""
    gc: str = ""
    smt: str = ""
    litho: str = ""
    location: str = ""
    obs: str = ""
    rows: List[FieldSpecimenRow] = field(default_factory=list)
    terminated: bool = True


def _tokenize(line: str) -> Tuple[List[str], str]:
    """Tokenise `line` en respectant les groupes entre guillemets (ex.
    "Hermite Island, West" - une virgule DANS les guillemets ne doit pas
    couper le champ) ; un "!" NON guillemete arrete la tokenisation (le
    reste de la ligne est un commentaire, retourne separement - equivalent
    du fait qu'une lecture Fortran list-directed s'arrete des que le
    nombre de variables demandees est atteint, ignorant tout ce qui suit)."""
    tokens: List[str] = []
    i, n = 0, len(line)
    comment = ""
    while i < n:
        while i < n and line[i].isspace():
            i += 1
        if i >= n:
            break
        if line[i] == '"':
            j = line.find('"', i + 1)
            if j == -1:
                j = n
            tokens.append(line[i + 1:j])
            i = j + 1
        elif line[i] == "!":
            comment = line[i + 1:].strip()
            break
        else:
            j = i
            while j < n and not line[j].isspace() and line[j] != "!":
                j += 1
            tokens.append(line[i:j])
            i = j
    return tokens, comment


def _to_float(token: str) -> Optional[float]:
    try:
        return float(token)
    except (TypeError, ValueError):
        return None


_LABEL_WORDS = {"site", "core"}


def _strip_comment_hash(stripped: str) -> Tuple[str, bool]:
    """Retire un "#" de tete (ex. "#EX02..." rencontre tel quel dans le
    fichier reel de l'utilisateur, colle au code site plutot que separe
    par un espace) - retourne (texte sans le "#", True si un "#" a ete
    retire). Le second element sert a decider si la ligne est une legende
    (voir _is_label_tokens) ou une vraie ligne de donnees."""
    if stripped.startswith("#"):
        return stripped[1:].strip(), True
    return stripped, False


def _is_label_tokens(tokens: List[str]) -> bool:
    return bool(tokens) and tokens[0].strip().lower() in _LABEL_WORDS


def _parse_site_line_simple(
    tokens: List[str], tool: str, ellipsoid: str, strati: str,
) -> Optional[FieldSite]:
    """"code_site lat_deg lat_min lon_deg lon_min annee utdif mois jour
    pendage_couche direction_ou_strike" (11 champs - voir docstring
    module). Tokens excedentaires au-dela du 11e sont ignores (tolerant,
    ex. valeur trainante accidentelle dans un fichier tape a la main).

    `strati` : convention utilisee pour la derniere colonne - "Dipdip"
    (direction de pendage, boussole donnant directement l'azimut de
    plongement - "some use a compass where they record bed dip
    direction") -> strike = valeur-90 (regle main droite) ; "Strikedip"
    (la valeur EST deja le strike) -> aucune conversion. Redevenue un
    reglage explicite (demande utilisateur : "we need the two options"),
    comme dans l'ancien fichier .ori, mais desormais GLOBAL a tout
    l'import (voir parse_orientation_file) plutot que reprecise par bloc
    "&year:" - ne varie pas d'un site a l'autre dans les fichiers reels
    vus jusqu'ici."""
    if len(tokens) < 11:
        return None
    vals = [_to_float(t) for t in tokens[1:11]]
    if any(v is None for v in vals):
        return None
    car_site = tokens[0]
    rlat, rlat_min, rlon, rlon_min, year, utdif, rmonth, riday, plane_dip, direction_value = vals

    lat = math.copysign(abs(rlat) + rlat_min / 60.0, rlat)
    rlong = math.copysign(abs(rlon) + rlon_min / 60.0, rlon)
    if ellipsoid.upper() != "WGS84":
        lat_min, lat_max, lon_min, lon_max = _DATUM_SHIFT_BOX
        if lat_min < lat < lat_max and lon_min < rlong < lon_max:
            lat += _DATUM_SHIFT[0]
            rlong += _DATUM_SHIFT[1]

    strike = direction_value - 90.0 if strati.strip().lower() == "dipdip" else direction_value

    return FieldSite(
        car_site=car_site, lat=lat, rlong=rlong, altitude=0.0,
        month=int(rmonth), day=int(riday), plane_dip=plane_dip, strike=strike,
        year=int(year), utdif=utdif, wgs=ellipsoid, outil=tool,
    )


def parse_orientation_file(
    lines: List[str], tool: str = "ASC", ellipsoid: str = "WGS84", strati: str = "Strikedip",
) -> Tuple[List[FieldSite], List[str]]:
    """Lit le fichier "simple" (voir docstring module) - `tool`/`ellipsoid`/
    `strati` s'appliquent a TOUS les sites du fichier (ne varient plus par
    site dans ce format, voir _parse_site_line_simple pour `strati`).
    Retourne (sites, warnings) - un site incomplet (jamais termine par une
    ligne "C..." avant la fin du fichier) est quand meme retourne (son
    `.prmag` reste exploitable), avec `terminated=False` et un
    avertissement."""
    sites: List[FieldSite] = []
    warnings: List[str] = []
    current: Optional[FieldSite] = None

    for lineno, raw in enumerate(lines, start=1):
        stripped = raw.strip()
        if not stripped or stripped.startswith("!"):
            continue
        working, had_hash = _strip_comment_hash(stripped)
        tokens, comment = _tokenize(working)
        if not tokens:
            continue
        if had_hash and _is_label_tokens(tokens):
            continue  # ligne de legende ("#site ..."/"#core ...")

        if current is None:
            site = _parse_site_line_simple(tokens, tool, ellipsoid, strati)
            if site is None:
                warnings.append(f"line {lineno}: malformed site line, ignored: {stripped!r}")
                continue
            current = site
            continue

        if len(tokens) < 6:
            warnings.append(f"line {lineno}: malformed core line, ignored: {stripped!r}")
            continue
        car_sample = tokens[0]
        if car_sample.strip().upper().startswith("C"):
            sites.append(current)
            current = None
            continue
        vals = [_to_float(t) for t in tokens[1:6]]
        if any(v is None for v in vals):
            warnings.append(f"line {lineno}: malformed core line, ignored: {stripped!r}")
            continue
        dip_core, azmag, azsun, hour, minute = vals
        current.rows.append(FieldSpecimenRow(
            car_sample=car_sample, dip_core=dip_core, azmag=azmag,
            azsun=azsun, hour=hour, minute=minute, comment=comment,
        ))

    if current is not None:
        current.terminated = False
        warnings.append(
            f"site « {current.car_site} » was not terminated by a 'C...' line "
            f"before the end of the file - included anyway.")
        sites.append(current)

    return sites, warnings


def parse_complement_file(lines: List[str]) -> Tuple[Dict[str, dict], List[str]]:
    """Lit un fichier "complement" (OPTIONNEL - voir docstring module) :
    une ligne 'code_site "Fm" "Age" "GC" "SMT" "Litho" "Location"
    ["Obs"]' fournit/complete la geologie d'un site deja present dans le
    fichier simple ; les lignes suivantes AU FORMAT CAROTTE (jusqu'au
    prochain code_site ou la fin du fichier) sont des carottes
    supplementaires a ajouter a ce meme site (oubliees/ajoutees apres le
    terrain). Retourne (dict code_site -> {fm,age,gc,smt,litho,location,
    obs,extra_rows}, warnings)."""
    result: Dict[str, dict] = {}
    warnings: List[str] = []
    current_site: Optional[str] = None

    for lineno, raw in enumerate(lines, start=1):
        stripped = raw.strip()
        if not stripped or stripped.startswith("!"):
            continue
        working, had_hash = _strip_comment_hash(stripped)
        tokens, comment = _tokenize(working)
        if not tokens:
            continue
        if had_hash and _is_label_tokens(tokens):
            continue

        # ligne de geologie : code_site suivi d'au moins une chaine entre
        # guillemets (ce qui la distingue d'une ligne carotte, purement
        # numerique apres le premier champ).
        if len(tokens) >= 2 and _to_float(tokens[1]) is None:
            car_site = tokens[0]
            rest = tokens[1:]
            if len(rest) < 6:
                warnings.append(
                    f"line {lineno}: incomplete geology line for site "
                    f"« {car_site} » (need Fm/Age/GC/SMT/Li/Loc), ignored: {stripped!r}")
                continue
            fm, age, gc, smt, litho, location = rest[:6]
            obs = rest[6] if len(rest) >= 7 else ""
            result[car_site] = {
                "fm": fm, "age": age, "gc": gc, "smt": smt,
                "litho": litho, "location": location, "obs": obs,
                "extra_rows": [],
            }
            current_site = car_site
            continue

        # sinon : ligne carotte supplementaire, rattachee au dernier
        # code_site de geologie rencontre.
        if current_site is None:
            warnings.append(
                f"line {lineno}: core line found before any site geology "
                f"line, ignored (no site to attach it to): {stripped!r}")
            continue
        if len(tokens) < 6:
            warnings.append(f"line {lineno}: malformed core line, ignored: {stripped!r}")
            continue
        vals = [_to_float(t) for t in tokens[1:6]]
        if any(v is None for v in vals):
            warnings.append(f"line {lineno}: malformed core line, ignored: {stripped!r}")
            continue
        dip_core, azmag, azsun, hour, minute = vals
        result[current_site]["extra_rows"].append(FieldSpecimenRow(
            car_sample=tokens[0], dip_core=dip_core, azmag=azmag,
            azsun=azsun, hour=hour, minute=minute, comment=comment,
        ))

    return result, warnings


def apply_complement(sites: List[FieldSite], complement: Dict[str, dict]) -> List[str]:
    """Applique `complement` (parse_complement_file) sur `sites` EN PLACE :
    geologie copiee sur le site correspondant, carottes supplementaires
    ajoutees a la suite de ses carottes existantes. Retourne les
    avertissements (code_site du complement sans site correspondant dans
    le fichier simple)."""
    warnings: List[str] = []
    by_site = {s.car_site: s for s in sites}
    for car_site, info in complement.items():
        site = by_site.get(car_site)
        if site is None:
            warnings.append(
                f"complement file: site « {car_site} » not found in the "
                f"field notes file, ignored.")
            continue
        site.fm = info["fm"] or site.fm
        site.age = info["age"] or site.age
        site.gc = info["gc"] or site.gc
        site.smt = info["smt"] or site.smt
        site.litho = info["litho"] or site.litho
        site.location = info["location"] or site.location
        site.obs = info["obs"] or site.obs
        site.rows.extend(info["extra_rows"])
    return warnings


def _wrap360(value: float) -> float:
    if value >= 360.0:
        return value - 360.0
    if value < 0.0:
        return value + 360.0
    return value


def compute_specimen_geometry(
    site: FieldSite, row: FieldSpecimenRow, decli_igrf: float, utdif: Optional[float] = None,
):
    """Port des lignes 227-334 de `orientation` pour UN specimen : retourne
    (caz, plane_strike, declin, comm). `decli_igrf` est calcule UNE FOIS
    par site (voir _site_decli_igrf), pas recalcule par specimen, comme
    dans le Fortran (le meme decli_igrf sert a toutes les carottes d'un
    site, avec ou sans visee solaire). `utdif` : `site.utdif` par defaut
    (None) - parametre EXPOSE (pas dans le Fortran) pour permettre a
    check_utdif de re-essayer d'autres decalages horaires sans dupliquer
    ce calcul."""
    if utdif is None:
        utdif = site.utdif
    ioutil = 1 if site.outil.strip().upper() == "ASC" else 2
    hour = row.hour + utdif
    minute = row.minute

    if row.azsun == 0.0 and row.hour == 0.0 and row.minute == 0.0:
        azgeo = row.azmag + decli_igrf
        declin = 0.0
        itest = 0
    else:
        day = float(site.day) + (hour + minute / 60.0) / 24.0 - 0.5
        dele, av = _sun_declination_and_eot(site.year, site.month, day)
        p = site.rlong + ((hour - 12.0) + (minute / 60.0)) * 15.0 + av
        if p > 360.0:
            p -= 360.0
        a = 90.0 - dele
        b = 90.0 - site.lat
        _dell, az = _spherical_azimuth(p, a, b)
        if 0.0 <= p < 180.0:
            az = 360.0 - az
        az += 180.0
        if az > 360.0:
            az -= 360.0
        azgeo = (row.azsun + az) if ioutil == 1 else (360.0 - row.azsun + az)
        azgeo = _wrap360(azgeo)
        declin = azgeo - row.azmag
        if declin < -180.0:
            declin += 360.0
        itest = 1

    # BUG Fortran confirme (d_err jamais assigne dans le source) - calcul
    # manifestement voulu, complete ici (voir docstring module).
    comm = ""
    if itest == 1:
        d_err = declin - decli_igrf
        if abs(d_err) > 10.0:
            comm = "err>10"
        elif abs(d_err) > 5.0:
            comm = "err>5"
        elif abs(d_err) > 3.0:
            comm = "err>3"

    caz = _wrap360(azgeo + 90.0)
    plane_strike = _wrap360(site.strike + decli_igrf)
    return caz, plane_strike, declin, comm


def _site_decli_igrf(site: FieldSite) -> float:
    return igrf_declination(site.lat, site.rlong, site.year, site.month, float(site.day))


def check_utdif(
    site: FieldSite,
    decli_igrf: Optional[float] = None,
    candidates: Optional[List[float]] = None,
    improvement_threshold: float = 3.0,
) -> Optional[str]:
    """Diagnostic (PAS dans le Fortran) : verifie que `site.utdif` (le
    decalage horaire UTC/fuseau saisi dans le fichier de terrain) est bien
    celui qui MINIMISE l'ecart entre la declinaison locale mesuree
    (`declin`, deduite de la visee solaire) et le modele IGRF - demande
    explicite utilisateur ("one of the most common error is with the UTM
    dif. Can you have a check that the difference between the local
    declination and IGRF is effectively the least important with the UTM
    difference in the file").

    La position du soleil est TRES sensible a l'heure (~15 deg d'azimut
    par heure d'ecart) alors que l'IGRF n'en depend pas du tout : un
    utdif errone se traduit generalement par un ecart declin/IGRF bien
    plus grand que necessaire - une erreur de saisie frequente ("one of
    the most common error"), detectable en balayant les decalages entiers
    voisins (-12 a +14 par defaut, la plage reelle des fuseaux horaires)
    et en comparant leur erreur moyenne a celle du utdif actuellement
    enregistre.

    Ne verifie que les carottes AVEC visee solaire (les autres ont
    `declin` force a 0, aucun signal exploitable). Retourne None si le
    site n'a aucune carotte avec visee solaire, ou si `site.utdif` est
    deja optimal (ou assez proche - `improvement_threshold`, en degres,
    evite de signaler un gain marginal/du bruit de mesure plutot qu'une
    vraie erreur de saisie)."""
    sun_rows = [r for r in site.rows if not (r.azsun == 0.0 and r.hour == 0.0 and r.minute == 0.0)]
    if not sun_rows:
        return None
    if decli_igrf is None:
        decli_igrf = _site_decli_igrf(site)
    if candidates is None:
        candidates = [float(u) for u in range(-12, 15)]
    if site.utdif not in candidates:
        candidates = list(candidates) + [site.utdif]

    def _mean_abs_error(utdif: float) -> float:
        errors = [
            abs(compute_specimen_geometry(site, r, decli_igrf, utdif=utdif)[2] - decli_igrf)
            for r in sun_rows
        ]
        return sum(errors) / len(errors)

    current_error = _mean_abs_error(site.utdif)
    best_utdif, best_error = site.utdif, current_error
    for cand in candidates:
        err = _mean_abs_error(cand)
        if err < best_error - 1e-9:
            best_utdif, best_error = cand, err

    if best_utdif == site.utdif or (current_error - best_error) < improvement_threshold:
        return None

    return (
        f"site « {site.car_site} »: utdif={site.utdif:g} gives a mean "
        f"|local declination - IGRF| of {current_error:.1f} deg over "
        f"{len(sun_rows)} sun-shot core(s); utdif={best_utdif:g} would give "
        f"{best_error:.1f} deg - check the UTM/timezone difference for this site."
    )


def check_utdif_for_sites(sites: List[FieldSite]) -> List[str]:
    """Applique check_utdif a chaque site, retourne la liste des
    avertissements (un par site suspect, aucun si tout est coherent)."""
    warnings = []
    for site in sites:
        w = check_utdif(site)
        if w:
            warnings.append(w)
    return warnings


def build_pmag_records(sites: List[FieldSite], ispec: int) -> List[Pmag]:
    """Construit les enregistrements Pmag (mesures VIDES) pour chaque
    site/carotte/specimen (A, B, ... jusqu'a `ispec` lettres) - equivalent
    des ecritures unit 32 (lignes 339-406 de `orientation`), en reutilisant
    le format MODERNE .prmag (voir docstring module).

    `comment` (champ QC logiciel) : uniquement "err>3/5/10" (voir
    compute_specimen_geometry) - JAMAIS l'annotation de terrain.
    `obs` (champ observation) : annotation de terrain PAR CAROTTE (ex.
    "TR"/"PMO" - Testigo Roto/Posiblemente mal orientado, voir le
    commentaire "!..." de la ligne carotte), combinee a l'obs du SITE si
    present - demande explicite utilisateur ("fill the obs like TR at the
    specimen level in the prmag file") : jusqu'ici l'annotation de
    carotte finissait dans `comment`, avec le flag "err>X", et `obs`
    restait fige a la valeur du site (identique pour toutes ses
    carottes) - desormais `obs` varie reellement par specimen."""
    ispec = max(1, min(3, int(ispec)))
    letters = _LETTERS[:ispec]
    records: List[Pmag] = []

    for site in sites:
        decli_igrf = _site_decli_igrf(site)
        car_year = str(site.year - 1900 if site.year < 2000 else site.year - 2000).zfill(2)
        for row in site.rows:
            caz, plane_strike, declin, comm = compute_specimen_geometry(site, row, decli_igrf)
            obs = " ".join(t for t in (site.obs, row.comment) if t)
            for letter in letters:
                specimen_id = f"{car_year}{site.car_site}{row.car_sample}{letter}"
                records.append(Pmag(
                    id=specimen_id, cin=row.dip_core, caz=caz,
                    dip=site.plane_dip, str_=plane_strike,
                    norme="v", vol=_DEFAULT_VOLUME, com=comm,
                    lat=site.lat, rlong=site.rlong, altitude=site.altitude,
                    year=site.year, month=site.month, day=site.day,
                    hour=int(row.hour), minute=int(row.minute),
                    azmag=row.azmag, azsun=row.azsun,
                    outilorient="A" if site.outil.strip().upper() == "ASC" else "P",
                    magic_fm=site.fm, magic_age=site.age, magic_gc=site.gc,
                    magic_smt=site.smt, magic_li=site.litho, magic_loc=site.location,
                    magic_obs=obs,
                ))
    return records


def write_prmag_from_field_notes(records: List[Pmag], out_path: str) -> int:
    """Ecrit `records` (mesures vides) au format .prmag moderne - un bloc
    de 4 lignes d'entete + l'entete de mesures (SANS ligne de mesure,
    aucun instrument n'a encore rien mesure), separes par une ligne vide.

    Cree aussi le .pmagani compagnon (vide, juste l'en-tete) s'il
    n'existe pas deja - demande explicite utilisateur ("lors de la
    creation d'un .prmag ce serait peut-etre bien de creer un pmagani
    meme si il est vide") : pret a etre archive (AMS_Py "Archive ASC
    into .pmagani...") sans devoir d'abord passer par un .asc/.ANI pour
    qu'il existe. N'ecrase JAMAIS un .pmagani deja present (voir
    calcul.create_empty_pmagani_if_missing)."""
    blocks = []
    for p in records:
        blocks.append(_sample_header_block(p) + "\n" + _MEAS_HEADER)
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n\n".join(blocks) + "\n")
        f.flush()
        os.fsync(f.fileno())
    create_empty_pmagani_if_missing(ani_path_for(out_path))
    return len(records)


def write_ged_file(sites: List[FieldSite], ispec: int, out_path: str) -> int:
    """Ecrit UN SEUL fichier .ged combinant tous les sites - meme
    comportement que le Fortran en mode combine (`unfichierged='Y'`,
    voir docstring module) - memes formats de ligne 200/666 que la
    branche par-site du source :
    - ligne d'entete (ECRITE UNE SEULE FOIS, a la fin) : nb_specimens+2,
      "????" (pas de site unique), 12 espaces, lat/lon DU DERNIER SITE
      traite, puis le suffixe fixe "xxx xxx xxx xxx 12 0  3  90 "
      (constante du format Fortran, jamais calculee).
    - une ligne par specimen (tous sites confondus, dans l'ordre du
      fichier d'entree) : id (12 car.), azimut(int), plongee(int),
      "B0  ", direction_couche(int), pendage_couche(int), puis le
      suffixe fixe "0   0   00  0   0   0   0     ".
    Retourne le nombre de specimens ecrits."""
    ispec = max(1, min(3, int(ispec)))
    letters = _LETTERS[:ispec]
    lines: List[str] = []
    last_lat, last_rlong = 0.0, 0.0

    for site in sites:
        decli_igrf = _site_decli_igrf(site)
        car_year = str(site.year - 1900 if site.year < 2000 else site.year - 2000).zfill(2)
        for row in site.rows:
            caz, plane_strike, _declin, _comm = compute_specimen_geometry(site, row, decli_igrf)
            for letter in letters:
                specimen_id = f"{car_year}{site.car_site}{row.car_sample}{letter}"
                lines.append(
                    f"{specimen_id:<12}{int(caz):4d}{int(row.dip_core):4d}B0  "
                    f"{int(plane_strike):4d}{int(site.plane_dip):4d}"
                    "0   0   00  0   0   0   0     "
                )
        last_lat, last_rlong = site.lat, site.rlong

    header = (
        f"{len(lines) + 2:4d}{'????':<4}{'':12s}"
        f"{last_lat:7.1f}{last_rlong:7.1f}xxx xxx xxx xxx 12 0  3  90 "
    )
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(header + "\n")
        f.write("\n".join(lines) + ("\n" if lines else ""))
        f.flush()
        os.fsync(f.fileno())
    return len(lines)


def parse_ged_file(path: str, encoding: str = "utf-8") -> List[Pmag]:
    """Lit un fichier .ged AGICO (colonnes fixes - voir write_ged_file
    ci-dessus pour le format exact, cette fonction en est l'inverse) et
    construit une liste de specimens .prmag SANS mesure - demande
    explicite utilisateur ("un menu de creation de prmag file a partir
    d'un fichier agico .ged. et ensuite archive data") : couvre le cas ou
    seul un .ged est disponible (pas le fichier "field notes" brut que
    `parse_orientation_file`/`build_pmag_records` attendent), par exemple
    un .ged recu d'un collegue ou reconstitue depuis le logiciel de
    l'instrument AGICO lui-meme.

    Un .ged ne porte QUE l'identifiant du specimen et sa geometrie (azimut/
    pendage carotte, direction/pendage de la strate) - contrairement aux
    field notes brutes, il ne connait ni le site, ni la date, ni une
    visee solaire distincte, ni la geologie : ces champs restent a leurs
    valeurs par defaut (0.0/"n.d") dans le `Pmag` retourne, a completer
    ensuite via `Complete sample information...`. lat/lon proviennent de
    la ligne d'entete (une seule paire pour tout le fichier, meme limite
    que `write_ged_file` en mode combine) - 0.0/0.0 si absents/illisibles.

    Le `.prmag` ainsi cree n'a PAS de mesures : l'etape suivante est
    d'utiliser `Archive new laboratory measurements` (import_new_data.
    archive_new_measurements) pour y attacher les mesures reelles de
    l'instrument, qui exige des specimens DEJA presents dans le .prmag
    cible - exactement ce que cette fonction prepare."""
    with open(path, "r", encoding=encoding, errors="replace") as f:
        lines = [raw.rstrip("\n") for raw in f if raw.strip()]
    if len(lines) < 2:
        return []

    header = lines[0]
    try:
        lat = float(header[20:27])
        rlong = float(header[27:34])
    except (ValueError, IndexError):
        lat = rlong = 0.0

    records: List[Pmag] = []
    for line in lines[1:]:
        specimen_id = line[0:12].strip()
        if not specimen_id:
            continue
        try:
            caz = float(line[12:16])
            cin = float(line[16:20])
            plane_strike = float(line[24:28])
            plane_dip = float(line[28:32])
        except (ValueError, IndexError):
            continue
        records.append(Pmag(
            id=specimen_id, cin=cin, caz=caz,
            dip=plane_dip, str_=plane_strike,
            norme="v", vol=_DEFAULT_VOLUME,
            lat=lat, rlong=rlong,
        ))
    return records


def write_diagnostics_report(sites: List[FieldSite], out_path: str) -> int:
    """Equivalent du fichier ".res" (unit 57, format 100) - une ligne de
    diagnostic par specimen : azimut geographique, declinaison locale
    mesuree, declinaison IGRF (comparaison/verification)."""
    lines = []
    n = 0
    for site in sites:
        decli_igrf = _site_decli_igrf(site)
        car_year = str(site.year - 1900 if site.year < 2000 else site.year - 2000).zfill(2)
        for row in site.rows:
            caz, _plane_strike, declin, comm = compute_specimen_geometry(site, row, decli_igrf)
            azgeo = _wrap360(caz - 90.0)
            specimen_id = f"{car_year}{site.car_site}{row.car_sample}A"
            note = f"  {comm}" if comm else ""
            lines.append(
                f"{specimen_id:<12s} heure:{int(row.hour):2d} {int(row.minute):2d}  "
                f"az geo:{azgeo:6.1f}  local Decli:{declin:6.1f}  igrf:{decli_igrf:6.1f}{note}"
            )
            n += 1
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + ("\n" if lines else ""))
        f.flush()
        os.fsync(f.fileno())
    return n
