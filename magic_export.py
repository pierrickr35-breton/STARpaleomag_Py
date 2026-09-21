"""
Port de `export2magic` ("export Rennes to Magic", fichiers_magic.f) - mode
"classique" (ichoixexport==1) : sites.txt, samples.txt, specimens.txt,
measurements.txt, locations.txt - PLUS les resultats de paleointensite deja
archives dans .pmagint (ichoixexport==2 cote Fortran, `export_int_2_magic` -
fusionnes ici dans le MEME specimens.txt plutot qu'un mode/fichier separe,
voir export_to_magic - demande explicite utilisateur "in export to Magic,
it is not asking for Paleointensity"). Magnetostratigraphie
(ichoixexport==3) reste HORS PERIMETRE (a ajouter plus tard si besoin).
L'export AMS (ichoixexport 4/5) est du code mort dans le Fortran (un
`return` le rend inatteignable) et n'est pas porte.

Ecarts deliberes par rapport au Fortran (tous discutes et valides) :
- Regroupement site/echantillon par TRI EXPLICITE sur magic_site/
  magic_sample (le Fortran detecte un "nouveau site" seulement par
  changement de valeur consecutive, ce qui suppose les donnees deja
  triees - un tri explicite est plus robuste, meme si les .ren sont
  souvent deja tries).
- Les 2 bugs reperes dans le Fortran (le flag "cooling rate" ne refletant
  que la DERNIERE mesure d'un echantillon au lieu de "au moins une", et un
  index de boucle perime pour la classification R/V d'un controle pTRM)
  sont CORRIGES ici, pas reproduits.
- Formatage numerique simple (`repr`/format Python standard), pas les
  conventions exactes du Fortran list-directed (`write(x,*)`) - le fichier
  reste un TSV MagIC valide, juste pas byte-identique a une sortie Fortran.
- `locations.txt` : continent/pays/region ne sont PAS des litteraux codes
  en dur ("Chili"/"Amerique du Sud" specifiques au jeu de donnees
  d'origine) - ce sont des parametres optionnels de `export_to_magic`
  (vide par defaut), a fournir une fois par export plutot qu'a completer
  a la main apres coup dans chaque fichier.

Champs Site/Sample/Fm/Age/GC/SMT/Li/Loc deja disponibles directement sur
chaque SelectedSample via testlect.decode_roche (magic_site, magic_sample,
magic_fm, magic_age, magic_gc, magic_smt, magic_li, magic_loc, magic_obs) -
pas besoin de re-parser la ligne roche ici.
"""

import math
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from selection import SelectedSample, Measurement, polere, corfor, corpen
from calcul import FitResult, AniTensor, AniMeanTensor
from complete_sample_info import _read_table_rows

# ---------------------------------------------------------------------------
# Age : equivalent de `checkage` (fichiers_magic.f:2225), etendu pour
# tolerer le format "low - high unit" produit par extract_magic.py (import
# MagIC, qui copie verbatim le champ age/geologic_age d'origine) en plus
# des delimiteurs Fortran natifs #/_/& des donnees Rennes brutes.
# ---------------------------------------------------------------------------

_AGE_UNIT_PATTERNS = [
    ("Ga", "Ga"),
    ("Ma", "Ma"),
    ("ka", "ka"),
    ("Cal AD", "Years Cal AD (+/-)"),
    ("Cal BP", "Years Cal BP"),
    ("AD", "Years AD (+/-)"),
    ("BP", "Years BP"),
]


def parse_age(age_str: str) -> Tuple[str, str, str, str, str]:
    """Retourne (age, age_sigma, age_low, age_high, age_unit)."""
    s = (age_str or "").strip()
    age_unit = ""
    for token, magic_unit in _AGE_UNIT_PATTERNS:
        if token in s:
            age_unit = magic_unit
            s = s.replace(token, " ").strip()
            break

    age = age_sigma = age_low = age_high = ""
    if "#" in s:
        a, b = s.split("#", 1)
        age, age_sigma = a.strip(), b.strip()
    elif "&" in s:
        parts = [p.strip() for p in s.split("&")]
        if len(parts) >= 3:
            age, age_low, age_high = parts[0], parts[1], parts[2]
        else:
            age = s
    elif "_" in s:
        a, b = s.split("_", 1)
        age_low, age_high = a.strip(), b.strip()
    elif re.search(r"\s-\s", s):
        a, b = re.split(r"\s-\s", s, maxsplit=1)
        age_low, age_high = a.strip(), b.strip()
    else:
        age = s

    return age, age_sigma, age_low, age_high, age_unit


# ---------------------------------------------------------------------------
# Petits helpers geometriques/formatage
# ---------------------------------------------------------------------------

def _wrap_lon(rlong: float) -> float:
    return 360.0 + rlong if rlong < 0.0 else rlong


def _bed_dip_direction(str_: float, dip: float) -> float:
    if dip == 0.0:
        return 0.0
    v = str_ + 90.0
    if v > 360.0:
        v -= 360.0
    return v


def _timestamp(ech: SelectedSample) -> str:
    if not (ech.year and ech.month and ech.day):
        return ""
    return (f"{int(ech.year):04d}-{int(ech.month):02d}-{int(ech.day):02d}T"
            f"{int(ech.hour):02d}:{int(ech.minute):02d}:00")


def _num(x: float, decimals: int = 5) -> str:
    return f"{x:.{decimals}f}"


def _num_sci(x: float, sig: int = 6) -> str:
    """Formatage scientifique (`sig` chiffres significatifs), comme les
    vrais fichiers MagIC pour ce type de grandeur (ex. magn_moment=
    "3.34e-05", volume="1e-06" - verifie sur magic_contribution_16645.txt/
    20340.txt) - demande explicite utilisateur ("in the export to magic,
    the volume and magnetic moment are wrong") : `_num()` (virgule fixe)
    tronquait silencieusement ces valeurs a "0.000000" pour un moment
    typique (~1e-8 Am2) - verifie sur donnees reelles (old_pmag.ren,
    17968 mesures) : 77 mesures perdaient COMPLETEMENT leur moment
    (arrondi exact a zero), 68% avaient moins de 3 chiffres significatifs
    a 6 decimales fixes. Utilise pour toute grandeur physique de magnitude
    tres variable (moment, volume, masse, susceptibilite) - PAS pour les
    angles/coordonnees (dec/inc/lat/lon/...), ou `_num()` (virgule fixe)
    reste approprie et est conserve."""
    return f"{x:.{sig}g}"


# ---------------------------------------------------------------------------
# sites.txt / locations.txt
# ---------------------------------------------------------------------------

_SITES_HEADER = [
    "site", "citations", "location", "geologic_classes", "lat", "lon",
    "elevation", "formation", "lithologies", "geologic_types",
    "bed_dip_direction", "bed_dip", "geographic_precision", "method_codes",
    "age", "age_sigma", "age_low", "age_high", "age_unit",
    "samples", "specimens", "result_quality", "dir_comp_name",
    "dir_tilt_correction", "dir_dec", "dir_inc", "dir_alpha95", "dir_r",
    "dir_k", "dir_n_specimens", "dir_n_specimens_lines",
    "dir_n_specimens_planes", "dir_polarity", "dir_nrm_origin",
    "vgp_lat", "vgp_lon", "vgp_dp", "vgp_dm",
    # AMS de site (moyenne DEJA calculee, voir anisotropy_site_magic_fields)
    # - demande explicite utilisateur ("l'exportation de l'AMS (niveau
    # specimen et sites)... vont aussi dans les fichiers sites.txt").
    "aniso_type", "aniso_tilt_correction", "aniso_v1", "aniso_v2", "aniso_v3",
    "aniso_p", "aniso_pp", "aniso_t", "aniso_l", "aniso_f",
    "aniso_perc", "aniso_total", "aniso_ll", "aniso_ff", "aniso_vg", "aniso_fl",
    # ecart avec AMS_Py.ams_stats.magic_site_aniso_fields comble ici (voir
    # anisotropy_site_magic_fields plus bas) : AMS_Py ecrit cette mise en
    # garde eta/zeta directement dans la colonne `description` de son
    # export sites.txt natif, ce que ce module ne faisait pas encore -
    # sans cette colonne, l'avertissement ne restait que dans le
    # docstring du code, invisible pour qui exporte depuis
    # STARpaleomag_Py plutot que depuis AMS_Py.
    "description",
]

_LOCATIONS_HEADER = [
    "location", "geologic_classes", "lithologies", "age", "age_sigma",
    "age_low", "age_high", "age_unit", "lat_n", "lat_s", "lon_e", "lon_w",
    "citations", "continent_ocean", "country", "description",
    "location_type", "region",
]

def _comp_name(component: str) -> str:
    """"Comp_A"/"Comp_B"/"Comp_C"... a partir de FitResult.component (A/B/C,
    etiquette de composante de magnetisation - voir son docstring) - PAS de
    `numcomp` (numero d'ajustement PCA 1-9, un concept different, deja
    signale a eviter ici : "I think it is best not to use the numcomp of
    individual samples for the mean"). Demande explicite utilisateur
    ("now that we also have component A or B, instead of Chrm_1 write
    Comp_A") : uniformise dir_comp/dir_comp_name (specimens.txt ET
    sites.txt) sur cette meme etiquette, plutot que Chrm_<numcomp>."""
    return f"Comp_{(component or 'A').strip().upper()}"


def load_site_metadata_table(path: str, encoding: str = "utf-8") -> Dict[str, Dict[str, str]]:
    """Table de metadonnees de site (formation/age/geologic_classes/
    geologic_types/lithologies/location), keyee par nom de site - demande
    explicite utilisateur ("can we also let the site-only path pull from
    a complement table... so those fields aren't just blank") : comble
    formation/lithologies/geologic_classes/geologic_types/age pour un
    site publie sur MagIC SANS specimen charge (voir build_sites_rows),
    qui n'a sinon aucune source pour ces champs (contrairement a un site
    avec specimens, qui les lit sur `ech.magic_fm`/etc.).

    Meme lecteur que complete_sample_info._read_table_rows (delimiteur
    tabulation OU virgule auto-detecte, colonnes normalisees espace->
    underscore) - accepte donc directement un fichier au format
    Complement_tect.txt (colonne 'site' REPETEE : ancien nom -> nouveau
    nom -> keyee ici par le DERNIER 'site' - le nouveau nom, celui sous
    lequel un site sans specimen est archive/nomme, voir build_site_mean_
    result) aussi bien qu'un fichier a une seule colonne 'site' (ex.
    site_locality_complement.txt). age_low/age_high/age_unit sont lus
    TELS QUELS s'ils sont presents (evite l'aller-retour par parse_age -
    inutile ici, sites.txt de MagIC a deja des colonnes age_low/age_high/
    age_unit separees)."""
    header, rows = _read_table_rows(path, encoding=encoding)
    if "site" not in header:
        raise ValueError("Table file must have a 'site' column header.")
    site_indices = [i for i, h in enumerate(header) if h == "site"]
    key_idx = site_indices[-1]
    table: Dict[str, Dict[str, str]] = {}
    for row in rows:
        if key_idx >= len(row):
            continue
        site = row[key_idx].strip()
        if not site:
            continue
        record = {
            col: row[i].strip()
            for i, col in enumerate(header)
            if i not in site_indices and i < len(row)
        }
        table[site] = record
    return table


def _apply_site_metadata(row: Dict[str, str], meta: Optional[Dict[str, str]]) -> None:
    """Comble UNIQUEMENT les champs encore vides de `row` (jamais n'ecrase
    une valeur deja connue, qu'elle vienne d'un specimen ou de `mean`) -
    voir load_site_metadata_table. `age` (chaine combinee) reconstruite
    depuis age_low/age_high/age_unit si fournis separement et qu'aucune
    colonne 'age' directe n'existe - meme convention que complete_sample_
    info._normalize_age, dans l'AUTRE sens (ici sites.txt veut le detail,
    pas la chaine combinee, donc age_low/age_high/age_unit priment)."""
    if not meta:
        return
    for field in ("formation", "lithologies", "geologic_classes", "geologic_types",
                  "location", "age_low", "age_high", "age_unit"):
        if not row.get(field) and meta.get(field):
            row[field] = meta[field]
    if not row.get("age_low") and not row.get("age_high") and meta.get("age") and not row.get("age"):
        row["age"] = meta["age"]
    # lat/lon - demande implicite (site sans specimen dont le VGP archive
    # est lui-meme a 0.0/0.0, voir build_sites_rows) : `_num()` valide
    # que la table fournit bien un nombre avant de l'utiliser, plutot que
    # de propager une valeur texte invalide dans sites.txt.
    if not row.get("lat") and meta.get("lat"):
        try:
            row["lat"] = _num(float(meta["lat"]))
        except ValueError:
            pass
    if not row.get("lon") and meta.get("lon"):
        try:
            row["lon"] = _num(_wrap_lon(float(meta["lon"])))
        except ValueError:
            pass


def _apply_site_aniso(row: Dict[str, str], mean: Optional[AniMeanTensor]) -> None:
    """Fusionne les colonnes AMS de site (anisotropy_site_magic_fields)
    dans `row`, et ajoute le code de methode du tenseur (ex. LP-AN-TRM)
    a `method_codes` s'il n'y est pas deja - meme principe que
    build_specimens_rows pour le niveau specimen (paleointensite/
    anisotropie), demande explicite utilisateur ("l'exportation de
    l'AMS (niveau specimen et sites)")."""
    if mean is None:
        return
    row.update(anisotropy_site_magic_fields(mean))
    _aniso_type, aniso_method = _ANISO_MAGIC_INFO.get(mean.code2, ("AMS", "LP-X"))
    codes = row.get("method_codes", "")
    if aniso_method not in codes.split(":"):
        row["method_codes"] = f"{codes}:{aniso_method}" if codes else aniso_method


_C_SUFFIX_RE = re.compile(r"_[a-z]{1,2}$")


def _strip_c_suffix(c: str) -> str:
    """Retire le suffixe anti-collision "_<lettre(s)>" d'un `c` (voir
    calcul._next_specimen_c : une lettre a-z, puis deux au-dela de 26
    collisions pour le meme specimen - jamais observe en pratique mais
    couvert quand meme) pour retrouver le nom de specimen d'origine SANS
    avoir besoin du resultat individuel correspondant charge en memoire -
    demande explicite utilisateur ("it should resolve also with the mean
    too as the codes are used to select the true specimen name"). Ne
    retire RIEN si `c` ne suit pas ce format exact (ex. un `c` numerique
    d'un ancien fichier .pmagres - voir _next_specimen_c, remplace un
    entier aleatoire 0-99999 dans les fichiers ecrits AVANT ce format)."""
    return _C_SUFFIX_RE.sub("", c)


def _site_mean_rows(site: str, results: List[FitResult]) -> List[Dict[str, str]]:
    """Equivalent de `magicmeanres` : TOUTES les moyennes de site (cat1=='F'
    via id "mean: <site>") correspondant a `site` (magic_site, ou son
    equivalent tronque a 6 caracteres comme la convention interne des id
    "mean: XXXXXX") - UNE ligne par resultat trouve (IS ET TC, PLUS
    plusieurs composantes A/B/C... si presentes), PAS seulement le
    premier - demande explicite utilisateur ("the export do not include
    the Tilt corrected results that should appear on a second line") :
    un site archive en IS ET en TC (voir _archive_fisher_mean, qui
    archive les deux automatiquement) n'exportait jusqu'ici QUE le
    premier trouve (le Fortran d'origine, magicmeanres, ne connaissait
    lui-meme qu'une seule orientation par appel - jamais deux archivees
    pour le meme site).

    BUGS REELS corriges ici (demande explicite utilisateur, "the export
    do not decode well the lines and planes") : `dir_k` lisait
    `r.par2_mean` (jamais renseigne pour une moyenne de site - toujours
    0.0) au lieu de `r.tx[0]` (ou k EST reellement stocke, voir
    build_site_mean_result: `tx=(stats.k, 0.0)`) ; `dir_n_specimens_lines`/
    `dir_n_specimens_planes` lisaient `r.tx[0]`/`r.tx[1]` (donc le K
    recopie comme "nombre de lignes", et 0.0 constant comme "nombre de
    plans") au lieu de `r.n_lines`/`r.n_planes` (les VRAIS comptes, voir
    _MEAN_FIELDS colonne "L/P")."""
    site6 = (site or "").strip()[:6]
    rows = []
    for r in results:
        if r.id[:5] != "mean:":
            continue
        if r.id[6:12].strip() != site6:
            continue
        tilt = "0" if r.par3_mean == 2.0 else "100"
        comp = _comp_name(r.component)
        # `r.liste` porte les `c` des resultats combines (ex.
        # "96CC0301B_a"), PAS leur nom de specimen - `c` est
        # l'identifiant anti-collision "<specimen>_<lettre>" attribue a
        # l'archivage (voir _next_specimen_c), jamais le nom de specimen
        # MagIC lui-meme - demande explicite utilisateur ("only the
        # number of the specimen is included. not with the extension _a
        # or b", puis "it should resolve also with the mean too as the
        # codes are used to select the true specimen name") : resout
        # chaque `c` en priorite vers le resultat individuel correspondant
        # (retrouve par egalite de `c` parmi `results` - le plus fiable,
        # marche meme si un specimen contient lui-meme un underscore) ;
        # si ce resultat individuel n'est PAS charge (ex. moyenne seule
        # selectionnee, mode 'm' plutot que 's'), retombe sur le retrait
        # du suffixe "_<lettre(s)>" du `c` lui-meme via _strip_c_suffix -
        # la convention de _next_specimen_c (1 lettre, puis 2 au-dela de
        # 26 collisions) suffit a elle seule a retrouver le nom. Deduplique
        # en gardant l'ordre (meme specimen combine deux fois - deux
        # composantes A/B distinctes du meme specimen physique, par
        # exemple - ne doit apparaitre qu'une fois dans la liste).
        by_c = {str(other.c): other.id for other in results if other.id[:5] != "mean:"}
        codes = [t.strip() for t in r.liste.replace("codes:", "").split(":") if t.strip()]
        specimens = ":".join(
            dict.fromkeys(by_c.get(c) or _strip_c_suffix(c) for c in codes)
        )
        rows.append({
            "specimens": specimens,
            "result_quality": "g",
            "dir_comp_name": comp,
            "dir_tilt_correction": tilt,
            "dir_dec": _num(r.dec, 1),
            "dir_inc": _num(r.inc, 1),
            "dir_alpha95": _num(r.mad, 1),
            "dir_k": _num(r.tx[0], 1),
            "dir_n_specimens": str(r.nb),
            "dir_n_specimens_lines": str(r.n_lines) if r.n_lines >= 0 else "",
            "dir_n_specimens_planes": str(r.n_planes) if r.n_planes >= 0 else "",
            "vgp_lat": _num(r.par4, 1),
            "vgp_lon": _num(r.par5, 1),
            "vgp_dp": _num(r.vgp_dp, 1) if r.vgp_dp else "",
            "vgp_dm": _num(r.vgp_dm, 1) if r.vgp_dm else "",
        })
    return rows


def build_sites_rows(
    samples: List[SelectedSample], results: List[FitResult],
    site_metadata: Optional[Dict[str, Dict[str, str]]] = None,
    aniso_mean_tensors: Optional[Dict[str, AniMeanTensor]] = None,
    include_site_only_means: bool = True,
) -> List[List[str]]:
    """`include_site_only_means` : False quand l'utilisateur n'exporte
    qu'une PARTIE du fichier - la seconde boucle plus bas (sites connus
    seulement par une moyenne "mean:", sans specimen dans `samples`)
    ajoutait sinon la moyenne de TOUS les sites du .pmagres, meme ceux
    hors de la selection (signale par l'utilisateur : "meme si on ne
    selectionne que qq echantillons, l'export de magic tente d'exporter
    tout le fichier").

    `site_metadata` (voir load_site_metadata_table) : ne COMBLE que les
    champs formation/lithologies/geologic_classes/geologic_types/age_low/
    age_high/age_unit/location encore VIDES (voir _apply_site_metadata) -
    utile pour un site SANS specimen (aucune autre source, voir plus bas
    dans cette fonction), mais applique aussi aux sites avec specimens
    dont l'un de ces champs n'a jamais ete renseigne (ex. .prmag jamais
    passe par "Complete sample information...").

    UNE LIGNE PAR RESULTAT de moyenne trouve pour le site (voir
    _site_mean_rows) - PAS une seule ligne meme quand le site porte une
    moyenne IS ET une moyenne TC (ou plusieurs composantes) - demande
    explicite utilisateur ("the export do not include the Tilt corrected
    results that should appear on a second line").

    `aniso_mean_tensors` (site -> AniMeanTensor 'A0' deja calcule, voir
    anisotropy_site_magic_fields) : demande explicite utilisateur
    ("l'exportation de l'AMS (niveau specimen et sites)"), meme principe
    que `aniso_tensors` niveau specimen (build_specimens_rows) - fusionne
    dans TOUTES les lignes du site (une moyenne AMS de site n'est pas
    liee a l'orientation IS/TC d'une moyenne directionnelle particuliere,
    contrairement a dir_dec/dir_inc)."""
    rows = []
    seen = []
    for ech in samples:
        site = ech.magic_site.strip()
        if not site or site in seen:
            continue
        seen.append(site)

        age, age_sigma, age_low, age_high, age_unit = parse_age(ech.magic_age)
        elevation = "" if ech.altitude <= 0.0 else _num(ech.altitude, 1)
        base_row = {
            "site": site,
            "citations": "This study",
            "location": ech.magic_loc,
            "geologic_classes": ech.magic_gc,
            "lat": _num(ech.lat),
            "lon": _num(_wrap_lon(ech.rlong)),
            "elevation": elevation,
            "formation": ech.magic_fm,
            "lithologies": ech.magic_li,
            "geologic_types": ech.magic_smt,
            "bed_dip_direction": _num(_bed_dip_direction(ech.str_, ech.dip), 1),
            "bed_dip": _num(ech.dip, 1),
            "geographic_precision": "0.0001",
            "method_codes": "FS-FD:GE-WGS84:FS-LOC-GPS",
            "age": age, "age_sigma": age_sigma, "age_low": age_low,
            "age_high": age_high, "age_unit": age_unit,
        }
        _apply_site_metadata(base_row, (site_metadata or {}).get(site))
        _apply_site_aniso(base_row, (aniso_mean_tensors or {}).get(site))
        means = _site_mean_rows(site, results)
        if not means:
            rows.append([base_row.get(col, "") for col in _SITES_HEADER])
        for mean in means:
            row = dict(base_row)
            row.update(mean)
            rows.append([row.get(col, "") for col in _SITES_HEADER])

    # Sites connus SEULEMENT via une moyenne archivee ("mean: <site>"),
    # SANS specimen charge dans `samples` - demande explicite utilisateur
    # ("in the case of legacy files, it might be interesting to publish
    # the site with its mean result (the published assuming that we lost
    # some data). can we have a site in prmag without data?") : un site
    # dont les mesures brutes sont perdues mais dont le resultat publie
    # (direction moyenne, VGP) est connu (ex. transcrit depuis une table
    # publiee - voir build_site_mean_result/le tableau du papier) reste
    # publiable sur MagIC - le format sites.txt de MagIC n'exige PAS de
    # specimens/measurements pour porter un resultat de site. lat/lon/
    # bed_dip proviennent alors du resultat lui-meme (r.lat/r.rlong/
    # r.dip/r.str_, voir build_site_mean_result) plutot que d'un
    # specimen ; formation/lithologies/age/geologic_classes n'ont ICI
    # AUCUNE source specimen - viennent UNIQUEMENT de `site_metadata` si
    # fourni (voir load_site_metadata_table/_apply_site_metadata),
    # restent vides sinon plutot que d'inventer une valeur.
    for r in (results if include_site_only_means else []):
        if r.id[:5] != "mean:":
            continue
        site = r.id[6:].strip()
        if not site or site in seen:
            continue
        seen.append(site)
        means = _site_mean_rows(site, results)
        if not means:
            continue
        base_row = {
            "site": site,
            "citations": "This study",
            "lat": _num(r.lat) if (r.lat or r.rlong) else "",
            "lon": _num(_wrap_lon(r.rlong)) if (r.lat or r.rlong) else "",
            "bed_dip_direction": _num(_bed_dip_direction(r.str_, r.dip), 1) if (r.dip or r.str_) else "",
            "bed_dip": _num(r.dip, 1) if (r.dip or r.str_) else "",
            "geographic_precision": "0.0001",
            "method_codes": "FS-FD:GE-WGS84:FS-LOC-GPS",
        }
        _apply_site_metadata(base_row, (site_metadata or {}).get(site))
        _apply_site_aniso(base_row, (aniso_mean_tensors or {}).get(site))
        for mean in means:
            row = dict(base_row)
            row.update(mean)
            rows.append([row.get(col, "") for col in _SITES_HEADER])
    return rows


_SITE_METADATA_PREVIEW_HEADER = [
    "site", "lat", "lon", "geologic_classes", "geologic_types", "lithologies", "formation",
]


def build_site_metadata_preview_rows(samples: List[SelectedSample]) -> List[Dict[str, str]]:
    """Une ligne PAR SITE (premier specimen rencontre pour ce site, MEME
    tri/convention que build_sites_rows - donc les memes valeurs que
    celles qui finiront dans sites.txt), avec les champs que MagIC exige
    au niveau site et que ce portage sait completer depuis un fichier
    complement (site/lat/lon/geologic_classes/geologic_types/
    lithologies/formation) - "" quand non renseigne sur le specimen.

    Sert de table d'apercu EXPORTABLE (voir
    app.ouvrir_export_magic_dialog) affichee/proposee AVANT le reste du
    dialogue d'export, plutot que de ne demander un fichier complement
    qu'a la toute fin sans que l'utilisateur ait pu voir ce qui manque -
    demande explicite utilisateur ("ce serait bien d'exporter d'abord un
    tableau avec site, lat lon Geologic classes, types lithology en
    invitant l'utilisateur d'avoir le fichier le plus adequat... plutot
    que de poser cette question a la fin"). Colonnes DELIBEREMENT
    identiques aux noms lus par load_site_metadata_table : le fichier
    ecrit ici (voir write_site_metadata_preview) peut etre complete a la
    main puis redonne tel quel comme fichier complement, sans renommer
    de colonne."""
    ordered = sorted(
        samples, key=lambda e: (e.magic_site.strip(), e.magic_sample.strip(), e.id.strip()))
    rows: List[Dict[str, str]] = []
    seen: set = set()
    for ech in ordered:
        site = ech.magic_site.strip()
        if not site or site in seen:
            continue
        seen.add(site)
        has_latlon = ech.lat != 0.0 or ech.rlong != 0.0
        rows.append({
            "site": site,
            "lat": _num(ech.lat) if has_latlon else "",
            "lon": _num(_wrap_lon(ech.rlong)) if has_latlon else "",
            "geologic_classes": ech.magic_gc,
            "geologic_types": ech.magic_smt,
            "lithologies": ech.magic_li,
            "formation": ech.magic_fm,
        })
    return rows


def write_site_metadata_preview(path: str, rows: List[Dict[str, str]]) -> None:
    """Ecrit la table d'apercu (build_site_metadata_preview_rows) en
    texte tabule - meme format (tabulation, en-tete en 1ere ligne) que
    ce que load_site_metadata_table/_read_table_rows sait relire, pour
    que ce fichier soit directement reutilisable comme fichier
    complement une fois complete a la main."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(_SITE_METADATA_PREVIEW_HEADER) + "\n")
        for row in rows:
            f.write("\t".join(row.get(col, "") for col in _SITE_METADATA_PREVIEW_HEADER) + "\n")


def build_locations_rows(
    samples: List[SelectedSample],
    continent_ocean: str = "", country: str = "", region: str = "",
    description: str = "", location_type: str = "Region",
) -> List[List[str]]:
    """Un `location` MagIC peut regrouper plusieurs sites ; ici on genere
    UNE ligne par valeur distincte de `magic_loc`, avec les bornes
    lat/lon agregees sur tous les sites qui y sont rattaches - le reste
    (continent/pays/region/description/type) est fourni par l'appelant
    (ex. saisi une fois au moment de l'export), PAS devine ou code en dur.

    `location_type` par defaut "Region" - demande explicite utilisateur
    ("in the export (file location) Location Type, location_type by
    default use Region"). Champ REQUIS par le modele de donnees MagIC
    (voir data_model3 : `validations: [cv("location_type"), required()]`,
    "Region" est le premier exemple donne) - laisser vide (comportement
    precedent) produisait un fichier locations.txt invalide au sens du
    modele MagIC."""
    groups: Dict[str, List[SelectedSample]] = {}
    for ech in samples:
        loc = ech.magic_loc.strip()
        if not loc:
            continue
        groups.setdefault(loc, []).append(ech)

    rows = []
    for loc in sorted(groups):
        members = groups[loc]
        lats = [m.lat for m in members]
        lons = [_wrap_lon(m.rlong) for m in members]
        gc = next((m.magic_gc for m in members if m.magic_gc), "")
        li = next((m.magic_li for m in members if m.magic_li), "")
        age, age_sigma, age_low, age_high, age_unit = parse_age(
            next((m.magic_age for m in members if m.magic_age), "")
        )
        row = {
            "location": loc,
            "geologic_classes": gc,
            "lithologies": li,
            "age": age, "age_sigma": age_sigma, "age_low": age_low,
            "age_high": age_high, "age_unit": age_unit,
            "lat_n": _num(max(lats)), "lat_s": _num(min(lats)),
            "lon_e": _num(max(lons)), "lon_w": _num(min(lons)),
            "citations": "This study",
            "continent_ocean": continent_ocean, "country": country,
            "description": description, "location_type": location_type,
            "region": region,
        }
        rows.append([row.get(col, "") for col in _LOCATIONS_HEADER])
    return rows


# ---------------------------------------------------------------------------
# samples.txt
# ---------------------------------------------------------------------------

# Vitesses de refroidissement de l'experience de vitesse de refroidissement
# (pas 'L' lent, 'Q' rapide) - valeurs fournies par l'utilisateur ("vitesse
# lente en laboratoire 1 K/min, rapide 10 K/min, pour la vitesse geologique
# mettre la meme vitesse que la lente au labo mais par Myr"). thellier_gui
# (PmagPy) ne calcule sa correction de vitesse de refroidissement que s'il
# trouve "<vitesse>:K/min" dans la description des mesures ET la vitesse
# ancienne (K/Ma, colonne `cooling_rate` de samples.txt).
_LAB_SLOW_COOLING_K_PER_MIN = 1.0
_LAB_FAST_COOLING_K_PER_MIN = 10.0
# Vitesse "geologique" = la vitesse LENTE de laboratoire exprimee en K/Ma
# (refroidissement de fours archeologiques, du meme ordre que le refroidissement
# lent de labo) : correction demandee par l'utilisateur ("transformer la
# vitesse lente en equivalent par /Ma" - la premiere version mettait 1 K/Ma,
# soit un rapport ~5e11 avec le labo). Meme conversion que thellier_gui
# (get_data : K/Ma / (1e6*365*24*60) = K/min), donc le rapport
# labo/geologique vaut exactement 1 pour la vitesse lente.
_GEOLOGICAL_COOLING_K_PER_MYR = _LAB_SLOW_COOLING_K_PER_MIN * (1.0e6 * 365.0 * 24.0 * 60.0)

_SAMPLES_HEADER = [
    "citations", "sample", "site", "geologic_classes", "lithologies",
    "geologic_types", "lat", "lon", "height", "timestamp", "orientation_quality",
    "azimuth", "dip", "bed_dip_direction", "bed_dip", "method_codes",
    "cooling_rate",
]


def build_samples_rows(samples: List[SelectedSample]) -> List[List[str]]:
    rows = []
    seen = []
    cooling_samples = {
        (e.magic_sample.strip() or e.id)
        for e in samples if any(m.cod1 in ("L", "Q") for m in e.mesures)
    }
    for ech in samples:
        # Bloc "specimen: n.d / sample: n.d" (testlect les lit comme id=""
        # - voir "specimen: n.d" dans le .prmag), sans aucune mesure : PAS
        # un vrai specimen, juste un porteur de metadonnees de SITE (lat/
        # lon/bed_dip/formation/age/geologie) pour un site dont les
        # donnees brutes sont perdues - demande explicite utilisateur
        # ("in the prmag... can we have a site without data?", "is this
        # OK?"). Contribue a build_sites_rows (qui l'utilise deja plus
        # haut) mais PAS ici : un nom d'echantillon vide est un champ cle
        # invalide pour MagIC (samples.txt exige un nom non vide et
        # unique) - verifie concretement (produisait une ligne avec
        # sample="").
        if not ech.id.strip() and not ech.mesures:
            continue
        sample = ech.magic_sample.strip() or ech.id
        if sample in seen:
            continue
        seen.append(sample)

        azimuth = ech.caz - 90.0
        if azimuth < 0.0:
            azimuth += 360.0
        if ech.hour == 0 and ech.minute == 0 and ech.azsun == 0.0:
            method_codes = "SO-POM:SO-CMD-NORTH"
        else:
            method_codes = "SO-POM:SO-SUN"

        height = "" if ech.stratigraphic_height is None else _num(ech.stratigraphic_height, 2)
        row = [
            "This study", sample, ech.magic_site, ech.magic_gc, ech.magic_li,
            ech.magic_smt, _num(ech.lat), _num(_wrap_lon(ech.rlong)), height,
            _timestamp(ech), "g", _num(azimuth, 1), _num(-ech.cin, 1),
            _num(_bed_dip_direction(ech.str_, ech.dip), 1), _num(ech.dip, 1),
            method_codes,
            f"{_GEOLOGICAL_COOLING_K_PER_MYR:.6g}" if sample in cooling_samples else "",
        ]
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# specimens.txt : ligne de base + lignes directionnelles (magicdirres)
# ---------------------------------------------------------------------------

_SPECIMENS_HEADER = [
    "specimen", "citations", "sample", "geologic_classes", "lithologies",
    "geologic_types", "azimuth", "dip", "volume", "weight", "method_codes",
    "dir_tilt_correction", "result_quality", "dir_nrm_origin",
    "meas_step_min", "meas_step_max", "meas_step_unit", "dir_comp",
    "dir_dec", "dir_inc", "dir_n_measurements", "dir_mad_anc", "dir_mad_free",
    "int_abs", "int_abs_sigma", "int_corr", "int_treat_dc_field",
    "int_n_measurements", "int_f", "int_g", "int_q", "int_mad_free",
    "int_dang", "int_b", "int_b_sigma", "int_b_beta", "int_k", "int_k_sse",
    "int_fvds", "int_frac", "int_gmax", "int_n_ptrm", "int_gamma",
    "int_corr_aniso", "int_corr_cooling_rate",
    "aniso_type", "aniso_s", "aniso_tilt_correction",
    "aniso_s_n_measurements", "aniso_s_sigma",
    "aniso_ftest", "aniso_ftest12", "aniso_ftest23", "aniso_ftest_quality",
    "description",
]

# code2 -> (aniso_type, method_codes) MagIC - MEME convention/memes 3
# entrees que AMS_Py.app._ANI_MAGIC_INFO (A0=ATRM, F0=AARM, N0=AMS -
# susceptibilite basse frequence, le cas le plus courant pour un .pmagani
# importe depuis AMS_Py/un ASC Agico, voir rn01_15_converted.pmagani reel
# - "prendre les A0" ne couvrait QUE le premier cas, code2 REELLEMENT
# rencontre etant souvent 'N0' - demande explicite utilisateur en
# verifiant "how to join the two exports from prmag & pmagres with the
# Anisotropy" sur un vrai fichier). Le defaut ("AMS","LP-X") ci-dessous
# (get avec fallback) couvrait deja N0 par coincidence ; explicite ici
# pour ne plus en dependre.
_ANISO_MAGIC_INFO = {
    "A0": ("ATRM", "LP-AN-TRM"),
    "F0": ("AARM", "LP-AN-ARM"),
    "N0": ("AMS", "LP-X"),
}


# aniso_s/aniso_tilt_correction/aniso_ftest* pour specimens.txt, depuis un
# AniTensor DEJA calcule (.pmagani, code2='A0') - demande explicite
# utilisateur ("il manque aussi les données d'anisotropie au niveau du
# fichier specimens (prendre les A0)"), MEME manque que celui deja corrige
# cote AMS_Py.exporter_magic_dialog (voir son docstring pour le detail des
# colonnes verifiees contre le data model reel) - ici pour le tenseur natif
# 'A0' de STARpaleomag_Py (calcul.AniTensor), pas le format AMS_Py.
# aniso_tilt_correction="-1" (coordonnees specimen) : aniso_s est exporte
# BRUT, sans reorientation - meme raisonnement/meme valeur que cote AMS_Py
# (deja l'hypothese par defaut de pmagpy.ipmag quand la colonne est
# absente).
def anisotropy_specimen_magic_fields(tensor: AniTensor) -> Dict[str, str]:
    trace = tensor.k11 + tensor.k22 + tensor.k33
    if not trace:
        return {}
    aniso_type, _method_codes = _ANISO_MAGIC_INFO.get(tensor.code2, ("AMS", "LP-X"))
    s1, s2, s3 = tensor.k11 / trace, tensor.k22 / trace, tensor.k33 / trace
    s4, s5, s6 = tensor.k12 / trace, tensor.k23 / trace, tensor.k13 / trace
    out = {
        "aniso_type": aniso_type,
        "aniso_s": f"{s1:.6f}:{s2:.6f}:{s3:.6f}:{s4:.6f}:{s5:.6f}:{s6:.6f}",
        "aniso_tilt_correction": "-1",
        "description": f"{aniso_type} tensor from STARpaleomag_Py, specimen coordinates",
    }
    n_pos = tensor.n_positions
    if n_pos is None and tensor.code2 == "A0":
        n_pos = 6  # ATRM natif = methode a 6 positions (detect_six_positions)
    if n_pos is not None:
        out["aniso_s_n_measurements"] = str(n_pos)
    if tensor.sigma is not None:
        out["aniso_s_sigma"] = f"{tensor.sigma:.6g}"
    if tensor.ftest is not None:
        out["aniso_ftest"] = f"{tensor.ftest:.6g}"
    if tensor.ftest12 is not None:
        out["aniso_ftest12"] = f"{tensor.ftest12:.6g}"
    if tensor.ftest23 is not None:
        out["aniso_ftest23"] = f"{tensor.ftest23:.6g}"
    if tensor.quality in ("g", "b"):
        out["aniso_ftest_quality"] = tensor.quality
    return out


# aniso_v1/v2/v3/aniso_p.../sites.txt (groupe "Anisotropy") depuis une
# moyenne de site DEJA calculee (calcul.AniMeanTensor, .pmagani section
# "#site mean tensor results") - demande explicite utilisateur ("je
# voudrais ajouter l'exportation de l'AMS (niveau specimen et sites).
# Les donnees d'AMS vont aussi dans les fichiers sites.txt et
# specimens.txt" - le niveau specimen existait deja, voir
# anisotropy_specimen_magic_fields ci-dessus). Port DIRECT de
# AMS_Py.ams_stats.magic_site_aniso_fields/shape_params (meme app,
# meme probleme deja resolu et verifie contre MagIC-data-model.txt
# reel - y compris le desaccord documente avec un bug repere dans
# pmagpy.ipmag ou aniso_v3 reutilise a tort l'angle e12 de v1 au lieu
# de e13) : PAS reimplemente a la main, seulement adapte aux champs
# deja disponibles sur AniMeanTensor (k1/k2/k3/dec/inc/alpha1_N/alpha2_N
# et P/T/L/F/Pprim DEJA fournis par l'appelant, contrairement a
# TensorialMeanResult ou tsmean/shape_params les calcule) - evite de
# redemontrer des formules deja etablies. Contrairement a AMS_Py (qui
# connait l'orientation Sa/IS/TC du resultat moyenne, STARpaleomag_Py
# n'a pas de champ AniMeanTensor dedie - mais AMS_Py l'ecrit deja en
# texte libre dans `info` (ex. 'tilt_correction: 0; site mean (AMS_Py)')
# a l'archivage (voir AMS_Py.ams_selection._mean_row_to_result, meme
# convention de lecture en sens inverse) : aniso_tilt_correction est
# RETROUVE depuis ce texte (voir _tilt_correction_from_info) plutot que
# fige a "-1" - BUG REEL corrige ici (demande explicite utilisateur :
# "la moyenne est en in situ et/ou TC rarement en SC. Il manque cette
# info dans le fichier pmagani" - l'info N'ETAIT PAS manquante dans le
# fichier, juste jamais lue cote export STARpaleomag_Py, qui annoncait
# systematiquement "-1"/coordonnees specimen meme pour une moyenne
# reellement IS ou TC). "-1" reste le defaut si `info` ne porte aucune
# mention exploitable (moyenne d'une origine plus ancienne/differente).
# Format "tau:dec:inc:eta/zeta:..." SANS espaces autour des ":"
# (contrairement a AMS_Py, " : ") - aligne sur pmagpy.convert_2_magic.py
# (source faisant foi, ':'.join(...)) et sur aniso_s ci-dessus, deja
# sans espaces dans ce module."""
_INFO_TILT_CORRECTION_RE = re.compile(r"tilt_correction:\s*(-?\d+)\b")
# code STOCKE dans `info` (convention AMS_Py, voir _FILE_CODE_TO_ORIENT
# cote AMS_Py : 1=Sa/2=IS/3=TC en interne, mais le TEXTE ecrit est deja
# le code fichier "1"/"0"/"100", PAS l'enum 1/2/3) -> code MagIC reel
# (voir AMS_Py.ams_stats.MAGIC_TILT_CORRECTION_CODE) : seul "1"
# (specimen/Sa) differe ("-1" en MagIC) ; "0"(IS) et "100"(TC) sont deja
# les codes MagIC tels quels.
_INFO_CODE_TO_MAGIC_TILT = {"1": "-1", "0": "0", "100": "100"}


def _tilt_correction_from_info(info: str) -> str:
    m = _INFO_TILT_CORRECTION_RE.search(info or "")
    if not m:
        return "-1"
    return _INFO_CODE_TO_MAGIC_TILT.get(m.group(1), "-1")


def _mean_tilt_correction(mean: AniMeanTensor) -> str:
    """`mean.tilt_correction` (colonne dediee, .pmagani ecrit apres
    l'ajout de ce champ - voir AniMeanTensor.__doc__) en priorite ; repli
    sur le texte libre `info` (`_tilt_correction_from_info`) pour un
    fichier plus ancien qui ne l'a jamais eue - demande explicite
    utilisateur ("oui ajouter une colonne avant info")."""
    if mean.tilt_correction:
        return _INFO_CODE_TO_MAGIC_TILT.get(mean.tilt_correction.strip(), "-1")
    return _tilt_correction_from_info(mean.info)


def anisotropy_site_magic_fields(mean: AniMeanTensor) -> Dict[str, str]:
    k1, k2, k3 = mean.k1, mean.k2, mean.k3
    aniso_type, _method_codes = _ANISO_MAGIC_INFO.get(mean.code2, ("AMS", "LP-X"))

    # aniso_perc/aniso_total/aniso_ll/aniso_ff/aniso_vg/aniso_fl : formules
    # officielles (MagIC-data-model.txt) - PAS deja sur AniMeanTensor
    # (seuls P/T/L/F/Pprim - Jelinek - y sont stockes), memes formules que
    # AMS_Py.ams_stats.magic_site_aniso_fields.
    total_k = k1 + k2 + k3
    aniso_perc = 100.0 * (k1 - k3) / total_k if total_k else 0.0
    mean_k = total_k / 3.0 if total_k else 0.0
    aniso_total = 100.0 * (k1 - k3) / mean_k if mean_k else 0.0
    aniso_ll = math.log(k1 / k2) if k2 else 0.0
    aniso_ff = math.log(k2 / k3) if k3 else 0.0
    aniso_vg = math.degrees(math.asin(math.sqrt((k2 - k3) / (k1 - k3)))) if (k1 - k3) else 0.0
    aniso_fl = (mean.F / mean.L) if (mean.L and mean.F is not None) else 0.0

    dec = (mean.dec1, mean.dec2, mean.dec3)
    inc = (mean.inc1, mean.inc2, mean.inc3)
    k = (mean.k1, mean.k2, mean.k3)
    # alpha_major(axis N)=alpha1_N (demi-grand axe, vers l'axe j - "eta"),
    # alpha_minor(axis N)=alpha2_N (demi-petit axe, vers l'axe k - "zeta") -
    # voir AniMeanTensor.__doc__ ; eta/zeta pointent vers les DEUX AUTRES
    # axes propres eux-memes (comme pmagpy.ipmag), pas une orientation
    # d'ellipse recalculee - meme simplification qu'AMS_Py, meme ecart
    # documente avec le formalisme par paire e12/e13/e23 de Hext (1963).
    alpha_major = (mean.alpha1_1, mean.alpha1_2, mean.alpha1_3)
    alpha_minor = (mean.alpha2_1, mean.alpha2_2, mean.alpha2_3)

    def aniso_v(i: int, j: int, kk: int) -> str:
        return ":".join(str(v) for v in [
            k[i], dec[i], inc[i], "eta/zeta",
            dec[j], inc[j], alpha_major[i],
            dec[kk], inc[kk], alpha_minor[i],
        ])

    fields = {
        "aniso_type": aniso_type,
        "aniso_tilt_correction": _mean_tilt_correction(mean),
        "aniso_v1": aniso_v(0, 1, 2),
        "aniso_v2": aniso_v(1, 0, 2),
        "aniso_v3": aniso_v(2, 0, 1),
        "aniso_perc": f"{aniso_perc:.4f}",
        "aniso_total": f"{aniso_total:.4f}",
        "aniso_ll": f"{aniso_ll:.6f}",
        "aniso_ff": f"{aniso_ff:.6f}",
        "aniso_vg": f"{aniso_vg:.4f}",
        "aniso_fl": f"{aniso_fl:.6f}",
    }
    if mean.P is not None:
        fields["aniso_p"] = f"{mean.P:.6f}"
    if mean.Pprim is not None:
        fields["aniso_pp"] = f"{mean.Pprim:.6f}"
    if mean.T is not None:
        fields["aniso_t"] = f"{mean.T:.6f}"
    if mean.L is not None:
        fields["aniso_l"] = f"{mean.L:.6f}"
    if mean.F is not None:
        fields["aniso_f"] = f"{mean.F:.6f}"
    # Meme colonne `description` que AMS_Py.ams_stats.magic_site_aniso_fields
    # (verifiee manquante ici, voir _SITES_HEADER) : documente dans le
    # fichier MagIC lui-meme (pas seulement dans ce docstring) la
    # simplification eta/zeta partagee par pmagpy.ipmag ET AMS_Py -
    # demande explicite utilisateur ("porter cette colonne description
    # manquante vers STARpaleomag_Py").
    fields["description"] = (
        "eta/zeta point toward the other two eigenvectors with Jelinek (1978) "
        "confidence semi-angles (alpha1_N/alpha2_N from the .pmagani mean "
        "tensor), not Hext (1963) pairwise e12/e13/e23 - see "
        "anisotropy_site_magic_fields docstring"
    )
    return fields

# Colonnes MagIC v3 (table "specimens", groupe paleointensite) construites
# depuis UNE ligne .pmagint (voir paleointensity.read_pmagint) - demande
# explicite utilisateur ("in export to Magic, it is not asking for
# Paleointensity") : export_to_magic ignorait entierement .pmagint jusqu'ici
# (voir le docstring de module - "Paleointensite... HORS PERIMETRE pour
# l'instant (a ajouter plus tard si besoin)"). Chaque colonne VERIFIEE
# contre le data model reel (data_model.json, table "specimens") plutot que
# devinee - notamment : int_abs/int_abs_sigma/int_treat_dc_field sont en
# TESLA (T) alors que .pmagint est en microteslas (Hlab/H/Hcorani/HcorCool)
# - conversion *1e-6 ici ; int_mad_free VERIFIE venir de pars["int_mad_free"]
# (voir paleointensity_magic.MagicPintResult.mad) - PAS une supposition ;
# int_gmax VERIFIE etre le nom reel de la colonne "gap_max" (grep direct du
# source pmagpy.pmag, gmax_key = 'int_gmax').
#
# int_abs (le resultat final) : precedence HcorCool > Hcorani > H, MEME
# ordre que `h_final` cote app.ouvrir_openfilepint_dialog (cooling puis
# anisotropie puis brut) - PAS Hcorani seul, qui ignorerait une correction
# de refroidissement deja appliquee.
#
# `ccr` (correlation Arai) n'est PAS mappe vers int_r2_corr : la colonne
# reelle est un COEFFICIENT DE DETERMINATION (r au carre, "r2"), alors que
# rien ne confirme si `ccr` cote STARpaleomag_Py est deja au carre ou non -
# risque de mismatch d'unite/definition pour un champ secondaire, laisse de
# cote plutot que devine.
def paleointensity_magic_fields(row: Dict[str, str]) -> Dict[str, str]:
    def f(key: str) -> Optional[float]:
        v = row.get(key, "n.d")
        if v in ("n.d", ""):
            return None
        try:
            return float(v)
        except ValueError:
            return None

    h_final = f("HcorCool")
    if h_final is None:
        h_final = f("Hcorani")
    if h_final is None:
        h_final = f("H")
    hlab = f("Hlab")
    b, sb = f("b"), f("sb")
    fcor, fcor_cool = f("fcor"), f("fcorCool")

    out: Dict[str, str] = {}
    # meas_step_min/max/unit MANQUAIENT ici (bug signale par l'utilisateur -
    # "le meas_temp_min et le meas_tem_max ne sont pas entres dans le
    # fichier specimens.txt") : deja renseignes pour un resultat
    # DIRECTIONNEL (_dir_rows_for_specimen, ci-dessus) mais jamais pour un
    # resultat de PALEOINTENSITE, seul cas ou `paleointensity_magic_fields`
    # est utilisee - `t1`/`t2` (colonnes .pmagint, voir _PMAGINT_HEADER)
    # sont la temperature du premier/dernier pas retenu dans l'intervalle
    # de fit (write_pmagint_line: `points[n1-1].temp`/`points[n2-1].temp`),
    # en degC comme `etape`/`FitResult.step_first/last` - meme conversion
    # +273/"K" que le cas thermique de `_dir_rows_for_specimen`.
    t1, t2 = f("t1"), f("t2")
    if t1 is not None and t2 is not None:
        out["meas_step_min"] = _num(min(t1, t2) + 273, 4)
        out["meas_step_max"] = _num(max(t1, t2) + 273, 4)
        out["meas_step_unit"] = "K"
    if h_final is not None:
        out["int_abs"] = _num_sci(h_final * 1.0e-6)
    if hlab is not None and b is not None and sb is not None:
        out["int_abs_sigma"] = _num_sci(abs(sb * hlab) * 1.0e-6)
    if hlab is not None:
        out["int_treat_dc_field"] = _num_sci(hlab * 1.0e-6)
    out["int_corr"] = "c" if (fcor is not None or fcor_cool is not None) else "u"
    if fcor is not None:
        out["int_corr_aniso"] = _num(fcor, 3)
    if fcor_cool is not None:
        out["int_corr_cooling_rate"] = _num(fcor_cool, 3)

    simple_map = [
        ("N", "int_n_measurements", 0), ("f", "int_f", 3), ("g", "int_g", 3),
        ("q", "int_q", 1), ("mad", "int_mad_free", 1), ("dang", "int_dang", 1),
        ("b", "int_b", 4), ("sb", "int_b_sigma", 4), ("sb_over_b", "int_b_beta", 4),
        ("k", "int_k", 4), ("k_sse", "int_k_sse", 5),
        ("fvds", "int_fvds", 3), ("frac", "int_frac", 3), ("gap_max", "int_gmax", 3),
        ("n_ptrm", "int_n_ptrm", 0), ("gamma", "int_gamma", 1),
    ]
    for src_key, dst_key, decimals in simple_map:
        v = f(src_key)
        if v is not None:
            out[dst_key] = str(int(v)) if decimals == 0 else _num(v, decimals)
    return out

_METHOD_CODE_BASE = {
    "A": "LT-AF-I", "I": "LT-IRM", "F": "LT-AF-Z", "N": "LT-NO",
    "D": "LT-T-Z", "S": "LT-T-Z", "R": "LT-T-I", "V": "LT-T-I",
    "P": "LT-PTRM-I", "X": "LT-T-I", "Y": "LT-T-I", "Z": "LT-T-I",
    "L": "LT-T-I", "Q": "LT-T-I",
}


def _specimen_method_codes(ech: SelectedSample) -> str:
    """Resume (pas dans une citation exacte du Fortran - la portion du
    rapport d'exploration couvrant cette boucle specifique etait
    paraphrasee) : union des codes de base rencontres parmi les mesures de
    l'echantillon, plus `:LP-CR-TRM` des qu'AU MOINS une mesure est un pas
    de vitesse de refroidissement (cod1 'L' ou 'Q') - CORRIGE par rapport
    au Fortran, qui ne regardait (bug de portee de variable) que la
    DERNIERE mesure de l'echantillon."""
    codes = []
    cooling = False
    for m in ech.mesures:
        code = _METHOD_CODE_BASE.get(m.cod1)
        if code and code not in codes:
            codes.append(code)
        if m.cod1 in ("L", "Q"):
            cooling = True
    if cooling:
        codes.append("LP-CR-TRM")
    return ":".join(codes) if codes else "LT-NO"


def _dir_rows_for_specimen(ech: SelectedSample, results: List[FitResult]) -> List[Dict[str, str]]:
    """Equivalent de `magicdirres` : pour chaque FitResult de cet
    echantillon, jusqu'a 3 lignes (une par etat de correction
    d'orientation - repere echantillon / geographique / tectonique, la
    tectonique etant sautee si dip==0)."""
    rows = []
    for r in [res for res in results if res.id.strip() == ech.id.strip()]:
        incr, decr = math.radians(r.inc), math.radians(r.dec)
        x = math.cos(incr) * math.cos(decr)
        y = math.cos(incr) * math.sin(decr)
        z = math.sin(incr)

        if r.demag == "F":
            step_min = r.step_first * 1.0e-4
            step_max = r.step_last * 1.0e-4
            step_unit = "T"
        else:
            step_min = r.step_first + 273
            step_max = r.step_last + 273
            step_unit = "K"

        if r.cat1 == "P":
            comp = "Plane"
        elif r.cat1 == "f":
            comp = "Fisher"
        elif r.cat1 == "s":
            comp = "Blanket"
        else:
            comp = _comp_name(r.component)

        base_codes = []
        if r.cat1 == "L" and r.orig == "o":
            base_codes.append("DE-BFL-A")
        elif r.cat1 == "L":
            base_codes.append("DE-BFL")
        elif r.cat1 == "P":
            base_codes.append("DE-BFP")
        elif r.cat1 == "f":
            base_codes.append("DE-FM")
        elif r.cat1 == "s":
            base_codes.append("DE-BLANKET")

        states = [(1, "-1", x, y, z)]
        xx, yy, zz = corfor(x, y, z, r.cin, r.caz)
        states.append((2, "0", xx, yy, zz))
        if r.dip != 0.0:
            xx3, yy3, zz3 = corpen(xx, yy, zz, r.dip, r.str_)
            states.append((3, "100", xx3, yy3, zz3))

        mad_anc = _num(r.mad, 1) if r.orig == "o" else ""
        mad_free = _num(r.mad, 1) if r.orig != "o" else ""

        for _iori, tilt, px, py, pz in states:
            _mag, dec, inc = polere(px, py, pz)
            rows.append({
                "dir_tilt_correction": tilt,
                "result_quality": "g",
                "meas_step_min": _num(step_min, 4),
                "meas_step_max": _num(step_max, 4),
                "meas_step_unit": step_unit,
                "dir_comp": comp,
                "dir_dec": _num(dec, 1),
                "dir_inc": _num(inc, 1),
                "dir_n_measurements": str(r.nb),
                "dir_mad_anc": mad_anc,
                "dir_mad_free": mad_free,
                "_method_codes_extra": ":".join(base_codes),
            })
    return rows


def build_specimens_rows(
    samples: List[SelectedSample], results: List[FitResult],
    pmagint_rows: Optional[Dict[str, Dict[str, str]]] = None,
    aniso_tensors: Optional[Dict[str, AniTensor]] = None,
) -> List[List[str]]:
    """`pmagint_rows` (specimen -> ligne .pmagint brute, voir
    paleointensity.read_pmagint) et `aniso_tensors` (specimen -> AniTensor
    'A0' deja calcule, voir calcul.read_ani_tensor) : quand presents pour
    un specimen, leurs champs (paleointensity_magic_fields/
    anisotropy_specimen_magic_fields) sont fusionnes sur la PREMIERE ligne
    de ce specimen SEULEMENT (bare row si aucun resultat directionnel,
    sinon la 1ere ligne directionnelle - coordonnees echantillon) - PAS
    repete sur les lignes IS/tectonique du meme specimen (ni le resultat
    de paleointensite ni le tenseur d'anisotropie ne dependent de l'etat
    de correction d'orientation directionnelle, le dupliquer serait juste
    du bruit) - demande explicite utilisateur ("in export to Magic, it is
    not asking for Paleointensity" puis "il manque aussi les données
    d'anisotropie au niveau du fichier specimens (prendre les A0)")."""
    rows = []
    for ech in samples:
        # Porteur de metadonnees de site sans vraie mesure - voir la meme
        # garde dans build_samples_rows (un nom de specimen vide est
        # invalide pour specimens.txt).
        if not ech.id.strip() and not ech.mesures:
            continue
        azimuth = ech.caz - 90.0
        if azimuth < 0.0:
            azimuth += 360.0
        if ech.norme == "v":
            volume, weight = _num_sci(ech.vol * 1.0e-6), ""
        else:
            volume, weight = "", _num_sci(ech.vol * 1.0e-3)

        base_method_codes = _specimen_method_codes(ech)
        sample_name = ech.magic_sample.strip() or ech.id

        base = {
            "specimen": ech.id, "citations": "This study", "sample": sample_name,
            "geologic_classes": ech.magic_gc, "lithologies": ech.magic_li,
            "geologic_types": ech.magic_smt, "azimuth": _num(azimuth, 1),
            "dip": _num(-ech.cin, 1), "volume": volume, "weight": weight,
        }

        extra_fields = {}
        extra_codes = []
        pint_row = (pmagint_rows or {}).get(ech.id)
        if pint_row is not None:
            extra_fields.update(paleointensity_magic_fields(pint_row))
            extra_codes.append("LP-PI-TRM")
        tensor = (aniso_tensors or {}).get(ech.id)
        if tensor is not None:
            aniso_fields = anisotropy_specimen_magic_fields(tensor)
            if aniso_fields:
                extra_fields.update(aniso_fields)
                _aniso_type, aniso_method = _ANISO_MAGIC_INFO.get(tensor.code2, ("AMS", "LP-X"))
                extra_codes.append(aniso_method)

        dir_rows = _dir_rows_for_specimen(ech, results)
        if not dir_rows:
            method_codes = base_method_codes
            for code in extra_codes:
                if code not in method_codes.split(":"):
                    method_codes += f":{code}"
            row = dict(base, method_codes=method_codes, **extra_fields)
            rows.append([row.get(col, "") for col in _SPECIMENS_HEADER])
            continue

        for i, dr in enumerate(dir_rows):
            extra = dr.pop("_method_codes_extra", "")
            method_codes = base_method_codes + (f":{extra}" if extra else "")
            row_fields = extra_fields if i == 0 else {}
            if row_fields:
                for code in extra_codes:
                    if code not in method_codes.split(":"):
                        method_codes += f":{code}"
            row = dict(base, method_codes=method_codes, **dr, **row_fields)
            rows.append([row.get(col, "") for col in _SPECIMENS_HEADER])
    return rows


# ---------------------------------------------------------------------------
# measurements.txt
# ---------------------------------------------------------------------------

_MEASUREMENTS_HEADER = [
    "citations", "analysts", "specimen", "experiment", "software_packages",
    "timestamp", "measurement", "quality", "standard", "treat_step_num",
    "treat_temp", "treat_ac_field", "treat_dc_field", "treat_dc_field_phi",
    "treat_dc_field_theta", "meas_temp", "dir_inc", "dir_dec",
    "magn_moment", "magn_volume", "magn_mass", "dir_csd", "susc_chi_volume",
    "susc_chi_mass", "method_codes", "instrument_codes", "description",
]

_INSTRUMENT_CODES = {
    "": ("2 positions", ""),
    "C1": ("2G magnetometer one position", "1"),
    "C4": ("2G magnetometer 4 positions", "4"),
    "JA": ("JR5 or JR6 automatic", "3"),
    "J2": ("JR5 or JR6 2 positions", "2"),
    # "S" : Schonstedt spinner (voir convert_ren_to_r.py, meme code deja
    # traite a part pour le champ "error", jeux de donnees de ~30 ans) -
    # demande explicite utilisateur ("Si S mettre Schonstedt").
    "S": ("Schonstedt", ""),
    # "J"/"J?" : modele Agico Spinner non precise dans les vieux fichiers
    # Rennes ; "J5"/"J6" : modele precise (JR5 vs JR6) - demande explicite
    # utilisateur ("si le code magnetometre est J, ou J? mettre Agico
    # Spinner si c'est J5 ou J6 : Agico Spinner Jr5 ou Agico Spinner
    # Jr6"), distinct de "JA"/"J2" ci-dessus (qui denotent le PROTOCOLE -
    # automatique/2 positions - pas le modele d'instrument).
    "J": ("Agico Spinner", ""),
    "J?": ("Agico Spinner", ""),
    "J5": ("Agico Spinner Jr5", ""),
    "J6": ("Agico Spinner Jr6", ""),
}


def _instrument(ins: str) -> str:
    ins = (ins or "").strip()
    if ins in _INSTRUMENT_CODES:
        return _INSTRUMENT_CODES[ins][0]
    if ins[:1] == "M":
        return "Molspin spinner"
    return "2 positions"


# Anisotropie (cod1 X/Y/Z, 6+ positions +/-) - demande explicite
# utilisateur ("anisotropy with code X+,Y+,Z+ X-,Y-,Z- is usually done by
# TRM acquisition and is found in files with paleointensity experiments,
# but there are some done with high field IRM, in that case, the strong
# field dc field is the step value... ask the user to confirm the kind of
# anisotropy experiment when it is dubious? and whether it wants to
# archive these data") : deux protocoles distincts partagent le meme
# cod1/cod2 - TRM (chauffe + champ labo, `etape` = temperature en degC,
# LP-AN-TRM) ou IRM fort champ (`etape` = champ fort en mT, LP-AN-IRM),
# indiscernables sans contexte. Auto-classifie quand le contexte est
# fiable, sinon laisse a l'appelant (app.py, seul endroit avec un moyen
# d'interroger l'utilisateur) le soin de demander - voir
# ouvrir_export_magic_dialog.
_ANISOTROPY_AXIS_COD1 = ("X", "Y", "Z")
_PALEOINT_COMPANION_COD1 = {"R", "V", "P", "S"}
_MAX_PLAUSIBLE_TEMP_C = 700.0  # au-dela, ne peut plus etre une temperature de chauffe usuelle


def classify_anisotropy_experiment(mesures: List[Measurement]) -> Optional[str]:
    """'trm' ou 'irm' si le contexte permet une classification fiable,
    None si DOUTEUX (l'appelant doit demander a l'utilisateur) - pour
    l'ENSEMBLE des pas d'anisotropie (cod1 X/Y/Z) d'un specimen :
    - 'trm' des qu'un pas de paleointensite genuine (R/V/P/S, un vrai
      protocole Thellier/IZZI sur CE specimen) est present - "found in
      files with paleointensity experiments".
    - 'irm' des qu'un `etape` d'un pas X/Y/Z depasse
      _MAX_PLAUSIBLE_TEMP_C (ne peut physiquement pas etre une
      temperature de chauffe - "the strong field dc field is the step
      value", un champ fort typique en mT depasse largement 700).
    - None (douteux) sinon - ex. `etape` dans une plage plausible pour
      les DEUX interpretations (temperature ordinaire OU champ IRM
      modere) sans compagnon paleointensite pour trancher. Verifie sur
      donnees reelles (Tibet_14_15_Pmag.txt) : plusieurs specimens
      genuinement ambigus (etape=520/510/150 sans R/V/P/S, ou etape=1100
      sans aucune ambiguite - largement au-dela de 700)."""
    axis_steps = [m for m in mesures if m.cod1 in _ANISOTROPY_AXIS_COD1]
    if not axis_steps:
        return None
    if any(m.cod1 in _PALEOINT_COMPANION_COD1 for m in mesures):
        return "trm"
    if any(m.etape > _MAX_PLAUSIBLE_TEMP_C for m in axis_steps):
        return "irm"
    return None


def _measurement_treatment(
    m: Measurement, prev: List[Measurement], ifield: float, anisotropy_kind: str = "trm",
    is_paleointensity: bool = False, is_pure_ii_protocol: bool = False,
) -> Tuple[str, float, float, float, float, float]:
    """Equivalent du `select case (mes(j).cod1)` de export2magic (voir §4
    du rapport d'exploration) : retourne (method_codes, treat_temp,
    treat_ac_field, treat_dc_field, treat_dc_field_phi, treat_dc_field_theta).
    `prev` = mesures precedentes DE CE MEME echantillon, dans l'ordre - sert
    au controle pTRM ('P'), avec le BON index (contrairement au bug du
    Fortran qui relisait une variable de boucle perimee).

    BUG CORRIGE ici (repere en construisant paleointensity_magic.py, PAS
    dans le Fortran - propre a ce portage) : R/V/P/X/Y/Z/L/Q ne
    renseignaient PAS `temp` (restait a 0.0, alors que ce sont tous des
    pas en temperature - meme ensemble que `convert_ren_to_r._TEMP_CODES`).
    Cela passait inapercu car aucun appelant existant (convert_ren_to_r.py,
    magic_export.build_measurements_rows) n'utilisait ce `temp` retourne
    par CETTE fonction pour ces cod1 - ils le recalculent independamment
    depuis `etape`. paleointensity_magic.py est le premier a en avoir
    besoin directement (pmag.sortarai exige un treat_temp correct sur
    CHAQUE pas, y compris les pas en champ R/V).

    R et V restent TOUS LES DEUX des pas EN CHAMP (theta=+90/-90), meme
    en Thellier sans 'S' - PAS de "role invers" selon le protocole (une
    version anterieure de ce fix l'avait suppose a tort, en lisant trop
    litteralement `paleointensity.detect_method_and_hlab`). Confirme par
    l'utilisateur sur des donnees reelles (specimen "02B") : la methode
    Thellier utilisee ici applique le champ dans DEUX orientations
    opposees (Z puis Z-, "there is no choice than (R+V)/2") - il n'existe
    PAS de pas reellement zero-field distinct de R ou V ; seule LA MOYENNE
    (R+V)/2 annule le pTRM et redonne le NRM restant. Ce fichier doit
    representer la mesure BRUTE (R/V sont reellement mesures en champ) -
    la reconstruction (R+V)/2 necessaire pour nourrir PmagPy est donc
    faite cote analyse, dans paleointensity_magic.py, PAS ici."""
    temp = af_field = dc_field = phi = theta = 0.0
    codes = "LT-NO"
    etape = m.etape

    if m.cod1 == "A":
        # `etape` = mT REEL directement (plus l'equivalent Oersted mT*10
        # d'avant - demande explicite utilisateur "convert all step
        # integer to float") : *0.001 (mT -> Tesla) plutot que l'ancien
        # *0.0001 qui compensait cette echelle *10 disparue.
        af_field = etape * 0.001
        dc_field = ifield * 1.0e-6
        codes = "LT-AF-I"
    elif m.cod1 == "I":
        dc_field = etape * 0.001
        codes = "LT-IRM"
    elif m.cod1 == "F":
        # Demande explicite utilisateur, en plusieurs temps :
        # 1) "during import F= just put LT-AF-Z ; we do not know the
        #    detail. FT add the code for AF tumbler".
        # 2) "all AF demagnetization were done [with] the 3 axis
        #    degausser if the instrument is C as the degausser is online
        #    with the magnetometer ; best to add the complement of the
        #    Magic code as defined before, especially for the FX,FY,FZ.
        #    For the other instrument, the AF degausser was a tumbler." -
        #    verifie sur donnees reelles (old_pmag.ren) : cod2 en 'X'/'Y'/
        #    'Z' pour cod1='F' n'apparait QUE pour l'instrument "C1"
        #    (32/32/46 occurrences), jamais pour "Mo"/"J2"/"S".
        # 3) CORRECTION IMPORTANTE (l'utilisateur les a lui-meme
        #    introduits dans le vocabulaire MagIC officiel - "it is
        #    official MagIC vocabulary? I introduce it to Magic!") :
        #    "LT-AF-Z-X"/"-Y"/"-Z" SONT des codes MagIC reels, verifie en
        #    RE-INTERROGEANT le vocabulaire EN LIGNE
        #    (https://www2.earthref.org/MagIC/method-codes.json, meme
        #    source que controlled_vocabularies3.py) plutot que le fichier
        #    LOCAL/PERIME embarque avec pmagpy 4.5.2
        #    (pmagpy/data_model/method_codes.json, qui ne les contient
        #    pas) - la premiere verification etait donc fausse (base sur
        #    une copie obsolete). Le vocabulaire en ligne contient aussi
        #    "LT-AF-Z-XZY"/"-YZX"/"-XYZ"/"-YXZ"/"-ZXY"/"-ZYX" (sequence
        #    complete des 3 axes, un code par ordre) - confirme par
        #    l'utilisateur : sur l'instrument "C1", cod2='+' correspond a
        #    la sequence Y,Z,X (LT-AF-Z-YZX) et cod2='-' a X,Z,Y
        #    (LT-AF-Z-XZY) - la tres large majorite des pas AF au C1 (cod2
        #    '+'/'-', ~4750 sur ~4780) portent donc CETTE information,
        #    X/Y/Z (un seul axe) restant l'exception rare.
        # etape = mT reel directement, voir commentaire du cas cod1=='A'
        # ci-dessus (*0.001 mT->Tesla, plus l'ancien *0.0001 Oersted-scale).
        af_field = etape * 0.001
        ins = (m.ins or "").strip().upper()
        if etape == 0:
            codes = "LT-NO"
        elif ins.startswith("C"):
            codes = {
                "X": "LT-AF-Z:LT-AF-Z-X", "Y": "LT-AF-Z:LT-AF-Z-Y", "Z": "LT-AF-Z:LT-AF-Z-Z",
                "+": "LT-AF-Z:LT-AF-Z-YZX", "-": "LT-AF-Z:LT-AF-Z-XZY",
            }.get(m.cod2, "LT-AF-Z")
        else:
            # Tout instrument hors "C*" : degausser tumbler (mecanique,
            # hors ligne) pour TOUT pas AF, quel que soit cod2 - un
            # tumbler n'a pas de notion d'ordre X/Y/Z discret.
            codes = "LT-AF-Z:LT-AF-Z-TUMB"
    elif m.cod1 == "N":
        temp = etape + 273
        # BUG CORRIGE ici (crash reel confirme : thellier_gui de PmagPy
        # plante avec "IndexError: list index out of range" dans get_data,
        # a la ligne `NRM = zijdblock[0][3]`, en ouvrant un export MagIC
        # produit par ce fichier). Cause: pour data_model==3, thellier_gui
        # filtre TOUTES les lignes de measurements.txt qui ne contiennent
        # aucun de LP-PI-TRM/LP-TRM/LP-PI-M/LP-AN/LP-CR-TRM AVANT de trier
        # les blocs par specimen - une ligne LT-NO seule (sans LP-PI-TRM)
        # est donc supprimee entierement, y compris pour son role de pas
        # zero-field initial (zijdblock), meme si sortarai lui-meme
        # n'exige aucun suffixe ZI/IZ/BT-IZZI particulier sur ce pas.
        # Reproduit et verifie directement en executant get_data() de
        # thellier_gui.py (tag v4.5.0 de PmagPy/PmagPy) sur un export reel
        # (specimen 15SO131604C) : son pas N s'exportait en "LT-NO" nu,
        # contrairement a un fichier MagIC de reference fonctionnel ou le
        # pas NRM porte toujours "LT-NO:LP-PI-TRM:...". Corrige seulement
        # pour les specimens ou un pas R/V/P/S (paleointensite genuine) est
        # present ailleurs - ne pas re-etiqueter un NRM d'un specimen
        # purement directionnel.
        codes = "LT-NO:LP-PI-TRM" if is_paleointensity else "LT-NO"
        if m.cod2 == "P" or (is_paleointensity and is_pure_ii_protocol):
            # Specimen R/V/P sans aucun pas 'S' (pas de vrai zero-field -
            # voir le commentaire au cas cod1=='R' ci-dessous : "R et V
            # sont TOUS LES DEUX en-champ ... seule LA MOYENNE (R+V)/2
            # annule le pTRM"). BUG CORRIGE ici (deuxieme, distinct du
            # "LT-NO nu" ci-dessus) : sans ce tag "LP-PI-II" sur le pas N,
            # thellier_gui de PmagPy (get_data) ne construit JAMAIS de
            # vrai Zijderveld multi-points pour ce protocole - son
            # scan ordinaire ne retient que les pas LT-NO/LT-T-Z/LT-M-Z/
            # LT-AF-Z (aucun ici, R/V sont tagues LT-T-I), et son
            # branchement special de reconstruction via araiblock (qui,
            # lui, gere bien le cas R/V par soustraction vectorielle)
            # n'est active QUE si `Data[s]['datablock'][0]` (le premier
            # pas du specimen ayant LP-PI-TRM/LP-PI-M, donc le pas N
            # lui-meme depuis le fix precedent) porte LP-PI-II/
            # LP-PI-M-II/LP-PI-T-II - jamais le cas jusqu'ici puisque
            # seul le pas V (jamais le premier) le portait. Verifie par
            # reproduction directe (specimen 15SO131604C, meme
            # get_data() extrait du tag v4.5.0 de PmagPy) : Zijderveld a
            # 1 seul point (le NRM) sans ce tag, correctement reconstruit
            # (1 point par temperature) avec.
            codes = "LT-NO:LP-PI-TRM:LP-PI-II:LP-PI-ALT"
    elif m.cod1 == "D":
        temp = etape + 273
        codes = "LT-NO" if etape < 21.0 else "LT-T-Z"
    elif m.cod1 == "S":
        temp = etape + 273
        codes = "LT-T-Z:LP-PI-TRM-IZ:LP-PI-TRM:LP-PI-ALT-PTRM:LP-PI-BT-IZZI"
    elif m.cod1 == "R":
        # REVERT (voir git history de cette session) : une version
        # anterieure de ce fix rendait R zero-field en THELLIER (izzi=
        # False), en lisant trop litteralement le commentaire de
        # detect_method_and_hlab ("R=zero-field"). Donnees reelles
        # (utilisateur, specimen "02B") : ni R ni V n'est un vrai pas
        # zero-field ici - c'est la methode Thellier "champ en Z puis Z-"
        # (confirmee par l'utilisateur : "field in Z and Z-, there is no
        # choice than (R+V)/2") - R et V sont TOUS LES DEUX en-champ,
        # orientations opposees (theta=+90/-90), et seule LA MOYENNE
        # (R+V)/2 annule le pTRM pour redonner le NRM restant. C'etait le
        # comportement D'ORIGINE de ce code (avant le fix errone) - remis
        # tel quel. La reconstruction (R+V)/2 necessaire pour PmagPy est
        # faite cote analyse (paleointensity_magic.py), PAS ici : ce
        # fichier doit representer la mesure BRUTE (R est reellement
        # mesuree en champ, ce serait mentir que d'ecrire "zero field").
        temp = etape + 273
        theta = 90.0
        dc_field = ifield * 1.0e-6
        if is_pure_ii_protocol:
            # Thellier classique (R/V/P sans 'S') : le champ est TOUJOURS
            # applique, il n'y a ni alternance zero-champ/en-champ ni
            # sequence IZZI - "LP-PI-TRM-ZI" et "LP-PI-BT-IZZI" (et
            # "LP-PI-ALT-PTRM", specifique au protocole IZZI) n'ont aucun
            # sens ici (demande explicite utilisateur : "LP-PI-TRM-ZI est
            # inutile en Thellier classic de meme que LP-PI-BT-IZZI").
            # Meme jeu que le pas V (deja "champ en Z puis Z-", protocole
            # LP-PI-II) - R et V sont les deux mesures en champ d'un meme
            # palier. Verifie que thellier_gui (sortarai) n'utilise ces
            # tags ZI/IZZI que pour apparier un pas zero-champ (LT-T-Z),
            # jamais present ici.
            codes = "LT-T-I:LP-PI-TRM:LP-PI-II:LP-PI-ALT"
        else:
            codes = "LT-T-I:LP-PI-TRM-ZI:LP-PI-TRM:LP-PI-ALT-PTRM:LP-PI-BT-IZZI"
    elif m.cod1 == "V":
        temp = etape + 273
        theta = -90.0
        dc_field = ifield * 1.0e-6
        codes = "LT-T-I:LP-PI-TRM:LP-PI-II:LP-PI-ALT"
    elif m.cod1 == "P":
        temp = etape + 273
        # Controle pTRM (bug corrige, signale par l'utilisateur sur des
        # donnees reelles - "the PTRM checks... are done in the field
        # direction of the R step. in this case it is not n.d but 90.0") :
        # le champ est reapplique dans l'orientation du pas R (theta=90),
        # jamais renseigne auparavant (restait a 0.0 -> "n.d" a l'ecriture,
        # cf. _THETA_CODES qui EXCLUT 'P').
        theta = 90.0
        last2 = [p.cod1 for p in prev[-2:]]
        if last2 == ["R", "V"] or last2 == ["V", "R"]:
            dc_field = ifield * 1.0e-6
            codes = "LT-PTRM-I:LP-PI-TRM:LP-PI-II:LP-PI-ALT"
        elif prev and prev[-1].cod1 == "R":
            # BUG CORRIGE (signale par l'utilisateur - "there is two ways
            # of doing this ... the second way is to do a PTRM check after
            # a R step, and then you do it in zero field") : ce controle
            # suit un pas 'R' (en champ) et est refait EN CHAMP NUL - c'est
            # le "pTRM tail check a temperature plus basse" officiel MagIC
            # LT-PTRM-Z ("After in laboratory field step, perform a zero
            # field cooling at a lower temperature", verifie sur le
            # vocabulaire en ligne, cf. www2.earthref.org/MagIC/method-
            # codes.json), PAS LT-PTRM-I ("After zero field step, perform
            # an in field cooling" - c'est l'AUTRE sens, tague dans le
            # `else` ci-dessous) : coder les deux LT-PTRM-I ici les
            # confondait dans `pmag.sortarai`, qui les repartit dans deux
            # listes DISTINCTES (ptrm_check vs zptrm_check) utilisees de
            # facon exclusive par PintPars pour DRAT/DRATS - un controle en
            # champ nul mal etiquete LT-PTRM-I y etait traite comme si son
            # propre moment (mesure SANS champ) etait un pTRM re-acquis,
            # faussant le calcul.
            dc_field = 0.0
            codes = "LT-PTRM-Z:LP-PI-TRM:LP-PI-ALT-PTRM:LP-PI-BT-IZZI"
        else:
            dc_field = ifield * 1.0e-6
            codes = "LT-PTRM-I:LP-PI-TRM:LP-PI-ALT-PTRM:LP-PI-BT-IZZI"
    elif m.cod1 in ("X", "Y", "Z"):
        theta = 90.0 if m.cod1 == "Z" else 0.0
        if m.cod1 == "X":
            phi = 180.0 if m.cod2 == "-" else 0.0
        elif m.cod1 == "Y":
            phi = 270.0 if m.cod2 == "-" else 90.0
        # LP-AN-TRM (anisotropie de TRM) - bug corrige, signale par
        # l'utilisateur ("the problem is with the anisotropy data.
        # :LP-AN-TRM should be added to the magic code during import") :
        # jamais tague auparavant, alors que convert_magic_to_r.py (sens
        # inverse) EXIGE deja "LP-AN-TRM" pour reconnaitre un pas ATRM a
        # l'import MagIC (voir _experiment_signature/_atrm_axis_sign) -
        # asymetrie corrigee ici, cote export/conversion.
        #
        # LP-AN-IRM (anisotropie d'IRM fort champ) AJOUTE ici - demande
        # explicite utilisateur (voir classify_anisotropy_experiment) :
        # `etape` represente alors un champ fort (mT), pas une
        # temperature - meme convention que cod1='I' (dc_field=etape*
        # 1e-3), PAS de treat_temp (contrairement au cas TRM).
        if anisotropy_kind == "irm":
            dc_field = etape * 0.001
            codes = "LT-IRM:LP-AN-IRM"
        else:
            temp = etape + 273
            # LT-T-I = pas EN CHAMP : le champ labo etait applique mais
            # treat_dc_field restait vide - thellier_gui/pmagpy lisent ce
            # champ (et phi/theta, ecrits meme a 0 des qu'un champ est
            # applique) pour ranger chaque pas ATRM dans l'une des 6
            # positions +/-X,Y,Z.
            dc_field = ifield * 1.0e-6
            codes = "LT-T-I:LP-AN-TRM"
    elif m.cod1 in ("L", "Q"):
        temp = etape + 273
        theta = 90.0
        dc_field = 0.0
        codes = "LT-T-I:LP-CR-TRM"
    else:
        temp = etape + 273
        codes = "LT-NO"

    return codes, temp, af_field, dc_field, phi, theta


def _izzi_order_tags(mesures: List[Measurement]) -> Dict[int, str]:
    """index de mesure -> "LP-PI-TRM-IZ"/"LP-PI-TRM-ZI" pour les pas R/S,
    deduit de l'ORDRE REEL de mesure (lequel des deux, R ou S, a la meme
    temperature, apparait en premier dans la sequence).

    BUG CORRIGE ici (repere en construisant paleointensity_magic.py, PAS
    dans le Fortran - propre a ce portage) : `_measurement_treatment`
    tague STATIQUEMENT cod1='S' avec "LP-PI-TRM-IZ" et cod1='R' avec
    "LP-PI-TRM-ZI", quel que soit l'ordre reel - or un protocole IZZI
    authentique ALTERNE les deux ordres d'une temperature a l'autre
    (verifie sur un vrai fichier IZZI, magic_contribution_19987.prmag,
    specimen kr01_01b3 : R-puis-S a 200degC, S-puis-R a 300degC). Ce tag
    est pourtant SPECIFIQUE PAR PAIRE dans le modele MagIC (lu par
    pmag.sortarai sur le pas en champ, par pmag.find_dmag_rec sur CHAQUE
    pas) - un tag toujours identique fausse silencieusement toute analyse
    downstream basee sur l'ordre IZ/ZI (ex. le test ZigZag de
    pmag.PintPars), y compris pour les fichiers deja exportes par
    `build_measurements_rows`."""
    tags: Dict[int, str] = {}
    for idx, m in enumerate(mesures):
        if m.cod1 not in ("R", "S"):
            continue
        partner_cod1 = "S" if m.cod1 == "R" else "R"
        partner_idx = next(
            (j for j, mm in enumerate(mesures) if mm.cod1 == partner_cod1 and mm.etape == m.etape),
            None,
        )
        if partner_idx is None:
            continue
        idx_R = idx if m.cod1 == "R" else partner_idx
        idx_S = partner_idx if m.cod1 == "R" else idx
        tags[idx] = "LP-PI-TRM-IZ" if idx_R < idx_S else "LP-PI-TRM-ZI"
    return tags


def build_measurements_rows(
    samples: List[SelectedSample], lab_analysts: str,
    anisotropy_kind_by_specimen: Optional[Dict[str, str]] = None,
    anisotropy_skip: Optional[set] = None,
) -> List[List[str]]:
    """`anisotropy_kind_by_specimen` ("trm"/"irm") et `anisotropy_skip`
    (id specimen -> ne pas archiver les pas X/Y/Z d'anisotropie) resolvent
    l'ambiguite TRM/IRM signalee par l'utilisateur pour les pas cod1 in
    (X,Y,Z) - voir classify_anisotropy_experiment. Un specimen absent de
    `anisotropy_kind_by_specimen` est traite en "trm" (comportement
    historique inchange)."""
    anisotropy_kind_by_specimen = anisotropy_kind_by_specimen or {}
    anisotropy_skip = anisotropy_skip or set()
    rows = []
    counter = 0
    for ech in samples:
        try:
            ifield = float(ech.com[:2])
        except (ValueError, TypeError):
            ifield = 0.0
        izzi_tags = _izzi_order_tags(ech.mesures)
        skip_anisotropy = ech.id in anisotropy_skip
        anisotropy_kind = anisotropy_kind_by_specimen.get(ech.id, "trm")
        is_paleointensity = any(m.cod1 in _PALEOINT_COMPANION_COD1 for m in ech.mesures)
        is_pure_ii_protocol = is_paleointensity and not any(m.cod1 == "S" for m in ech.mesures)
        # Thellier classique : au palier ou l'ATRM est mesuree (pas X/Y),
        # les pas R et V SONT les positions Z+ et Z- (voir
        # calcul.detect_six_positions, qui les utilise deja comme
        # substituts) - demande explicite utilisateur : "ajouter
        # LP-AN-TRM a la meme etape pour R et V". Limite au protocole
        # classique (R et V existent tous deux) et a l'ATRM (pas d'IRM,
        # pas d'anisotropie ecartee par l'utilisateur).
        # thellier_gui exige, des 3 lignes LP-CR-TRM ou plus, une ligne
        # LT-PTRM-I (controle d'alteration = refroidissement rapide
        # repete) sinon IndexError sur alteration_check[0] - la DERNIERE
        # ligne L/Q d'un specimen qui en a >= 3 joue ce role.
        cooling_idx = [j for j, m in enumerate(ech.mesures) if m.cod1 in ("L", "Q")]
        alteration_check_idx = cooling_idx[-1] if len(cooling_idx) >= 3 else None
        atrm_etapes = set()
        if is_pure_ii_protocol and not skip_anisotropy and anisotropy_kind == "trm":
            atrm_etapes = {m.etape for m in ech.mesures if m.cod1 in ("X", "Y")}

        for j, m in enumerate(ech.mesures):
            if skip_anisotropy and m.cod1 in ("X", "Y", "Z"):
                continue
            counter += 1
            mag, dec, inc = polere(m.x, m.y, m.z)
            if ech.norme == "m":
                magn_mass = mag * 1.0e3 / ech.vol if ech.vol else 0.0
                magn_volume = ""
                chi_mass = m.s * 1.0e-7 / ech.vol if (ech.vol and m.s) else ""
                chi_volume = ""
            else:
                magn_volume = mag * 1.0e6 / ech.vol if ech.vol else 0.0
                magn_mass = ""
                chi_volume = m.s * 10.0 * 1.0e-5 / ech.vol if (ech.vol and m.s) else ""
                chi_mass = ""

            codes, temp, af_field, dc_field, phi, theta = _measurement_treatment(
                m, ech.mesures[:j], ifield, anisotropy_kind, is_paleointensity, is_pure_ii_protocol)
            if j in izzi_tags:
                parts = [p for p in codes.split(":") if p not in ("LP-PI-TRM-IZ", "LP-PI-TRM-ZI")]
                parts.append(izzi_tags[j])
                codes = ":".join(parts)
            if m.cod1 in ("R", "V") and m.etape in atrm_etapes and "LP-AN-TRM" not in codes:
                codes += ":LP-AN-TRM"

            # dir_csd derive de m.q ("error") - demande explicite
            # utilisateur ("lors de l'importation, lorsque l'instrument
            # est S, mettre n.d lors de l'import et ne pas la prendre en
            # compte lors de l'exportation vers magic") : la definition
            # de ce facteur pour l'instrument "S" (spinner, donnees
            # anciennes) n'est pas comparable a celle des instruments 2G
            # cryo modernes (C1/C4/JA/J2) - ne PAS le calculer/exporter
            # pour ces mesures, plutot que de deriver un dir_csd sur une
            # base non comparable.
            ins_short = (m.ins or "").strip()
            csd = None if ins_short == "S" else 0.1 + math.degrees(math.atan2(m.q / 100.0, 1.0))
            ins_desc = _instrument(m.ins)
            # "<vitesse>:K/min" : format lu par thellier_gui (get_data, bloc
            # cooling rate : `description.split(":")`, valeur juste avant
            # "K/min") ; le texte qualitatif reste en tete.
            description = ""
            if m.cod1 == "L":
                description = f"slow cooling:{_LAB_SLOW_COOLING_K_PER_MIN:g}:K/min"
            elif m.cod1 == "Q":
                description = f"fast cooling:{_LAB_FAST_COOLING_K_PER_MIN:g}:K/min"
            if j == alteration_check_idx:
                codes = "LT-PTRM-I:LP-CR-TRM"

            row = [
                "This study", lab_analysts, ech.id, f"{ech.id}_Rem_Mag",
                "STARpaleomag_Py", "", f"{counter}-{ech.id}", "g", "u",
                str(j + 1),
                _num(temp, 1) if temp else "",
                _num(af_field, 6) if af_field else "",
                _num(dc_field, 6) if dc_field else "",
                # phi/theta = 0 est une VRAIE valeur des qu'un champ est
                # applique (dc_field != 0) - pas "inconnu". Ecrit "" pour 0,
                # _write_tsv supprimait la colonne entiere quand aucune
                # ligne n'avait phi != 0 (aucun pas d'anisotropie X/Y dans
                # l'export) et thellier_gui plantait avec KeyError
                # 'treatment_dc_field_phi' (get_data lit phi/theta sur CHAQUE
                # pas en champ LT-T-I/LT-PTRM-I) - signale par l'utilisateur.
                _num(phi, 1) if (phi or dc_field) else "",
                _num(theta, 1) if (theta or dc_field) else "",
                "293",
                _num(inc, 1), _num(dec, 1),
                _num_sci(mag),
                _num_sci(magn_volume) if magn_volume != "" else "",
                _num_sci(magn_mass) if magn_mass != "" else "",
                _num(csd, 2) if csd is not None else "",
                _num_sci(chi_volume) if chi_volume != "" else "",
                _num_sci(chi_mass) if chi_mass != "" else "",
                codes, ins_desc, description,
            ]
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

@dataclass
class MagicExportResult:
    paths: Dict[str, str] = field(default_factory=dict)
    counts: Dict[str, int] = field(default_factory=dict)


_THELLIER_GUI_ANISO_COLUMNS = (
    "aniso_s", "aniso_ftest", "aniso_ftest12", "aniso_s_n_measurements",
    "aniso_s_sigma", "aniso_type", "description",
)


def _write_tsv(path: str, table_name: str, header: List[str], rows: List[List[str]]) -> None:
    # Colonnes ENTIEREMENT vides sur TOUTES les lignes (ex. aucune donnee
    # de paleointensite dans tout l'export) supprimees avant ecriture -
    # demande explicite utilisateur ("est-ce possible d'éviter les
    # colonnes vides, par exemple si il n'y a pas de paleointensité,
    # éviter ces colonnes") : les en-tetes (_SITES_HEADER/_SPECIMENS_
    # HEADER/...) listent le SUR-ENSEMBLE de colonnes possibles (int_*/
    # aniso_*/vgp_*...), la plupart des exports n'en utilisant qu'un
    # sous-ensemble - une colonne ecrite alors qu'AUCUNE ligne ne la
    # renseigne n'ajoute aucune information, juste du bruit. Gardee des
    # qu'AU MOINS une ligne la renseigne (une valeur ponctuelle, meme
    # rare, reste une vraie donnee - jamais retiree).
    if rows:
        keep = [i for i in range(len(header)) if any(str(row[i]).strip() for row in rows)]
        # thellier_gui n'utilise un tenseur de specimens.txt que si TOUTES
        # ces colonnes existent (sinon "Incomplete anisotropy data ...
        # Ignoring anisotropy data") - gardees ensemble des qu'un tenseur
        # est exporte, meme si certaines valeurs (sigma, F-test) sont vides.
        if table_name == "specimens" and "aniso_s" in [header[i] for i in keep]:
            keep = sorted(set(keep) | {header.index(c) for c in _THELLIER_GUI_ANISO_COLUMNS if c in header})
    else:
        keep = list(range(len(header)))
    if len(keep) != len(header):
        header = [header[i] for i in keep]
        rows = [[row[i] for i in keep] for row in rows]
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(f"tab delimited\t{table_name}\n")
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(str(c) for c in row) + "\n")


def export_to_magic(
    samples: List[SelectedSample],
    results: List[FitResult],
    out_dir: str,
    lab_analysts: str = "",
    continent_ocean: str = "", country: str = "", region: str = "",
    anisotropy_kind_by_specimen: Optional[Dict[str, str]] = None,
    anisotropy_skip: Optional[set] = None,
    pmagint_rows: Optional[Dict[str, Dict[str, str]]] = None,
    aniso_tensors: Optional[Dict[str, AniTensor]] = None,
    aniso_mean_tensors: Optional[Dict[str, AniMeanTensor]] = None,
    site_metadata: Optional[Dict[str, Dict[str, str]]] = None,
    include_site_only_means: bool = True,
) -> MagicExportResult:
    """Equivalent (mode classique, ichoixexport==1) de `export2magic`, PLUS
    les resultats de paleointensite deja archives dans .pmagint
    (`pmagint_rows`, voir build_specimens_rows/paleointensity_magic_fields)
    - le Fortran d'origine traitait ca dans un mode SEPARE (ichoixexport==2,
    export_int_2_magic), hors perimetre ici jusqu'a cette demande explicite
    utilisateur ("in export to Magic, it is not asking for Paleointensity") ;
    fusionne desormais dans le MEME specimens.txt plutot qu'un fichier a
    part, coherent avec la convention MagIC v3 (une ligne specimen peut
    porter a la fois dir_*/int_*/aniso_* quand ils viennent de la meme
    sequence Thellier/IZZI). `aniso_tensors` (specimen -> AniTensor 'A0'
    deja calcule, voir calcul.read_ani_tensor) : meme principe, colonnes
    aniso_s/aniso_ftest* - demande explicite utilisateur ("il manque aussi
    les données d'anisotropie au niveau du fichier specimens (prendre les
    A0)"). `aniso_mean_tensors` (site -> AniMeanTensor 'A0' deja calcule,
    voir calcul.read_ani_mean_tensor) : meme principe au niveau SITE
    plutot que specimen, colonnes aniso_v1/v2/v3/aniso_p... dans
    sites.txt - demande explicite utilisateur ("l'exportation de l'AMS
    (niveau specimen et sites). Les données d'AMS vont aussi dans les
    fichiers sites.txt et specimens.txt").
    `samples` est trie par (magic_site, magic_sample, id) avant traitement
    - voir l'ecart documente en tete de module."""
    ordered = sorted(
        samples,
        key=lambda e: (e.magic_site.strip(), e.magic_sample.strip(), e.id.strip()),
    )

    os.makedirs(out_dir, exist_ok=True)

    sites_rows = build_sites_rows(
        ordered, results, site_metadata=site_metadata, aniso_mean_tensors=aniso_mean_tensors,
        include_site_only_means=include_site_only_means)
    locations_rows = build_locations_rows(ordered, continent_ocean, country, region)
    samples_rows = build_samples_rows(ordered)
    specimens_rows = build_specimens_rows(
        ordered, results, pmagint_rows=pmagint_rows, aniso_tensors=aniso_tensors)
    measurements_rows = build_measurements_rows(
        ordered, lab_analysts, anisotropy_kind_by_specimen, anisotropy_skip)

    files = [
        ("sites.txt", "sites", _SITES_HEADER, sites_rows),
        ("locations.txt", "locations", _LOCATIONS_HEADER, locations_rows),
        ("samples.txt", "samples", _SAMPLES_HEADER, samples_rows),
        ("specimens.txt", "specimens", _SPECIMENS_HEADER, specimens_rows),
        ("measurements.txt", "measurements", _MEASUREMENTS_HEADER, measurements_rows),
    ]

    result = MagicExportResult()
    for filename, table_name, header, rows in files:
        path = os.path.join(out_dir, filename)
        _write_tsv(path, table_name, header, rows)
        result.paths[table_name] = path
        result.counts[table_name] = len(rows)
    return result
