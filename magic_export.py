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
from calcul import FitResult, AniTensor

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
    "dir_k", "dir_n_samples", "dir_n_specimens_lines",
    "dir_n_specimens_planes", "dir_polarity", "dir_nrm_origin",
    "vgp_lat", "vgp_lon", "vgp_dp", "vgp_dm",
]

_LOCATIONS_HEADER = [
    "location", "geologic_classes", "lithologies", "age", "age_sigma",
    "age_low", "age_high", "age_unit", "lat_n", "lat_s", "lon_e", "lon_w",
    "citations", "continent_ocean", "country", "description",
    "location_type", "region",
]

_COMP_NAMES = {0: "Comp_A", 1: "Comp_A", 2: "Comp_B", 3: "Comp_C"}


def _site_mean_row(site: str, results: List[FitResult]) -> Optional[Dict[str, str]]:
    """Equivalent de `magicmeanres` : cherche une moyenne de site
    (cat1=='F' via id "mean: <site>") correspondant a `site` (magic_site,
    ou son equivalent tronque a 6 caracteres comme le convention interne
    des id "mean: XXXXXX")."""
    site6 = (site or "").strip()[:6]
    for r in results:
        if r.id[:5] != "mean:":
            continue
        if r.id[6:12].strip() != site6:
            continue
        tilt = "0" if r.par3_mean == 2.0 else "100"
        comp = _COMP_NAMES.get(r.numcomp, "Comp_A")
        specimens = ":".join(
            t.strip() for t in r.liste.replace("codes:", "").split(":") if t.strip()
        )
        return {
            "specimens": specimens,
            "result_quality": "g",
            "dir_comp_name": comp,
            "dir_tilt_correction": tilt,
            "dir_dec": _num(r.dec, 1),
            "dir_inc": _num(r.inc, 1),
            "dir_alpha95": _num(r.mad, 1),
            "dir_k": _num(r.par2_mean, 1),
            "dir_n_samples": str(r.nb),
            "dir_n_specimens_lines": str(int(r.tx[0])),
            "dir_n_specimens_planes": str(int(r.tx[1])),
            "vgp_lat": _num(r.par4, 1),
            "vgp_lon": _num(r.par5, 1),
            "vgp_dp": _num(r.vgp_dp, 1) if r.vgp_dp else "",
            "vgp_dm": _num(r.vgp_dm, 1) if r.vgp_dm else "",
        }
    return None


def build_sites_rows(samples: List[SelectedSample], results: List[FitResult]) -> List[List[str]]:
    rows = []
    seen = []
    for ech in samples:
        site = ech.magic_site.strip()
        if not site or site in seen:
            continue
        seen.append(site)

        age, age_sigma, age_low, age_high, age_unit = parse_age(ech.magic_age)
        elevation = "" if ech.altitude <= 0.0 else _num(ech.altitude, 1)
        row = {
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
        mean = _site_mean_row(site, results)
        if mean:
            row.update(mean)
        rows.append([row.get(col, "") for col in _SITES_HEADER])
    return rows


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

_SAMPLES_HEADER = [
    "citations", "sample", "site", "geologic_classes", "lithologies",
    "geologic_types", "lat", "lon", "height", "timestamp", "orientation_quality",
    "azimuth", "dip", "bed_dip_direction", "bed_dip", "method_codes",
]


def build_samples_rows(samples: List[SelectedSample]) -> List[List[str]]:
    rows = []
    seen = []
    for ech in samples:
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

# code2 ('A0' seul demande explicitement - "prendre les A0") -> (aniso_type,
# method_codes) MagIC - MEME convention que AMS_Py._ANI_MAGIC_INFO (A0=TRM).
_ANISO_MAGIC_INFO = {"A0": ("ATRM", "LP-AN-TRM")}


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
    }
    if tensor.n_positions is not None:
        out["aniso_s_n_measurements"] = str(tensor.n_positions)
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
            comp = f"ChRM_{r.numcomp}"

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
    "C1": ("2G_Cryo 1 position", "1"),
    "C4": ("2G Cryo 4 positions", "4"),
    "JA": ("JR5 or JR6 automatic", "3"),
    "J2": ("JR5 or JR6 2 positions", "2"),
}


def _instrument(ins: str) -> str:
    ins = (ins or "").strip()
    if ins in _INSTRUMENT_CODES:
        return _INSTRUMENT_CODES[ins][0]
    if ins[:1] == "M":
        return "Molspin 6 positions"
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
        codes = "LT-NO"
        if m.cod2 == "P":
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
                m, ech.mesures[:j], ifield, anisotropy_kind)
            if j in izzi_tags:
                parts = [p for p in codes.split(":") if p not in ("LP-PI-TRM-IZ", "LP-PI-TRM-ZI")]
                parts.append(izzi_tags[j])
                codes = ":".join(parts)

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
            description = "slow cooling" if m.cod1 == "L" else ("fast cooling" if m.cod1 == "Q" else "")

            row = [
                "This study", lab_analysts, ech.id, f"{ech.id}_Rem_Mag",
                "Rennes_Pmag_Starmac", "", f"{counter}-{ech.id}", "g", "u",
                str(j + 1),
                _num(temp, 1) if temp else "",
                _num(af_field, 6) if af_field else "",
                _num(dc_field, 6) if dc_field else "",
                _num(phi, 1) if phi else "",
                _num(theta, 1) if theta else "",
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


def _write_tsv(path: str, table_name: str, header: List[str], rows: List[List[str]]) -> None:
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
    A0)").
    `samples` est trie par (magic_site, magic_sample, id) avant traitement
    - voir l'ecart documente en tete de module."""
    ordered = sorted(
        samples,
        key=lambda e: (e.magic_site.strip(), e.magic_sample.strip(), e.id.strip()),
    )

    os.makedirs(out_dir, exist_ok=True)

    sites_rows = build_sites_rows(ordered, results)
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
