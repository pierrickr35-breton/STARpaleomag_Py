#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import math
import os
import re
import sys
import tkinter as tk
from tkinter import filedialog
from typing import Dict, List, Optional, Tuple
import pandas as pd

# Protocoles MagIC reels mais HORS PERIMETRE STARpaleomag_Py (rock-magnetisme, pas
# demagnetisation/paleointensite) - demande explicite utilisateur ("Magic
# file also contain Anisotropy of magnetic susceptibility and other
# magnetic experiments that are not processed by STARpaleomag_Py"). Definie ICI
# (pas dans convert_magic_to_r.py, qui l'utilisait a l'origine) pour eviter
# un import circulaire : convert_magic_to_r.py importe aussi depuis ce
# module (magic_results_to_redo_lines/magic_site_means, voir plus bas) -
# ce module-ci ne doit rien importer en retour de convert_magic_to_r.py.
_OUT_OF_SCOPE_PROTOCOLS = {
    "LP-AN-MS": "anisotropy of magnetic susceptibility (AMS)",
    "LP-AN-ARM": "anisotropy of ARM",
    "LP-AN-IRM": "anisotropy of IRM",
    "LP-HYS": "hysteresis loop",
    # LP-BCR (coercivity of remanence/backfield) et LP-IRM (courbe
    # d'acquisition IRM) NE SONT PAS hors perimetre - demande explicite
    # utilisateur ("LP-BCR, LP-IRM are measurements that can be processed
    # by STARpaleomag_Py") : leurs pas utilisent le meme code de traitement
    # LT-IRM deja gere (cod1='I') que le controle IRM ponctuel d'une
    # experience de paleointensite - aucun mapping supplementaire requis,
    # juste ne plus les ecarter en amont.
    "LP-X": "susceptibility characterization (vs field/freq/temp)",
    # MPMS (Quantum Design, instrument_codes typiquement "IRM-MPMS3") -
    # cycles thermiques haut-champ SANS RAPPORT avec la NRM - demande
    # explicite utilisateur ("in the magic file there are MPMS high field
    # data not related to the natural remanent magnetization. Do not
    # import these data"). Verifie sur donnees reelles
    # (magic_contribution_20340.txt) : 1184 mesures MPMS/29425, aucune ne
    # matchait un prefixe LT-* reconnu par convert_magic_to_r._derive_cod
    # avant ce correctif - elles tombaient dans le repli "non reconnu"
    # (cod1='?'), au lieu d'etre ecartees en amont comme les autres
    # protocoles hors perimetre.
    "LP-ZFC": "MPMS zero-field-cooled low-temperature cycling",
    "LP-FC": "MPMS field-cooled low-temperature cycling",
    "LP-CW-SIRM": "MPMS SIRM cooling/warming (e.g. Verwey transition)",
}


def clean_str(val):
    """Nettoie une chaîne de tous les espaces superflus et espaces insecables."""
    if pd.isna(val) or val is None:
        return ""
    return str(val).replace("\xa0", " ").strip()


def read_magic_file(filepath):
    """Lit un fichier MAGIC individuel (sites.txt/samples.txt/...) en
    sautant systématiquement la 1ère ligne d'en-tête de format."""
    try:
        df = pd.read_csv(
            filepath, sep="\t", skiprows=1, header=0, dtype=str
        )
        return _clean_columns(df)
    except Exception as e:
        print(f"❌ Read error ({os.path.basename(filepath)}): {e}")
        return None


def _clean_columns(df):
    df.columns = df.columns.astype(str).str.strip().str.lower()
    df.columns = (
        df.columns.str.replace(" ", "_")
        .str.replace("#", "")
        .str.replace('"', "")
    )
    for col in df.columns:
        df[col] = df[col].apply(clean_str)
    return df


def split_combined_magic_file(filepath):
    """Scinde un fichier de contribution MagIC combiné (le format
    téléchargé depuis earthref.org/MagIC : une seule table par bloc,
    chaque bloc précédé d'une ligne "tab delimited\\t<nom>" et séparé du
    suivant par une ligne ">>>>>>>>>>") en un dict {nom_de_table:
    DataFrame} - mêmes colonnes nettoyées que read_magic_file() sur un
    fichier individuel, pour rester utilisable de façon interchangeable."""
    with open(filepath, "r", encoding="utf-8-sig", errors="replace") as f:
        content = f.read()

    tables = {}
    for block in content.split(">>>>>>>>>>"):
        block_lines = block.strip("\n").splitlines()
        if not block_lines or not block_lines[0].lower().startswith("tab delimited"):
            continue
        header_parts = block_lines[0].split("\t", 1)
        table_name = header_parts[1].strip().lower() if len(header_parts) > 1 else ""
        rest = "\n".join(block_lines[1:]).strip("\n")
        if not table_name or not rest:
            continue
        try:
            df = pd.read_csv(io.StringIO(rest), sep="\t", header=0, dtype=str)
        except Exception as e:
            print(f"⚠️ Could not read table « {table_name} »: {e}")
            continue
        tables[table_name] = _clean_columns(df)
    return tables


def parse_float_val(val_str, default=0.0):
    """Convertit une chaîne en float de manière sécurisée."""
    try:
        return float(val_str)
    except (ValueError, TypeError):
        return default


def calc_xyz(dec_val, inc_val, mag_val):
    """Calcule les composantes cartésiennes (X, Y, Z)."""
    try:
        dec_rad = math.radians(dec_val)
        inc_rad = math.radians(inc_val)

        x = mag_val * math.cos(inc_rad) * math.cos(dec_rad)
        y = mag_val * math.cos(inc_rad) * math.sin(dec_rad)
        z = mag_val * math.sin(inc_rad)

        return x, y, z
    except (ValueError, TypeError):
        return 0.0, 0.0, 0.0


def format_1pe10_3(val):
    """Formate un nombre flottant selon la convention Fortran 1PE10.3."""
    s = f"{val:.3E}"
    if "E+" in s:
        base, exp = s.split("E+")
        s = f"{base}E+{int(exp):02d}"
    elif "E-" in s:
        base, exp = s.split("E-")
        s = f"{base}E-{int(exp):02d}"
    return f"{s:>10}"


def build_entity_dict(df, key_candidates):
    """Construit un dictionnaire d'entités indexé par leur identifiant (en minuscules)."""
    if df is None or df.empty:
        return {}

    key_col = None
    for candidate in key_candidates:
        if candidate in df.columns:
            key_col = candidate
            break

    if not key_col:
        return {}

    result_dict = {}
    for _, row in df.iterrows():
        key_val = clean_str(row[key_col]).lower()
        if key_val:
            result_dict[key_val] = row.to_dict()

    return result_dict


def get_value_from_dicts(dict_list, candidates, default="n.d"):
    """Cherche la valeur dans une liste de dictionnaires selon une liste de clés possibles."""
    for d in dict_list:
        if not d:
            continue
        for key in candidates:
            if key in d:
                val = clean_str(d[key])
                if val and val.lower() != "nan":
                    return val
    return default


def extract_volume(spec_dict, sample_dict, row_dict, default_vol=9.99):
    """Extrait le volume avec repli prioritaire : Specimen -> Sample -> Mesure -> Défaut."""
    vol_str = ""

    for d in [spec_dict, sample_dict, row_dict]:
        if not d:
            continue
        for k, v in d.items():
            if "vol" in k:
                val = clean_str(v)
                if val and val.lower() != "nan":
                    vol_str = val
                    break
        if vol_str:
            break

    vol_val = parse_float_val(vol_str, default=0.0)

    if vol_val == 0.0:
        return default_vol
    elif 0.0 < vol_val < 1.0:
        return vol_val * 1e6
    else:
        return vol_val


def format_custom_output(
    sites_file=None,
    samples_file=None,
    specimens_file=None,
    measurements_file=None,
    combined_file=None,
    output_file="extracted_data_custom.ren",
    default_volume=9.99,
):
    """`combined_file` : fichier de contribution MagIC unique (format
    téléchargé depuis earthref.org, voir split_combined_magic_file) - mode
    alternatif aux 4 fichiers séparés (sites_file/samples_file/
    specimens_file/measurements_file), fourni en priorité si présent."""
    print("📂 Loading and processing MAGIC files...")

    if combined_file:
        tables = split_combined_magic_file(combined_file)
        sites_df = tables.get("sites")
        samples_df = tables.get("samples")
        specimens_df = tables.get("specimens")
        measurements_df = tables.get("measurements")
    else:
        sites_df = read_magic_file(sites_file)
        samples_df = read_magic_file(samples_file)
        specimens_df = read_magic_file(specimens_file)
        measurements_df = read_magic_file(measurements_file)

    if measurements_df is None or measurements_df.empty:
        print("❌ Measurements table not found or empty.")
        return

    specimens_dict = build_entity_dict(
        specimens_df, ["specimen", "specimen_name", "specimen_id"]
    )
    samples_dict = build_entity_dict(
        samples_df, ["sample", "sample_name", "sample_id"]
    )
    sites_dict = build_entity_dict(
        sites_df, ["site", "site_name", "site_id"]
    )

    meas_spec_col = None
    for cand in ["specimen", "specimen_name", "specimen_id"]:
        if cand in measurements_df.columns:
            meas_spec_col = cand
            break

    if not meas_spec_col:
        print(f"❌ Specimen column not found. Available columns: {list(measurements_df.columns)}")
        return

    specimen_groups = measurements_df.groupby(meas_spec_col, sort=False)
    output_lines = []
    n_excluded_rows = 0
    n_excluded_specimens = 0

    for specimen_name, group in specimen_groups:
        # Ecarte les mesures MPMS/rock-magnetisme HORS PERIMETRE (memes
        # protocoles que convert_magic_to_r.py - demande explicite
        # utilisateur : "in the magic file there are MPMS high field data
        # not related to the natural remanent magnetization. Do not
        # import these data"). Cette fonction (contrairement a
        # convert_magic_to_r.py) classe cod1/cod2 UNIQUEMENT depuis
        # treat_ac_field/treat_temp BRUTS (pas method_codes) - sans ce
        # filtre, un cycle thermique MPMS (treat_temp non nul) etait
        # silencieusement pris pour un pas de demagnetisation thermique
        # ordinaire (cod1='D'), pas juste rejete comme "non reconnu".
        if "method_codes" in group.columns:
            codes_col = group["method_codes"].fillna("")
            mask = ~codes_col.apply(lambda c: any(p in c for p in _OUT_OF_SCOPE_PROTOCOLS))
            n_excluded_rows += int((~mask).sum())
            if not mask.all():
                group = group[mask]
            if group.empty:
                n_excluded_specimens += 1
                continue

        first_row = group.iloc[0].to_dict()

        spec_key = clean_str(specimen_name).lower()
        spec_dict = specimens_dict.get(spec_key, {})

        sample_name = spec_dict.get(
            "sample", first_row.get("sample", "n.d")
        )
        sample_dict = samples_dict.get(clean_str(sample_name).lower(), {})

        site_name = spec_dict.get("site") or sample_dict.get(
            "site", first_row.get("site", "n.d")
        )
        site_dict = sites_dict.get(clean_str(site_name).lower(), {})

        dict_chain = [spec_dict, sample_dict, site_dict, first_row]

        specimen_id = clean_str(specimen_name)[:12]

        core_dip_raw = parse_float_val(
            get_value_from_dicts(
                dict_chain,
                ["dip", "core_dip", "orientation_dip", "dip_sample"],
                "0.0",
            )
        )
        core_az_raw = parse_float_val(
            get_value_from_dicts(
                dict_chain,
                [
                    "azimuth",
                    "core_azimuth",
                    "orientation_azimuth",
                    "azimuth_sample",
                ],
                "0.0",
            )
        )
        bed_dip_raw = parse_float_val(
            get_value_from_dicts(dict_chain, ["bed_dip"], "0.0")
        )
        # bug corrige (signale par l'utilisateur, meme erreur trouvee et
        # corrigee dans StereoUtils_Py) : MagIC ne fournit QU'UNE SEULE
        # variable de ce type, "bed_dip_direction" (pas de "bed_strike" -
        # confirme par l'utilisateur) - l'ancien code cherchait les deux
        # comme interchangeables puis appliquait +90 SANS CONDITION plus
        # bas (ligne `str_trans = bed_str_raw + 90`), alors que la bonne
        # conversion dip_direction -> strike est un DECALAGE DE -90
        # (dip_direction = strike + 90, verifie/corrige avec la meme
        # convention que StereoUtils_Py).
        bed_str_raw = (
            parse_float_val(
                get_value_from_dicts(dict_chain, ["bed_dip_direction"], "0.0")
            )
            - 90.0
        ) % 360.0

        volume_cm3 = extract_volume(
            spec_dict, sample_dict, first_row, default_vol=default_volume
        )
        vol_str = (
            f"{volume_cm3:5.2f}"
            if abs(volume_cm3 - default_volume) < 1e-3
            else f"{volume_cm3:5.1f}"
        )

        lat_raw = parse_float_val(
            get_value_from_dicts(
                dict_chain, ["latitude", "lat", "site_lat"], "0.0"
            )
        )
        lon_raw = parse_float_val(
            get_value_from_dicts(
                dict_chain, ["longitude", "lon", "site_lon"], "0.0"
            )
        )

        fm = get_value_from_dicts(
            dict_chain, ["geologic_formations", "formation"]
        )
        age = get_value_from_dicts(dict_chain, ["geologic_age", "age"])
        
        # MAPPINGS DE MÉTADONNÉES MIS À JOUR
        gc = get_value_from_dicts(dict_chain, ["geologic_classes", "geologic_class"])
        smt = get_value_from_dicts(dict_chain, ["geologic_types", "geologic_type"])
        li = get_value_from_dicts(dict_chain, ["lithology", "lithologies"])
        loc = get_value_from_dicts(dict_chain, ["locality"])

        inc_trans = -core_dip_raw
        az_trans = (core_az_raw + 90.0) % 360.0
        # bed_str_raw est deja la vraie valeur de strike (conversion faite
        # plus haut si la source etait bed_dip_direction) - pas de decalage
        # supplementaire ici (l'ancien +90 double-comptait la conversion).
        str_trans = bed_str_raw % 360.0

        header_1 = (
            f"Id:{specimen_id:<12} "
            f"in:{inc_trans:5.1f} "
            f"az:{az_trans:5.1f} "
            f"dip:{bed_dip_raw:5.1f} "
            f"str:{str_trans:5.1f} "
            f"v:{vol_str} "
            f"com:            "
        )

        # H/T/azm/azs : certaines variables ne sont pas toujours
        # disponibles selon le laboratoire d'origine (ex. les fichiers
        # samples.txt/specimens.txt reels contiennent souvent "height"
        # mais rarement un timestamp de terrain ou des azimuts bruts
        # separes) - on les extrait quand elles existent, "n.d"/0 sinon,
        # plutot que de toujours ecrire des placeholders comme avant.
        height_raw = get_value_from_dicts(
            dict_chain, ["height", "elevation", "site_elevation"], ""
        )
        try:
            height_val = float(height_raw)
        except (ValueError, TypeError):
            height_val = None
        h_str = f"{height_val:.1f}" if height_val is not None else "0"

        timestamp_raw = get_value_from_dicts(
            dict_chain, ["timestamp", "date_collected", "date"], ""
        )
        year = month = day = hour = minute = 0
        m = re.match(r"(\d{4})-(\d{2})-(\d{2})[T ](\d{2}):(\d{2})", timestamp_raw or "")
        if m:
            year, month, day, hour, minute = (int(g) for g in m.groups())

        azm_raw = get_value_from_dicts(
            dict_chain, ["azimuth_mag", "azm", "magnetic_azimuth"], ""
        )
        azs_raw = get_value_from_dicts(
            dict_chain, ["azimuth_sun", "azs", "sun_azimuth"], ""
        )
        azm_str = azm_raw if azm_raw else "n.d"
        azs_str = azs_raw if azs_raw else "n.d"

        header_2 = (
            f"L: {lat_raw:8.5f} G: {lon_raw:8.5f}  "
            f"H: {h_str:>3}  T:{year:4d}  {month:2d} {day:2d} {hour:2d} {minute:2d}   "
            f"azm:{azm_str:>5} azs:{azs_str:>5}  Or:A13_0_3_90"
        )

        header_3 = (
            f'Site: "{site_name}" Sample: "{sample_name}" Fm: "{fm}" Age: "{age}" '
            f'GC: "{gc}" SMT: "{smt}" Li: "{li}" Loc: "{loc}" Obs: "none"'
        )

        output_lines.append(header_1)
        output_lines.append(header_2)
        output_lines.append(header_3)

        for _, row_series in group.iterrows():
            row = row_series.to_dict()

            tac_raw = parse_float_val(
                get_value_from_dicts(
                    [row], ["treat_ac_field_meas", "treat_ac_field"], "0.0"
                )
            )
            ttemp_raw = parse_float_val(
                get_value_from_dicts(
                    [row], ["treat_temp_meas", "treat_temp"], "0.0"
                )
            )
            mag = parse_float_val(
                get_value_from_dicts(
                    [row], ["magn_moment_meas", "magn_moment"], "1.0"
                )
            )
            dec = parse_float_val(
                get_value_from_dicts(
                    [row], ["dir_dec_meas", "dir_dec"], "0.0"
                )
            )
            inc = parse_float_val(
                get_value_from_dicts(
                    [row], ["dir_inc_meas", "dir_inc"], "0.0"
                )
            )

            ttemp_c = ttemp_raw - 273.0 if ttemp_raw > 200 else ttemp_raw

            if tac_raw > 0:
                val_mt = tac_raw * 1000.0 if tac_raw < 1.0 else tac_raw
                ietape = int(round(val_mt * 10.0))
                code1 = "F"
                code2 = "="
            elif ttemp_c > 0:
                ietape = int(round(ttemp_c))
                code1 = "D"
                code2 = "="
            else:
                ietape = 0
                code1 = "N"
                code2 = "O"

            x, y, z = calc_xyz(dec, inc, mag)

            q_val = 0
            nd_str = "nd"
            susc_val = 0.0

            line_str = (
                f" {ietape:>4}{code1}{code2}"
                f" {format_1pe10_3(x)} {format_1pe10_3(y)} {format_1pe10_3(z)}"
                f" {q_val:>3} {nd_str:>2}{susc_val:7.1f}"
            )
            output_lines.append(line_str)

        output_lines.append("")

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(output_lines))

    if n_excluded_rows:
        print(f"⚠️  Excluded {n_excluded_rows} out-of-scope measurement(s) "
              f"(MPMS/rock-magnetism protocols, {n_excluded_specimens} specimen(s) fully excluded).")
    print(f"✅ Export finished at: {output_file}")
    return output_file


# ---------------------------------------------------------------------------
# Import des résultats déjà interprétés (specimens.txt) vers un fichier
# "redo" (voir calcul.fit_from_redo_file) - permet de rejouer dans STARpaleomag_Py
# les meilleurs ajustements (lignes/plans) déjà calculés par pmagpy/demag_gui.
# ---------------------------------------------------------------------------

def _convert_step_range(row):
    """Convertit meas_step_min/max vers la convention d'étape STARpaleomag_Py selon
    meas_step_unit : 'T'/'Tesla' -> mT (etape F, ex. 0.018T -> 18.0 mT -
    valeur physique reelle directement, plus d'echelle Oersted-equivalente
    depuis "convert all step integer to float"), 'K'/'Kelvin' -> degC
    (etape D, ex. 573K -> 300). Le K->degC est fait par SOUSTRACTION de
    273 (comme `ttemp_c = ttemp_raw - 273.0` déjà utilisé plus haut dans
    ce fichier pour la même conversion sur les mesures elles-mêmes) - pas
    par addition. `meas_step_unit` est un champ texte LIBRE du data model
    3 (pas de vocabulaire controle) : comparaison par PREFIXE (startswith
    "K"/"T") plutot qu'egalite stricte - bug reel corrige, une vraie
    contribution MagIC (magic_contribution_19491.txt) ecrit "Kelvin" en
    toutes lettres, jamais reconnu par l'egalite stricte a "K"
    (magic_results_to_redo_lines retournait silencieusement une liste vide
    sur ce fichier). Retourne (None, None) si non convertible."""
    unit = clean_str(row.get("meas_step_unit", "")).upper()
    smin = parse_float_val(row.get("meas_step_min", ""), None)
    smax = parse_float_val(row.get("meas_step_max", ""), None)
    if smin is None or smax is None:
        return None, None
    if unit.startswith("T"):
        return smin * 1000.0, smax * 1000.0
    if unit.startswith("K"):
        return smin - 273.0, smax - 273.0
    return smin, smax


def _classify_fit(method_codes):
    """(bestfit, ancr, comp) a partir de method_codes, ou (None,None,None)
    si aucun ajustement de type ligne/plan n'est present :
    - DE-BFL-A ou DE-BFL-O -> ('L','o',1) ; DE-BFL (seul) -> ('L','n',2)
    - DE-BFP-G -> ('P','o',1) ; DE-BFP (seul) -> ('P','n',2)"""
    codes = {c.strip() for c in (method_codes or "").split(":") if c.strip()}
    if "DE-BFL-A" in codes or "DE-BFL-O" in codes:
        return "L", "o", 1
    if "DE-BFL" in codes:
        return "L", "n", 2
    if "DE-BFP-G" in codes:
        return "P", "o", 1
    if "DE-BFP" in codes:
        return "P", "n", 2
    return None, None, None


def magic_results_to_redo_lines(specimens_source, combined=False) -> list:
    """Meme extraction que import_specimens_results_to_redo (lit
    specimens.txt ou la table "specimens" d'une contribution combinee,
    filtre result_quality=='g' + method_codes de type ligne/plan, deduplique
    par (specimen, ajustement)), mais retourne les lignes "redo" EN MEMOIRE
    plutot que de les ecrire dans un fichier - reutilisee a la fois par
    l'ecriture d'un fichier redo (toujours utile pour en preparer un a la
    main/le revoir avant de lancer les ajustements) ET par l'archivage
    DIRECT vers .pmagres, sans fichier intermediaire - demande explicite
    utilisateur ("with new format file for the results, i think it is
    possible to fill the pmagres file during import and not going through
    a redo file"). Retourne None si la table specimens est introuvable/vide
    (distinct d'une liste vide, qui signifie "aucun resultat exploitable")."""
    if combined:
        tables = split_combined_magic_file(specimens_source)
        df = tables.get("specimens")
    else:
        df = read_magic_file(specimens_source)

    if df is None or df.empty:
        print("❌ Specimens table not found or empty.")
        return None

    seen = set()
    lines = []
    for _, row_series in df.iterrows():
        row = row_series.to_dict()
        if clean_str(row.get("result_quality", "")).lower() != "g":
            continue
        bestfit, ancr, comp = _classify_fit(row.get("method_codes", ""))
        if bestfit is None:
            continue
        smin, smax = _convert_step_range(row)
        if smin is None:
            continue
        specimen = clean_str(row.get("specimen", ""))[:12]
        if not specimen:
            continue
        key = (specimen, bestfit, ancr, smin, smax, comp)
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"{specimen:<12} {bestfit}  {ancr}    {smin} {smax} {comp}")
    return lines


# Variantes LP-PI-* (vocabulaire MagIC, www2.earthref.org/MagIC/
# method-codes.json) qui ne sont PAS un protocole Thellier/IZZI
# thermique classique - AF, micro-ondes, multi-specimen, ou un proxy
# alternatif (ARM/IRM/susceptibilite) au lieu d'une vraie sequence
# TRM/four - jamais rejouables par le natif STARpaleomag_Py (afficher_arai
# attend des paliers thermiques R/V/P) meme si la ligne porte aussi le
# tag generique "LP-PI" - voir magic_pint_results_to_redo_lines.
_PI_NON_TRM_CODES = {
    "LP-PI-AFAF", "LP-PI-ARM", "LP-PI-IRM", "LP-PI-X", "LP-PI-PARM",
    "LP-PI-REL", "LP-PI-REL-PT", "LP-PI-SXTAL", "LP-PI-CHEM",
    "LP-PI-TRIAXE", "LP-PI-MULT", "LP-PI-MULT-DB", "LP-PI-MULT-FL",
    "LP-PI-M", "LP-PI-M-II", "LP-PI-M-IZ", "LP-PI-M-IZZI",
    "LP-PI-M-PERP", "LP-PI-M-QP", "LP-PI-M-ZI",
}


def _is_accepted_pi_codes(codes: set) -> bool:
    """Un ensemble de method_codes designe-t-il une determination
    Thellier/IZZI thermique rejouable par le natif STARpaleomag_Py
    (afficher_arai) ? "LP-PI-TRM"(-*) ou "LP-PI" nu (sans code
    _PI_NON_TRM_CODES) ou "LP-PI-BT-IZZI" - ce dernier code de vocabulaire
    controle MagIC reel (voir pmagpy/data_model/er_methods.txt,
    "lab_protocol") designe un protocole IZZI thermique evalue par
    Bootstrap, PAS dans _PI_NON_TRM_CODES (rejouable de la meme facon
    qu'un LP-PI-TRM/LP-PI nu) - BUG REEL corrige (demande explicite
    utilisateur "il n'est pas propose de redo"), confirme sur
    magic_contribution_20522.txt/20480.txt (specimens tagues UNIQUEMENT
    "LP-PI-BT-IZZI", jamais "LP-PI-TRM" ni "LP-PI" nu au niveau
    specimens.txt, jusqu'ici tous silencieusement rejetes)."""
    is_trm_variant = any(c == "LP-PI-TRM" or c.startswith("LP-PI-TRM-") for c in codes)
    is_generic_pi = "LP-PI" in codes and not (codes & _PI_NON_TRM_CODES)
    is_bt_izzi = "LP-PI-BT-IZZI" in codes
    return is_trm_variant or is_generic_pi or is_bt_izzi


def magic_pint_results_to_redo_lines(specimens_source, combined=False) -> list:
    """Meme principe que `magic_results_to_redo_lines` (specimens.txt deja
    interprete -> lignes "redo"), pour les determinations de PALEOINTENSITE
    plutot que les ajustements directionnels ligne/plan - demande explicite
    utilisateur ("during import of Magic paleoint data, can you write a
    redo_pint_contribution_num.txt file with the specimen number and step1
    and step2 found in specimens.txt file. This redo file can be used in
    view Paleointensity data").

    Format de sortie : "specimen tmin tmax cooling aniso_corr" (une ligne
    par specimen), le format attendu par `ouvrir_openfilepint_dialog`/
    `openfilepint` (visi_Paleoint.f) - PAS le format ligne/plan de
    `magic_results_to_redo_lines` (colonnes bestfit/ancr/comp n'ont pas de
    sens pour une interpretation de paleointensite). `cooling`/
    `aniso_corr` (colonnes MagIC `int_corr_cooling_rate`/`int_corr_aniso`
    - facteurs multiplicatifs SANS dimension deja calcules par le
    logiciel source de la contribution, cf. le modele de donnees MagIC :
    "Cooling Rate/Anisotropy Correction Factor for intensity") - demande
    explicite utilisateur ("while writing the redo file for
    paleointensity, add the int_corr_cooling_rate, and int_corr_aniso in
    the redo file") - `0` (sentinelle "absent", jamais une vraie valeur
    physique pour un facteur multiplicatif sur une intensite) quand la
    colonne MagIC est vide/absente, comme le fait deja le lecteur
    existant pour son 4e champ optionnel de taux de refroidissement.
    `ouvrir_openfilepint_dialog` applique `cooling` inconditionnellement
    (comme avant) mais `aniso_corr` UNIQUEMENT en repli, si STARpaleomag_Py n'a
    PAS de tenseur .pmagani natif pour ce specimen - demande explicite
    utilisateur, pour eviter de corriger deux fois le meme effet
    (tenseur natif ET facteur MagIC deja applique par le logiciel
    source).

    Selection des lignes : method_codes contient "LP-PI-TRM" (protocole
    Thellier/IZZI thermique explicite) OU le tag generique "LP-PI" seul
    (protocole paleointensite non davantage precise) ET result_quality==
    'g'. BUG REEL corrige (signale par l'utilisateur : "during import of
    Magic data with paleointensity, it is possible to create a redo
    file? can find it" - le fichier n'etait pas cree du tout) : la
    premiere version n'acceptait QUE "LP-PI-TRM" litteralement, alors
    qu'une contribution reelle (Miriam_Magic/magic_contribution_20536.txt,
    276 lignes result_quality='g' exploitables) tague ses determinations
    Thellier classiques (pTRM check + fit anchore - method_codes
    "LP-NO:LP-PI:LP-PI-ALT-PTRM:DE-BFL-A") avec le tag generique "LP-PI"
    SANS jamais ecrire "LP-PI-TRM" explicitement - toutes silencieusement
    ecartees (0 ligne, aucun fichier ecrit, aucun avertissement). La
    presence de "LP-PI-ALT-PTRM" (pTRM check, un concept qui n'existe que
    pour un protocole TRM/thermique - voir vocabulaire MagIC) confirme
    qu'un "LP-PI" nu designe ici bien un Thellier/IZZI thermique classique,
    rejouable par le natif STARpaleomag_Py (paleointensity.py/afficher_arai) au
    meme titre qu'un "LP-PI-TRM" explicite. Exclut toujours les variantes
    clairement NON rejouables de cette facon (`_PI_NON_TRM_CODES`
    ci-dessous - AF/micro-ondes/multi-specimen/proxy alternatif) meme si
    elles portent aussi le tag "LP-PI", et exclut deliberement les lignes
    "IE-BICEP" (methode Bootstrapped Common Enclosing Prism, calcul
    multi-specimens sans borne de temperature rejouable individuellement
    de la meme facon, verifie sur une vraie contribution -
    magic_contribution_19491.txt - ou ces lignes existent en DOUBLON du
    meme int_abs/meas_step_min/max qu'une ligne LP-PI-TRM sœur, mais ne
    representent pas un second ajustement independant a rejouer).
    Deduplique par specimen (une determination retenue par specimen).
    Retourne None si la table specimens est introuvable/vide.

    Specimen SANS meas_step_min/meas_step_max exploitable dans
    specimens.txt -> AUCUNE ligne redo pour lui, jamais de repli sur
    measurements.txt - demande explicite utilisateur : une premiere
    version comblait ce cas par l'etendue min/max de TOUTES les mesures
    Thellier/IZZI du specimen, mais "les bornes des temperatures de
    mesures ne servent pas vraiment" (l'etendue complete des mesures
    n'a rien a voir avec la fenetre reellement retenue pour l'ajustement
    - un point de depart trompeur, pas honnete, pour le fichier redo)."""
    if combined:
        tables = split_combined_magic_file(specimens_source)
        df = tables.get("specimens")
    else:
        df = read_magic_file(specimens_source)

    if df is None or df.empty:
        print("❌ Specimens table not found or empty.")
        return None

    # Les facteurs de correction ne sont PAS forcement portes par la MEME
    # ligne specimens.txt que celle qui fournit smin/smax - BUG REEL
    # corrige (signale par l'utilisateur : "did not read the two
    # parameters in the file and wrote 0.0"), confirme sur une vraie
    # contribution (Miriam_Magic/magic_contribution_20536.txt, specimen
    # 15ZA0701B) : la ligne "g"/"LP-NO:LP-PI:LP-PI-ALT-PTRM:DE-BFL-A" (le
    # fit de direction, qui fournit smin/smax) a int_corr_cooling_rate/
    # int_corr_aniso VIDES, alors qu'une ligne SOEUR du MEME specimen
    # (`int_corr`='c', method_codes SANS le tag de fit "DE-BFL-A",
    # result_quality VIDE - donc jamais selectionnee par le filtre
    # quality=='g' ci-dessous) porte les valeurs reelles (0.893/0.836).
    # Repere ici INDEPENDAMMENT du filtre quality/method_codes utilise
    # pour smin/smax, sur TOUTE ligne du specimen qui porte l'une des
    # deux valeurs - 138/138 specimens de ce fichier reel n'ont jamais
    # plus d'une telle ligne (verifie).
    corr_by_specimen: Dict[str, Tuple[float, float]] = {}
    for _, row_series in df.iterrows():
        row = row_series.to_dict()
        cool_raw = clean_str(row.get("int_corr_cooling_rate", ""))
        ani_raw = clean_str(row.get("int_corr_aniso", ""))
        if not cool_raw and not ani_raw:
            continue
        specimen = clean_str(row.get("specimen", ""))[:12]
        if not specimen:
            continue
        corr_by_specimen[specimen] = (
            parse_float_val(cool_raw, 0.0), parse_float_val(ani_raw, 0.0))

    # `result_quality` est une colonne OPTIONNELLE du data model MagIC -
    # BUG REEL corrige (signale par l'utilisateur : "dans l'importation
    # de certaines contributions avec des donnees de paleointensite, il
    # n'est pas propose de redo") : quand la colonne est ENTIEREMENT
    # ABSENTE de specimens.txt (verifie sur une vraie contribution,
    # magic_contribution_20545.txt - 120 lignes LP-PI-TRM, colonne
    # result_quality absente), `row.get("result_quality", "")` renvoyait
    # silencieusement "" pour CHAQUE ligne, donc `!= "g"` rejetait TOUT -
    # 0 ligne, aucun fichier redo, sans qu'aucune ligne n'ait ete
    # explicitement marquee mauvaise. Distinct du cas ou la colonne
    # EXISTE mais une ligne particuliere est vide/differente de "g" (ce
    # filtre reste applique tel quel dans ce cas - une contribution qui
    # renseigne bien result_quality l'a fait pour distinguer bon/mauvais,
    # a respecter).
    has_quality_col = "result_quality" in df.columns
    seen = set()
    lines = []
    for _, row_series in df.iterrows():
        row = row_series.to_dict()
        if has_quality_col and clean_str(row.get("result_quality", "")).lower() != "g":
            continue
        codes = {c.strip() for c in clean_str(row.get("method_codes", "")).split(":") if c.strip()}
        if not _is_accepted_pi_codes(codes):
            continue
        smin, smax = _convert_step_range(row)
        if smin is None:
            continue
        specimen = clean_str(row.get("specimen", ""))[:12]
        if not specimen or specimen in seen:
            continue
        seen.add(specimen)
        cooling, aniso_corr = corr_by_specimen.get(specimen, (0.0, 0.0))
        lines.append(f"{specimen} {smin} {smax} {cooling} {aniso_corr}")
    return lines


def import_specimens_results_to_redo(specimens_source, output_path, combined=False):
    """Lit specimens.txt (ou la table "specimens" d'un fichier de
    contribution combiné si `combined=True`) et écrit un fichier "redo"
    (`specimen bestfit ancr step_min step_max comp`, voir
    calcul.fit_from_redo_file) pour chaque résultat result_quality=='g'
    dont method_codes désigne un ajustement de ligne/plan. Une seule ligne
    par (specimen, ajustement) même si le même résultat apparaît plusieurs
    fois pour différentes orientations (DA-DIR/DA-DIR-GEO/DA-DIR-TILT) -
    le fichier redo est indépendant de l'orientation."""
    lines = magic_results_to_redo_lines(specimens_source, combined=combined)
    if lines is None:
        return None

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + ("\n" if lines else ""))

    print(f"✅ Redo file written: {output_path} ({len(lines)} line(s))")
    return output_path


# Protocole d'anisotropie (method_codes) -> code2 .pmagani - demande
# explicite utilisateur ("lors de l'importation des fichiers Magic,
# archiver dans le fichier pmagani les donnees de tenseurs d'anisotropie
# avec le code A0 si c'est un tenseur de TRM, N0 pour l'AMS, AA pour
# l'ARM"). 'AA' est DELIBEREMENT distinct du 'F0' utilise par le calcul
# natif STARpaleomag_Py/l'export AMS_Py pour un tenseur d'ARM - choix confirme
# par l'utilisateur (question posee : un tenseur ARM deja calcule par une
# contribution MagIC externe, potentiellement par un tout autre logiciel/
# laboratoire, reste distingable d'un F0 recalcule nativement a partir
# des 6 mesures ; Inverse_ANI_correction... est etendu avec ce meme code
# 'AA' comme option supplementaire - voir calcul._ANI_CODE2).
_ANISO_METHOD_TO_CODE2 = {
    "LP-AN-TRM": "A0",
    "LP-AN-ARM": "AA",
    "LP-AN-MS": "N0",
}


def magic_anisotropy_rows(specimens_source, combined=False) -> List[dict]:
    """Lit specimens.txt (ou la table "specimens" d'une contribution
    combinee si `combined=True`) et retourne un tenseur d'anisotropie
    DEJA CALCULE par cette contribution (colonne `aniso_s`, convention
    MagIC trace-normalisee "s1:s2:s3:s4:s5:s6" = k11:k22:k33:k12:k23:k13,
    memes conventions que celles utilisees par AMS_Py pour l'export
    inverse - .pmagani -> MagIC) pour chaque ligne qui en porte une -
    demande explicite utilisateur ("lors de l'importation des fichiers
    Magic, archiver dans le fichier pmagani les donnees de tenseurs
    d'anisotropie").

    Beaucoup de lignes specimens.txt ne portent PAS d'`aniso_s` du tout
    (resultats de direction/paleointensite plutot que d'anisotropie) -
    ecartees silencieusement, ce n'est pas une anomalie. `code2` est
    assigne depuis `method_codes` (voir _ANISO_METHOD_TO_CODE2) plutot
    que depuis `aniso_type` (present dans certaines contributions mais
    pas garanti par le modele de donnees) - une ligne dont le protocole
    d'anisotropie n'est pas reconnu (aucun des trois codes method_codes
    ci-dessus) est ECARTEE plutot que devinee.

    Chaque element retourne est un dict avec les cles: specimen, code2,
    k11, k22, k33, k12, k23, k13 (floats, trace-normalises - PAS en
    unites physiques absolues, comme le sont deja les 15 variantes du
    calcul natif STARpaleomag_Py : la correction d'inversion ne se sert que des
    RATIOS entre composantes, l'echelle absolue n'entre pas en jeu),
    n_positions, sigma, ftest, ftest12, ftest23 (float ou None si la
    colonne MagIC correspondante est absente/vide) et quality ('g'/'b'
    ou None). Retourne [] si specimens.txt est introuvable/vide ou ne
    contient aucune ligne d'anisotropie exploitable."""
    if combined:
        tables = split_combined_magic_file(specimens_source)
        df = tables.get("specimens")
    else:
        df = read_magic_file(specimens_source)
    if df is None or df.empty:
        return []

    rows: List[dict] = []
    for _, row_series in df.iterrows():
        row = row_series.to_dict()
        aniso_s = clean_str(row.get("aniso_s", ""))
        if not aniso_s:
            continue
        parts = aniso_s.split(":")
        if len(parts) != 6:
            continue
        try:
            k11, k22, k33, k12, k23, k13 = (float(p) for p in parts)
        except ValueError:
            continue

        codes = {c.strip() for c in clean_str(row.get("method_codes", "")).split(":") if c.strip()}
        code2 = next((v for k, v in _ANISO_METHOD_TO_CODE2.items() if k in codes), None)
        if code2 is None:
            continue

        specimen = clean_str(row.get("specimen", ""))[:12]
        if not specimen:
            continue

        n_positions = None
        raw_n = clean_str(row.get("aniso_n_pos", ""))
        if raw_n:
            try:
                n_positions = int(float(raw_n))
            except ValueError:
                n_positions = None

        quality = clean_str(row.get("aniso_ftest_quality", "")).lower() or None
        if quality not in ("g", "b"):
            quality = None

        rows.append({
            "specimen": specimen, "code2": code2,
            "k11": k11, "k22": k22, "k33": k33, "k12": k12, "k23": k23, "k13": k13,
            "n_positions": n_positions,
            "sigma": parse_float_val(row.get("aniso_s_sigma", ""), None) if row.get("aniso_s_sigma") else None,
            "ftest": parse_float_val(row.get("aniso_ftest", ""), None) if row.get("aniso_ftest") else None,
            "ftest12": parse_float_val(row.get("aniso_ftest12", ""), None) if row.get("aniso_ftest12") else None,
            "ftest23": parse_float_val(row.get("aniso_ftest23", ""), None) if row.get("aniso_ftest23") else None,
            "quality": quality,
        })
    return rows


# ---------------------------------------------------------------------------
# Moyennes de site DEJA calculees (sites.txt, colonnes dir_*) -> "mean:"
# dans .pmagres - demande explicite utilisateur ("in Magic file site.txt,
# the site mean direction is sometimes provided, is it possible to archive
# these results in pmagres").
# ---------------------------------------------------------------------------

# Convention MagIC dir_tilt_correction (data_model3, "Percentage tilt
# correction applied... geographic (0%) [to] stratigraphic (100%)",
# heritee de pmag_sites.site_tilt_correction en DM2.5 - -1 y designait le
# repere specimen/echantillon, non corrige) -> code d'orientation STARpaleomag_Py
# (par3_mean, meme convention que self.orientation : 1=echantillon,
# 2=in-situ, 3=apres pendage).
#
# "-1" (coordonnees specimen/echantillon) est volontairement ABSENT de
# cette table, tout comme toute valeur manquante ou non reconnue - demande
# explicite utilisateur ("eviter les moyennes en coordonnees echantillon") :
# une moyenne de Fisher sur plusieurs specimens n'a de sens que si tous les
# specimens ont deja ete ramenes a un repere COMMUN (au moins geographique),
# chaque specimen ayant son propre repere de carotte/echantillon, mutuellement
# non alignes - une moyenne "en coordonnees echantillon" ne peut donc jamais
# etre une moyenne de plusieurs specimens. Une ligne sites.txt dont
# dir_tilt_correction vaut "-1" ou est absent/non reconnu est ECARTEE par
# `magic_site_means` (voir plus bas), jamais reetiquetee sous un autre code
# d'orientation par defaut.
_TILT_CORRECTION_TO_ORIENTATION = {"0": 2.0, "100": 3.0}


def _resolve_tilt_orientation(tilt: str) -> Optional[float]:
    """Resout dir_tilt_correction en code d'orientation STARpaleomag_Py par
    comparaison NUMERIQUE contre _TILT_CORRECTION_TO_ORIENTATION, pas par
    correspondance de chaine exacte - le format varie d'une contribution a
    l'autre ("0"/"100" dans certaines, "0.0"/"100.0" dans d'autres, ex.
    magic_contribution_20340.txt) et une correspondance de chaine ratait
    silencieusement "0.0", faisant passer TOUTES les moyennes de site
    d'une telle contribution pour "coordonnees echantillon" et les
    ecartant a tort (bug reel signale par l'utilisateur : "there is no
    site mean data in the pmagres file during import magic" sur la
    contribution 20340, ou dir_tilt_correction vaut "0.0" pour 200/208
    sites). Retourne None si absent/non reconnu (dont "-1"/"-1.0",
    coordonnees specimen explicites)."""
    try:
        value = float(tilt)
    except ValueError:
        return None
    for key, orientation in _TILT_CORRECTION_TO_ORIENTATION.items():
        if abs(value - float(key)) < 1e-6:
            return orientation
    return None


def magic_site_means(sites_source, combined=True) -> Optional[List[dict]]:
    """Lit sites.txt (ou la table "sites" d'une contribution combinee) et
    retourne un dict par (site, orientation) portant une direction
    moyenne DEJA calculee (dir_dec/dir_inc renseignes) - une contribution
    fournit souvent 2 lignes par site (geographique ET apres pendage,
    dir_tilt_correction='0'/'100'). Chaque dict : site, orientation (2/3,
    voir _TILT_CORRECTION_TO_ORIENTATION), dec, inc, alpha95, k, n,
    specimens (liste de noms, depuis la colonne "specimens" - "samples" en
    repli si absente), lat, lon (coordonnees du SITE, pas du specimen),
    vgp_lat, vgp_lon, component (dir_comp_name, "A" si absent - voir
    calcul.FitResult.component). Une ligne dont dir_tilt_correction vaut "-1" ou est
    absent/non reconnu est ECARTEE (jamais rangee par defaut en
    1/echantillon - voir le commentaire pres de
    _TILT_CORRECTION_TO_ORIENTATION), avec un compte-rendu imprime.
    Retourne None si la table sites est introuvable/vide (distinct d'une
    liste vide = aucune moyenne directionnelle presente/exploitable)."""
    if combined:
        tables = split_combined_magic_file(sites_source)
        df = tables.get("sites")
    else:
        df = read_magic_file(sites_source)
    if df is None or df.empty:
        print("❌ Sites table not found or empty.")
        return None
    if "dir_dec" not in df.columns or "dir_inc" not in df.columns:
        return []

    rows = []
    n_skipped_sample_coords = 0
    for _, row_series in df.iterrows():
        row = row_series.to_dict()
        dec = parse_float_val(row.get("dir_dec", ""), None)
        inc = parse_float_val(row.get("dir_inc", ""), None)
        if dec is None or inc is None:
            continue
        site = clean_str(row.get("site", ""))
        if not site:
            continue
        tilt = clean_str(row.get("dir_tilt_correction", ""))
        orientation = _resolve_tilt_orientation(tilt)
        if orientation is None:
            n_skipped_sample_coords += 1
            continue

        specimens_str = clean_str(row.get("specimens", "")) or clean_str(row.get("samples", ""))
        specimens = [s for s in specimens_str.split(":") if s]

        n = parse_float_val(row.get("dir_n_specimens", ""), None)
        if n is None:
            n = parse_float_val(row.get("dir_n_samples", ""), len(specimens))

        rows.append({
            "site": site,
            "orientation": orientation,
            "dec": dec,
            "inc": inc,
            "alpha95": parse_float_val(row.get("dir_alpha95", ""), 0.0),
            "k": parse_float_val(row.get("dir_k", ""), 0.0),
            "n": int(round(n)),
            "specimens": specimens,
            "lat": parse_float_val(row.get("lat", ""), 0.0),
            "lon": parse_float_val(row.get("lon", ""), 0.0),
            "vgp_lat": parse_float_val(row.get("vgp_lat", ""), 0.0),
            "vgp_lon": parse_float_val(row.get("vgp_lon", ""), 0.0),
            # Ovale de confiance du VGP - present directement dans CERTAINES
            # contributions reelles (verifie : magic_contribution_20340.txt)
            # ; 0.0 si absent, calcule alors depuis alpha95/inc au moment de
            # l'archivage (voir convert_magic_to_r._convert_magic_results,
            # calcul.dp_dm_from_a95) - demande explicite utilisateur.
            "vgp_dp": parse_float_val(row.get("vgp_dp", ""), 0.0),
            "vgp_dm": parse_float_val(row.get("vgp_dm", ""), 0.0),
            # Etiquette de composante de magnetisation (voir
            # calcul.FitResult.component) - "Site direction component
            # name" dans le modele MagIC (verifie live sur
            # earthref.org/MagIC/data-models/3.0.json), permet a une
            # contribution reelle de distinguer plusieurs moyennes du meme
            # site (ex. "A"/"B") des l'import. "A" si absent (ancien
            # champ, valeur par defaut de FitResult.component).
            "component": clean_str(row.get("dir_comp_name", "")) or "A",
        })
    if n_skipped_sample_coords:
        print(f"⚠️  {n_skipped_sample_coords} site mean(s) skipped: dir_tilt_correction "
              f"missing or '-1' (specimen coordinates) - a multi-specimen mean is not "
              f"meaningful in sample coordinates.")
    return rows


if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)

    data_dir = filedialog.askdirectory(
        title="Select the folder containing the 4 MagIC files"
    )

    root.destroy()

    if not data_dir:
        print("❌ No selection made. Operation cancelled.")
        sys.exit(0)

    sites_path = os.path.join(data_dir, "sites.txt")
    samples_path = os.path.join(data_dir, "samples.txt")
    specimens_path = os.path.join(data_dir, "specimens.txt")
    measurements_path = os.path.join(data_dir, "measurements.txt")
    output_path = os.path.join(data_dir, "extracted_data_custom.ren")

    format_custom_output(
        sites_file=sites_path,
        samples_file=samples_path,
        specimens_file=specimens_path,
        measurements_file=measurements_path,
        output_file=output_path,
        default_volume=9.99,
    )
