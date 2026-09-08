"""
Convertit un fichier .ren (format Rennes historique, voir testlect.read_ren_file)
vers le nouveau format .r discute avec l'utilisateur (cf. exemple.txt) :

- le "code" de mesure historique (ex. "100RA", "50F+") a d'abord ete
  conserve tel quel en une seule colonne (decision initiale : "keep the
  codes from the legacy Rennes data... but add also their meaning in
  magic format"), PUIS retire une fois les colonnes step/cod1/cod2
  ajoutees separement (devenu totalement redondant - demande explicite
  utilisateur).
- les champs MagIC (method_codes, instrument_codes, treat_temp,
  treat_ac_field, treat_dc_field, treat_dc_field_phi, treat_dc_field_theta)
  sont calcules UNE FOIS ici, a la conversion, en reutilisant les tables
  DEJA EPROUVEES de magic_export.py (_measurement_treatment, _instrument)
  plutot que d'etre re-derives a chaque export - objectif explicite de
  l'utilisateur ("make the new format more interoperable with Magic").
- bed_dip_strike est stocke DIRECTEMENT (valeur native p.str_, telle que
  mesuree sur le terrain, convention Rennes) - decision finale de
  l'utilisateur (revient sur un choix intermediaire de cette session qui
  stockait bed_dip_direction ; voir DIFF_WITH_MAGIC pour la conversion
  vers MagIC : bed_dip_direction = bed_dip_strike + 90, meme convention
  etablie/verifiee cette session, voir memoire
  project_stereoutils_py_port.md / project_starmac_fortran_gotchas.md).

2 bugs reperes en relisant magic_export.py._measurement_treatment pour ce
portage sont CORRIGES ici (pas reproduits) :
  1) `ifield = float(ech.com[:2])` ne lisait que les 2 premiers caracteres
     du commentaire (ex. "50_V_R1_1711" -> 50, mais silencieusement faux
     des que le champ fait >=100 ou <10 uT). Remplace par une regex sur
     le prefixe numerique complet.
  2) `treat_dc_field` n'etait jamais renseigne pour les pas X/Y/Z (theta/
     phi calcules mais champ applique perdu, alors que ces pas sont
     precisement les pas d'anisotropie de TRM discutes avec l'utilisateur).
     Corrige ici : meme valeur que R/V (ifield*1e-6).

`quality` (lettre g/b, NOUVEAU champ introduit par l'utilisateur, absent
des .ren existants) est force a "g" par cette conversion (demande
explicite utilisateur : aucune mesure n'etait auparavant marquee
mauvaise, MagIC porte deja cette colonne dans de nombreuses archives -
"g" par defaut a l'import, correction manuelle en "b" ensuite si
l'utilisateur identifie un probleme sur une mesure precise). `age` reste
un champ texte libre (MagIC
fournit "1835 - 1835 AD" deja pre-forme via geologic_age - aucune donnee
source ne permet de le re-decouper fiablement en age_low/age_high/
age_unit distincts a ce stade)."""

import argparse
import os
import re
from typing import List, Optional, Tuple

from testlect import Pmag, read_ren_file
from magic_export import _measurement_treatment, _instrument, _specimen_method_codes
from calcul import results_path_for, convert_legacy_results_file

FORMAT_HEADER = (
    "#STARpaleomag_Py .prmag v1  angles=deg  fields in milliTesla (mT) for strong "
    "fields AF or IRM and in microTesla (uT) for low field paleointensity "
    "or ARM  temperatures in degC  date=ISO8601"
)
DIFF_WITH_MAGIC = (
    "#difference with Magic: azimuth is azimuth of X +90, dip (as measured "
    "with ASC tool) is - dip of sample's x direction from the horizontal; "
    "bed_dip_direction = bed_dip_strike + 90"
)
S_FIELD_HEADER = "#s = magnetic susceptibility in 1e-5 SI"

# cod1 pour lesquels chaque champ treat_* est effectivement calcule par
# _measurement_treatment (voir magic_export.py:479-556) - au dela de ce
# perimetre le champ reste a sa valeur d'initialisation (0.0), qui ne
# distingue pas "non applicable" de "vraiment zero" : on choisit ici de
# n'ecrire une valeur que la ou la fonction source calcule reellement
# quelque chose, "n.d" partout ailleurs.
_AC_FIELD_CODES = {"A", "F"}
_DC_FIELD_CODES = {"A", "I", "R", "V", "P", "X", "Y", "Z", "L", "Q"}
# 'P' (controle pTRM) ajoute a _THETA_CODES - bug corrige, signale par
# l'utilisateur ("the PTRM checks... are done in the field direction of
# the R step. in this case it is not n.d but 90.0") : _measurement_treatment
# renseigne desormais theta=90.0 pour 'P', mais ce set decidait seul de
# l'ecrire ou non dans le fichier - sans 'P' ici, la valeur calculee
# aurait ete silencieusement remplacee par "n.d" a l'ecriture.
_THETA_CODES = {"R", "V", "P", "X", "Y", "Z", "L", "Q"}
_PHI_CODES = {"R", "V", "X", "Y", "Z"}
# thermal demag (N/D/S) ET pas d'in-field paleointensite/anisotropie
# (R/V/P/X/Y/Z) ET vitesse de refroidissement (L/Q) sont tous des pas en
# temperature (`etape` en degC) - demande explicite utilisateur ("for the
# paleointensity and LT-T experiment treat_temp is the step"). Seuls A/F
# (Oersted, AF/ARM) et I (IRM) restent a temperature ambiante (n.d).
_TEMP_CODES = {"N", "D", "S", "R", "V", "P", "X", "Y", "Z", "L", "Q"}

_NUM_RE = re.compile(r"^\s*(-?\d+(?:\.\d+)?)")


def _iso_date(p: Pmag) -> str:
    if not p.year:
        return "n.d"
    return f"{p.year:04d}-{p.month:02d}-{p.day:02d}T{p.hour:02d}:{p.minute:02d}"


def _guess_sample_name(p: Pmag) -> str:
    """Nom d'echantillon (champ MagIC 'sample') : magic_sample si connu
    (import MagIC), sinon heuristique historique Rennes (specimen sans sa
    derniere lettre, ex. '11CL7801A' -> '11CL7801')."""
    if p.magic_sample:
        return p.magic_sample
    if p.id and p.id[-1:].isalpha():
        return p.id[:-1]
    return p.id


def _parse_ifield_uT(com: str) -> Optional[float]:
    """Champ applique (uT) code en tete du commentaire (ex.
    '50_V_R1_1711' -> 50.0). Remplace `float(ech.com[:2])` de
    magic_export.py - voir bug #1 dans le docstring du module."""
    m = _NUM_RE.match(com or "")
    return float(m.group(1)) if m else None


def _nd(value, fmt: Optional[str] = None) -> str:
    """"n.d" quand `value` est absente, sinon `format(value, fmt)`. Si
    `fmt` porte une largeur (ex. "7.1f"), "n.d" est aussi cale a cette
    largeur (justifie a droite) pour garder la colonne alignee meme sur
    les lignes ou le champ ne s'applique pas - demande explicite
    utilisateur ("format... to keep aligned the columns")."""
    if value is None:
        if fmt:
            width_digits = "".join(ch for ch in fmt.split(".")[0] if ch.isdigit())
            if width_digits:
                return "n.d".rjust(int(width_digits))
        return "n.d"
    return format(value, fmt) if fmt is not None else str(value)


def _sample_header_block(p: Pmag) -> str:
    sample = _guess_sample_name(p)
    site = p.magic_site or (p.id[:6] if len(p.id) >= 6 else p.id)
    volume, mass = (None, p.vol) if p.norme == "m" else (p.vol, None)

    # bed_dip_strike : valeur native Rennes (p.str_), aucune conversion -
    # demande explicite utilisateur (change bed_dip_direction -> bed_dip_
    # strike). Pour MagIC (bed_dip_direction) : bed_dip_direction =
    # bed_dip_strike + 90 (voir DIFF_WITH_MAGIC).
    bed_dip_strike = p.str_
    bed_dip = p.dip

    method_summary = _specimen_method_codes(p)

    line_a = (
        f"specimen: {p.id}\tsample: {sample}\tsite: {site}\t"
        f"volume: {_nd(volume, '.2f')}\tmass: {_nd(mass, '.2f')}\t"
        f"lat: {p.lat:.5f}\tlon: {p.rlong:.5f}\televation: {p.altitude:.1f}\t"
        f"stratigraphic_height: {_nd(p.stratigraphic_height, '.2f')}\t"
        f"comment: {p.com.strip() or 'n.d'}"
    )
    line_b = (
        f"azimuth: {p.caz:.1f}\tdip: {p.cin:.1f}\tdate: {_iso_date(p)}\t"
        f"magnetic_azimuth: {_nd(p.azmag or None, '.1f')}\t"
        f"solar_azimuth: {_nd(p.azsun or None, '.1f')}\t"
        f"orient_tool: {p.outilorient or 'n.d'}"
    )
    line_c = f"bed_dip_strike: {bed_dip_strike:.1f}\tbed_dip: {bed_dip:.1f}"
    line_d = (
        f"formation: {p.magic_fm or 'n.d'}\tage: {p.magic_age.strip() or 'n.d'}\t"
        f"geologic_classes: {p.magic_gc or 'n.d'}\tgeologic_types: {p.magic_smt or 'n.d'}\t"
        f"lithologies: {p.magic_li or 'n.d'}\tlocation: {p.magic_loc or 'n.d'}\t"
        f"obs: {p.magic_obs or 'n.d'}\tmethod_codes: {method_summary}"
    )
    return "\n".join([line_a, line_b, line_c, line_d])


_MEAS_HEADER = (
    "step\tcod1\tcod2\tx\ty\tz\terror\tquality\tinstrument\ts\t"
    "treat_temp\ttreat_ac_field\ttreat_dc_strongfield\ttreat_dc_lowfield\t"
    "treat_dc_field_phi\ttreat_dc_field_theta\t"
    "method_codes\tinstrument_codes\ttreat_step_num"
)

# treat_dc_field (MagIC) est scinde en 2 colonnes ici, une par regime
# d'unite (demande explicite utilisateur, resout l'ambiguite "l'unite
# depend de la ligne" signalee precedemment) : treat_dc_strongfield (mT,
# IRM) et treat_dc_lowfield (uT, ARM/paleointensite) - mutuellement
# exclusives par cod1, une seule des deux non-"n.d" par ligne. Recombiner
# en treat_dc_field (Tesla) au moment de l'export MagIC.
_DC_STRONGFIELD_CODES = {"I"}
_DC_LOWFIELD_CODES = _DC_FIELD_CODES - _DC_STRONGFIELD_CODES

def _step_value(m) -> float:
    """`m.etape` est deja la valeur physique reelle (mT pour A/F, degC
    sinon) - testlect.parse_measure_line convertit l'Oersted historique
    (mT*10) DES LA LECTURE d'un .ren, plus a la conversion (ni a une
    echelle dediee) ici - demande explicite utilisateur ("on importing
    legacy files from Rennes, we can divide by 10 from oersted to
    mT... convert all step integer to float")."""
    return float(m.etape)


def _measurement_rows(p: Pmag) -> List[str]:
    ifield = _parse_ifield_uT(p.com)
    ifield_val = ifield if ifield is not None else 0.0
    rows = [_MEAS_HEADER]
    for j, m in enumerate(p.mesures):
        codes, temp, af_field, dc_field, phi, theta = _measurement_treatment(
            m, p.mesures[:j], ifield_val)

        # bug #2 (voir docstring module) : dc_field jamais renseigne pour
        # X/Y/Z dans le Fortran/magic_export.py d'origine - corrige ici.
        if m.cod1 in ("X", "Y", "Z") and ifield is not None:
            dc_field = ifield_val * 1.0e-6

        # Unites lisibles (demande explicite utilisateur) : AF et IRM
        # (champs forts) en mT ; ARM/paleointensite (champs faibles) en
        # uT ; temperature en degC (pas Kelvin - `temp` retourne par
        # _measurement_treatment est deja etape+273, on reprend
        # directement etape). IRM (mT) et ARM/paleointensite (uT)
        # occupent maintenant 2 colonnes distinctes plutot qu'une seule
        # colonne treat_dc_field dont l'unite dependait de la ligne.
        af_field_mT = af_field * 1.0e3 if m.cod1 in _AC_FIELD_CODES else None
        dc_strongfield = dc_field * 1.0e3 if m.cod1 in _DC_STRONGFIELD_CODES else None
        dc_lowfield = dc_field * 1.0e6 if m.cod1 in _DC_LOWFIELD_CODES else 0.0
        # AF/ARM/IRM (A/F/I) : pas de pas thermique, temperature ambiante
        # - demande explicite utilisateur : 0.0 plutot que "n.d".
        temp_c = float(m.etape) if m.cod1 in _TEMP_CODES else 0.0

        step_str = f"{_step_value(m):6.1f}"
        ins_short = (m.ins or "").strip()
        ins_full = _instrument(m.ins)

        # "error" (m.q) - demande explicite utilisateur ("lors de
        # l'importation, lorsque l'instrument est S, mettre n.d lors de
        # l'import et ne pas la prendre en compte lors de
        # l'exportation vers magic") : pour l'instrument "S" (spinner,
        # jeux de donnees de ~30 ans), la definition de ce facteur etait
        # differente de la convention moderne (2G cryo C1/C2) - pas
        # comparable, donc pas importee comme si elle l'etait.
        error_str = "n.d" if ins_short == "S" else str(m.q)

        # Largeurs fixes (demande explicite utilisateur : "format... to
        # keep aligned the columns") - s en f7.1 (equivalent Fortran),
        # treat_* alignes de la meme facon ; method_codes justifie a
        # gauche sur une largeur fixe (texte de longueur variable - les
        # codes les plus longs depasseront quand meme, inevitable pour
        # une colonne texte).
        row = "\t".join([
            step_str, m.cod1, m.cod2,
            f"{m.x:11.3E}", f"{m.y:11.3E}", f"{m.z:11.3E}",
            error_str, "g", ins_short, f"{m.s:7.1f}",
            _nd(temp_c, "6.1f"),
            _nd(af_field_mT, "7.2f"),
            _nd(dc_strongfield, "7.2f"),
            _nd(dc_lowfield, "7.2f"),
            _nd(phi if m.cod1 in _PHI_CODES else None, "6.1f"),
            _nd(theta if m.cod1 in _THETA_CODES else None, "6.1f"),
            f"{codes:<40}", ins_full, str(j + 1),
        ])
        rows.append(row)
    return rows


def convert_sample(p: Pmag) -> str:
    return _sample_header_block(p) + "\n" + "\n".join(_measurement_rows(p))


def convert_file(
    path_in: str, path_out: str,
    legacy_results_path: Optional[str] = None,
) -> Tuple[int, int, List[str]]:
    """Retourne (nb echantillons convertis, nb resultats convertis depuis
    un ANCIEN fichier .r compagnon - 0 si aucun n'existe, avertissements
    de renumerotation - voir calcul.convert_legacy_results_file). Le
    fichier .r (positions de colonnes fixes) est converti vers .pmagres a
    cote de `path_out` si trouve - demande explicite utilisateur ("during
    the convert legacy file to the new format, is it possible to convert
    the .r file with the results").

    `legacy_results_path` est OPTIONNEL : par defaut, meme base que
    `path_in` (.ren -> .r), ce qui suppose que le .r legacy est a cote du
    .ren en cours de conversion - vrai pour un usage direct de "Convert
    .ren to new format .r...". Mais le pipeline unifie "Import legacy
    files..." (convert_legacy_ren.convert_legacy_auto) produit un .ren
    intermediaire nomme "<base>_converted.ren", DIFFERENT de la base du
    fichier legacy d'origine "<base>.txt" a cote duquel se trouve le vrai
    "<base>.r" - sans ce parametre explicite, la recherche automatique
    "<base>_converted.r" ne trouverait jamais ce fichier - demande
    explicite utilisateur ("if there is a .r result file, can it be
    converted to the .pmagres format")."""
    samples = read_ren_file(path_in)
    schemes = {p.outilorient for p in samples if p.outilorient}
    lines = [
        FORMAT_HEADER,
        f"#converted from {os.path.basename(path_in)} by convert_ren_to_r.py",
    ]
    if len(schemes) == 1:
        lines.append(f"#Agico orientation scheme {schemes.pop()}")
    elif len(schemes) > 1:
        lines.append(f"#WARNING: multiple orientation schemes in source file: {sorted(schemes)}")
    lines.append(DIFF_WITH_MAGIC)
    lines.append(S_FIELD_HEADER)
    lines.append("")

    blocks = [convert_sample(p) for p in samples]
    with open(path_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        f.write("\n\n".join(blocks) + "\n")

    old_results_path = legacy_results_path or (os.path.splitext(path_in)[0] + ".r")
    new_results_path = results_path_for(path_out)
    n_results, warnings = convert_legacy_results_file(old_results_path, new_results_path)
    return len(samples), n_results, warnings


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a .ren file to the new .r format")
    parser.add_argument("ren_file")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    out = args.output or os.path.splitext(args.ren_file)[0] + ".prmag"
    n, n_results, warnings = convert_file(args.ren_file, out)
    print(f"{n} sample(s) converted -> {out}")
    if n_results:
        print(f"{n_results} result(s) converted -> {results_path_for(out)}")
    for w in warnings:
        print(f"WARNING: {w}")
