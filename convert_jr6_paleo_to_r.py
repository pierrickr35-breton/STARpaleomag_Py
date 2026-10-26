"""
Convertit un fichier AGICO JR6 (".jr6", colonnes fixes) contenant une
experience de PALEOINTENSITE (Thellier-Thellier) vers .prmag - demande
explicite utilisateur ("ajouter l'importation de fichiers JR6 pour
specifiquement de la paleointensite. Les etapes de temperature sont
souvent par groupe : temp (etape en Z+, code R), temp+1 (en Z-, code V),
temp+2 (etape de PTRM test). Demander des l'ouverture la valeur du champ
applique et demander si elle s'applique a tous les echantillons. Creer un
fichier .prmag").

Distinct de import_new_data.parse_jr6_file (archivage de nouvelles mesures
de DESAIMANTATION dans un .prmag existant) : ici le fichier cree le .prmag
et les paliers TDnnn sont classes en N/R/V/P.

Colonnes du fichier JR6 (largeurs fixes - mappage confirme par
pmagpy.convert_2_magic.jr6_jr6, `char_count`, et verifie sur Isla Mala.jr6) :

    [0:10]  specimen          [40:44] azimuth        [64:67] param1 \
    [10:18] etape (NRM/TDnnn) [44:48] dip            [67:70] param2  } schema
    [18:24] x                 [48:52] foliation az   [70:73] param3  } Agico
    [24:30] y                 [52:56] foliation dip  [73:76] param4 /
    [30:36] z                 [56:64] lineation      [76:80] csd
    [36:40] exposant

Conventions reprises du pipeline JR6 historique (ImportJR6data.f95,
`importinoldren`) et de import_new_data.parse_jr6_file - PAS reinventees :

  - moment x = valeur * 10**(exposant - 5), avec volume 10.00 cm3 ecrit
    dans l'en-tete : MAG(A/m) = x / 10 cm3 retombe exactement sur
    valeur * 10**exposant, l'aimantation (A/m) rapportee par le JR6.
    (Verifie sur 25AN_RGonzalez : "5.49 ... -2" -> 5.490E-07 dans le .ren
    produit par le Fortran.) La paleointensite est un RAPPORT de moments
    (pente d'Arai) et ne depend pas de ce choix de volume.
  - azimut/pendage du carottage = colonnes azimuth/dip tels quels
    (`az:`/`in:` du Fortran) -> `azimuth:`/`dip:` du .prmag.
  - instrument "J6".

Classification des paliers TDnnn (ordre du fichier, PAR SPECIMEN - le JR6
entrelace les specimens d'un meme tour de mesure) :

  - NRM                 -> N, cod2 "0", etape 0
  - T            (neuf) -> R (aimantation dans +Z), etape T, nouveau groupe
  - T+1 apres R(T)      -> V (champ inverse -Z), meme etape T et meme lettre
  - T+2 apres R(T)      -> P (controle pTRM de T)
  - T identique a un R anterieur -> P (controle pTRM de T) : l'operateur a
    nomme le controle "TD540" au lieu de "TD542" dans Isla Mala.jr6.

cod2 = lettre de groupe A, B, C... incrementee a chaque nouveau R ; V porte
la lettre de son R, P celle du groupe COURANT (le groupe qui vient d'etre
mesure) et NON de la temperature qu'il re-teste - meme convention que les
.prmag Thellier existants (SanJuan_Pmag.prmag : "150 P C" apres le groupe
C = 250). `etape` d'un P = temperature re-testee (voir
paleointensity.py : slot courant pour le trace, etape pour le calcul).

Champ labo : ecrit en prefixe numerique du `comment:` du specimen (meme
convention que SanJuan_Pmag.prmag, `comment: 50`) - c'est de la que
magic_export._measurement_treatment tire treat_dc_lowfield /
phi / theta de R, V et P a l'ecriture.

Site : "n.d" (le JR6 n'en porte pas, comme Montpellier/Utrecht), a raffiner
ensuite via la complete_sample_info. Seul le protocole R/V/P est traite
(pas IZZI/S). Les paliers non reconnus sont ECARTES et journalises, jamais
devines.
"""

import argparse
import os
import re
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple, Union

from testlect import Measurement, Pmag
from convert_ren_to_r import _measurement_rows
from magic_export import _specimen_method_codes

FORMAT_HEADER = (
    "#STARpaleomag_Py .prmag v1  angles=deg  fields in milliTesla (mT) for strong "
    "fields AF or IRM and in microTesla (uT) for low field paleointensity "
    "or ARM  temperatures in degC  date=ISO8601"
)
S_FIELD_HEADER = "#s = magnetic susceptibility in 1e-5 SI"

# Volume nominal du pipeline JR6 historique (v:10.00) - voir docstring.
JR6_NOMINAL_VOLUME_CM3 = 10.0

# (debut, fin) des champs a largeur fixe - voir docstring du module.
_COLS = {
    "specimen": (0, 10), "step": (10, 18), "x": (18, 24), "y": (24, 30),
    "z": (30, 36), "expon": (36, 40), "azimuth": (40, 44), "dip": (44, 48),
    "fol_az": (48, 52), "fol_dip": (52, 56),
    "p1": (64, 67), "p2": (67, 70), "p3": (70, 73), "p4": (73, 76),
}

_TD_RE = re.compile(r"^TD\s*(\d+(?:\.\d+)?)$")
_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _field(line: str, name: str) -> str:
    a, b = _COLS[name]
    return line[a:b].strip()


def _num(line: str, name: str, default: float = 0.0) -> float:
    s = _field(line, name)
    try:
        return float(s)
    except ValueError:
        return default


def read_jr6_paleo_file(path: str, encoding: str = "latin-1") -> "OrderedDict[str, dict]":
    """Parse un .jr6 -> OrderedDict specimen -> {caz, cin, bed_strike,
    bed_dip, scheme, steps:[(label, x, y, z, expon)]} (ordre de premiere
    apparition preserve, mesures dans l'ordre du fichier). Les lignes
    illisibles sont ignorees (comptees dans `skipped`)."""
    specimens: "OrderedDict[str, dict]" = OrderedDict()
    with open(path, "r", encoding=encoding, errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if len(line) < 40 or not line[0:10].strip():
                continue
            try:
                x, y, z = _num(line, "x"), _num(line, "y"), _num(line, "z")
                expon = int(_field(line, "expon"))
            except ValueError:
                continue
            name = _field(line, "specimen").upper().replace(" ", "_")
            sp = specimens.get(name)
            if sp is None:
                fol_dip = _num(line, "fol_dip")
                # bed_dip_direction = bed_dip_strike + 90 (voir
                # convert_ren_to_r.DIFF_WITH_MAGIC) ; sans pendage de
                # foliation, azimut de pendage 0 = "pas de litage".
                strike = (_num(line, "fol_az") - 90.0) % 360.0 if fol_dip else 0.0
                p = [_field(line, k) or "0" for k in ("p1", "p2", "p3", "p4")]
                sp = specimens[name] = {
                    "specimen": name,
                    "caz": _num(line, "azimuth"),
                    "cin": _num(line, "dip"),
                    "bed_strike": strike,
                    "bed_dip": fol_dip,
                    "scheme": "A" + "_".join(p),
                    "steps": [],
                }
            sp["steps"].append((_field(line, "step").upper(), x, y, z, expon))
    return specimens


def classify_steps(steps: List[tuple]) -> Tuple[List[dict], List[str]]:
    """Classe les paliers d'UN specimen (ordre du fichier) en N/R/V/P - voir
    docstring du module. Retourne (paliers classes, avertissements). Chaque
    palier classe : {cod1, cod2, etape, x, y, z, label} ; x/y/z en Am2."""
    out: List[dict] = []
    warnings: List[str] = []
    r_temps: List[float] = []      # temperatures des R deja vus
    last_r: Optional[float] = None
    v_pending = False              # dernier R sans V pour l'instant
    group = -1

    def letter() -> str:
        return _LETTERS[group % len(_LETTERS)]

    for label, x, y, z, expon in steps:
        scale = 10.0 ** (expon - 5)
        base = {"x": x * scale, "y": y * scale, "z": z * scale, "label": label}

        if label.startswith("NRM"):
            out.append({**base, "cod1": "N", "cod2": "0", "etape": 0.0})
            continue
        m = _TD_RE.match(label)
        if not m:
            warnings.append(f"step {label!r} not recognised (only NRM and TDnnn) - skipped")
            continue
        t = float(m.group(1))

        if last_r is not None and v_pending and t == last_r + 1:
            out.append({**base, "cod1": "V", "cod2": letter(), "etape": last_r})
            v_pending = False
        elif any(t == r + 2 for r in r_temps):
            target = next(r for r in reversed(r_temps) if t == r + 2)
            out.append({**base, "cod1": "P", "cod2": letter(), "etape": target})
        elif t in r_temps:
            out.append({**base, "cod1": "P", "cod2": letter(), "etape": t})
        else:
            if v_pending:
                warnings.append(f"R at {last_r:g} has no V step (next step {label})")
            group += 1
            out.append({**base, "cod1": "R", "cod2": letter(), "etape": t})
            r_temps.append(t)
            last_r, v_pending = t, True

    if v_pending:
        warnings.append(f"R at {last_r:g} has no V step (end of specimen)")
    return out, warnings


def _fmt(x: float) -> str:
    return f"{x:g}"


def _header_block(sp: dict, p: Pmag, field_uT: float) -> str:
    specimen = sp["specimen"]
    sample = specimen[:-1] if specimen[-1:].isalpha() else specimen
    line_a = (
        f"specimen: {specimen}\tsample: {sample}\tsite: n.d\t"
        f"volume: {JR6_NOMINAL_VOLUME_CM3:.2f}\tmass: n.d\t"
        f"lat: n.d\tlon: n.d\televation: n.d\t"
        f"stratigraphic_height: n.d\t"
        f"comment: {_fmt(field_uT)}"
    )
    line_b = (
        f"azimuth: {sp['caz']:.1f}\tdip: {sp['cin']:.1f}\t"
        f"date: n.d\tmagnetic_azimuth: n.d\tsolar_azimuth: n.d\t"
        f"orient_tool: {sp['scheme']}"
    )
    line_c = f"bed_dip_strike: {sp['bed_strike']:.1f}\tbed_dip: {sp['bed_dip']:.1f}"
    line_d = (
        "formation: n.d\tage: n.d\tage_low: n.d\tage_high: n.d\tage_unit: n.d\t"
        "geologic_classes: n.d\tgeologic_types: n.d\t"
        "lithologies: n.d\tlocation: n.d\t"
        f"obs: n.d\tmethod_codes: {_specimen_method_codes(p)}"
    )
    return "\n".join([line_a, line_b, line_c, line_d])


def _pmag_for(sp: dict, steps: List[dict], field_uT: float) -> Pmag:
    mesures = [
        Measurement(etape=s["etape"], cod1=s["cod1"], cod2=s["cod2"],
                    x=s["x"], y=s["y"], z=s["z"], q=0, ins="J6", s=0.0)
        for s in steps
    ]
    return Pmag(
        id=sp["specimen"], cin=sp["cin"], caz=sp["caz"], dip=sp["bed_dip"],
        str_=sp["bed_strike"], norme="v", vol=JR6_NOMINAL_VOLUME_CM3,
        com=_fmt(field_uT), outilorient=sp["scheme"], mesures=mesures,
    )


def convert_file(
    path_in: str, path_out: str,
    field_uT: Union[float, Dict[str, float]],
) -> Tuple[int, List[str]]:
    """Convertit un .jr6 de paleointensite -> .prmag. `field_uT` = un
    nombre (meme champ pour tous les specimens) ou un dict specimen ->
    champ (uT). Retourne (nb specimens ecrits, lignes de rapport)."""
    specimens = read_jr6_paleo_file(path_in)
    if not specimens:
        raise ValueError("no JR6 measurement lines found in this file")

    def field_of(name: str) -> float:
        value = field_uT.get(name) if isinstance(field_uT, dict) else field_uT
        if value is None or value <= 0:
            raise ValueError(f"no valid applied field (uT) for specimen {name}")
        return float(value)

    blocks: List[str] = []
    report: List[str] = []
    for name, sp in specimens.items():
        steps, warnings = classify_steps(sp["steps"])
        if not steps:
            report.append(f"{name}: no usable step - skipped")
            continue
        field = field_of(name)
        p = _pmag_for(sp, steps, field)
        counts = {c: sum(1 for s in steps if s["cod1"] == c) for c in "NRVP"}
        report.append(
            f"{name}: field {_fmt(field)} uT, {counts['R']} R / {counts['V']} V / "
            f"{counts['P']} P (+{counts['N']} NRM)"
            + "".join(f"\n    WARNING: {w}" for w in warnings)
        )
        blocks.append(_header_block(sp, p, field) + "\n" + "\n".join(_measurement_rows(p)))

    if not blocks:
        raise ValueError(
            "no usable NRM/TDnnn step found - is this a thermal (paleointensity) JR6 file?")

    lines = [
        FORMAT_HEADER,
        f"#converted from {os.path.basename(path_in)} by convert_jr6_paleo_to_r.py "
        f"(AGICO JR6 paleointensity R/V/P -> .prmag, nominal volume="
        f"{JR6_NOMINAL_VOLUME_CM3:.2f} cm3)",
        S_FIELD_HEADER,
        "",
    ]
    with open(path_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        f.write("\n\n".join(blocks) + "\n")
    return len(blocks), report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a JR6 paleointensity file to .prmag")
    parser.add_argument("jr6_file")
    parser.add_argument("--field", type=float, required=True, help="Applied lab field (uT)")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    out = args.output or os.path.splitext(args.jr6_file)[0] + ".prmag"
    nb, rep = convert_file(args.jr6_file, out, args.field)
    print(f"{nb} specimen(s) written to {out}")
    print("\n".join(rep))
