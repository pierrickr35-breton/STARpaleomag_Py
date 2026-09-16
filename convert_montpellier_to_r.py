"""
Convertit le format "Montpellier" (fichier ".montpellier" - demande
explicite utilisateur "convert also the format from Montpellier") vers
.prmag. Table plate auto-descriptive, UNE ligne par mesure :

    Sample Step Intensity Dec Inc Azimuth Plunge Strike Dip
    BN4.1A 20 0.0033090577239462235 143.83182599744745 1.3952186578458248 52 59 182 44
    BN4.2A 0.0 0.013459457344537823 203.44795925963203 -3.9884234888185004 338 44 182 44
    ...

Colonnes VERIFIEES (pas supposees) contre les donnees Utrecht/PMAG2 deja
converties et validees pour CE MEME jeu de donnees (BN4.*, meme
collection - voir Presentation_15Septembre/Utrecht_Import/Utrecht/BN.col,
et convert_utrecht_to_r.py) : `coreAzimuth`/`coreDip`/`beddingStrike`/
`beddingDip` de BN.col pour "BN4.1A" valent 52/59/182/44, IDENTIQUES aux
colonnes Azimuth/Plunge/Strike/Dip du fichier Montpellier pour ce meme
specimen - confirme qu'il s'agit de la MEME collection, reexportee par un
pipeline logiciel different (Montpellier plutot que Utrecht/PMAG2) :

  - Dec/Inc = direction REPERE SPECIMEN BRUT (PAS geographique) :
    confirme a PLEINE PRECISION FLOTTANTE (143.83182599744745 /
    1.3952186578458248) en recalculant selection.polere(x,y,z) depuis
    les x/y/z BRUTS de BN.col pour le meme palier (BN4.1A, palier 20).
  - Intensity = magnitude du VECTEUR SPECIMEN BRUT, DEJA en A/m
    (magnetisation, PAS moment) : confirme a pleine precision flottante
    contre |x,y,z| de BN.col (uA/m) * 1e-6.
  - Azimuth/Plunge = coreAzimuth/coreDip BRUTS PMAG2 (PAS caz/cin
    STARpaleomag_Py directement) : MEME transformation que
    convert_utrecht_to_r.py, DEJA verifiee sur 156 interpretations
    reelles - `caz = Azimuth + 90`, `cin = 90 - Plunge`. Volontairement
    DIFFERENT de convert_ipgp_to_r.py (qui utilise b= directement) :
    malgre une syntaxe d'en-tete a=/b=/... superficiellement proche dans
    d'autres formats de cette session, Montpellier suit la convention
    PMAG2/Utrecht, pas celle de l'export Johan Guyodo/CryoMag - confirme
    numeriquement ici, pas suppose par analogie.
  - Strike/Dip = bedding BRUTS, ecrits DIRECTEMENT (identiques a
    beddingStrike/beddingDip de BN.col) - meme convention que Utrecht.

Volume : ABSENT de ce format (aucune colonne, contrairement a
Utrecht/IPGP/CryoMag qui portent tous un volume). x/y/z ne peuvent donc
etre reconstruits en MOMENT (Am2, attendu par .prmag) sans supposer un
volume : defaut 10.80 cm3 (MEME valeur par defaut que
ams_prmag.create_prmag_from_legacy_ani - un plug de carotte standard),
personnalisable via `default_volume_cm3`. La direction (dec/inc) et la
forme des ajustements NE DEPENDENT PAS de ce choix (le volume n'affecte
que l'affichage MAG(A/m) absolu) - a corriger manuellement si le volume
reel differe par specimen.

Site : absent egalement (comme Utrecht) - "n.d", a raffiner manuellement.

Type de palier (AF/thermique) : ce format ne porte AUCUN marqueur
explicite (paliers nus, comme la variante "Johan Guyodo" de
convert_ipgp_to_r.py) - MEME heuristique utilisateur reutilisee ici
("si la plupart des etapes sont a moins de 150, c'est de l'AF en mT",
voir convert_ipgp_to_r._is_af). Le NRM est detecte par la VALEUR du
palier (== 0.0 exactement), PAS par sa position dans le bloc : BN4.1A et
BN4.5 commencent tous deux directement a un palier 20 (thermique), SANS
aucun palier 0 dans les donnees - alors que BN4.2A/BN4.3A/BN4.4A
commencent bien a 0.0. Un simple "1er palier du bloc = NRM" (utilise a
tort dans une version anterieure de convert_ipgp_to_r.py) aurait ete
faux pour BN4.1A/BN4.5 - corrige ici des le depart.
"""

import argparse
import math
import os
from collections import OrderedDict
from typing import List, Optional

from magic_export import _measurement_treatment
from convert_ipgp_to_r import _is_af, _nd, _FakeMeasurement

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

DEFAULT_VOLUME_CM3 = 10.80


def read_montpellier_file(path: str) -> List[dict]:
    """Parse un fichier .montpellier -> liste de dict {specimen, caz,
    cin, bed_dip_strike, bed_dip, steps:[{step,intensity,dec,inc}]}, un
    dict par specimen (ordre de premiere apparition preserve)."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    header = lines[0].split()
    idx = {name: i for i, name in enumerate(header)}
    required = ("Sample", "Step", "Intensity", "Dec", "Inc", "Azimuth", "Plunge", "Strike", "Dip")
    missing = [k for k in required if k not in idx]
    if missing:
        raise ValueError(f"missing expected column(s) in .montpellier header: {missing}")

    specimens: "OrderedDict[str, dict]" = OrderedDict()
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) < len(header):
            continue
        sample = fields[idx["Sample"]]
        az = float(fields[idx["Azimuth"]])
        pl = float(fields[idx["Plunge"]])
        sp = specimens.setdefault(sample, {
            "specimen": sample,
            "caz": (az + 90.0) % 360.0,
            "cin": 90.0 - pl,
            "bed_dip_strike": float(fields[idx["Strike"]]),
            "bed_dip": float(fields[idx["Dip"]]),
            "steps": [],
        })
        sp["steps"].append({
            "step": float(fields[idx["Step"]]),
            "intensity": float(fields[idx["Intensity"]]),
            "dec": float(fields[idx["Dec"]]),
            "inc": float(fields[idx["Inc"]]),
        })
    return list(specimens.values())


def _sample_header_block(sp: dict, volume_cm3: float) -> str:
    specimen = sp["specimen"]
    sample = specimen[:-1] if specimen[-1:].isalpha() else specimen
    line_a = (
        f"specimen: {specimen}\tsample: {sample}\tsite: n.d\t"
        f"volume: {volume_cm3:.2f}\tmass: n.d\t"
        f"lat: n.d\tlon: n.d\televation: n.d\t"
        f"stratigraphic_height: n.d\t"
        f"comment: n.d"
    )
    line_b = (
        f"azimuth: {sp['caz']:.1f}\tdip: {sp['cin']:.1f}\t"
        f"date: n.d\tmagnetic_azimuth: n.d\tsolar_azimuth: n.d\torient_tool: n.d"
    )
    line_c = f"bed_dip_strike: {sp['bed_dip_strike']:.1f}\tbed_dip: {sp['bed_dip']:.1f}"
    line_d = (
        "formation: n.d\tage: n.d\t"
        "geologic_classes: n.d\tgeologic_types: n.d\t"
        "lithologies: n.d\tlocation: n.d\t"
        "obs: n.d\tmethod_codes: n.d"
    )
    return "\n".join([line_a, line_b, line_c, line_d])


def _measurement_rows(sp: dict, volume_cm3: float) -> List[str]:
    cod1_demag_fallback = "F" if _is_af([st["step"] for st in sp["steps"]]) else "D"
    to_moment = volume_cm3 * 1.0e-6  # magnetization (A/m) * volume(cm3->m3) = moment (Am2)

    rows = [_MEAS_HEADER]
    for j, st in enumerate(sp["steps"]):
        if st["step"] == 0.0:
            cod1, etape = "N", 0.0
        else:
            cod1, etape = cod1_demag_fallback, st["step"]

        dec_r, inc_r = math.radians(st["dec"]), math.radians(st["inc"])
        x = st["intensity"] * math.cos(inc_r) * math.cos(dec_r) * to_moment
        y = st["intensity"] * math.cos(inc_r) * math.sin(dec_r) * to_moment
        z = st["intensity"] * math.sin(inc_r) * to_moment

        fake_m = _FakeMeasurement(etape=etape, cod1=cod1, cod2="0")
        codes, _temp_k, af_field, _dc, _phi, _theta = _measurement_treatment(fake_m, [], 0.0)

        af_field_mT = af_field * 1.0e3 if cod1 == "F" else None
        temp_c = etape if cod1 in ("N", "D") else 0.0

        row = "\t".join([
            f"{etape:6.1f}", cod1, "0",
            f"{x:11.3E}", f"{y:11.3E}", f"{z:11.3E}",
            "n.d", "g", "n.d", "n.d",
            _nd(temp_c, "6.1f"),
            _nd(af_field_mT, "7.2f"),
            _nd(None, "7.2f"), _nd(0.0, "7.2f"),
            _nd(None, "6.1f"), _nd(None, "6.1f"),
            f"{codes:<40}", "n.d", str(j + 1),
        ])
        rows.append(row)
    return rows


def convert_file(path_in: str, path_out: str, default_volume_cm3: float = DEFAULT_VOLUME_CM3) -> int:
    """Convertit un fichier .montpellier -> .prmag. Retourne le nombre de
    specimens ecrits. Voir docstring module pour le volume par defaut
    (absent de ce format)."""
    specimens = read_montpellier_file(path_in)

    lines = [
        FORMAT_HEADER,
        f"#converted from {os.path.basename(path_in)} by convert_montpellier_to_r.py "
        f"(Montpellier flat format -> .prmag, assumed volume={default_volume_cm3:.2f} cm3)",
        S_FIELD_HEADER,
        "",
    ]
    blocks = [
        _sample_header_block(sp, default_volume_cm3) + "\n"
        + "\n".join(_measurement_rows(sp, default_volume_cm3))
        for sp in specimens
    ]

    with open(path_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        f.write("\n\n".join(blocks) + "\n")

    return len(blocks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a Montpellier .montpellier file to .prmag")
    parser.add_argument("montpellier_file")
    parser.add_argument("-o", "--output")
    parser.add_argument("--volume", type=float, default=DEFAULT_VOLUME_CM3,
                         help="Assumed specimen volume in cm3 (not carried by this format)")
    args = parser.parse_args()
    out = args.output or os.path.splitext(args.montpellier_file)[0] + ".prmag"
    nb = convert_file(args.montpellier_file, out, args.volume)
    print(f"{nb} specimen(s) written to {out} (assumed volume={args.volume:.2f} cm3)")
