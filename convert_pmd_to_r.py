"""
Convertit un DOSSIER de fichiers ".pmd" (un fichier par specimen) vers UN
seul .prmag - demande explicite utilisateur ("ajouter une importation de
donnees a partir de fichiers .PMD et importer tout un dossier"). Construit
et verifie sur le jeu de donnees reel Lu et al. (2026), DatasetS1/
thermal_demag/PMD (334 specimens, 6263 mesures, section de Sigou, Shanxi) :

    <ligne de titre, vide ici>
    JC100-1A  a=264     b= 76     s= 77     d= 43     v=10.0E-6m3
    PAL    Xc    (Am2)    Yc    (Am2)    Zc    (Am2)    MAG(Am2)    Dg    Ig    Ds    Is
    NRM    -8.36E-10    -1.34E-08    7.54E-09    1.54E-08    202.0    9.9    205.6    -25.1
    100C    3.45E-09    -1.36E-08    4.64E-09    1.48E-08    195.4    -8.7    208.4    -44.7

Variante "colonnes separees par des blancs" du format .pmd (PaleoMac/
Enkin) : paliers SUFFIXES "C" ("100C") ou nus ("640"), PAS prefixes
("T020") comme l'export CryoMag deja gere par convert_ipgp_to_r.py - dont
le parseur ne reconnait aucune ligne de ces fichiers (teste), d'ou ce
module separe. Les prefixes T/A/M et suffixes C/mT sont geres aussi.

Conventions VERIFIEES NUMERIQUEMENT (pas supposees) sur les 6263 mesures,
en repassant Xc/Yc/Zc par selection.corfor()/corpen() et en comparant aux
Dg/Ig et Ds/Is du fichier (mediane 0.05 deg, max 0.25 deg = l'arrondi des
3 chiffres significatifs de Xc/Yc/Zc) ; toute autre combinaison testee
(a, a+-90, a+180 x b, 90-b, -b, b-90 ; s, s+-90, s+180 x d, -d) donne
une mediane > 18 deg :

  - azimuth (caz) = a + 90, dip (cin) = b - MEME regle que
    convert_ipgp_to_r.py, PAS celle d'Utrecht/Montpellier (cin = 90 - b).
  - bed_dip_strike = s, bed_dip = d, directs.
  - Xc/Yc/Zc = MOMENT en Am2 (repere specimen), ecrits tels quels ;
    volume = v= (m3 -> cm3). NRM ~1.5e-8 Am2 pour 10 cm3 = 1.5e-3 A/m,
    normal pour des red beds.

Specimen = NOM DU FICHIER, pas le nom de la ligne d'en-tete : celui-ci est
tronque dans ~1 fichier sur 3 ("JC100" dans JC100-2.pmd) ou colle au a=
("JC117.5-1Aa=280"). Site = partie avant le dernier "-" (JC100-1A et
JC100-2 -> site JC100 : un niveau stratigraphique, "site" au sens du
readme des auteurs), "n.d" si le nom n'a pas de "-".

Mesures REPETEES (178 fichiers sur 334, 376 lignes) - readme des auteurs :
"In cases where a single temperature step was measured multiple times,
only the final measurement was used for the principal component
analysis". Trois cas reels, tous traites comme "remplace par une mesure
ulterieure" :
  - meme palier present plusieurs fois (adjacent ou non) : seule la
    DERNIERE occurrence compte ;
  - plusieurs "NRM" en tete de fichier : seule la derniere compte ;
  - un "NRM" au MILIEU de la sequence (7 fichiers) : mesure dont le
    palier n'a pas ete saisi, refaite aussitot sous le bon nom (les 7 cas
    sont quasi identiques a la ligne suivante) - jamais une vraie NRM.
`keep_last=True` (defaut) ne garde que la mesure finale ; False garde
tout, les mesures remplacees etant marquees quality 'b' (voir "Remove
bad-quality (b) steps..."). L'ordre du fichier est conserve dans les deux
cas : 42 fichiers ont des temperatures non croissantes (reprise d'un
palier apres coup) - a regarder specimen par specimen, signale dans le
rapport.

Lignes ignorees : palier non numerique autre que NRM ("Hold", 1 ligne -
mesure du porte-echantillon, 1.8e-10 Am2). Suffixe inconnu ("689t"/"500t",
2 lignes) : lu comme le palier numerique, signale dans le rapport.
"""

import argparse
import os
import re
from typing import Dict, List, Optional, Tuple

from magic_export import _measurement_treatment
from convert_ipgp_to_r import _is_af, _nd, _FakeMeasurement
from convert_montpellier_to_r import FORMAT_HEADER, S_FIELD_HEADER, _MEAS_HEADER

DEFAULT_VOLUME_CM3 = 10.0

_NUM = r"-?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?"
_ORIENT_RE = re.compile(
    r"a=\s*(" + _NUM + r")\s+b=\s*(" + _NUM + r")\s+s=\s*(" + _NUM + r")\s+d=\s*(" + _NUM + r")")
_VOLUME_RE = re.compile(r"v=\s*(" + _NUM + r")\s*m3", re.IGNORECASE)
_STEP_RE = re.compile(r"^([A-Za-z]*)(\d+(?:\.\d*)?)([A-Za-z°]*)$")


def list_pmd_files(folder: str) -> List[str]:
    """Tous les .pmd du dossier ET de ses sous-dossiers (le dossier fourni
    par l'utilisateur, thermal_demag/, contient les fichiers dans un
    sous-dossier PMD/), en ordre naturel (JC2 avant JC10, JC117.5 apres
    JC117)."""
    found = []
    for root, _dirs, files in os.walk(folder):
        for name in files:
            if name.lower().endswith(".pmd") and not name.startswith("."):
                found.append(os.path.join(root, name))

    def natural_key(path):
        stem = os.path.splitext(os.path.basename(path))[0]
        return [(0, float(tok), "") if tok[:1].isdigit() else (1, 0.0, tok.lower())
                for tok in re.findall(r"\d+(?:\.\d+)?|\D+", stem)]

    return sorted(found, key=lambda p: (natural_key(p), p))


def _decode_step_label(label: str) -> Optional[Tuple[Optional[str], float, bool]]:
    """-> (cod1 | None si le type n'est pas dans le label, valeur,
    suffixe/prefixe non reconnu ?), ou None si le label n'est pas un palier
    (ex. "Hold")."""
    label = label.strip()
    if label.upper() == "NRM":
        return "N", 0.0, False
    m = _STEP_RE.match(label)
    if not m:
        return None
    prefix, value, suffix = m.group(1).upper(), float(m.group(2)), m.group(3).upper()
    if prefix == "T" or suffix in ("C", "°C") and not prefix:
        return "D", value, False
    if prefix in ("A", "M", "AF") and not suffix or suffix == "MT" and not prefix:
        return "F", value, False
    return None, value, bool(prefix or suffix)


def read_pmd_file(path: str) -> Optional[dict]:
    """Parse un fichier .pmd -> dict {specimen, site, caz, cin,
    bed_dip_strike, bed_dip, volume_cm3, steps:[{label, cod1, etape, x, y,
    z, superseded}], skipped:[labels], odd:[labels]}, ou None si le
    fichier n'a pas de ligne d'en-tete a=/b=/s=/d=."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln.rstrip("\r\n") for ln in f]

    header_idx = next((i for i, ln in enumerate(lines) if _ORIENT_RE.search(ln)), None)
    if header_idx is None:
        return None
    header = lines[header_idx]
    a, b, s, d = (float(v) for v in _ORIENT_RE.search(header).groups())
    vol = _VOLUME_RE.search(header)
    volume_cm3 = float(vol.group(1)) * 1.0e6 if vol else DEFAULT_VOLUME_CM3

    raw, skipped, odd = [], [], []
    for ln in lines[header_idx + 1:]:
        fields = ln.split()
        if len(fields) < 4:
            continue
        try:
            x, y, z = float(fields[1]), float(fields[2]), float(fields[3])
        except ValueError:
            continue  # ligne des noms de colonnes ("PAL  Xc  (Am2) ...")
        decoded = _decode_step_label(fields[0])
        if decoded is None:
            skipped.append(fields[0])
            continue
        cod1, value, is_odd = decoded
        if is_odd:
            odd.append(fields[0])
        raw.append({"label": fields[0], "cod1": cod1, "value": value, "x": x, "y": y, "z": z})

    # Type des paliers sans marqueur (nus, ou suffixe inconnu) : celui des
    # paliers explicites du meme fichier s'il y en a, sinon l'heuristique
    # utilisateur de convert_ipgp_to_r._is_af.
    explicit = {r["cod1"] for r in raw if r["cod1"] in ("D", "F")}
    if len(explicit) == 1:
        fallback = explicit.pop()
    else:
        fallback = "F" if _is_af([r["value"] for r in raw if r["cod1"] is None]) else "D"
    steps = []
    for r in raw:
        cod1 = r["cod1"]
        if cod1 is None:
            cod1 = "N" if r["value"] == 0.0 else fallback
        steps.append({
            "label": r["label"], "cod1": cod1, "etape": 0.0 if cod1 == "N" else r["value"],
            "x": r["x"], "y": r["y"], "z": r["z"], "superseded": False,
        })

    # Mesures remplacees par une mesure ulterieure (voir docstring module).
    n_lead = 0
    while n_lead < len(steps) and steps[n_lead]["cod1"] == "N":
        n_lead += 1
    last_seen: Dict[Tuple[str, float], int] = {}
    for i, st in enumerate(steps):
        if st["cod1"] == "N":
            st["superseded"] = i != n_lead - 1
            continue
        key = (st["cod1"], st["etape"])
        if key in last_seen:
            steps[last_seen[key]]["superseded"] = True
        last_seen[key] = i

    specimen = os.path.splitext(os.path.basename(path))[0]
    return {
        "specimen": specimen,
        "site": specimen.rsplit("-", 1)[0] if "-" in specimen else "n.d",
        "caz": (a + 90.0) % 360.0, "cin": b,
        "bed_dip_strike": s, "bed_dip": d,
        "volume_cm3": volume_cm3,
        "steps": steps, "skipped": skipped, "odd": odd,
    }


def _sample_header_block(sp: dict) -> str:
    specimen = sp["specimen"]
    sample = specimen[:-1] if specimen[-1:].isalpha() else specimen
    line_a = (
        f"specimen: {specimen}\tsample: {sample}\tsite: {sp['site']}\t"
        f"volume: {sp['volume_cm3']:.2f}\tmass: n.d\t"
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


def _measurement_rows(sp: dict, keep_last: bool) -> List[str]:
    rows = [_MEAS_HEADER]
    steps = [st for st in sp["steps"] if not (keep_last and st["superseded"])]
    for j, st in enumerate(steps):
        cod1, etape = st["cod1"], st["etape"]
        fake_m = _FakeMeasurement(etape=etape, cod1=cod1, cod2="0")
        codes, _temp_k, af_field, _dc, _phi, _theta = _measurement_treatment(fake_m, [], 0.0)
        af_field_mT = af_field * 1.0e3 if cod1 == "F" else None
        temp_c = etape if cod1 in ("N", "D") else 0.0
        rows.append("\t".join([
            f"{etape:6.1f}", cod1, "0",
            f"{st['x']:11.3E}", f"{st['y']:11.3E}", f"{st['z']:11.3E}",
            "n.d", "b" if st["superseded"] else "g", "n.d", "n.d",
            _nd(temp_c, "6.1f"),
            _nd(af_field_mT, "7.2f"),
            _nd(None, "7.2f"), _nd(0.0, "7.2f"),
            _nd(None, "6.1f"), _nd(None, "6.1f"),
            f"{codes:<40}", "n.d", str(j + 1),
        ]))
    return rows


def _short_list(names: List[str], limit: int = 12) -> str:
    return ", ".join(names[:limit]) + (f", ... (+{len(names) - limit})" if len(names) > limit else "")


def convert_files(paths_in: List[str], path_out: str, keep_last: bool = True) -> Tuple[int, List[str]]:
    """Convertit une liste de fichiers .pmd -> UN .prmag. Retourne (nombre
    de specimens ecrits, lignes de rapport)."""
    specimens, unreadable, duplicates = [], [], []
    seen = set()
    for path in paths_in:
        sp = read_pmd_file(path)
        if sp is None or not sp["steps"]:
            unreadable.append(os.path.basename(path))
            continue
        if sp["specimen"] in seen:
            duplicates.append(sp["specimen"])
            continue
        seen.add(sp["specimen"])
        specimens.append(sp)
    if not specimens:
        raise ValueError("no readable .pmd file (no 'a= b= s= d=' header line found)")

    source = os.path.basename(os.path.commonpath([os.path.abspath(p) for p in paths_in]))
    lines = [
        FORMAT_HEADER,
        f"#converted from {len(specimens)} .pmd file(s) in {source} by convert_pmd_to_r.py "
        f"(PMD -> .prmag, repeated steps: "
        f"{'last measurement only' if keep_last else 'all kept, superseded ones quality b'})",
        S_FIELD_HEADER,
        "",
    ]
    blocks = [
        _sample_header_block(sp) + "\n" + "\n".join(_measurement_rows(sp, keep_last))
        for sp in specimens
    ]
    with open(path_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        f.write("\n\n".join(blocks) + "\n")

    n_meas = sum(len(sp["steps"]) for sp in specimens)
    repeated = [sp["specimen"] for sp in specimens if any(st["superseded"] for st in sp["steps"])]
    n_superseded = sum(1 for sp in specimens for st in sp["steps"] if st["superseded"])
    mid_nrm = [sp["specimen"] for sp in specimens
               if any(st["cod1"] == "N" and i > 0 and sp["steps"][i - 1]["cod1"] != "N"
                      for i, st in enumerate(sp["steps"]))]
    unordered = []
    for sp in specimens:
        kept = [st["etape"] for st in sp["steps"] if not st["superseded"]]
        if any(t2 < t1 for t1, t2 in zip(kept, kept[1:])):
            unordered.append(sp["specimen"])
    sites = {sp["site"] for sp in specimens}
    types = {st["cod1"] for sp in specimens for st in sp["steps"]} - {"N"}

    report = [
        f"{len(specimens)} specimen(s), {len(sites)} site(s), "
        f"{n_meas - (n_superseded if keep_last else 0)} measurement(s) written "
        f"({'thermal' if types == {'D'} else 'AF' if types == {'F'} else 'thermal and AF'}).",
    ]
    if repeated:
        report.append(
            f"Repeated steps in {len(repeated)} specimen(s): {n_superseded} earlier measurement(s) "
            + ("dropped, only the last one of each step kept." if keep_last
               else "kept and marked quality 'b' (the last one of each step is 'g')."))
    if mid_nrm:
        report.append(
            f"'NRM' label in the middle of the sequence (unlabelled measurement, superseded "
            f"by the next step) in {len(mid_nrm)}: {_short_list(mid_nrm)}")
    if unordered:
        report.append(
            f"Steps not in increasing order (a step redone later; file order kept) in "
            f"{len(unordered)}: {_short_list(unordered)}")
    odd = [f"{sp['specimen']} ({', '.join(sp['odd'])})" for sp in specimens if sp["odd"]]
    if odd:
        report.append(f"Unknown step suffix, read as the plain step value: {_short_list(odd)}")
    skipped = [f"{sp['specimen']} ({', '.join(sp['skipped'])})" for sp in specimens if sp["skipped"]]
    if skipped:
        report.append(f"Non-step line(s) ignored: {_short_list(skipped)}")
    if duplicates:
        report.append(f"Same specimen name in two files, second one ignored: {_short_list(duplicates)}")
    if unreadable:
        report.append(f"Not a readable .pmd, ignored: {_short_list(unreadable)}")
    return len(specimens), report


def convert_folder(folder: str, path_out: str, keep_last: bool = True) -> Tuple[int, List[str]]:
    paths = list_pmd_files(folder)
    if not paths:
        raise ValueError(f"no .pmd file found in {folder} (sub-folders included)")
    return convert_files(paths, path_out, keep_last)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a folder of .pmd files to one .prmag")
    parser.add_argument("folder")
    parser.add_argument("-o", "--output")
    parser.add_argument("--keep-all", action="store_true",
                         help="Keep repeated measurements (superseded ones marked quality b)")
    args = parser.parse_args()
    out = args.output or os.path.join(args.folder, os.path.basename(os.path.normpath(args.folder)) + ".prmag")
    nb, rep = convert_folder(args.folder, out, keep_last=not args.keep_all)
    print(f"{nb} specimen(s) written to {out}")
    print("\n".join(rep))
