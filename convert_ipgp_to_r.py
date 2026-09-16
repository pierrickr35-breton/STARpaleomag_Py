"""
Convertit les fichiers du format "IPGP" (Institut de Physique du Globe de
Paris - demande explicite utilisateur, fichier d'exemple fourni
explicitement signale incomplet : "voici un nouveau format (IPGP) que je
souhaite convertir. Ce n'est pas un exemple complet.") vers .prmag/.pmagres.
Deux fichiers INDEPENDANTS a convertir (confirme par l'utilisateur -
"1" = les deux) :

1) Fichier de mesures brutes ("Remanence measurement file", attribue a
   Johan Guyodo) -> .prmag. Un bloc par specimen :

       <site>
       <specimen>  a=<az>  b=<dip>  s=<strike>  d=<dip>  v=<vol>m3  <date> <heure> [PAL <en-tete colonnes>]
       [PAL <en-tete colonnes>]                      (si pas deja sur la ligne precedente)
       <palier> <Xc> <Yc> <Zc> <MAG> <Dg> <Ig> <Ds> <Is> <a95>
       ...

   Colonnes (Xc/Yc/Zc en Am2 DEJA, pas de conversion d'unite necessaire -
   contrairement a Utrecht/uA/m, voir convert_utrecht_to_r.py) : palier,
   Xc, Yc, Zc, MAG(A/m), Dg/Ig (declinaison/inclinaison GEOGRAPHIQUE,
   apres correction carotte), Ds/Is (declinaison/inclinaison
   STRUCTURALE, apres correction pendage en plus), a95 (precision
   instrument par mesure, non reutilisee - .prmag n'a pas de colonne
   pour cette valeur au niveau mesure).

   Ligne "site" : un token seul juste avant l'en-tete specimen (ex.
   "Hz12"/"hz12" dans l'exemple fourni) - hypothese retenue faute
   d'autre source de site dans ce format (contrairement a l'ANI legacy
   Rennes ou au nom de fichier Utrecht) : CHAQUE bloc specimen est
   precede d'un token seul qui reapparait identique (aux capitales pres)
   sur les 3 specimens de l'exemple fourni, cense denoter un
   regroupement (site/etude) - A CONFIRMER par l'utilisateur des qu'un
   exemple plus complet sera disponible (fichier explicitement signale
   incomplet).

   a=/b=/s=/d= -> caz/cin/bed_dip_strike/bed_dip : transformation
   VERIFIEE NUMERIQUEMENT (pas une hypothese) sur les 4 mesures "step 0"
   des 3 specimens de l'exemple fourni, en comparant Dg/Ig/Ds/Is donnes
   par le fichier au resultat de selection.corfor/corpen appliques a
   Xc/Yc/Zc :

       caz = (a + 90.0) % 360.0     (meme decalage +90 que Utrecht/MagIC,
                                      voir convert_utrecht_to_r.py)
       cin = b                       (directement, PAS de 90-b)
       bed_dip_strike = s            (directement, PAS de +90 ici)
       bed_dip = d                   (directement)

   Concordance Dg/Ig et Ds/Is a 0.1 pres sur les 4 points testes
   (0975A/0765B/0490D, plusieurs paliers) - transformation consideree
   fiable pour ce format.

   AF vs thermique : ce format ne porte AUCUN marqueur explicite du type
   de desaimantation. Heuristique donnee par l'utilisateur ("si la
   plupart des etapes sont a moins de 150, c'est de l'AF en mT") :
   applique par specimen, sur les paliers non nuls de son propre bloc -
   voir `_is_af`.

2) Fichier "Results" IPGP (CSV, developpe par Dragomir Dragomirov,
   dragomirov@ipgp.fr) -> .pmagres. Ajustements PCA DEJA calcules par
   leur logiciel, donnes dans 3 reperes (sample/geo/tect) x 2 versions
   (libre/ancree a l'origine) par ligne. Repere retenu pour
   FitResult.dec/inc : decl_sample/incl_sample (repere specimen BRUT,
   MEME convention que le reste de .pmagres - correction affichee a la
   demande via calcul._correct_dec_inc, pas stockee corrigee - voir
   convert_utrecht_to_r.py). Chaque ligne produit DEUX FitResult (le
   fichier IPGP donne toujours les deux solutions) :
     - orig='n' (libre)  : decl_sample/incl_sample, mad
     - orig='o' (ancree) : decl_sample_anchored/incl_sample_anchored, amad
   alph95/aalph95/mdf/demag_loss_* n'ont pas de colonne equivalente dans
   FitResult (deja le cas pour MAD seul avec Utrecht) - non reportes.
   cat1 force a 'L' (ligne) : aucune colonne de ce format ne distingue
   ligne/plan (pas de marqueur equivalent au "TAU1"/"TAU3" d'Utrecht).
   numcomp = component_id (reutilise tel quel).

   site : ce format n'a PAS de colonne site, et .pmagres ne stocke de
   toute facon jamais le site par resultat (seul le .prmag compagnon le
   porte, par specimen) - rien de particulier a faire ici, le site se
   resout normalement des que le .prmag et le .pmagres sont charges
   ensemble dans l'appli.

   Doublons EXACTS detectes (l'exemple fourni contient deux lignes
   "111A" strictement identiques, lignes 71/72 du fichier - signale a
   l'utilisateur, probable artefact de l'exemple plutot qu'une
   caracteristique reelle du format) : ecartes silencieusement au sens
   ou ils ne sont PAS ecrits deux fois (contrairement aux variantes
   jackknife d'un .ANI, qui sont un phenomene reel documente - voir
   ams_prmag.create_prmag_from_legacy_ani), mais comptes et rapportes
   par `convert_results_file`.
"""

import argparse
import ast
import csv
import os
import re
from typing import Dict, List, Optional, Tuple

from calcul import FitResult, archivres, results_path_for
from magic_export import _measurement_treatment

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

_NUM = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
_HEADER_LINE_RE = re.compile(
    r"^(?P<specimen>\S+)\s+a=\s*(?P<a>" + _NUM + r")\s+b=\s*(?P<b>" + _NUM + r")\s+"
    r"s=\s*(?P<s>" + _NUM + r")\s+d=\s*(?P<d>" + _NUM + r")\s+v=\s*(?P<v>" + _NUM + r")\s*m3\s*"
    r"(?P<date>\d{2}-\d{2}-\d{4})\s+(?P<time>\d{2}:\d{2})"
)


def _nd(value, fmt: Optional[str] = None) -> str:
    if value is None:
        if fmt:
            width_digits = "".join(ch for ch in fmt.split(".")[0] if ch.isdigit())
            if width_digits:
                return "n.d".rjust(int(width_digits))
        return "n.d"
    return format(value, fmt) if fmt is not None else str(value)


def _is_af(step_values: List[float]) -> bool:
    """Heuristique donnee par l'utilisateur : "si la plupart des etapes
    sont a moins de 150, c'est de l'AF en mT" (sinon thermique en degC).
    Applique aux paliers non nuls (le NRM, palier 0, ne renseigne pas sur
    le type de desaimantation)."""
    nonzero = [v for v in step_values if v != 0.0]
    if not nonzero:
        return True
    below = sum(1 for v in nonzero if v < 150.0)
    return below >= len(nonzero) / 2.0


def _iso_date(date_str: Optional[str], time_str: Optional[str]) -> str:
    """"08-04-2008" -> suppose DD-MM-YYYY (convention europeenne, source
    IPGP/Paris - PAS verifie sur un exemple complet, a confirmer)."""
    if not date_str:
        return "n.d"
    try:
        dd, mm, yyyy = date_str.split("-")
        iso = f"{yyyy}-{mm}-{dd}"
    except ValueError:
        return "n.d"
    return f"{iso}T{time_str}" if time_str else iso


class _FakeMeasurement:
    """Juste assez de champs pour reutiliser magic_export._measurement_treatment
    (voir convert_utrecht_to_r.py, meme usage)."""
    def __init__(self, etape: float, cod1: str, cod2: str = "0"):
        self.etape = etape
        self.cod1 = cod1
        self.cod2 = cod2
        self.ins = None


def read_ipgp_measurements(path: str) -> List[dict]:
    """Parse un fichier "Remanence measurement file" IPGP -> liste de
    dict {specimen, site, caz, cin, bed_dip_strike, bed_dip, volume_cm3,
    date, time, steps:[{step,x,y,z,mag,dg,ig,ds,is_,a95}]}.

    Classification ligne par ligne (voir docstring module pour le detail
    des 3 hypotheses verifiees/a confirmer) :
      - une ligne contenant "a=" ET "v=" demarre un nouveau bloc specimen ;
      - une ligne commencant par "PAL" (en-tete de colonnes, sur la meme
        ligne que le specimen ou sur la suivante selon les blocs de
        l'exemple fourni) est ignoree (colonnes fixes, connues d'avance) ;
      - une ligne dont le premier token est numerique, a l'interieur d'un
        bloc specimen actif, est une ligne de mesure ;
      - toute autre ligne non vide est retenue comme "site en attente"
        (ecrase par la ligne suivante du meme type tant qu'aucun bloc
        specimen ne l'a consommee - absorbe donc sans effet de bord un
        titre de fichier avant la premiere vraie ligne de site)."""
    specimens: List[dict] = []
    current: Optional[dict] = None
    pending_site: Optional[str] = None

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        raw_lines = f.readlines()

    for raw in raw_lines:
        line = raw.strip()
        if not line:
            continue

        m = _HEADER_LINE_RE.match(line)
        if m:
            if current is not None:
                specimens.append(current)
            a, b, s, d = (float(m.group(k)) for k in ("a", "b", "s", "d"))
            current = {
                "specimen": m.group("specimen"),
                "site": pending_site or "n.d",
                "caz": (a + 90.0) % 360.0,
                "cin": b,
                "bed_dip_strike": s,
                "bed_dip": d,
                "volume_cm3": float(m.group("v")) * 1.0e6,  # m3 -> cm3
                "date": m.group("date"),
                "time": m.group("time"),
                "steps": [],
            }
            continue

        if line.upper().startswith("PAL"):
            # en-tete de colonnes (deja connu, colonnes fixes) - ignore
            continue

        first_token = line.split()[0]
        try:
            step_val = float(first_token)
            is_data_row = True
        except ValueError:
            is_data_row = False

        if is_data_row and current is not None:
            fields = line.split()
            if len(fields) < 9:
                continue
            step, xc, yc, zc, mag, dg, ig, ds, is_val = (float(v) for v in fields[:9])
            a95 = float(fields[9]) if len(fields) > 9 else None
            current["steps"].append({
                "step": step, "x": xc, "y": yc, "z": zc, "mag": mag,
                "dg": dg, "ig": ig, "ds": ds, "is_": is_val, "a95": a95,
            })
            continue

        # ni en-tete specimen, ni ligne PAL, ni ligne de mesure -> site
        pending_site = line

    if current is not None:
        specimens.append(current)
    return specimens


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
        f"date: {_iso_date(sp.get('date'), sp.get('time'))}\t"
        f"magnetic_azimuth: n.d\tsolar_azimuth: n.d\torient_tool: n.d"
    )
    line_c = f"bed_dip_strike: {sp['bed_dip_strike']:.1f}\tbed_dip: {sp['bed_dip']:.1f}"
    line_d = (
        "formation: n.d\tage: n.d\t"
        "geologic_classes: n.d\tgeologic_types: n.d\t"
        "lithologies: n.d\tlocation: n.d\t"
        "obs: n.d\tmethod_codes: n.d"
    )
    return "\n".join([line_a, line_b, line_c, line_d])


def _measurement_rows(sp: dict) -> List[str]:
    is_af = _is_af([st["step"] for st in sp["steps"]])
    cod1_demag = "F" if is_af else "D"
    rows = [_MEAS_HEADER]
    for j, st in enumerate(sp["steps"]):
        cod1 = "N" if j == 0 else cod1_demag
        step_val = st["step"]
        etape = 0.0 if cod1 == "N" else step_val
        fake_m = _FakeMeasurement(etape=etape, cod1=cod1, cod2="0")
        codes, _temp_k, af_field, _dc, _phi, _theta = _measurement_treatment(fake_m, [], 0.0)

        af_field_mT = af_field * 1.0e3 if cod1 == "F" else None
        temp_c = etape if cod1 in ("N", "D") else 0.0
        step_str = 0.0 if cod1 == "N" else step_val

        row = "\t".join([
            f"{step_str:6.1f}", cod1, "0",
            f"{st['x']:11.3E}", f"{st['y']:11.3E}", f"{st['z']:11.3E}",
            "n.d", "g", "n.d", "n.d",
            _nd(temp_c, "6.1f"),
            _nd(af_field_mT, "7.2f"),
            _nd(None, "7.2f"), _nd(0.0, "7.2f"),
            _nd(None, "6.1f"), _nd(None, "6.1f"),
            f"{codes:<40}", "n.d", str(j + 1),
        ])
        rows.append(row)
    return rows


def convert_measurements_file(path_in: str, path_out: str) -> int:
    """Convertit un fichier de mesures brutes IPGP -> .prmag. Retourne le
    nombre de specimens ecrits."""
    specimens = read_ipgp_measurements(path_in)

    lines = [
        FORMAT_HEADER,
        f"#converted from {os.path.basename(path_in)} by convert_ipgp_to_r.py "
        "(IPGP raw measurement format -> .prmag)",
        S_FIELD_HEADER,
        "",
    ]
    blocks = [
        _sample_header_block(sp) + "\n" + "\n".join(_measurement_rows(sp))
        for sp in specimens
    ]

    with open(path_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        f.write("\n\n".join(blocks) + "\n")

    return len(blocks)


def _parse_steps_list(raw: str) -> List[float]:
    try:
        return [float(v) for v in ast.literal_eval(raw)]
    except (ValueError, SyntaxError):
        return []


def _fit_results_for_row(row: dict) -> List[FitResult]:
    specimen = row["sample"]
    steps = _parse_steps_list(row.get("steps", "[]"))
    demag = "F" if _is_af(steps) else "D"
    step_first = int(round(min(steps))) if steps else 0
    step_last = int(round(max(steps))) if steps else 0
    numcomp = int(float(row.get("component_id") or 1))
    nb = int(float(row.get("number_points") or len(steps)))

    def _f(key: str) -> float:
        try:
            return float(row[key])
        except (KeyError, ValueError):
            return 0.0

    common = dict(
        id=specimen, cat1="L", cat2="", demag=demag, numcomp=numcomp, nb=nb,
        step_first=step_first, step_last=step_last,
    )
    results = [
        FitResult(orig="n", dec=_f("decl_sample"), inc=_f("incl_sample"),
                   mad=_f("mad"), **common),
        FitResult(orig="o", dec=_f("decl_sample_anchored"), inc=_f("incl_sample_anchored"),
                   mad=_f("amad"), **common),
    ]
    return results


def convert_results_file(path_in: str, path_out: str) -> Tuple[int, int]:
    """Convertit le fichier "Results" IPGP (CSV) -> .pmagres. Retourne
    (nb resultats ecrits, nb lignes doublons exactes ecartees). Pas de
    notion de site (.pmagres ne stocke jamais le site par resultat, voir
    docstring module - un specimen commun au fichier de mesures resout
    deja son site via le .prmag compagnon une fois les deux charges
    ensemble dans l'appli)."""
    if os.path.exists(path_out):
        # meme raison que convert_utrecht_to_r.convert_files : eviter de
        # dupliquer les interpretations a une re-conversion vers le meme
        # chemin de sortie (archivres ajoute, il n'ecrase pas).
        os.remove(path_out)

    seen_rows = set()
    nb_written = 0
    nb_duplicates = 0
    existing_ids: Optional[set] = None

    with open(path_in, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = tuple(row.get(k, "") for k in reader.fieldnames)
            if key in seen_rows:
                nb_duplicates += 1
                continue
            seen_rows.add(key)
            for res in _fit_results_for_row(row):
                _c, existing_ids = archivres(res, path_out, existing_ids)
                nb_written += 1

    return nb_written, nb_duplicates


def convert_files(
    measurements_path: Optional[str] = None,
    results_path: Optional[str] = None,
    prmag_out: Optional[str] = None,
) -> Tuple[int, int, int]:
    """Convertit un fichier de mesures et/ou un fichier de resultats IPGP.
    `prmag_out` : chemin du .prmag de sortie (deduit du fichier de mesures
    si non fourni) - le .pmagres compagnon est deduit via results_path_for.
    Retourne (nb specimens, nb resultats, nb doublons ecartes)."""
    if not measurements_path and not results_path:
        raise ValueError("at least one of measurements_path/results_path is required")

    if not prmag_out:
        base = measurements_path or results_path
        prmag_out = os.path.splitext(base)[0] + ".prmag"

    nb_specimens = 0
    if measurements_path:
        nb_specimens = convert_measurements_file(measurements_path, prmag_out)

    nb_results = nb_duplicates = 0
    if results_path:
        nb_results, nb_duplicates = convert_results_file(results_path, results_path_for(prmag_out))

    return nb_specimens, nb_results, nb_duplicates


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert IPGP measurement/results files to .prmag/.pmagres")
    parser.add_argument("-m", "--measurements", help="Raw measurement file (Remanence measurement file)")
    parser.add_argument("-r", "--results", help="Results CSV file (component_id,sample,... header)")
    parser.add_argument("-o", "--output", help="Output .prmag path")
    args = parser.parse_args()
    nb_sp, nb_res, nb_dup = convert_files(args.measurements, args.results, args.output)
    out = args.output or os.path.splitext(args.measurements or args.results)[0] + ".prmag"
    print(f"{nb_sp} specimen(s) written to {out}")
    print(f"{nb_res} result(s) written to {results_path_for(out)} ({nb_dup} exact duplicate row(s) skipped)")
