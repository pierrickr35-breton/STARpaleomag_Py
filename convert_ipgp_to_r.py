"""
Convertit les fichiers du format "IPGP" (Institut de Physique du Globe de
Paris - demande explicite utilisateur, fichier d'exemple fourni
explicitement signale incomplet : "voici un nouveau format (IPGP) que je
souhaite convertir. Ce n'est pas un exemple complet.") vers .prmag/.pmagres.
Deux fichiers INDEPENDANTS a convertir (confirme par l'utilisateur -
"1" = les deux) :

1) Fichier de mesures brutes -> .prmag. Un bloc par specimen. DEUX
   variantes reelles vues jusqu'ici, gerees par le MEME parseur
   (`read_ipgp_measurements`) :

   a) "Remanence measurement file" (attribue a Johan Guyodo) - paliers
      NUS (nombre seul, ex. "0", "5", ... "200"), et un token seul
      (ex. "Hz12"/"hz12") juste avant chaque en-tete specimen :

          Hz12
          0975A     a=340.0   b=-12.0   s=340.0   d=  1.0   v= 8.0E-6m3  08-04-2008 20:50 PAL  Xc (Am2) ...
          PAL  Xc (Am2)  Yc (Am2)  Zc (Am2)  MAG(A/m)   Dg    Ig    Ds    Is   a95
          0   4.10E-08 -2.52E-09  3.27E-08  6.56E-03 335.7  50.5 336.9  50.6  1.5

   b) Export CryoMag/.pmd (PaleoMac - "CryoMag 2.0c-User modified data
      file exported to PaleoMac", demande explicite utilisateur "pour
      les fichiers .pmd") - pas de ligne "site" (la 1ere ligne du
      fichier est un titre logiciel fixe, PAS un site), et paliers
      PREFIXES d'une lettre indiquant directement le type (ex. "T020" =
      thermique 20 degC) :

          K15371    a= 25.0   b= 39.0   s=314.0   d= 42.0   v=11.0E-6m3  01-03-2016 00:05
          STEP  Xc (Am2)  Yc (Am2)  Zc (Am2)  MAG(A/m)   Dg    Ig    Ds    Is   a95
          T020  3.07E-08 -1.25E-08  5.85E-08  6.11E-03  13.4  22.9  15.3 -13.3  0.1

   Colonnes (Xc/Yc/Zc en Am2 DEJA, pas de conversion d'unite necessaire -
   contrairement a Utrecht/uA/m, voir convert_utrecht_to_r.py) : palier,
   Xc, Yc, Zc, MAG(A/m), Dg/Ig (declinaison/inclinaison GEOGRAPHIQUE,
   apres correction carotte), Ds/Is (declinaison/inclinaison
   STRUCTURALE, apres correction pendage en plus), a95 (precision
   instrument par mesure, non reutilisee - .prmag n'a pas de colonne
   pour cette valeur au niveau mesure).

   Ligne "site" : reservee aux lignes a UN SEUL token (ex. "Hz12") -
   distingue une vraie ligne de site d'un titre logiciel multi-mots
   (variante b, qui n'en a pas) ; reste une hypothese A CONFIRMER pour
   la variante (a) (aucune autre source de site dans ce format,
   contrairement a l'ANI legacy Rennes ou au nom de fichier Utrecht).

   Type de palier (cod1/etape), voir `_decode_step_label` : un label
   "NRM" -> N/0.0 ; un label prefixe d'une lettre + un nombre -> le
   PREFIXE donne le type directement ('T...' -> thermique/D, 'A...' ->
   AF/F, sans ambiguite, pas d'heuristique necessaire pour la variante
   CryoMag/.pmd) ; un label NUMERIQUE NU (variante Johan Guyodo, sans
   prefixe) -> position 0 du bloc = N, sinon heuristique donnee par
   l'utilisateur ("si la plupart des etapes sont a moins de 150, c'est
   de l'AF en mT") appliquee aux seuls paliers nus du bloc - voir
   `_is_af`.

   a=/b=/s=/d= -> caz/cin/bed_dip_strike/bed_dip : transformation
   VERIFIEE NUMERIQUEMENT (pas une hypothese), sur les DEUX variantes
   (4 mesures de l'exemple Johan Guyodo ET 5 mesures de l'exemple
   CryoMag/.pmd K15371, sources differentes), en comparant Dg/Ig/Ds/Is
   donnes par chaque fichier au resultat de selection.corfor/corpen
   appliques a Xc/Yc/Zc :

       caz = (a + 90.0) % 360.0     (meme decalage +90 que Utrecht/MagIC,
                                      voir convert_utrecht_to_r.py)
       cin = b                       (directement, PAS de 90-b)
       bed_dip_strike = s            (directement, PAS de +90 ici)
       bed_dip = d                   (directement)

   Concordance Dg/Ig a 0.1 pres sur les 9 points testes des deux
   fichiers. Ds/Is concorde exactement sur l'exemple Johan Guyodo, avec
   un residu de 0.3 a 1.3 degre sur l'exemple CryoMag/.pmd (K15371) -
   coherent avec l'arrondi a 3 chiffres significatifs de Xc/Yc/Zc
   (residu qui grandit avec les mesures aux composantes les plus
   faibles, et se propage/amplifie a la 2e rotation corpen alors que
   Dg/Ig - une seule rotation corfor - reste quasi exact) : PAS un
   signe d'erreur de convention, mais a garder en tete si un ecart plus
   grand apparait sur des donnees reelles completes.

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
    """Parse un fichier de mesures brutes IPGP/CryoMag -> liste de dict
    {specimen, site, caz, cin, bed_dip_strike, bed_dip, volume_cm3, date,
    time, steps:[{step_label,x,y,z,mag,dg,ig,ds,is_,a95}]} - gere les
    deux variantes documentees dans la docstring du module (paliers nus
    vs prefixes, avec/sans ligne de site).

    Classification ligne par ligne :
      - une ligne contenant "a=" ET "v=" demarre un nouveau bloc specimen ;
      - une ligne contenant "Xc" et "(Am2)" est un en-tete de colonnes
        (que le 1er mot soit "PAL" ou "STEP") - ignoree (colonnes fixes,
        connues d'avance) ;
      - une ligne d'au moins 9 champs dont les champs 2 a 9 sont
        numeriques, a l'interieur d'un bloc specimen actif, est une
        ligne de mesure (le 1er champ, le "palier", est garde tel quel -
        `step_label` - decode plus tard par `_decode_step_label`) ;
      - une ligne restante A UN SEUL token (pas de sequence de mesure
        valide, pas d'en-tete) est retenue comme "site en attente"
        (ecrasee par la ligne suivante du meme type tant qu'aucun bloc
        specimen ne l'a consommee) - une ligne a PLUSIEURS mots (ex. un
        titre logiciel comme "CryoMag 2.0c-User modified data file
        exported to PaleoMac") n'est PAS retenue comme site."""
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

        if "Xc" in line and "(Am2)" in line:
            # en-tete de colonnes ("PAL ..." ou "STEP ...") - ignore
            continue

        fields = line.split()
        is_data_row = False
        if current is not None and len(fields) >= 9:
            try:
                xc, yc, zc, mag, dg, ig, ds, is_val = (float(v) for v in fields[1:9])
                is_data_row = True
            except ValueError:
                is_data_row = False

        if is_data_row:
            a95 = float(fields[9]) if len(fields) > 9 else None
            current["steps"].append({
                "step_label": fields[0], "x": xc, "y": yc, "z": zc, "mag": mag,
                "dg": dg, "ig": ig, "ds": ds, "is_": is_val, "a95": a95,
            })
            continue

        if len(fields) == 1:
            pending_site = line

    if current is not None:
        specimens.append(current)
    return specimens


def _decode_step_label(label: str) -> Optional[Tuple[str, float]]:
    """"NRM" -> ('N', 0.0) ; un label prefixe d'une ou plusieurs lettres
    suivies d'un nombre -> type donne par la 1ere lettre ('T...' =
    thermique/'D', 'A...' = AF/'F') ; sinon (label numerique nu ou non
    reconnu) -> None, l'appelant retombe sur l'heuristique de position/
    majorite (voir docstring module)."""
    label_u = label.strip().upper()
    if label_u == "NRM":
        return "N", 0.0
    m = re.match(r"^([A-Za-z]+)(" + _NUM + r")$", label_u)
    if not m:
        return None
    prefix, num = m.group(1), float(m.group(2))
    if prefix.startswith("T"):
        return "D", num
    if prefix.startswith("A"):
        return "F", num
    return None


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
    # Heuristique de secours (voir _decode_step_label / docstring module) :
    # seulement pour les paliers NUS (sans prefixe de lettre) - un bloc
    # entierement prefixe (variante CryoMag/.pmd) n'a jamais besoin de ce
    # repli, chaque palier etant deja type par son propre label.
    bare_numeric_steps = []
    for st in sp["steps"]:
        if _decode_step_label(st["step_label"]) is None:
            try:
                bare_numeric_steps.append(float(st["step_label"]))
            except ValueError:
                pass
    cod1_demag_fallback = "F" if _is_af(bare_numeric_steps) else "D"

    rows = [_MEAS_HEADER]
    for j, st in enumerate(sp["steps"]):
        decoded = _decode_step_label(st["step_label"])
        if decoded is not None:
            cod1, etape = decoded
        elif j == 0:
            cod1, etape = "N", 0.0
        else:
            try:
                step_val = float(st["step_label"])
            except ValueError:
                step_val = 0.0
            cod1, etape = cod1_demag_fallback, step_val

        fake_m = _FakeMeasurement(etape=etape, cod1=cod1, cod2="0")
        codes, _temp_k, af_field, _dc, _phi, _theta = _measurement_treatment(fake_m, [], 0.0)

        af_field_mT = af_field * 1.0e3 if cod1 == "F" else None
        temp_c = etape if cod1 in ("N", "D") else 0.0

        row = "\t".join([
            f"{etape:6.1f}", cod1, "0",
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
