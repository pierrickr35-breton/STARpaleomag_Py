"""
Convertit un fichier ".col" du groupe d'Utrecht (paleomagnetism.org /
PMAG2, https://github.com/Jollyfant/PMAG2 - malgre l'extension, c'est du
JSON) vers le nouveau format .prmag - demande explicite utilisateur
("import data from the Utrecht group... Will it be possible to write a
converter for these data?"). Chaque specimen porte deja des pas de
demagnetisation resolus en x/y/z (repere "specimen", pas de calibration
Rennes a inverser) ET, le cas echeant, une ou plusieurs interpretations
(ajustements de droite/plan deja calcules par leur logiciel) - importees
ici vers un fichier .pmagres compagnon ("import the interpretation in
Pmagres").

Notes :

1) Site : ces fichiers ne portent AUCUNE notion de site (le groupe
   d'Utrecht a abandonne cette notion : "one sample one site" - position
   du groupe, PAS un choix de ce portage), et le nom du site ne peut pas
   etre extrait du nom de specimen (confirme par l'utilisateur : "it is
   not possible to extract a site name from the specimen name" - noms du
   type "BN4.1A"/"SR11.1I" sans separateur fiable). Convention retenue,
   demande explicite utilisateur ("assume that the site is the name of
   the col file") : `site` = nom du fichier .col source (sans extension),
   voir `convert_files`. Approximation assumee (un fichier peut regrouper
   plusieurs sites reels, ex. "SR1-3-4-10-11.col") - a raffiner
   manuellement dans STARpaleomag_Py au besoin (meme mecanisme que
   `magic_site`/`Site: "..."`, deja decouple du nom de specimen pour les
   imports MagIC "case 1").

2) Orientation carotte (coreAzimuth/coreDip) : une premiere version de ce
   convertisseur ecrivait ces valeurs telles quelles, faute de
   transformation validee (une recherche exhaustive - 6 permutations
   d'axes x 8 signes x 9 variantes d'angle - n'avait trouve aucune
   combinaison universelle sur les 36 specimens interpretes de BN.col).
   Convention donnee par l'utilisateur ("the coredip is usually
   90-dip_rennes and the azimuth is azimuth of the X axis as in Magic") :

       caz (STARpaleomag_Py) = coreAzimuth + 90   (meme regle que MagIC, voir
                                            DIFF_WITH_MAGIC)
       cin (STARpaleomag_Py) = 90 - coreDip

   Verifiee numeriquement sur les 156 interpretations reelles des 4
   fichiers fournis (BN/SR1-3-4-10-11/SR2/SR5-6-7-8-9) : `corfor(x,y,z,
   cin,caz)` applique au vecteur specimen de chaque interpretation
   reproduit le vecteur geographique deja calcule par PMAG2 avec un
   produit scalaire normalise > 0.9999 dans 155 cas sur 156 - le seul
   ecart est un antipode EXACT (dot=-1.0, pas une valeur intermediaire),
   symptome connu d'une ambiguite de signe sur l'eigenvecteur d'un
   ajustement de droite (PCA), pas un defaut de la formule. Les vues
   in-situ/apres pendage (orientation 2/3) sont donc desormais fiables
   pour les echantillons importes par ce convertisseur.

   bed_dip_strike/bed_dip (pendage de la couche, PAS l'orientation de la
   carotte) restent, comme avant, ecrits directement : PMAG2 utilise deja
   la meme convention `dipDirection = strike + 90` que celle documentee
   dans DIFF_WITH_MAGIC/le header .prmag.

3) Unite de x/y/z ("check the unit intensity of the magnetization... in
   prmag x,y,z are in Am2") : Utrecht stocke x/y/z DEJA divise par le
   volume (uA/m, une magnetisation), pas en Am2 (un moment) comme .prmag
   l'attend - confirme dans le code source de PMAG2
   (interpretation/js/importing.js:importUtrecht, "Step is in pico Am^2
   .. divide by sample volume to get uAm/m!") et par coherence physique
   (sans reconversion, la magnetisation recalculee par STARpaleomag_Py serait de
   l'ordre de 1e8 A/m, impossible pour une roche). Reconverti en Am2 ici
   (voir _measurement_rows : moment = uA/m * volume(cm3) * 1e-12) avant
   ecriture dans .prmag.
"""

import argparse
import json
import os
from typing import Dict, List, Optional, Tuple

from calcul import FitResult, archivres, results_path_for
from magic_export import _measurement_treatment
from selection import polere

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

# "alternating" -> 'F' (AF, meme code que convert_ren_to_r.py/Rennes) ;
# "thermal" -> 'D' (demag thermique generique, PAS 'S' qui est reserve au
# pas zero-field d'un protocole IZZI - absent de ces fichiers, demag
# simple uniquement, aucune interpretation de type paleointensite vue
# dans les 4 fichiers fournis).
_COD1_BY_DEMAG_TYPE = {"alternating": "F", "thermal": "D"}


def _nd(value, fmt: Optional[str] = None) -> str:
    if value is None:
        if fmt:
            width_digits = "".join(ch for ch in fmt.split(".")[0] if ch.isdigit())
            if width_digits:
                return "n.d".rjust(int(width_digits))
        return "n.d"
    return format(value, fmt) if fmt is not None else str(value)


class _FakeMeasurement:
    """Juste assez de champs pour reutiliser magic_export._measurement_treatment
    sans construire un vrai testlect.Measurement (pas de x/y/z/q/s
    pertinents a ce stade). `ins` DOIT rester present (meme None) : le
    cas cod1=='F' de _measurement_treatment lit `m.ins` inconditionnellement
    (`(m.ins or "").strip().upper()`) - AttributeError sinon des qu'un
    fichier Utrecht a demagnetizationType=="alternating" (bug reel trouve
    en portant convert_ipgp_to_r.py, qui reutilise ce meme patron)."""
    def __init__(self, etape: float, cod1: str, cod2: str = "0"):
        self.etape = etape
        self.cod1 = cod1
        self.cod2 = cod2
        self.ins = None


def _sample_header_block(sp: dict, site: str = "n.d") -> str:
    name = sp["name"]
    lat, lon = sp.get("latitude"), sp.get("longitude")
    age_min, age_max = sp.get("ageMin"), sp.get("ageMax")
    age = f"{_nd(age_min)} - {_nd(age_max)} Ma" if (age_min is not None or age_max is not None) else "n.d"

    line_a = (
        f"specimen: {name}\tsample: {name}\tsite: {site}\t"
        f"volume: {_nd(sp.get('volume'), '.2f')}\tmass: n.d\t"
        f"lat: {_nd(lat, '.5f')}\tlon: {_nd(lon, '.5f')}\televation: n.d\t"
        f"stratigraphic_height: n.d\t"
        f"comment: n.d"
    )
    # azimuth/dip : convention Utrecht -> STARpaleomag_Py donnee par l'utilisateur
    # et verifiee numeriquement sur les 156 interpretations reelles
    # fournies (voir docstring module) : caz=coreAzimuth+90 (meme regle
    # que MagIC), cin=90-coreDip.
    core_az, core_dip = sp.get("coreAzimuth"), sp.get("coreDip")
    caz = (core_az + 90.0) % 360.0 if core_az is not None else None
    cin = (90.0 - core_dip) if core_dip is not None else None
    line_b = (
        f"azimuth: {_nd(caz, '.1f')}\tdip: {_nd(cin, '.1f')}\t"
        f"date: n.d\tmagnetic_azimuth: n.d\tsolar_azimuth: n.d\torient_tool: n.d"
    )
    line_c = f"bed_dip_strike: {_nd(sp.get('beddingStrike'), '.1f')}\tbed_dip: {_nd(sp.get('beddingDip'), '.1f')}"
    line_d = (
        f"formation: n.d\tage: {age}\t"
        f"geologic_classes: {sp.get('geology') or 'n.d'}\tgeologic_types: n.d\t"
        f"lithologies: {sp.get('lithology') or 'n.d'}\tlocation: n.d\t"
        f"obs: n.d\tmethod_codes: n.d"
    )
    return "\n".join([line_a, line_b, line_c, line_d])


def _measurement_rows(sp: dict) -> List[str]:
    cod1_demag = _COD1_BY_DEMAG_TYPE.get(sp.get("demagnetizationType"), "D")
    # Utrecht stocke x/y/z DEJA divise par le volume de l'echantillon (voir
    # PMAG2 interpretation/js/importing.js:importUtrecht - "Step is in pico
    # Am^2 .. divide by sample volume to get uAm/m!"), donc en uA/m
    # (magnetisation), PAS en Am2 (moment) comme .prmag l'attend (voir
    # selection._mag_values : `mtot*1e6/vol(cm3)` = A/m, coherent
    # uniquement si x/y/z EST le moment Am2, pas deja la magnetisation).
    # Verification physique : sans cette conversion, une magnetisation
    # M=moment/volume appliquee UNE SECONDE FOIS sur des valeurs deja en
    # uA/m donnerait des M de l'ordre de 1e8 A/m (impossible pour une
    # roche - les valeurs reelles vont de ~1e-4 a quelques dizaines de
    # A/m) ; en repartant du moment Am2 correctement reconstitue on
    # retrouve un NRM de l'ordre de 1e-8 Am2, plausible pour un
    # specimen de ~10 cm3. Reconversion : moment(Am2) = M(uA/m) * 1e-6
    # (uA/m -> A/m) * volume(cm3) * 1e-6 (cm3 -> m3) = M * volume * 1e-12.
    volume_cm3 = sp.get("volume") or 0.0
    to_moment = volume_cm3 * 1.0e-12
    rows = [_MEAS_HEADER]
    for j, step in enumerate(sp["steps"]):
        cod1 = "N" if j == 0 else cod1_demag
        step_val = float(step["step"])
        # etape : valeur physique reelle directement (mT pour A/F, degC
        # sinon - plus d'echelle Oersted-equivalente, voir testlect.
        # Measurement/magic_export._measurement_treatment). Le NRM (j==0)
        # est toujours ecrit avec etape=0, quelle que soit la valeur brute
        # Utrecht (ex. thermal commence a "20", pas "0").
        etape = 0.0 if cod1 == "N" else step_val
        fake_m = _FakeMeasurement(etape=etape, cod1=cod1, cod2="0")
        codes, _temp_k, af_field, _dc, _phi, _theta = _measurement_treatment(fake_m, [], 0.0)

        af_field_mT = af_field * 1.0e3 if cod1 == "F" else None
        temp_c = etape if cod1 in ("N", "D") else 0.0
        step_str = 0.0 if cod1 == "N" else step_val
        x_moment, y_moment, z_moment = step["x"] * to_moment, step["y"] * to_moment, step["z"] * to_moment

        row = "\t".join([
            f"{step_str:6.1f}", cod1, "0",
            f"{x_moment:11.3E}", f"{y_moment:11.3E}", f"{z_moment:11.3E}",
            "n.d", "g", "n.d", "n.d",
            _nd(temp_c, "6.1f"),
            _nd(af_field_mT, "7.2f"),
            _nd(None, "7.2f"), _nd(0.0, "7.2f"),
            _nd(None, "6.1f"), _nd(None, "6.1f"),
            f"{codes:<40}", "n.d", str(j + 1),
        ])
        rows.append(row)
    return rows


def _fit_results_for_specimen(sp: dict) -> List[FitResult]:
    """Une interpretation Utrecht ("type": TAU1=droite, TAU3=plan) ->
    un FitResult, dec/inc pris directement du vecteur "specimen" (repere
    SPECIMEN, deja resolu par PMAG2) - PAS de corfor() ici, meme convention
    que le reste de .pmagres (dec/inc toujours stockees en repere
    specimen BRUT, corrigees seulement a l'affichage via
    calcul._correct_dec_inc + recompute_fit_geometry, qui utilise
    desormais le azimuth/dip ecrits par _sample_header_block - voir
    docstring module pour la conversion coreAzimuth/coreDip -> caz/cin).
    step_first/step_last dans la MEME convention interne que
    _measurement_rows (Oersted-equivalent pour l'AF), pour que
    calcul.recompute_fit_geometry (comparaison sur `m.etape`) retrouve
    bien les mesures correspondantes apres relecture du .prmag."""
    cod1_demag = _COD1_BY_DEMAG_TYPE.get(sp.get("demagnetizationType"), "D")
    demag = "F" if cod1_demag == "F" else "D"
    results = []
    for numcomp, it in enumerate(sp.get("interpretations", []), start=1):
        cat1 = "L" if it.get("type", "").startswith("TAU1") else "P"
        coords = it.get("specimen", {}).get("coordinates", {})
        if not coords:
            continue
        _mag, dec, inc = polere(coords.get("x", 0.0), coords.get("y", 0.0), coords.get("z", 0.0))

        step_labels = it.get("steps", [])
        if not step_labels:
            continue
        # step_vals (mT/degC reels, tels que fournis par Utrecht/PMAG2) se
        # comparent desormais directement a Measurement.etape (lui aussi
        # en unite reelle - plus d'echelle Oersted-equivalente a
        # appliquer ici pour un demag AF, voir testlect.Measurement).
        step_vals = [float(s) for s in step_labels]
        step_first_raw, step_last_raw = min(step_vals), max(step_vals)

        results.append(FitResult(
            id=sp["name"], cat1=cat1, cat2="",
            orig="o" if it.get("anchored") else "n",
            demag=demag, numcomp=numcomp, nb=len(step_labels),
            dec=dec, inc=inc, mad=it.get("MAD", 0.0),
            step_first=int(round(step_first_raw)),
            step_last=int(round(step_last_raw)),
        ))
    return results


def convert_files(paths_in: List[str], path_out: str) -> Tuple[int, int]:
    """Convertit UN OU PLUSIEURS fichiers .col vers un seul .prmag (+ un
    seul .pmagres compagnon) - demande explicite utilisateur ("import all
    the col files in a single .prmag file and assume that the site is the
    name of the col file"). Le site de chaque specimen est le nom du
    fichier .col source (sans extension) : ces fichiers n'ont aucune
    notion de site propre (voir docstring module), mais un fichier .col
    correspond en pratique a une collecte/collection reelle - un
    groupement plus utile que "n.d" partout, quitte a le raffiner
    manuellement ensuite (mecanisme `Site:`/`magic_site`, deja decouple
    du nom de specimen)."""
    all_specimens: List[Tuple[dict, str]] = []  # (specimen, site)
    for path_in in paths_in:
        with open(path_in, "r", encoding="utf-8") as f:
            data = json.load(f)
        site = os.path.splitext(os.path.basename(path_in))[0]
        all_specimens.extend((sp, site) for sp in data.get("specimens", []))

    source_names = ", ".join(os.path.basename(p) for p in paths_in)
    lines = [
        FORMAT_HEADER,
        f"#converted from {source_names} by convert_utrecht_to_r.py "
        "(Utrecht paleomagnetism.org/PMAG2 .col -> .prmag)",
        S_FIELD_HEADER,
        "",
    ]
    blocks = [
        _sample_header_block(sp, site) + "\n" + "\n".join(_measurement_rows(sp))
        for sp, site in all_specimens
    ]

    with open(path_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        f.write("\n\n".join(blocks) + "\n")

    path_results = results_path_for(path_out)
    if os.path.exists(path_results):
        # archivres() ajoute a un fichier existant plutot que de l'ecraser
        # (usage normal : accumulation interactive depuis l'appli) - une
        # re-conversion vers le MEME chemin de sortie dupliquerait donc
        # les interpretations si on ne repartait pas d'un fichier vide,
        # contrairement a .prmag (toujours ecrit "w", donc deja neuf a
        # chaque appel).
        os.remove(path_results)
    existing_ids: Optional[set] = None
    nb_results = 0
    for sp, _site in all_specimens:
        for res in _fit_results_for_specimen(sp):
            _c, existing_ids = archivres(res, path_results, existing_ids)
            nb_results += 1

    return len(blocks), nb_results


def convert_file(path_in: str, path_out: str) -> Tuple[int, int]:
    """Cas particulier a un seul fichier de `convert_files` - conserve
    pour compatibilite (usage CLI direct)."""
    return convert_files([path_in], path_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert one or more .col files (Utrecht/PMAG2) to a .prmag file")
    parser.add_argument("utrecht_files", nargs="+")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    out = args.output or os.path.splitext(args.utrecht_files[0])[0] + ".prmag"
    nb, nb_res = convert_files(args.utrecht_files, out)
    print(f"{nb} specimen(s) written to {out}")
    print(f"{nb_res} interpretation(s) written to {results_path_for(out)}")
