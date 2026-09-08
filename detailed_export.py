"""
Port de `exportpmagren` ("export detailed Rennes", dataselect.f) : fichier
texte unique avec, pour chaque echantillon, un bloc de parametres complet
puis le tableau de mesures.

`exporttolatex` (l'export Latex compagnon dans le Fortran d'origine, que
`exportpmagren` enchainait automatiquement) a ete PORTE puis RETIRE -
demande explicite utilisateur ("we can remove the Latex export... a
paleomagnetist who want to see the data can do it through
STARpaleomag_Py") : necessitait une chaine LaTeX installee juste pour
produire un PDF, un cout injustifie des lors que l'application elle-meme
reste le bon outil de consultation. export_detailed_txt genere desormais
AUTOMATIQUEMENT une copie PDF monospace du meme texte a la place (voir
_render_text_as_pdf) - demande explicite utilisateur ("a PDF copy of the
text export: do it automatically with the export text").

Hors perimetre (documente, pas un oubli) :
- Declinaison IGRF (`decli_igrf`/`declin`, via `orient_sample`/
  IGRFstarmac.f) : IGRF n'est pas porte ailleurs dans ce projet
  (selection.py le note deja) - ecrit "n.d" comme le fait le Fortran
  lui-meme pour ses propres cas de donnees manquantes.

Les champs Site/Sample/Fm/Age/GC/SMT/Li/Loc viennent directement des
champs magic_* deja decodes depuis la ligne roche (testlect.decode_roche)
- pas besoin de la reparser ici.
"""

import os
from typing import Dict, List, Optional

from selection import SelectedSample, Measurement, polere, corfor, corpen

_HEADER_COMMENT = [
    "! Paleomagnetic laboratory - Geosciences Rennes",
    "! codes for the Equipments when available 2G magnetometer - Molspin spinner - Agico JR6",
    "! 2G one measurement between two zeros  : code C1",
    "! 2G four measurements between two zeros : code C4",
    "! Molspin spinner 6 positions : code Mo",
    "! Jr6a spinner automated mode : code Ja",
    "! Jr6/Jr5 spinner 2 positions : code J2",
    "! thermal demagnetization: D+ and D- = sample orientation changed in the furnace between two steps along +Z and -Z",
    "! AF 2G Online three axis F+ = sequence coils  X,Z,Y, F- = sequence coils  Y,Z,X",
    "! FX = along X; FY = along Y, FZ = along Z;  FG = combine FX,FY and FZ to remove GRM",
    "! Paleointensities with the Coe/Tauxe or Thellier method",
    "! R: indicates field along Z axis; V: indicates field along -Z ; P: Ptrm check",
]


def _header_lines(source_file: Optional[str]) -> List[str]:
    """`_HEADER_COMMENT` + une ligne "data from file: ..." quand
    `source_file` est fourni - demande explicite utilisateur ("as we are
    now importing data from other files, we should adapt the header and
    indicate the name of the file") : le bloc d'origine ("Equipments 2G
    magnetometer - Molspin spinner - Agico JR6", "MMTD Furnace...")
    supposait implicitement une acquisition faite a Rennes ; les donnees
    exportees peuvent maintenant venir d'une contribution MagIC ou d'un
    autre laboratoire (import Utrecht, etc.) sans ces equipements -
    "codes for the Equipments when available" et le nom du fichier source
    rendent l'en-tete correct dans les deux cas plutot que de sous-
    entendre a tort une acquisition Rennes."""
    lines = list(_HEADER_COMMENT)
    if source_file:
        lines.append(f"! data from file: {source_file}")
    return lines


def _sample_display_name(ech: SelectedSample) -> str:
    """Nom d'echantillon = specimen id sans le dernier caractere A-E s'il
    en a un (convention specimen = sample + lettre de sous-carotte)."""
    sid = ech.id.strip()
    if sid and sid[-1] in "ABCDE":
        return sid[:-1]
    return sid


def _dc_field_string(m: Measurement, prev: List[Measurement], rfield: float, thellier: bool) -> str:
    """Equivalent du bloc de 8 `if` construisant `dc_field` (dataselect.f,
    juste avant l'ecriture de chaque ligne de mesure).

    Deux corrections par rapport a la premiere version du port :

    1) `dc_field` est declare `CHARACTER*9` en Fortran (dataselect.f:1853)
       - CHAQUE `write(dc_field, ...)`, quelle que soit la branche, aboutit
       donc TOUJOURS a exactement 9 caracteres (complete/tronque
       automatiquement). Les gabarits Python precedents ("  0:0:%d " etc.)
       n'imposaient PAS cette largeur fixe - une valeur a 3 chiffres, par
       exemple, produisait une chaine plus longue que prevu, decalant
       toutes les colonnes suivantes (Mag, Dsc, Isc...) sur cette ligne
       precise - demande explicite utilisateur ("is it also possible to
       align the text for the data in the listing"). Construit maintenant
       le contenu "x:y:z" puis le force a 9 caracteres (`:>9.9s`), quelle
       que soit la branche.

    2) Utilise PRIORITAIREMENT m.treat_dc_field (colonne treat_dc_lowfield
       des .prmag - par MESURE, voir testlect.Measurement/read_prmag_file)
       quand elle est renseignee, plutot que le seul `rfield` reconstruit
       depuis le commentaire (ech.com[:2] - convention Fortran d'origine,
       qui ne lit jamais de champ reel par mesure et reste le seul repli
       pour les fichiers .ren historiques sans cette colonne) - demande
       explicite utilisateur ("the dc field value is not exported") :
       pour un fichier issu d'une contribution MagIC, ech.com n'encode
       rien d'utilisable et rfield restait a 0.0 sur toutes les lignes,
       alors que la vraie valeur EST deja lue par read_prmag_file, juste
       jamais utilisee dans cet export."""
    effective = m.treat_dc_field if m.treat_dc_field is not None else rfield

    def axis(x: str = "0", y: str = "0", z: str = "0") -> str:
        return f"{x}:{y}:{z}"

    if m.cod1 == "R":
        content = axis(z=str(int(effective)))
    elif m.cod1 == "V":
        content = axis(z=str(int(-effective)))
    elif m.cod1 == "P":
        if thellier:
            content = axis(z=str(int(effective)))
        elif prev and prev[-1].cod1 == "S":
            content = axis(z=str(int(effective)))
        else:
            content = axis()
    elif m.cod1 == "Z" and m.cod2 == "+":
        content = axis(z=str(int(effective)))
    elif m.cod1 == "Z" and m.cod2 == "-":
        content = axis(z=str(int(-effective)))
    elif m.cod1 == "Y" and m.cod2 == "+":
        content = axis(y=str(int(effective)))
    elif m.cod1 == "Y" and m.cod2 == "-":
        content = axis(y=str(int(-effective)))
    elif m.cod1 == "X" and m.cod2 == "+":
        content = axis(x=str(int(effective)))
    elif m.cod1 == "X" and m.cod2 == "-":
        content = axis(x=str(int(-effective)))
    else:
        content = axis()
    return f"{content:>9.9s}"


def _sample_param_lines(ech: SelectedSample, depthsam: Optional[float]) -> List[str]:
    """Bloc de parametres (str2..str20 du Fortran)."""
    lines = [" --------------  Parameters sample & data   ---------------- ", ""]
    lines.append(f"Site     :  {ech.magic_site}")
    lines.append(f"Sample   :  {_sample_display_name(ech)}")
    lines.append(f"Specimen :  {ech.id}")
    if ech.norme == "v":
        lines.append(f"Volume   :{ech.vol:8.3f}     masse :   n.d ")
    else:
        lines.append(f"Volume   :   n.d      masse :  {ech.vol:8.3f}")

    if depthsam is not None:
        lines.append(f"Lat :{ech.lat:10.5f}    Long :{ech.rlong:12.5f}"
                      f"  height magnetostratigraphy :{depthsam:7.1f}")
    else:
        lines.append(f"Lat :{ech.lat:10.5f}    Long :{ech.rlong:12.5f}"
                      f"     Elevation :{ech.altitude:7.1f}")

    lines.append(f"Sampling date     :  {int(ech.year):4d}   {int(ech.month):2d}   {int(ech.day):2d}")
    lines.append(f"Sampling time UTM :  {int(ech.hour):2d}  {int(ech.minute):2d}")

    lines.append(f"azimuth mag :{ech.azmag:6.1f}  IGRF  Declination :   n.d")
    if ech.azsun == 0.0 and ech.hour == 0.0:
        lines.append(f"azimuth sun :{ech.azsun:6.1f}  Local Declination :    n.d ")
    else:
        lines.append(f"azimuth sun :{ech.azsun:6.1f}  Local Declination :   n.d")

    lines.append('Orientation :  "use AGICO code A12_0_3_90"')
    lines.append(f"core azimuth   :  {ech.caz:6.1f}")
    lines.append(f"core dip       :  {ech.cin:6.1f}")
    lines.append(f"Strike bedding :  {ech.str_:6.1f}")
    lines.append(f"Dip bedding    :  {ech.dip:6.1f}")

    lines.append(f'Formation   :   "{ech.magic_fm}"')
    lines.append(f'Age         :   "{ech.magic_age}"')
    lines.append(f'Geology     :   "{ech.magic_gc} : {ech.magic_smt} : {ech.magic_li}"')
    lines.append(f'Locality    :   "{ech.magic_loc}"')
    lines.append(f'Observation :   "{ech.magic_obs}"')

    rfield = 0.0
    try:
        rfield = float(ech.com[:2])
    except (ValueError, TypeError):
        pass
    if rfield != 0.0:
        lines.append(f"dc applied magnetic field : {rfield:4.1f}  µT")
    else:
        lines.append(f"dc applied magnetic field : {rfield:4.1f}  µT")
    return lines


def _measurement_table_lines(ech: SelectedSample) -> List[str]:
    """Equivalent du tableau `Step code dc_field Mag ... Dsc Isc Dis Iis
    Dtc Itc q ins K` (dataselect.f)."""
    if ech.norme == "m":
        lines = [" Step code  dc_Field    Mag(Am2)     Am2/kg    Dsc   Isc     "
                 "Dis   Iis     Dtc   Itc    q  Mag     K"]
    else:
        lines = [" Step code  dc_field    Mag(Am2)      A/m      Dsc   Isc     "
                 "Dis   Iis     Dtc   Itc    q  Mag     K"]

    rfield = 0.0
    try:
        rfield = float(ech.com[:2])
    except (ValueError, TypeError):
        pass
    thellier = any(m.cod1 == "V" for m in ech.mesures)

    for j, m in enumerate(ech.mesures):
        dc_field = _dc_field_string(m, ech.mesures[:j], rfield, thellier)

        mag1, dec1, inc1 = polere(m.x, m.y, m.z)
        x2, y2, z2 = corfor(m.x, m.y, m.z, ech.cin, ech.caz)
        _mag2, dec2, inc2 = polere(x2, y2, z2)
        x3, y3, z3 = corpen(x2, y2, z2, ech.dip, ech.str_)
        _mag3, dec3, inc3 = polere(x3, y3, z3)

        # m.etape est deja la valeur physique reelle (mT pour AF, degC
        # sinon) - plus d'echelle Oersted a compenser ici.
        step_mag = float(m.etape)
        if ech.norme == "m":
            rxx = mag1 * 1.0e3 / ech.vol if ech.vol else 0.0
            k = m.s * 1.0e-7 / ech.vol if ech.vol else 0.0
        else:
            rxx = mag1 * 1.0e6 / ech.vol if ech.vol else 0.0
            k = m.s * 10.0 * 1.0e-5 / ech.vol if ech.vol else 0.0

        lines.append(
            f"{step_mag:6.1f} {m.cod1}{m.cod2}  {dc_field}  {mag1:10.3E}  {rxx:10.3E}"
            f"  {dec1:6.1f}{inc1:6.1f}  {dec2:6.1f}{inc2:6.1f}  {dec3:6.1f}{inc3:6.1f}"
            f"{m.q:4d}  {m.ins:<2s}  {k:10.3E}"
        )
    return lines


# ---------------------------------------------------------------------------
# export detailed Rennes (exportpmagren)
# ---------------------------------------------------------------------------

def _render_text_as_pdf(lines: List[str], out_path: str) -> None:
    """Copie PDF du meme texte, police MONOSPACE - demande explicite
    utilisateur ("a PDF copy of the text export: do it automatically
    with the export text"), suite a "perhaps a PDF copy of the text file
    might be useful as text editor and fonts may change with different
    computers" : un .txt ouvert dans un editeur quelconque peut perdre
    l'alignement en colonnes si la police par defaut n'est pas a chasse
    fixe - le PDF fige une police monospace, l'alignement reste garanti
    partout. Remplace l'export Latex (retire - necessitait une chaine
    LaTeX installee juste pour obtenir un PDF, un cout que
    STARpaleomag_Py lui-meme rend inutile pour consulter les donnees).

    matplotlib (deja une dependance du projet, utilisee partout ailleurs
    pour les graphiques) plutot qu'une nouvelle dependance PDF dediee.
    Page US Letter PAYSAGE (les lignes de donnees font ~100-130
    caracteres, plus larges que hautes) ; taille de police calculee pour
    que la ligne la plus longue tienne sur la largeur de page, puis
    paginee en consequence."""
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.figure import Figure

    # chr(12) (form feed) prefixe certaines lignes dans le .txt - une
    # convention d'imprimante (saut de page) reprise du Fortran, qui n'a
    # pas de glyphe dans une police normale (avertissement matplotlib
    # "Glyph 12 missing") : la pagination du PDF est deja geree ci-dessous
    # independamment, ce caractere ne sert plus a rien ici.
    lines = [l.replace("\x0c", "") for l in lines]

    page_w, page_h = 11.0, 8.5  # pouces, paysage
    margin = 0.4
    usable_w_pt = (page_w - 2 * margin) * 72.0
    usable_h_pt = (page_h - 2 * margin) * 72.0

    max_len = max((len(l) for l in lines), default=1) or 1
    # largeur moyenne d'un caractere monospace ~ 0.60 * taille de police
    font_size = max(5.0, min(9.0, usable_w_pt / (max_len * 0.60)))
    linespacing = 1.15
    line_height_pt = font_size * linespacing
    lines_per_page = max(1, int(usable_h_pt / line_height_pt))

    with PdfPages(out_path) as pdf:
        for start in range(0, len(lines), lines_per_page):
            page_lines = lines[start:start + lines_per_page]
            fig = Figure(figsize=(page_w, page_h))
            fig.text(
                margin / page_w, 1.0 - margin / page_h, "\n".join(page_lines),
                family="monospace", fontsize=font_size, va="top", ha="left",
                linespacing=linespacing,
            )
            pdf.savefig(fig)


def export_detailed_txt(
    samples: List[SelectedSample],
    location: str,
    out_path: str,
    heights: Optional[Dict[str, float]] = None,
    source_file: Optional[str] = None,
) -> None:
    lines = _header_lines(source_file)
    lines += ["", "", f"Location :{location}", ""]

    for ech in samples:
        if len(ech.mesures) < 2:
            continue
        depthsam = heights.get(ech.id.strip()) if heights else None
        lines.append("")
        lines.append(chr(12) + " --------------  Parameters sample & data   ---------------- ")
        lines.append("")
        lines.extend(_sample_param_lines(ech, depthsam)[1:])  # sans re-repeter le titre
        lines.append("")
        lines.extend(_measurement_table_lines(ech))

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    pdf_path = os.path.splitext(out_path)[0] + ".pdf"
    _render_text_as_pdf(lines, pdf_path)
