"""Export Google Earth (.kml) des sites du fichier .prmag charge - demande
explicite utilisateur ("in Starpaleomag, is it possible to add a menu
create kml file from prmag", sous "export results to Stereo_Py").

Une Placemark par site (`magic_site`, repli sur les 6 premiers caracteres
de l'id - meme convention que calcul.site_of_result), au format de
description demande : Nb pmag cores / date / strike dip / Nb specimens
demagnetized / Lithology / Geologic type / Age / Location. Lit
directement `donnees` (List[Pmag], le fichier ENTIER charge) plutot que la
selection courante - une carte de sites n'a pas de raison de dependre de
ce qui est selectionne pour les graphiques."""

from typing import Dict, List
from xml.sax.saxutils import escape as _xml_escape

from testlect import Pmag

_DEFAULT_STYLE_URL = "root://styleMaps#default+nicon=0x304+hicon=0x314"
# Codes de desaimantation AF / thermique, y compris les etapes Thellier
# (R/V/P) - demande explicite utilisateur ; pas N/I (NRM seule, IRM).
_DEMAG_CODES = {"A", "F", "D", "S", "T", "K", "R", "V", "P"}


def _sample_name(p: Pmag) -> str:
    """Nom de carotte (= echantillon) : `magic_sample` si renseigne, sinon
    l'id de specimen sans sa lettre finale A-E (convention specimen =
    echantillon + lettre de sous-carotte, voir detailed_export)."""
    if p.magic_sample.strip():
        return p.magic_sample.strip()
    sid = p.id.strip()
    return sid[:-1] if sid and sid[-1] in "ABCDE" else sid


def _is_demagnetized(p: Pmag) -> bool:
    return sum(1 for m in p.mesures if m.cod1 in _DEMAG_CODES) >= 2


def _description_html(specimens: List[Pmag]) -> str:
    first = specimens[0]
    cores = {_sample_name(p) for p in specimens}
    n_demag = sum(1 for p in specimens if _is_demagnetized(p))
    parts = [f"<B>Nb pmag cores : </B>{len(cores)}<BR>"]
    if first.year:
        parts.append(f"<B>date : </B>{first.year:04d}-{first.month:02d}-{first.day:02d}<BR>")
    parts.append(f"<B>Strike : </B>{first.str_:.1f} <B>dip : </B>{first.dip:.1f}<BR>")
    parts.append(f"<B>Nb specimens demagnetized : </B>{n_demag}<BR>")
    for label, value in (
        ("Lithology", first.magic_li), ("Geologic type", first.magic_smt),
        ("Age", first.magic_age), ("Location", first.magic_loc),
    ):
        if value.strip():
            parts.append(f"<B>{label} : </B>{_xml_escape(value.strip())}<BR>")
    return "".join(parts)


def export_kml(donnees: List[Pmag], filepath: str, doc_name: str = "Sites") -> Dict[str, int]:
    """Ecrit `filepath` (.kml). Retourne {"sites": n ecrits, "no_coords":
    n sites ignores faute de latitude/longitude (0/0 = non renseigne)}."""
    by_site: Dict[str, List[Pmag]] = {}
    for p in donnees:
        by_site.setdefault(p.magic_site.strip() or p.id.strip()[:6], []).append(p)

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<kml xmlns="http://www.opengis.net/kml/2.2">',
        "<Document>",
        f"<name>{_xml_escape(doc_name)}</name>",
    ]
    written = no_coords = 0
    for site in sorted(by_site):
        specimens = by_site[site]
        lat = sum(p.lat for p in specimens) / len(specimens)
        lon = sum(p.rlong for p in specimens) / len(specimens)
        if lat == 0.0 and lon == 0.0:
            no_coords += 1
            continue
        lines += [
            "<Placemark>",
            f"<name>{_xml_escape(site)}</name>",
            f"<description><![CDATA[{_description_html(specimens)}]]></description>",
            "<open>1</open>",
            f"<styleUrl>{_DEFAULT_STYLE_URL}</styleUrl>",
            f"<Point><coordinates>{lon:.6f},{lat:.6f},0</coordinates></Point>",
            "</Placemark>",
        ]
        written += 1
    lines += ["</Document>", "</kml>"]
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return {"sites": written, "no_coords": no_coords}
