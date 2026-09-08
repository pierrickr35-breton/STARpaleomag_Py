"""
Archive de NOUVELLES mesures (acquises APRES la creation d'un .prmag, sur
le meme instrument de labo) dans ce .prmag deja existant - demande
explicite utilisateur ("After the creation of the .prmag file, it is
possible that some new data will be acquired on the magnetometer in the
lab. At this stage, I will not be able to rewrite the data acquisition
software... So I need the possibility to archive new data acquired in
the legacy files. I need also to upload those acquire with the JR6
magnetometer").

Deux sources, meme logique de fusion :
- `parse_legacy_new_measurements` : nouvelle acquisition dans l'ancien
  format Rennes ("Id:"/mesures, 1/2/3 lignes d'entete) - reutilise
  testlect.read_ren_file (deja tolerant a ces variantes).
- `parse_jr6_file` : fichier JR6 brut, colonnes fixes - port de la partie
  parsing de reference/ImportJR6/ImportJR6data.f95 (branches
  `importinoldren`/`importinoldtxt`, PAS `createnewren` - cod2='=' pour
  un pas "AD", pas 'T').

ECART DELIBERE par rapport a reference/ImportJR6/importinpmagren.f
(confirme par l'utilisateur, "we assume that data will be archived only
for specimens already defined in the .prmag file") : un specimen de la
nouvelle acquisition SANS correspondance dans le .prmag est ECARTE et
journalise, PAS cree comme le faisait le Fortran d'origine
(`if(ntkk==0) pmag(nb_ech+1)=pmagmini(it)`).

Deduplication (comme le Fortran, meme triplet de comparaison) : une
nouvelle mesure dont (etape, cod1, cod2) correspond DEJA a une mesure
existante du specimen est IGNOREE - permet de relancer l'import sur un
fichier qui recouvre partiellement des donnees deja archivees sans
dupliquer, y compris a l'INTERIEUR du nouveau lot lui-meme.

Patch le .prmag EN PLACE, en INSERANT les nouvelles lignes de mesure a la
fin du bloc du specimen concerne (les 4 lignes d'entete et les mesures
DEJA presentes restent OCTET PRES identiques - seules des lignes sont
AJOUTEES, jamais reecrites) - meme philosophie que
complete_sample_info.py. Les nouvelles lignes sont formattees par
convert_ren_to_r._measurement_rows (memes conventions treat_*/
method_codes/instrument_codes que le reste du fichier - PAS reimplementees
ici), appliquee a un Pmag temporaire qui inclut aussi les mesures DEJA
presentes (necessaire pour un contexte correct sur les pas R/V/P, dont la
classification depend de la mesure PRECEDENTE) - seules les lignes
correspondant aux nouvelles mesures sont effectivement inserees. Une
sauvegarde `.bak` est ecrite avant toute modification, une seule fois."""

import os
from typing import Dict, List, Tuple

from testlect import Measurement, Pmag, _prmag_kv_line, read_ren_file, read_prmag_file
from convert_ren_to_r import _measurement_rows

_JR6_UNRECOGNIZED = ("?", "?", 999)


def _jr6_step_code(truc1: str) -> Tuple[str, str, float]:
    """Port de la classification de code d'etape JR6 (ImportJR6data.f95,
    branches `importinoldren`/`importinoldtxt`) : `truc1` (8 caracteres)
    ex. "NRM     ", "AD10.0  ", "TD300   ", "ARM100  ", "IRM1000 ".
    "AD"/"ARM" (champs AF/ARM) sont deja en mT reel dans le fichier JR6 -
    plus multiplies par 10 (l'ancienne echelle Oersted-equivalente de
    Measurement.etape a ete supprimee, demande explicite utilisateur
    "convert all step integer to float")."""
    t = truc1.strip().upper()
    if t.startswith("NRM"):
        return "N", "O", 0.0
    if t.startswith("AD"):
        try:
            step = float(t[2:])
        except ValueError:
            return _JR6_UNRECOGNIZED
        return "F", "=", step
    if t.startswith("TD"):
        try:
            step = float(t[2:])
        except ValueError:
            return _JR6_UNRECOGNIZED
        return "D", "=", round(step)
    if t.startswith("ARM"):
        try:
            step = float(t[3:])
        except ValueError:
            return _JR6_UNRECOGNIZED
        return "A", "Z", step
    if t.startswith("IRM"):
        try:
            step = float(t[3:])
        except ValueError:
            return _JR6_UNRECOGNIZED
        return "I", "Z", round(step)
    return _JR6_UNRECOGNIZED


def parse_jr6_file(path: str, encoding: str = "latin-1") -> Dict[str, List[Measurement]]:
    """Parse un fichier JR6 brut (colonnes fixes : specimen[0:10],
    code d'etape[10:18] (8 car.), x/y/z (6 car. chacun), exposant
    (4 car.)) en dict specimen -> liste de Measurement. "q"
    (qualite/erreur) et "s" (susceptibilite) ne sont jamais reellement
    renseignes par le Fortran d'origine pour ce chemin d'import - mis a 0
    ici plutot que de reproduire un comportement indetermine."""
    by_specimen: Dict[str, List[Measurement]] = {}
    with open(path, "r", encoding=encoding, errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if len(line) < 40:
                continue
            specimen = line[0:10].strip().upper().replace(" ", "_")
            if not specimen:
                continue
            truc1 = line[10:18]
            try:
                x = float(line[18:24])
                y = float(line[24:30])
                z = float(line[30:36])
                iexpo = int(line[36:40])
            except ValueError:
                continue
            cod1, cod2, etape = _jr6_step_code(truc1)
            r = 10.0 ** (iexpo - 5)
            meas = Measurement(
                etape=etape, cod1=cod1, cod2=cod2,
                x=x * r, y=y * r, z=z * r, q=0, ins="J6", s=0.0,
            )
            by_specimen.setdefault(specimen, []).append(meas)
    return by_specimen


def parse_legacy_new_measurements(path: str, encoding: str = "latin-1") -> Dict[str, List[Measurement]]:
    """Nouvelle acquisition dans l'ancien format Rennes (meme instrument
    de labo, logiciel d'acquisition non reecrit - demande explicite
    utilisateur) - reutilise read_ren_file (tolerant a 1/2/3 lignes
    d'entete par specimen). Seules les mesures sont conservees, le reste
    (site/sample/lat/lon...) est ignore : on archive UNIQUEMENT dans des
    specimens DEJA existants du .prmag cible."""
    samples = read_ren_file(path, encoding=encoding)
    return {p.id: p.mesures for p in samples if p.mesures}


def _measurement_key(m: Measurement) -> tuple:
    return (m.etape, m.cod1, m.cod2)


def archive_new_measurements(
    prmag_path: str, new_by_specimen: Dict[str, List[Measurement]],
) -> Tuple[int, int, int, List[str]]:
    """Fusionne `new_by_specimen` (specimen -> nouvelles mesures) dans un
    .prmag EXISTANT - voir docstring module pour le detail complet
    (insertion en fin de bloc, deduplication, specimens non trouves
    ecartes plutot que crees).

    Retourne (nb_specimens_completes, nb_mesures_ajoutees,
    nb_mesures_deja_presentes_ignorees, liste_specimens_sans_correspondance)."""
    existing_by_id = {p.id: p for p in read_prmag_file(prmag_path)}

    with open(prmag_path, "r", encoding="utf-8") as f:
        lines = [raw.rstrip("\n") for raw in f]
    original_text = "\n".join(lines) + "\n"
    n = len(lines)
    i = 0

    def skip_blank():
        nonlocal i
        while i < n and (lines[i].strip() == "" or lines[i].lstrip().startswith("#")):
            i += 1

    insertions: Dict[int, List[str]] = {}
    n_specimens_done = 0
    n_added = 0
    n_dup = 0
    seen_specimens = set()

    skip_blank()
    while i < n:
        if i + 4 >= n:
            break
        line_a = _prmag_kv_line(lines[i])
        specimen = line_a.get("specimen", "").strip()
        i += 4  # line_a/b/c/d
        i += 1  # ligne d'en-tete des mesures
        block_end = i
        while block_end < n and lines[block_end].strip() != "":
            block_end += 1
        i = block_end

        seen_specimens.add(specimen)
        new_meas = new_by_specimen.get(specimen)
        if new_meas:
            existing_p = existing_by_id.get(specimen)
            existing_mesures = existing_p.mesures if existing_p else []
            existing_keys = {_measurement_key(m) for m in existing_mesures}

            deduped: List[Measurement] = []
            seen_new = set()
            for m in new_meas:
                key = _measurement_key(m)
                if key in existing_keys or key in seen_new:
                    n_dup += 1
                    continue
                seen_new.add(key)
                deduped.append(m)

            if deduped:
                temp = Pmag(id=specimen, com=(existing_p.com if existing_p else ""),
                            mesures=list(existing_mesures) + deduped)
                rows = _measurement_rows(temp)[1:]  # [0] = ligne d'en-tete des colonnes
                insertions[block_end] = rows[-len(deduped):]
                n_added += len(deduped)
                n_specimens_done += 1
        skip_blank()

    unmatched = sorted(sp for sp in new_by_specimen if sp not in seen_specimens)

    if insertions:
        out_lines: List[str] = []
        for idx, line in enumerate(lines):
            if idx in insertions:
                out_lines.extend(insertions[idx])
            out_lines.append(line)
        if n in insertions:
            out_lines.extend(insertions[n])

        backup_path = prmag_path + ".bak"
        if not os.path.exists(backup_path):
            with open(backup_path, "w", encoding="utf-8") as f:
                f.write(original_text)
        with open(prmag_path, "w", encoding="utf-8") as f:
            f.write("\n".join(out_lines) + "\n")

    return n_specimens_done, n_added, n_dup, unmatched
