"""
Complete les metadonnees Site/Formation/Age/GC/SMT/Li/Loc/Obs (et lat/lon,
et en mode specimen, sample/site) d'un fichier .prmag DEJA converti,
depuis une table externe - demande explicite utilisateur, en 2 temps :
d'abord "we can build a routine to complete the sample information later
on" (pendant la discussion sur l'import "as-is" des fichiers legacy), puis
la specification concrete : "can we put a menu like complete sample
information with data in a table".

Trois modes, chacun avec ses propres colonnes de table (`#` optionnel en
tete de la ligne d'en-tete, comme les fichiers complement existants de
convert_legacy_ren.py) :

- Mode SITE (`complete_site_info`) : suppose que specimen/sample sont
  DEJA corrects dans le .prmag (typiquement le cas apres un import MagIC
  ou une conversion legacy qui a bien derive le site des 6 premiers
  caracteres du specimen) - la table ne fournit QUE les metadonnees
  complementaires, UNE LIGNE PAR SITE, appliquee a TOUS les specimens de
  ce site :
      #site	lat	lon	formation	age	geologic_classes	geologic_types	lithologies	location	obs
  La colonne 'site' EST la cle de recherche (nom deja present dans le
  .prmag) - mais peut etre RENOMMEE en repetant l'en-tete 'site' une
  seconde fois, cette seconde colonne fournissant alors le nouveau nom -
  demande explicite utilisateur ("i cannot select site and replace it
  by an other site"), format exact d'un fichier reel de l'utilisateur
  (Complement_tect.txt : ancien nom de campagne -> nouveau nom
  consolide) :
      #site	site	lat	lon	...
  Ce meme fichier fournit aussi l'age deja separe en age_low/age_high/
  age_unit plutot qu'une seule colonne 'age' (voir `_normalize_age`), et
  utilise des VIRGULES comme separateur plutot qu'une tabulation (voir
  `_sniff_delimiter` - les deux styles sont acceptes, detectes
  automatiquement depuis la ligne d'en-tete).

- Mode SPECIMEN (`complete_specimen_info`) : ne suppose RIEN de fiable sur
  sample/site dans le .prmag - la table fournit AUSSI sample/site
  (ecrases si presents), UNE LIGNE PAR SPECIMEN. Accepte aussi une colonne
  `stratigraphic_height` (metres, positif vers le haut - meme champ/unite
  que `samples.height` du modele MagIC) - demande explicite utilisateur
  ("add an additional variable: stratigraphic_position (or the Magic
  equivalent)... this field will replace the need to load a file for
  magnetostratigraphic studies") :
      #specimen	sample	site	lat	lon	stratigraphic_height	formation	age	geologic_classes	geologic_types	lithologies	location	obs

- Mode SAMPLE (`complete_sample_height`) : cas particulier, UNIQUEMENT
  pour `stratigraphic_height` - demande explicite utilisateur ("add a
  specific case for stratigraphic_height filled from sample and height
  as all specimens have the same height") : plusieurs specimens (A/B/C)
  decoupes du MEME echantillon physique partagent necessairement la
  meme position stratigraphique - repeter cette valeur par SPECIMEN
  (mode ci-dessus) est inutilement verbeux des qu'un sample porte
  plusieurs specimens ; ce mode applique UNE LIGNE PAR SAMPLE a tous les
  specimens qui le partagent :
      #sample	stratigraphic_height

`stratigraphic_height` varie par SAMPLE (voire par specimen si vraiment
necessaire), contrairement a lat/lon/formation qui ne varient pas au sein
d'un site - c'est pourquoi le mode SITE ne l'accepte jamais (un site
combine plusieurs samples a des hauteurs differentes).

Dans les trois modes, une colonne peut etre omise/laissee vide sans
consequence (seules les cles reellement fournies, non vides, sont
appliquees) - completer PROGRESSIVEMENT au fil de ce qui devient connu est
le but explicite de cette routine, pas un remplissage en un seul coup.

Patch le fichier .prmag EN PLACE au niveau du TEXTE BRUT des 4 lignes
d'entete par specimen (pas une reserialisation depuis les objets Pmag deja
parses par testlect.read_prmag_file) : seules les cles explicitement
fournies par la table sont remplacees, tout le reste (volume/mass/
elevation/comment, azimuth/dip/date, bed_dip, method_codes, et TOUTES les
lignes de mesure) reste OCTET PRES identique. Une reserialisation depuis
Pmag perdrait `method_codes` (jamais stocke dans le dataclass Pmag - voir
testlect.py) et risquerait de reformater des valeurs numeriques
differemment de l'original.

Une sauvegarde `.bak` du fichier .prmag est ecrite AVANT toute
modification, une seule fois (jamais ecrasee si elle existe deja - garde
le tout premier etat, avant la toute premiere completion)."""

import os
from typing import Dict, List, Optional, Tuple

_ROCHE_TABLE_FIELDS = (
    "formation", "age", "geologic_classes", "geologic_types", "lithologies", "location", "obs",
)


def _split_line(line: str, delimiter: str = "\t") -> List[str]:
    return [f.strip() for f in line.rstrip("\n").split(delimiter)]


def _sniff_delimiter(header_line: str) -> str:
    """Les tables de completion existantes (complement_Tibet.txt) sont
    tabulees, mais Complement_tect.txt (demande explicite utilisateur -
    correspondance ancien/nouveau nom de site) est separee par des
    virgules, avec des espaces de remplissage autour de chaque champ
    (ex. '96CC02     ,   C04     ,   -25.989   , ...') - deja absorbes
    par le .strip() de _split_line. Priorite a la tabulation si les deux
    sont presentes (jamais le cas en pratique, mais une valeur pourrait
    contenir une virgule)."""
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    return "\t"


def _read_table_rows(path: str, encoding: str) -> Tuple[List[str], List[List[str]]]:
    """Bugs reels corriges, decouverts sur une vraie table (complement_
    Tibet.txt, 100 sites, testee contre Tibet_14_15_Pmag_converted.prmag) :
    1) Les colonnes du header sont normalisees (espace -> underscore) - la
       table reelle ecrit "geologic classes"/"geologic types" (espace),
       jamais reconnues par _ROCHE_TABLE_FIELDS ("geologic_classes"/
       "geologic_types", underscore) : ces deux champs restaient
       silencieusement jamais appliques, meme pour les sites par ailleurs
       correctement mis a jour (formation/age/lithologies/location/obs).
    2) Les lignes plus courtes que le header (derniere(s) colonne(s) vide(s)
       sans tabulation finale - frequent pour "obs", souvent vide) sont
       COMPLETEES avec des chaines vides plutot que rejetees : avant ce
       fix, le site etait alors totalement ABSENT de `table`, donc
       rapporte a tort comme "sans correspondance" - verifie sur la table
       reelle : 72 des 100 lignes etaient dans ce cas, faisant echouer la
       mise a jour de 1662 des 2378 specimens du fichier de test."""
    with open(path, "r", encoding=encoding, errors="replace") as f:
        raw_lines = [l for l in f if l.strip()]
    if not raw_lines:
        return [], []
    delimiter = _sniff_delimiter(raw_lines[0])
    header = [h.lstrip("#").strip().lower().replace(" ", "_") for h in _split_line(raw_lines[0], delimiter)]
    rows = []
    for l in raw_lines[1:]:
        row = _split_line(l, delimiter)
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        rows.append(row)
    return header, rows


def _prmag_kv_line(line: str) -> dict:
    """Meme parsing que testlect._prmag_kv_line (non importe directement -
    module de bas niveau, pas de dependance croisee necessaire ici)."""
    result = {}
    for chunk in line.split("\t"):
        if ":" not in chunk:
            continue
        k, _sep, v = chunk.partition(":")
        result[k.strip()] = v.strip()
    return result


def _rewrite_kv_line(original_line: str, updates: Dict[str, str]) -> str:
    """Reconstruit une ligne 'cle: valeur\\tcle2: valeur2...' en ne
    changeant QUE les cles presentes dans `updates`, toutes les autres
    cles/valeurs de `original_line` restant identiques - y compris leur
    ORDRE et leur formatage d'origine (voir docstring module : jamais
    reconstruite depuis des champs Pmag deja parses)."""
    chunks = original_line.split("\t")
    new_chunks = []
    seen = set()
    for chunk in chunks:
        if ":" not in chunk:
            new_chunks.append(chunk)
            continue
        k, _sep, _v = chunk.partition(":")
        key = k.strip()
        if key in updates:
            new_chunks.append(f"{key}: {updates[key]}")
            seen.add(key)
        else:
            new_chunks.append(chunk)
    for key, value in updates.items():
        if key not in seen:
            new_chunks.append(f"{key}: {value}")
    return "\t".join(new_chunks)


_INFO_FIELDS = ("lat", "lon") + _ROCHE_TABLE_FIELDS


def _missing_fields(line_a: dict, line_d: dict, check_height: bool) -> List[str]:
    """Champs encore a "n.d" (placeholder ecrit par TOUS les convertisseurs
    .prmag quand une valeur est inconnue - convert_ren_to_r/
    convert_magic_to_r/convert_utrecht_to_r - voir leurs `_nd`/formattage
    de formation/age/.../lat/lon/stratigraphic_height) APRES application de
    la table de completion - demande explicite utilisateur ("it will be
    good to know what are the samples with missing information, when we
    complete the file"). `stratigraphic_height` seulement en mode
    specimen/sample (le mode site ne le touche jamais, cf. docstring
    module - le signaler comme "manquant" y serait du bruit puisque cette
    routine n'a jamais pour but de le remplir dans ce mode)."""
    missing = []
    for field in ("lat", "lon"):
        if line_a.get(field, "").strip().lower() in ("", "n.d"):
            missing.append(field)
    if check_height and line_a.get("stratigraphic_height", "").strip().lower() in ("", "n.d"):
        missing.append("stratigraphic_height")
    for field in _ROCHE_TABLE_FIELDS:
        if line_d.get(field, "").strip().lower() in ("", "n.d"):
            missing.append(field)
    return missing


def _normalize_age(record: dict) -> None:
    """Complement_tect.txt (et master_site_list_with_geology.csv dont il
    reprend la forme) fournit l'age SEPARE en colonnes age_low/age_high/
    age_unit (deja au vocabulaire Magic reel, ex. 'Ma') plutot qu'une
    seule colonne 'age' - reconstruit la chaine 'age_low - age_high
    age_unit' en place (meme forme que magic_export.parse_age sait
    re-decouper via son branchement '\\s-\\s' pour age_low/age_high, et
    la table _AGE_UNIT_PATTERNS pour age_unit) UNIQUEMENT si 'age' n'est
    pas deja fournie directement - ne change rien pour les tables
    existantes qui fournissent 'age' toute faite (mode specimen/sample)."""
    if record.get("age"):
        return
    low = record.get("age_low", "").strip()
    high = record.get("age_high", "").strip()
    if not (low and high):
        return
    unit = record.get("age_unit", "").strip()
    record["age"] = f"{low} - {high} {unit}".strip()


def _apply_table(
    prmag_path: str, table: Dict[str, dict], key_field: str,
) -> Tuple[int, List[str], List[Tuple[str, List[str]]]]:
    """Parcourt le .prmag bloc par specimen (meme detection que
    testlect.read_prmag_file), patche les lignes 'specimen:'/'formation:'
    de chaque bloc dont la cle (`site`, `specimen` ou `sample` selon
    `key_field`) est trouvee dans `table`. Retourne (nb_specimens_mis_a_
    jour, liste_specimens_sans_correspondance, liste_specimens_encore_
    incomplets) - ce dernier couvre AUSSI les specimens matches dont la
    table ne fournissait qu'une partie des champs (pas seulement les
    specimens sans aucune correspondance)."""
    with open(prmag_path, "r", encoding="utf-8") as f:
        lines = [raw.rstrip("\n") for raw in f]
    original_text = "\n".join(lines) + "\n"
    n = len(lines)
    i = 0
    n_updated = 0
    unmatched: List[str] = []
    still_missing: List[Tuple[str, List[str]]] = []

    def skip_blank():
        # meme comportement que testlect.read_prmag_file._skip_blank_and_comments
        # (lignes vides ET commentaires "#...") - un premier essai qui ne
        # sautait que les lignes vides desalignait la lecture des le
        # premier bloc (les 5 lignes d'en-tete "#..." du fichier .prmag
        # etaient alors lues a tort comme un faux bloc specimen).
        nonlocal i
        while i < n and (lines[i].strip() == "" or lines[i].lstrip().startswith("#")):
            i += 1

    skip_blank()
    while i < n:
        if i + 4 >= n:
            break
        idx_a = i
        line_a = _prmag_kv_line(lines[i]); i += 1
        i += 1  # line_b - jamais touchee
        i += 1  # line_c - jamais touchee
        idx_d = i
        line_d = _prmag_kv_line(lines[i]); i += 1
        i += 1  # ligne d'en-tete des mesures
        while i < n and lines[i].strip() != "":
            i += 1

        specimen = line_a.get("specimen", "").strip()
        site = line_a.get("site", "").strip()
        sample = line_a.get("sample", "").strip()
        lookup_key = {"specimen": specimen, "sample": sample}.get(key_field, site)
        record = table.get(lookup_key)
        if record is None:
            unmatched.append(specimen or f"(line {idx_a + 1})")
            missing = _missing_fields(
                _prmag_kv_line(lines[idx_a]), _prmag_kv_line(lines[idx_d]),
                check_height=(key_field in ("specimen", "sample")),
            )
            if missing:
                still_missing.append((specimen or f"(line {idx_a + 1})", missing))
            skip_blank()
            continue

        _normalize_age(record)
        a_updates = {}
        if record.get("lat"):
            try:
                a_updates["lat"] = f"{float(record['lat']):.5f}"
            except ValueError:
                pass
        if record.get("lon"):
            try:
                a_updates["lon"] = f"{float(record['lon']):.5f}"
            except ValueError:
                pass
        # stratigraphic_height varie specimen par specimen au sein d'une
        # meme section/site (c'est tout l'interet d'une etude
        # magnetostratigraphique) - contrairement a lat/lon/formation,
        # appliquer UNE valeur de table a tous les specimens d'un site
        # (mode site) ecraserait cette variation ; le mode specimen
        # l'applique (une ligne par specimen), et le mode sample aussi -
        # cas particulier explicite utilisateur ("add a specific case for
        # stratigraphic_height filled from sample and height as all
        # specimens have the same height") : plusieurs specimens (A/B/C)
        # decoupes du MEME echantillon physique partagent necessairement
        # la meme position - c'est d'ailleurs un champ SAMPLE, pas
        # specimen, dans le modele MagIC reel (samples.height, voir
        # convert_magic_to_r._sample_header_block) - seul le mode site
        # NE l'applique jamais (un site combine plusieurs samples a des
        # hauteurs differentes).
        if key_field in ("specimen", "sample") and record.get("stratigraphic_height"):
            try:
                a_updates["stratigraphic_height"] = f"{float(record['stratigraphic_height']):.2f}"
            except ValueError:
                pass
        if key_field == "specimen":
            if record.get("sample"):
                a_updates["sample"] = record["sample"]
            if record.get("site"):
                a_updates["site"] = record["site"]
        # Mode site : une SECONDE colonne "site" dans la table (deux
        # colonnes de meme nom d'en-tete, ex. Complement_tect.txt :
        # "site,site,Lat,Lon,...") est le nom vers lequel RENOMMER - la
        # premiere est la cle de recherche (nom deja present dans le
        # .prmag), demande explicite utilisateur ("i cannot select site
        # and replace it by an other site"). Pas de collision avec le
        # bloc specimen ci-dessus (mutuellement exclusif sur key_field).
        if key_field == "site" and record.get("site"):
            new_site = record["site"].strip()
            if new_site:
                a_updates["site"] = new_site
        if a_updates:
            lines[idx_a] = _rewrite_kv_line(lines[idx_a], a_updates)

        d_updates = {field: record[field] for field in _ROCHE_TABLE_FIELDS if record.get(field)}
        if d_updates:
            lines[idx_d] = _rewrite_kv_line(lines[idx_d], d_updates)

        if a_updates or d_updates:
            n_updated += 1

        missing = _missing_fields(
            _prmag_kv_line(lines[idx_a]), _prmag_kv_line(lines[idx_d]),
            check_height=(key_field in ("specimen", "sample")),
        )
        if missing:
            still_missing.append((specimen or f"(line {idx_a + 1})", missing))
        skip_blank()

    if n_updated:
        backup_path = prmag_path + ".bak"
        if not os.path.exists(backup_path):
            with open(backup_path, "w", encoding="utf-8") as f:
                f.write(original_text)
        with open(prmag_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    return n_updated, unmatched, still_missing


def complete_site_info(
    prmag_path: str, table_path: str, encoding: str = "utf-8",
) -> Tuple[int, List[str], List[Tuple[str, List[str]]]]:
    """Table indexee par SITE (colonne 'site' obligatoire) - voir
    docstring module. Retourne (nb_specimens_mis_a_jour,
    liste_specimens_sans_site_correspondant_dans_la_table,
    liste_specimens_encore_incomplets).

    Une SECONDE colonne 'site' (meme nom d'en-tete repete, ex.
    "site\\tsite\\tLat\\tLon\\t..." comme Complement_tect.txt : ancien
    nom de campagne -> nouveau nom consolide) est acceptee pour RENOMMER
    le site - demande explicite utilisateur ("i cannot select site and
    replace it by an other site"). La PREMIERE colonne 'site' reste
    toujours la cle de recherche (nom deja present dans le .prmag) ;
    construire `record` par position (pas par dict(zip(header, row)),
    qui ne garderait que la DERNIERE colonne 'site' et perdrait la cle
    de recherche pour les tables a une seule colonne 'site' - non, en
    fait l'inverse : cela ferait passer la valeur de renommage comme cle
    de recherche par erreur) fait naturellement atterrir la seconde
    colonne 'site', si presente, dans record['site'] (lu par
    _apply_table comme cible de renommage) sans rien changer pour les
    tables existantes a une seule colonne 'site' (record ne contient
    alors aucune cle 'site', comportement inchange)."""
    header, rows = _read_table_rows(table_path, encoding=encoding)
    if "site" not in header:
        raise ValueError("Table file must have a 'site' column header (site-level mode).")
    lookup_idx = header.index("site")
    table: Dict[str, dict] = {}
    for row in rows:
        if lookup_idx >= len(row):
            continue
        site = row[lookup_idx].strip()
        if not site:
            continue
        record = {
            col: row[i].strip()
            for i, col in enumerate(header)
            if i != lookup_idx and i < len(row)
        }
        table[site] = record
    return _apply_table(prmag_path, table, key_field="site")


def complete_specimen_info(
    prmag_path: str, table_path: str, encoding: str = "utf-8",
) -> Tuple[int, List[str], List[Tuple[str, List[str]]]]:
    """Table indexee par SPECIMEN complet (colonne 'specimen' obligatoire) -
    voir docstring module. Retourne (nb_specimens_mis_a_jour,
    liste_specimens_sans_correspondance_dans_la_table,
    liste_specimens_encore_incomplets)."""
    header, rows = _read_table_rows(table_path, encoding=encoding)
    if "specimen" not in header:
        raise ValueError("Table file must have a 'specimen' column header (specimen-level mode).")
    table: Dict[str, dict] = {}
    for row in rows:
        record = dict(zip(header, row))
        specimen = record.get("specimen", "").strip()
        if specimen:
            table[specimen] = record
    return _apply_table(prmag_path, table, key_field="specimen")


def complete_sample_height(
    prmag_path: str, table_path: str, encoding: str = "utf-8",
) -> Tuple[int, List[str], List[Tuple[str, List[str]]]]:
    """Table indexee par SAMPLE (colonnes 'sample' ET 'stratigraphic_height'
    obligatoires) - cas particulier explicite utilisateur ("add a
    specific case for stratigraphic_height filled from sample and height
    as all specimens have the same height") : une ligne par SAMPLE
    physique (pas par specimen ni par site) - la valeur s'applique a
    TOUS les specimens de ce fichier partageant ce `sample` (typiquement
    plusieurs lettres A/B/C decoupees du meme echantillon), evitant de
    repeter la meme hauteur autant de fois qu'il y a de specimens. Ne
    touche QUE `stratigraphic_height` - pour lat/lon/geologie, voir
    complete_site_info/complete_specimen_info. Retourne (nb_specimens_
    mis_a_jour, liste_specimens_sans_sample_correspondant_dans_la_table,
    liste_specimens_encore_incomplets - couvre TOUS les champs, pas
    seulement stratigraphic_height, pour une vue d'ensemble utile meme
    dans ce mode restreint)."""
    header, rows = _read_table_rows(table_path, encoding=encoding)
    if "sample" not in header:
        raise ValueError("Table file must have a 'sample' column header (sample-level mode).")
    if "stratigraphic_height" not in header:
        raise ValueError("Table file must have a 'stratigraphic_height' column header (sample-level mode).")
    table: Dict[str, dict] = {}
    for row in rows:
        record = dict(zip(header, row))
        sample = record.get("sample", "").strip()
        if sample:
            table[sample] = record
    return _apply_table(prmag_path, table, key_field="sample")
