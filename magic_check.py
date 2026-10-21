"""
Controle detaille d'une contribution MagIC (fichier magic_contribution_*.txt)
- demande explicite utilisateur ("lors de l'importation des contributions
Magic, est-ce possible de voir si il y a des doubles de sites, samples et
specimens. Est-ce possible de faire un rapport detaille sur les problemes
rencontres lors de l'importation des contributions. Rechercher si possible
tous les problemes possibles et imaginables").

Deux familles de controles :

1) Regles du modele de donnees MagIC 3.0 lui-meme (data_model.json livre
   avec PmagPy), appliquees colonne par colonne : tables et colonnes
   inconnues, colonnes obligatoires (required, requiredUnless, requiredIf,
   requiredUnlessSynthetic), types (Number/Integer), bornes min/max (valeurs
   fixes ou autre colonne), vocabulaires controles (cv), codes de methode,
   references entre tables (in("samples.sample")...), unicite.

2) Problemes reels rencontres dans des contributions publiees, que le
   modele ne voit pas (voir l'audit des 40 contributions Lu et al. 2026) :
   doublons de noms (lignes identiques / metadonnees contradictoires / noms
   presque identiques), noms suspects (guillemets, "0.0"), hierarchie
   orpheline, paliers mesures plusieurs fois, lignes sans code de
   traitement LT- (illisibles par Demag GUI), directions de measurements
   hors repere specimen, incoherences magn_x/y/z <-> dir/moment, unites
   suspectes (cm3 au lieu de m3, g au lieu de kg, degC au lieu de K, mT au
   lieu de T), dip positifs, coordonnees nulles ou grossieres, classe
   geologique extraterrestre pour des echantillons orientes, etc.

Le controle ne modifie rien : il produit un rapport texte trie par gravite
(ERROR / WARNING / INFO), chaque probleme avec son nombre d'occurrences et
quelques exemples (table, ligne du fichier, nom).
"""

import collections
import difflib
import json
import math
import os
import re
from typing import Dict, List, Optional, Tuple

_SEVERITIES = ("ERROR", "WARNING", "INFO")
# Categories placed first in the report - user request ("start with error
# on Orientation"): an orientation error makes every in situ or
# tilt-corrected direction of the specimen wrong.
_ORIENTATION_CATEGORIES = ("Orientation", "Directions")
_ORIENTATION_COLUMNS = ("azimuth", "dip", "bed_dip", "bed_dip_direction")
_MAX_EXAMPLES = 6
_NAME_TABLES = (("locations", "location"), ("sites", "site"), ("samples", "sample"), ("specimens", "specimen"))
_PARENT = {"sites": "location", "samples": "site", "specimens": "sample", "measurements": "specimen"}
_DEMAG_GUI_LT = ("LT-NO", "LT-AF-Z", "LT-T-Z", "LT-M-Z", "LT-LT-Z")

# Colonnes qui decrivent l'objet lui-meme (et non un resultat) : deux lignes
# du meme nom qui different sur l'une d'elles se contredisent.
_IDENTITY_COLUMNS = {
    "locations": ("location_type", "lat_s", "lat_n", "lon_w", "lon_e", "geologic_classes", "lithologies"),
    "sites": ("location", "lat", "lon", "geologic_classes", "geologic_types", "lithologies", "age", "age_low",
              "age_high", "age_unit", "height"),
    "samples": ("site", "lat", "lon", "azimuth", "dip", "bed_dip", "bed_dip_direction", "geologic_classes",
                "geologic_types", "lithologies", "height"),
    "specimens": ("sample", "volume", "weight", "geologic_classes", "lithologies"),
}


# ----------------------------------------------------------------------
# Modele de donnees
# ----------------------------------------------------------------------

def _data_model_dir() -> Optional[str]:
    try:
        import pmagpy
    except ImportError:
        return None
    path = os.path.join(os.path.dirname(pmagpy.__file__), "data_model")
    return path if os.path.isdir(path) else None


DATA_MODEL_URL = "https://magic.earthref.org/MagIC/data-models/3.0.json"
_REFRESH_DAYS = 30


def _cache_path() -> str:
    return os.path.join(os.path.expanduser("~"), ".STARpaleomag_Py", "magic_data_model_3.0.json")


def _bundled_path() -> str:
    import sys
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, "magic_model", "data_model_3.0.json")


def _official_data_model() -> Tuple[Optional[dict], str]:
    """Modele officiel MagIC 3.0 tel que publie par EarthRef - demande
    explicite utilisateur ("est-ce possible de verifier aupres de
    https://magic.earthref.org/MagIC/data-models/3.0") : la copie livree
    avec PmagPy date de 2019 alors que le modele officiel continue d'evoluer
    (colonnes ajoutees, bornes elargies - specimens.dip passe a -180/180,
    samples.bed_dip n'a plus de borne, measurements.sequence n'est plus
    obligatoire...). Ordre : copie locale rafraichie depuis EarthRef au plus
    tous les 30 jours (delai de 5 s, sans bloquer si hors ligne), sinon la
    copie livree avec STARpaleomag_Py (magic_model/). Retourne (modele,
    provenance)."""
    import time
    import urllib.request
    cache = _cache_path()
    fresh = os.path.exists(cache) and time.time() - os.path.getmtime(cache) < _REFRESH_DAYS * 86400
    if not fresh:
        try:
            import ssl
            try:
                # Python de python.org sous macOS : pas de certificats systeme
                # tant que "Install Certificates" n'a pas ete lance - certifi
                # les fournit (meme solution que StereoUtils_Py).
                import certifi
                context = ssl.create_default_context(cafile=certifi.where())
            except ImportError:
                context = ssl.create_default_context()
            with urllib.request.urlopen(DATA_MODEL_URL, timeout=5, context=context) as response:
                text = response.read().decode("utf-8-sig")
            if "tables" in json.loads(text):
                os.makedirs(os.path.dirname(cache), exist_ok=True)
                with open(cache, "w", encoding="utf-8") as f:
                    f.write(text)
        except Exception:
            pass
    for path, origin in ((cache, f"EarthRef ({DATA_MODEL_URL})"), (_bundled_path(), "copy shipped with "
                                                                    "STARpaleomag_Py")):
        try:
            with open(path, encoding="utf-8-sig") as f:
                return json.load(f), origin
        except (OSError, ValueError):
            continue
    return None, ""


class DataModel:
    def __init__(self):
        folder = _data_model_dir()
        if folder is None:
            raise RuntimeError("PmagPy's MagIC data model files were not found")

        def load(name):
            with open(os.path.join(folder, name), encoding="utf-8-sig") as f:
                return json.load(f)

        model, self.origin = _official_data_model()
        if model is None:
            model, self.origin = load("data_model.json"), "copy shipped with PmagPy"
        self.version = model.get("data_model_version", "?")
        self.updated = model.get("updated_day", "?")
        self.tables = model["tables"]
        # Codes de methode officiels (magic_model/method_codes.json, telecharge
        # par l'utilisateur depuis magic.earthref.org/MagIC/method-codes le
        # 2026-10-21, 673 codes) plutot que la liste de PmagPy (2019, 630) :
        # 47 codes ajoutes depuis (LP-T, LP-AF, LP-CA-*, LT-AF-Z-XYZ...), 4
        # retires (LP-PI-CA-*, remplaces par LP-CA-*).
        official_codes = os.path.join(os.path.dirname(_bundled_path()), "method_codes.json")
        try:
            with open(official_codes, encoding="utf-8-sig") as f:
                code_tables = json.load(f)
            self.codes_origin = "official list (October 2026)"
        except (OSError, ValueError):
            code_tables = load("method_codes.json")
            self.codes_origin = "PmagPy list (2019)"
        self.codes = {c["code"] for cat in code_tables.values() for c in cat["codes"]}
        self.codes_upper = {c.upper(): c for c in self.codes}
        vocab_files = sorted(f for f in os.listdir(folder) if f.startswith("controlled_vocabularies_"))
        vocab = load(vocab_files[-1]) if vocab_files else {}
        self.vocab = {
            name: {str(item["item"]) for item in entry.get("items", []) if "item" in item}
            for name, entry in vocab.items() if isinstance(entry, dict)
        }
        # Listes officielles plus recentes livrees avec STARpaleomag_Py
        # (magic_model/vocabulary_<nom>.txt, une valeur par ligne) : elles
        # REMPLACENT celles de PmagPy (2019) - EarthRef ne publie plus ces
        # listes en fichier ; "lithology" (1003 termes, oct. 2026, fournie
        # par l'utilisateur) contient par exemple "Granitoid", absent de la
        # liste de PmagPy et signale a tort auparavant.
        self.updated_vocabularies = []
        folder_model = os.path.dirname(_bundled_path())
        try:
            names = sorted(os.listdir(folder_model))
        except OSError:
            names = []
        for fname in names:
            m = re.match(r"vocabulary_(.+)\.txt$", fname)
            if not m:
                continue
            with open(os.path.join(folder_model, fname), encoding="utf-8") as f:
                terms = {line.strip() for line in f if line.strip() and not line.startswith("#")}
            if terms:
                self.vocab[m.group(1)] = terms
                self.updated_vocabularies.append(m.group(1))
        self.vocab_upper = {name: {v.upper(): v for v in values} for name, values in self.vocab.items()}

    def columns(self, table: str) -> Dict[str, dict]:
        return self.tables.get(table, {}).get("columns", {})


_MODEL: Optional[DataModel] = None


def data_model() -> DataModel:
    global _MODEL
    if _MODEL is None:
        _MODEL = DataModel()
    return _MODEL


# ----------------------------------------------------------------------
# Lecture brute (garde les doublons de tables, les numeros de ligne, la
# longueur reelle de chaque ligne - ce que parse_magic_contribution efface)
# ----------------------------------------------------------------------

class Table:
    def __init__(self, name: str, header_line: int, columns: List[str]):
        self.name = name
        self.header_line = header_line
        self.columns = columns
        self.rows: List[Dict[str, str]] = []
        self.lines: List[int] = []
        self.bad_length: List[Tuple[int, int]] = []


def read_tables(path: str) -> List[Table]:
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        lines = f.read().replace("\r\n", "\n").replace("\r", "\n").split("\n")
    tables: List[Table] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip() or line.startswith(">>>>"):
            i += 1
            continue
        parts = line.split("\t")
        name = parts[1].strip() if len(parts) > 1 else parts[0].strip()
        if i + 1 >= len(lines):
            break
        columns = [c.strip() for c in lines[i + 1].split("\t")]
        table = Table(name, i + 2, columns)
        i += 2
        while i < len(lines) and not lines[i].startswith(">>>>"):
            if lines[i].strip():
                values = lines[i].split("\t")
                if len(values) != len(columns):
                    table.bad_length.append((i + 1, len(values)))
                table.rows.append({c: v.strip() for c, v in zip(columns, values)})
                table.lines.append(i + 1)
            i += 1
        tables.append(table)
    return tables


# ----------------------------------------------------------------------
# Rapport
# ----------------------------------------------------------------------

class Report:
    def __init__(self, path: str):
        self.path = path
        self.issues: "collections.OrderedDict[Tuple[str, str, str], dict]" = collections.OrderedDict()
        self.summary: List[str] = []
        self.positive_dip: Tuple[int, int] = (0, 0)  # (samples with dip > 0, samples with a dip)

    def add(self, severity: str, category: str, message: str, example: Optional[str] = None, count: int = 1,
            aggregate: bool = False):
        """`aggregate` : l'exemple est deja la liste (tronquee) de tous les
        cas - pas de ligne "... (n more)" dans le rapport."""
        key = (severity, category, message)
        issue = self.issues.setdefault(key, {"count": 0, "examples": [], "aggregate": aggregate})
        issue["count"] += count
        if example is not None and len(issue["examples"]) < _MAX_EXAMPLES:
            issue["examples"].append(example)

    def counts(self) -> Dict[str, int]:
        out = {s: 0 for s in _SEVERITIES}
        for (severity, _c, _m) in self.issues:
            out[severity] += 1
        return out

    def format(self, error_mark: str = "") -> str:
        """Texte du rapport. Les problemes d'orientation (categories
        _ORIENTATION_CATEGORIES) viennent en tete, toutes gravites
        confondues (erreurs d'abord), puis le reste par gravite. Avec
        `error_mark` (marqueur de la console de l'appli), chaque ligne
        d'erreur en est entouree pour etre affichee en rouge et en gras -
        demande explicite utilisateur ("is it possible to write it in red
        and bold") ; sans marqueur (fichier texte), le texte est le meme."""
        counts = self.counts()
        n_orient = sum(1 for k in self.issues if k[1] in _ORIENTATION_CATEGORIES and k[0] == "ERROR")
        out = [
            f"MagIC contribution check - {os.path.basename(self.path)}",
            f"{counts['ERROR']} error(s) (orientation: {n_orient}), {counts['WARNING']} warning(s), "
            f"{counts['INFO']} note(s) (each line below is one kind of problem, with its number of occurrences).",
            "",
        ]
        out.extend(self.summary)

        def emit(items, show_severity):
            current = None
            for (sev, category, message), issue in items:
                if category != current:
                    out.append(f"--- {category}")
                    current = category
                mark = error_mark if sev == "ERROR" else ""
                label = f"{sev}: " if show_severity else ""
                out.append(f"{mark}  [{issue['count']}] {label}{message}{mark}")
                for example in issue["examples"]:
                    out.append(f"        e.g. {example}")
                if issue["count"] > len(issue["examples"]) and issue["examples"] and not issue["aggregate"]:
                    out.append(f"        ... ({issue['count'] - len(issue['examples'])} more)")

        orientation = [(k, v) for k, v in self.issues.items() if k[1] in _ORIENTATION_CATEGORIES]
        if orientation:
            out.append("")
            title = "===== ORIENTATION ====="
            out.append(f"{error_mark}{title}{error_mark}" if n_orient else title)
            emit(sorted(orientation, key=lambda kv: (_SEVERITIES.index(kv[0][0]), kv[0][1])), True)
        for severity in _SEVERITIES:
            items = [(k, v) for k, v in self.issues.items()
                     if k[0] == severity and k[1] not in _ORIENTATION_CATEGORIES]
            if not items:
                continue
            out.append("")
            title = f"===== {severity}S ====="
            out.append(f"{error_mark}{title}{error_mark}" if severity == "ERROR" else title)
            emit(sorted(items, key=lambda kv: kv[0][1]), False)
        if not self.issues:
            out.append("No problem found.")
        return "\n".join(out) + "\n"


def _where(table: Table, k: int, name_col: Optional[str] = None) -> str:
    row = table.rows[k]
    name = row.get(name_col or _name_column(table.name) or "", "")
    return f"{table.name} line {table.lines[k]}" + (f" ({name})" if name else "")


def _name_column(table: str) -> Optional[str]:
    return dict(_NAME_TABLES).get(table) or ("measurement" if table == "measurements" else None)


def _num(value: str) -> Optional[float]:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _tokens(value: str) -> List[str]:
    return [t.strip() for t in value.split(":") if t.strip()]


def _normalized_name(name: str) -> str:
    return re.sub(r"\s+", "", name.replace('"', "").replace("”", "").replace("“", "")
                  .replace("'", "")).upper()


# ----------------------------------------------------------------------
# 1) Regles du modele MagIC
# ----------------------------------------------------------------------

def _check_structure(tables: List[Table], model: DataModel, report: Report):
    seen = collections.Counter(t.name for t in tables)
    for name, n in seen.items():
        if name not in model.tables:
            report.add("ERROR", "Structure", f"unknown table '{name}' (not in the MagIC data model)")
        if n > 1:
            report.add("WARNING", "Structure", f"table '{name}' appears {n} times in the file")
    for table in tables:
        cols = model.columns(table.name)
        for col, n in collections.Counter(table.columns).items():
            if n > 1:
                report.add("ERROR", "Structure", f"column {table.name}.{col} appears {n} times in the header")
        if not cols:
            continue
        for col in table.columns:
            if col and col not in cols:
                filled = sum(1 for r in table.rows if r.get(col))
                elsewhere = [t for t in model.tables if col in model.columns(t)]
                hint = f" - exists in table(s) {', '.join(elsewhere)}" if elsewhere else ""
                report.add("WARNING", "Unknown column",
                           f"{table.name}.{col} is not a MagIC 3.0 column (filled in {filled}/{len(table.rows)} "
                           f"rows){hint}")
        for line, n in table.bad_length:
            report.add("ERROR", "Structure",
                       f"{table.name}: row with a number of fields different from the header "
                       f"({len(table.columns)} columns)", f"line {line}: {n} fields")
        if table.rows and table.columns and all(not c for c in table.columns):
            report.add("ERROR", "Structure", f"{table.name}: empty header line")


def _check_columns(tables: List[Table], model: DataModel, report: Report, index: Dict[str, set]):
    by_name = {t.name: t for t in tables}
    sample_rows = {r.get("sample", ""): r for r in by_name["samples"].rows} if "samples" in by_name else {}
    site_rows = {r.get("site", ""): r for r in by_name["sites"].rows} if "sites" in by_name else {}

    def inherited_value(row, col):
        """Valeur donnee a un niveau superieur : pour un specimen, son
        sample puis le site de ce sample ; pour un sample, son site - la
        pratique MagIC (remarque utilisateur : "when the lat lon and
        geologic information are given at the site level, they usually
        apply to the sample and specimen level") ; l'import fait d'ailleurs
        de meme (convert_magic_to_r._sample_header_block)."""
        sample = sample_rows.get(row.get("sample", ""), {})
        if sample.get(col):
            return sample.get(col)
        site = site_rows.get(row.get("site") or sample.get("site", ""), {})
        return site.get(col, "")

    synthetic_samples = set()
    for table in tables:
        if table.name == "samples":
            for r in table.rows:
                if "synthetic" in (r.get("material_type", "") + r.get("method_codes", "")).lower():
                    synthetic_samples.add(r.get("sample", ""))

    for table in tables:
        cols = model.columns(table.name)
        present = set(table.columns)
        for col, spec in cols.items():
            validations = spec.get("validations", [])
            ctype = spec.get("type")
            for v in validations:
                if v == "required()":
                    if col not in present:
                        if table.rows:
                            report.add("ERROR", "Required column",
                                       f"{table.name}.{col} is required but the column is absent")
                    else:
                        empty = [k for k, r in enumerate(table.rows) if not r.get(col)]
                        for k in empty:
                            report.add("ERROR", "Required column", f"{table.name}.{col} is required but empty",
                                       _where(table, k))
                elif v == "requiredUnlessSynthetic()" and table.rows:
                    for k, r in enumerate(table.rows):
                        if not r.get(col) and r.get("sample", "") not in synthetic_samples \
                                and not inherited_value(r, col):
                            report.add("WARNING", "Required column",
                                       f"{table.name}.{col} is empty at every level (specimen, sample and "
                                       "site)" if table.name == "specimens" else
                                       f"{table.name}.{col} is empty at both the sample and site level",
                                       _where(table, k))
                elif v.startswith("requiredUnless(") and table.rows:
                    others = re.findall(r'"([^"]+)"', v)
                    for k, r in enumerate(table.rows):
                        if not r.get(col) and not any(r.get(o) for o in others):
                            report.add("WARNING", "Required column",
                                       f"{table.name}: neither {col} nor {'/'.join(others)} is given",
                                       _where(table, k))
                elif v.startswith("requiredIf(") and table.rows:
                    other = re.findall(r'"([^"]+)"', v)[0]
                    for k, r in enumerate(table.rows):
                        if r.get(other) and not r.get(col):
                            report.add("ERROR", "Required column",
                                       f"{table.name}.{other} is given without {col}", _where(table, k))
            if col not in present:
                continue
            for k, r in enumerate(table.rows):
                value = r.get(col, "")
                if not value:
                    continue
                if ctype in ("Number", "Integer"):
                    x = _num(value)
                    if x is None:
                        report.add("ERROR", "Value type", f"{table.name}.{col} is not a number",
                                   f"{_where(table, k)}: '{value}'")
                        continue
                    if ctype == "Integer" and x != int(x):
                        report.add("WARNING", "Value type", f"{table.name}.{col} should be an integer",
                                   f"{_where(table, k)}: {value}")
                    for v in validations:
                        m = re.match(r'(min|max)\((.+)\)$', v)
                        if not m:
                            continue
                        bound_s = m.group(2).strip()
                        if bound_s.startswith('"'):
                            bound = _num(r.get(bound_s.strip('"'), ""))
                            label = bound_s.strip('"')
                        else:
                            bound, label = _num(bound_s), bound_s
                        if bound is None:
                            continue
                        if (m.group(1) == "min" and x < bound) or (m.group(1) == "max" and x > bound):
                            hint = ""
                            if col in ("lon", "lon_w", "lon_e") and x < 0:
                                hint = " - negative longitude, MagIC expects 0-360"
                            report.add("ERROR" if not hint else "WARNING",
                                       "Orientation" if col in _ORIENTATION_COLUMNS else "Out of range",
                                       f"{table.name}.{col} {'below' if m.group(1) == 'min' else 'above'} "
                                       f"{m.group(1)} {label}{hint}", f"{_where(table, k)}: {value}")
                for v in validations:
                    if v.startswith("cv("):
                        vocab_name = v[4:-2]
                        allowed = model.vocab.get(vocab_name)
                        if not allowed:
                            continue
                        for token in (_tokens(value) if ctype == "List" else [value]):
                            if token in allowed:
                                continue
                            standard = model.vocab_upper[vocab_name].get(token.upper())
                            if standard:
                                report.add("INFO", "Controlled vocabulary",
                                           f"{table.name}.{col}: '{token}' written with a different case "
                                           f"from the standard '{standard}'", _where(table, k))
                            else:
                                close = difflib.get_close_matches(token, sorted(allowed), n=3, cutoff=0.7)
                                hint = f" (closest standard terms: {', '.join(close)})" if close else ""
                                report.add("WARNING", "Controlled vocabulary",
                                           f"{table.name}.{col}: '{token}' is not in the MagIC vocabulary "
                                           f"'{vocab_name}'{hint}", _where(table, k))
                    elif v == 'type("method_codes")':
                        for code in _tokens(value):
                            if code in model.codes:
                                continue
                            standard = model.codes_upper.get(code.upper())
                            if standard:
                                report.add("WARNING", "Method codes",
                                           f"{table.name}: method code '{code}' should be written '{standard}'",
                                           _where(table, k))
                            else:
                                report.add("WARNING", "Method codes",
                                           f"{table.name}: '{code}' is not a MagIC method code",
                                           _where(table, k))
                    elif v.startswith("in("):
                        target = v[4:-2]
                        if target not in index:
                            continue
                        for token in (_tokens(value) if ctype == "List" else [value]):
                            if token not in index[target]:
                                report.add("ERROR", "Broken link between tables",
                                           f"{table.name}.{col} refers to a {target.split('.')[1]} absent "
                                           f"from the {target.split('.')[0]} table",
                                           f"{_where(table, k)}: '{token}'")
                    elif v == "unique()":
                        pass  # voir _check_measurements (colonne sequence)


# ----------------------------------------------------------------------
# 2) Problemes reels hors modele
# ----------------------------------------------------------------------

def _check_duplicates(tables: Dict[str, Table], report: Report):
    for table_name, name_col in _NAME_TABLES:
        table = tables.get(table_name)
        if table is None or name_col not in table.columns:
            continue
        by_name: Dict[str, List[int]] = collections.defaultdict(list)
        for k, r in enumerate(table.rows):
            by_name[r.get(name_col, "")].append(k)
        n_names = len(by_name)
        for name, ks in by_name.items():
            if len(ks) < 2 or not name:
                continue
            rows = [table.rows[k] for k in ks]
            by_signature: Dict[tuple, List[int]] = collections.defaultdict(list)
            for k in ks:
                by_signature[tuple(sorted(table.rows[k].items()))].append(k)
            signatures = by_signature
            for same in by_signature.values():
                if len(same) > 1:
                    report.add("WARNING", "Duplicates",
                               f"{table_name}: the same {name_col} appears on identical rows (pure duplicates)",
                               f"{name_col} '{name}': {len(same)} identical rows (lines "
                               f"{', '.join(str(table.lines[k]) for k in same[:6])})")
            conflicts, orient_conflicts = [], []
            for col in _IDENTITY_COLUMNS.get(table_name, ()):
                values = {r.get(col, "") for r in rows if r.get(col, "")}
                if len(values) > 1:
                    (orient_conflicts if col in _ORIENTATION_COLUMNS else conflicts).append(
                        f"{col}={'/'.join(sorted(values)[:4])}")
            if orient_conflicts:
                report.add("ERROR", "Orientation",
                           f"{table_name}: the same {name_col} is given contradictory orientations "
                           "(azimuth, dip or bedding)",
                           f"{name_col} '{name}' ({len(ks)} rows): {'; '.join(orient_conflicts)}")
            if conflicts:
                report.add("ERROR", "Duplicates",
                           f"{table_name}: the same {name_col} is described with contradictory values",
                           f"{name_col} '{name}' ({len(ks)} rows): {'; '.join(conflicts)}")
            elif len(signatures) > 1:
                # Plusieurs lignes de resultats d'un meme objet : normal
                # quand elles sont dans des reperes differents (specimen /
                # in situ / tectonique, dir_tilt_correction -1/0/100) -
                # remarque utilisateur ("lorsqu'il y a des resultats en
                # coordonnee specimen et IS, c'est normal d'avoir 2 lignes").
                # Signale seulement deux resultats du MEME repere pour la MEME
                # composante (et le meme type de resultat).
                def result_key(r):
                    kind = "dir" if r.get("dir_dec") else ("int" if r.get("int_abs") else "other")
                    return (kind, r.get("dir_tilt_correction", ""),
                            r.get("dir_comp") or r.get("dir_comp_name") or "")
                keys = collections.Counter(result_key(table.rows[k]) for k in ks)
                same = [key for key, n in keys.items() if n > 1 and key[0] != "other"]
                if same:
                    kind, tilt, comp = same[0]
                    report.add("WARNING", "Duplicates",
                               f"{table_name}: several {'direction' if kind == 'dir' else 'intensity'} results "
                               f"for the same {name_col} in the same coordinate system and component",
                               f"{name_col} '{name}': {keys[same[0]]} rows with dir_tilt_correction="
                               f"{tilt or 'n.d'}, component {comp or 'n.d'}")
        report.summary.append(f"  {table_name}: {len(table.rows)} row(s), {n_names} distinct {name_col}(s)")

        # Noms presque identiques : meme nom une fois retires guillemets,
        # espaces et casse.
        groups: Dict[str, set] = collections.defaultdict(set)
        for name in by_name:
            if name:
                groups[_normalized_name(name)].add(name)
        for variants in groups.values():
            if len(variants) > 1:
                report.add("WARNING", "Duplicates",
                           f"{table_name}: {name_col} names identical except for quotes, spaces or case "
                           "(probably the same object entered twice)",
                           " / ".join(f"'{v}'" for v in sorted(variants)))
        for name in by_name:
            if not name:
                continue
            if any(q in name for q in ('"', "“", "”", "'")):
                report.add("WARNING", "Suspicious names", f"{table_name}: {name_col} name containing quotes",
                           f"'{name}'")
            elif name != name.strip() or "  " in name:
                report.add("INFO", "Suspicious names",
                           f"{table_name}: {name_col} name with leading/trailing or double spaces", f"'{name}'")
            if re.fullmatch(r"[0.\-]+", name):
                report.add("WARNING", "Suspicious names",
                           f"{table_name}: {name_col} name that looks like a placeholder or a number", f"'{name}'")
            elif name.endswith(("-", ".", "_")):
                report.add("INFO", "Suspicious names", f"{table_name}: {name_col} name ending with a separator",
                           f"'{name}'")


def _check_hierarchy(tables: Dict[str, Table], report: Report):
    children = {"locations": ("sites", "location"), "sites": ("samples", "site"),
                "samples": ("specimens", "sample"), "specimens": ("measurements", "specimen")}
    for parent, (child, col) in children.items():
        pt, ct = tables.get(parent), tables.get(child)
        if pt is None or ct is None:
            continue
        used = {r.get(col, "") for r in ct.rows}
        name_col = _name_column(parent)
        orphans = sorted({r.get(name_col, "") for r in pt.rows} - used - {""})
        if orphans:
            report.add("INFO", "Hierarchy", f"{parent} with no {child} attached ({len(orphans)})",
                       ", ".join(orphans[:12]) + (" ..." if len(orphans) > 12 else ""), count=len(orphans), aggregate=True)
    for child, col in (("samples", "site"), ("specimens", "sample")):
        ct = tables.get(child)
        if ct is None:
            continue
        name_col = _name_column(child)
        parents: Dict[str, set] = collections.defaultdict(set)
        for r in ct.rows:
            if r.get(col):
                parents[r.get(name_col, "")].add(r.get(col))
        for name, ps in parents.items():
            if len(ps) > 1:
                report.add("ERROR", "Hierarchy", f"{child}: the same {name_col} belongs to several {col}s",
                           f"'{name}' -> {', '.join(sorted(ps))}")
    # Un site par echantillon et sites partageant la meme coordonnee
    sites, samples = tables.get("sites"), tables.get("samples")
    if sites is not None and samples is not None:
        n_sites = len({r.get("site") for r in sites.rows})
        n_samples = len({r.get("sample") for r in samples.rows})
        coords = collections.Counter((r.get("lat"), r.get("lon")) for r in
                                     {r.get("site"): r for r in sites.rows}.values() if r.get("lat") and r.get("lon"))
        shared = sum(1 for n in coords.values() if n > 1)
        if n_sites and n_sites >= 0.9 * n_samples and shared:
            report.add("INFO", "Hierarchy",
                       f"about one site per sample ({n_sites} sites for {n_samples} samples) and {shared} "
                       "coordinate(s) shared by several sites: the real sites may only survive as coordinates "
                       "(see 'Group the sites with the same lat/lon' at import)")


def _check_samples_sites(tables: Dict[str, Table], report: Report):
    samples, sites, specimens = tables.get("samples"), tables.get("sites"), tables.get("specimens")
    site_rows = {r.get("site"): r for r in sites.rows} if sites else {}
    if samples is not None:
        # dip au niveau sample, ou a defaut au niveau specimen (MagIC permet
        # d'orienter le specimen lui-meme ; vu dans une vraie contribution,
        # 20476, ou samples.dip est vide et specimens.dip renseigne).
        dips = [(_num(r.get("dip", "")), samples, k) for k, r in enumerate(samples.rows)]
        dips = [(d, t, k) for d, t, k in dips if d is not None]
        if specimens is not None and "dip" in specimens.columns:
            with_dip = {samples.rows[k].get("sample", "") for _d, _t, k in dips}
            dips += [(d, specimens, k) for k, r in enumerate(specimens.rows)
                     for d in [_num(r.get("dip", ""))]
                     if d is not None and r.get("sample", "") not in with_dip]
        positive = [(t, k) for d, t, k in dips if d > 0]
        report.positive_dip = (len(positive), len(dips))
        if dips and len(positive) > len(dips) / 2:
            # ERROR (plus WARNING) - remarque utilisateur ("the incorrect
            # positive dip") : un dip de convention erronee fausse toutes les
            # directions in situ et tectoniques si on l'importe tel quel.
            report.add("ERROR", "Orientation",
                       f"dip positive for {len(positive)}/{len(dips)} samples/specimens - MagIC dip (plunge of "
                       "the X axis) is normally negative or zero; typical of the Utrecht convention "
                       "(dip = 90 - plunge) or of dip entered with the wrong sign",
                       _where(*positive[0]))
        for k, r in enumerate(samples.rows):
            az, dip = r.get("azimuth", ""), r.get("dip", "")
            if not az or not dip:
                report.add("WARNING", "Orientation", "samples: no azimuth/dip (unoriented sample: directions "
                           "only usable in specimen coordinates)", _where(samples, k))
            if r.get("bed_dip_direction") and not r.get("bed_dip"):
                report.add("WARNING", "Orientation", "samples: bed_dip_direction without bed_dip",
                           _where(samples, k))
            if r.get("bed_dip") and not r.get("bed_dip_direction"):
                report.add("WARNING", "Orientation", "samples: bed_dip without bed_dip_direction",
                           _where(samples, k))
            site = site_rows.get(r.get("site"))
            lat, lon = _num(r.get("lat", "")), _num(r.get("lon", ""))
            if site and lat is not None and lon is not None:
                slat, slon = _num(site.get("lat", "")), _num(site.get("lon", ""))
                if slat is not None and slon is not None and (
                        abs(lat - slat) > 0.05 or abs(((lon - slon) + 180) % 360 - 180) > 0.05):
                    report.add("WARNING", "Coordinates",
                               "samples: lat/lon more than 0.05 deg away from those of their site",
                               f"{_where(samples, k)}: {lat}/{lon} vs site {slat}/{slon}")
            classes = (r.get("geologic_classes", "") or (site or {}).get("geologic_classes", "")).lower()
            if any(c in classes for c in ("extraterrestrial", "meteorite", "lunar", "martian")) and az and dip:
                report.add("WARNING", "Geology",
                           "samples classed extraterrestrial/meteorite/lunar/martian but oriented in the field "
                           "(azimuth/dip given) - check geologic_classes (often a slip for 'Extrusive')",
                           f"{_where(samples, k)}: {classes}")
    for table in (sites, samples):
        if table is None:
            continue
        for k, r in enumerate(table.rows):
            lat, lon = _num(r.get("lat", "")), _num(r.get("lon", ""))
            if lat == 0 and lon == 0:
                report.add("WARNING", "Coordinates", f"{table.name}: lat = lon = 0 (missing coordinates?)",
                           _where(table, k))
            for col in ("lat", "lon"):
                v = r.get(col, "")
                if v and _num(v) is not None and "." in v and len(v.split(".")[1]) <= 1 and table is sites:
                    report.add("INFO", "Coordinates",
                               "sites: coordinate given with one decimal or less (~10 km)",
                               f"{_where(table, k)}: {col}={v}")
                    break
            if table is sites and r.get("dir_dec") and not r.get("dir_tilt_correction"):
                report.add("WARNING", "Results", "sites: direction without dir_tilt_correction (coordinate "
                           "system unknown)", _where(table, k))
            if table is sites and r.get("vgp_lat") and not r.get("dir_dec"):
                report.add("INFO", "Results", "sites: VGP given without the site direction", _where(table, k))
            lo, hi = _num(r.get("age_low", "")), _num(r.get("age_high", ""))
            if lo is not None and hi is not None and lo > hi:
                report.add("ERROR", "Age", f"{table.name}: age_low greater than age_high",
                           f"{_where(table, k)}: {lo} > {hi}")
            if (r.get("age") or r.get("age_low") or r.get("age_high")) and not r.get("age_unit") \
                    and "age_unit" in table.columns:
                report.add("WARNING", "Age", f"{table.name}: age given without age_unit", _where(table, k))
    if specimens is not None:
        for k, r in enumerate(specimens.rows):
            vol, wt = _num(r.get("volume", "")), _num(r.get("weight", ""))
            if vol is not None and vol > 1e-3:
                report.add("WARNING", "Units", "specimens.volume larger than 1 litre - probably cm3 written in "
                           "a column in m3 (the import divides it by 1e6)", f"{_where(specimens, k)}: {vol}")
            if vol is not None and vol == 0:
                report.add("WARNING", "Units", "specimens.volume = 0", _where(specimens, k))
            if wt is not None and wt > 1.0:
                report.add("WARNING", "Units", "specimens.weight larger than 1 kg - probably grams written in a "
                           "column in kg", f"{_where(specimens, k)}: {wt}")
            if r.get("dir_dec") and not r.get("dir_tilt_correction"):
                report.add("WARNING", "Results", "specimens: direction without dir_tilt_correction (coordinate "
                           "system unknown)", _where(specimens, k))
        n_none = sum(1 for r in specimens.rows if not r.get("volume") and not r.get("weight"))
        if n_none:
            report.add("INFO", "Units", f"specimens without volume nor weight ({n_none}): magnetization per "
                       "unit volume/mass cannot be computed (the import assumes 10.8 cm3)", count=n_none, aggregate=True)


def _measurement_vector_issues(r: Dict[str, str], check_direction: bool = True) -> List[str]:
    x, y, z = (_num(r.get(c, "")) for c in ("magn_x", "magn_y", "magn_z"))
    dec, inc, mag = (_num(r.get(c, "")) for c in ("dir_dec", "dir_inc", "magn_moment"))
    if None in (x, y, z):
        return []
    issues = []
    norm = math.sqrt(x * x + y * y + z * z)
    if check_direction and dec is not None and inc is not None and norm > 0:
        dr, ir = math.radians(dec), math.radians(inc)
        c = (x * math.cos(ir) * math.cos(dr) + y * math.cos(ir) * math.sin(dr) + z * math.sin(ir)) / norm
        if math.degrees(math.acos(max(-1.0, min(1.0, c)))) > 1.0:
            issues.append("measurements: dir_dec/dir_inc more than 1 deg from the magn_x/y/z direction "
                          "(the import uses magn_x/y/z)")
    if mag is not None and abs(mag - norm) > 0.02 * max(abs(mag), norm):
        issues.append("measurements: magn_moment differs by more than 2% from |magn_x, magn_y, magn_z|")
    return issues


def _angle(d1, i1, d2, i2) -> float:
    r = math.radians
    c = math.sin(r(i1)) * math.sin(r(i2)) + math.cos(r(i1)) * math.cos(r(i2)) * math.cos(r(d1 - d2))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _check_direction_frame(tables: Dict[str, Table], report: Report):
    """Quand measurements porte a la fois magn_x/y/z (repere specimen) et
    dir_dec/dir_inc, retrouve dans QUEL repere sont ces directions en les
    recalculant avec l'orientation de l'echantillon (pmag.dogeo) et le
    pendage (pmag.dotilt) : specimen, geographique, tectonique - et, si
    rien ne colle, tectonique en prenant bed_dip_direction -/+ 90 (le cas
    d'une direction de couche (strike) saisie a la place de la direction
    de pendage, vu dans une contribution reelle)."""
    meas, specimens, samples = tables.get("measurements"), tables.get("specimens"), tables.get("samples")
    if meas is None or samples is None or "magn_x" not in meas.columns or "dir_dec" not in meas.columns:
        return
    try:
        from pmagpy import pmag
    except ImportError:
        return
    sample_of = {r.get("specimen", ""): r.get("sample", "") for r in specimens.rows} if specimens else {}
    sample_rows = {r.get("sample", ""): r for r in samples.rows}
    frames = collections.Counter()
    examples: Dict[str, str] = {}
    for k, r in enumerate(meas.rows):
        x, y, z = (_num(r.get(c, "")) for c in ("magn_x", "magn_y", "magn_z"))
        dec, inc = _num(r.get("dir_dec", "")), _num(r.get("dir_inc", ""))
        if None in (x, y, z, dec, inc) or x == y == z == 0:
            continue
        s_row = sample_rows.get(sample_of.get(r.get("specimen", ""), r.get("specimen", "")), {})
        az, dip = _num(s_row.get("azimuth", "")), _num(s_row.get("dip", ""))
        dsp = math.degrees(math.atan2(y, x)) % 360.0
        isp = math.degrees(math.atan2(z, math.hypot(x, y)))
        candidates = [("specimen", dsp, isp)]
        if az is not None and dip is not None:
            dg, ig = pmag.dogeo(dsp, isp, az, dip)
            candidates.append(("geographic", dg, ig))
            ddir, bdip = _num(s_row.get("bed_dip_direction", "")), _num(s_row.get("bed_dip", ""))
            if ddir is not None and bdip is not None:
                for label, d in (("tilt-corrected", ddir), ("tilt-corrected with bed_dip_direction - 90",
                                                           ddir - 90.0),
                                 ("tilt-corrected with bed_dip_direction + 90", ddir + 90.0)):
                    dt, it = pmag.dotilt(dg, ig, d % 360.0, bdip)
                    candidates.append((label, dt, it))
        best = min(candidates, key=lambda c: _angle(dec, inc, c[1], c[2]))
        frame = best[0] if _angle(dec, inc, best[1], best[2]) <= 2.0 else "none of them"
        frames[frame] += 1
        examples.setdefault(frame, _where(meas, k))
    total = sum(frames.values())
    if not total or frames.get("specimen", 0) == total:
        return
    for frame, n in frames.most_common():
        if frame == "specimen":
            continue
        if frame.startswith("tilt-corrected with"):
            msg = (f"measurements.dir_dec/dir_inc are tilt-corrected directions, but only match when the "
                   f"bedding is taken as {frame.split('with ')[1]}: samples.bed_dip_direction probably holds "
                   "the strike of the beds, not their dip direction")
            severity = "ERROR"
        elif frame == "none of them":
            msg = ("measurements.dir_dec/dir_inc match neither the magn_x/y/z direction nor its geographic or "
                   "tilt-corrected version (within 2 deg)")
            severity = "WARNING"
        else:
            msg = (f"measurements.dir_dec/dir_inc are in {frame} coordinates, not in specimen coordinates as "
                   "MagIC expects (magn_x/y/z are used by the import)")
            severity = "WARNING"
        report.add(severity, "Directions", msg + f" - {n}/{total} rows", examples[frame], count=n)


def _check_measurements(tables: Dict[str, Table], report: Report):
    meas = tables.get("measurements")
    if meas is None:
        report.add("WARNING", "Measurements", "no measurements table (results only)")
        return
    seq = collections.Counter(r.get("sequence") for r in meas.rows if r.get("sequence"))
    for value, n in seq.items():
        if n > 1:
            report.add("ERROR", "Measurements", "measurements.sequence must be unique", f"{value}: {n} rows")
    names = collections.Counter(r.get("measurement") for r in meas.rows if r.get("measurement"))
    for value, n in names.items():
        if n > 1:
            report.add("WARNING", "Measurements", "measurement name used several times", f"'{value}': {n} rows")

    steps: Dict[tuple, List[int]] = collections.defaultdict(list)
    nrm_specimens = set()
    specimens_seen = set()
    for k, r in enumerate(meas.rows):
        codes = _tokens(r.get("method_codes", ""))
        spec = r.get("specimen", "")
        specimens_seen.add(spec)
        if not any(c.startswith("LT-") for c in codes):
            report.add("WARNING", "Treatment codes",
                       "measurements without any LT- treatment code (only LP-/other codes): the type of step "
                       "must be guessed from treat_temp/treat_ac_field, and PmagPy's Demag GUI ignores these rows",
                       _where(meas, k))
        elif not any(c in _DEMAG_GUI_LT for c in codes) and any(c.startswith("LP-DIR") for c in codes):
            report.add("WARNING", "Treatment codes",
                       "directional (LP-DIR) measurement whose LT- code is none of LT-NO/LT-AF-Z/LT-T-Z/LT-M-Z/"
                       "LT-LT-Z: PmagPy's Demag GUI ignores it", f"{_where(meas, k)}: {r.get('method_codes')}")
        if "LT-HT-Z" in codes:
            report.add("WARNING", "Treatment codes", "LT-HT-Z used for zero-field heating steps - thermal "
                       "demagnetization is coded LT-T-Z (and the NRM LT-NO); PmagPy's Demag GUI and the "
                       "STARpaleomag_Py import do not recognize these steps", _where(meas, k))
        for c in codes:
            if c.startswith("DA-DIR"):
                report.add("WARNING", "Directions", f"measurements carrying {c}: dir_dec/dir_inc are not in "
                           "specimen coordinates as MagIC expects for measurements", _where(meas, k))
        # Avec une table samples, la direction est controlee (et son repere
        # identifie) par _check_direction_frame ; ici seulement le moment.
        for issue in _measurement_vector_issues(r, check_direction="samples" not in tables):
            report.add("WARNING", "Directions", issue, _where(meas, k))
        if not any(r.get(c) for c in ("magn_moment", "magn_volume", "magn_mass", "magn_x", "magn_uncal",
                                      "susc_chi_volume", "susc_chi_mass")):
            report.add("WARNING", "Measurements", "measurement with no magnetization value", _where(meas, k))
        if r.get("dir_dec") and not r.get("dir_inc") or r.get("dir_inc") and not r.get("dir_dec"):
            report.add("ERROR", "Measurements", "dir_dec without dir_inc (or the reverse)", _where(meas, k))
        temp, af, dc = _num(r.get("treat_temp", "")), _num(r.get("treat_ac_field", "")), _num(r.get("treat_dc_field", ""))
        if temp is not None and 0 < temp < 250 and "LT-LT-Z" not in codes:
            report.add("WARNING", "Units", "treat_temp below 250 K - probably degrees Celsius written in a "
                       "column in kelvin", f"{_where(meas, k)}: {temp}")
        if temp is not None and temp > 1300:
            report.add("WARNING", "Units", "treat_temp above 1300 K (1027 degC)", f"{_where(meas, k)}: {temp}")
        if af is not None and af > 1.0:
            report.add("WARNING", "Units", "treat_ac_field above 1 T - probably mT written in a column in "
                       "tesla", f"{_where(meas, k)}: {af}")
        if dc is not None and abs(dc) > 10.0:
            report.add("WARNING", "Units", "treat_dc_field above 10 T - probably mT or uT written in a column "
                       "in tesla", f"{_where(meas, k)}: {dc}")
        if dc is not None and any(c in codes for c in ("LT-T-I", "LT-PTRM-I")) and abs(dc) > 1e-3:
            report.add("WARNING", "Units", "paleointensity laboratory field above 1000 uT - unit error?",
                       f"{_where(meas, k)}: {dc} T")
        is_nrm = "LT-NO" in codes or (not temp or temp <= 300) and not af and not dc and \
            not any(c.startswith(("LT-T", "LT-AF", "LT-IRM", "LT-M")) for c in codes)
        if is_nrm:
            nrm_specimens.add(spec)
        lt = ",".join(sorted(c for c in codes if c.startswith("LT-")))
        key = (spec, r.get("experiment", ""), lt, r.get("treat_temp", ""), r.get("treat_ac_field", ""),
               r.get("treat_dc_field", ""), r.get("treat_dc_field_phi", ""), r.get("treat_dc_field_theta", ""),
               r.get("treat_step_num", "") if lt == "LT-NO" else "")
        steps[key].append(k)
    for key, ks in steps.items():
        if len(ks) < 2:
            continue
        rows = [meas.rows[k] for k in ks]
        vectors = {tuple(r.get(c, "") for c in ("magn_moment", "magn_x", "magn_y", "magn_z", "dir_dec",
                                                    "dir_inc")) for r in rows}
        same = len(vectors) == 1
        report.add("WARNING" if same else "INFO", "Repeated steps",
                   "the same treatment step measured several times on one specimen with identical values "
                   "(duplicated rows)" if same else
                   "the same treatment step measured several times on one specimen (repeated measurement - "
                   "the import keeps them all, in file order)",
                   f"{key[0]} {key[2] or 'no LT code'} temp={key[3] or '-'} ac={key[4] or '-'} "
                   f"dc={key[5] or '-'}: {len(ks)} rows (lines {', '.join(str(meas.lines[k]) for k in ks[:5])})")
    without_nrm = sorted(specimens_seen - nrm_specimens - {""})
    if without_nrm:
        report.add("INFO", "Measurements", f"specimens with no NRM / untreated measurement ({len(without_nrm)})",
                   ", ".join(without_nrm[:12]) + (" ..." if len(without_nrm) > 12 else ""), count=len(without_nrm), aggregate=True)


def check_contribution(path: str) -> Report:
    model = data_model()
    report = Report(path)
    tables_list = read_tables(path)
    tables: Dict[str, Table] = {}
    for t in tables_list:
        if t.name in tables:  # table repetee : fusionnee pour les controles de contenu
            tables[t.name].rows.extend(t.rows)
            tables[t.name].lines.extend(t.lines)
        else:
            tables[t.name] = t
    with open(path, "rb") as f:
        if f.read(3) == b"\xef\xbb\xbf":
            report.add("INFO", "Structure", "the file starts with a UTF-8 byte-order mark (BOM) - invisible, "
                       "but some programs then fail to recognize the first table")
    contrib = tables.get("contribution")
    if contrib and contrib.rows:
        r = contrib.rows[0]
        report.summary.append(
            f"  contribution {r.get('id', '?')} version {r.get('version', '?')}, {r.get('contributor', '?')}, "
            f"{r.get('timestamp', '?')[:10]}, reference {r.get('reference', '?')}, data model "
            f"{r.get('data_model_version', '?')}")
        if r.get("data_model_version", "3.0") not in ("3.0", "3"):
            report.add("WARNING", "Structure", f"data model version {r.get('data_model_version')} (checks below "
                       "use MagIC 3.0)")
    report.summary.append(f"  reference used for the checks: MagIC data model {model.version} updated "
                          f"{model.updated}, {model.origin}; method codes: {model.codes_origin}; controlled "
                          "vocabularies: "
                          + (f"{', '.join(model.updated_vocabularies)} from the official lists (October 2026), "
                             "the others" if model.updated_vocabularies else "all")
                          + " as shipped with PmagPy (2019)")

    index: Dict[str, set] = {}
    for t in tables.values():
        for col in t.columns:
            index[f"{t.name}.{col}"] = {r.get(col, "") for r in t.rows}

    _check_structure(tables_list, model, report)
    _check_columns(list(tables.values()), model, report, index)
    _check_duplicates(tables, report)
    _check_hierarchy(tables, report)
    _check_samples_sites(tables, report)
    _check_measurements(tables, report)
    _check_direction_frame(tables, report)
    return report


def report_path_for(path_in: str, path_out: Optional[str] = None) -> str:
    base = os.path.splitext(path_out or path_in)[0]
    return base + "_MagIC_check.txt"


def write_check_report(path_in: str, path_out: Optional[str] = None) -> Tuple[str, Dict[str, int]]:
    """Controle `path_in` et ecrit le rapport a cote de `path_out` (ou du
    fichier d'entree). Retourne (chemin du rapport, nombre de types de
    problemes par gravite)."""
    report = check_contribution(path_in)
    out = report_path_for(path_in, path_out)
    with open(out, "w", encoding="utf-8") as f:
        f.write(report.format())
    return out, report.counts()


def is_magic_contribution(path: str) -> bool:
    try:
        with open(path, encoding="utf-8-sig", errors="replace") as f:
            first = f.readline()
    except OSError:
        return False
    return first.lower().startswith("tab") and "\t" in first


def check_folder(folder: str, out_dir: Optional[str] = None) -> Tuple[str, List[str]]:
    """Controle TOUTES les contributions MagIC d'un dossier (fichiers .txt
    dont la 1re ligne est "tab delimited...", sous-dossiers non compris) -
    demande explicite utilisateur ("can you run the check for all the
    contributions in the folder"). Un rapport par contribution + un tableau
    recapitulatif dans `out_dir` (defaut <dossier>/MagIC_checks), et un
    controle ENTRE contributions : meme reference (DOI) ou memes noms de
    specimens dans deux contributions (le meme article depose deux fois,
    vu dans un vrai dossier). Retourne (chemin du recapitulatif, lignes de
    resume)."""
    out_dir = out_dir or os.path.join(folder, "MagIC_checks")
    os.makedirs(out_dir, exist_ok=True)
    paths = sorted(os.path.join(folder, f) for f in os.listdir(folder)
                   if f.lower().endswith(".txt") and not f.endswith("_MagIC_check.txt")
                   and is_magic_contribution(os.path.join(folder, f)))
    table = ["contribution\terrors\twarnings\tnotes\tpositive dip (samples)\treference\terror kinds"]
    references: Dict[str, List[str]] = collections.defaultdict(list)
    specimens_of: Dict[str, set] = {}
    lines = []
    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        try:
            report = check_contribution(path)
        except Exception as e:
            lines.append(f"{name}: check failed ({e})")
            continue
        with open(os.path.join(out_dir, name + "_MagIC_check.txt"), "w", encoding="utf-8") as f:
            f.write(report.format())
        tables = {t.name: t for t in read_tables(path)}
        contrib = tables.get("contribution")
        ref = (contrib.rows[0].get("reference", "") if contrib and contrib.rows else "").upper()
        if ref:
            references[ref].append(name)
        if "specimens" in tables:
            specimens_of[name] = {r.get("specimen", "") for r in tables["specimens"].rows} - {""}
        c = report.counts()
        errors = "; ".join(k[2] for k in report.issues if k[0] == "ERROR")
        n_pos, n_dip = report.positive_dip
        dip_col = f"{n_pos}/{n_dip}" if n_dip else "no dip"
        table.append(f"{name}\t{c['ERROR']}\t{c['WARNING']}\t{c['INFO']}\t{dip_col}\t{ref}\t{errors}")
        lines.append(f"{name}: {c['ERROR']} error(s), {c['WARNING']} warning(s), {c['INFO']} note(s)")

    across = []
    for ref, names in references.items():
        if len(names) > 1:
            across.append(f"same reference {ref} in: {', '.join(names)}")
    names = sorted(specimens_of)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            common = specimens_of[a] & specimens_of[b]
            if common and len(common) >= 0.5 * min(len(specimens_of[a]), len(specimens_of[b])):
                across.append(f"{a} and {b} share {len(common)} specimen names "
                              f"({len(specimens_of[a])} and {len(specimens_of[b])} specimens)")
    summary = os.path.join(out_dir, "00_summary.tsv")
    with open(summary, "w", encoding="utf-8") as f:
        f.write("\n".join(table) + "\n")
        if across:
            f.write("\n# Across contributions (probably the same data deposited twice)\n")
            f.write("\n".join("# " + a for a in across) + "\n")
    n_err = sum(1 for row in table[1:] if row.split("\t")[1] != "0")
    head = [f"{len(paths)} MagIC contribution(s) checked in {folder}: {n_err} with at least one error.",
            f"Reports and summary table -> {out_dir}"]
    if across:
        head.append("Across contributions (probably the same data deposited twice):")
        head.extend("  " + a for a in across)
    return summary, head + [""] + lines


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        print(check_contribution(p).format())
