"""
Age d'un specimen dans le .prmag : quatre champs separes, comme dans MagIC -
demande explicite utilisateur ("pour lever les ambiguites, ne faudrait-il
pas mettre age, age_low, age_high et age_unit") :

    formation: Cerro Toro Fm.	age: 80	age_low: 70	age_high: 90	age_unit: Ma	geologic_classes: ...

Avant, un seul champ texte `age:` portait tour a tour "80", "80 ± 10 Ma",
"70 - 90 Ma" ou "Turonian", que l'export devait deviner. Les anciens
fichiers restent lisibles : quand age_low/age_high/age_unit sont absents,
ils sont deduits du texte de `age:` (split_age_text), et ecrits au prochain
enregistrement.
"""

import re
from typing import Tuple

_AGE_UNIT_PATTERNS = [
    # Noms MagIC COMPLETS d'abord : une chaine importee de MagIC ecrit
    # "665 # 20 AD Years Cal AD (+/-)" ; avec seul "Cal AD" reconnu, il
    # restait "AD Years ... (+/-)" colle a l'ecart-type.
    ("Years Cal AD (+/-)", "Years Cal AD (+/-)"),
    ("Years Cal BP", "Years Cal BP"),
    ("Years AD (+/-)", "Years AD (+/-)"),
    ("Years BP", "Years BP"),
    ("Ga", "Ga"),
    ("Ma", "Ma"),
    ("ka", "ka"),
    ("Cal AD", "Years Cal AD (+/-)"),
    ("Cal BP", "Years Cal BP"),
    ("AD", "Years AD (+/-)"),
    ("BP", "Years BP"),
]

_ND = ("", "n.d", "N.D", "nd")


def _clean(value) -> str:
    v = ("" if value is None else str(value)).strip()
    return "" if v in _ND else v


def parse_age(age_str: str) -> Tuple[str, str, str, str, str]:
    """Texte libre -> (age, age_sigma, age_low, age_high, age_unit)."""
    s = (age_str or "").strip()
    age_unit = ""
    for token, magic_unit in _AGE_UNIT_PATTERNS:
        if token in s:
            age_unit = magic_unit
            s = s.replace(token, " ").strip()
            if magic_unit.startswith("Years"):
                s = re.sub(r"\b(?:AD|BP|Cal)\b", " ", s)
                s = re.sub(r"\s+", " ", s).strip()
            break

    age = age_sigma = age_low = age_high = ""
    # "60 ± 5 Ma" / "60 +/- 5 Ma" : age et incertitude, comme "60#5".
    pm = re.match(r"^\s*([-+]?\d+(?:\.\d+)?)\s*(?:±|\+/-|\+-)\s*(\d+(?:\.\d+)?)\s*$", s)
    if pm:
        age, age_sigma = pm.group(1), pm.group(2)
    elif "#" in s:
        a, b = s.split("#", 1)
        age, age_sigma = a.strip(), b.strip()
    elif "&" in s:
        parts = [p.strip() for p in s.split("&")]
        if len(parts) >= 3:
            age, age_low, age_high = parts[0], parts[1], parts[2]
        else:
            age = s
    elif "_" in s:
        a, b = s.split("_", 1)
        age_low, age_high = a.strip(), b.strip()
    elif re.search(r"\s-\s", s):
        a, b = re.split(r"\s-\s", s, maxsplit=1)
        age_low, age_high = a.strip(), b.strip()
    else:
        age = s

    return age, age_sigma, age_low, age_high, age_unit


def _num(s: str):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _fmt(x: float) -> str:
    return f"{x:g}"


def split_age_text(age_str: str) -> Tuple[str, str, str, str]:
    """Texte libre d'un ancien `age:` -> (age, age_low, age_high, age_unit).
    "80 ± 10 Ma" -> ("80", "70", "90", "Ma") ; "70 - 90 Ma" -> ("", "70",
    "90", "Ma") ; "80" -> ("80", "", "", "") ; "Turonian" reste dans age."""
    age, sigma, low, high, unit = parse_age(_clean(age_str))
    a, s = _num(age), _num(sigma)
    if a is not None and s is not None and not (low or high):
        low, high = _fmt(a - s), _fmt(a + s)
    return age, low, high, unit


def normalize_age_fields(age="", age_low="", age_high="", age_unit="") -> Tuple[str, str, str, str]:
    """Quatre champs propres. Si age_low/age_high/age_unit sont tous vides,
    ils sont deduits du texte de `age` (ancien format). Si seule l'unite
    est collee a l'age ("80 Ma"), elle en est separee."""
    age, age_low, age_high, age_unit = (_clean(v) for v in (age, age_low, age_high, age_unit))
    if not (age_low or age_high or age_unit):
        return split_age_text(age)
    if age and _num(age) is None:
        a2, l2, h2, u2 = split_age_text(age)
        if _num(a2) is not None or l2 or h2:
            age = a2
            age_low, age_high, age_unit = age_low or l2, age_high or h2, age_unit or u2
    return age, age_low, age_high, age_unit


def age_fields_chunks(age="", age_low="", age_high="", age_unit="") -> str:
    """'age: ..\\tage_low: ..\\tage_high: ..\\tage_unit: ..' pour une
    ligne d'entete .prmag (n.d pour un champ vide)."""
    vals = normalize_age_fields(age, age_low, age_high, age_unit)
    return "\t".join(f"{k}: {v or 'n.d'}" for k, v in zip(("age", "age_low", "age_high", "age_unit"), vals))


def set_age_fields_in_line(line: str, age="", age_low="", age_high="", age_unit="") -> str:
    """Remplace dans une ligne 'formation: ...\\tage: ...' les champs d'age
    par les quatre champs, places juste apres `age` (ou apres `formation`
    s'il n'y a pas de champ age) - le reste de la ligne est inchange."""
    chunks = [c for c in line.split("\t")
              if c.partition(":")[0].strip() not in ("age", "age_low", "age_high", "age_unit")]
    new = age_fields_chunks(age, age_low, age_high, age_unit).split("\t")
    pos = 1 if chunks and chunks[0].partition(":")[0].strip() == "formation" else 0
    return "\t".join(chunks[:pos] + new + chunks[pos:])


def age_display(age="", age_low="", age_high="", age_unit="") -> str:
    """Texte lisible : "80 (70-90) Ma", "70-90 Ma", "80 Ma"..."""
    age, low, high, unit = normalize_age_fields(age, age_low, age_high, age_unit)
    rng = f"{low}-{high}" if (low or high) else ""
    if age and rng:
        txt = f"{age} ({rng})"
    else:
        txt = age or rng
    return f"{txt} {unit}".strip() if txt else ""
