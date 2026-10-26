"""
Utilitaires de reduction/conversion de donnees du menu Graphics de
StarmacOSX_x.f95 (convertthelli, removestep, elimineGRM, convzmoins,
exportthellier) - toutes MUTENT `ech.mesures` en place, comme leurs
equivalents Fortran (aucune ne dessine quoi que ce soit ; elles vivent dans
le menu Graphics de l'original mais sont en realite des operations sur les
donnees, generalement suivies d'un nouveau trace Zijderveld manuel).

`testinduite` (correction du champ induit parasite, calcul.f:3687-3774) a
ete porte puis RETIRE a la demande explicite de l'utilisateur ("on peut la
retirer, personne ne va penser a faire ces tests pour tester le champ
residuel au niveau de la mesure du Cryo 2G") - fonctionnalite jugee sans
usage reel pour ce laboratoire, pas seulement deplacee de menu.
"""

import math
import os
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Tuple

from selection import SelectedSample, polere
from calcul import _fortran_e


# ---------------------------------------------------------------------------
# convertthelli / thelli : sequence Thellier (D/N/T/K/S/R/P/L/X/Y/Z/V) ->
# sequence NRM/TRM simple (cod1='D', cod2='p') pour le trace Zijderveld -
# plotorthog.f:780-921
# ---------------------------------------------------------------------------

def convert_thellier_to_nrm(ech: SelectedSample) -> None:
    """Port exact du `select case (mes(i).cod1)` de `convertthelli`
    (plotorthog.f:800-887). D/N/T/K/S copies telles quelles (cod1='D',
    cod2='p'). R/P/L supprimees. Pour X/Y/Z avec cod2 '+'/'-' : fusionnee
    (moyenne) avec la ligne PRECEDENTE seulement si celle-ci a le MEME
    cod1 (2e passage du meme axe) ; sinon supprimee - ce n'est jamais la
    ligne SUIVANTE qui declenche la fusion. Pour V : fusionnee TOUJOURS
    avec la ligne precedente, QUEL QUE SOIT son cod1 - typiquement un 'R'
    juste avant (meme si ce 'R' est par ailleurs supprime de la sortie) :
    c'est ainsi que R et V sont combines pour extraire la NRM (TRM
    partielle appliquee/mesuree en deux positions Z/-Z ; R porte le pas
    thermique, V le controle pTRM, leur moyenne annule le biais directionnel
    entre les deux). Toute autre cod1 : supprimee (case default)."""
    mesures = ech.mesures
    result = []

    def merge(prev, cur):
        x = (prev.x + cur.x) / 2.0
        y = (prev.y + cur.y) / 2.0
        z = (prev.z + cur.z) / 2.0
        q = (prev.q + cur.q) // 2
        s = prev.s if (cur.s == 0 or prev.s == 0) else (cur.s + prev.s) / 2.0
        result.append(replace(cur, x=x, y=y, z=z, q=q, s=s, cod1="D", cod2="p"))

    for i, m in enumerate(mesures):
        if m.cod1 in ("D", "N", "T", "K", "S"):
            result.append(replace(m, cod1="D", cod2="p"))
        elif m.cod1 in ("R", "P", "L"):
            continue
        elif m.cod1 in ("X", "Y", "Z"):
            if m.cod2 in ("+", "-") and i > 0 and mesures[i - 1].cod1 == m.cod1:
                merge(mesures[i - 1], m)
        elif m.cod1 == "V":
            if i > 0:
                merge(mesures[i - 1], m)
    ech.mesures = result


# ---------------------------------------------------------------------------
# removestep : supprime toutes les lignes d'un etape donne, renomme les
# cod2 de R/V/P sequentiellement - plotpaleoint2.f:1918-1979
# ---------------------------------------------------------------------------

_CLEANPALEO_CODES = ("X", "Y", "Z", "Q", "L", "F")


def remove_step(ech: SelectedSample, etape: int) -> int:
    """D'abord `cleanpaleo` sur X/Y/Z/Q/L/F (supprime les lignes GRM/qualite
    non pertinentes pour l'interpretation paleointensite), puis supprime
    toutes les lignes de l'etape donne. Renumerotation par ETAPE - reecrite
    proprement en Python a la demande explicite de l'utilisateur ("the menu
    remove step is supposed to be dedicated to paleointensity which are
    numbered RA,VA,...RI,VI, etc ... I tried to do that in Fortran but it
    is possible to do it cleanly in Python").

    Une lettre par etape DISTINCTE (ordre d'apparition parmi R/V,
    mesures[0] exclue), partagee par R et V de cette etape - remplace les
    DEUX compteurs INDEPENDANTS du `removestep` Fortran original
    (plotpaleoint2.f:1958-1975 : A,B,C... pour R, separement A,B,C... pour
    V), qui ne restent synchronises que si R et V sont en nombre egal et
    dans le meme ordre.

    'P' (verification pTRM) reprend la lettre du DERNIER R/V renumerote
    rencontre (equivalent, pour ce format de donnees ou X/Y/Z/Q/L/F ont
    deja ete retires, au `mes(i).cod2=mes(i-1).cod2` du Fortran) - PAS la
    lettre de sa propre etape revisitee (`m.etape`) : verifie sur donnees
    reelles (06A, SanJuan_Pmag.prmag) que le cod2 d'un P a l'origine
    correspond au point de reference COURANT du diagramme d'Arai (pour
    xt, cf. paleointensity.py `check.xt = vecdiff(ji, i)` apparie par
    cod2==cod2 avec un V), independant de l'etape re-testee - `m.etape`
    sert separement a retrouver le R de l'etape ciblee (`check.yt`/
    `xtptrm`, apparie par etape==etape). Retourne le nombre de lignes
    supprimees (cleanpaleo + etape)."""
    before = len(ech.mesures)
    ech.mesures = [m for m in ech.mesures if m.cod1 not in _CLEANPALEO_CODES]
    ech.mesures = [m for m in ech.mesures if m.etape != etape]
    removed = before - len(ech.mesures)

    mesures = ech.mesures
    etape_to_letter: Dict[float, str] = {}
    last_cod2 = None
    for m in mesures[1:]:
        if m.cod1 in ("R", "V"):
            letter = etape_to_letter.get(m.etape)
            if letter is None:
                ic = len(etape_to_letter)
                letter = chr(65 + ic) if ic < 26 else "?"
                etape_to_letter[m.etape] = letter
            m.cod2 = letter
            last_cod2 = letter
        elif m.cod1 == "P" and last_cod2 is not None:
            m.cod2 = last_cod2
    return removed


def remove_bad_quality_steps(ech: SelectedSample) -> Tuple[int, List[float]]:
    """Detecte les etapes marquees qualite 'b' (mauvais) dans ech.mesures -
    colonne "quality" du .prmag, voir testlect.Measurement.quality - et
    les supprime UNE PAR UNE via remove_step (meme renumerotation des
    cod2 R/V/P que la suppression manuelle d'etape, meme suppression en
    cascade d'une eventuelle verification pTRM qui se referait a l'etape
    supprimee - remove_step filtre deja sur `m.etape`, qui est justement
    l'etape CIBLEE par un 'P', pas son propre ordre d'execution) - demande
    explicite utilisateur ("je voudrais utiliser le critere de qualite
    b/g dans les donnees pour ne pas prendre en compte cette etape...
    un des problemes en paleointensite est eventuellement qu'une serie
    d'echantillons ne soient pas mis correctement dans le four... il faut
    alors renumeroter les etapes automatiquement et aussi supprimer la
    PTRM check si elle se refere a l'etape rejetee"). Retourne (nombre
    total de lignes supprimees, etapes retirees, triees) - liste vide si
    aucune etape 'b' trouvee (aucune mutation dans ce cas)."""
    bad_etapes = sorted({
        m.etape for m in ech.mesures if (m.quality or "g").strip().lower() == "b"
    })
    total_removed = 0
    for etape in bad_etapes:
        total_removed += remove_step(ech, etape)
    return total_removed, bad_etapes


# ---------------------------------------------------------------------------
# elimineGRM (elimineGRM_DZ / elimineGRM_ZD) : combine les triplets
# cod2='X','Y','Z' consecutifs (apres un pas 'F') en un seul point corrige -
# plotorthog.f:1536-1846
# ---------------------------------------------------------------------------

def eliminate_grm(ech: SelectedSample, method: int = 1) -> int:
    """Port exact de `elimineGRM_DZ` (method=1, plotorthog.f:1692-1846) /
    `elimineGRM_ZD` (method=2, plotorthog.f:1536-1691), selectionnees par
    le dispatcher `elimineGRM` (StarmacOSX_x.f95:975) selon la commune
    `/prefstarmac/grm`. Seules les lignes cod1='F' sont concernees : celles
    dont cod2 est '+'/'-'/'=' ou dont l'etape vaut 0 sont copiees telles
    quelles ; un triplet consecutif cod2='X','Y','Z' est reduit en UN
    point (method=1 : substitution axe-par-axe x du point X, y du point Y,
    z du point Z ; method=2 : moyenne simple des 3 points) ; toute autre
    ligne 'F' isolee (ni copiee ni fusionnable) est supprimee, comme dans
    le Fortran (aucune des deux conditions du `select case` ne s'applique).
    Toutes les lignes d'un autre cod1 sont copiees inchangees. Retourne le
    nombre de triplets reduits."""
    mesures = ech.mesures
    n = len(mesures)
    result = []
    i = 0
    count = 0
    while i < n:
        m = mesures[i]
        if m.cod1 != "F":
            result.append(m)
            i += 1
            continue
        if m.cod2 in ("+", "-", "=") or m.etape == 0:
            result.append(m)
            i += 1
        elif (i + 2 < n and m.cod2 == "X" and mesures[i + 1].cod2 == "Y"
                and mesures[i + 2].cod2 == "Z"):
            mx, my, mz = mesures[i], mesures[i + 1], mesures[i + 2]
            if method == 1:
                x, y, z = mx.x, my.y, mz.z
            else:
                x = (mx.x + my.x + mz.x) / 3.0
                y = (mx.y + my.y + mz.y) / 3.0
                z = (mx.z + my.z + mz.z) / 3.0
            result.append(replace(mx, x=x, y=y, z=z, cod2="="))
            count += 1
            i += 3
        else:
            i += 1  # ligne 'F' isolee : ni copiee ni fusionnee (comme le Fortran)
    ech.mesures = result
    return count


def has_grm_triplets(ech: SelectedSample) -> bool:
    """Vrai si le specimen contient au moins un triplet AF consecutif
    cod1='F', cod2 'X','Y','Z' (celui que eliminate_grm reduit)."""
    ms = ech.mesures
    return any(
        ms[i].cod1 == "F" and ms[i].cod2 == "X" and ms[i + 1].cod2 == "Y" and ms[i + 2].cod2 == "Z"
        for i in range(len(ms) - 2)
    )


def reduce_grm_before_fit(samples, method: int = 1) -> List[Tuple[str, int, int]]:
    """Reduit les triplets AF X/Y/Z de chaque specimen qui en contient,
    AVANT un ajustement (best fit) - demande explicite utilisateur ("au
    cours de l'analyse en best fit, reduire les GRM d'abord si il y a des
    mesures en FX, FY et FZ"). Methode 1 (substitution axe par axe) par
    defaut : c'est elle qui redonne les ajustements des anciens .r du
    Fortran (Fuego_Pmag.r). Retourne [(id, triplets reduits, mesures
    restantes)] pour les specimens modifies."""
    done = []
    for ech in samples:
        if has_grm_triplets(ech):
            n = eliminate_grm(ech, method=method)
            done.append((ech.id, n, len(ech.mesures)))
    return done


# ---------------------------------------------------------------------------
# detect_grm : diagnostic (pas de correction) de contamination par GRM sur
# un degausser 3 axes EN LIGNE avec le magnetometre - demande explicite
# utilisateur ("Most laboratories using 3 axis degausser online with the
# Cryogenic magnetometer... A GRM is a laboratory magnetization, orthogonal
# of the axis of the coil of the AF demagnetization that increase with AF
# field... Can you write a test to detect such behavior especially for the
# Magic database contribution").
#
# Physique (Stephenson 1993 ; Dankers & Zijderveld 1981) : la GRM acquise
# lors d'un pas de demagnetisation AF est orientee PERPENDICULAIREMENT au
# DERNIER axe de la bobine utilise a ce pas, et sa magnitude CROIT avec le
# champ AF applique - surtout marquee sur les grains mono-domaine avec une
# legere anisotropie magnetique. La plupart des laboratoires utilisent
# TOUJOURS la meme sequence d'axes (donc le meme axe "dernier"), ce qui rend
# la GRM invisible en isolation (elle se confond avec la direction NRM sur
# tout le trajet de desaimantation) - contrairement a Rennes, qui alterne
# deliberement 2 sequences opposees (cod2 '+'=Y,Z,X et '-'=X,Z,Y sur
# l'instrument "C1"/"C2", voir magic_export._measurement_treatment) d'un
# palier au suivant PRECISEMENT pour reveler la GRM (son signe/sa direction
# s'inverse avec la sequence, alors qu'une desaimantation "propre" ne
# devrait pas en dependre) - c'est le meme mecanisme que corrige
# `eliminate_grm` ci-dessus (mais celui-ci ne traite QUE les triplets X/Y/Z
# consecutifs a un instrument sans degausser en ligne, jamais le cas '+'/'-'
# d'un C1/C2 - AUCUNE correction n'existe aujourd'hui pour ce cas, d'ou
# l'utilite de ce diagnostic).
#
# Methode (revisee - la premiere version, comparant simplement le vecteur
# difference entre paliers consecutifs a l'axe de bobine, ne distinguait pas
# une GRM d'une simple decroissance NATURELLEMENT oblique par rapport aux
# axes de bobine (verifie : un cas synthetique PROPRE, sans aucune GRM mais
# de direction non alignee sur X/Y/Z, produisait un faux positif tout aussi
# fort qu'un cas contamine) :
#
# 1) Ajuste une droite (ACP) a travers TOUS les points (x,y,z) des paliers
#    AF exploitables - la tendance de decroissance "propre" attendue.
# 2) Pour chaque palier, le RESIDU (ecart au point projete sur cette droite)
#    capture tout ce que la decroissance simple n'explique pas.
# 3) Decompose ce residu en une composante ALONG le dernier axe de bobine
#    utilise a ce palier et une composante PERP (orthogonale a cet axe) - la
#    GRM, PAR DEFINITION, n'a aucune composante le long de l'axe qui l'a
#    generee ; `perp_fraction` (part du residu qui est bien perpendiculaire
#    a l'axe) verifie cette signature GEOMETRIQUE SPECIFIQUE, ce qui exclut
#    la plupart des causes alternatives de non-linearite (bruit, 2e
#    composante) qui n'ont pas de raison de s'aligner sur les axes de bobine.
# 4) La magnitude de cette composante perp (normalisee par la NRM initiale
#    `arno`, meme convention que paleointensity.py) doit CROITRE avec le
#    champ AF (correlation de Pearson) - signature temporelle attendue de la
#    GRM - et atteindre une fraction non negligeable de la NRM au palier le
#    plus fort.
#
# Fonctionne QUE le laboratoire alterne (Rennes, cod2 '+'/'-') ou utilise
# TOUJOURS la meme sequence (la plupart des laboratoires, dixit
# l'utilisateur) : aucune des deux etapes ci-dessus ne suppose une
# alternance, contrairement a une comparaison directe entre groupes '+' et
# '-' (qui ne fonctionnerait QUE pour Rennes) - c'est precisement le point
# souleve par l'utilisateur ("especially for the Magic database
# contribution", des donnees d'autres laboratoires qui n'alternent pas).
# Seuils (`correlation_threshold`/`ratio_threshold`/`perp_fraction_min`)
# HEURISTIQUES, pas un critere physique absolu : ce diagnostic vise a
# ATTIRER L'ATTENTION pour inspection visuelle (Zijderveld), pas a rejeter
# automatiquement des donnees.
# ---------------------------------------------------------------------------

# Vecteur unitaire du DERNIER axe de la bobine utilise a ce palier - X/Y/Z
# pour un pas mono-axe, '+'/'-' pour la sequence complete 3 axes d'un
# instrument "C1"/"C2" (voir magic_export.py : '+' = Y,Z,X (dernier axe X),
# '-' = X,Z,Y (dernier axe Y), confirme par l'utilisateur cette session).
_AXIS_UNIT: Dict[str, tuple] = {
    "X": (1.0, 0.0, 0.0),
    "Y": (0.0, 1.0, 0.0),
    "Z": (0.0, 0.0, 1.0),
    "+": (1.0, 0.0, 0.0),
    "-": (0.0, 1.0, 0.0),
}


def _pearson(xs: List[float], ys: List[float]) -> float:
    n = len(xs)
    mean_x, mean_y = sum(xs) / n, sum(ys) / n
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    return cov / math.sqrt(var_x * var_y) if var_x > 0 and var_y > 0 else 0.0


@dataclass
class GrmCheckResult:
    specimen: str
    n_steps: int            # nombre de paliers AF exploitables
    correlation: float      # correlation de Pearson (champ AF, composante perp/NRM)
    max_ratio: float        # composante perp/NRM au palier le plus fort
    mean_perp_fraction: float  # part moyenne du residu (vs droite ACP) qui est perp a l'axe
    suspected: bool         # les 3 grandeurs ci-dessus au-dessus des seuils


def detect_grm(
    ech: SelectedSample,
    min_steps: int = 4,
    correlation_threshold: float = 0.6,
    ratio_threshold: float = 0.02,
    perp_fraction_min: float = 0.7,
) -> Optional[GrmCheckResult]:
    """Retourne None si le specimen n'a pas assez de paliers AF exploitables
    (moins de `min_steps` sur un instrument "C1"/"C2" - un degausser tumbler
    hors ligne n'a pas d'axe "dernier" bien defini, voir docstring
    ci-dessus, et est donc exclu). Sinon un `GrmCheckResult` dont
    `suspected` indique si CE specimen merite une inspection visuelle
    (Zijderveld) pour contamination GRM probable."""
    mesures = ech.mesures
    af_steps = [
        m for m in mesures
        if m.cod1 == "F" and m.cod2 in _AXIS_UNIT
        and (m.ins or "").strip().upper().startswith("C")
    ]
    n = len(af_steps)
    if n < min_steps:
        return None

    ideb = next((m for m in mesures if m.cod1 == "N"), mesures[0] if mesures else None)
    if ideb is None:
        return None
    arno = math.sqrt(ideb.x ** 2 + ideb.y ** 2 + ideb.z ** 2)
    if arno == 0.0:
        return None

    # La droite ACP de reference est ajustee sur la PREMIERE MOITIE (bas
    # champ) des paliers exploitables SEULEMENT, pas sur tous - bug de
    # conception corrige pendant le developpement de cette fonction (verifie
    # sur donnees synthetiques) : une ACP sur l'ENSEMBLE des points laisse
    # la GRM elle-meme (grande a haut champ) DEFORMER la droite ajustee, ce
    # qui absorbe une partie de la deviation qu'on cherche justement a
    # detecter et rend le residu haut-champ artificiellement petit/non
    # monotone. En n'ajustant que sur le bas champ (ou la GRM est encore
    # negligeable), la droite reflete la VRAIE direction de decroissance,
    # et le residu des paliers plus forts (mesure contre CETTE droite)
    # grandit alors proprement avec le champ si une GRM est presente.
    ref_n = max(3, n // 2)
    ref_pts = [(m.x, m.y, m.z) for m in af_steps[:ref_n]]
    cx = sum(p[0] for p in ref_pts) / ref_n
    cy = sum(p[1] for p in ref_pts) / ref_n
    cz = sum(p[2] for p in ref_pts) / ref_n
    ref_centered = [(p[0] - cx, p[1] - cy, p[2] - cz) for p in ref_pts]

    # Direction principale (plus grande variance) via la matrice de
    # covariance 3x3 - evite une dependance a numpy pour ce seul calcul,
    # coherent avec le reste de datatools.py (aucune autre fonction du
    # fichier n'importe numpy).
    sxx = sum(p[0] * p[0] for p in ref_centered)
    syy = sum(p[1] * p[1] for p in ref_centered)
    szz = sum(p[2] * p[2] for p in ref_centered)
    sxy = sum(p[0] * p[1] for p in ref_centered)
    sxz = sum(p[0] * p[2] for p in ref_centered)
    syz = sum(p[1] * p[2] for p in ref_centered)
    cov_mat = [[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]]
    direction = _dominant_eigenvector(cov_mat)
    if direction is None:
        return None

    fields: List[float] = []
    perp_ratios: List[float] = []
    perp_fractions: List[float] = []
    for m in af_steps:
        p = (m.x - cx, m.y - cy, m.z - cz)
        dot_dir = p[0] * direction[0] + p[1] * direction[1] + p[2] * direction[2]
        resid = (
            p[0] - dot_dir * direction[0],
            p[1] - dot_dir * direction[1],
            p[2] - dot_dir * direction[2],
        )
        resid_mag = math.sqrt(sum(c * c for c in resid))
        if resid_mag == 0.0:
            continue
        ax, ay, az = _AXIS_UNIT[m.cod2]
        along = resid[0] * ax + resid[1] * ay + resid[2] * az
        perp = (resid[0] - along * ax, resid[1] - along * ay, resid[2] - along * az)
        perp_mag = math.sqrt(sum(c * c for c in perp))
        fields.append(float(m.etape))
        perp_ratios.append(perp_mag / arno)
        perp_fractions.append(perp_mag / resid_mag)

    if len(fields) < min_steps:
        return None

    correlation = _pearson(fields, perp_ratios)
    max_ratio = max(perp_ratios)
    mean_perp_fraction = sum(perp_fractions) / len(perp_fractions)
    suspected = (
        correlation >= correlation_threshold
        and max_ratio >= ratio_threshold
        and mean_perp_fraction >= perp_fraction_min
    )

    return GrmCheckResult(
        specimen=ech.id, n_steps=len(fields), correlation=correlation,
        max_ratio=max_ratio, mean_perp_fraction=mean_perp_fraction,
        suspected=suspected,
    )


def _dominant_eigenvector(mat: List[List[float]], iterations: int = 100) -> Optional[tuple]:
    """Vecteur propre de plus grande valeur propre d'une matrice 3x3
    symetrique semi-definie positive (covariance), par iteration de la
    puissance - suffisant ici (pas besoin de la valeur propre elle-meme,
    juste de la direction), evite une dependance numpy."""
    v = (1.0, 1.0, 1.0)
    for _ in range(iterations):
        nv = (
            mat[0][0] * v[0] + mat[0][1] * v[1] + mat[0][2] * v[2],
            mat[1][0] * v[0] + mat[1][1] * v[1] + mat[1][2] * v[2],
            mat[2][0] * v[0] + mat[2][1] * v[1] + mat[2][2] * v[2],
        )
        norm = math.sqrt(sum(c * c for c in nv))
        if norm == 0.0:
            return None
        v = (nv[0] / norm, nv[1] / norm, nv[2] / norm)
    return v


# ---------------------------------------------------------------------------
# convzmoins : inverse y/z et recode en 'R' les mesures cod1='Z' cod2='-' -
# calcul.f:3929-3951 (citee integralement, algorithme trivial)
# ---------------------------------------------------------------------------

def convert_z_minus(ech: SelectedSample) -> int:
    """Retourne le nombre de mesures converties."""
    count = 0
    for m in ech.mesures:
        if m.cod1 == "Z" and m.cod2 == "-":
            m.y = -m.y
            m.z = -m.z
            m.cod2 = "R"
            count += 1
    return count


# ---------------------------------------------------------------------------
# exportthellier : export au format .tdt (ThellierTool) - fichiers.f:1033-1091
# ---------------------------------------------------------------------------

_THELLIER_STEP_OFFSET = {"R": 0.1, "P": 0.2, "V": 0.5}


def export_thellier_tdt(ech: SelectedSample, out_dir: str) -> str:
    """Port de `exportthellier` (fichiers.f:1033-1091 - version compilee du
    projet, PAS fichiers_mod_magic.f qui n'est pas dans le makefile).
    Applique d'abord `cleanpaleo` sur X/Y/Z/Q/L/F (lignes GRM/qualite non
    pertinentes) SUR UNE COPIE des mesures (ne mute pas `ech`, contrairement
    au Fortran qui mute les donnees globales - l'export ne doit pas avoir
    d'effet de bord sur la selection en memoire). Ecrit `<out_dir>/<id>.tdt`.
    Retourne le chemin ecrit."""
    mesures = [m for m in ech.mesures if m.cod1 not in _CLEANPALEO_CODES]
    path = os.path.join(out_dir, f"{ech.id}.tdt")
    lines = ["Thellier-tdt", f"{(ech.com or '')[:2]}\t0.0\t0.0\t0.0\t0.0"]
    vol = ech.vol or 1.0
    for m in mesures:
        mag, dec, inc = polere(m.x, m.y, m.z)
        rxx = mag * 1.0e3 / vol
        step = m.etape + _THELLIER_STEP_OFFSET.get(m.cod1, 0.0)
        lines.append(f"{ech.id}\t{step:.1f}\t{_fortran_e(rxx, 12, 5)}\t{dec:.2f}\t{inc:.2f}")
    with open(path, "w", encoding="iso-8859-1", errors="replace") as f:
        f.write("\n".join(lines) + "\n")
    return path
