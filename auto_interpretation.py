"""
Interpretation AUTOMATIQUE (suggestion) de diagrammes de desaimantation -
demande explicite utilisateur ("une routine intelligente qui ferait des
interpretations des diagrammes de desaimantation"), deuxieme des deux
routines demandees (la premiere, l'evaluation d'interpretations deja
faites, est dans interpretation_quality.py).

Principe (choisi explicitement plutot qu'un modele IA/ML - discute avec
l'utilisateur : "je recommande une approche heuristique... plus
verifiable, coherent avec l'esprit du projet") : la sequence de
desaimantation est decoupee en segments rectilignes par une segmentation
OPTIMALE (programmation dynamique, voir _segment_sequence), les points de
transition entre deux composantes pouvant rester hors de tout segment -
comme dans une interpretation manuelle. Les ajustements eux-memes passent
EXCLUSIVEMENT par calcul.linear_fit/planar_fit (aucune nouvelle formule
statistique).

1) Composante PRINCIPALE (primary) : le segment le plus haut en
   temperature/champ. Ancree a l'origine si l'ajustement ancre explique
   mieux les points que le libre (decide dans l'optimisation) ; des
   points de fin peuvent etre ecartes (magnetisation parasite en toute fin
   de sequence - demande explicite utilisateur : "often, the last steps
   at high temperature is a spurious magnetization that departs from
   zero... it is best to discard the last points and keep the vector
   anchored to the origin"). Plan (grand cercle) si la droite est
   mediocre (site 14NQ04).
2) Composante SECONDAIRE (secondary) : le segment precedent, toujours
   libre (une surimpression secondaire ne vise pas necessairement
   l'origine).

Remplace (2026-11) l'extension gloutonne d'origine (fenetre agrandie point
par point tant que le MAD ne sautait pas de plus de 5 deg), qui fusionnait
deux composantes successives des que la courbure entre elles etait
progressive. Compare aux interpretations manuelles retenues dans les
moyennes de site (premier choix a moins de 10 deg) : Tibet 85.7% -> 96.0%,
Beagle 67% -> 80%, Tierra del Fuego 88% -> 95%, Peninsule Antarctique 98%
-> 98%, CARF 76% -> 94%. Quelle composante entre dans la moyenne d'un site
(caracteristique haute temperature ou reaimantation) reste une decision
geologique : l'utilisateur choisit p/s/b a chaque specimen.

Les memes points que le Zijderveld affiche (zijderveld.draw_zijderveld) -
runs IRM/thermique-d'IRM exclus (selection.split_experiments/
experiment_kind) - pour que la composante suggeree corresponde a ce que
l'utilisateur voit reellement dans ce diagramme.

Notation reutilisee de interpretation_quality.py (grade_mad, seuils
MAD/fraction NRM) pour une echelle de confiance COHERENTE entre "evaluer
une interpretation existante" et "noter une suggestion automatique".

LIMITE ASSUMEE (avertissement donne des la premiere discussion, a
rappeler a l'utilisateur) : fonctionne bien sur une decroissance propre a
une seule composante par segment ; reste faillible sur des composantes
qui se chevauchent (transition progressive plutot que franche). Une
suggestion est TOUJOURS a valider par l'utilisateur, jamais archivee
automatiquement - `propose_components` ne fait AUCUNE ecriture dans
self.results/.pmagres, c'est a l'appelant (app.py) de proposer
l'archivage via les mecanismes existants (fit_line + archivres) si
l'utilisateur accepte une suggestion.
"""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from calcul import (
    AUTO_OUTLIER_MIN_DEG, FitResult, build_site_mean_result, combine_antipodal_groups,
    fisher_combine_lines_planes, fisher_from_results, fisher_mean, fit_line, fit_plane, linear_fit,
    planar_fit, site_lat_lon_from_donnees, split_outliers, _correct_dec_inc, _next_specimen_c,
    _unit_vector,
)
from interpretation_quality import (
    grade_mad, _nrm_fraction_grade, _GRADE_ORDER, _LINE_MAD_THRESHOLD,
    _PLANE_MAD_THRESHOLD, _vector_diff,
)
from selection import apply_orientation, zijderveld_measurements, polere

MIN_POINTS = 4          # points minimum d'un segment libre
MIN_POINTS_ANCHORED = 3  # un segment ancre a l'origine : l'origine sert de 4e contrainte
_TRIM_MAD_RATIO = 0.5    # voir propose_components : debut de la principale releve
_TRIM_MIN_POINTS = 5
MAX_SEGMENTS = 3
# Couts de la segmentation (voir _segment_sequence), en deg^2 : un point
# non attribue coute autant qu'un point a MAD=10 deg, un segment de plus
# autant que ~17 deg^2 par point sur 18 points. Regles sur les
# interpretations manuelles de Tibet, verifies sans re-reglage sur Beagle,
# Tierra del Fuego, Peninsule Antarctique et CARF : toute la plage
# 50-100 / 150-600 donne le meme resultat a ~1% pres.
GAP_COST = 100.0
SEGMENT_COST = 300.0
# La composante principale ne commence jamais en dessous de ces paliers
# (demande explicite utilisateur : "pour les composantes principales,
# mettre un seuil d'etapes a ne pas prendre en dessous de 10mT ou 130 degC
# en thermique") - les premiers paliers portent les composantes visqueuses.
PRIMARY_MIN_AF = 10.0        # mT
PRIMARY_MIN_THERMAL = 130.0  # degC
_AF_CODES = ("A", "F")
# composante haute temperature courte en fin de sequence (voir _terminal_component)
TERMINAL_DROP = 0.5          # le dernier point garde moins de la moitie de l'aimantation restante
TERMINAL_MIN_ANGLE = 20.0    # deg entre la composante terminale et le segment libre
TERMINAL_MAX_POINTS = 5
_THERMAL_CODES = ("D", "S", "T", "K", "R", "V", "P")


@dataclass
class ComponentSuggestion:
    label: str            # "primary" / "secondary"
    step_first: float
    step_last: float
    anchored: bool
    dec: float
    inc: float
    mad: float
    nb: int
    nrm_fraction: Optional[float]
    grade: str
    notes: List[str] = field(default_factory=list)
    kind: str = "line"    # "line" ou "plane" (grand cercle) - voir propose_components
    # indices (0-based) du segment dans ech.mesures - les mesures du Zijderveld
    # en sont le debut ; retrouver le segment par la valeur du palier echouait
    # sur les paliers mesures deux fois
    index_first: int = 0
    index_last: int = 0
    # lettre de composante imposee par l'analyse du site (voir
    # _site_uniform_fit) ; None : A pour primary, B pour secondary
    component: Optional[str] = None
    # composante haute temperature portee par les derniers paliers seuls
    # (voir _terminal_component) - exclue du controle de recouvrement du site
    terminal: bool = False
    # raison d'exclure ce specimen de la moyenne de site (voir
    # _site_great_circles : aimantation remplacee par la foudre)
    site_outlier: Optional[str] = None


_zijderveld_measurements = zijderveld_measurements  # voir selection.py (partage avec interpretation_quality.py)


def _nrm_fraction(mesures: list, start: int, end: int) -> Optional[float]:
    """Meme calcul VECTORIEL (VDS, voir interpretation_quality._nrm_fraction
    et _vector_diff) que le module d'evaluation - notation reutilisee pour
    une echelle de confiance COHERENTE entre suggestion automatique et
    interpretation deja archivee (voir docstring module). Reimplemente ici
    (plutot que d'appeler directement interpretation_quality._nrm_fraction)
    sur la liste de mesures deja filtree (indices 0-based [start,end]
    inclus) plutot que sur ech.mesures/des indices 1-based - evite
    d'exiger que la selection Zijderveld soit un prefixe strict de
    ech.mesures (vrai aujourd'hui mais pas une garantie a figer ici)."""
    if not mesures or len(mesures) < 2:
        return None
    vds_total = sum(_vector_diff(mesures[i], mesures[i + 1]) for i in range(len(mesures) - 1))
    m_last = mesures[-1]
    vds_total += math.sqrt(m_last.x ** 2 + m_last.y ** 2 + m_last.z ** 2)
    if vds_total == 0:
        return None
    interval_vd = sum(_vector_diff(mesures[i], mesures[i + 1]) for i in range(start, end))
    return interval_vd / vds_total


def _mad(pts: np.ndarray, anchored: bool) -> float:
    """MAD (Kirschvink 1980) sans les tests de calcul.linear_fit - appele
    pour CHAQUE intervalle possible de la sequence (voir
    _segment_sequence), d'ou cette version minimale."""
    c = pts if anchored else pts - pts.mean(axis=0)
    l3, l2, l1 = np.linalg.eigvalsh(c.T @ c)
    if l1 <= 0:
        return 90.0
    return math.degrees(math.atan(math.sqrt(abs((l2 + l3) / l1))))


def _treatment_breaks(mesures: list) -> set:
    """Indices k ou la sequence change de traitement entre k et k+1 (code
    different, ou retour a un palier plus bas : 14NQ1708B, thermique
    140-290 degC puis AF 10-120 mT). Un ajustement ne les franchit pas :
    un resultat .pmagres n'a qu'un code de desaimantation."""
    return {k for k in range(len(mesures) - 1)
            if mesures[k + 1].cod1 != mesures[k].cod1 or mesures[k + 1].etape < mesures[k].etape}


def _segment_sequence(points_xyz: list, first: int, min_last_start: Optional[int] = None,
                      max_prev_end: Optional[int] = None,
                      breaks: Optional[set] = None) -> Tuple[List[Tuple[int, int]], bool]:
    """Decoupage OPTIMAL de points_xyz[first:] en 1 a MAX_SEGMENTS
    segments rectilignes (indices inclus, dans l'ordre des traitements),
    les points restants etant laisses hors de toute composante - demande
    explicite utilisateur ("ameliorer la routine de calcul automatique des
    directions et en faire une comparaison avec les calculs manuels du
    fichier Tibet").

    Remplace l'ancienne extension gloutonne (une fenetre agrandie point par
    point tant que le MAD ne sautait pas de plus de 5 deg) : sur les 833
    specimens de Tibet, quand deux composantes se succedent avec une
    courbure progressive, le MAD montait doucement sans jamais sauter, et
    la fenetre avalait les deux composantes (ex. 15XL2301A : 140-500 ancre,
    entre les deux composantes manuelles 140-260 et 400-500). Ici on
    minimise sur TOUS les decoupages possibles (programmation dynamique) :
        somme des n*MAD^2 des segments + SEGMENT_COST par segment
        + GAP_COST par point non attribue
    ce qui reproduit la facon de faire manuelle : des segments droits, et
    la zone de transition entre deux composantes laissee de cote.

    Le DERNIER segment (le plus haut en traitement) peut etre ancre a
    l'origine ; son ancrage est decide DANS l'optimisation (MAD ancre
    contre MAD libre), et il peut n'avoir que MIN_POINTS_ANCHORED points
    (composante caracteristique courte, ex. 15NQ5714A 630-675 degC, n=3
    en manuel). Des points de fin non attribues restent possibles
    (magnetisation parasite en toute fin de sequence, voir l'ancienne
    docstring de module).

    `min_last_start` (indice) : le dernier segment ne peut pas commencer
    avant ce point - contrainte venue du SITE (voir _site_overlap_window :
    zone ou les spectres de deblocage des deux composantes se recouvrent).
    `max_prev_end` (indice) : les segments qui precedent le dernier
    s'arretent au plus tard a ce point (meme origine).
    Ignorees si elles ne laissent aucun segment possible.

    Retourne (segments, dernier_segment_ancre) ; ([], False) si moins de
    MIN_POINTS_ANCHORED points."""
    pts = np.asarray(points_xyz, dtype=float)
    n = len(pts)
    steps = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    norms = np.linalg.norm(pts, axis=1)
    cost_free, cost_anc = {}, {}
    breaks = breaks or set()
    for i in range(first, n):
        for j in range(i + 1, n):
            if j - 1 in breaks:
                break  # `breaks` : voir _treatment_breaks
            k = j - i + 1
            path = float(steps[i:j].sum())
            # meme test de linearite que calcul.linear_fit (rapport < 0.5
            # refuse) : un segment que linear_fit rejetterait ensuite ne
            # doit pas etre retenu ici.
            if k >= MIN_POINTS and np.linalg.norm(pts[i] - pts[j]) >= 0.5 * path:
                cost_free[(i, j)] = k * _mad(pts[i:j + 1], False) ** 2
            if k >= MIN_POINTS_ANCHORED and norms[i] >= 0.5 * (path + norms[j]):
                cost_anc[(i, j)] = k * _mad(pts[i:j + 1], True) ** 2
    inf = float("inf")
    # best[k][s] : cout minimal des points first..k-1 avec s segments LIBRES
    # (le dernier segment, eventuellement ancre, est ajoute ensuite).
    best = [[inf] * MAX_SEGMENTS for _ in range(n + 1)]
    back = [[None] * MAX_SEGMENTS for _ in range(n + 1)]
    best[first][0] = 0.0
    for k in range(first, n):
        for s in range(MAX_SEGMENTS):
            b = best[k][s]
            if b == inf:
                continue
            if b + GAP_COST < best[k + 1][s]:
                best[k + 1][s], back[k + 1][s] = b + GAP_COST, (k, s, None)
            if s + 1 < MAX_SEGMENTS:
                for j in range(k + MIN_POINTS - 1, n):
                    if max_prev_end is not None and j > max_prev_end:
                        break
                    if (k, j) not in cost_free:
                        continue
                    c = b + cost_free[(k, j)] + SEGMENT_COST
                    if c < best[j + 1][s + 1]:
                        best[j + 1][s + 1], back[j + 1][s + 1] = c, (k, s, (k, j))

    choice = None
    for anchored, costs in ((False, cost_free), (True, cost_anc)):
        for (i, j), c in costs.items():
            if min_last_start is not None and i < min_last_start:
                continue
            for s in range(MAX_SEGMENTS):
                if best[i][s] == inf:
                    continue
                total = best[i][s] + c + SEGMENT_COST + GAP_COST * (n - 1 - j)
                if choice is None or total < choice[0]:
                    choice = (total, s, i, j, anchored)
    if choice is None:
        if min_last_start is not None or max_prev_end is not None:
            return _segment_sequence(points_xyz, first, breaks=breaks)
        return [], False
    _, s, i, j, anchored = choice
    segments = [(i, j)]
    k = i
    while k > first:
        k, s, seg = back[k][s]
        if seg:
            segments.append(seg)
    return segments[::-1], anchored


# palier aberrant isole - demande explicite utilisateur (14NQ1708B : "en
# fait l'etape a 350 est une erreur de mesure. Test facile a faire")
_SPIKE_NEIGHBOUR_RATIO = 0.2   # voisins entre eux / palier a ses voisins
_SPIKE_MIN_JUMP = 0.5          # saut minimal, en fraction du moment des voisins
_SPIKE_MIN_NRM_FRACTION = 0.1  # voisins d'au moins 10 % de la NRM (pas la fin bruitee)
_SPIKE_CODES = ("D", "F")      # paliers de desaimantation seulement (pas Thellier V/R/Z...)


def spurious_steps(ech) -> List[Tuple[object, str]]:
    """Paliers isoles aberrants de la sequence Zijderveld de `ech` : le
    vecteur du palier i est loin de ses DEUX voisins alors que ces voisins
    sont proches l'un de l'autre - |v(i+1) - v(i-1)| < _SPIKE_NEIGHBOUR_RATIO
    x min(|v(i) - v(i-1)|, |v(i+1) - v(i)|), le plus petit de ces sauts
    depassant _SPIKE_MIN_JUMP x le moment des voisins. Une desaimantation
    normale a l'inverse : deux paliers d'ecart enlevent plus qu'un seul.
    14NQ1708B : 350 degC, moment x1.8 et 45 deg de 290 degC, puis l'AF
    10 mT reprend la direction et l'intensite de 290 degC. Le premier et le
    dernier palier ne sont pas testables, ni un palier 0 ; le palier et ses
    deux voisins doivent etre des paliers de desaimantation (_SPIKE_CODES),
    et les voisins garder au moins _SPIKE_MIN_NRM_FRACTION de la NRM.
    Sur les 5 jeux de donnees de reference : 40 paliers sur ~5600
    specimens, dont 36 hors de tout ajustement manuel. Retourne
    [(mesure, raison)]."""
    mesures = _zijderveld_measurements(ech)
    out = []
    if len(mesures) < 3:
        return out
    nrm = max(np.linalg.norm([m.x, m.y, m.z]) for m in mesures[:2])
    for i in range(1, len(mesures) - 1):
        if mesures[i].etape == 0 or any(m.cod1 not in _SPIKE_CODES for m in mesures[i - 1:i + 2]):
            continue
        p, c, n = (np.array([m.x, m.y, m.z]) for m in mesures[i - 1:i + 2])
        a, b, ab = np.linalg.norm(c - p), np.linalg.norm(n - c), np.linalg.norm(n - p)
        scale = max(np.linalg.norm(p), np.linalg.norm(n))
        if (scale <= 0 or min(np.linalg.norm(p), np.linalg.norm(n)) < _SPIKE_MIN_NRM_FRACTION * nrm
                or min(a, b) < _SPIKE_MIN_JUMP * scale or ab >= _SPIKE_NEIGHBOUR_RATIO * min(a, b)):
            continue
        cosang = float(c @ (p + n)) / (np.linalg.norm(c) * np.linalg.norm(p + n) or 1.0)
        ang = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
        out.append((mesures[i], f"moment x{np.linalg.norm(c) / (0.5 * (np.linalg.norm(p) + np.linalg.norm(n))):.1f} "
                                f"and {ang:.0f} deg from its neighbours, which agree with each other"))
    return out


def remove_spurious_steps(samples: list) -> List[Tuple[str, str, str]]:
    """Retire de ech.mesures (en memoire, comme le menu Remove step) les
    paliers de spurious_steps. Retourne [(specimen, palier, raison)]."""
    removed = []
    for ech in samples:
        for m, why in spurious_steps(ech):
            ech.mesures.remove(m)
            removed.append((ech.id.strip(), f"{m.cod1} {m.etape:g}", why))
    return removed


def _first_demag_index(mesures: list) -> int:
    """Premier palier de traitement : saute le NRM brut ET ses eventuelles
    repetitions (plusieurs mesures a l'etape 0 - ex. 14NQ0701A, "0N 0N",
    ou l'ancienne regle min_start=1 gardait le second NRM dans le fit)."""
    i = 0
    while i < len(mesures) - 1 and mesures[i].etape == 0:
        i += 1
    return i


_SUSCEPTIBILITY_JUMP_RATIO = 3.0  # voir _susceptibility_jump_note


def _susceptibility_jump_note(mesures: list, start: int, end: int) -> Optional[str]:
    """Detecte un saut brutal de susceptibilite (colonne `s`, mesuree a
    chaque palier - voir le rappel de cette pratique dans le guide
    utilisateur, section "A way of working") dans la fenetre [start,end]
    (en incluant le point juste AVANT `start`, pour detecter un saut qui
    demarre pile au debut de la fenetre) - signe frequent d'alteration
    thermique (creation d'une nouvelle phase magnetique pendant le
    chauffage) qui peut fausser une direction/un plan sans que le MAD du
    fit lui-meme ne le revele (un ajustement sur peu de points, surtout
    un PLAN a 4 points, absorbe
    facilement le bruit d'une alteration sans que ca degrade son MAD) -
    demande explicite utilisateur ("voici des exemples de calcul de
    composantes secondaire (site 14NQ03) : est-ce possible d'ameliorer
    l'option auto-interpret pour se rapprocher des interpretations
    faites manuellement").

    Valeurs `s` NULLES ignorees (convention deja observee sur donnees
    reelles Tibet : `s` n'est mesure qu'une etape sur deux, 0.0 sur les
    etapes intermediaires - PAS une susceptibilite reellement nulle).

    NOTE SEULEMENT - ne modifie NI le grade NI la decision de recherche.
    Volontairement PAS un seuil dur qui exclurait/degraderait
    automatiquement : verifie sur donnees reelles (meme fichier Tibet)
    qu'un saut du meme ordre de grandeur (x4) apparait aussi DANS une
    fenetre que l'interpretation manuelle garde deliberement telle
    quelle (14NQ0104A, 530-630 degC, ancree, MAD=0.8 - un test qui
    rejetait/degradait automatiquement sur ce seul critere aurait
    degrade cette interpretation, pourtant bonne). Un vrai saut
    d'alteration (14NQ0301A, 575-650 degC : x2.9 puis x3.5 puis x7.3
    puis x3.5, jusqu'a x194 par rapport a la valeur de depart) reste
    bien plus extreme que ce faux-positif possible - le seuil x3 attrape
    les deux, d'ou une simple note plutot qu'un rejet."""
    prev_s = None
    for m in mesures[max(start - 1, 0):end + 1]:
        if not m.s:
            continue
        if prev_s and m.s / prev_s >= _SUSCEPTIBILITY_JUMP_RATIO:
            return (
                f"susceptibility jumps x{m.s / prev_s:.1f} within/around this interval "
                f"(step {m.etape:.0f}) - possible thermal alteration, worth checking manually"
            )
        prev_s = m.s
    return None


def _make_suggestion(
    label: str, mesures: list, start: int, end: int, anchored: bool, fit: dict,
    n_trimmed: int = 0, kind: str = "line", n_trimmed_leading: int = 0,
) -> ComponentSuggestion:
    """`fit` est le dict retourne par calcul.linear_fit (cle "direction")
    pour kind="line", ou calcul.planar_fit (cle "pole") pour kind="plane" -
    demande explicite utilisateur ("check the logic implemented in the
    manual interpretation of site 14NQ04") : dec/inc d'un plan sont ceux
    du POLE du grand cercle (meme convention que calcul.fit_plane), PAS
    une direction caracteristique directement comparable a celle d'une
    ligne - a afficher/interpreter comme tel."""
    dx, dy, dz = fit["pole"] if kind == "plane" else fit["direction"]
    _, dec, inc = polere(dx, dy, dz)
    nrm_fraction = _nrm_fraction(mesures, start, end)
    threshold = _PLANE_MAD_THRESHOLD if kind == "plane" else _LINE_MAD_THRESHOLD
    grade = grade_mad(fit["mad"], threshold)
    notes = []
    if kind == "plane":
        notes.append("plane fit (great circle) - dec/inc shown are the POLE, not a characteristic direction")
    if n_trimmed:
        notes.append(
            f"discarded {n_trimmed} trailing high-treatment point(s) that departed from "
            f"the origin-ward trend, to keep the fit anchored"
        )
    if n_trimmed_leading:
        notes.append(
            f"discarded {n_trimmed_leading} leading low-treatment point(s) with noisy/viscous "
            f"behavior that broke the fit, to find a stable trend"
        )
    frac_grade = _nrm_fraction_grade(nrm_fraction)
    if frac_grade is not None and _GRADE_ORDER[frac_grade] < _GRADE_ORDER.get(grade, 99):
        notes.append(f"only {nrm_fraction * 100:.0f}% of initial NRM described by this interval")
        grade = frac_grade
    # `<= MIN_POINTS` : signale explicitement un ajustement au tout minimum de points acceptes,
    # notamment un PLAN a exactement 4 points (un plan a peu de degres
    # de liberte - 3 points le definissent deja presque completement -
    # une "bonne" MAD a ce nombre minimal est donc un test faible, voir
    # _susceptibility_jump_note pour un exemple reel concret).
    if fit["nb"] <= MIN_POINTS:
        notes.append(f"only {fit['nb']} points - statistically fragile, especially for a {kind} fit")
    susc_note = _susceptibility_jump_note(mesures, start, end)
    if susc_note:
        notes.append(susc_note)
    return ComponentSuggestion(
        label=label, step_first=mesures[start].etape, step_last=mesures[end].etape,
        anchored=anchored, dec=dec, inc=inc, mad=fit["mad"], nb=fit["nb"],
        nrm_fraction=nrm_fraction, grade=grade, notes=notes, kind=kind,
        index_first=start, index_last=end,
    )


def _primary_min_index(mesures: list, first: int) -> int:
    """Premier point ou la composante principale peut commencer : 10 mT en
    AF, 130 degC en thermique (PRIMARY_MIN_AF/PRIMARY_MIN_THERMAL)."""
    for k in range(first, len(mesures)):
        m = mesures[k]
        if m.cod1 in _AF_CODES and m.etape < PRIMARY_MIN_AF:
            continue
        if m.cod1 in _THERMAL_CODES and m.etape < PRIMARY_MIN_THERMAL:
            continue
        return k
    return first


def _af_without_grm_control(mesures: list, first: int) -> bool:
    """Desaimantation AF sans controle de la GRM : ni triplets X/Y/Z
    (reduits en '=' par datatools.reduce_grm_before_fit), ni alternance des
    sequences d'axes '+'/'-' (Rennes, instrument C1/C2) - typiquement une
    seule sequence d'axes (cod2 '0', jeux MagIC de Mongolie/Oslo).

    Demande explicite utilisateur ("voir aussi mon commentaire sur
    Gerritsen et al pour les analyses, en particulier en AF sans controle
    de GRM") : si la GRM n'est ni detectee ni corrigee, la ChRM est
    biaisee des qu'elle n'est pas ancree a l'origine (site YD10, Turquie) ;
    la composante principale est alors ancree des que l'ajustement ancre du
    meme segment reste bon (MAD <= 10 deg). Un ancrage OBLIGATOIRE a ete
    essaye et abandonne : sur la Mongolie (MagIC 19428), il remplacait
    parfois un ajustement libre 10-120 mT (MAD 1.2) par une queue de 3
    points dans une autre direction (TT 31B)."""
    af = [m for m in mesures[first:] if m.cod1 in _AF_CODES]
    if len(af) < len(mesures[first:]) / 2:
        return False
    codes = {m.cod2 for m in af}
    return "=" not in codes and not {"+", "-"} <= codes


def _terminal_component(points_xyz: list, segment: Tuple[int, int]) -> Optional[Tuple[int, dict]]:
    """Composante haute temperature portee par les derniers paliers seuls -
    demande explicite utilisateur ("maintenant le programme n'evalue plus
    la composante haute temperature qui etait bien evaluee auparavant sur
    certains sites comme le NQ04").

    A 14NQ04, l'aimantation restante en fin de segment libre (la
    reaimantation, 140-665 degC) est importante (5-110% du NRM), d'une
    tout autre direction, et disparait d'un coup aux 2-3 derniers paliers
    (14NQ0407B : 110% a 665 degC, 5% a 675 degC). La segmentation optimale
    prefere prolonger le segment libre et laisser le dernier palier hors
    segment plutot que de payer un segment ancre de 3 points ; l'ancienne
    extension gloutonne trouvait cette composante (650-680 ancre, comme
    les interpretations manuelles).

    Cherche, parmi les derniers points du segment libre, le debut k d'un
    segment ancre k..fin (au moins MIN_POINTS_ANCHORED points) tel que :
    le dernier point garde moins de TERMINAL_DROP de l'aimantation au
    point k (vraie chute finale, pas une aimantation parasite en fin de
    sequence), la direction ancree s'ecarte de plus de TERMINAL_MIN_ANGLE
    du segment libre, et l'ajustement ancre est bon (MAD <= 10 deg). Garde
    le k le plus bas qui reste bon. Retourne (k, ajustement) ou None."""
    pts = np.asarray(points_xyz, dtype=float)
    n = len(pts)
    i, j = segment
    free = linear_fit(points_xyz[i:j + 1], anchored=False, mad_threshold=180.0)
    if free is None:
        return None
    # directions ORIENTEES : sens de l'aimantation retiree le long du segment
    # libre, sens du vecteur restant pour la composante terminale - une
    # composante terminale de polarite inverse est presque antiparallele a la
    # reaimantation (14NQ0401A, 14NQ0403B) et un angle entre droites (non
    # oriente) la confondait avec elle
    removed = np.asarray(free["direction"], dtype=float)
    if removed @ (pts[i] - pts[j]) < 0:
        removed = -removed
    best = None
    for k in range(max(i + MIN_POINTS - 1, n - TERMINAL_MAX_POINTS), n - MIN_POINTS_ANCHORED + 1):
        if k > j or np.linalg.norm(pts[-1]) >= TERMINAL_DROP * np.linalg.norm(pts[k]):
            continue
        fit = linear_fit(points_xyz[k:], anchored=True, mad_threshold=_LINE_MAD_THRESHOLD)
        if fit is None or grade_mad(fit["mad"], _LINE_MAD_THRESHOLD) not in ("excellent", "good"):
            continue
        remaining = np.asarray(fit["direction"], dtype=float)
        if remaining @ pts[k] < 0:
            remaining = -remaining
        cosang = max(-1.0, min(1.0, float(remaining @ removed)))
        if math.degrees(math.acos(cosang)) <= TERMINAL_MIN_ANGLE:
            continue
        best = (k, fit)
        break
    return best


def propose_components(ech, clean_from: Optional[float] = None,
                       low_until: Optional[float] = None) -> List[ComponentSuggestion]:
    """Propose 0, 1 ou 2 composantes (primary, et secondary si la
    segmentation trouve un segment avant elle) pour ce specimen. Ne modifie rien,
    n'ecrit rien - voir docstring module.

    Le NRM BRUT (mesure(s) a l'etape 0 en debut de sequence, non
    traitees - voir _first_demag_index) est TOUJOURS exclu, INCONDITIONNELLEMENT (pas
    seulement si cod1=='N' - demande explicite utilisateur : "by default,
    it is best to not include the NRM0... start the interpretation at
    best at the first demag step") - convention paleomagnetique standard,
    PAS un critere statistique : verifie sur donnees reelles
    (old_pmag.ren) que le MAD change a peine (0.54 -> 0.45 deg) selon
    qu'on l'inclue ou non, alors que l'interpretation humaine l'exclut
    systematiquement (le NRM brut porte souvent une composante visqueuse
    non representative, retiree des le premier palier de traitement,
    meme s'il "tombe" pres de la droite ajustee par coincidence
    statistique).

    `clean_from` (temperature/champ, fourni par propose_components_for_site) :
    premier palier ou le site ne retire plus que la composante principale ;
    la composante principale commence au plus tot au point qui le precede.
    `low_until` : dernier palier ou le site ne retire encore que la
    composante secondaire ; les segments qui precedent la principale
    s'arretent au plus tard a ce palier."""
    mesures = _zijderveld_measurements(ech)
    if len(mesures) < MIN_POINTS:
        return []
    points_xyz = [(m.x, m.y, m.z) for m in mesures]
    first = _first_demag_index(mesures)
    if len(mesures) - first < MIN_POINTS:
        first = max(0, len(mesures) - MIN_POINTS)

    min_last_start = None
    if clean_from is not None:
        before = [k for k in range(first, len(mesures)) if mesures[k].etape < clean_from]
        if before and before[-1] > first:
            min_last_start = before[-1]
    max_prev_end = None
    if low_until is not None:
        upto = [k for k in range(first, len(mesures)) if mesures[k].etape <= low_until]
        if upto:
            max_prev_end = upto[-1]
    lowest = _primary_min_index(mesures, first)
    af_no_grm = _af_without_grm_control(mesures, first)
    segments, last_anchored = _segment_sequence(points_xyz, first, min_last_start, max_prev_end,
                                                _treatment_breaks(mesures))
    terminal = None
    if segments:
        terminal = _terminal_component(points_xyz, segments[-1])
        if terminal is not None:
            k = terminal[0]
            segments = segments[:-1] + [(segments[-1][0], k), (k, len(points_xyz) - 1)]
            last_anchored = True
    suggestions: List[ComponentSuggestion] = []
    # les deux segments les plus hauts en traitement : la composante
    # caracteristique ("primary") et celle qui la precede ("secondary") ;
    # un 3e segment eventuel, plus bas, est le plus souvent visqueux.
    for rank, (i, j) in enumerate(reversed(segments[-2:])):
        label = "primary" if rank == 0 else "secondary"
        anchored = last_anchored if rank == 0 else False
        # paliers sous 10 mT / 130 degC retires de la composante principale
        # APRES la segmentation : imposer le seuil a la segmentation elle-meme
        # la poussait a redecouper une composante unique (14NQ0604B, 210-555
        # degC en manuel, devenait 100-460 + une fausse composante 500-665)
        starts = [i]
        if rank == 0 and i < lowest and j - lowest + 1 >= (MIN_POINTS_ANCHORED if anchored else MIN_POINTS):
            starts.insert(0, lowest)
        fit = None
        for i in starts:  # sans les paliers bas, sinon (test de linearite refuse) le segment entier
            fit = linear_fit(points_xyz[i:j + 1], anchored=anchored, mad_threshold=180.0)
            if fit is None:  # test de linearite de linear_fit
                fit = linear_fit(points_xyz[i:j + 1], anchored=not anchored, mad_threshold=180.0)
                if fit is not None:
                    anchored = not anchored
            if fit is not None:
                break
        if fit is None:
            continue
        n_low = i - starts[-1]
        if rank == 0 and af_no_grm and not anchored:
            anc = linear_fit(points_xyz[i:j + 1], anchored=True, mad_threshold=180.0)
            if anc is not None and grade_mad(anc["mad"], _LINE_MAD_THRESHOLD) in ("excellent", "good"):
                fit, anchored = anc, True
        # debut de la composante principale releve quand le MAD y est
        # nettement meilleur - demande explicite utilisateur ("pour la
        # composante principale si le MAD est meilleur en prenant une
        # etape superieure; faire ce choix", 15NQ5612A : 400-675 MAD 3.8,
        # choix manuel 550-675 MAD 1.5) : premier depart dont le MAD
        # tombe sous _TRIM_MAD_RATIO fois celui du segment, en gardant au
        # moins _TRIM_MIN_POINTS paliers.
        n_raised = 0
        if rank == 0 and _TRIM_MAD_RATIO:
            for k in range(i + 1, j - _TRIM_MIN_POINTS + 2):
                t = linear_fit(points_xyz[k:j + 1], anchored=anchored, mad_threshold=180.0)
                if t is not None and t["mad"] <= _TRIM_MAD_RATIO * fit["mad"]:
                    n_raised, mad_before, i, fit = k - i, fit["mad"], k, t
                    break
        # Plan (grand cercle) si la droite de la composante principale est
        # mediocre (MAD > 10 deg, "marginal") - demande explicite
        # utilisateur (site 14NQ04, Tibet), meme regle qu'avant la
        # segmentation : une droite "good"/"excellent" reste toujours
        # preferee. Etendre au "marginal" de Evaluate interpretations
        # (droite ancree a plus de 15 deg de la droite libre) a ete teste
        # ("lorsque l'interpretation par une ligne donne un calcul
        # marginal: prendre un plan si celui ci est bon") : +181 plans
        # sur les 5 jeux de reference, moyennes un peu moins proches des
        # manuelles (13CB35 8 -> 17 deg, P29 21 -> 31) ; meme en
        # n'exigeant un plan "excellent" que si la droite libre est aussi
        # mauvaise : +67 plans, 3 moyennes meilleures (15NQ55 5.9 -> 3.6),
        # 2 moins bonnes. Non retenu.
        if rank == 0 and grade_mad(fit["mad"], _LINE_MAD_THRESHOLD) not in ("excellent", "good"):
            plane = planar_fit(points_xyz[i:j + 1], normalize=True)
            if plane is not None and grade_mad(plane["mad"], _PLANE_MAD_THRESHOLD) in ("excellent", "good"):
                suggestions.append(_make_suggestion(label, mesures, i, j, True, plane, kind="plane"))
                continue
        n_trailing = len(mesures) - 1 - j if (rank == 0 and anchored) else 0
        suggestions.append(_make_suggestion(label, mesures, i, j, anchored, fit, n_trimmed=n_trailing))
        if n_raised:
            suggestions[-1].notes.append(
                f"start raised by {n_raised} step(s): MAD {fit['mad']:.1f} instead of "
                f"{mad_before:.1f} from the segment start")
        if n_low:
            suggestions[-1].notes.append(
                f"{n_low} step(s) below {PRIMARY_MIN_AF:g} mT / {PRIMARY_MIN_THERMAL:g} degC left out")
    if terminal is not None and suggestions and suggestions[0].label == "primary" \
            and suggestions[0].index_first == terminal[0]:
        suggestions[0].terminal = True
        suggestions[0].notes.append(
            "high-temperature component carried by the last steps only (sharp final drop, "
            "direction different from the component before it)")
    if af_no_grm and suggestions and suggestions[0].label == "primary":
        top = suggestions[0]
        if top.kind == "plane":
            top.notes.append("AF without GRM control: this great circle may come from GRM "
                             "rather than from a lightning overprint")
        elif top.anchored:
            top.notes.append("AF without GRM control: anchored to the origin, an unanchored "
                             "fit would be biased by GRM")
    return suggestions


def format_suggestions(ech_id: str, suggestions: List[ComponentSuggestion]) -> str:
    """Rendu compact - une ligne courte par suggestion (demande explicite
    utilisateur : "the line is too long and difficult to read"), dec/inc
    omis ici (visibles directement sur le Zijderveld superpose)."""
    if not suggestions:
        return f"{ech_id}: no component could be automatically proposed (too few usable points, or too noisy)."
    lines = [f"{ech_id} - auto-interpretation suggestions:"]
    for s in suggestions:
        nrm_txt = f"  NRM={s.nrm_fraction * 100:.0f}%" if s.nrm_fraction is not None else ""
        anc = "plane" if s.kind == "plane" else ("anc" if s.anchored else "free")
        lines.append(
            f"  {s.label:<9s} {s.step_first:.0f}-{s.step_last:.0f}  {anc}  n={s.nb}"
            f"  MAD={s.mad:.1f}{nrm_txt}  [{s.grade}]"
        )
        for note in s.notes:
            lines.append(f"      -> {note}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Consensus au niveau du SITE - demande explicite utilisateur ("you miss the
# point at a site level. if some samples have a well defined secondary
# magnetization in a wide temperature range (as site NQ03) then it is
# likely that this behavior is the same for all samples") suite a
# l'exemple reel 14NQ03 (Tibet) : sur 17 specimens, 17 partagent en
# realite le MEME palier bien defini (~100/140 a 530-650 degC) - mais il
# n'est en tete ("primary") que pour les specimens ou il a pu s'ancrer a
# l'origine ; pour les 8 autres, un petit ajustement (4 points, souvent un
# plan) pris dans la queue alteree du chauffage (voir
# _susceptibility_jump_note) usurpe la place de "primary" alors que le
# MEME palier large existe bien chez eux aussi, releque en "secondary".
# Une suggestion par specimen ISOLE ne peut pas le savoir - propose_
# components_for_site regarde TOUS les specimens du site ensemble pour
# le reperer.
#
# RESULTAT FINAL - purement consultatif (voir docstring de
# propose_components_for_site pour l'historique complet) : deux
# contre-exemples reels rencontres en verifiant sur TOUT le fichier
# Tibet (pas seulement 14NQ03) ont ecarte toute idee d'ECHANGER
# automatiquement les etiquettes primary/secondary - (1) site 14NQ04 :
# le palier large libre y est au contraire la composante SECONDAIRE (le
# vrai "primary" utilise pour la moyenne de site y est le palier ANCRE
# haute temperature - l'exact INVERSE de 14NQ03) ; (2) specimen
# 14NQ0403B (meme site) : meme apres avoir corrige (1), un simple
# recouvrement de FENETRE ne garantit pas la bonne DIRECTION - aurait
# remplace une interpretation deja quasi correcte (marginale, mais
# juste) par une franchement fausse. propose_components_for_site se
# limite donc a AJOUTER UNE NOTE des deux cotes en cas de desaccord entre
# le "primary" d'un specimen et le consensus du site, sans jamais
# reordonner ni relabelliser - la decision reste entierement a
# l'utilisateur, qui voit les deux options et le Zijderveld.
# ---------------------------------------------------------------------------

def _is_well_graded(s: ComponentSuggestion) -> bool:
    """Une suggestion sans reserve majeure - ni fragile (nombre de points
    minimal), ni suspecte d'alteration (voir _susceptibility_jump_note) -
    utilisable comme reference pour degager un consensus de site. Le
    grade seul ("good"/"excellent") ne suffit pas : c'est precisement le
    cas du plan a 4 points sur la queue alteree qui motive cette
    fonction - il peut etre grade "good" tout en portant deja la note
    d'alteration."""
    if s.grade not in ("excellent", "good"):
        return False
    if any("statistically fragile" in n or "possible thermal alteration" in n for n in s.notes):
        return False
    return True


def _site_consensus_window(per_specimen: Dict[str, List[ComponentSuggestion]]) -> Optional[Tuple[float, float]]:
    """Fenetre de traitement consensuelle du site : mediane des step_first/
    step_last, construite UNIQUEMENT a partir des suggestions DEJA
    etiquetees "primary" et bien notees - jamais un melange primary+
    secondary.

    Correction critique suite a un contre-exemple reel (site 14NQ04,
    Tibet - demande explicite utilisateur : "still not good. Check the
    results in .pmagres for site NQ04") : a ce site, comp=1 (le vrai
    "primary", celui qui entre dans la moyenne de site, voir .pmagres) EST
    la fenetre haute temperature ANCREE (500-685 degC, ligne ou plan) ;
    comp=2 (secondaire, une surimpression) EST la fenetre large basse/
    moyenne temperature LIBRE (180-650 degC) pour LES 11 specimens - soit
    l'INVERSE exact de 14NQ03, ou c'etait la fenetre large qui etait le
    "primary". Une premiere version de cette fonction melangeait les
    suggestions primary ET secondary de tout le site dans UNE seule
    mediane ; a 14NQ04, la recherche "secondary" (libre, large) reussit
    plus souvent (bien notee sur 9 specimens sur 11) que la recherche
    "primary" (ancree, souvent fragile/4 points, bien notee sur seulement
    7) - le melange faisait donc pencher le consensus vers la fenetre
    large et PROMOUVAIT A TORT la surimpression au rang de "primary" pour
    les specimens dont le vrai "primary" ancre etait fragile (402A, 403B,
    405B, 408B), inversant l'interpretation. Se limiter aux seules
    suggestions DEJA "primary" evite ce retournement : le role
    geologiquement significatif (primary vs secondary) varie d'un site a
    l'autre et ne peut pas se deviner par la seule geometrie/frequence de
    reussite - seul ce qui est deja convenu comme "primary" ailleurs sur
    le site sert de reference pour repecher un specimen dont LE SIEN a
    echoue.

    None si moins de 3 suggestions "primary" bien notees (pas assez pour
    degager un consensus fiable plutot que de deviner)."""
    candidates = [
        s for suggestions in per_specimen.values() for s in suggestions
        if s.label == "primary" and _is_well_graded(s)
    ]
    if len(candidates) < 3:
        return None

    def median(vals: List[float]) -> float:
        vals = sorted(vals)
        n = len(vals)
        mid = n // 2
        return vals[mid] if n % 2 else (vals[mid - 1] + vals[mid]) / 2.0

    return median([s.step_first for s in candidates]), median([s.step_last for s in candidates])


def _overlap_fraction(s: ComponentSuggestion, window: Tuple[float, float]) -> float:
    """Fraction de la fenetre consensuelle couverte par [s.step_first,
    s.step_last] (1.0 = couvre la fenetre entiere ou plus)."""
    w0, w1 = window
    lo, hi = min(s.step_first, s.step_last), max(s.step_first, s.step_last)
    overlap = max(0.0, min(hi, w1) - max(lo, w0))
    span = max(w1 - w0, 1e-9)
    return overlap / span


_SITE_OVERLAP_THRESHOLD = 0.7  # voir propose_components_for_site


_OVERLAP_CLEAN_DEG = 10.0     # voir _site_overlap_window
_OVERLAP_MIN_SEPARATION = 15.0
_OVERLAP_MIN_SPECIMENS = 3
_OVERLAP_MAX_WIDTH = 0.25


def _axial_mean(vectors: List[np.ndarray]) -> np.ndarray:
    """Axe principal (vecteurs unitaires), oriente comme la majorite."""
    t = sum(np.outer(v, v) for v in vectors)
    w, e = np.linalg.eigh(t)
    axis = e[:, 2]
    if sum(float(v @ axis) for v in vectors) < 0:
        axis = -axis
    return axis


def _site_overlap_window(specimens: list, per_specimen: Dict[str, List[ComponentSuggestion]]
                         ) -> Tuple[Optional[float], Optional[float]]:
    """Zone ou les spectres de deblocage des deux composantes du site se
    recouvrent : (dernier palier ou le site ne retire encore QUE la
    composante secondaire, premier palier a partir duquel le SITE ne retire plus que la composante
    principale - demande explicite utilisateur (site 14XL19, Tibet :
    "autour de 300 degC il y a un melange des temperatures de blocage des
    deux composantes... il faut analyser les composantes non seulement par
    echantillon mais en ayant la vision du site").

    Sur un specimen seul, la zone de recouvrement des spectres de
    deblocage est une courbure progressive que la segmentation peut
    absorber dans la composante principale si le MAD reste faible
    (14XL1912A : 210-600 ancre, MAD 4.3, inclinaison tiree de 6 deg vers
    la composante basse temperature). Au niveau du site elle est
    evidente : en in situ, on mesure pour chaque vecteur difference (ce
    qui est retire entre deux paliers) son ecart a la direction moyenne
    des "primary" du site, vers la direction moyenne des "secondary" ; la
    mediane de cet ecart par palier vaut 20-60 deg dans la zone de
    melange et retombe a +/-5 deg ensuite (14XL19 : encore 13 deg a
    370 degC, 3 deg a 410 degC -> composante principale a partir de 370,
    comme les interpretations manuelles).

    Meme mesure dans l'autre sens pour la composante secondaire (ecart a
    la direction "secondary" du site vers la "primary") - demande
    explicite utilisateur ("OUI" a la proposition de traiter aussi la
    composante basse temperature) : a 14XL19 elle ne doit pas depasser
    290 degC (13 deg vers la composante HT des 330 degC), comme les
    interpretations manuelles.

    (None, None) si les deux bornes ne se degagent pas (pas assez de
    specimens, composantes trop proches, aucun melange visible, ou un
    seul cote visible) ; la
    composante principale commence au plus tot au point qui precede la
    seconde borne, la secondaire s'arrete au plus tard a la premiere."""
    by_id = {ech.id: ech for ech in specimens}
    prim, sec = {}, []
    for ech_id, suggestions in per_specimen.items():
        ech = by_id[ech_id]
        for sg in suggestions:
            if sg.kind != "line" or sg.terminal:
                continue
            d, i = math.radians(sg.dec), math.radians(sg.inc)
            v = np.array(apply_orientation(math.cos(i) * math.cos(d), math.cos(i) * math.sin(d),
                                           math.sin(i), ech, 2))
            if sg.label == "primary":
                prim[ech_id] = (v, sg)
            else:
                sec.append(v)
    if len(prim) < _OVERLAP_MIN_SPECIMENS or len(sec) < _OVERLAP_MIN_SPECIMENS:
        return None, None
    h_site = _axial_mean([v for v, _ in prim.values()])
    l_site = _axial_mean(sec)
    if math.degrees(math.acos(min(1.0, abs(float(h_site @ l_site))))) < _OVERLAP_MIN_SEPARATION:
        return None, None

    # ecart (deg) de chaque vecteur difference a une composante, vers l'autre
    toward_low: Dict[float, List[float]] = {}
    toward_high: Dict[float, List[float]] = {}
    for ech_id, (v, sg) in prim.items():
        ech = by_id[ech_id]
        h = h_site if v @ h_site >= 0 else -h_site
        lp = l_site - (l_site @ h) * h
        lp /= np.linalg.norm(lp)
        hp = h - (h @ l_site) * l_site
        hp /= np.linalg.norm(hp)
        mesures = _zijderveld_measurements(ech)
        for k in range(_first_demag_index(mesures) + 1, len(mesures)):
            if mesures[k].etape > sg.step_last:
                break
            a, b = mesures[k - 1], mesures[k]
            dv = np.array(apply_orientation(a.x - b.x, a.y - b.y, a.z - b.z, ech, 2))
            if dv @ h > 0:
                toward_low.setdefault(b.etape, []).append(math.degrees(math.atan2(dv @ lp, dv @ h)))
            if dv @ l_site > 0:
                toward_high.setdefault(b.etape, []).append(math.degrees(math.atan2(dv @ hp, dv @ l_site)))

    def medians(deviations: Dict[float, List[float]]) -> List[Tuple[float, float]]:
        return sorted((t, float(np.median(v))) for t, v in deviations.items()
                      if len(v) >= _OVERLAP_MIN_SPECIMENS)

    high = medians(toward_low)
    onset = None
    for t, med in reversed(high):
        if med > _OVERLAP_CLEAN_DEG:
            break
        onset = t
    if onset is not None and not any(med > 2 * _OVERLAP_CLEAN_DEG for t, med in high if t < onset):
        onset = None

    low = medians(toward_high)
    lt_end = None
    for t, med in low:
        if med > _OVERLAP_CLEAN_DEG:
            break
        lt_end = t
    if lt_end is not None and not any(med > 2 * _OVERLAP_CLEAN_DEG for t, med in low if t > lt_end):
        lt_end = None
    # un vrai recouvrement se voit des DEUX cotes ; une seule borne trahit
    # plutot un site heterogene (15NQ56 : "primary" dispersees, la borne
    # basse seule y deplacait la composante principale de 15NQ5616A de
    # 630-675 a 140-490 degC)
    if lt_end is None or onset is None or lt_end >= onset:
        return None, None
    # zone de melange large (plus du quart de la gamme de traitement) : seule
    # la borne haute est gardee - 14NQ07 donnait 260-655 degC, toute la
    # reaimantation, que la borne basse coupait a 260 degC ; les vrais
    # recouvrements sont etroits (14XL19 290-410, 15XL23 300-370)
    top = max(t for t, _ in high + low)
    if onset - lt_end > _OVERLAP_MAX_WIDTH * top:
        return None, onset
    return lt_end, onset


_FREE_DANG_DEG = 10.0       # voir _site_uniform_fit
_FREE_MAD_RATIO = 2.0
_FREE_MIN_SPECIMENS = 3
_FREE_SITE_FRACTION = 0.25
_FREE_SAME_COMPONENT_DEG = 20.0


def _site_uniform_fit(specimens: list, per_specimen: Dict[str, List[ComponentSuggestion]]) -> Optional[str]:
    """Composante qui ne va PAS a l'origine sur une partie du site -
    demande explicite utilisateur ("parfois sur un site comme le 14NQ03,
    la composante ne passe pas par l'origine pour certains echantillons et
    la meme composante semble passer par l'origine sur certains
    echantillons. Dans ce cas, tous les echantillons doivent etre analyses
    de la meme facon avec un fit libre et un numero de composante B").

    Critere : sur le segment de chaque "primary" (droite), un specimen "ne
    va pas a l'origine" si l'ajustement ancre s'ecarte de plus de 10 deg du
    libre ET double au moins son MAD ; si au moins 3 specimens ET au moins
    25% du site sont dans ce cas, toutes les "primary" du site deviennent
    libres, composante B. Regle sur les choix manuels des 5 jeux de
    reference : detecte 14NQ03, 14NQ09, 15NQ43, 15NQ52 (Tibet) et J25
    (CARF), ajustes tout en libre a la main, et AUCUN des 141 sites ajustes
    tout en ancre. L'angle seul (> 10 deg sur 25% du site) declenchait sur
    19 de ces sites ancres, surtout en AF (Beagle, Terre de Feu) ou un
    segment court et bruite donne des ajustements libre et ancre differents
    sans que la composante rate l'origine. Les sites libres ou libre et
    ancre donnent la meme direction (14NQ15, 15NQ49) ne sont pas concernes.

    Il doit s'agir de la MEME composante : en in situ, la moyenne des
    specimens qui ratent l'origine et celle des autres sont a moins de
    20 deg (14NQ03 2 deg, 15NQ52 13 deg). Sinon le site porte deux
    composantes - une reaimantation libre et une composante haute
    temperature courte et ancree (15NQ57, 15NQ50, 14NQ06, 14NQ13 : 27-58
    deg) - que les interpretations manuelles traitent chacune a sa facon.

    Sinon, si le site porte une seule composante mais que certaines
    "primary" sont ancrees et d'autres libres, toutes sont ancrees
    (composante A) - demande explicite utilisateur ("pourquoi la logique
    se perd au fil des sites. Elle a ete bien appliquee au site NQ03") :
    sur 14NQ01-28, la regle libre n'agissait que sur 14NQ03/14NQ09 et
    laissait 14NQ07 (13 ancres, 6 libres), 14NQ12, 14NQ24... melanges ;
    tous les echantillons d'un site doivent etre traites de la meme facon.
    Un specimen dont l'ajustement ancre echoue (test de linearite, MAD >
    15 deg) reste libre.

    Modifie per_specimen EN PLACE ; retourne "free", "anchored" ou None."""
    by_id = {ech.id: ech for ech in specimens}

    def in_situ(ech_id: str, direction) -> np.ndarray:
        return np.array(apply_orientation(*direction, by_id[ech_id], 2))

    lines = []  # (id, rate l'origine, ajustement libre, ajustement ancre, direction in situ actuelle)
    for ech_id, suggestions in per_specimen.items():
        if not suggestions or suggestions[0].label != "primary" or suggestions[0].kind != "line":
            continue
        sg = suggestions[0]
        pts = [(m.x, m.y, m.z) for m in by_id[ech_id].mesures[sg.index_first:sg.index_last + 1]]
        free = linear_fit(pts, anchored=False, mad_threshold=180.0)
        anchored = linear_fit(pts, anchored=True, mad_threshold=180.0)
        if free is None:
            continue
        # un ajustement ancre refuse (test de linearite) rate l'origine
        off = anchored is None
        if anchored is not None:
            cosang = min(1.0, abs(float(np.dot(free["direction"], anchored["direction"]))))
            off = (math.degrees(math.acos(cosang)) > _FREE_DANG_DEG
                   and anchored["mad"] >= _FREE_MAD_RATIO * free["mad"])
            if anchored["mad"] > _LINE_MAD_THRESHOLD:
                anchored = None  # inutilisable pour un ajustement ancre
        current = (anchored if sg.anchored else free) or free
        lines.append((ech_id, off, free, anchored, in_situ(ech_id, current["direction"])))
    if len(lines) < _FREE_MIN_SPECIMENS:
        return None

    def same_component(group_a: List[np.ndarray], group_b: List[np.ndarray]) -> bool:
        if not group_a or not group_b:
            return True
        cosang = min(1.0, abs(float(_axial_mean(group_a) @ _axial_mean(group_b))))
        return math.degrees(math.acos(cosang)) <= _FREE_SAME_COMPONENT_DEG

    n_off = sum(1 for line in lines if line[1])
    if n_off >= _FREE_MIN_SPECIMENS and n_off >= _FREE_SITE_FRACTION * len(lines):
        if not same_component([in_situ(e, f["direction"]) for e, off, f, _, _ in lines if off],
                              [v for _, off, _, _, v in lines if not off]):
            return None
        mode, component = "free", "B"
        note = (f"site-level check: this component misses the origin on {n_off} of {len(lines)} "
                f"specimens - fitted free on every specimen of the site, component B")
    else:
        if all(per_specimen[e][0].anchored for e, *_ in lines):
            return None
        if not same_component([v for e, *_, v in lines if per_specimen[e][0].anchored],
                              [v for e, *_, v in lines if not per_specimen[e][0].anchored]):
            return None
        mode, component = "anchored", None
        note = ("site-level check: the component goes to the origin on this site - "
                "anchored on every specimen of the site")

    if mode == "free":
        # une principale qui ne s'ajuste pas en libre (moins de MIN_POINTS
        # paliers : une composante finale ancree a 3 points ; MAD > 15 deg)
        # est retiree et la composante suivante devient la composante du
        # site - demande explicite utilisateur (14NQ03 : "pourquoi des
        # composantes C et une B en haute temperatures ?" - 14NQ0317A
        # 650-690 degC libre a 3 points notee B, sa composante 140-650
        # notee C ; 14NQ0315A 630-685 MAD 15.4 non ecrite, 140-630 notee C)
        dropped = set()
        for ech_id, _, free, _, _ in lines:
            suggestions = per_specimen[ech_id]
            sg = suggestions[0]
            if (sg.index_last - sg.index_first + 1 >= MIN_POINTS and free is not None
                    and free["mad"] <= _LINE_MAD_THRESHOLD) or len(suggestions) < 2:
                continue
            removed = suggestions.pop(0)
            nxt = suggestions[0]
            nxt.label = "primary"
            nxt.notes.append(
                f"site-level check: the high-treatment segment {removed.step_first:g}-"
                f"{removed.step_last:g} ({removed.nb} points) cannot be fitted free - "
                f"left out, this component is the site's component")
            dropped.add(ech_id)
        lines = [line for line in lines if line[0] not in dropped]
    for ech_id, _, free, anchored, _ in lines:
        sg = per_specimen[ech_id][0]
        fit = free if mode == "free" else anchored
        if fit is None or sg.anchored == (mode == "anchored"):
            continue
        mesures = _zijderveld_measurements(by_id[ech_id])
        new = _make_suggestion("primary", mesures, sg.index_first, sg.index_last, mode == "anchored", fit)
        new.notes += [n for n in sg.notes if n.startswith("site-level")]
        per_specimen[ech_id][0] = new
    for suggestions in per_specimen.values():
        for sg in suggestions:
            if sg.label == "primary" and sg.kind == "line":
                sg.notes.append(note)
            if component:
                sg.component = component if sg.label == "primary" else "C"
    return mode


def _site_component_letters(specimens: list, per_specimen: Dict[str, List[ComponentSuggestion]]) -> None:
    """Lettre de composante attribuee par DIRECTION, a l'echelle du site -
    demande explicite utilisateur (14NQ04 : "maintenant le programme
    n'evalue plus la composante haute temperature..."). Sur les specimens a
    deux composantes, la "primary" (A) et la "secondary" (B) donnent les
    directions de reference du site (in situ). Un specimen qui n'a qu'une
    composante (14NQ0402A : la composante haute temperature n'y est pas
    isolee) la voit numerotee B si elle est a moins de
    _FREE_SAME_COMPONENT_DEG de la reference B et plus proche d'elle que de
    la reference A ; sinon elle melangeait la reaimantation a la moyenne A.
    Ne touche pas un site deja numerote par _site_uniform_fit."""
    by_id = {ech.id: ech for ech in specimens}

    def in_situ(ech_id: str, sg: ComponentSuggestion) -> np.ndarray:
        d, i = math.radians(sg.dec), math.radians(sg.inc)
        return np.array(apply_orientation(math.cos(i) * math.cos(d), math.cos(i) * math.sin(d),
                                          math.sin(i), by_id[ech_id], 2))

    pairs = [(e, s) for e, s in per_specimen.items()
             if len(s) == 2 and all(sg.kind == "line" and sg.component is None for sg in s)]
    if len(pairs) < _FREE_MIN_SPECIMENS:
        return
    a_vecs = [in_situ(e, s[0]) for e, s in pairs]
    b_vecs = [in_situ(e, s[1]) for e, s in pairs]
    ref_a = _axial_mean(a_vecs)
    ref_b = _axial_mean(b_vecs)
    # intervalles typiques (medianes) des deux composantes du site
    window_a = (float(np.median([s[0].step_first for _, s in pairs])),
                float(np.median([s[0].step_last for _, s in pairs])))
    window_b = (float(np.median([s[1].step_first for _, s in pairs])),
                float(np.median([s[1].step_last for _, s in pairs])))

    def oriented(ref, vectors):
        n_pos = sum(1 for v in vectors if float(v @ ref) >= 0)
        if 3 * max(n_pos, len(vectors) - n_pos) < 2 * len(vectors):
            return None  # pas de polarite majoritaire nette
        return ref if 2 * n_pos >= len(vectors) else -ref

    ref_a_or, ref_b_or = oriented(ref_a, a_vecs), oriented(ref_b, b_vecs)
    for ech_id, suggestions in per_specimen.items():
        if len(suggestions) != 1 or suggestions[0].kind != "line" or suggestions[0].component:
            continue
        sg = suggestions[0]
        v = in_situ(ech_id, sg)
        to_a = math.degrees(math.acos(min(1.0, abs(float(v @ ref_a)))))
        to_b = math.degrees(math.acos(min(1.0, abs(float(v @ ref_b)))))
        if to_b < to_a and to_b <= _FREE_SAME_COMPONENT_DEG:
            sg.component = "B"
            sg.notes.append(
                f"site-level check: same direction as this site's second component "
                f"({to_b:.0f} deg from it, {to_a:.0f} deg from the first) - numbered B")
            continue
        # regle de l'intervalle - demande explicite utilisateur ("la regle
        # c'est qu'au sein d'un site tous les vecteurs dans le ~ meme
        # intervalle de desaimantation devrait avoir la meme categorie de
        # composante") : 14NQ0404B, droite unique 170-675 degC comme
        # 14NQ0402A (B), mais a 11 deg de l'antipode de A et 25 de B, etait
        # numerotee A. Une droite unique qui recouvre surtout l'intervalle
        # B du site est B, sauf si elle a la polarite des droites A et pas
        # celle des droites B (specimen sans reaimantation).
        span = max(sg.step_last - sg.step_first, 1e-9)

        def overlap(w):
            return max(0.0, min(sg.step_last, w[1]) - max(sg.step_first, w[0])) / span

        if overlap(window_b) <= overlap(window_a):
            continue
        if (ref_a_or is not None and ref_b_or is not None and float(ref_a_or @ ref_b_or) < 0
                and not (float(v @ ref_b_or) > 0 and float(v @ ref_a_or) < 0)):
            continue
        sg.component = "B"
        sg.notes.append(
            f"site-level check: same demagnetization interval as this site's second component "
            f"({window_b[0]:g}-{window_b[1]:g}; first component {window_a[0]:g}-{window_a[1]:g}) and "
            f"same polarity - numbered B ({to_b:.0f} deg from it, {to_a:.0f} deg from the first)")


_GC_MAX_LINE_K = 10.0       # voir _site_great_circles
_GC_A95_RATIO = 0.5
_GC_MAX_A95 = 12.0
_GC_PLANE_MAX_MAD = 10.0
_GC_MIN_PLANES = 4
_GC_LINE_KEEP_DEG = 10.0
_LIGHTNING_NRM_RATIO = 10.0
_GC_MAX_TWO_COMP_FRACTION = 0.5


def _sequence_end(mesures: list, first: int) -> int:
    """Dernier palier de la premiere serie de traitement croissante (15NQ3105B
    recommence a 0 mT apres 80 mT)."""
    j = first
    while j + 1 < len(mesures) and mesures[j + 1].etape >= mesures[j].etape:
        j += 1
    return j


def _best_plane(pts: list) -> Tuple[Optional[dict], int, int]:
    """Grand cercle d'un specimen pour _site_great_circles : l'intervalle
    le plus long (au moins MIN_POINTS paliers) de MAD <= _GC_PLANE_MAX_MAD,
    a egalite le MAD le plus faible. Toute la desaimantation ne convient
    pas toujours - demande explicite utilisateur (15NQ35 : "the scatter of
    magnetization tend to be reduced in the last steps") : sur 15NQ35 les
    grands cercles de toute la serie ont un MAD de 16 a 23 deg, les
    interpretations manuelles prennent 0-80, 0-70, 10-60, 15-90 mT."""
    best, best_key, best_range = None, None, (0, 0)
    for i in range(len(pts)):
        for j in range(i + MIN_POINTS - 1, len(pts)):
            fit = planar_fit(pts[i:j + 1], normalize=True)
            if fit is None or fit["mad"] > _GC_PLANE_MAX_MAD:
                continue
            key = (j - i, -fit["mad"])
            if best_key is None or key > best_key:
                best, best_key, best_range = fit, key, (i, j)
    return best, best_range[0], best_range[1]


def _site_great_circles(specimens: list, per_specimen: Dict[str, List[ComponentSuggestion]]) -> bool:
    """Site analyse par intersection de grands cercles - demande explicite
    utilisateur (15NQ31 : "NRM de forte intensite, dispersees en
    directions. Au cours de la desaimantation AF les directions vont vers
    une direction commune, sauf le plus fortement aimante qui a vu son
    aimantation naturelle remplacee par l'effet de foudre [...] ici
    l'intersection de grands cercles est la base").

    Chaque specimen recoit un grand cercle sur le plus long intervalle de
    sa premiere serie de traitement (NRM compris) de MAD <= 10 deg (voir
    _best_plane). Le site passe en grands cercles si les
    droites "primary" sont dispersees (k < _GC_MAX_LINE_K, polarites
    regroupees) et si l'intersection des grands cercles (au moins 4,
    McFadden & McElhinny) est au moins deux fois plus precise
    (a95) que la meilleure composante en droites (A ou B : a 15NQ52 les
    grands cercles ne font que retrouver la reaimantation B, deja bien
    groupee) et sous _GC_MAX_A95. Regle calee sur les 5 jeux de reference :
    declenche sur 15NQ31, 15NQ52, 15NQ43, 14NQ09, P04, P38, P39, M07...
    (moyenne de site 5 a 50 deg plus proche des moyennes manuelles), aucun
    site eloigne ; les sites a deux composantes ou les grands cercles
    melangent les deux (15NQ55, 15NQ57, 15NQ54, J05) ne passent pas le
    critere a95.

    Un specimen garde sa droite si elle est a moins de _GC_LINE_KEEP_DEG
    de l'intersection (point stable atteint) ou si son grand cercle est
    mal defini ; sinon son grand cercle devient la composante A.

    Foudre : un specimen dont le NRM depasse _LIGHTNING_NRM_RATIO fois la
    mediane du site (15NQ3103B, 26 fois) et dont la droite n'est pas
    proche de l'intersection (85 deg) est exclu de la moyenne
    (site_outlier) : son aimantation naturelle a ete remplacee et l'AF ne
    la fait pas converger. Il n'entre pas dans le calcul de
    l'intersection.

    Modifie per_specimen EN PLACE ; retourne True si le site est passe en
    grands cercles."""
    by_id = {ech.id: ech for ech in specimens}

    def in_situ(ech, v) -> np.ndarray:
        return np.array(apply_orientation(*v, ech, 2))

    def unit(dec, inc) -> np.ndarray:
        d, i = math.radians(dec), math.radians(inc)
        return np.array([math.cos(i) * math.cos(d), math.cos(i) * math.sin(d), math.sin(i)])

    def dec_inc(v) -> Tuple[float, float]:
        _, dec, inc = polere(*v)
        return dec, inc

    rows = []  # (ech, mesures, intensite NRM, droite in situ, grand cercle, pole in situ, fin)
    secondaries = []
    for ech in specimens:
        mesures = _zijderveld_measurements(ech)
        suggestions = per_specimen.get(ech.id) or []
        if len(mesures) < MIN_POINTS or not suggestions:
            continue
        last = _sequence_end(mesures, _first_demag_index(mesures))
        pts = [(m.x, m.y, m.z) for m in mesures[:last + 1]]
        plane, start, end = _best_plane(pts)
        primary = suggestions[0]
        line = (in_situ(ech, unit(primary.dec, primary.inc))
                if primary.label == "primary" and primary.kind == "line" else None)
        pole = in_situ(ech, plane["pole"]) if plane is not None else None
        rows.append((ech, mesures, float(np.linalg.norm(pts[0])), line, plane, pole, (start, end)))
        secondaries += [dec_inc(in_situ(ech, unit(sg.dec, sg.inc))) for sg in suggestions[1:]
                        if sg.label == "secondary" and sg.kind == "line"]
    lines = [dec_inc(r[3]) for r in rows if r[3] is not None]
    if len(rows) < _GC_MIN_PLANES or len(lines) < 3:
        return False
    # site a deux composantes : les grands cercles contiennent tous la
    # composante basse temperature, leur intersection seule retombe sur
    # elle (voir _site_mixed_ht) - 14NQ13 : 12 specimens sur 17 ; P46
    # (3 sur 10) garde l'intersection
    two_comp = sum(1 for ech in specimens
                   if any(sg.label == "secondary" and sg.kind == "line"
                          for sg in per_specimen.get(ech.id) or []))
    if two_comp > _GC_MAX_TWO_COMP_FRACTION * len(rows):
        return False
    median_nrm = float(np.median([r[2] for r in rows]))
    lightning = {r[0].id for r in rows if median_nrm > 0 and r[2] >= _LIGHTNING_NRM_RATIO * median_nrm}
    poles = [dec_inc(r[5]) for r in rows
             if r[4] is not None and r[4]["mad"] <= _GC_PLANE_MAX_MAD and r[0].id not in lightning]
    if len(poles) < _GC_MIN_PLANES:
        return False
    line_stats = fisher_mean(combine_antipodal_groups(lines)[0])
    best_line_a95 = line_stats.a95
    if len(secondaries) >= 3:
        best_line_a95 = min(best_line_a95, fisher_mean(combine_antipodal_groups(secondaries)[0]).a95)
    gc = fisher_combine_lines_planes([], poles)
    if not (line_stats.k < _GC_MAX_LINE_K and gc.a95 <= _GC_A95_RATIO * best_line_a95
            and gc.a95 <= _GC_MAX_A95):
        return False

    target = unit(gc.dec, gc.inc)
    note = (f"site-level check: directions scattered (k={line_stats.k:.1f}) but great circles "
            f"converge (intersection {gc.dec:.0f}/{gc.inc:.0f} in situ, a95={gc.a95:.1f}) - "
            f"site analysed by great-circle intersection")
    for ech, mesures, nrm, line, plane, pole, (start, end) in rows:
        suggestions = per_specimen[ech.id]
        primary = suggestions[0]
        line_angle = (math.degrees(math.acos(max(-1.0, min(1.0, float(line @ target)))))
                      if line is not None else None)
        if ech.id in lightning and (line_angle is None or line_angle > _GC_LINE_KEEP_DEG):
            primary.site_outlier = (f"probable lightning: NRM {nrm / median_nrm:.0f} times the site "
                                    f"median, direction not on the great-circle convergence")
            primary.notes.append(f"site-level check: {primary.site_outlier} - left out of the site mean")
        elif (line_angle is None or line_angle > _GC_LINE_KEEP_DEG) and plane is not None \
                and plane["mad"] <= _GC_PLANE_MAX_MAD:
            primary = _make_suggestion("primary", mesures, start, end, True, plane, kind="plane")
            suggestions[0] = primary
            primary.notes.append(note)
        else:
            primary.notes.append(note + " - this specimen keeps its line (stable end point)")
        primary.component = "A"
        for sg in suggestions[1:]:
            sg.component = "B"
    return True


_MIX_MIN_LINES = 3          # voir _site_mixed_ht
_MIX_MAX_PLANES_PER_LINE = 2  # au moins 1/3 de droites, au plus 2/3 de plans
_MIX_LINE_OFF_DEG = 20.0    # droite A a plus de 20 deg de la reference A : HT non isolee
_MIX_PLANE_NEAR_DEG = 25.0  # le grand cercle doit passer a moins de 25 deg de la reference A
_MIX_PLANE_BACK = 99         # le plan peut couvrir toute la composante B
_MIX_MIN_AB_DEG = 30.0      # A et B distinctes (angle axial)
_MIX_TERMINAL_DEG = 30.0    # droite ancree de fin de desaimantation (voir _site_mixed_ht)


def _ht_plane(mesures: list, start: int, last: int, earliest: int) -> Tuple[Optional[dict], int, int]:
    """Grand cercle HT de _site_mixed_ht : il doit aller jusqu'au bout de la
    desaimantation (dernier ou avant-dernier palier) et commencer dans la
    fin de la composante basse temperature (entre `earliest` et `start`),
    comme les plans manuels de Tibet (420-665 degC a 14NQ13) ; le plus
    long de MAD <= _GC_PLANE_MAX_MAD, a egalite le MAD le plus faible."""
    best, best_key, best_range = None, None, (0, 0)
    for j in (last, last - 1):
        for i in range(earliest, start + 1):
            if j - i + 1 < MIN_POINTS:
                continue
            fit = planar_fit([(m.x, m.y, m.z) for m in mesures[i:j + 1]], normalize=True)
            if fit is None or fit["mad"] > _GC_PLANE_MAX_MAD:
                continue
            key = (j - i, -fit["mad"])
            if best_key is None or key > best_key:
                best, best_key, best_range = fit, key, (i, j)
    return best, best_range[0], best_range[1]


def _site_mixed_ht(specimens: list, per_specimen: Dict[str, List[ComponentSuggestion]]) -> int:
    """Composante haute temperature d'un site a deux composantes : droites
    ET grands cercles combines (McFadden & McElhinny) - demande explicite
    utilisateur ("les grands cercles sont aussi a prendre avec precaution
    car sauf pour la foudre, ils sont aussi controles par la composante
    basse temperature" ; "l'ideal, c'est environ au moins 1/3 lignes et
    2/3 plans maximum pour forcer l'intersection sur une direction qui
    n'est pas antipodale a celle en basse temperature").

    Les grands cercles d'un site a deux composantes contiennent tous la
    direction basse temperature (B) : leur intersection seule retombe sur
    B ou sur son antipode. Ce sont les droites HT isolees (references A
    des specimens a deux composantes) qui fixent la direction ; un grand
    cercle n'y apporte que son point le plus proche, comme dans les
    choix manuels de Tibet (14NQ04, 14NQ13, 15NQ50 : plans 240-500 degC
    jusqu'a la fin, combines a autant ou plus de droites ancrees
    610-690 degC).

    Un specimen recoit un grand cercle A quand sa composante HT n'est pas
    isolee : il n'a qu'une droite B, ou sa droite A est a plus de
    _MIX_LINE_OFF_DEG de la reference A (droite melangee). Le plan va de
    la fin de sa composante B (ou du debut de sa droite melangee) a la
    fin de la premiere serie de traitement (_best_plane, MAD <= 10 deg)
    et doit passer a moins de _MIX_PLANE_NEAR_DEG de la reference A. Au
    plus _MIX_MAX_PLANES_PER_LINE plans par droite A retenue (les MAD
    les plus faibles). Ne s'applique pas a un site deja analyse en
    grands cercles (_site_great_circles) ni si A et B ne sont pas
    distinctes. Modifie per_specimen EN PLACE ; retourne le nombre de
    plans ajoutes."""
    by_id = {ech.id: ech for ech in specimens}

    def in_situ(ech, v) -> np.ndarray:
        return np.array(apply_orientation(*v, ech, 2))

    def unit(dec, inc) -> np.ndarray:
        d, i = math.radians(dec), math.radians(inc)
        return np.array([math.cos(i) * math.cos(d), math.cos(i) * math.sin(d), math.sin(i)])

    def axial(u, v) -> float:
        return math.degrees(math.acos(min(1.0, abs(float(u @ v)))))

    def letter(sg: ComponentSuggestion) -> str:
        return sg.component or ("A" if sg.label == "primary" else "B")

    if any(sg.kind == "plane" for s in per_specimen.values() for sg in s):
        return 0
    pairs = [(e, s) for e, s in per_specimen.items()
             if e in by_id and len(s) == 2 and letter(s[0]) == "A" and letter(s[1]) == "B"]
    if len(pairs) < _FREE_MIN_SPECIMENS:
        return 0
    ref_a = _axial_mean([in_situ(by_id[e], unit(s[0].dec, s[0].inc)) for e, s in pairs])
    ref_b = _axial_mean([in_situ(by_id[e], unit(s[1].dec, s[1].inc)) for e, s in pairs])
    a_vecs = [in_situ(by_id[e], unit(s[0].dec, s[0].inc)) for e, s in pairs]
    b_vecs = [in_situ(by_id[e], unit(s[1].dec, s[1].inc)) for e, s in pairs]

    def oriented(ref, vectors):
        """ref orientee vers la polarite majoritaire ; None sans polarite
        nette (moins de 2/3 des directions)."""
        n_pos = sum(1 for v in vectors if float(v @ ref) >= 0)
        if 3 * max(n_pos, len(vectors) - n_pos) < 2 * len(vectors):
            return None
        return ref if 2 * n_pos >= len(vectors) else -ref

    ref_a_or, ref_b_or = oriented(ref_a, a_vecs), oriented(ref_b, b_vecs)
    # A et B de polarites opposees sont distinctes meme presque antipodales
    # (14NQ06 : A ~165/-15, B ~345/35, 30 deg axiaux mais 150 orientes)
    opposite = ref_a_or is not None and ref_b_or is not None and float(ref_a_or @ ref_b_or) < 0
    if axial(ref_a, ref_b) < _MIX_MIN_AB_DEG and not opposite:
        return 0

    def terminal_line(ech, mesures, start, last):
        """Droite ANCREE de fin de desaimantation, apres la composante B,
        comme les choix manuels de 14NQ06 (3 a 5 paliers 635-677 degC,
        MAD jusqu'a 13 deg) : la plus longue de MAD <= _LINE_MAD_THRESHOLD
        a moins de _MIX_TERMINAL_DEG de la reference A (orientee si la
        polarite du site est nette)."""
        pts = [(m.x, m.y, m.z) for m in mesures]
        for i in range(max(start, 0), last - MIN_POINTS_ANCHORED + 2):
            fit = linear_fit(pts[i:last + 1], anchored=True, mad_threshold=180.0)
            if fit is None or fit["mad"] > _LINE_MAD_THRESHOLD:
                continue
            v = in_situ(ech, np.array(fit["direction"]))
            ang = (math.degrees(math.acos(max(-1.0, min(1.0, float(v @ ref_a_or)))))
                   if ref_a_or is not None else axial(v, ref_a))
            if ang <= _MIX_TERMINAL_DEG:
                return i, fit, ang
        return None

    def thermal(ech) -> bool:
        return any(m.cod1 in _THERMAL_CODES for m in _zijderveld_measurements(ech))

    n_lines, candidates, line_thermal, terminal = 0, [], [], []
    for ech in specimens:
        suggestions = per_specimen.get(ech.id) or []
        if not suggestions or any(sg.site_outlier for sg in suggestions):
            continue
        a_lines = [sg for sg in suggestions if letter(sg) == "A" and sg.kind == "line"]
        b_lines = [sg for sg in suggestions if letter(sg) == "B" and sg.kind == "line"]
        if a_lines and axial(in_situ(ech, unit(a_lines[0].dec, a_lines[0].inc)), ref_a) <= _MIX_LINE_OFF_DEG:
            n_lines += 1
            line_thermal.append(thermal(ech))
            continue
        mesures = _zijderveld_measurements(ech)
        if len(mesures) < MIN_POINTS:
            continue
        last = _sequence_end(mesures, _first_demag_index(mesures))
        if not a_lines and b_lines:
            found = terminal_line(ech, mesures, b_lines[0].index_last, last)
            if found is not None:
                terminal.append((ech, mesures, found))
                n_lines += 1
                line_thermal.append(thermal(ech))
                continue
        if a_lines:
            start = a_lines[0].index_first
        elif b_lines:
            start = b_lines[0].index_last
        else:
            continue
        plane, i, j = _ht_plane(mesures, start, last,
                                max(b_lines[0].index_first, start - _MIX_PLANE_BACK) if b_lines else start)
        if plane is None:
            continue
        near = 90.0 - axial(in_situ(ech, plane["pole"]), ref_a)
        if near > _MIX_PLANE_NEAR_DEG:
            continue
        candidates.append((plane["mad"], ech, mesures, i, j, plane, near))
    if n_lines < _MIX_MIN_LINES or not (candidates or terminal):
        return 0
    for ech, mesures, (i, fit, ang) in terminal:
        suggestions = per_specimen[ech.id]
        last = _sequence_end(mesures, _first_demag_index(mesures))
        sg = _make_suggestion("primary", mesures, i, last, True, fit)
        sg.component = "A"
        sg.notes.append(
            f"site-level check: high-temperature component of the site found in the last steps "
            f"({sg.step_first:g}-{sg.step_last:g}, anchored, {ang:.0f} deg from the site's HT lines)")
        for s_ in suggestions:
            s_.component = s_.component or "B"
            if s_.label == "primary":
                s_.label = "secondary"
        suggestions.insert(0, sg)
    # meme traitement que la majorite des droites HT (15XL2309B, seul
    # specimen AF d'un site thermique)
    site_thermal = 2 * sum(line_thermal) >= len(line_thermal)
    candidates = [c for c in candidates if thermal(c[1]) == site_thermal]
    candidates.sort(key=lambda c: c[0])
    kept = candidates[:_MIX_MAX_PLANES_PER_LINE * n_lines]
    for mad, ech, mesures, i, j, plane, near in kept:
        suggestions = per_specimen[ech.id]
        sg = _make_suggestion("primary", mesures, i, j, True, plane, kind="plane")
        sg.component = "A"
        sg.notes.append(
            f"site-level check: high-temperature component not isolated on this specimen - "
            f"great circle ({sg.step_first:g}-{sg.step_last:g}, passes {near:.0f} deg from the "
            f"site's HT lines) combined with the {n_lines} HT lines of the site")
        for s in [s for s in suggestions if letter(s) == "A" and s.kind == "line"]:
            suggestions.remove(s)
        for s in suggestions:
            s.component = s.component or "B"
            if s.label == "primary":
                s.label = "secondary"
        suggestions.insert(0, sg)
    return len(kept)


def propose_components_for_site(specimens: list) -> Dict[str, List[ComponentSuggestion]]:
    """Comme propose_components, specimen par specimen, PUIS ajoute une
    note (SEULEMENT une note - voir plus bas pourquoi) quand la
    suggestion en tete ("primary") d'un specimen s'ecarte de la fenetre
    de traitement typique des "primary" bien notes des AUTRES specimens
    du meme site (recouvrement < 70%, voir _overlap_fraction) alors
    qu'une autre suggestion BIEN NOTEE (_is_well_graded) du meme specimen
    la couvre a 70% ou plus - demande explicite utilisateur ("you miss
    the point at a site level...", voir le commentaire de section
    ci-dessus pour l'exemple reel 14NQ03 qui a motive cette fonction).

    N'ECHANGE PLUS automatiquement les etiquettes primary/secondary (une
    premiere version le faisait). Abandonne apres un DEUXIEME
    contre-exemple reel, plus grave que le premier (voir
    _site_consensus_window pour le premier, 14NQ04 vs 14NQ03) : meme en
    comparant desormais primary-contre-primary uniquement, specimen
    14NQ0403B (site 14NQ04) montre qu'un simple recouvrement de FENETRE
    (temperature/champ) ne garantit PAS que la DIRECTION du candidat soit
    la bonne - son "primary" existant (650-680 degC, ancre, MAD=2.4,
    dec=71.3/inc=-61.4) correspond en realite presque exactement a
    l'interpretation manuelle (650-680, dec=71.0/inc=-62.1) mais n'est
    note que "marginal" a cause de sa fragilite statistique (4 points,
    comme le fit manuel n'en utilise que 3) ; son "secondary" (fenetre
    large 120-625, libre) recouvre la fenetre consensuelle a 71% - juste
    au-dessus du seuil - mais pointe dans une direction totalement
    differente (dec=236/inc=76, une toute autre composante). Un echange
    automatique aurait donc REMPLACE une interpretation deja quasi
    correcte par une franchement fausse. Comparer des directions (et pas
    seulement des fenetres) demanderait de combiner correctement lignes
    ET plans au niveau du site - probleme explicitement identifie comme
    hors de portee ici ("la combinaison de
    plusieurs plans/lignes de PLUSIEURS specimens ... laisse pour un
    travail separe"). Plutot que d'echanger a l'aveugle sur un seul
    critere insuffisant, cette fonction se contente donc de SIGNALER le
    desaccord et de laisser le choix a l'utilisateur, qui voit les DEUX
    suggestions, leurs dec/inc, et le Zijderveld superpose - meme posture
    que _susceptibility_jump_note (note seule, jamais une decision
    automatique).

    Ne recalcule et ne modifie AUCUN chiffre ni aucune etiquette
    (dec/inc/MAD/n/label inchanges, toujours ceux de propose_components,
    meme ordre) - ajoute uniquement des notes explicatives. Si aucune
    fenetre consensuelle ne se degage (moins de 3 "primary" bien notes
    sur tout le site), renvoie exactement ce que propose_components
    aurait donne specimen par specimen, sans y toucher."""
    per_specimen = {ech.id: propose_components(ech) for ech in specimens}
    first_pass = dict(per_specimen)
    # recouvrement des spectres de deblocage vu au niveau du site (voir
    # _site_overlap_window) : chaque specimen est re-segmente pour que la
    # composante secondaire s'arrete avant la zone de melange et que la
    # principale commence apres ; deux iterations suffisent (les
    # directions du site se resserrent).
    limits = (None, None)
    for _ in range(3):
        new_limits = _site_overlap_window(specimens, per_specimen)
        if new_limits == (None, None) or new_limits == limits:
            break
        limits = new_limits
        # un specimen a composante terminale (voir _terminal_component) garde
        # son premier decoupage : le recouvrement de site concerne deux
        # composantes etalees sur une large gamme, pas les derniers paliers
        per_specimen = {ech.id: first_pass[ech.id] if any(sg.terminal for sg in first_pass[ech.id])
                        else propose_components(ech, clean_from=limits[1], low_until=limits[0])
                        for ech in specimens}
    low_until, clean_from = limits
    for suggestions in per_specimen.values():
        for sg in suggestions:
            if sg.label == "primary" and clean_from is not None:
                sg.notes.append(
                    f"site-level check: both components still unblock together below "
                    f"{clean_from:g} - this component starts after that overlap")
            if sg.label == "secondary" and low_until is not None:
                sg.notes.append(
                    f"site-level check: the next component starts unblocking above "
                    f"{low_until:g} - this component stops before that overlap")
    _site_uniform_fit(specimens, per_specimen)
    _site_component_letters(specimens, per_specimen)
    if _site_great_circles(specimens, per_specimen):
        # site en grands cercles : la fenetre de consensus des droites n'a
        # plus de sens (15NQ3105B invitait a revenir a sa droite)
        return per_specimen
    _site_mixed_ht(specimens, per_specimen)
    window = _site_consensus_window(per_specimen)
    if window is None:
        return per_specimen

    for suggestions in per_specimen.values():
        if not suggestions:
            continue
        top = suggestions[0]
        if _overlap_fraction(top, window) >= _SITE_OVERLAP_THRESHOLD:
            continue  # deja coherent avec le site, rien a signaler

        best_alt, best_overlap = None, _SITE_OVERLAP_THRESHOLD
        for s in suggestions[1:]:
            if not _is_well_graded(s):
                continue
            frac = _overlap_fraction(s, window)
            if frac >= best_overlap:
                best_alt, best_overlap = s, frac
        if best_alt is None:
            continue

        top.notes.append(
            f"this site's other well-graded primaries mostly fall in "
            f"{window[0]:.0f}-{window[1]:.0f} - this specimen's own "
            f"'{best_alt.label}' suggestion ({best_alt.step_first:.0f}-{best_alt.step_last:.0f}) "
            f"overlaps that window instead; worth comparing both directions manually"
        )
        best_alt.notes.append(
            f"overlaps this site's consensus primary window "
            f"({window[0]:.0f}-{window[1]:.0f}, built from other specimens) better than "
            f"the standalone top pick does - worth checking manually, but NOT auto-promoted "
            f"(window overlap alone does not guarantee this is the same component)"
        )

    return per_specimen


AUTO_MIN_SPECIMENS = 3  # specimens minimum pour une moyenne de site automatique
# k minimal (in situ) d'une moyenne de site automatique - demande explicite
# utilisateur (15NQ34 : "highly scattered and with k < 5, these means should
# not be reported")
AUTO_MIN_K = 5.0


# droite a plus de _FAR_LINE_DEG de la moyenne de site : grand cercle ou
# rejet - demande explicite utilisateur ("faire aussi un test lorsque la
# ligne est a plus de 45 deg de la moyenne du site, tenter un grand cercle
# ou sinon rejeter le specimen")
_FAR_LINE_DEG = 45.0
_FAR_PLANE_NEAR_DEG = 25.0


def _far_lines_to_planes(results: List[FitResult], ech_of: dict, sg_of: dict, out) -> List[FitResult]:
    """Voir _FAR_LINE_DEG. Moyenne in situ (polarites regroupees) des
    resultats d'un site/composante ; chaque DROITE a plus de _FAR_LINE_DEG
    d'elle est remplacee par le grand cercle le plus long du specimen (MAD
    <= _GC_PLANE_MAX_MAD, au moins MIN_POINTS paliers) qui commence au plus
    tard au debut de la droite, finit au plus tot a ce debut, et passe a
    moins de _FAR_PLANE_NEAR_DEG de la moyenne calculee SANS les droites
    eloignees ; sinon elle est retiree de la moyenne (outlier, le resultat
    specimen reste ecrit). Modifie out.specimens ; retourne la liste pour
    la moyenne."""
    def dev(res, mean_vec):
        v = _unit_vector(*_correct_dec_inc(res, 2))
        ang = math.degrees(math.acos(max(-1.0, min(1.0, float(v @ mean_vec)))))
        return abs(90.0 - ang) if res.cat1 == "P" else min(ang, 180.0 - ang)

    if len(results) < AUTO_MIN_SPECIMENS:
        return results
    stats = fisher_from_results(results, orientation=2, combine_antipodal=True)
    mean_vec = _unit_vector(stats.dec, stats.inc)
    far = [r for r in results if r.cat1 != "P" and dev(r, mean_vec) > _FAR_LINE_DEG]
    if not far:
        return results
    near = [r for r in results if r not in far]
    if len(near) < AUTO_MIN_SPECIMENS:
        return results
    stats = fisher_from_results(near, orientation=2, combine_antipodal=True)
    mean_vec = _unit_vector(stats.dec, stats.inc)
    kept = list(near)
    for r in far:
        ech, sg = ech_of.get(r.id), sg_of.get(id(r))
        best = None
        if ech is not None and sg is not None:
            mesures = _zijderveld_measurements(ech)
            pts = [(m.x, m.y, m.z) for m in mesures]
            first = _first_demag_index(mesures)
            last = _sequence_end(mesures, first)
            for i in range(first, sg.index_first + 1):
                for j in range(max(sg.index_first, i + MIN_POINTS - 1), last + 1):
                    fit = planar_fit(pts[i:j + 1], normalize=True)
                    if fit is None or fit["mad"] > _GC_PLANE_MAX_MAD:
                        continue
                    key = (j - i, -fit["mad"])
                    if best is None or key > best[0]:
                        best = (key, i, j)
        plane = None
        if best is not None:
            _, i, j = best
            plane = fit_plane(ech, i + 1, j + 1, numcomp=r.numcomp)
            if plane is not None and dev(plane, mean_vec) > _FAR_PLANE_NEAR_DEG:
                plane = None
        d = dev(r, mean_vec)
        site = (ech.magic_site or "").strip() if ech is not None else ""
        if plane is None:
            out.outliers.append((site, r.component, r.id, d))
            continue
        plane.component, plane.c = r.component, r.c
        out.specimens[out.specimens.index(r)] = plane
        out.far_to_plane.append((site, r.component, r.id, d, plane.step_first, plane.step_last))
        kept.append(plane)
    return kept


@dataclass
class AutomaticResults:
    specimens: List[FitResult]
    means: List[FitResult]
    not_written: List[Tuple[str, str]]               # (specimen, raison)
    outliers: List[Tuple[str, str, str, float]]      # (site, composante, specimen, ecart)
    excluded: List[Tuple[str, str, str]] = field(default_factory=list)  # (site, specimen, raison)
    scattered: List[Tuple[str, str, int, float]] = field(default_factory=list)  # (site, composante, N, k)
    # (site, composante, specimen, ecart de la droite, palier debut, palier fin du grand cercle)
    far_to_plane: List[Tuple[str, str, str, float, float, float]] = field(default_factory=list)
    spurious: List[Tuple[str, str, str]] = field(default_factory=list)  # (specimen, palier, raison)


def automatic_results(samples: list, donnees, existing_ids: Optional[set] = None) -> AutomaticResults:
    """Mode TOUT AUTOMATIQUE - demande explicite utilisateur ("OUI" a la
    proposition d'ajouter a l'application ce qui produisait le
    _autotest.pmagres de Tibet) : archive les suggestions de
    propose_components_for_site comme des resultats (composante A = 1re
    suggestion, B = 2e), puis calcule pour chaque site et chaque composante
    la moyenne de Fisher in situ ET apres pendage (lignes et plans combines,
    polarites opposees regroupees), sans les valeurs aberrantes
    (split_outliers, plus de 2 ecarts-types angulaires ET plus de
    AUTO_OUTLIER_MIN_DEG - "en mode tout automatique mettre un seuil a
    10 deg" - et specimens en desaccord avec leurs soeurs de la meme
    carotte, magic_sample), avec au moins AUTO_MIN_SPECIMENS specimens.
    Une moyenne de k < AUTO_MIN_K (in situ) n'est pas ecrite (scattered).

    `existing_ids` : identifiants `c` deja presents dans le fichier de
    destination, pour ne pas en reprendre un.

    Les specimens sans site (magic_site vide) sont interpretes seuls et
    n'entrent dans aucune moyenne. N'ecrit rien : voir
    calcul.write_results_file. Les ajustements passent par fit_line/
    fit_plane, comme un ajustement manuel (seuils MAD 15/25 deg compris)."""
    # paliers aberrants isoles retires d'abord (en memoire, comme Remove step)
    spurious = remove_spurious_steps(samples)
    by_site: Dict[str, list] = {}
    suggestions: Dict[str, List[ComponentSuggestion]] = {}
    for ech in samples:
        site = (ech.magic_site or "").strip()
        if site:
            by_site.setdefault(site, []).append(ech)
        else:
            suggestions[ech.id] = propose_components(ech)
    for echs in by_site.values():
        suggestions.update(propose_components_for_site(echs))

    out = AutomaticResults([], [], [], [])
    out.spurious = spurious
    # `c` deja pris dans le .pmagres de destination (mode site par site,
    # qui ecrit dans le .pmagres de l'utilisateur)
    existing: set = set(existing_ids or ())
    per_site: Dict[Tuple[str, str], List[FitResult]] = {}
    sg_of: dict = {}
    ech_of = {ech.id: ech for ech in samples}
    for ech in samples:
        for sg in suggestions.get(ech.id, []):
            component, numcomp = ("A", 1) if sg.label == "primary" else ("B", 2)
            if sg.component:
                component = sg.component
            if sg.kind == "plane":
                res = fit_plane(ech, sg.index_first + 1, sg.index_last + 1, numcomp=numcomp)
            else:
                res = fit_line(ech, sg.index_first + 1, sg.index_last + 1,
                               anchored=sg.anchored, numcomp=numcomp)
            if res is None:
                out.not_written.append((ech.id, f"{sg.label} {sg.kind}: MAD {sg.mad:.1f} above the limit"))
                continue
            res.component = component
            sg_of[id(res)] = sg
            res.c = _next_specimen_c(ech.id, existing)
            existing.add(res.c)
            out.specimens.append(res)
            site = (ech.magic_site or "").strip()
            if site and sg.site_outlier:
                out.excluded.append((site, ech.id, sg.site_outlier))
            elif site:
                per_site.setdefault((site, component), []).append(res)

    sample_of = {ech.id: (ech.magic_sample or "").strip() for ech in samples}
    for (site, component), results in sorted(per_site.items()):
        if len(results) < AUTO_MIN_SPECIMENS:
            continue
        results = _far_lines_to_planes(results, ech_of, sg_of, out)
        kept, rejected, _ = split_outliers(results, orientation=2, combine_antipodal=True,
                                           min_deviation=AUTO_OUTLIER_MIN_DEG, sample_of=sample_of)
        out.outliers += [(site, component, r.id, d) for r, d in rejected]
        in_situ = fisher_from_results(kept, orientation=2, combine_antipodal=True)
        if in_situ.k < AUTO_MIN_K:
            out.scattered.append((site, component, in_situ.n, in_situ.k))
            continue
        lat, lon = site_lat_lon_from_donnees(site, donnees) or (0.0, 0.0)
        for orientation in (2, 3):
            stats = fisher_from_results(kept, orientation=orientation, combine_antipodal=True)
            out.means.append(build_site_mean_result(stats, kept, site, orientation, lat, lon,
                                                    component=component))
    return out
