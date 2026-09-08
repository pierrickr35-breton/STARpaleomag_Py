"""
Interpretation AUTOMATIQUE (suggestion) de diagrammes de desaimantation -
demande explicite utilisateur ("une routine intelligente qui ferait des
interpretations des diagrammes de desaimantation"), deuxieme des deux
routines demandees (la premiere, l'evaluation d'interpretations deja
faites, est dans interpretation_quality.py).

Principe (choisi explicitement plutot qu'un modele IA/ML - discute avec
l'utilisateur : "je recommande une approche heuristique... plus
verifiable, coherent avec l'esprit du projet") : PCA en extension gloutonne
(greedy), reutilisant EXCLUSIVEMENT calcul.linear_fit (aucune nouvelle
formule statistique) :

1) Composante PRINCIPALE (primary) : part des DERNIERS points (haute
   temperature/champ fort - la ou la composante caracteristique est
   generalement isolee) et etend la fenetre VERS L'ARRIERE tant que le
   MAD reste acceptable et n'augmente pas brusquement (signe qu'une autre
   composante a ete atteinte). Ancree a l'origine par defaut (hypothese
   standard pour une composante caracteristique qui doit decroitre vers
   zero). Si l'ancrage echoue DES le dernier point (magnetisation
   parasite en toute fin de sequence - demande explicite utilisateur :
   "often, the last steps at high temperature is a spurious
   magnetization that departs from zero... it is best to discard the
   last points and keep the vector anchored to the origin"), retire
   d'abord 1 a 3 points du cote haute temperature et reessaie ANCRE
   (voir _primary_search) - repli en libre seulement si meme cela
   echoue.
2) Composante SECONDAIRE (secondary) : meme extension gloutonne mais
   depuis le PREMIER point, VERS L'AVANT, limitee aux points non deja
   couverts par la composante principale. Libre par defaut (une
   surimpression secondaire ne vise pas necessairement l'origine), repli
   en ancree si le libre echoue d'emblee. Absente si aucun point ne
   reste disponible.

Convention "ancre = primary" (demande explicite utilisateur : "by
default, if the component is anchored to the origin, use primary") : si
la recherche secondaire (avant) finit ancree alors que la principale
(arriere) ne l'est pas, les etiquettes sont echangees - l'ancrage a
l'origine prime sur le sens de recherche pour decider quelle composante
est "primary".

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

from calcul import FitResult, linear_fit, planar_fit
from interpretation_quality import (
    grade_mad, _nrm_fraction_grade, _GRADE_ORDER, _LINE_MAD_THRESHOLD,
    _PLANE_MAD_THRESHOLD, _vector_diff,
)
from selection import zijderveld_measurements, polere

MIN_POINTS = 4
MAX_MAD_JUMP = 5.0  # deg - au-dela, on considere avoir atteint une autre composante


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
    kind: str = "line"    # "line" ou "plane" (grand cercle) - voir _primary_search


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


def _greedy_backward(points_xyz: list, anchored: bool, min_start: int = 0) -> Optional[dict]:
    """Meilleure fenetre [start, n-1] trouvee en etendant vers l'arriere
    depuis les derniers MIN_POINTS points - voir docstring module, etape
    1. `min_start` : ne descend pas en dessous de cet indice (voir
    propose_components - exclut le NRM brut, index 0). Retourne None si
    meme la fenetre minimale echoue."""
    n = len(points_xyz)
    if n - min_start < MIN_POINTS:
        return None
    best = None
    prev_mad = None
    for start in range(n - MIN_POINTS, min_start - 1, -1):
        fit = linear_fit(points_xyz[start:n], anchored=anchored, mad_threshold=180.0)
        if fit is None:
            break
        mad = fit["mad"]
        if mad > _LINE_MAD_THRESHOLD:
            break
        if prev_mad is not None and mad - prev_mad > MAX_MAD_JUMP:
            break
        best = (start, n - 1, fit)
        prev_mad = mad
    return best


def _greedy_backward_plane(points_xyz: list, min_start: int = 0) -> Optional[tuple]:
    """Variante PLAN de _greedy_backward - demande explicite utilisateur,
    apres verification sur des donnees reelles (site 14NQ04, fichier
    Tibet) que _primary_search (LIGNE seule) manque ou force un mauvais
    ajustement sur des specimens ou l'interpretation manuelle utilisait un
    plan (grand cercle) pour la composante haute temperature/champ fort -
    cas frequent quand cette composante n'est pas isolee (chevauchement
    avec une composante non totalement retiree). `planar_fit` n'a pas de
    notion d'ancrage (toujours "libre" au sens geometrique - un plan
    contraint une direction sans en choisir une seule ; convention
    Fortran/calcul.fit_plane : orig='o' pour un plan, quoi qu'il en soit -
    voir _make_plane_suggestion) et rejette lui-meme un MAD>25 (voir
    calcul.planar_fit) - pas besoin de reverifier le seuil ici comme pour
    une ligne."""
    n = len(points_xyz)
    if n - min_start < MIN_POINTS:
        return None
    best = None
    prev_mad = None
    for start in range(n - MIN_POINTS, min_start - 1, -1):
        fit = planar_fit(points_xyz[start:n], normalize=True)
        if fit is None:
            break
        mad = fit["mad"]
        if prev_mad is not None and mad - prev_mad > MAX_MAD_JUMP:
            break
        best = (start, n - 1, fit)
        prev_mad = mad
    return best


def _greedy_forward(points_xyz: list, anchored: bool, min_start: int, limit: int) -> Optional[dict]:
    """Symetrique de _greedy_backward : etend vers l'avant depuis les
    premiers MIN_POINTS points a partir de `min_start` (voir
    propose_components - exclut le NRM brut, index 0), sans depasser
    l'indice `limit` (exclu - les points deja pris par la composante
    principale)."""
    if limit - min_start < MIN_POINTS:
        return None
    best = None
    prev_mad = None
    for end in range(min_start + MIN_POINTS - 1, limit):
        fit = linear_fit(points_xyz[min_start:end + 1], anchored=anchored, mad_threshold=180.0)
        if fit is None:
            break
        mad = fit["mad"]
        if mad > _LINE_MAD_THRESHOLD:
            break
        if prev_mad is not None and mad - prev_mad > MAX_MAD_JUMP:
            break
        best = (min_start, end, fit)
        prev_mad = mad
    return best


def _secondary_search(points_xyz: list, min_start: int, limit: int, max_start_trim: int = 3):
    """Symetrique de _primary_search, mais sur le bord AVANT (bas
    temperature/champ faible) de la recherche secondaire - trouve en
    reponse au meme signalement utilisateur que _susceptibility_jump_note
    (site 14NQ03, specimen 14NQ0301A) : la toute premiere fenetre libre
    (min_start..min_start+MIN_POINTS-1) peut echouer d'emblee (MAD trop
    haut des le depart) si les tout premiers paliers portent une
    surimpression visqueuse/incoherente non representative - meme raison
    que le repli ancre de _primary_search sur le bord haute temperature,
    en miroir. Verifie sur donnees reelles : sans ce recul, 14NQ0301A
    tombait sur un ajustement ancre 0-555 MAD=14.2 ("marginal") au lieu
    du 230-555 libre MAD=2.1 ("excellent") que retrouve ce recul (2
    points ecartes) - et qui correspond a l'interpretation manuelle
    (230-555, comp=2, libre, MAD=2.1).

    Essaie D'ABORD libre (comportement par defaut d'une composante
    secondaire, cf. module docstring), en reculant le debut de 0 a
    `max_start_trim` points des que la fenetre minimale echoue
    completement ; PUIS seulement ancre si meme ce recul echoue
    (meme ordre libre-puis-ancre que l'appelant avant ce changement).

    Retourne (start, end, fit, anchored, n_trimmed) ou None."""
    for anchored in (False, True):
        for start_trim in range(0, max_start_trim + 1):
            start = min_start + start_trim
            if limit - start < MIN_POINTS:
                break
            candidate = _greedy_forward(points_xyz, anchored=anchored, min_start=start, limit=limit)
            if candidate is not None:
                cstart, cend, fit = candidate
                return cstart, cend, fit, anchored, start_trim
    return None


def _primary_search(points_xyz: list, min_start: int, max_end_trim: int = 3):
    """Recherche de la composante PRINCIPALE - demande explicite
    utilisateur ("often, the last steps at high temperature is a spurious
    magnetization that departs from zero, if the main component was
    going to the origin, it is best to discard the last points and keep
    the vector anchored to the origin") : essaie D'ABORD un fit ANCRE se
    terminant au tout DERNIER point ; si ce dernier point (ou les
    derniers) est une magnetisation parasite qui fait echouer l'ancrage
    des le depart (MAD/linearite), retire 1 a `max_end_trim` points du
    cote haute temperature/champ fort et reessaie ANCRE avant de se
    replier sur un fit LIBRE - preference EXPLICITE pour "ancre mais plus
    court" plutot que "libre mais complet".

    PUIS compare au meilleur ajustement de PLAN (grand cercle) sur la
    meme fenetre haute temperature/champ fort - demande explicite
    utilisateur, apres verification sur des specimens reels (site 14NQ04,
    Tibet) ou l'interpretation manuelle utilisait un plan pour cette
    composante alors que la recherche LIGNE seule echouait entierement ou
    ne trouvait qu'un ajustement mediocre. Le plan n'est retenu que
    lorsqu'il fait clairement mieux qu'une ligne inexistante ou de qualite
    "marginal"/"poor" (voir grade_mad) - une ligne "good"/"excellent"
    reste toujours preferee (une direction unique est plus utile qu'un
    grand cercle des qu'elle est bien resolue). PORTEE VOLONTAIREMENT
    LIMITEE au niveau du specimen individuel - la combinaison de plusieurs
    plans/lignes de PLUSIEURS specimens en une direction de site est un
    probleme distinct et non trivial (l'intersection de grands cercles
    peut converger sur la direction d'une composante secondaire bien
    groupee plutot que sur la primaire - demande explicite utilisateur,
    laisse pour un travail separe), pas traite ici.

    Retourne (start, end, fit, anchored, n_trimmed, kind) ou None -
    kind = "line" ou "plane"."""
    n = len(points_xyz)
    line_result = None
    for end_trim in range(0, max_end_trim + 1):
        end = n - 1 - end_trim
        if end - min_start + 1 < MIN_POINTS:
            break
        candidate = _greedy_backward(points_xyz[:end + 1], anchored=True, min_start=min_start)
        if candidate is not None:
            start, cend, fit = candidate
            line_result = (start, cend, fit, True, end_trim)
            break
    if line_result is None:
        candidate = _greedy_backward(points_xyz, anchored=False, min_start=min_start)
        if candidate is not None:
            start, end, fit = candidate
            line_result = (start, end, fit, False, 0)

    line_grade = grade_mad(line_result[2]["mad"], _LINE_MAD_THRESHOLD) if line_result is not None else None
    if line_result is not None and line_grade in ("excellent", "good"):
        start, end, fit, anchored, n_trimmed = line_result
        return start, end, fit, anchored, n_trimmed, "line"

    plane_result = _greedy_backward_plane(points_xyz, min_start=min_start)
    if plane_result is not None:
        start, end, fit = plane_result
        plane_grade = grade_mad(fit["mad"], _PLANE_MAD_THRESHOLD)
        if plane_grade in ("excellent", "good"):
            return start, end, fit, True, 0, "plane"

    if line_result is not None:
        start, end, fit, anchored, n_trimmed = line_result
        return start, end, fit, anchored, n_trimmed, "line"
    if plane_result is not None:
        start, end, fit = plane_result
        return start, end, fit, True, 0, "plane"
    return None


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
    un PLAN a 4 points comme dans `_greedy_backward_plane`, absorbe
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
    # `<= MIN_POINTS` (pas `< MIN_POINTS`, code mort - une recherche
    # gloutonne ne retourne JAMAIS moins de MIN_POINTS points, la
    # condition originale ne se declenchait donc jamais) : signale
    # explicitement un ajustement au tout minimum de points acceptes,
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
    )


def propose_components(ech) -> List[ComponentSuggestion]:
    """Propose 0, 1 ou 2 composantes (primary, et secondary si des points
    restent disponibles avant elle) pour ce specimen. Ne modifie rien,
    n'ecrit rien - voir docstring module.

    Le NRM BRUT (TOUTE PREMIERE mesure de la sequence, non traitee) est
    TOUJOURS exclu des deux recherches, INCONDITIONNELLEMENT (pas
    seulement si cod1=='N' - demande explicite utilisateur : "by default,
    it is best to not include the NRM0... start the interpretation at
    best at the first demag step") - convention paleomagnetique standard,
    PAS un critere statistique : verifie sur donnees reelles
    (old_pmag.ren) que le MAD change a peine (0.54 -> 0.45 deg) selon
    qu'on l'inclue ou non, alors que l'interpretation humaine l'exclut
    systematiquement (le NRM brut porte souvent une composante visqueuse
    non representative, retiree des le premier palier de traitement,
    meme s'il "tombe" pres de la droite ajustee par coincidence
    statistique)."""
    mesures = _zijderveld_measurements(ech)
    if len(mesures) < MIN_POINTS:
        return []
    points_xyz = [(m.x, m.y, m.z) for m in mesures]
    min_start = 1 if len(mesures) > MIN_POINTS else 0

    primary_result = _primary_search(points_xyz, min_start=min_start)

    primary_start = len(mesures)
    primary_suggestion = None
    if primary_result is not None:
        start, end, fit, anchored, n_trimmed, kind = primary_result
        primary_start = start
        primary_suggestion = _make_suggestion("primary", mesures, start, end, anchored, fit, n_trimmed, kind)

    secondary_result = _secondary_search(points_xyz, min_start=min_start, limit=primary_start)
    secondary_suggestion = None
    if secondary_result is not None:
        start, end, fit, secondary_anchored, n_trimmed_leading = secondary_result
        secondary_suggestion = _make_suggestion(
            "secondary", mesures, start, end, secondary_anchored, fit,
            n_trimmed_leading=n_trimmed_leading,
        )

    # "by default, if the component is anchored to the origin, use
    # primary" (demande explicite utilisateur) : une composante ancree a
    # l'origine EST par convention la caracteristique/terminale - si la
    # recherche "secondaire" (avant) finit ancree alors que la
    # "principale" (arriere) ne l'est pas, on echange simplement les
    # etiquettes plutot que le sens de recherche (start/end restent
    # corrects, seul le label et le numcomp implicite changent).
    suggestions: List[ComponentSuggestion] = [s for s in (primary_suggestion, secondary_suggestion) if s is not None]
    if (
        primary_suggestion is not None and secondary_suggestion is not None
        and secondary_suggestion.anchored and not primary_suggestion.anchored
    ):
        primary_suggestion.label, secondary_suggestion.label = "secondary", "primary"
        suggestions = [secondary_suggestion, primary_suggestion]

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
    hors de portee ici (voir _primary_search : "la combinaison de
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
