"""
Equivalent Python d'une premiere tranche du menu Fortran "Calcul" :
ajustement de droites par ACP (subroutines `linear`/`eigen` de linesplans.f,
pilotees interactivement par `ajuslig`) et statistiques de Fisher
(subroutine `fisher`/`cpolar` de calcul.f, utilisees par `fishmes`/`fishres`).

Travaille sur les SelectedSample/Measurement de selection.py.
"""

import math
import os
import re
import string
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Tuple

import numpy as np

from selection import (
    Measurement, SelectedSample, apply_orientation, corpen, polere, select_samples,
    normalized_intensity, _fmt_step,
)
from testlect import Pmag

# Voir app.HEADER_MARK / selection._HEADER_MARK - marque une ligne de
# titres de colonnes pour un affichage en gras cote console (STARpaleomagApp.
# _afficher), meme valeur, redefinie localement dans chaque module plutot
# qu'importee (evite tout couplage sur un nom prive).
_HEADER_MARK = "\x01"


# ---------------------------------------------------------------------------
# Equivalent de la structure /Resultats/ (starmac_OSX.inc), limite aux
# champs utilises par un ajustement de droite (cat1='L').
# ---------------------------------------------------------------------------

@dataclass
class FitResult:
    """Un ajustement de droite (equivalent d'un enregistrement `res`/`tr`)."""
    # identifiant anti-collision (equivalent res.c), rempli a l'archivage -
    # "<specimen>_<lettre>" (voir _next_specimen_c), PAS un entier aleatoire
    # 0-99999 (l'ancien schema, remplace suite a un fichier legacy reel
    # corrompu - voir _next_specimen_c pour le detail).
    c: str = ""
    id: str = ""
    cin: float = 0.0
    caz: float = 0.0
    dip: float = 0.0
    str_: float = 0.0
    cat1: str = "L"
    cat2: str = " "
    orig: str = "n"        # 'o' = ancre a l'origine, 'n' = non ancre
    demag: str = ""
    numcomp: int = 1
    # Etiquette de COMPOSANTE de magnetisation (A/B/C...) - PAS dans le
    # Fortran, demande explicite utilisateur ("when there is different
    # components of magnetizations within the same site, we may have two
    # or three means by site... We need to add a column component A,B,C").
    # Distinct de `numcomp` (numero d'ajustement PCA 1-9, deja existant) :
    # `numcomp` designe QUEL segment de mesures a ete ajuste sur un
    # specimen donne, `component` designe la nature physique du signal
    # retenu (ex. surimpression basse temperature vs ChRM haute
    # temperature) - assigne explicitement, jamais derive de `numcomp`
    # ("I think it is best not to use the numcomp of individual samples
    # for the mean"). Meme concept que sites.dir_comp_name/
    # specimens.dir_comp du modele MagIC (verifie live). Defaut "A" (la
    # composante principale/unique quand un seul jeu existe par site).
    component: str = "A"
    nb: int = 0
    dec: float = 0.0       # declinaison/inclinaison BRUTES (repere echantillon,
    inc: float = 0.0       # non corrigees - comme res.dec/res.inc en Fortran)
    mad: float = 0.0       # maximum angular deviation (deg)
    step_first: int = 0
    step_last: int = 0
    tx: Tuple[float, float] = (0.0, 0.0)  # extremites du segment (pour tracer
    ty: Tuple[float, float] = (0.0, 0.0)  # la droite ajustee sur un Zijderveld)
    tz: Tuple[float, float] = (0.0, 0.0)
    # Champs specifiques aux resultats agreges "mean:" (moyenne de site,
    # equivalent id="mean: <site>", cat1='F', cat2='i' - Fisher). Inutilises
    # pour un resultat normal (L/P/f/s). par2_mean/par3_mean : memes
    # positions de colonnes que par2/par3 (step_first/step_last), mais des
    # statistiques continues (pas des numeros d'etape entiers) - champs
    # separes pour ne pas perdre leur precision decimale au passage par
    # step_first/step_last (entiers, arrondis) lors d'un aller-retour
    # lecture/ecriture. par4/par5 = 2e paire lat/long ou VGP selon le
    # contexte, `liste` = "codes:c1:c2:..." (les `c` des resultats
    # individuels combines dans la moyenne).
    lat: float = 0.0
    rlong: float = 0.0
    par2_mean: float = 0.0
    par3_mean: float = 0.0
    par4: float = 0.0
    par5: float = 0.0
    # Ovale de confiance du VGP (demi-axes, degres) - voir dp_dm_from_a95.
    # Uniquement pour un "mean:" (comme par4/par5 juste au-dessus).
    vgp_dp: float = 0.0
    vgp_dm: float = 0.0
    # Nombre de resultats ligne/plan combines dans la moyenne - PAS dans le
    # Fortran (qui stockait un compte equivalent directement dans tr(i).
    # tx(1)/tx(2), repurpose ici pour k - voir _format_mean_line), calcule
    # jusqu'ici seulement A L'AFFICHAGE par recoupement fragile de `liste`
    # contre les resultats actuellement charges (list_results/
    # _mean_line_plane_counts - "?" si les specimens ne sont pas dans la
    # selection courante). Demande explicite utilisateur ("add a column
    # l/p with the number of lines and planes used in the mean Fisher for
    # the site... useful in the table for the publication") : desormais
    # PERSISTE avec la moyenne des l'archivage (voir
    # build_site_mean_result), fiable independamment de ce qui est charge
    # au moment de l'affichage. -1 = non renseigne (colonne reservee mais
    # pas peuplee - moyennes importees d'un autre format, ex. MagIC :
    # "This column will only be populated from now on. But during import
    # of other files keep the place for this column" - jamais deduit
    # apres coup pour ces imports, meme si l'information serait
    # techniquement calculable).
    n_lines: int = -1
    n_planes: int = -1
    liste: str = ""


@dataclass
class FisherStats:
    """Statistiques de Fisher (1953) sur un jeu de directions."""
    n: int
    dec: float
    inc: float
    r: float
    k: float
    a95: float
    csd: float
    # Nombre de directions inversees par combine_antipodal_groups avant le
    # calcul (0 si l'option n'a pas ete demandee, ou si un seul mode a ete
    # detecte - rien a inverser).
    n_flipped: int = 0


# ---------------------------------------------------------------------------
# ajuslig / linear / eigen : ajustement de droite par ACP (Kirschvink 1980)
# ---------------------------------------------------------------------------

def linear_fit(
    points: List[Tuple[float, float, float]], anchored: bool, mad_threshold: float = 15.0,
) -> Optional[dict]:
    """Equivalent de la subroutine `linear` : ajustement par ACP d'une droite
    passant (si `anchored`) par l'origine, sur les points (x,y,z) fournis
    (repere echantillon, dans l'ordre des etapes).

    Retourne None si le test de linearite ou le MAD (>`mad_threshold` deg,
    15 par defaut pour ajuslig/Zijderveld) echouent (equivalent de
    itestlin<0 en Fortran), sinon un dict avec direction (vecteur propre
    principal, unitaire), mad, nb, et les extremites (point1, point2) du
    segment ajuste pour affichage. `linearpal` (paleointensite, meme
    algorithme) utilise un seuil de 35 - voir fit_arai_direction.
    """
    pts = np.asarray(points, dtype=float)
    n = len(pts)
    if n == 0:
        return None

    anchor = np.zeros(3) if anchored else pts[-1]
    # Le centre utilise pour l'ACP est l'origine si ancre, sinon le centroide
    # (c'est le comportement d'origine de `linear` : sx/sy/sz ne sont
    # accumules que si `ancr .ne. 'o'`).
    center = np.zeros(3) if anchored else pts.mean(axis=0)

    # Test de linearite : distance directe (premier point -> ancre) comparee
    # a la somme des longueurs des segments consecutifs le long du chemin.
    sl = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
    if anchored:
        sl += float(np.linalg.norm(pts[-1]))
        mom = float(np.linalg.norm(pts[0]))
    else:
        mom = float(np.linalg.norm(pts[0] - anchor))
    linearity_ratio = (mom / sl) if sl > 0 else 1.0
    if sl > 0 and linearity_ratio < 0.5:
        return None
    if mom == 0:
        return None

    diffs = pts - center
    cov = diffs.T @ diffs

    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    l1, l2, l3 = eigvals[order]
    direction = eigvecs[:, order[0]]

    if l1 <= 0:
        return None
    mad = math.degrees(math.atan(math.sqrt(abs((l2 + l3) / l1))))
    if mad > mad_threshold:
        return None

    first_delta = pts[0] - anchor
    if np.dot(direction, first_delta) < 0:
        direction = -direction

    dot = float(np.dot(direction, center - anchor))
    point1 = center - dot * first_delta / mom
    point2 = point1 + mom * direction

    return {
        "direction": direction,
        "mad": mad,
        "nb": n,
        "linearity_ratio": linearity_ratio,
        "point1": tuple(point1),
        "point2": tuple(point2),
    }


def fit_line(
    ech: SelectedSample,
    jdeb: int,
    jfin: int,
    anchored: bool = True,
    numcomp: int = 1,
) -> Optional[FitResult]:
    """Equivalent non-interactif de `ajuslig` pour un seul segment : ajuste
    une droite sur ech.mesures[jdeb-1:jfin] (indices 1-based, comme affiches
    par list_measurements). Retourne None si le fit echoue (MAD>15 ou test
    de linearite, comme itestlin<0 en Fortran)."""
    subset = ech.mesures[jdeb - 1:jfin]
    fit = linear_fit([(m.x, m.y, m.z) for m in subset], anchored)
    if fit is None:
        return None

    dx, dy, dz = fit["direction"]
    _, dec, inc = polere(dx, dy, dz)

    return FitResult(
        id=ech.id,
        cin=ech.cin,
        caz=ech.caz,
        dip=ech.dip,
        str_=ech.str_,
        cat1="L",
        cat2=" ",
        orig="o" if anchored else "n",
        demag=subset[-1].cod1,
        numcomp=numcomp,
        nb=fit["nb"],
        dec=dec,
        inc=inc,
        mad=fit["mad"],
        step_first=subset[0].etape,
        step_last=subset[-1].etape,
        tx=(float(fit["point1"][0]), float(fit["point2"][0])),
        ty=(float(fit["point1"][1]), float(fit["point2"][1])),
        tz=(float(fit["point1"][2]), float(fit["point2"][2])),
    )


def _read_int(prompt: str, default: Optional[int] = None) -> Optional[int]:
    chaine = input(prompt).strip()
    if not chaine:
        return default
    try:
        return int(chaine)
    except ValueError:
        return default


def fit_line_interactive(selected: List[SelectedSample]) -> List[FitResult]:
    """Equivalent interactif de `ajuslig` : pour chaque echantillon
    selectionne, affiche la liste des etapes, demande le premier/dernier
    step a utiliser (0 = passer a l'echantillon suivant), l'ancrage a
    l'origine et le numero de composante, puis affiche le resultat."""
    results: List[FitResult] = []

    for ech in selected:
        if len(ech.mesures) < 2:
            continue

        steps = [m.etape for m in ech.mesures]
        print(f"\n{ech.id} :")
        print("  " + "  ".join(f"{i + 1}:{s}" for i, s in enumerate(steps)))

        jdeb = _read_int(" First step (1-based, 0=skip): ", 0)
        if not jdeb:
            continue
        jfin = _read_int(" Final step: ", None)
        if jfin is None or jfin < jdeb or jfin > len(steps):
            continue

        ancr = input(" anchored to the origin (Y/n)? ").strip().lower()
        anchored = ancr != "n"

        numcomp = _read_int(" component number (1-9)? ", 1) or 1

        fit = fit_line(ech, jdeb, jfin, anchored=anchored, numcomp=numcomp)
        if fit is None:
            print(" negative linear test (MAD too high or non-linear trend)")
            continue

        print(
            f"  dec={fit.dec:5.1f}  inc={fit.inc:5.1f}  mad={fit.mad:5.1f}  "
            f"nb points: {fit.nb}"
        )
        results.append(fit)

    return results


# ---------------------------------------------------------------------------
# fisher / cpolar : statistiques de Fisher (1953)
# ---------------------------------------------------------------------------

def fisher_mean(directions: List[Tuple[float, float]]) -> FisherStats:
    """Equivalent de `fisher`+`cpolar` : moyenne de Fisher d'une liste de
    directions (dec, inc) en degres."""
    n = len(directions)
    if n == 0:
        raise ValueError("empty list of directions")

    xsum = ysum = zsum = 0.0
    for dec, inc in directions:
        decr, incr = math.radians(dec), math.radians(inc)
        xsum += math.cos(incr) * math.cos(decr)
        ysum += math.cos(incr) * math.sin(decr)
        zsum += math.sin(incr)

    r = math.sqrt(xsum ** 2 + ysum ** 2 + zsum ** 2)
    if r == 0:
        raise ValueError("zero resultant - directions undefined")
    x, y, z = xsum / r, ysum / r, zsum / r

    mean_inc = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    mean_dec = math.degrees(math.atan2(y, x)) % 360.0

    if n > 1 and (n - r) > 0:
        k = (n - 1) / (n - r)
    else:
        k = 1.0

    if n > 1:
        fac = 1.0 / (n - 1)
        p = 20 ** fac - 1
        a = 1 - ((n - r) / r) * p
        a95 = math.degrees(math.acos(max(-1.0, min(1.0, a)))) if abs(a) < 1.0 else 90.0
    else:
        a95 = 180.0

    if n > 1:
        tet = 0.0
        for dec, inc in directions:
            decr, incr = math.radians(dec), math.radians(inc)
            xi = math.cos(incr) * math.cos(decr)
            yi = math.cos(incr) * math.sin(decr)
            zi = math.sin(incr)
            chord = math.sqrt((xi - x) ** 2 + (yi - y) ** 2 + (zi - z) ** 2)
            tet += (2 * math.asin(min(1.0, chord / 2))) ** 2
        csd = math.degrees(math.sqrt(tet / (n - 1)))
    else:
        csd = 0.0

    return FisherStats(n=n, dec=mean_dec, inc=mean_inc, r=r, k=k, a95=a95, csd=csd)


def combine_antipodal_groups(
    directions: List[Tuple[float, float]],
) -> Tuple[List[Tuple[float, float]], int]:
    """Detecte deux groupes de directions quasi-antipodaux (typiquement les
    intervalles de polarite normale/inverse d'une courte section
    magnetostratigraphique) et inverse (dec+180, -inc) le groupe le MOINS
    nombreux, pour permettre le calcul d'une seule moyenne de Fisher
    combinee - demande explicite utilisateur ("when in the results there
    are two nearly antipodal groups, best to invert the group with less
    results. This is the case in a short magnetostratigraphic section").

    Meme principe que `pmag.flip` (pmagpy, code disponible dans
    StereoUtils_Py) : l'axe de reference est le vecteur propre dominant de
    la matrice d'orientation (bidirectionnelle - construite a partir de
    x.x^T, donc insensible a la polarite de chaque direction). Un jeu de
    directions unimodal (un seul mode) donne un axe proche de sa propre
    direction moyenne : toutes les directions tombent alors du meme cote
    (angle <= 90 deg) et rien n'est inverse, ce qui evite de retourner par
    erreur un simple point isole/aberrant d'un jeu par ailleurs unimodal.

    Retourne (directions ajustees, nombre de directions inversees) - la
    liste d'entree n'est jamais modifiee en place."""
    n = len(directions)
    if n < 2:
        return list(directions), 0

    vecs = np.array([
        (
            math.cos(math.radians(inc)) * math.cos(math.radians(dec)),
            math.cos(math.radians(inc)) * math.sin(math.radians(dec)),
            math.sin(math.radians(inc)),
        )
        for dec, inc in directions
    ])
    orient = vecs.T @ vecs
    eigvals, eigvecs = np.linalg.eigh(orient)
    axis = eigvecs[:, int(np.argmax(eigvals))]

    angles = np.degrees(np.arccos(np.clip(vecs @ axis, -1.0, 1.0)))
    group_a = [k for k in range(n) if angles[k] <= 90.0]
    group_b = [k for k in range(n) if angles[k] > 90.0]
    if not group_a or not group_b:
        return list(directions), 0

    minority = group_b if len(group_b) <= len(group_a) else group_a
    adjusted = list(directions)
    for k in minority:
        dec, inc = adjusted[k]
        adjusted[k] = ((dec + 180.0) % 360.0, -inc)
    return adjusted, len(minority)


# cat1 traites comme des LIGNES par `fishres` (calcul.f:84-142, meme
# regroupement pour 'L','f' ET 's' - seul 'P' est traite a part comme un
# pole de grand cercle) - voir fisher_from_results/fisher_combine_lines_planes.
_LINE_LIKE_CAT1 = ("L", "f", "s")


def fisher_combine_lines_planes(
    lines: List[Tuple[float, float]], planes: List[Tuple[float, float]],
) -> FisherStats:
    """Port de `fishgc` (calcul.f:169-515, verifie contre le source Fortran
    - demande explicite utilisateur "check the fortran") : moyenne de
    Fisher combinant des directions de ligne (`lines`) et des poles de
    grand cercle (`planes`, ex. le pole retourne par `fit_plane`) selon
    McFadden & McElhinny (EPSL, 1988, v87) - message imprime par le
    Fortran lui-meme ("Mean direction following McFadden et McElhinny...
    combining great circles and directions"). Chaque plan est
    iterativement reprojete sur son grand cercle, au point le plus proche
    de la moyenne courante (Gauss-Seidel : le point du plan en cours de
    mise a jour est d'abord RETIRE de la somme courante avant d'etre
    recalcule, puis rajoute - meme sequence que le Fortran, PAS une simple
    recomputation en parallele), jusqu'a 100 iterations ou convergence
    (`r1 < r0*(1+1e-6)` apres au moins 16 iterations - memes constantes
    que le Fortran).

    Fonctionnalite "secteur" du Fortran (restreindre la reprojection a un
    ARC plutot qu'au grand cercle complet, via ad1/ai1/ad2/ai2) : PAS
    portee - entierement DESACTIVEE dans le source lui-meme (`iyes='n'`
    code en dur, toute la saisie/logique de secteur est en commentaires)
    - code mort, jamais atteignable meme en Fortran.

    `csd` n'est pas calcule par `fishgc` dans ce cas combine (contrairement
    a `fisher_mean` seule) - retourne toujours 0.0 ici. Si `planes` est
    vide, les formules k/a95 de McFadden&McElhinny se reduisent
    algebriquement a celles de Fisher standard (verifie) - mais cette
    fonction ne prend PAS ce raccourci elle-meme, c'est a l'appelant de
    preferer `fisher_mean` directement dans ce cas (voir
    fisher_from_results) pour eviter 100 iterations pour rien."""
    m, imm = len(lines), len(planes)
    if m + imm < 2:
        raise ValueError("need at least 2 combined line/plane results for a Fisher mean")

    def _unit(dec: float, inc: float) -> Tuple[float, float, float]:
        decr, incr = math.radians(dec), math.radians(inc)
        return math.cos(incr) * math.cos(decr), math.cos(incr) * math.sin(decr), math.sin(incr)

    def _project(u0: float, v0: float, w0: float, pqr: Tuple[float, float, float]):
        p, q, r = pqr
        tau = u0 * p + v0 * q + w0 * r
        # garde anti-crash (PAS dans le Fortran, ou une division IEEE par 0
        # produirait Inf plutot qu'une exception Python) - n'affecte que le
        # cas degenere ou la moyenne courante coincide EXACTEMENT avec le
        # pole d'un plan (mathematiquement indefini : tout point du grand
        # cercle est alors equidistant).
        ro = math.sqrt(max(1e-12, 1.0 - tau * tau))
        return (u0 - tau * p) / ro, (v0 - tau * q) / ro, (w0 - tau * r) / ro

    sx = sy = sz = 0.0
    if m == 0:
        # Fortran : demande une direction de depart (defaut (0,0) si
        # l'utilisateur ne saisit rien) quand aucune ligne n'est disponible -
        # equivalent headless : toujours le defaut (0,0).
        u0, v0, w0 = 1.0, 0.0, 0.0
    else:
        for dec, inc in lines:
            x, y, z = _unit(dec, inc)
            sx += x; sy += y; sz += z
        norm = math.sqrt(sx * sx + sy * sy + sz * sz)
        u0, v0, w0 = sx / norm, sy / norm, sz / norm

    pqr_list = [_unit(dec, inc) for dec, inc in planes]
    xg = [0.0] * imm
    yg = [0.0] * imm
    zg = [0.0] * imm

    for i in range(imm):
        xg[i], yg[i], zg[i] = _project(u0, v0, w0, pqr_list[i])
        sx += xg[i]; sy += yg[i]; sz += zg[i]
        norm = math.sqrt(sx * sx + sy * sy + sz * sz)
        u0, v0, w0 = sx / norm, sy / norm, sz / norm
    r0 = math.sqrt(sx * sx + sy * sy + sz * sz)

    r1 = r0
    for itr in range(1, 101):
        for j in range(imm):
            sx -= xg[j]; sy -= yg[j]; sz -= zg[j]
            norm = math.sqrt(sx * sx + sy * sy + sz * sz)
            u0, v0, w0 = sx / norm, sy / norm, sz / norm
            xg[j], yg[j], zg[j] = _project(u0, v0, w0, pqr_list[j])
            sx += xg[j]; sy += yg[j]; sz += zg[j]
            norm = math.sqrt(sx * sx + sy * sy + sz * sz)
            u0, v0, w0 = sx / norm, sy / norm, sz / norm
        r1 = math.sqrt(sx * sx + sy * sy + sz * sz)
        if r1 < r0 * (1.0 + 1e-6) and itr > 15:
            break
        r0 = r1

    _, dec, inc = polere(sx, sy, sz)
    # Cas structurellement 0/0 (m==0 et imm==2 : deux grands cercles SANS
    # ligne convergent generiquement sur un unique point d'intersection
    # exact, r1==2==m+imm) : le Fortran (IEEE) obtiendrait NaN silencieusement,
    # Python leve ZeroDivisionError sur une division flottante par 0 - meme
    # resultat (NaN) obtenu ici explicitement plutot que de planter.
    denom = 2 * (m + imm - r1)
    rk = (2 * m + imm - 2) / denom if abs(denom) > 1e-12 else float("nan")
    # nprim = m + imm/2 (division ENTIERE, comme en Fortran - chaque plan
    # compte pour une demi-donnee dans le calcul des degres de liberte).
    nprim = m + imm // 2
    if nprim <= 1:
        # Fortran : `1/(float(nprim)-1)` -> IEEE Inf, `(nprim-1)*Inf` -> NaN
        # (nprim-1 vaut 0 comme entier), NaN se propage jusqu'au test
        # `abs(alpha95)<1.0` (toujours FAUX pour NaN) -> alpha95=90.0. Meme
        # resultat obtenu ici explicitement (Python leve ZeroDivisionError
        # sur une division flottante par 0, contrairement au Fortran/IEEE).
        alpha95 = 90.0
    else:
        fac = 1.0 / (nprim - 1)
        p_ = 20 ** fac - 1
        a = 1 - (nprim - 1) * p_ / (rk * r1)
        alpha95 = math.degrees(math.acos(a)) if abs(a) < 1.0 else 90.0

    return FisherStats(n=m + imm, dec=dec, inc=inc, r=r1, k=rk, a95=alpha95, csd=0.0)


def fisher_from_measurements(
    selected: List[SelectedSample], orientation: int = 1, combine_antipodal: bool = False,
) -> FisherStats:
    """Equivalent de `fishmes` : moyenne de Fisher calculee directement sur
    TOUTES les mesures de la selection (chaque etape de demagnetisation
    compte comme une direction independante), apres correction d'orientation.
    `combine_antipodal` : voir combine_antipodal_groups (option, defaut
    False - ne change pas le comportement existant par defaut)."""
    directions = []
    for ech in selected:
        for m in ech.mesures:
            xx, yy, zz = apply_orientation(m.x, m.y, m.z, ech, orientation)
            _, dec, inc = polere(xx, yy, zz)
            directions.append((dec, inc))
    n_flipped = 0
    if combine_antipodal:
        directions, n_flipped = combine_antipodal_groups(directions)
    stats = fisher_mean(directions)
    return replace(stats, n_flipped=n_flipped)


def fisher_from_results(
    results: List[FitResult], orientation: int = 1, combine_antipodal: bool = False,
) -> FisherStats:
    """Equivalent de `fishres`+`fishgc` (calcul.f) - demande explicite
    utilisateur ("when lines and planes are selected for a fisher result,
    use the combined L & P", verifiee contre le source Fortran : "check
    the fortran"). Sur les resultats de type ligne (cat1 in 'L'/'f'/'s' -
    meme regroupement que `fishres`, calcul.f:84-142 - PAS uniquement
    'L') ET de type plan (cat1=='P', combines via McFadden & McElhinny -
    voir fisher_combine_lines_planes), en reappliquant la correction
    d'orientation de chaque resultat (son propre cin/caz/dip/str, via
    `_correct_dec_inc` - la meme transformation que le Fortran applique
    UNIFORMEMENT a tr(j) avant de classer en ligne/plan). Auparavant
    limitee aux seuls resultats 'L' - les plans (et 'f'/'s') etaient
    silencieusement ignores.

    `combine_antipodal` : voir combine_antipodal_groups (option, defaut
    False) - applique UNIQUEMENT aux directions de type ligne (un pole de
    grand cercle n'a pas de "polarite" a inverser)."""
    lines: List[Tuple[float, float]] = []
    planes: List[Tuple[float, float]] = []
    for res in results:
        if res.cat1 in _LINE_LIKE_CAT1:
            lines.append(_correct_dec_inc(res, orientation))
        elif res.cat1 == "P":
            planes.append(_correct_dec_inc(res, orientation))
    if not lines and not planes:
        raise ValueError("no line or plane result ('L'/'P'/'f'/'s') in the list")

    n_flipped = 0
    if combine_antipodal and lines:
        lines, n_flipped = combine_antipodal_groups(lines)

    stats = fisher_mean(lines) if not planes else fisher_combine_lines_planes(lines, planes)
    return replace(stats, n_flipped=n_flipped)


# ---------------------------------------------------------------------------
# ajusplans / planar : ajustement de plan (grand cercle) par ACP - meme
# principe que linear_fit, mais le pole du plan est le vecteur propre de la
# PLUS PETITE valeur propre (linear_fit prend la plus grande). Pas de test de
# linearite (n'a pas de sens pour un plan) ; seuil MAD 25 deg (pas 15).
# ---------------------------------------------------------------------------

def planar_fit(points: List[Tuple[float, float, float]], normalize: bool = True) -> Optional[dict]:
    """Equivalent de la subroutine `planar` : ajustement par ACP d'un plan
    (grand cercle) passant par les points fournis. `normalize` : equivalent
    de `norm=='o'` (points ramenes sur la sphere unite avant l'ACP - option
    par defaut de `ajusplans`). Retourne None si MAD>25 (equivalent
    itestplan=-1).

    BUG REEL corrige ici (demande explicite utilisateur : "there is a
    problem in the calculation of the great circles... seems a problem in
    calculation of the fit") : la matrice de dispersion etait calculee
    apres avoir SOUSTRAIT LA MOYENNE des points (`diffs = pts - center`),
    ce qui ajuste le plan le plus proche du NUAGE DE POINTS LUI-MEME (son
    propre centroide) - correct pour un plan 3D generique, mais FAUX pour
    un ajustement de grand cercle a la Kirschvink [1980] : un grand
    cercle est par definition un plan passant par le CENTRE de la sphere
    (l'origine), pas par le centroide des points de desaimantation
    (qui n'a aucune raison de coincider avec l'origine, surtout sur un
    arc partiel/asymetrique du cercle - cas frequent en donnees reelles).
    Verifie sur 102 plans reels (pmagCARF.r, Arriagada et al. 2006) :
    l'ancienne formule donnait un pole a plus de 15 deg. du veritable
    grand cercle (ecart median au pole de 90 deg. : 19.7 deg., vs 1.9 deg.
    pour le meme plan tel qu'archive par le Fortran d'origine) pour 50 des
    102 plans ; retirer le centrage (utiliser directement `pts.T @ pts`,
    sans soustraire `pts.mean(axis=0)`) ramene l'ecart median a 1.7 deg.,
    conforme a l'archive Fortran."""
    pts = np.asarray(points, dtype=float)
    n = len(pts)
    if n == 0:
        return None
    if normalize:
        norms = np.linalg.norm(pts, axis=1)
        norms[norms == 0] = 1.0
        pts = pts / norms[:, None]

    cov = pts.T @ pts

    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    l1, l2, l3 = eigvals[order]
    if l1 <= 0 or l2 <= 0:
        return None
    pole = eigvecs[:, order[2]]  # plus petite valeur propre = pole du plan

    mad = math.degrees(math.atan(math.sqrt(abs(l3 / l2 + l3 / l1))))
    if mad > 25.0:
        return None

    return {"pole": pole, "mad": mad, "nb": n}


def fit_plane(
    ech: SelectedSample,
    jdeb: int,
    jfin: int,
    normalize: bool = True,
    numcomp: int = 1,
) -> Optional[FitResult]:
    """Equivalent non-interactif de `ajusplans` pour un seul segment (indices
    1-based, comme fit_line). `ancr` est toujours 'o' cote Fortran (hardcode,
    pas de prompt) - reproduit ici (orig="o"), sauf si l'appelant l'ecrase
    explicitement (cas de ajusligredo, ou le fichier redo peut demander 'n')."""
    subset = ech.mesures[jdeb - 1:jfin]
    fit = planar_fit([(m.x, m.y, m.z) for m in subset], normalize=normalize)
    if fit is None:
        return None

    px, py, pz = fit["pole"]
    _, dec, inc = polere(px, py, pz)

    return FitResult(
        id=ech.id,
        cin=ech.cin,
        caz=ech.caz,
        dip=ech.dip,
        str_=ech.str_,
        cat1="P",
        cat2=" ",
        orig="o",
        demag=subset[-1].cod1,
        numcomp=numcomp,
        nb=fit["nb"],
        dec=dec,
        inc=inc,
        mad=fit["mad"],
        step_first=subset[0].etape,
        step_last=subset[-1].etape,
    )


# ---------------------------------------------------------------------------
# ajusfisher : direction moyenne de Fisher sur une plage d'etapes (au lieu
# d'un ajustement PCA) - reutilise fisher_mean() (deja verifiee identique a
# `fisher`/`cpolar`). Repli sur une direction brute unique (label 735 du
# Fortran) pour les echantillons a moins de 3 mesures.
# ---------------------------------------------------------------------------

def fit_single_direction(ech: SelectedSample, index: int = -1) -> FitResult:
    """Equivalent de la branche 'direction unique' de `ajusfisher` (label
    735) : utilise directement le vecteur d'UNE mesure (par defaut la
    derniere, comme le Fortran `jdeb=ech(i).ifin`) comme direction, sans
    moyenne de Fisher."""
    m = ech.mesures[index]
    _, dec, inc = polere(m.x, m.y, m.z)
    return FitResult(
        id=ech.id, cin=ech.cin, caz=ech.caz, dip=ech.dip, str_=ech.str_,
        cat1="s", cat2="s", orig="n", demag=m.cod1, numcomp=1,
        nb=1, dec=dec, inc=inc, mad=0.0,
        step_first=m.etape, step_last=m.etape,
    )


def fit_fisher_direction(
    ech: SelectedSample, jdeb: int, jfin: int, numcomp: int = 1
) -> FitResult:
    """Equivalent de `ajusfisher` (branche moyenne de Fisher) sur
    ech.mesures[jdeb-1:jfin]. `tx` porte (kappa, csd) - meme detournement de
    champ que le Fortran (res.tx(1)=rk, res.tx(2)=stv), pas des coordonnees
    de segment."""
    subset = ech.mesures[jdeb - 1:jfin]
    directions = []
    for m in subset:
        _, dec, inc = polere(m.x, m.y, m.z)
        directions.append((dec, inc))
    stats = fisher_mean(directions)

    return FitResult(
        id=ech.id, cin=ech.cin, caz=ech.caz, dip=ech.dip, str_=ech.str_,
        cat1="f", cat2=" ", orig="n",
        demag=subset[-1].cod1, numcomp=numcomp, nb=stats.n,
        dec=stats.dec, inc=stats.inc, mad=stats.a95,
        step_first=subset[0].etape, step_last=subset[-1].etape,
        tx=(stats.k, stats.csd),
    )


# ---------------------------------------------------------------------------
# ajusligauto : version automatisee de ajuslig - une droite par echantillon
# sur la totalite de ses mesures, meme ancrage pour tous, sauvegarde sans
# confirmation (echantillons a moins de 2 mesures ignores).
# ---------------------------------------------------------------------------

def fit_lines_auto(selected: List[SelectedSample], anchored: bool = True) -> List[FitResult]:
    """Equivalent de `ajusligauto`. `numcomp` vaut 1 si ancre, 2 sinon (meme
    regle que le Fortran : `res.numcomp=2 ; if(ancr.ne."n") res.numcomp=1`)."""
    numcomp = 1 if anchored else 2
    results = []
    for ech in selected:
        if len(ech.mesures) < 2:
            continue
        fit = fit_line(ech, 1, len(ech.mesures), anchored=anchored, numcomp=numcomp)
        if fit is not None:
            results.append(fit)
    return results


# ---------------------------------------------------------------------------
# ajusligredo : rejoue des ajustements (droite ou plan) a partir d'un fichier
# texte "redo" (un ajustement par ligne). Selection d'echantillon par id EXACT
# (equivalent simplifie de `seloneech` - pas de joker '*'/'?' ici, contrai-
# rement a selection.select_samples qui les gere deja pour le menu Selection
# donnees). Plage de steps resolue par VALEUR d'etape (pas position) : jdeb/
# jfin = dernier index dont l'etape est <= tempmin/tempmax (balayage lineaire
# "dernier qui correspond", pas une recherche exacte - fidele au Fortran).
# ---------------------------------------------------------------------------

def _resolve_step_range_by_value(
    ech: SelectedSample, tempmin: float, tempmax: float
) -> Optional[Tuple[int, int]]:
    jdeb = jfin = None
    for idx, m in enumerate(ech.mesures, start=1):
        if m.etape <= tempmin:
            jdeb = idx
        if m.etape <= tempmax:
            jfin = idx
    if jdeb is None or jfin is None:
        return None
    return jdeb, jfin


def _parse_redo_numcomp(token: str) -> Tuple[int, str]:
    """Derniere colonne d'une ligne de redo file (voir fit_from_redo_file) :
    accepte soit un ANCIEN numero (1-9, converti en lettre via
    _component_from_numcomp - MEME principe que l'import legacy, "change
    the component number to letter during import of legacy"), soit DEJA
    une lettre de composante (A/B/C...) - demande explicite utilisateur
    ("quand on utilise un redo file, transformer le numero de composante
    1 en A; 2 en B, 3 en C etc; si le redo file a deja les composantes en
    A,B etc, garder les lettres"). `numcomp` (l'entier, usage INTERNE
    seulement - voir FitResult.numcomp, ex. couleur du trace Zijderveld)
    est deduit de la lettre dans ce 2e cas (A=1, B=2, ...) plutot que
    fige a 1 pour toutes les composantes. Retourne (1, "A") si le token
    n'est ni un nombre ni une lettre seule."""
    token = token.strip()
    if token.isdigit():
        numcomp = int(token)
        return numcomp, _component_from_numcomp(numcomp)
    if len(token) == 1 and token.isalpha():
        letter = token.upper()
        return ord(letter) - ord("A") + 1, letter
    return 1, "A"


def fit_from_redo_file(
    pmag_list: List[Pmag], redo_lines: List[str]
) -> List[FitResult]:
    """Equivalent de `ajusligredo`. Chaque ligne : 'sample L|P o|n tempmin
    tempmax numcomp' (ex: 'NWD1-11A L n 400 560 2'). Lignes consecutives
    identiques ignorees (la seconde), comme le Fortran (chaine2==chaine).
    Echantillons a moins de 3 mesures ignores (comme `nbmes<3` cote Fortran).

    `pmag_list` : la TOTALITE des donnees chargees (equivalent self.donnees/
    pmag(:)), PAS une selection prealable - chaque ligne recharge elle-meme
    l'echantillon qu'elle designe via `call seloneech(samnum)` (nbech=0;
    nbmes=0; recherche dans pmag(:) par id, etapmin=0/etapmax=9999/
    demag='*'), independamment de toute selection en cours. Un fichier redo
    peut donc rejouer des ajustements sur des specimens jamais selectionnes
    manuellement au prealable.

    La derniere colonne (numcomp) est convertie en lettre de composante
    (voir _parse_redo_numcomp) et assignee a `FitResult.component` -
    demande explicite utilisateur (voir _parse_redo_numcomp)."""
    results: List[FitResult] = []
    last_line: Optional[str] = None

    for raw in redo_lines:
        line = raw.strip()
        if not line or line == last_line:
            last_line = line
            continue
        last_line = line

        parts = line.split()
        if len(parts) < 6:
            continue
        samnum, lp, ancr, tempmin_s, tempmax_s, numcomp_s = parts[:6]
        if lp not in ("L", "P"):
            continue
        matches = select_samples(pmag_list, samnum, verbose=False)
        ech = matches[0] if matches else None
        if ech is None or len(ech.mesures) < 3:
            continue
        numcomp, component = _parse_redo_numcomp(numcomp_s)
        try:
            tempmin, tempmax = float(tempmin_s), float(tempmax_s)
        except ValueError:
            continue

        rng = _resolve_step_range_by_value(ech, tempmin, tempmax)
        if rng is None:
            continue
        jdeb, jfin = rng
        anchored = ancr.lower() != "n"

        if lp == "L":
            fit = fit_line(ech, jdeb, jfin, anchored=anchored, numcomp=numcomp)
        else:
            fit = fit_plane(ech, jdeb, jfin, normalize=True, numcomp=numcomp)
            if fit is not None:
                # contrairement a ajusplans (ancr toujours 'o'), ajusligredo
                # derive orig de l'ancr du fichier meme pour un plan (le
                # calcul PCA lui-meme ignore ancr, seul le champ change).
                fit.orig = "o" if anchored else "n"
        if fit is not None:
            fit.component = component
            results.append(fit)

    return results


# ---------------------------------------------------------------------------
# Resultats : init / listage (equivalents de initres / lisres, limites aux
# lignes - lisres gere aussi les plans et n'est donc que partiellement porte)
# ---------------------------------------------------------------------------

def init_results() -> List[FitResult]:
    """Equivalent de `initres` : reinitialise la liste des resultats."""
    print("\n-- Init the results list -- ")
    return []


def _correct_dec_inc(res: FitResult, orientation: int) -> Tuple[float, float]:
    """Equivalent de `correct(rd,ri,cdip,caz,dip,str)` (StarUtil.f:66),
    reapplique a la direction BRUTE stockee dans un resultat (res.dec/
    res.inc, repere echantillon d'origine) - duplique de
    stereo.py:_correct_dec_inc pour eviter un import circulaire (stereo.py
    importe deja FitResult depuis ce module).

    Un resultat "mean:" est DEJA fige dans UNE orientation precise
    (par3_mean : 1=echantillon/2=in-situ/3=apres pendage) - retourne
    TOUJOURS la valeur stockee TELLE QUELLE pour un "mean:", quelle que
    soit `orientation` demandee. Une tentative anterieure de reprojeter
    une moyenne IS en TC (et reciproquement) via le pendage/strike du
    site a ete EXPLICITEMENT ECARTEE - demande utilisateur ("this is not
    exactly what is needed... only the means in IS are plotted if you
    are in IS and only the selected TC means if you are in TC") : le
    comportement voulu n'est PAS une reprojection mais un FILTRAGE - une
    moyenne archivee dans une orientation ne doit etre selectionnee/
    tracee QUE quand cette meme orientation est active (voir
    _load_mean_results et draw_stereo_results), jamais recalculee dans
    une autre."""
    if res.id[:5] == "mean:":
        return res.dec, res.inc
    incr, decr = math.radians(res.inc), math.radians(res.dec)
    x = math.cos(incr) * math.cos(decr)
    y = math.cos(incr) * math.sin(decr)
    z = math.sin(incr)
    xx, yy, zz = apply_orientation(x, y, z, res, orientation)
    _, dec, inc = polere(xx, yy, zz)
    return dec, inc


_ORIENT_HEADER = {
    1: "results in sample coordinates",
    2: "results in in-situ coordinates",
    3: "results after tilt correction",
}


_ORIENT_MODE_TAG = {1: "Sa", 2: "IS", 3: "TC"}


def dp_dm_from_a95(a95: float, inc: float) -> Tuple[float, float]:
    """dP,dM (demi-axes, en degres, de l'ovale de confiance du VGP derive
    d'une direction moyenne + a95) - dP le long du grand cercle site-pole,
    dM perpendiculaire. Meme formule que StereoUtils_Py
    (`stereo_pmagutils._dp_dm_dir_to_vgp`, elle-meme portee de
    pmagoutils.f:1758-1761, VGPC1/VGPC3) et que `_e95` ci-dessous (qui n'en
    prenait que la moyenne pour l'affichage `lisres`) - exposee separement
    pour l'archivage .pmagres (demande explicite utilisateur : "during the
    archive of the mean direction in pmagres, the VGP is calculated. Is it
    possible to add the dp dm ellipse of confidence derived from the a95,
    available in Stereo_Py"). Verifiee contre une VRAIE contribution MagIC
    publiant deja vgp_dp/vgp_dm (magic_contribution_20340.txt, site "AynC" :
    a95=6.4, inc=27.3 -> dp=3.8/dm=7.0 publies, dp=3.80/dm=6.97 calcules
    ici). `inc` est l'inclinaison de la DIRECTION moyenne (repere site),
    pas celle du VGP lui-meme."""
    if a95 <= 0:
        return 0.0, 0.0
    ai = math.radians(inc)
    pl = 1.57095 - math.atan(0.5 * math.tan(ai))
    dm = a95 * math.sin(pl) / math.cos(ai)
    dp = 2 * a95 * (1 / (1 + 3 * math.cos(ai) ** 2))
    return dp, dm


def _e95(a95: float, inc: float) -> float:
    """Equivalent exact de la formule e95 de `lisres` (dataselect.f, juste
    avant le format 201) : deduite de a95 et de l'inclinaison de la
    moyenne (paleolatitude via VGP), PAS d'une donnee stockee - calculable
    a l'affichage, comme dans le Fortran."""
    dp, dm = dp_dm_from_a95(a95, inc)
    return (dm + dp) / 2


def dir_to_vgp(dec: float, inc: float, site_lat: float, site_lon: float) -> Tuple[float, float]:
    """Direction (dec, inc) mesuree au site (site_lat, site_lon) -> pole VGP
    (vgp_lat, vgp_lon), tout en degres. Meme formule que StereoUtils_Py
    (`stereo_pmagutils.dodi_vgp`, portee de `pmagoutils.f:4001-4074`) -
    demande explicite utilisateur ("in fisher results, a mean is
    calculated for the site, if the user want to archive this mean, this
    mean needs to be recorded in the .pmagres file") : jusqu'ici seul
    l'import MagIC archivait un VGP (deja calcule par la contribution
    source, jamais recalcule ici - voir convert_magic_to_r), STARpaleomag_Py
    n'avait pas sa propre transformation direction->VGP."""
    dec_r, inc_r = math.radians(dec), math.radians(inc)
    slat_r, slon_r = math.radians(site_lat), math.radians(site_lon)
    p = math.atan2(2.0, math.tan(inc_r))
    plat = math.asin(
        math.sin(slat_r) * math.cos(p) + math.cos(slat_r) * math.sin(p) * math.cos(dec_r)
    )
    cos_plat = math.cos(plat)
    beta = (math.sin(p) * math.sin(dec_r)) / cos_plat if cos_plat else 0.0
    beta = math.asin(max(-1.0, min(1.0, beta)))
    if math.cos(p) >= math.sin(slat_r) * math.sin(plat):
        plon = slon_r + beta
    else:
        plon = slon_r + math.pi - beta
    plon %= 2.0 * math.pi
    return math.degrees(plat), math.degrees(plon)


def _find_site_specimen(site: str, donnees):
    """N'IMPORTE QUEL specimen de `donnees` appartenant a `site` - d'abord
    par `magic_site` (le champ 'site:' du .prmag, EDITABLE independamment
    de l'id via "Complete sample information" en mode site avec
    renommage - voir complete_sample_info.py), et SEULEMENT en repli par
    id.startswith(site) (ancienne convention "6 premiers caracteres de
    l'id = site", plus fiable qu'un magic_site jamais renseigne).

    BUG REEL corrige ici (demande explicite utilisateur, apres avoir
    renomme des sites - "in the calculation of a mean for a site, it does
    not use the information from the site") : magic_site est la source de
    verite du nom de site (peut differer de l'id, delibarement JAMAIS
    renomme lui-meme, voir complete_sample_info.py - "specimen ID kept
    untouched (site is just metadata, not derived from the ID)") ; ne
    chercher que par `id.startswith(site)` faisait echouer TOUTE recherche
    de lat/lon/pendage de site des qu'un site avait ete renomme (ex. ancien
    nom de campagne "95CH43" renomme en "P10" : aucun id de specimen ne
    commence plus par "P10", l'id garde son ancien prefixe) - le VGP
    archive retombait alors silencieusement sur 0.0/0.0."""
    if not donnees:
        return None
    ech = next((s for s in donnees if getattr(s, "magic_site", "") == site), None)
    if ech is not None:
        return ech
    return next((s for s in donnees if s.id.startswith(site)), None)


def site_of_result(res: FitResult, donnees) -> str:
    """Nom de site d'un resultat SPECIMEN (pas "mean:") : `magic_site` du
    specimen correspondant dans `donnees` (source de verite, voir
    _find_site_specimen) si trouve, sinon repli sur les 6 premiers
    caracteres de `res.id` (ancienne convention, comportement inchange si
    `donnees` n'est pas fourni ou si le specimen n'y est pas trouve)."""
    ech = next((s for s in (donnees or []) if s.id == res.id), None)
    site = getattr(ech, "magic_site", "") if ech is not None else ""
    return site or res.id[:6]


def site_lat_lon_from_donnees(site: str, donnees) -> Optional[Tuple[float, float]]:
    """Latitude/longitude du site, lues sur N'IMPORTE QUEL specimen de
    `donnees` appartenant a ce site (voir _find_site_specimen - une
    propriete de site, identique pour tous ses specimens). Retourne None
    si `donnees` n'est pas fourni ou si aucun specimen du site n'y est
    trouve."""
    ech = _find_site_specimen(site, donnees)
    if ech is None:
        return None
    return ech.lat, ech.rlong


def build_site_mean_result(
    stats: FisherStats,
    contributing: List[FitResult],
    site: str,
    orientation: int,
    site_lat: float = 0.0,
    site_lon: float = 0.0,
    component: str = "A",
) -> FitResult:
    """Construit un resultat "mean:" (moyenne de site) archivable, a partir
    des statistiques de Fisher (fisher_from_results) et des resultats
    individuels qui y contribuent - demande explicite utilisateur ("in
    fisher results, a mean is calculated for the site, if the user want to
    archive this mean, this mean needs to be recorded in the .pmagres
    file"). Jusqu'ici, seul l'import MagIC (convert_magic_to_r) produisait
    des lignes "mean:" ; aucune fonction n'existait pour archiver une
    moyenne calculee localement par l'application elle-meme.

    `site_lat`/`site_lon` : necessaires au calcul du VGP (dir_to_vgp) -
    0.0/0.0 si non disponibles (voir site_lat_lon_from_donnees), le VGP est
    alors archive mais denue de sens (0,0 n'est jamais un vrai site) -
    laisse a l'appelant le soin de prevenir l'utilisateur plutot que de
    bloquer l'archivage. `component` : etiquette de composante de
    magnetisation (A/B/C..., voir FitResult.component) - demandee
    explicitement a l'archivage, jamais deduite de `numcomp` des resultats
    individuels ("I think it is best not to use the numcomp of individual
    samples for the mean").

    n_lines/n_planes (colonne "L/P" du .pmagres) : comptes des resultats de
    `contributing` par type - cat1 in 'L'/'f'/'s' pour n_lines, cat1=='P'
    pour n_planes (meme regroupement que `fisher_from_results`/`fishres`,
    PAS uniquement 'L'/'P') - demande explicite utilisateur ("add a column
    l/p with the number of lines and planes used in the mean Fisher for
    the site... useful in the table for the publication"). Toujours
    renseigne ici (une moyenne archivee par cette fonction sait forcement
    ce qu'elle combine), au contraire d'une moyenne importee d'un autre
    format (voir FitResult.n_lines)."""
    vgp_lat, vgp_lon = dir_to_vgp(stats.dec, stats.inc, site_lat, site_lon)
    vgp_dp, vgp_dm = dp_dm_from_a95(stats.a95, stats.inc)
    codes = ":".join(str(r.c) for r in contributing if r.c)
    n_lines = sum(1 for r in contributing if r.cat1 in _LINE_LIKE_CAT1)
    n_planes = sum(1 for r in contributing if r.cat1 == "P")
    return FitResult(
        id=f"mean: {site}", cat1="F", cat2="i",
        nb=stats.n, dec=stats.dec, inc=stats.inc, mad=stats.a95,
        tx=(stats.k, 0.0), par3_mean=float(orientation),
        lat=site_lat, rlong=site_lon,
        par4=vgp_lat, par5=vgp_lon,
        vgp_dp=vgp_dp, vgp_dm=vgp_dm,
        component=component or "A",
        n_lines=n_lines, n_planes=n_planes,
        liste="codes:" + codes,
    )


def _mean_line_plane_counts(r: FitResult, results: List[FitResult]) -> Optional[Tuple[int, int]]:
    """Equivalent de "L nnn  P nnn" dans `lisres` : compte, parmi les
    resultats specimen PRESENTS dans `results` (meme liste passee a
    list_results) dont le `c` figure dans `r.liste`, combien sont des
    droites (cat1=='L') et des plans (cat1=='P'). Contrairement au
    Fortran (qui stocke ces 2 comptes directement dans tr(i).tx(1)/tx(2)
    a l'archivage), le format .pmagres ne les stocke pas - tx(1) y est
    deja pris par k (voir _format_mean_line) - donc on les retrouve ici
    par recoupement via `liste`. Retourne None si `results` ne contient
    aucun des specimens references (ex: `results` == uniquement les
    moyennes, chargees avec carselect='m') - dans ce cas le compte est
    inconnu plutot que faussement 0/0."""
    matched = mean_components(r, results)
    if not matched:
        return None
    nlines = sum(1 for m in matched if m.cat1 == "L")
    nplanes = sum(1 for m in matched if m.cat1 == "P")
    return nlines, nplanes


def mean_components(r: FitResult, results: List[FitResult]) -> List[FitResult]:
    """Retourne les resultats specimen (parmi `results`) references dans
    `r.liste` ("codes:c1:c2:...") - memes matching que
    _mean_line_plane_counts, mais renvoie les FitResult eux-memes (pas
    seulement un compte lignes/plans) - demande explicite utilisateur
    ("when the results is a mean, can you plot the stereo with the mean
    and individual results") : `data+interpretation` (afficher_visres,
    app.py) ignorait jusqu'ici toute ligne "mean:" (comme le Fortran
    d'origine), meme quand ses composants individuels etaient bien
    charges a cote (mode "s" de Select results..., voir
    _load_site_results) - desormais utilise pour tracer moyenne + ses
    composants ensemble sur un seul stereo (build_stereo_results_figure
    accepte deja une liste melant moyenne(s) et resultats individuels,
    voir "Stereo Results"). Liste vide si `r.liste` ne reference rien de
    present dans `results` (ex. `results` ne contient que des moyennes,
    carselect='m')."""
    ids = [tok for tok in r.liste.replace("codes:", "").split(":") if tok.strip()]
    if not ids:
        return []
    by_c = {str(other.c): other for other in results if other.id[:5] != "mean:"}
    return [by_c[i] for i in ids if i in by_c]


def _mean_site_strike_dip(r: FitResult, donnees) -> Optional[Tuple[float, float]]:
    """Equivalent de la boucle `do inech=1,nb_ech ... if(tr(i).id(7:12)==
    pmag(inech).id(1:6))` : strike/pendage du site, lus sur N'IMPORTE
    QUEL specimen de `donnees` appartenant a ce site (voir
    _find_site_specimen - le pendage est une propriete du site, identique
    pour tous ses specimens). Retourne None si `donnees` n'est pas fourni
    ou si aucun specimen du site n'y est trouve."""
    site = _mean_site_name(r.id)
    ech = _find_site_specimen(site, donnees)
    if ech is None:
        return None
    return ech.str_, ech.dip


def list_results(results: List[FitResult], orientation: int = 1, donnees=None) -> str:
    """Equivalent de `lisres` (dataselect.f:690-820) : tableau des
    ajustements. Comme le Fortran (bloc `select case (iorient)` lignes
    722-784), dec/inc AFFICHES sont recalcules a partir de la direction
    BRUTE stockee (repere echantillon d'origine, jamais modifiee) en
    appliquant la correction d'orientation COURANTE - la meme direction
    stockee s'affiche donc differemment selon iorient (échantillon/in
    situ/apres pendage), au lieu de rester figee dans le repere
    echantillon (voir recompute_fit_geometry pour le pre-requis : cin/
    caz/dip/str_ doivent avoir ete completes depuis le specimen d'origine
    pour un resultat charge depuis .pmagres, sinon cette correction est
    fausse).

    Ligne "mean:" (moyenne de site) : dec/inc affiches NE SONT PAS
    recalcules (voir _correct_dec_inc) - une moyenne est deja figee dans
    l'orientation `par3_mean`, affichee avec SA PROPRE etiquette
    (Sa)/(IS)/(TC) plutot que celle demandee en tete de tableau, pour ne
    pas laisser croire a une reprojection qui n'a pas lieu - demande
    explicite utilisateur sur des moyennes importees de MagIC, deja en
    In-Situ/apres pendage : "les moyennes de site ne se calculent pas en
    coordonnees echantillons mais en In situ et/ou bedding correction".
    Memes autres colonnes que `lisres`, y compris (contrairement a une
    version anterieure de cette fonction) :
    - e95 : calcule ici (`_e95`), pure fonction de a95/inclinaison, comme
      dans le Fortran - pas une donnee stockee.
    - "L nnn  P nnn" (compte lignes/plans composant la moyenne) : la valeur
      PERSISTEE (r.n_lines/r.n_planes, voir build_site_mean_result) si
      disponible - fiable quel que soit ce qui est charge par ailleurs.
      Sinon (moyenne importee d'un autre format, colonne reservee mais pas
      peuplee - voir FitResult.n_lines), repli sur l'ancien recoupement de
      `r.liste` (les `c` des specimens composants) CONTRE `results`
      lui-meme - fiable seulement si `results` contient aussi ces
      specimens (typiquement `load_results(..., carselect='s')`, qui les
      inclut expres pour cet usage) ; sinon affiche "?" plutot qu'un faux
      "0/0" (voir `_mean_line_plane_counts`).
    - strike/dip du site : lus sur un specimen de `donnees` (parametre
      optionnel, ex. self.donnees) dont l'id commence par le nom du
      site ; affiche "?" si `donnees` n'est pas fourni.
    Le Fortran stocke lignes/plans directement dans tr(i).tx(1)/tx(2) a
    l'archivage - le format .pmagres a repurpose ce slot pour k (voir
    _format_mean_line, k=tx[0]), d'ou le recoupement via `liste` ici au
    lieu d'une simple lecture."""
    tag = _ORIENT_MODE_TAG.get(orientation, "Sa")
    lines = [
        _ORIENT_HEADER.get(orientation, _ORIENT_HEADER[1]),
        f"{_HEADER_MARK}     Sample          comp  cat  orig  demag   step1  stepn   nb   dec    inc   mad   ({tag}){_HEADER_MARK}",
    ]
    for i, r in enumerate(results, start=1):
        dec, inc = _correct_dec_inc(r, orientation)
        if r.id[:5] == "mean:":
            # "2/5 (line/plane)" - demande explicite utilisateur ("?/?
            # 2/5 replace by 2/5 (line/plane)"), remplace l'ancien
            # "L  2  P  5" (illisible une fois le vrai bug de lecture
            # corrige - voir _format_mean_line).
            if r.n_lines >= 0:
                lp_txt = f"{r.n_lines}/{r.n_planes} (line/plane)"
            else:
                counts = _mean_line_plane_counts(r, results)
                lp_txt = f"{counts[0]}/{counts[1]} (line/plane)" if counts else "?/? (line/plane)"
            e95 = _e95(r.mad, inc)
            # dp/dm : la valeur ARCHIVEE (voir dp_dm_from_a95) si presente,
            # sinon calculee a la volee depuis a95/inc (meme repli que e95
            # ci-dessus) - couvre les moyennes archivees avant l'ajout de
            # ces colonnes (voir _parse_mean_line, vgp_dp/vgp_dm=0.0 par
            # defaut pour un ancien fichier).
            dp_show, dm_show = (
                (r.vgp_dp, r.vgp_dm) if (r.vgp_dp or r.vgp_dm) else dp_dm_from_a95(r.mad, inc)
            )
            strdip = _mean_site_strike_dip(r, donnees)
            strdip_txt = f"str={strdip[0]:5.1f} dip={strdip[1]:4.1f}" if strdip else "str=?  dip=?"
            # (Sa)/(IS)/(TC) AFFICHEE : celle de l'orientation ARCHIVEE
            # (r.par3_mean), PAS celle demandee en tete de tableau (voir
            # _correct_dec_inc - une moyenne n'est jamais reprojetee, la
            # valeur stockee est toujours celle affichee).
            own_tag = _ORIENT_MODE_TAG.get(int(r.par3_mean), "?")
            # Ligne separee pour [codes:...] (demande explicite utilisateur,
            # "when you list the results with a mean... a linefeed might be
            # useful") : `r.liste` porte un `c` par specimen combine dans la
            # moyenne, sans borne de longueur - accolee en fin de la ligne
            # de statistiques deja longue, elle rendait certaines moyennes
            # (site a beaucoup de specimens) illisibles sur une seule ligne
            # sans retour a la ligne, la fenetre Text ayant wrap="none".
            lines.append(
                f"{i:4d}: {r.id:<13s}[{r.component or 'A'}]    {r.cat1}{r.cat2}    {r.orig}     {r.demag:<3s}"
                f"  {lp_txt}  {r.nb:4d}  {dec:6.1f} {inc:6.1f} ({own_tag})  a95={r.mad:5.1f} e95={e95:5.1f}"
                f"  k={r.tx[0]:8.1f}  lat={r.lat:9.5f} lon={r.rlong:9.5f}"
                f"  VGP=({r.par4:6.1f},{r.par5:6.1f})  dp/dm=({dp_show:.1f}/{dm_show:.1f})  {strdip_txt}\n"
                f"          [{r.liste}]"
            )
        else:
            lines.append(
                f"{i:4d}: {r.id:<13s}[{r.component or 'A'}]    {r.cat1}{r.cat2}    {r.orig}     {r.demag:<3s}"
                f"  {r.step_first:5.0f}  {r.step_last:5.0f}  {r.nb:4d}  {dec:6.1f} {inc:6.1f}  {r.mad:5.1f}"
            )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# archivres / selres : sauvegarde et lecture du fichier ".r" (equivalent de
# fichiers.f:archivres et dataselect.f:selres). Format fixe (positions de
# colonnes exactes, pas un simple split) pour rester compatible avec les
# fichiers .r reels produits par Starmac_OSX - structure /Resultats/ complete
# (starmac_OSX.inc lignes 37-47) : c, id, cin, caz, dip, str, cat1, cat2,
# orig, demag, numcomp, nb, dec, inc, par1(mad), par2(step1), par3(stepn),
# tx(2), ty(2), tz(2). Le champ `liste` (cas special "mean:", agregats
# multi-echantillons) n'est pas gere ici - hors perimetre de ce portage.
# ---------------------------------------------------------------------------

def _fortran_e(value: float, width: int = 10, decimals: int = 3) -> str:
    """Equivalent d'un champ Fortran ExW.d (ex E10.3) : mantisse "0.ddd"
    (toujours <1 en valeur absolue), exposant signe sur 2 chiffres - PAS le
    style "d.ddE+xx" de Python/C (`%.3e`). Utilise par datatools.py (export
    detaille) - PLUS par le format .r (qui ne stocke plus tx/ty/tz, voir
    plus bas), garde ici comme utilitaire generique."""
    if value == 0.0:
        mantissa, exp = 0.0, 0
    else:
        sign = -1.0 if value < 0 else 1.0
        av = abs(value)
        exp = math.floor(math.log10(av)) + 1
        mantissa = av / (10.0 ** exp)
        if round(mantissa, decimals) >= 1.0:
            mantissa /= 10.0
            exp += 1
        mantissa *= sign
    sign_char = "-" if mantissa < 0 else " "
    digits = f"{abs(mantissa):.{decimals}f}"[2:]  # "0.123" -> "123"
    exp_sign = "+" if exp >= 0 else "-"
    text = f"{sign_char}0.{digits}E{exp_sign}{abs(exp):02d}"
    return text.rjust(width)[:width]


def results_path_for(data_path: str) -> str:
    """Equivalent de `filr=fil1(1:(nlen-4))//'.r'` : nom du fichier
    resultats derive du fichier de donnees charge (meme dossier, extension
    '.pmagres' - anciennement '.r', renomme a la demande explicite de
    l'utilisateur pour eviter toute ambiguite avec d'autres usages de
    '.r')."""
    base, _ext = os.path.splitext(data_path)
    return base + ".pmagres"


def ani_path_for(data_path: str) -> str:
    """Equivalent de results_path_for, pour le fichier d'anisotropie
    .pmagani (compagnon du fichier de donnees .prmag, MEME nom de base) -
    demande explicite utilisateur ("perhaps, we should have a new .ANI
    file that will be associated to the Prmag file") : formalise en une
    fonction UNIQUE, nommee, ce qui etait deja fait de facon dispersee/
    dupliquee dans app.py (4 sites derivant tour a tour de
    self.results_path ou de results_path_for applique a self.results_path,
    sans jamais passer par une fonction dediee comme results_path_for) -
    meme convention que results_path_for (extension REMPLACEE, pas
    simplement ajoutee au nom complet).

    Extension ".pmagani" (et non plus ".ANI") - demande explicite
    utilisateur ("can we write the name of the ani extension as
    .pmagani") : nom explicitement lie a .prmag/.pmagres, alors que
    ".ANI" restait ambigu (aussi utilise, format different, par
    AMS_OSX_AWE/AMS_Py). Les anciens fichiers ".ANI" restent lisibles
    directement (voir read_ani_tensor, dispatch par extension) ; leur
    CONVERSION explicite vers .pmagani se fait cote AMS_Py (ams_selection.
    import_legacy_ani) - demande explicite utilisateur ("in starmac_Py,
    you can delete the import legacy .ANI file, it is better to handle
    this in AMS_Py")."""
    base, _ext = os.path.splitext(data_path)
    return base + ".pmagani"


def pmagint_path_for(data_path: str) -> str:
    """Equivalent de results_path_for/ani_path_for, pour le fichier de
    RESULTATS de paleointensite .pmagint (compagnon du fichier de donnees
    .prmag, MEME nom de base) - demande explicite utilisateur ("pour
    garder la meme logique, je souhaite avoir un fichier supplementaire
    avec la ligne de resultats de paleointensite... qu'on pourrait
    appeler nomfichier.pmagint comme pmagres ou pmagani. il contiendrait
    les donnees comme dans tempint.txt"). tempint.txt (visi_Paleoint.f:
    40,1472,1487) est le nom fige, non associe au fichier de donnees, du
    fichier d'archivage des resultats de paleointensite que produisait
    `openfilepint` en Fortran (jamais porte en Python jusqu'ici - voir
    write_pmagint_line) ; .pmagint le remplace avec la convention de
    nommage deja etablie pour .pmagres/.pmagani."""
    base, _ext = os.path.splitext(data_path)
    return base + ".pmagint"


# ---------------------------------------------------------------------------
# Format .pmagres (anciennement .r, renomme a la demande de l'utilisateur)
# reorganise (demande explicite utilisateur, remplace l'ancien
# format Fortran a positions de colonnes fixes/E10.3) : deux sections
# separees dans le MEME fichier - "specimen results" (un ajustement par
# ligne) puis "site mean results" (une moyenne de site par ligne),
# tab-separees, cles auto-descriptives. Tel que specifie par l'utilisateur
# (fichier exemple_Pmag.r) :
# - PAS de tx/ty/tz (extremites du segment ajuste, pour le trace sur un
#   Zijderveld) : "to draw the line on the zijderveld plot we can redo
#   the calculation" - voir recompute_fit_geometry, appelee au chargement
#   plutot que stockee.
# - dec/inc restent en repere ECHANTILLON (pas de correction d'orientation
#   stockee : "as we have the orientation of the samples in the main
#   file, they may not be needed in the result file").
# - `random` (l'ancien `c:NNNNN`) EST conserve - PAS decoratif : sert de
#   cle de reference stable, un specimen pouvant porter plusieurs
#   ajustements candidats (ex. deux droites a des plages d'etapes
#   differentes) et la moyenne de site devant pouvoir dire PRECISEMENT
#   lesquels elle combine (`included_samples` = "codes:<random>:<random>:...").
# - cat2 abandonne (redondant - jamais filtre nulle part, cat1 suffit).
# - Insertion "au bon endroit" (pas toujours en fin de fichier comme
#   avant) : un nouveau resultat specimen est insere a la fin de la
#   SECTION 1 (juste avant la section "site mean", pas apres) ; une
#   nouvelle moyenne va en fin de fichier (section 2, toujours derniere).
#   Les lignes "!" (annotations manuelles de l'utilisateur, ex.
#   "! primary magnetization") et les lignes vides sont preservees telles
#   quelles - demande explicite ("In fortran I was writing at the end of
#   the file" -> ne plus faire ca aveuglement).
# ---------------------------------------------------------------------------

_SPECIMEN_HEADER = "#specimen results in sample coordinates"
_MEAN_HEADER = "#site mean results"

# (libelle, largeur) - espaces-alignees comme le fichier exemple fourni par
# l'utilisateur (exemple_Pmag.r), PAS des tabulations : plus proche de
# l'habitude de lecture Fortran a colonnes fixes (voir aussi .prmag). La
# derniere colonne (random / included_samples) n'est pas paddee (longueur
# variable, rien apres elle sur la ligne).
_SPECIMEN_FIELDS = [
    ("specimen", 12), ("step1", 7), ("step2", 7), ("L/P", 5), ("anc/not", 8),
    ("demag", 7), ("n", 6), ("dec", 7), ("inc", 7),
    ("mad", 7), ("component", 10), ("random", 8),
]
_MEAN_FIELDS = [
    ("site", 9), ("type", 6), ("n", 5), ("dec", 7), ("inc", 7), ("a95", 7),
    ("k", 9), ("IS/TC", 7), ("lat_site", 12), ("long_site", 12),
    ("VGP_lat", 8), ("VGP_lon", 8), ("VGP_dp", 7), ("VGP_dm", 7),
    ("component", 10), ("L/P", 8),
    ("included_samples", 0),
]
# "comp" (numcomp, 1-9 - numero d'ajustement PCA, sans aucun role dans le
# calcul - voir fit_line/fit_plane, purement une etiquette) RETIRE de
# _SPECIMEN_FIELDS/du fichier ecrit - demande explicite utilisateur
# ("autant partir sur un format plus etabli et retirer le numcomp"),
# suite a "there is an ambiguity between the component (number and
# letter)... it will be confusing for new user otherwise". "magcomp"
# renomme "component" (au lieu d'une abreviation cryptique) - plus besoin
# de le distinguer visuellement d'un "comp" numerique qui n'existe plus.
# `FitResult.numcomp` reste un champ INTERNE (fit_lines_auto : 1=ancre/
# 2=non-ancre, redo files, diagnostics) - simplement plus ecrit dans
# .pmagres. Retro-compatibilite en LECTURE sur 3 formats possibles (voir
# _parse_specimen_line) : le plus recent (12 colonnes, sans comp, avec
# component) est ecrit desormais ; deux anciens formats restent lisibles.

# "anc"/"not" (anchored/not-anchored) au lieu de 'o'/'n' (francais
# origine/non) dans la colonne "anc/not" du fichier .pmagres - demande
# explicite utilisateur. Le champ interne `FitResult.orig` garde 'o'/'n'
# (utilise ainsi partout ailleurs dans calcul.py/app.py/stereo.py) - seule
# la serialisation fichier change. La lecture reste retro-compatible avec
# d'anciens fichiers .pmagres deja ecrits avec 'o'/'n'.
_ANCHOR_TO_FILE_CODE = {"o": "anc", "n": "not"}
_FILE_CODE_TO_ANCHOR = {"anc": "o", "not": "n", "o": "o", "n": "n"}

# 0/100 (convention MagIC dir_tilt_correction : geographique/apres pendage
# complet) au lieu de 2.0/3.0 dans la colonne "IS/TC" du fichier .pmagres -
# demande explicite utilisateur ("for the mean, can we use 0 and 100 for
# IS and full TC like in magic"). Le champ interne `par3_mean` garde la
# convention STARpaleomag_Py existante (1/2/3, partagee avec `self.orientation`
# ailleurs dans l'appli, voir _ORIENT_MODE_TAG) - seule la serialisation
# fichier change, comme _ANCHOR_TO_FILE_CODE ci-dessus. "1" (coordonnees
# echantillon) n'est normalement plus produit pour une moyenne (voir
# extract_magic.magic_site_means) mais reste gere en lecture/ecriture pour
# ne pas casser d'anciens fichiers .pmagres.
_ORIENT_TO_FILE_CODE = {1.0: "1", 2.0: "0", 3.0: "100"}
_FILE_CODE_TO_ORIENT = {"1": 1.0, "0": 2.0, "100": 3.0}


def _parse_anchor_token(token: str) -> str:
    t = token.strip().lower()
    return _FILE_CODE_TO_ANCHOR.get(t, t[:1] if t else "n")


def _row(fields: List[Tuple[str, int]], values: List[str]) -> str:
    """BUG REEL corrige ici (trouve en testant import_published_means avec
    un nom de site de 9 caracteres, la largeur nominale de la colonne
    'site' - voir _MEAN_FIELDS) : `value.ljust(width)` ne separe PAS deux
    colonnes des que `value` atteint ou depasse `width` (ex. site="J0123456"
    a 8 caracteres, colonne 'site' large de 9 -> aucun padding ajoute,
    "J0123456" colle directement a la colonne suivante) - la ligne ecrite
    a alors un token de MOINS une fois relue par whitespace-split (meme
    categorie de bug que la colonne 'liste' vide - voir _format_mean_line).
    Garantit desormais au moins UN espace separateur apres chaque colonne
    a largeur fixe, meme si la valeur deborde de sa largeur nominale."""
    parts = []
    for (_label, width), value in zip(fields, values):
        if not width:
            parts.append(value)
        elif len(value) < width:
            parts.append(value.ljust(width))
        else:
            parts.append(value + " ")
    return "".join(parts).rstrip()


def _cols_header(fields: List[Tuple[str, int]]) -> str:
    """En-tete de colonnes, prefixee "#" - le "#" REMPLACE un caractere de
    remplissage (pas un ajout) pour que les colonnes suivantes restent
    alignees avec les lignes de donnees en dessous (sans prefixe) - meme
    convention que exemple_Pmag.r ("#specimen" = 9 caracteres, comme
    "14NQ0101A")."""
    labels = [label for label, _w in fields]
    labels[0] = "#" + labels[0]
    return _row(fields, labels)


def _fmt1(value: float, decimals: int = 1) -> str:
    return f"{float(value):.{decimals}f}"


def _mean_site_name(res_id: str) -> str:
    return res_id[6:] if res_id[:6].lower() == "mean: " else res_id


def _format_specimen_line(res: FitResult) -> str:
    return _row(_SPECIMEN_FIELDS, [
        res.id, _fmt1(res.step_first), _fmt1(res.step_last),
        res.cat1, _ANCHOR_TO_FILE_CODE.get(res.orig, res.orig), res.demag,
        str(res.nb),
        _fmt1(res.dec), _fmt1(res.inc), _fmt1(res.mad),
        (res.component or "A").strip(), str(res.c),
    ])


def _format_mean_line(res: FitResult) -> str:
    tilt_code = _ORIENT_TO_FILE_CODE.get(res.par3_mean, _fmt1(res.par3_mean, 0))
    lp_token = f"{res.n_lines}/{res.n_planes}" if res.n_lines >= 0 else "?/?"
    # BUG REEL corrige ici (signale par l'utilisateur, "?/?     2/5 replace
    # by 2/5" - le L/P se retrouvait faux/manquant a la relecture) : `liste`
    # (derniere colonne, largeur 0 - voir _MEAN_FIELDS) etait ecrite TELLE
    # QUELLE, y compris vide ("" si aucun `contributing` connu - ex. une
    # moyenne construite depuis une table externe, sans resultats
    # individuels a combiner). `_row` fait un `.rstrip()` final sur la
    # ligne entiere : une derniere colonne vide efface alors AUSSI le
    # remplissage (largeur 8) de `lp_token` juste avant - la ligne ecrite
    # a donc UN token de MOINS que prevu. `_parse_mean_line` (retro-
    # compatibilite sur le NOMBRE de colonnes) retombe alors sur la
    # branche a 16 colonnes au lieu de 17 : lp_token y est FIGE a "?/?"
    # (jamais relu depuis le fichier dans cette branche) et le VRAI
    # lp_token ecrit ("6/0" par ex.) est mal interprete comme `liste`.
    # Corrige a la source en ne laissant JAMAIS `liste` vide en ecriture -
    # "codes:" (prefixe vide, deja ce que mean_components() interprete
    # comme "aucun composant individuel connu") garantit une derniere
    # colonne toujours non-vide, donc jamais absorbee par le rstrip.
    liste = res.liste or "codes:"
    return _row(_MEAN_FIELDS, [
        _mean_site_name(res.id), "Fi", str(res.nb),
        _fmt1(res.dec), _fmt1(res.inc), _fmt1(res.mad), _fmt1(res.tx[0]),
        tilt_code, f"{res.lat:.5f}", f"{res.rlong:.5f}",
        _fmt1(res.par4), _fmt1(res.par5),
        _fmt1(res.vgp_dp), _fmt1(res.vgp_dm),
        (res.component or "A").strip(), lp_token, liste,
    ])


def _format_result_line(res: FitResult) -> str:
    """Equivalent de `archivres` : choisit le format specimen ou "mean:"
    (si `res.id[:5]=="mean:"`) - meme test que le Fortran d'origine."""
    if res.id[:5] == "mean:":
        return _format_mean_line(res)
    return _format_specimen_line(res)


_KNOWN_CAT1_LETTERS = ("L", "P", "f", "F", "s", "S")

# Variante encore plus ancienne de .r decouverte sur un fichier reel
# (pmagCARF.r, campagnes 94-98/03) : la colonne "anc/not" (index 4, voir
# _SPECIMEN_FIELDS) est carrement ABSENTE de la ligne quand le fit n'est
# PAS ancre (pas meme un espace-caractere, "o" est le SEUL cas ecrit) -
# demande explicite utilisateur ("I just need a file with specimen temp1
# temp2 ancr (o or white sspace) comp"). Ce fichier a en plus DEUX
# colonnes numeriques (correction de pendage, jamais utilisees ailleurs -
# le fichier .prmag porte deja bed_dip_strike) inserees avant la lettre
# L/P, elle-meme TOUJOURS collee sans espace au nombre qui la precede
# (ex. "11.2L", "1.3P" - artefact de largeur de colonne Fortran figee).
# Consequence : `.split()` decale TOUS les champs a partir de l'index 3
# d'une position variable (0, 1 ou 2 selon anc/not present ou non et la
# largeur des colonnes de pendage) - la lecture "moderne" ci-dessus
# retombait donc sur des valeurs n'importe quoi pour cat1/orig/demag
# (ex. cat1='64.1' au lieu de 'L'/'P'), silencieusement (aucune
# ValueError - ces chaines sont toutes des valeurs de type valide) : TOUS
# les resultats du fichier reel se retrouvaient avec un cat1 hors de
# ('L','P','f','F') et disparaissaient de `data+interpretation`
# (afficher_visres filtre sur cat1) sans aucun message d'erreur. Detecte
# ici en verifiant que parts[3] (cat1 attendu) est bien une des lettres
# connues - sinon, retrouve la position reelle de L/P en cherchant le
# token colle (`_OLD_VARIANT_LP_RE`) a partir de l'index 3."""
_OLD_VARIANT_LP_RE = re.compile(r"^-?\d+\.\d+([LP])$")
_OLD_VARIANT_DEMAG_LETTERS = set("DdSsFfNnIi")


def _component_from_numcomp(numcomp: int) -> str:
    """"A"/"B"/"C"... depuis un ANCIEN `numcomp` (1-9) - demande explicite
    utilisateur ("I think we should update everything to the letter and
    change the component number to letter during import of legacy") :
    les tres vieux formats de resultats (voir _parse_specimen_line_old_
    variant/_parse_legacy_result_line ci-dessous) n'ont JAMAIS de lettre
    A/B/C dediee - seul `numcomp` (1-9) y distingue deja, de facto, les
    differentes interpretations d'un meme specimen (numero d'ajustement
    PCA, seule information disponible a l'epoque) - AVANT l'ajout du
    champ `component` (demande explicite utilisateur separee, "when
    there is different components of magnetizations within the same
    site... We need to add a column component A,B,C"), tout import
    legacy figeait `component` a "A" par defaut, perdant cette
    distinction de fait au passage au nouveau format. UNIQUEMENT pour un
    resultat SPECIMEN (jamais une moyenne de site - deja explicitement
    exclu ailleurs : "I think it is best not to use the numcomp of
    individual samples for the mean"). 1->A, 2->B, ... 9->I ; toute
    valeur hors 1-9 retombe sur "A" plutot que produire un caractere
    non-alphabetique."""
    if 1 <= numcomp <= 9:
        return chr(ord("A") + numcomp - 1)
    return "A"


def _parse_specimen_line_old_variant(parts: List[str]) -> Optional[FitResult]:
    """Voir docstring de `_KNOWN_CAT1_LETTERS` ci-dessus pour le contexte -
    retrouve L/P par recherche du token colle plutot qu'une position fixe,
    puis avance un pointeur colonne par colonne (anc/not optionnel, demag
    optionnel) au lieu de supposer que toutes les colonnes sont presentes.
    `component` PAS trouvable comme telle dans ce tres vieux format (pas
    de lettre A/B/C/... nulle part dans le fichier reel de test - cette
    distinction n'existait probablement pas encore comme colonne dediee)
    - DERIVEE de `numcomp` (voir _component_from_numcomp) plutot que
    figee a "A" pour tous les resultats : demande explicite utilisateur
    ("change the component number to letter during import of legacy").
    `c` (juste un disambiguateur d'id, pas une donnee scientifique)
    reprend le token numerique restant s'il y en a un, sinon vide."""
    j = None
    for i in range(3, len(parts)):
        tok = parts[i]
        if tok in ("L", "P") or _OLD_VARIANT_LP_RE.match(tok):
            j = i
            break
    if j is None:
        return None
    lp = parts[j][-1]

    ptr = j + 1
    ancr = "n"
    if ptr < len(parts) and parts[ptr] == "o":
        ancr = "o"
        ptr += 1
    demag = ""
    if ptr < len(parts) and len(parts[ptr]) == 1 and parts[ptr] in _OLD_VARIANT_DEMAG_LETTERS:
        demag = parts[ptr]
        ptr += 1

    try:
        numcomp = int(parts[ptr]); ptr += 1
        nb = int(parts[ptr]); ptr += 1
        dec = float(parts[ptr]); ptr += 1
        inc = float(parts[ptr]); ptr += 1
        mad = float(parts[ptr]); ptr += 1
        return FitResult(
            id=parts[0].strip(),
            step_first=int(round(float(parts[1]))), step_last=int(round(float(parts[2]))),
            cat1=lp, orig=ancr, demag=demag, numcomp=numcomp, nb=nb,
            dec=dec, inc=inc, mad=mad,
            component=_component_from_numcomp(numcomp),
            c=(parts[ptr].strip() if ptr < len(parts) else ""),
        )
    except (ValueError, IndexError):
        return None


def _parse_specimen_line(parts: List[str]) -> Optional[FitResult]:
    if len(parts) < 12:
        return None
    if parts[3].strip() not in _KNOWN_CAT1_LETTERS:
        return _parse_specimen_line_old_variant(parts)
    try:
        # 3 formats possibles, du plus recent au plus ancien (voir
        # _SPECIMEN_FIELDS - "comp"/numcomp retire, "magcomp" renomme
        # "component" - demande explicite utilisateur "autant partir sur
        # un format plus etabli et retirer le numcomp") :
        #   - 13 colonnes : comp ET magcomp presents (format juste avant
        #     ce changement).
        #   - 12 colonnes, RECENT : sans comp, avec component - parts[10]
        #     (juste avant `c`) est une LETTRE.
        #   - 12 colonnes, TRES ANCIEN : avec comp, sans component/
        #     magcomp du tout - parts[10] est mad, un NOMBRE.
        # Les deux variantes a 12 colonnes ne se distinguent PAS par leur
        # longueur (identique) mais par le CONTENU de parts[10] - meme
        # principe que _parse_mean_line (verifie le nombre de colonnes),
        # pousse un cran plus loin ici car la longueur seule ne suffit
        # plus a lever l'ambiguite.
        if len(parts) >= 13:
            numcomp, nb_s = int(parts[6]), parts[7]
            dec_s, inc_s, mad_s = parts[8], parts[9], parts[10]
            component, c = parts[11].strip() or "A", parts[12].strip()
        else:
            try:
                float(parts[10])
                ancient_with_comp = True
            except ValueError:
                ancient_with_comp = False
            if ancient_with_comp:
                numcomp, nb_s = int(parts[6]), parts[7]
                dec_s, inc_s, mad_s = parts[8], parts[9], parts[10]
                component, c = "A", parts[11].strip()
            else:
                numcomp, nb_s = 1, parts[6]
                dec_s, inc_s, mad_s = parts[7], parts[8], parts[9]
                component, c = parts[10].strip() or "A", parts[11].strip()
        return FitResult(
            id=parts[0].strip(),
            step_first=int(round(float(parts[1]))), step_last=int(round(float(parts[2]))),
            cat1=parts[3].strip(), orig=_parse_anchor_token(parts[4]), demag=parts[5].strip(),
            numcomp=numcomp, nb=int(nb_s),
            dec=float(dec_s), inc=float(inc_s), mad=float(mad_s),
            component=component, c=c,
        )
    except ValueError:
        return None


def _parse_mean_line(parts: List[str]) -> Optional[FitResult]:
    if len(parts) < 13:
        return None
    try:
        tilt_token = parts[7].strip()
        orientation = _FILE_CODE_TO_ORIENT.get(tilt_token)
        if orientation is None:
            orientation = float(tilt_token)  # retro-compat anciens fichiers (1.0/2.0/3.0)
        # VGP_dp/VGP_dm ajoutes apres coup (demande explicite utilisateur) -
        # colonnes optionnelles inserees AVANT `liste` : un fichier .pmagres
        # ecrit avant ce changement n'a que 13 colonnes (liste en dernier),
        # un fichier ecrit apres en a 15 (dp/dm avant liste) - retro-
        # compatible en verifiant le nombre de colonnes plutot qu'un index
        # fixe pour `liste`.
        # "magcomp" (component A/B/C) insere apres vgp_dp/vgp_dm, "L/P"
        # (nb lignes/plans combines, ex. "10/2" - "?/?" si non renseigne)
        # insere apres magcomp, avant `liste` (toujours en dernier) -
        # demandes explicites utilisateur successives. 4 longueurs
        # possibles selon l'anciennete du fichier : 13 (avant vgp_dp/
        # vgp_dm), 15 (vgp_dp/vgp_dm), 16 (+ magcomp), 17 (+ L/P) -
        # verifie la plus longue en premier.
        if len(parts) >= 17:
            vgp_dp, vgp_dm = float(parts[12]), float(parts[13])
            component = parts[14].strip() or "A"
            lp_token, liste = parts[15].strip(), parts[16].strip()
        elif len(parts) >= 16:
            vgp_dp, vgp_dm = float(parts[12]), float(parts[13])
            component, liste = parts[14].strip() or "A", parts[15].strip()
            lp_token = "?/?"
        elif len(parts) >= 15:
            vgp_dp, vgp_dm, liste = float(parts[12]), float(parts[13]), parts[14].strip()
            component = "A"
            lp_token = "?/?"
        else:
            vgp_dp, vgp_dm, liste = 0.0, 0.0, parts[12].strip()
            component = "A"
            lp_token = "?/?"
        n_lines, n_planes = -1, -1
        if lp_token and lp_token != "?/?" and "/" in lp_token:
            try:
                nl_s, np_s = lp_token.split("/", 1)
                n_lines, n_planes = int(nl_s), int(np_s)
            except ValueError:
                n_lines, n_planes = -1, -1
        return FitResult(
            id="mean: " + parts[0].strip(), cat1="F", cat2="i",
            nb=int(parts[2]), dec=float(parts[3]), inc=float(parts[4]), mad=float(parts[5]),
            tx=(float(parts[6]), 0.0), par3_mean=orientation,
            lat=float(parts[8]), rlong=float(parts[9]),
            par4=float(parts[10]), par5=float(parts[11]),
            vgp_dp=vgp_dp, vgp_dm=vgp_dm, component=component,
            n_lines=n_lines, n_planes=n_planes, liste=liste,
        )
    except ValueError:
        return None


def _iter_result_lines(path: str):
    """Lit les deux sections du fichier .r - bascule de la section
    "specimen" vers "mean" des que la ligne d'en-tete "#site mean results"
    est rencontree. Les lignes vides et les annotations "!..." (ajoutees a
    la main par l'utilisateur) sont silencieusement ignorees ici (mais
    preservees telles quelles par archivres, qui ne fait jamais table
    rase du fichier)."""
    if not os.path.exists(path):
        return
    in_mean_section = False
    with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("!"):
                continue
            if stripped.startswith("#"):
                if "site mean" in stripped.lower():
                    in_mean_section = True
                continue
            parts = stripped.split()
            res = _parse_mean_line(parts) if in_mean_section else _parse_specimen_line(parts)
            if res is not None:
                yield res


def _next_specimen_c(specimen_id: str, existing_ids: set) -> str:
    """Identifiant "<specimen>_<lettre>" (a, b, c... puis aa, ab... au-dela
    de 26 interpretations du meme specimen, jamais observe en pratique) -
    REMPLACE un entier aleatoire 0-99999 (l'ancien schema) suite a un
    fichier reel corrompu (Tibet_14_15_Pmag_.r, deux campagnes de terrain
    fusionnees dans un seul .r) : 91 valeurs "c" dupliquees sur 1213,
    certaines associant meme un resultat specimen et une moyenne "mean:"
    totalement sans rapport - demande explicite utilisateur ("the most
    simple might be specimen_a, specimen_b instead of a random number as
    it is unlikely to have so many interpretations by specimen"). Un
    entier purement aleatoire peut coincider entre deux specimens SANS
    RAPPORT des que deux fichiers generes independamment sont fusionnes ;
    cet identifiant porte deja le nom du specimen - une collision
    necessiterait desormais le MEME specimen ET la MEME lettre (une
    veritable duplication de la meme interpretation, facilement reperee),
    pas un accident purement numerique."""
    for letter in string.ascii_lowercase:
        candidate = f"{specimen_id}_{letter}"
        if candidate not in existing_ids:
            return candidate
    for l1 in string.ascii_lowercase:
        for l2 in string.ascii_lowercase:
            candidate = f"{specimen_id}_{l1}{l2}"
            if candidate not in existing_ids:
                return candidate
    raise RuntimeError(f"too many interpretations stored for specimen {specimen_id}")


def archivres(
    res: FitResult, path: str, existing_ids: Optional[set] = None
) -> Tuple[str, set]:
    """Equivalent de `archivres`, mais insere CHAQUE resultat dans la
    bonne section plutot que de toujours ajouter en fin de fichier brut
    (demande explicite utilisateur : "can we split in two the file and
    write at the right place. In fortran I was writing at the end of the
    file") : un resultat specimen va a la fin de la section 1 (juste
    avant l'en-tete "#site mean results" s'il existe deja, en sautant les
    lignes vides qui la precedent pour ne pas casser l'espacement) ; une
    moyenne de site va toujours en fin de fichier (section 2, creee au
    besoin avec son propre en-tete). Genere un id anti-collision
    "<specimen>_<lettre>" (voir _next_specimen_c - REMPLACE l'ancien
    entier aleatoire 0-99999) - `existing_ids` : ensemble des `c`/`random`
    deja connus (evite de rouvrir le fichier a chaque appel ; passer None
    la premiere fois puis reutiliser l'ensemble retourne). Retourne `(c
    attribue, existing_ids mis a jour)`."""
    if existing_ids is None:
        existing_ids = {r.c for r in _iter_result_lines(path)}

    c = _next_specimen_c(res.id, existing_ids)
    existing_ids.add(c)
    res.c = c

    is_mean = res.id[:5] == "mean:"
    new_line = _format_result_line(res)

    lines: List[str] = []
    if os.path.exists(path):
        with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
            lines = f.read().splitlines()

    mean_header_idx = next(
        (i for i, l in enumerate(lines) if "site mean" in l.strip().lower() and l.strip().startswith("#")),
        None,
    )

    if is_mean:
        if mean_header_idx is None:
            if lines and lines[-1].strip():
                lines.append("")
            lines.append(_MEAN_HEADER)
            lines.append(_cols_header(_MEAN_FIELDS))
        lines.append(new_line)
    else:
        if not lines:
            lines = [_SPECIMEN_HEADER, _cols_header(_SPECIMEN_FIELDS)]
        if mean_header_idx is None:
            lines.append(new_line)
        else:
            insert_idx = mean_header_idx
            while insert_idx > 0 and not lines[insert_idx - 1].strip():
                insert_idx -= 1
            lines.insert(insert_idx, new_line)

    with open(path, "w", encoding="iso-8859-1", errors="replace", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
        # flush explicite (buffer Python) + fsync (cache OS -> disque) - demande
        # explicite utilisateur ("is it possible to flush the buffer after the
        # result is written, to be sure that the text is written") : chaque
        # resultat est archive individuellement au fil de la session (voir
        # docstring ci-dessus), un crash/force-quit juste apres un calcul ne
        # doit pas pouvoir perdre ce resultat parce qu'il ne serait encore
        # que dans un buffer (Python et/ou OS), pas sur le disque.
        f.flush()
        os.fsync(f.fileno())
    return c, existing_ids


def import_published_means(
    pmagres_path: str, table_path: str, encoding: str = "utf-8",
) -> Tuple[int, List[str]]:
    """Archive des moyennes de site PUBLIEES (table externe, ex. Table 1
    d'un papier) comme resultats "mean:" - demande explicite utilisateur
    ("how do we fill the prmag in that case", suite a "can we have a site
    in prmag without data?") : pour un site dont les mesures BRUTES sont
    perdues (fichier legacy), le SEUL moyen d'obtenir un resultat "mean:"
    exploitable (par data+interpretation, l'export MagIC - voir
    magic_export.build_sites_rows - etc.) etait jusqu'ici soit de calculer
    une moyenne depuis des specimens reellement charges (build_site_mean_
    result), soit un import MagIC deja publie, soit un script ad hoc
    (voir la conversion manuelle faite cette session pour Table 1
    d'Arriagada et al. 2006) - generalise ici en fonctionnalite reutilisable.

    Table (colonnes obligatoires 'site','orientation','dec','inc','a95',
    'k','n' - meme lecteur que complete_sample_info._read_table_rows,
    tabulation OU virgule auto-detectee) :
        #site   orientation   dec   inc   a95   k   n
        J01     IS            46.1  -55.4 7.5   80  6
    'orientation' : IS/TC (insensible a la casse, ou directement 0/100 -
    meme convention fichier que _FILE_CODE_TO_ORIENT). Colonnes
    optionnelles : 'n_lines','n_planes' (colonne L/P du .pmagres, vide si
    omises - voir _format_mean_line) ; 'component' (A/B/C..., 'A' par
    defaut) ; 'lat','lon' (site, necessaires au calcul du VGP - VGP
    calcule UNIQUEMENT si l'orientation est TC, meme raison que partout
    ailleurs dans ce module : un VGP n'a de sens que depuis une direction
    deja en repere geographique/apres correction de pendage, jamais
    depuis de l'in-situ).

    Retourne (nb_moyennes_archivees, liste_erreurs_par_ligne) - une ligne
    en erreur (orientation non reconnue, dec/inc/a95/k/n non numeriques)
    est SIGNALEE et IGNOREE, n'interrompt jamais l'import des autres
    lignes."""
    from complete_sample_info import _read_table_rows

    header, rows = _read_table_rows(table_path, encoding=encoding)
    required = ("site", "orientation", "dec", "inc", "a95", "k", "n")
    missing = [c for c in required if c not in header]
    if missing:
        raise ValueError(
            f"Table file must have columns: {', '.join(required)} (missing: {', '.join(missing)})."
        )

    orient_aliases = {
        "is": 2.0, "insitu": 2.0, "in-situ": 2.0, "in_situ": 2.0, "0": 2.0, "2": 2.0,
        "tc": 3.0, "tilt": 3.0, "tiltcorrected": 3.0, "tilt-corrected": 3.0,
        "tilt_corrected": 3.0, "100": 3.0, "3": 3.0,
    }

    n_imported = 0
    errors: List[str] = []
    existing_ids: Optional[set] = None
    for i, row in enumerate(rows, start=2):
        record = dict(zip(header, row))

        def col(name: str, default: str = "") -> str:
            return record.get(name, default).strip()

        site = col("site")
        if not site:
            continue
        orient_key = col("orientation").lower().replace(" ", "").replace("_", "").replace("-", "")
        orientation = orient_aliases.get(orient_key)
        if orientation is None:
            errors.append(f"line {i} (site {site}): unrecognized orientation "
                           f"'{col('orientation')}' (use IS or TC)")
            continue
        try:
            dec, inc = float(col("dec")), float(col("inc"))
            a95, k = float(col("a95")), float(col("k"))
            n = int(float(col("n")))
        except ValueError:
            errors.append(f"line {i} (site {site}): dec/inc/a95/k/n must be numeric")
            continue

        n_lines = int(float(col("n_lines"))) if col("n_lines") else -1
        n_planes = int(float(col("n_planes"))) if col("n_planes") else -1
        component = col("component") or "A"
        try:
            lat = float(col("lat")) if col("lat") else 0.0
            lon = float(col("lon")) if col("lon") else 0.0
        except ValueError:
            lat = lon = 0.0

        vgp_lat = vgp_lon = vgp_dp = vgp_dm = 0.0
        if orientation == 3.0 and (lat or lon):
            vgp_lat, vgp_lon = dir_to_vgp(dec, inc, lat, lon)
            vgp_dp, vgp_dm = dp_dm_from_a95(a95, inc)

        res = FitResult(
            id=f"mean: {site}", cat1="F", cat2="i",
            dec=dec, inc=inc, mad=a95, tx=(k, 0.0), nb=n,
            par3_mean=orientation, component=component,
            lat=lat, rlong=lon, par4=vgp_lat, par5=vgp_lon,
            vgp_dp=vgp_dp, vgp_dm=vgp_dm, n_lines=n_lines, n_planes=n_planes,
        )
        _, existing_ids = archivres(res, pmagres_path, existing_ids)
        n_imported += 1

    return n_imported, errors


def recompute_fit_geometry(res: FitResult, donnees) -> FitResult:
    """Complete un FitResult CHARGE DEPUIS .pmagres avec ce que ce format
    ne stocke plus, en le recalculant/le relisant depuis le specimen
    d'origine (`donnees`, ex. self.donnees) :

    1) cin/caz/dip/str_ (orientation du specimen - core+pendage) : le
       nouveau format .pmagres ne les stocke PAS ("as we have the
       orientation of the samples in the main file, they may not be
       needed in the result file"), donc un FitResult rechargee les a
       tous a 0.0 (valeurs par defaut du dataclass). BUG REEL corrige ici
       (signale par l'utilisateur : "the core correction is not correct
       for the results in the stereo... it seems that there is a -90
       applied instead of the core correction") : `corfor(x,y,z,cin=0,
       caz=0)` ne renvoie PAS l'identite, c'est une vraie rotation dont
       le residu ressemble a un decalage de 90 - `_correct_dec_inc`
       (calcul.py/stereo.py) appliquait donc une fausse correction des
       qu'un resultat charge depuis le fichier etait affiche en
       orientation 2/3 (in-situ/pendage), alors que le Zijderveld (qui
       lit cin/caz directement sur le SelectedSample vivant, jamais sur
       le FitResult) n'etait pas touche - d'ou l'ecart stereo/Zijderveld
       constate. Applique a TOUT resultat specimen (pas seulement L/P).
    2) tx/ty/tz (segment ajuste, pour le trace sur un Zijderveld) : pas
       stockes non plus ("to draw the line on the zijderveld plot we can
       redo the calculation") - recalcule uniquement pour L/P (seuls
       types avec un segment a tracer).

    Les resultats "mean:" (moyennes de site, pas un seul specimen) sont
    retournes inchanges. Si le specimen est introuvable dans `donnees`,
    les champs concernes restent a leur defaut (comportement inchange)."""
    if res.id[:5] == "mean:":
        return res
    ech = next((s for s in donnees if s.id == res.id), None)
    if ech is None:
        return res

    res = replace(res, cin=ech.cin, caz=ech.caz, dip=ech.dip, str_=ech.str_)

    if res.cat1 not in ("L", "P") or not getattr(ech, "mesures", None):
        return res
    jdeb = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == res.step_first), None)
    jfin = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == res.step_last), None)
    if jdeb is None or jfin is None or jfin < jdeb:
        return res

    anchored = res.orig == "o"
    if res.cat1 == "L":
        recomputed = fit_line(ech, jdeb, jfin, anchored=anchored, numcomp=res.numcomp)
    else:
        recomputed = fit_plane(ech, jdeb, jfin, normalize=True, numcomp=res.numcomp)
    if recomputed is None:
        return res
    return replace(res, tx=recomputed.tx, ty=recomputed.ty, tz=recomputed.tz)


# ---------------------------------------------------------------------------
# Conversion d'un ANCIEN fichier .r (positions de colonnes fixes, format
# Fortran d'origine) vers le nouveau .pmagres - demande explicite
# utilisateur ("during the convert legacy file to the new format, is it
# possible to convert the .r file with the results"), typiquement appelee
# en meme temps que la conversion .ren -> .prmag (voir convert_ren_to_r.py).
# Positions figees ici UNIQUEMENT pour lire un vieux fichier existant - le
# lecteur "en direct" (_iter_result_lines) n'utilise plus que le nouveau
# format.
# ---------------------------------------------------------------------------

_LEGACY_R_COLS = {
    "prefix": (0, 2), "c": (2, 7), "id": (9, 21),
    "cin": (23, 28), "caz": (29, 34), "dip": (35, 40), "str_": (41, 46),
    "cat1": (48, 49), "cat2": (49, 50), "orig": (53, 54), "demag": (59, 60),
    "numcomp": (65, 66), "nb": (68, 71), "dec": (72, 77), "inc": (78, 83),
    "par1": (85, 89), "par2": (91, 97), "par3": (98, 104),
    "tx1": (105, 115), "tx2": (116, 126), "ty1": (127, 137),
    "ty2": (138, 148), "tz1": (149, 159), "tz2": (160, 170),
}
_LEGACY_R_COLS_MEAN = {
    "tx1": (105, 115), "tx2": (116, 126),
    "lat": (128, 138), "rlong": (140, 150),
    "par4": (152, 158), "par5": (160, 166),
    "liste": (169, None),
}


def _parse_legacy_result_line(line: str) -> Optional[FitResult]:
    """Un champ "orig" (anchored/not, colonne 53:54) vide/blanc dans
    l'ancien fichier .r est traite comme "n" (non ancre) - demande
    explicite utilisateur ("during the import of legacy results, the
    missing information for anchored or not, means not")."""
    if not line.startswith("c:") or len(line) < _LEGACY_R_COLS["par3"][1]:
        return None

    def seg(key, cols=_LEGACY_R_COLS):
        a, b = cols[key]
        return line[a:b if b is not None else len(line)]

    try:
        prefix_kwargs = dict(
            c=seg("c").strip(), id=seg("id").strip(),
            cin=float(seg("cin")), caz=float(seg("caz")),
            dip=float(seg("dip")), str_=float(seg("str_")),
            cat1=seg("cat1"), cat2=seg("cat2"),
            orig=seg("orig").strip() or "n",
            demag=seg("demag").strip(), numcomp=int(seg("numcomp")),
            nb=int(seg("nb")), dec=float(seg("dec")), inc=float(seg("inc")),
            mad=float(seg("par1")),
        )
        par2_raw, par3_raw = float(seg("par2")), float(seg("par3"))
    except ValueError:
        return None

    if prefix_kwargs["id"][:5] == "mean:":
        try:
            return FitResult(
                **prefix_kwargs, par2_mean=par2_raw, par3_mean=par3_raw,
                tx=(float(seg("tx1", _LEGACY_R_COLS_MEAN)), float(seg("tx2", _LEGACY_R_COLS_MEAN))),
                lat=float(seg("lat", _LEGACY_R_COLS_MEAN)), rlong=float(seg("rlong", _LEGACY_R_COLS_MEAN)),
                par4=float(seg("par4", _LEGACY_R_COLS_MEAN)), par5=float(seg("par5", _LEGACY_R_COLS_MEAN)),
                liste=seg("liste", _LEGACY_R_COLS_MEAN).strip(),
            )
        except ValueError:
            return None

    prefix_kwargs["step_first"] = int(round(par2_raw))
    prefix_kwargs["step_last"] = int(round(par3_raw))
    try:
        return FitResult(
            **prefix_kwargs,
            # component DERIVEE de numcomp - specimen SEULEMENT, jamais la
            # branche "mean:" ci-dessus (voir _component_from_numcomp -
            # demande explicite utilisateur "change the component number
            # to letter during import of legacy").
            component=_component_from_numcomp(prefix_kwargs["numcomp"]),
            tx=(float(seg("tx1")), float(seg("tx2"))),
            ty=(float(seg("ty1")), float(seg("ty2"))),
            tz=(float(seg("tz1")), float(seg("tz2"))),
        )
    except ValueError:
        return None


def convert_legacy_results_file(old_path: str, new_path: str) -> Tuple[int, List[str]]:
    """Convertit un ANCIEN fichier .r (colonnes fixes) vers le nouveau
    .pmagres (sections specimen/site mean, colonnes auto-descriptives).
    tx/ty/tz de l'ancien fichier sont IGNORES (pas ecrits dans le nouveau
    format - voir recompute_fit_geometry, qui les recalcule a la demande).

    RENUMEROTE au passage chaque `c`/"random" en "<specimen>_<lettre>"
    (voir _next_specimen_c) plutot que de recopier l'ancien entier
    0-99999 tel quel - demande explicite utilisateur, suite a un fichier
    reel corrompu (Tibet_14_15_Pmag_.r, deux campagnes de terrain
    fusionnees dans un seul .r ancien format) : l'ancien entier ALEATOIRE
    entrait en collision entre deux specimens SANS RAPPORT, corrompant
    silencieusement les "codes:" des moyennes de site ("can you check the
    random number... perhaps rename it site_rannum... the most simple
    might be specimen_a, specimen_b instead of a random number... can you
    implement this in the import legacy file").

    La reconstruction des "codes:" d'une moyenne relit chaque ANCIEN
    numero de sa liste et le retrouve dans une table {ancien numero ->
    nouvel id}, mise a jour au fur et a mesure que les lignes specimen
    sont lues - PAS une fenetre positionnelle "les N dernieres lignes
    lues" (un premier essai de cette fonction faisait cette hypothese ;
    fausse sur donnees reelles : un site peut lister TOUS ses resultats
    specimen de PLUSIEURS composantes d'affilee, puis PLUSIEURS paires de
    moyennes a la fin - chacune resumant un SOUS-ENSEMBLE different, pas
    forcement le plus recent - verifie sur le site "14NQ04" du fichier
    reel : 22 resultats specimen (2 groupes de 11) precedent 4 lignes
    "mean:", les 2 premieres resumant le 1er groupe, les 2 suivantes le
    2e). La table est indexee par l'ANCIEN numero justement PARCE QU'il
    reste unique LOCALEMENT (les collisions reelles n'existent qu'ENTRE
    sites distants du fichier, jamais entre deux resultats d'un meme
    site - verifie sur le fichier reel) - une entree peut donc etre
    ecrasee par un doublon plus loin dans le fichier SANS consequence, la
    moyenne concernee ayant deja consomme la bonne valeur avant que ce
    doublon distant n'apparaisse. Un ancien numero reference par une
    moyenne mais introuvable dans la table (jamais vu, ou deja ecrase par
    un doublon survenu ENTRE temps - non observe sur le fichier reel mais
    possible en theorie) est signale dans les avertissements retournes
    plutot que de deviner une correspondance fausse.

    Retourne (nombre de resultats convertis, avertissements) - (0, [])
    si `old_path` n'existe pas ou ne contient aucune ligne reconnue
    (`new_path` n'est alors pas cree)."""
    if not os.path.exists(old_path):
        return 0, []

    letter_counts: Dict[str, int] = {}
    old_c_to_new: Dict[str, str] = {}
    warnings: List[str] = []
    specimen_lines, mean_lines = [], []
    with open(old_path, "r", encoding="iso-8859-1", errors="replace") as f:
        for raw in f:
            res = _parse_legacy_result_line(raw.rstrip("\n"))
            if res is None:
                continue
            if res.id[:5] == "mean:":
                old_codes = [t for t in res.liste.replace("codes:", "").split(":") if t.strip()]
                new_codes = []
                missing = []
                for old_c in old_codes:
                    new_c = old_c_to_new.get(old_c)
                    if new_c is None:
                        # Repli OBSERVE sur un fichier reel (Tibet_14_15_Pmag_.r,
                        # sites 14NQ13/14NQ15/14NQ23) : un "6" de tete manquant
                        # de facon DETERMINISTE (ex. "8624" au lieu de "68624") -
                        # confirme sur 3 cas independants, dont 2 ou le meme
                        # specimen apparait CORRECTEMENT (avec le "6") dans une
                        # AUTRE moyenne du meme fichier - probable troncature
                        # d'un champ largeur fixe cote Fortran d'origine, pas
                        # une simple coincidence. Recupere ce cas precis (SEUL
                        # candidat "6"+old_c connu) plutot que de le signaler
                        # manquant, mais le note explicitement (pas silencieux)
                        # pour verification par l'utilisateur.
                        recovered = old_c_to_new.get("6" + old_c)
                        if recovered is not None:
                            new_codes.append(recovered)
                            warnings.append(
                                f"{res.id}: old id '{old_c}' not found as-is, recovered as "
                                f"'6{old_c}' -> {recovered} (likely dropped leading digit "
                                f"'6' - please verify)."
                            )
                        else:
                            missing.append(old_c)
                    else:
                        new_codes.append(new_c)
                if missing:
                    warnings.append(
                        f"{res.id}: {len(missing)} of {len(old_codes)} referenced specimen "
                        f"result(s) not found (old id(s) {missing}) - codes list left incomplete."
                    )
                res.liste = "codes:" + ":".join(new_codes)
                mean_lines.append(_format_mean_line(res))
            else:
                old_c = res.c
                n = letter_counts.get(res.id, 0)
                letter = string.ascii_lowercase[n] if n < 26 else f"{string.ascii_lowercase[n // 26 - 1]}{string.ascii_lowercase[n % 26]}"
                letter_counts[res.id] = n + 1
                res.c = f"{res.id}_{letter}"
                old_c_to_new[old_c] = res.c
                specimen_lines.append(_format_specimen_line(res))

    if not specimen_lines and not mean_lines:
        return 0, warnings

    out_lines: List[str] = []
    if specimen_lines:
        out_lines += [_SPECIMEN_HEADER, _cols_header(_SPECIMEN_FIELDS), *specimen_lines]
    if mean_lines:
        if out_lines:
            out_lines.append("")
        out_lines += [_MEAN_HEADER, _cols_header(_MEAN_FIELDS), *mean_lines]

    with open(new_path, "w", encoding="iso-8859-1", errors="replace", newline="\n") as f:
        f.write("\n".join(out_lines) + "\n")
        f.flush()
        os.fsync(f.fileno())
    return len(specimen_lines) + len(mean_lines), warnings


def _load_data_results(
    path: str, pattern: str, cat1: str, numcomp: Optional[int], component: str = "*"
) -> List[FitResult]:
    """Equivalent de la branche `carselect=="d"` (par defaut) de `selres` :
    resultats normaux (L/P/f/s...), les lignes "mean:" sont exclues -
    fidele au Fortran, ou une ligne mean ne matche NI la branche mean
    (carselect!='m') NI la branche normale (verifie `chaineres(10:15).ne.
    "mean: "`). `component` : filtre PAS dans le Fortran (voir
    FitResult.component) - "*" (defaut) = pas de filtre."""
    pattern = (pattern or "*").upper().strip()
    nlen = len(pattern) if pattern != "*" else 0
    component = (component or "*").strip()
    results = []
    for res in _iter_result_lines(path):
        if res.id[:5] == "mean:":
            continue
        if pattern != "*" and res.id[:nlen].upper() != pattern[:nlen]:
            continue
        if cat1 != "*" and res.cat1 != cat1:
            continue
        if numcomp is not None and res.numcomp != numcomp:
            continue
        if component != "*" and (res.component or "A").upper() != component.upper():
            continue
        results.append(res)
    return results


def _load_mean_results(
    path: str, pattern: str, iorient: int, component: str = "*"
) -> List[FitResult]:
    """Equivalent de la branche `carselect=="m"` de `selres` : UNIQUEMENT les
    moyennes de site ("mean: <site>"), filtrees par `res.par3==float(iorient)`
    (le code d'orientation enregistre avec la moyenne) et par `pattern`
    (compare apres le prefixe "mean: ", ajoute automatiquement - equivalent
    de `enteteres="mean: "//entete` puis `numero=enteteres//chaine`).
    `component` : filtre PAS dans le Fortran (voir FitResult.component) -
    "*" (defaut) = pas de filtre. Compare la propre etiquette de LA
    MOYENNE, jamais celle de ses resultats constitutifs ("I think it is
    best not to use the numcomp de samples individuels pour la moyenne").

    Filtre d'orientation EXACT (`res.par3_mean == iorient`) - demande
    explicite utilisateur ("If you are in in situ, you select only the
    mean that are in in situ. If you want to select in TC, the user
    should decide to switch to TC") : le comportement voulu n'est PAS une
    reprojection automatique mais un filtrage strict - basculer
    l'orientation courante de l'appli (self.orientation) est la SEULE
    facon de rendre selectionnable une moyenne archivee dans une autre
    orientation."""
    site = (pattern or "*").strip()
    full_pattern = "*" if site in ("", "*") else f"MEAN: {site.upper()}"
    nlen = len(full_pattern) if full_pattern != "*" else 0
    component = (component or "*").strip()

    results = []
    for res in _iter_result_lines(path):
        if res.id[:5] != "mean:":
            continue
        if res.par3_mean != float(iorient):
            continue
        if full_pattern != "*" and res.id[:nlen].upper() != full_pattern[:nlen]:
            continue
        if component != "*" and (res.component or "A").upper() != component.upper():
            continue
        results.append(res)
    return results


def available_mean_orientations(path: str, pattern: str = "*") -> List[int]:
    """Diagnostic (PAS dans le Fortran) : orientations (`par3_mean`, meme
    convention que `self.orientation` - 1/2/3) sous lesquelles il existe au
    moins une moyenne "mean:" correspondant a `pattern`, IGNORE le filtre
    d'orientation - sert a expliquer un "0 resultat" en mode m/s de
    `selres` quand l'orientation courante ne correspond a aucune moyenne
    enregistree (`_load_mean_results` exige une correspondance EXACTE)."""
    site = (pattern or "*").strip()
    full_pattern = "*" if site in ("", "*") else f"MEAN: {site.upper()}"
    nlen = len(full_pattern) if full_pattern != "*" else 0

    found = set()
    for res in _iter_result_lines(path):
        if res.id[:5] != "mean:":
            continue
        if full_pattern != "*" and res.id[:nlen].upper() != full_pattern[:nlen]:
            continue
        if res.par3_mean == int(res.par3_mean):
            found.add(int(res.par3_mean))
    return sorted(found)


def _load_site_results(
    path: str, pattern: str, iorient: int, component: str = "*"
) -> List[FitResult]:
    """Equivalent de la branche `carselect=="s"` de `selres` (label 444) :
    les moyennes matchees par `_load_mean_results`, PLUS pour chacune les
    resultats individuels references dans son champ `liste`
    ("codes:c1:c2:...", les `c` des resultats combines) - equivalent de
    `decodelisteres` (recherche par `c`, ajout sans autre filtre).
    `component` filtre le CHOIX de la/les moyenne(s) (voir
    _load_mean_results) - les resultats individuels associes suivent
    ensuite sans filtre supplementaire, quelle que soit leur propre
    etiquette component.

    Ordre du resultat : POUR CHAQUE moyenne, ses lignes/plans individuels
    D'ABORD, puis la moyenne elle-meme, repete pour chaque moyenne
    selectionnee - demande explicite utilisateur ("after the selection
    site+best fits, is it possible to list the lines involved in the
    mean, then the mean, for all selected means"). PAS l'ordre du Fortran
    d'origine (qui listait toutes les moyennes d'abord, puis tous les
    specimens ensuite en un seul bloc separe) - purement l'ordre
    d'affichage/de self.results, n'affecte aucun calcul (fisher_from_results
    etc. ne dependent pas de l'ordre)."""
    means = _load_mean_results(path, pattern, iorient, component=component)
    if not means:
        return means

    results: List[FitResult] = []
    for mean_res in means:
        wanted = set()
        for token in mean_res.liste.replace("codes:", "").split(":"):
            token = token.strip()
            if token:
                wanted.add(token)
        if wanted:
            for res in _iter_result_lines(path):
                if res.id[:5] != "mean:" and res.c in wanted:
                    results.append(res)
        results.append(mean_res)
    return results


def load_results(
    path: str,
    pattern: str = "*",
    carselect: str = "d",
    cat1: str = "*",
    numcomp: Optional[int] = None,
    iorient: int = 1,
    component: str = "*",
) -> List[FitResult]:
    """Equivalent de `selres` (dataselect.f), les 3 modes `carselect` :

    - 'd' (Data, par defaut) : resultats normaux (L/P/f/s), filtres par
      `pattern` (id), `cat1`, `numcomp`. Les moyennes "mean:" sont exclues.
    - 'm' (Mean) : uniquement les moyennes de site, `pattern` = le nom du
      site SANS le prefixe "mean: " (ajoute automatiquement), filtrees par
      orientation (`iorient`, doit correspondre a `res.par3`). `cat1`/
      `numcomp` ignores (pas prompted par le Fortran en mode Mean).
    - 's' (Site = Data+Mean) : une moyenne matchee (meme filtre que 'm')
      PLUS les resultats individuels qui la composent (son champ `liste`).

    `component` : filtre PAS dans le Fortran, sur FitResult.component
    (A/B/C..., voir sa docstring) - demande explicite utilisateur pour
    distinguer plusieurs moyennes d'un meme site portant des composantes
    de magnetisation differentes ("when there is different components of
    magnetizations within the same site, we may have two or three means
    by site... We need to add a column component"). "*" (defaut) = pas de
    filtre, applicable aux 3 modes. En mode 'm'/'s', compare l'etiquette
    de LA MOYENNE elle-meme, jamais celle de ses resultats individuels
    constitutifs.
    """
    if not os.path.exists(path):
        return []
    carselect = (carselect or "d").strip().lower()
    if carselect == "m":
        return _load_mean_results(path, pattern, iorient, component=component)
    if carselect == "s":
        return _load_site_results(path, pattern, iorient, component=component)
    return _load_data_results(path, pattern, cat1, numcomp, component=component)


# ---------------------------------------------------------------------------
# mdf : step de demi-desaimantation (MDF si AF, MDT si thermique), par deux
# methodes (norme simple, distance restante le long du chemin) - calcul.f:519-678
# ---------------------------------------------------------------------------

@dataclass
class MdfResult:
    id: str
    demagcod: str   # 'M.D.T.' / 'M.D.F.'
    unite: str
    xmdf_path: float   # xmdf2 : methode "distance restante sur le chemin"
    xmdf_norm: float   # xmdf1 : methode "norme simple du vecteur"
    rapport: float


_MDF_CODES = {
    "D": ("M.D.T.", "deg C"), "d": ("M.D.T.", "deg C"),
    "S": ("M.D.T.", "deg C"), "s": ("M.D.T.", "deg C"),
    "F": ("M.D.F.", "mT"), "f": ("M.D.F.", "mT"),
    "N": ("A.R.N.", "......."), "n": ("A.R.N.", "......."),
    "I": ("A.R.I.", "......."), "i": ("A.R.i.", "......."),
}


def compute_mdf(ech: SelectedSample) -> Optional[MdfResult]:
    """Equivalent de `mdf` : interpolation lineaire du step ou la remanence
    tombe a 50% de sa valeur initiale, par deux methodes. None pour les
    echantillons ARM/AR (dernier cod1 in N/I/n/i), comme le Fortran (qui
    n'imprime alors rien pour cet echantillon)."""
    mesures = ech.mesures
    n = len(mesures)
    if n < 2:
        return None

    steps = [m.etape for m in mesures]
    xs = [m.x for m in mesures]
    ys = [m.y for m in mesures]
    zs = [m.z for m in mesures]

    # steps[j]==steps[j-1] est possible avec des donnees paleointensite
    # (etapes R/V a la meme temperature) - non prevu par le Fortran (division
    # par zero), on saute simplement le croisement dans ce cas (mdf n'a pas
    # de sens pour ce type de protocole de toute facon).
    rint1 = [math.sqrt(xs[i] ** 2 + ys[i] ** 2 + zs[i] ** 2) for i in range(n)]
    rmax1 = rint1[0]
    xmdf1 = 0.0
    if rmax1 != 0:
        for j in range(1, n):
            if steps[j] == steps[j - 1]:
                continue
            if (rint1[j] / rmax1) < 0.5 and (rint1[j - 1] / rmax1) >= 0.5:
                amdf1 = (rint1[j] - rint1[j - 1]) / ((steps[j] - steps[j - 1]) * rmax1)
                bmdf1 = rint1[j] / rmax1 - steps[j] * amdf1
                if amdf1 != 0:
                    xmdf1 = (0.5 - bmdf1) / amdf1

    rint2 = [0.0] * n
    rint2[n - 1] = rint1[n - 1]
    for j in range(n - 2, -1, -1):
        seg = math.sqrt((xs[j] - xs[j + 1]) ** 2 + (ys[j] - ys[j + 1]) ** 2 + (zs[j] - zs[j + 1]) ** 2)
        rint2[j] = rint2[j + 1] + seg
    rmax2 = rint2[0]
    xmdf2 = 0.0
    if rmax2 != 0:
        for j in range(1, n):
            if steps[j] == steps[j - 1]:
                continue
            if (rint2[j] / rmax2) < 0.5 and (rint2[j - 1] / rmax2) >= 0.5:
                amdf2 = (rint2[j] - rint2[j - 1]) / ((steps[j] - steps[j - 1]) * rmax2)
                bmdf2 = rint2[j] / rmax2 - steps[j] * amdf2
                if amdf2 != 0:
                    xmdf2 = (0.5 - bmdf2) / amdf2

    demagcod, unite = _MDF_CODES.get(mesures[-1].cod1, ("A.R.N.", "......."))
    if demagcod in ("A.R.N.", "A.R.I.", "A.R.i."):
        return None

    # `steps` (etape) est directement la valeur physique reelle (mT pour
    # un pas AF, degC pour thermique - plus d'echelle Oersted-equivalente,
    # voir testlect.Measurement) : xmdf1/xmdf2, interpolation LINEAIRE de
    # ces steps, sortent donc DEJA dans la bonne unite - AUCUNE mise a
    # l'echelle supplementaire ici (l'ancienne division par 10, necessaire
    # quand `etape` restait Oersted-scale, ferait desormais doublement
    # petit le resultat).

    rapport = xmdf2 / xmdf1 if xmdf1 != 0 else 0.0
    return MdfResult(id=ech.id, demagcod=demagcod, unite=unite,
                      xmdf_path=xmdf2, xmdf_norm=xmdf1, rapport=rapport)


def list_mdf(selected: List[SelectedSample]) -> str:
    lines = []
    for ech in selected:
        r = compute_mdf(ech)
        if r is None:
            continue
        lines.append(
            f"{r.id} {r.demagcod} norm {r.xmdf_path:7.2f} {r.unite}  "
            f"A.R.N. {r.xmdf_norm:7.2f} {r.unite}  Ratio: {r.rapport:7.2f}"
        )
    return "\n".join(lines) if lines else "(no result - ARM/AR samples excluded, or no 50% crossing)"


# ---------------------------------------------------------------------------
# mdi / mds : moyenne arithmetique/geometrique de l'intensite et de la
# susceptibilite sur toute la selection - calcul.f:683-865
# ---------------------------------------------------------------------------

@dataclass
class MeanIntensityResult:
    id6: str
    n: int
    unit: str
    arith_mean: float
    arith_sd: float
    arith_se: float
    geom_mean: float
    geom_inf: float
    geom_sup: float


@dataclass
class MeanSuscResult:
    n: int
    arith_mean: float
    arith_sd: float
    geom_mean: float
    geom_inf: float
    geom_sup: float


def _arith_geom_stats(values: List[float]) -> tuple:
    n = len(values)
    amoy = sum(values) / n
    somcar = sum((v - amoy) ** 2 for v in values)
    ecart = math.sqrt(somcar / (n - 1)) if n > 1 else 0.0
    erreur = ecart / math.sqrt(n) if n > 1 else 0.0

    logs = [math.log10(v) for v in values if v > 0]
    if len(logs) == n and n > 1:
        amoyg = sum(logs) / n
        somcarg = sum((v - amoyg) ** 2 for v in logs)
        ecartg = math.sqrt(somcarg / (n - 1))
        rgeom, rinf, rsup = 10 ** amoyg, 10 ** (amoyg - ecartg), 10 ** (amoyg + ecartg)
    else:
        rgeom = rinf = rsup = 0.0
    return amoy, ecart, erreur, rgeom, rinf, rsup


def compute_mean_intensity(selected: List[SelectedSample]) -> Optional[MeanIntensityResult]:
    """Equivalent de `mdi` : moyennes arithmetique/geometrique de l'intensite
    (norme du moment) sur TOUTES les mesures de la selection. None si la
    selection melange des normalisations masse/volume (comme le Fortran :
    "on ne peut pas melanger les choux et les carottes"), ou si <=1 mesure.

    None aussi si la selection melange des specimens AVEC et SANS volume/
    masse renseigne (meme categorie d'incompatibilite que masse/volume
    melanges - une moyenne ne peut pas melanger des Am2 bruts avec des A/m
    normalises). Si AUCUN specimen n'a de volume/masse, calcule la moyenne
    du moment total BRUT (Am2) plutot que d'appliquer quand meme le
    facteur 1e3/1e6 a un volume absent (bug Fortran confirme - le Fortran
    d'origine divisait toujours par ech.vol sans le verifier) - demande
    explicite utilisateur ("when the volume or the mass of the sample is
    not given, best to show the data in total moment as done in the
    list")."""
    if not selected:
        return None
    norme0 = selected[0].norme
    if any(ech.norme != norme0 for ech in selected[1:]):
        return None
    has_vol0 = bool(selected[0].vol)
    if any(bool(ech.vol) != has_vol0 for ech in selected[1:]):
        return None

    if not has_vol0:
        unit = "Am2"
        rint = [
            math.sqrt(m.x ** 2 + m.y ** 2 + m.z ** 2)
            for ech in selected for m in ech.mesures
        ]
    else:
        factor = 1.0e3 if norme0 == "m" else 1.0e6
        unit = "Am2/kg" if norme0 == "m" else "A/m"
        rint = [
            math.sqrt(m.x ** 2 + m.y ** 2 + m.z ** 2) * factor / ech.vol
            for ech in selected for m in ech.mesures
        ]
    if len(rint) <= 1:
        return None

    amoy, ecart, erreur, rgeom, rinf, rsup = _arith_geom_stats(rint)
    return MeanIntensityResult(
        id6=(getattr(selected[0], "magic_site", "") or selected[0].id[:6]), n=len(rint), unit=unit,
        arith_mean=amoy, arith_sd=ecart, arith_se=erreur,
        geom_mean=rgeom, geom_inf=rinf, geom_sup=rsup,
    )


def compute_mean_susceptibility(selected: List[SelectedSample]) -> Optional[MeanSuscResult]:
    """Equivalent de `mds`, appele par `mdi`. None si <=2 mesures avec
    susceptibilite non nulle (comme le Fortran, `if(itot==1) return` /
    `if(itot==2) return`).

    Specimens SANS volume/masse (0/None) EXCLUS de cette moyenne - PAS de
    "total brut" honnete pour une susceptibilite (contrairement au moment,
    qui a un equivalent Am2 non normalise significatif - voir
    compute_mean_intensity) : mieux vaut ne pas inclure une valeur
    faussee que d'appliquer quand meme le facteur 1e-4 a un volume absent
    comme le faisait le Fortran d'origine (demande explicite utilisateur,
    "the previous Fortran was systematically dividing by mass and volume
    and was not expecting no vol and no mass")."""
    suscint = [
        m.s * 1e-4 / ech.vol
        for ech in selected for m in ech.mesures if m.s != 0.0 and ech.vol
    ]
    if len(suscint) <= 2:
        return None
    amoy, ecart, _erreur, rgeom, rinf, rsup = _arith_geom_stats(suscint)
    return MeanSuscResult(n=len(suscint), arith_mean=amoy, arith_sd=ecart,
                           geom_mean=rgeom, geom_inf=rinf, geom_sup=rsup)


def format_mean_intensity(
    mi: MeanIntensityResult, ms: Optional[MeanSuscResult] = None,
    lat: float = 0.0, rlong: float = 0.0,
) -> str:
    lines = [
        f"{mi.id6}  Nb: {mi.n:4d}  Mean Arith. Int.: {mi.arith_mean:.3e}"
        f"        sd: {mi.arith_sd:.3e}   error: {mi.arith_se:.3e}",
        f"{mi.id6}  Nb: {mi.n:4d}  Mean Geom. Int.: {mi.geom_mean:.3e}"
        f"  ± sd.inf: {mi.geom_inf:.3e} ± sd.sup: {mi.geom_sup:.3e}",
    ]
    if ms is not None:
        lines.append(
            f"{mi.id6}   {lat:10.5f}  {rlong:10.5f}   "
            f"{mi.geom_mean:.3e}   {mi.geom_mean - mi.geom_inf:.3e}   {mi.geom_sup - mi.geom_mean:.3e}   "
            f"{ms.geom_mean:.3e}   {ms.geom_mean - ms.geom_inf:.3e}   {ms.geom_sup - ms.geom_mean:.3e}"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# meaninc : inclinaison moyenne de McFadden & Reid (1982), estimateur du
# maximum de vraisemblance par iteration de Newton-Raphson - calcul.f:1011-1129
# ---------------------------------------------------------------------------

@dataclass
class MeanInclinationResult:
    n: int
    arith_mean: float
    biased_inc: float     # inclinaison biaisee (Arason)
    unbiased_inc: float   # inclinaison non biaisee (McFadden)
    precision: float
    alpha95: float
    alpha_pos: float
    alpha_neg: float


def compute_mean_inclination(
    selected: List[SelectedSample], orientation: int = 1
) -> Optional[MeanInclinationResult]:
    """Equivalent de `meaninc`. None si <2 mesures ou si l'iteration ne
    converge pas (equivalent implicite d'une boucle Fortran sans limite -
    ici bornee a 1000 iterations par securite)."""
    dirs = []
    for ech in selected:
        for m in ech.mesures:
            xx, yy, zz = apply_orientation(m.x, m.y, m.z, ech, orientation)
            _, _dec, inc = polere(xx, yy, zz)
            dirs.append(inc)

    k = len(dirs)
    if k < 2:
        return None

    eps = 1e-4
    conv = math.pi / 180.0
    som = sum(dirs)
    som1 = som2 = som3 = 0.0
    for d in dirs:
        dd = (90.0 - d) * conv
        som1 += math.cos(dd)
        som2 += math.sin(dd)
        som3 += dd

    amar = som / k
    par = som3 / k

    converged = False
    for _ in range(1000):
        a1 = math.cos(par) ** 2
        a2 = math.sin(par) ** 2
        ff = k * math.cos(par) + (a2 - a1) * som1 - 2 * math.sin(par) * math.cos(par) * som2
        ffp = -k * math.sin(par) + 2 * som1 * math.sin(2 * par) - 2 * som2 * math.cos(2 * par)
        if ffp == 0:
            return None
        xx = par - ff / ffp
        converged = abs(xx - par) < eps
        par = xx
        if converged:
            break
    if not converged:
        return None

    c = math.cos(par) * som1 + math.sin(par) * som2
    s = math.sin(par) * som1 - math.cos(par) * som2
    biais = 180.0 * s / (math.pi * c) if c != 0 else 0.0
    prec = (k - 1.0) / (2.0 * (k - c)) if (k - c) != 0 else 0.0
    ainc = 90.0 - math.degrees(par) + biais
    aai = 90.0 - math.degrees(par)
    resu = ((prec - 1.0) * k + 1.0) / prec if prec != 0 else 0.0
    toto = (0.05 ** (-1.0 / (k - 1.0)) - 1.0) * (k - resu) / resu if resu != 0 else 0.0
    toto = max(0.0, min(2.0, toto))
    alp95 = math.degrees(math.acos(1.0 - toto))

    return MeanInclinationResult(
        n=k, arith_mean=amar, biased_inc=aai, unbiased_inc=ainc,
        precision=prec, alpha95=alp95, alpha_pos=alp95 + biais, alpha_neg=alp95 - biais,
    )


def format_mean_inclination(r: MeanInclinationResult) -> str:
    return (
        f"arithmetic mean: inclination = {r.arith_mean:.2f}\n"
        f"biased inclination (Arason)           : {r.biased_inc:.2f}\n"
        f"unbiased inclination (McFadden)        : {r.unbiased_inc:.2f}\n"
        f"precision parameter                    : {r.precision:.2f}\n"
        f"alpha95 estimate                       : {r.alpha95:.2f}\n\n"
        f"mean inclination: {r.unbiased_inc:.2f} +/- {r.alpha95:.2f}\n"
        f"                    | {r.biased_inc:.2f} +{r.alpha_pos:.2f} -{r.alpha_neg:.2f}"
    )


# ---------------------------------------------------------------------------
# Koenigs : rapport de Koenigsberger Qn = aimantation remanente / aimantation
# induite dans un champ de reference - calcul.f:3869-3926
# ---------------------------------------------------------------------------

@dataclass
class KoenigsbergerRow:
    id: str
    etape: int
    cod1: str
    cod2: str
    mag: float
    ind_mag: float
    ratio: float


def compute_koenigsberger(selected: List[SelectedSample], valk: float) -> List[KoenigsbergerRow]:
    """Equivalent de `Koenigs` : `valk` = champ de reference en microteslas.
    `ratio` (Qn, le resultat analytique principal) est INVARIANT a un
    volume/masse absent (le meme facteur divise rxx et rind, s'annule
    dans le quotient) ; `mag`/`ind_mag` (valeurs intermediaires affichees)
    ne recoivent pas le facteur 1e3/1e6/1e-7 sans volume/masse (0/None) -
    demande explicite utilisateur ("the previous Fortran was
    systematically dividing by mass and volume and was not expecting no
    vol and no mass"), meme principe que compute_mean_intensity."""
    valc = valk / (4 * math.pi / 10.0)
    rows = []
    for ech in selected:
        vol = ech.vol or 1.0
        for m in ech.mesures:
            if m.s == 0.0:
                continue
            rxx = math.sqrt(m.x ** 2 + m.y ** 2 + m.z ** 2)
            if not ech.vol:
                rind = valc * m.s
            elif ech.norme == "m":
                rxx = rxx * 1.0e3 / vol
                rind = valc * m.s * 1e-7 / vol
            else:
                rxx = rxx * 1.0e6 / vol
                rind = valc * m.s * 10.0 * 1e-5 / vol
            ratio = rxx / rind if rind != 0 else 0.0
            rows.append(KoenigsbergerRow(
                id=ech.id, etape=m.etape, cod1=m.cod1, cod2=m.cod2,
                mag=rxx, ind_mag=rind, ratio=ratio,
            ))
    return rows


def format_koenigsberger(rows: List[KoenigsbergerRow]) -> str:
    lines = [f"{_HEADER_MARK}Sample        step       Mag           K            Koenigsberger ratio{_HEADER_MARK}"]
    for r in rows:
        fake_m = Measurement(etape=r.etape, cod1=r.cod1, cod2=r.cod2, x=0, y=0, z=0, q=0, ins="", s=0)
        lines.append(
            f"{r.id:<12s}  {_fmt_step(fake_m)}  {r.mag:.3e}  {r.ind_mag:.3e}   {r.ratio:7.3f}"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# diffmes : difference vectorielle entre mesures consecutives (etape j moins
# etape j+1) - dataselect.f:1156-1236
# ---------------------------------------------------------------------------

def list_diff_measurements(selected: List[SelectedSample], orientation: int = 1) -> str:
    """Equivalent de `diffmes`."""
    lines = []
    ij = 1
    for ech in selected:
        for j in range(len(ech.mesures) - 1):
            a, b = ech.mesures[j], ech.mesures[j + 1]
            dx, dy, dz = a.x - b.x, a.y - b.y, a.z - b.z
            xx, yy, zz = apply_orientation(dx, dy, dz, ech, orientation)
            mag, dec, inc = polere(xx, yy, zz)
            # sans volume/masse (0/None) : moment total brut (Am2), pas un
            # faux nombre issu d'une division par un volume absent (bug
            # Fortran confirme) - demande explicite utilisateur.
            mag, _unit = normalized_intensity(ech, mag)
            lines.append(
                f"{ij:4d}: {ech.id:<12s}  {_fmt_step(a)}  "
                f"{mag:10.3e} {dec:6.1f} {inc:6.1f}  q={a.q:<4d} {a.ins:<2s} s={a.s:.1f}"
            )
            ij += 1
    return "\n".join(lines) if lines else "(aucune mesure)"


# ---------------------------------------------------------------------------
# viscos : indice de viscosite sur des paires N+/N- - viscos.f
# MUTE selected en place (comme le Fortran).
# ---------------------------------------------------------------------------

def apply_viscosity_test(selected: List[SelectedSample]) -> List[str]:
    """Equivalent de `viscos` : pour chaque echantillon dont les deux
    premieres mesures sont codees N+/N-, remplace la premiere par le
    vecteur MOYEN et la seconde par le vecteur DIFFERENCE (cod2 mis a
    '0'/'v'), .q de la seconde recevant l'indice de viscosite (%). Retourne
    les avertissements (echantillons dont les 2 premieres mesures ne sont
    pas N+/N-, comme les 2 messages `write` du Fortran)."""
    warnings = []
    for ech in selected:
        if len(ech.mesures) < 2:
            continue
        m0, m1 = ech.mesures[0], ech.mesures[1]
        if not (m0.cod1 == "N" and m1.cod1 == "N"):
            warnings.append(f"{ech.id} : attention code different de N et N")
            continue
        if not (m0.cod2 == "+" and m1.cod2 == "-"):
            warnings.append(f"{ech.id} : attention code different de + et -")
            continue

        xm, ym, zm = (m0.x + m1.x) / 2.0, (m0.y + m1.y) / 2.0, (m0.z + m1.z) / 2.0
        xd, yd, zd = (m0.x - m1.x) / 2.0, (m0.y - m1.y) / 2.0, (m0.z - m1.z) / 2.0
        iq1 = (m0.q + m1.q) // 2
        num = math.sqrt(xd ** 2 + yd ** 2 + zd ** 2)
        den = math.sqrt(xm ** 2 + ym ** 2 + zm ** 2)
        visco = 100.0 * num / den if den != 0 else 0.0

        m0.x, m0.y, m0.z = xm, ym, zm
        m0.q = iq1
        m0.cod2 = "0"
        m1.x, m1.y, m1.z = xd, yd, zd
        m1.q = int(visco)
        m1.cod2 = "v"
    return warnings


# ---------------------------------------------------------------------------
# soustra : soustrait le vecteur brut d'une mesure de toutes les autres de
# l'echantillon, puis supprime cette ligne - dataselect.f:1239-1283
# MUTE ech.mesures en place.
# ---------------------------------------------------------------------------

def apply_subtraction(ech: SelectedSample, row_index: int) -> Optional[str]:
    """Equivalent de `soustra`. `row_index` : 1-based (comme les autres
    prompts "numero de ligne"). Retourne un message d'erreur si invalide,
    sinon None."""
    idx = row_index - 1
    if not (0 <= idx < len(ech.mesures)):
        return "invalid line number"
    xs, ys, zs = ech.mesures[idx].x, ech.mesures[idx].y, ech.mesures[idx].z
    for m in ech.mesures:
        m.x -= xs
        m.y -= ys
        m.z -= zs
    del ech.mesures[idx]
    return None


# ---------------------------------------------------------------------------
# holderarm : enregistre les 6 mesures ARM d'un porte-echantillon vide
# (x+,x-,y+,y-,z+,z-), pour une future correction dans `anisot` (non porte) -
# calcul.f:2337-2414
# ---------------------------------------------------------------------------

@dataclass
class ArmHolderBackground:
    x: Tuple[float, float, float, float, float, float]
    y: Tuple[float, float, float, float, float, float]
    z: Tuple[float, float, float, float, float, float]


def record_arm_holder(
    ech: SelectedSample, ixp: int, ixm: int, iyp: int, iym: int, izp: int, izm: int
) -> Optional[ArmHolderBackground]:
    """Equivalent de `holderarm`. Indices 1-based dans ech.mesures. None si
    un indice est invalide."""
    idx = [ixp, ixm, iyp, iym, izp, izm]
    if any(not (1 <= i <= len(ech.mesures)) for i in idx):
        return None
    xs = tuple(ech.mesures[i - 1].x for i in idx)
    ys = tuple(ech.mesures[i - 1].y for i in idx)
    zs = tuple(ech.mesures[i - 1].z for i in idx)
    return ArmHolderBackground(x=xs, y=ys, z=zs)


# ---------------------------------------------------------------------------
# vitref : correction de vitesse de refroidissement (TRM rapide labo vs TRM
# lente en four), a partir d'un motif de mesures L/Q (etape lente) precede
# de R/V (etapes rapides, meme cod2) - calcul.f:3461-3625 (branche "live",
# le reste apres `go to 555` est mort/inutilise cote Fortran).
# ---------------------------------------------------------------------------

@dataclass
class CoolingRateResult:
    id: str
    correction_avant: float
    correction_apres: float
    correction_mean: float
    evolution_pct: float
    evolution_lente_avant_pct: float
    evolution_lente_apres_pct: float


def detect_cooling_rate_rows(ech: SelectedSample) -> Optional[Tuple[int, int, int, int, int]]:
    """Equivalent du mode automatise de `vitref` : cherche un pas 'L' suivi
    d'un pas 'Q', puis dans les 8 mesures precedentes les pas 'R' et 'V' de
    MEME cod2. Retourne (irr,irv,ilr,iravant,irapres) en indices 1-based
    (memes conventions que les autres prompts "numero de ligne"), ou None."""
    mesures = ech.mesures
    for j in range(len(mesures) - 1):
        if mesures[j].cod1 == "L" and mesures[j + 1].cod1 == "Q":
            ilr, irapres = j, j + 1
            irr = irv = iravant = None
            for ki in range(1, 9):
                idx = j - ki
                if idx < 0:
                    break
                if mesures[idx].cod1 == "R" and mesures[idx].cod2 == mesures[j].cod2:
                    irr = idx
                if mesures[idx].cod1 == "V" and mesures[idx].cod2 == mesures[j].cod2:
                    irv = idx
                    iravant = idx
            if irr is not None and irv is not None:
                return irr + 1, irv + 1, ilr + 1, iravant + 1, irapres + 1
    return None


def compute_cooling_rate(
    ech: SelectedSample, irr: int, irv: int, ilr: int, iravant: int, irapres: int
) -> CoolingRateResult:
    """Equivalent du calcul de `vitref` (partie "live" avant le `go to 555`
    mort). Indices 1-based."""
    m = ech.mesures
    vol = ech.vol or 1.0

    def mag(x, y, z):
        mg, _dec, _inc = polere(x, y, z)
        return mg * 1.0e6 / vol

    r, v, l, av, ap = m[irr - 1], m[irv - 1], m[ilr - 1], m[iravant - 1], m[irapres - 1]

    arnx, arny, arnz = (r.x + v.x) / 2.0, (r.y + v.y) / 2.0, (r.z + v.z) / 2.0
    atr_rapide_mag = mag((r.x - v.x) / 2.0, (r.y - v.y) / 2.0, (r.z - v.z) / 2.0)  # noqa: F841 (echo Fortran)
    atrlente = mag(l.x - arnx, l.y - arny, l.z - arnz)
    atrrapavant = mag(av.x - arnx, av.y - arny, av.z - arnz)
    atrrapapres = mag(ap.x - arnx, ap.y - arny, ap.z - arnz)

    corr_avant = atrrapavant / atrlente if atrlente else 0.0
    corr_apres = atrrapapres / atrlente if atrlente else 0.0
    evol = ((atrrapapres / atrrapavant) - 1.0) * 100.0 if atrrapavant else 0.0
    evol_la = ((atrlente / atrrapavant) - 1.0) * 100.0 if atrrapavant else 0.0
    evol_lp = ((atrlente / atrrapapres) - 1.0) * 100.0 if atrrapapres else 0.0

    return CoolingRateResult(
        id=ech.id, correction_avant=corr_avant, correction_apres=corr_apres,
        correction_mean=(corr_avant + corr_apres) / 2.0,
        evolution_pct=evol, evolution_lente_avant_pct=evol_la, evolution_lente_apres_pct=evol_lp,
    )


def format_cooling_rate(r: CoolingRateResult) -> str:
    return (
        f"Sample: {r.id}  % Evolution: {r.evolution_pct:.1f}%"
        f"  slow/Fast before: {r.evolution_lente_avant_pct:.1f}%"
        f"  slow/Fast after: {r.evolution_lente_apres_pct:.1f}%\n"
        f"Sample: {r.id}  cooling correction (1st fast): {r.correction_avant:.3f}\n"
        f"Sample: {r.id}  cooling correction (last fast): {r.correction_apres:.3f}\n"
        f"Sample: {r.id}  mean cooling correction: {r.correction_mean:.3f}"
    )


# ---------------------------------------------------------------------------
# anisot / anisoauto : calcul du tenseur d'anisotropie (ARM/TRM 6 positions)
# - calcul.f:2418-3310 (anisot, saisie manuelle) / 3951-4900+ (anisoauto,
# detection automatique). Porte desormais LES 15 VARIANTES ecrites par le
# Fortran (A0/A+/A-/A1/B1/A2/B2/A3/B3/A4/B4/A5/B5/A6/B6 - demande explicite
# utilisateur "implement the 15 outputs for the anisotropy tensor with the
# A0 code being the major one" - "Etape 2" de la demande initiale, qui
# avait deliberement differe ces 14 variantes "jackknife"). Le bloc de
# correction d'alteration/evolution (testevol) reste hors perimetre.
#
# Principe des 14 variantes au-dela de A0 (calcul.f:2691-3256) : A0 =
# demi-difference +/- par paire d'axes (s'annule pour tout offset commun
# aux deux mesures, ex. NRM residuelle - LE tenseur "officiel"). Les
# variantes A+/A- reprennent les MEMES composantes que A0 mais SANS les
# symetriser (A+ = triangle "ligne i", A- = triangle "colonne i", i.e.
# "transpose"). Les paires (A1,B1)/(A2,B2) soustraient la MOYENNE (pas la
# demi-difference) de la paire X (`anrmx*`, une estimation du "bruit"/
# remanence parasite commune a X+ et X-) aux 3 mesures "+" (A1) ou "-"
# (A2, signe inverse) de CHAQUE paire d'axes - (A3,B3)/(A4,B4) et
# (A5,B5)/(A6,B6) refont le meme calcul en soustrayant la moyenne de la
# paire Y, puis de la paire Z, respectivement. Comparer ces 14 variantes
# a A0 sert de diagnostic de robustesse (une forte divergence signale que
# la remanence parasite/NRM residuelle contamine significativement le
# tenseur mesure) - aucune n'est "la" mesure a utiliser en aval (voir
# compute_anicor_factor, qui continue a n'utiliser QUE A0).
#
# VERIFIE OCTET-PRES contre un vrai fichier .ANI fourni par l'utilisateur
# (Miriam_2025b/SanJuan_Pmag.ANI, echantillon "02A", etape 460) : calcul a
# la main des 6 valeurs k11/k22/k33/k12/k23/k13 a partir des mesures
# reelles du .ren (460X+/460X-/460Y+/460Y-/460RH/460VH) - identiques a la
# ligne 'A0' reelle du .ANI (aucune ligne de base porte-echantillon dans
# ce cas, holder=None).
# ---------------------------------------------------------------------------

@dataclass
class AniTensor:
    id: str
    code2: str  # 'A0'=TRM, 'F0'=ARM (calcul natif/AMS_Py), 'N0'=susceptibilite,
    # 'AA'=ARM DEJA calcule par une contribution MagIC importee (voir
    # extract_magic.magic_anisotropy_rows) - deliberement distinct de 'F0'
    k11: float
    k22: float
    k33: float
    k12: float
    k23: float
    k13: float
    # Champs ajoutes pour le format .pmagani (voir write_ani_tensors/
    # read_ani_tensor) - demande explicite utilisateur ("What will be the
    # most complete .pmagani file") : nombre de positions ayant servi au
    # calcul et statistiques de Hext (sigma, test F global/F12/F23), issues
    # UNIQUEMENT du chemin PmagPy (voir anisotropy_magic.compute_aarm_
    # pmagpy) - le calcul natif STARpaleomag_Py (jackknife geometrique, calcul.f)
    # ne les calcule pas ; None/absent si non disponible.
    n_positions: Optional[int] = None
    sigma: Optional[float] = None
    ftest: Optional[float] = None
    ftest12: Optional[float] = None
    ftest23: Optional[float] = None
    # f_crit/quality : valeur F critique (95%, Hext 1963) et verdict 'g'
    # (F>F_crit, anisotropie significative)/'b' (non significative), issus
    # de pmagpy (params["aniso_ftest_quality"]) - PAS des colonnes .pmagani
    # dediees (le fichier reste au meme schema), seulement utilises pour
    # composer la note "satisfactory/not satisfactory" dans la colonne
    # info (voir _format_pmagani_line) - demande explicite utilisateur
    # ("ajouter le calcul de PmagPy et les erreurs dans pmagani, ainsi que
    # l'estimation dans le comment satisfactory or not satisfactory").
    f_crit: Optional[float] = None
    quality: Optional[str] = None
    # "Y"/"N" - colonne dediee ajoutee cote AMS_Py (ams_selection.
    # AMSMeasurement.export/mark_pmagani_export, MEME fichier .pmagani
    # partage) - demande explicite utilisateur ("is it possible to export
    # from AMS_py only the data and mean tensors that we want to export
    # and only these selected data will be taken into account in the
    # main export from Starpaleomag") : "N" signifie que ce tenseur doit
    # etre IGNORE par ouvrir_export_magic_dialog. "Y" par defaut (ancien
    # fichier jamais marque cote AMS_Py = tout exporte, comme avant
    # l'ajout de cette colonne).
    export: str = "Y"


@dataclass
class AniMeanTensor:
    """Une moyenne tensorielle au niveau SITE (PAS specimen) dans le
    fichier .pmagani - demande explicite utilisateur ("can we have the
    pmagani as pmagres with two blocks one at a specimen level, the next
    one at site level with both blocks having their own set of
    parameters"), meme principe que .pmagres (voir _SPECIMEN_HEADER/
    _MEAN_HEADER plus haut dans ce module) : deux sections dans le MEME
    fichier, chacune avec ses PROPRES colonnes plutot que de forcer les
    deux niveaux dans un seul schema commun. Les 3 axes propres (k1>=
    k2>=k3, convention Jelinek/AMS_Py - PAS k11/k22/k33 bruts comme le
    niveau specimen, qui n'a pas de notion d'axes propres) avec leur
    dec/inc et l'ellipse de confiance a 95% (alpha1/alpha2 = demi-grand/
    demi-petit axe, meme convention que AMS_Py.ams_stats.AxisResult -
    voir sa docstring pour l'ecart documente avec le formalisme Hext de
    pmagpy). `id` porte le nom du SITE (pas un specimen) - `n` : nombre
    de specimens effectivement utilises dans la moyenne (apres exclusion
    des tenseurs isotropes, meme convention que TensorialMeanResult.n
    cote AMS_Py)."""
    id: str  # nom du site
    code2: str
    n: int
    k1: float
    k2: float
    k3: float
    dec1: float
    inc1: float
    dec2: float
    inc2: float
    dec3: float
    inc3: float
    alpha1_1: float = 0.0
    alpha2_1: float = 0.0
    alpha1_2: float = 0.0
    alpha2_2: float = 0.0
    alpha1_3: float = 0.0
    alpha2_3: float = 0.0
    # Parametres de forme (Jelinek 1978, meme convention/formules que
    # AMS_Py.ams_stats.shape_params - PAS recalcules ici, fournis par
    # l'appelant) - "n.d" possible si non fournis (ex. tenseur moyen
    # isotrope, ou provenance qui ne les calcule pas).
    P: Optional[float] = None
    T: Optional[float] = None
    L: Optional[float] = None
    F: Optional[float] = None
    Pprim: Optional[float] = None
    # Orientation de la moyenne (Sa/IS/TC) - colonne DEDIEE (demande
    # explicite utilisateur "ajouter une colonne avant info", suite a
    # "in the file pmagani, there is no information about the IS/TC") :
    # meme convention fichier "1"/"0"/"100" que la colonne "IS/TC" de
    # .pmagres (voir _ORIENT_TO_FILE_CODE/_FILE_CODE_TO_ORIENT) - PAS le
    # code MagIC (qui utiliserait "-1" au lieu de "1"). "" = inconnue
    # (fichier ecrit avant l'ajout de cette colonne - voir magic_export.
    # _tilt_correction_from_info pour le repli sur `info`, seule source
    # pour un tel fichier).
    tilt_correction: str = ""
    # Meme colonne "export" (voir AniTensor.export ci-dessus), niveau
    # site mean - MEME cle (site, code2, tilt_correction) que
    # ams_selection.mark_pmagani_export utilise cote AMS_Py pour cibler
    # une ligne mean precise.
    export: str = "Y"
    info: str = ""


_ANI_CODE2 = {1: "A0", 2: "F0", 3: "N0", 4: "AA"}
_ANI_POSITION_KEYS = ("X+", "X-", "Y+", "Y-", "Z+", "Z-")


def detect_six_positions(
    ech: "SelectedSample", force_zplus_label: Optional[str] = None,
) -> Optional[Dict[str, Measurement]]:
    """Equivalent de la detection automatique des 6 positions dans
    anisoauto (calcul.f:4005-4059) : trouve l'etape des mesures 'X'
    (cod1='X', X+/X-), puis a cette MEME etape les mesures codees 'R' et
    'V' (cod1='R'/'V') servent de substituts pour Z+/Z- quand
    l'echantillon n'a pas de vraies mesures Z+/Z- - le "label" de
    substitution est cod1+cod2 de la mesure trouvee (ex. 'RH'/'VH',
    verifie sur donnees reelles). `force_zplus_label` : remplace le label
    Z+ detecte (ex. 'ZB' - voir _check_trm_evolution) plutot que le
    substitut R normal, pour la reclassification apres detection d'une
    evolution/alteration significative (calcul.f:4107 `Zplus="ZB"`).
    Retourne None si l'etape X est absente ou si les 6 positions ne sont
    pas toutes identifiees (equivalent "decoding incomplete" du Fortran).

    BUG Fortran CONFIRME et CORRIGE ici (pas silencieusement reproduit -
    demande explicite utilisateur : "the determination of anisotropy does
    not work; I tried on the magic contribution 20595") : le Fortran
    d'origine (calcul.f:4011-4013) cherche un substitut 'R'/'V' a l'etape
    X SANS JAMAIS VERIFIER qu'une vraie mesure Z+/Z- (cod1='Z') existe
    deja a cette meme etape - un 'R'/'V' present pour une AUTRE raison (ex.
    une experience de paleointensite sur le meme specimen dont un pas
    partage accidentellement le meme numero d'etape que l'anisotropie)
    ecrase alors la vraie mesure Z+/Z- deja trouvee (dernier match gagne,
    dans le Fortran comme dans le port). Confirme sur une contribution
    MagIC reelle (20595, specimen HP01-01) : etape 600 porte a la fois
    Z+/Z- reels ET des mesures 'R0'/'RN' issues d'un protocole Thellier
    separe sur le meme specimen - le tenseur resultant etait corrompu
    (k33 proche de zero, incoherent avec k11/k22). Desormais, le substitut
    R/V n'est cherche QUE pour l'axe (Z+ et/ou Z-) qui n'a PAS deja de
    vraie mesure a cette etape."""
    item_rv = None
    for m in ech.mesures:
        if m.cod1 == "X":
            item_rv = m.etape
    if item_rv is None:
        return None

    has_real_zplus = any(m.etape == item_rv and m.cod1 == "Z" and m.cod2 == "+" for m in ech.mesures)
    has_real_zminus = any(m.etape == item_rv and m.cod1 == "Z" and m.cod2 == "-" for m in ech.mesures)

    zplus_label = zminus_label = None
    for m in ech.mesures:
        if m.etape != item_rv:
            continue
        if m.cod1 == "R" and not has_real_zplus:
            zplus_label = m.cod1 + m.cod2
        elif m.cod1 == "V" and not has_real_zminus:
            zminus_label = m.cod1 + m.cod2
    if force_zplus_label:
        zplus_label = force_zplus_label

    positions: Dict[str, Measurement] = {}
    for m in ech.mesures:
        testcode = m.cod1 + m.cod2
        if zplus_label and testcode == zplus_label:
            testcode = "Z+"
        elif zminus_label and testcode == zminus_label:
            testcode = "Z-"
        if testcode in _ANI_POSITION_KEYS:
            positions[testcode] = m

    if len(positions) != 6:
        return None
    return positions


def _position_vector(
    positions: Dict[str, Measurement], key: str, idx: int,
    holder: Optional[ArmHolderBackground],
) -> Tuple[float, float, float]:
    """Composantes (x,y,z) de la position `key` (idx 0..5 dans l'ordre
    X+,X-,Y+,Y-,Z+,Z- - meme ordre que holder.x/y/z), ligne de base
    porte-echantillon deja soustraite si fournie."""
    m = positions[key]
    x, y, z = m.x, m.y, m.z
    if holder is not None:
        x -= holder.x[idx]
        y -= holder.y[idx]
        z -= holder.z[idx]
    return x, y, z


@dataclass
class AnisotropyPositionDiag:
    """Une ligne de diagnostic 'TRM values' (calcul.f:4340-4477) : mesure
    d'une position apres soustraction de la NRM residuelle moyenne
    (nrm_mean), en unites physiques (A/m ou Am2/kg) + dec/inc."""
    key: str
    measurement: Measurement
    intensity: float
    dec: float
    inc: float


@dataclass
class AnisotropyPairNrm:
    """Une ligne 'ARN X+X-'/'ARN Y+Y-'/'ARN Z+Z-' (calcul.f:4206-4328) :
    estimation de la NRM residuelle a partir d'UNE SEULE paire d'axes
    (moyenne (position+ + position-)/2), AVANT la moyenne des 3 qui donne
    `nrm_mean`/"Mean NRM" - demande explicite utilisateur ("I want also
    the listing of the three NRM determined from each pair... and not
    only the mean") : les 3 valeurs prises separement font apparaitre si
    UNE paire specifique porte une NRM residuelle nettement differente
    des deux autres (probleme localise a cette paire), ce que la seule
    moyenne des 3 ne montre pas."""
    key: str  # "X+X-" / "Y+Y-" / "Z+Z-"
    intensity: float
    dec: float
    inc: float


def _nrm_mean_and_diag(
    positions: Dict[str, Measurement], holder: Optional[ArmHolderBackground],
    norme: str, vol: float,
) -> Tuple[Tuple[float, float, float], List[AnisotropyPositionDiag], float,
           Dict[str, Tuple[float, float, float]], List[AnisotropyPairNrm], AnisotropyPairNrm]:
    """Equivalent de calcul.f:4206-4479 : moyenne des 3 paires (X+/X-,
    Y+/Y-, Z+/Z-) -> `nrm_mean` (la "NRM residuelle" commune aux 6
    positions - une anisotropie basee sur des TRM PARTIELLES laisse
    souvent une part de NRM non remplacee, qui s'annule dans le tenseur
    final par difference mais fausse visuellement chaque position prise
    isolement - demande explicite de l'utilisateur). Retourne (nrm_mean,
    liste de diagnostics par position, deviation_pct (`deviatTRM`),
    composantes brutes NRM-soustraites par position - utilisees par
    _check_position_inversion -, liste des 3 estimations PAR PAIRE de la
    NRM residuelle - voir AnisotropyPairNrm -, et cette MEME nrm_mean
    exprimee en moment/dec/inc plutot qu'en x/y/z brut - demande
    explicite utilisateur ("provide Mean NRM... in moment dec inc; it is
    easier for the user to compare with the three NRM above"), meme
    format que les 3 estimations par paire pour une comparaison directe)."""
    idx_of = {"X+": 0, "X-": 1, "Y+": 2, "Y-": 3, "Z+": 4, "Z-": 5}
    vecs = {k: _position_vector(positions, k, idx_of[k], holder) for k in _ANI_POSITION_KEYS}

    pair_mean_x = tuple((vecs["X+"][c] + vecs["X-"][c]) / 2 for c in range(3))
    pair_mean_y = tuple((vecs["Y+"][c] + vecs["Y-"][c]) / 2 for c in range(3))
    pair_mean_z = tuple((vecs["Z+"][c] + vecs["Z-"][c]) / 2 for c in range(3))
    nrm_mean = tuple((pair_mean_x[c] + pair_mean_y[c] + pair_mean_z[c]) / 3 for c in range(3))

    def scale(intensity: float) -> float:
        return intensity * 1.0e3 / vol if norme == "m" else intensity * 1.0e6 / vol

    def pair_nrm(key: str, pair_mean: Tuple[float, float, float]) -> AnisotropyPairNrm:
        mag, dec, inc = polere(*pair_mean)
        return AnisotropyPairNrm(key=key, intensity=scale(mag), dec=dec, inc=inc)

    pair_nrms = [
        pair_nrm("X+X-", pair_mean_x), pair_nrm("Y+Y-", pair_mean_y), pair_nrm("Z+Z-", pair_mean_z),
    ]
    nrm_mean_diag = pair_nrm("Mean", nrm_mean)

    diags: List[AnisotropyPositionDiag] = []
    raw_components: Dict[str, Tuple[float, float, float]] = {}
    sommetrm = 0.0
    trmmean = 0.0
    prev_raw = None
    for key in _ANI_POSITION_KEYS:
        raw = tuple(vecs[key][c] - nrm_mean[c] for c in range(3))
        raw_components[key] = raw
        mag, dec, inc = polere(*raw)
        trmmean += mag
        diags.append(AnisotropyPositionDiag(
            key=key, measurement=positions[key], intensity=scale(mag), dec=dec, inc=inc,
        ))
        if key in ("X-", "Y-", "Z-"):
            sommetrm += math.sqrt(sum((raw[c] + prev_raw[c]) ** 2 for c in range(3)))
        prev_raw = raw
    sommetrm /= 3.0
    trmmean /= 6.0
    deviation_pct = 100.0 * sommetrm / trmmean if trmmean else 0.0

    return nrm_mean, diags, deviation_pct, raw_components, pair_nrms, nrm_mean_diag


def replace_position_by_symmetry(
    positions: Dict[str, Measurement], holder: Optional[ArmHolderBackground], replace_key: str,
) -> Dict[str, Measurement]:
    """Reconstruit UNE position suspecte a partir de SA PROPRE PAIRE et de
    la moyenne des DEUX AUTRES paires - port du bloc DESACTIVE de
    calcul.f:4481-4552 (le meme "c ... voulez-vous supprimer une des 6
    mesures ? laquelle:" que le prompt lui-meme, jamais active dans la
    version livree) - demande explicite utilisateur ("when we remove one
    line, it is replaced by the antiparallel... also easier to write the
    number of the line from 1 to 6").

    PAS une nouvelle mesure a choisir : `new_vector = -partner_vector +
    2 * avg(other_two_pair_means)` - ou `partner_vector` est l'AUTRE
    membre de la MEME paire (X- pour X+, etc.) et `other_two_pair_means`
    les moyennes (position+ + position-)/2 des deux paires NON concernees
    (dans lesquelles la mesure suspecte n'intervient pas). Le raisonnement :
    en l'absence de tout probleme, chaque paire devrait donner la MEME
    estimation de la NRM residuelle commune aux 6 positions ; si UNE
    position est suspecte, les deux AUTRES paires restent une estimation
    fiable de cette NRM commune, et permettent de reconstruire la position
    manquante en supposant sa propre paire parfaitement antiparallele une
    fois recalee sur ce niveau de NRM partage - PAS un simple `-partner`
    (qui ignorerait completement la NRM residuelle commune).

    Retourne une COPIE de `positions` avec seulement `replace_key` change
    (vers une Measurement SYNTHETIQUE, meme etape/cod1/cod2/q/ins/s que
    l'originale, x/y/z recalcules) - la mesure BRUTE d'origine dans
    ech.mesures n'est JAMAIS modifiee, contrairement au Fortran qui
    mutait `mes(ixp)` directement (aurait corrompu la mesure source pour
    le reste de la session)."""
    idx_of = {"X+": 0, "X-": 1, "Y+": 2, "Y-": 3, "Z+": 4, "Z-": 5}
    vecs = {k: _position_vector(positions, k, idx_of[k], holder) for k in _ANI_POSITION_KEYS}
    pair_mean = {
        axis: tuple((vecs[f"{axis}+"][c] + vecs[f"{axis}-"][c]) / 2 for c in range(3))
        for axis in ("X", "Y", "Z")
    }

    axis, sign = replace_key[0], replace_key[1]
    partner_key = axis + ("-" if sign == "+" else "+")
    other_axes = [a for a in ("X", "Y", "Z") if a != axis]
    anrm_other = tuple(
        (pair_mean[other_axes[0]][c] + pair_mean[other_axes[1]][c]) / 2 for c in range(3))
    partner_vec = vecs[partner_key]
    new_vec = tuple(-partner_vec[c] + 2 * anrm_other[c] for c in range(3))

    # _position_vector soustrait la ligne de base porte-echantillon - la
    # rajouter ici pour que la Measurement synthetique, une fois repassee
    # dans le pipeline normal (qui la soustrait de nouveau), reconstruise
    # exactement `new_vec`.
    if holder is not None:
        idx = idx_of[replace_key]
        new_vec = (new_vec[0] + holder.x[idx], new_vec[1] + holder.y[idx], new_vec[2] + holder.z[idx])

    original = positions[replace_key]
    synthetic = replace(original, x=new_vec[0], y=new_vec[1], z=new_vec[2])
    new_positions = dict(positions)
    new_positions[replace_key] = synthetic
    return new_positions


# Direction TRM ATTENDUE pour chaque position (repere specimen, meme
# repere que Measurement.x/y/z) - PAS dans le Fortran (celui-ci ne teste
# que le SIGNE de X/Y, voir _check_position_inversion ci-dessous) :
# demande explicite utilisateur ("parfois par erreur, les échantillons
# ont été mal orientés dans le four... est-ce possible de mettre un
# warning lorsque le code ne correspond pas à la vraie direction
# enregistrée"). Un champ fort applique selon +X (ou -X, +Y...) impose
# une TRM le long de cet axe ; l'anisotropie devie cette direction du
# champ applique, mais reste dans un cone raisonnable (l'utilisateur
# fixe 45 deg - une anisotropie meme forte, quelques dizaines de %,
# n'ecarte jamais la TRM a plus de la moitie d'un angle droit du champ
# applique) : au-dela, c'est le signe d'une position reellement inversee
# lors de l'acquisition (four/porte-echantillon), pas un effet
# d'anisotropie.
_EXPECTED_POSITION_DIRECTION = {
    "X+": (1.0, 0.0, 0.0), "X-": (-1.0, 0.0, 0.0),
    "Y+": (0.0, 1.0, 0.0), "Y-": (0.0, -1.0, 0.0),
    "Z+": (0.0, 0.0, 1.0), "Z-": (0.0, 0.0, -1.0),
}
_POSITION_DEVIATION_WARNING_DEG = 45.0


def _position_deviation_deg(vec: Tuple[float, float, float], expected: Tuple[float, float, float]) -> float:
    """Angle (deg) entre `vec` (mesure, ligne de base deja soustraite) et
    `expected` (axe theorique de cette position, voir
    _EXPECTED_POSITION_DIRECTION) - 0.0 si `vec` est nul (rien a comparer,
    pas une deviation infinie)."""
    mag = math.sqrt(sum(c * c for c in vec))
    if mag == 0.0:
        return 0.0
    dot = sum((c / mag) * e for c, e in zip(vec, expected))
    dot = max(-1.0, min(1.0, dot))
    return math.degrees(math.acos(dot))


def _find_misoriented_positions(
    positions: Dict[str, Measurement], holder: Optional["ArmHolderBackground"],
) -> List[Tuple[str, float]]:
    """Verifie les 6 positions (pas seulement X/Y comme _check_position_
    inversion) contre leur direction attendue - retourne [(position,
    deviation_deg), ...] pour chaque position au-dela de
    _POSITION_DEVIATION_WARNING_DEG, triees par deviation decroissante -
    [] si aucune. Demande explicite utilisateur ("mettre un warning
    lorsque le code ne correspond pas à la vraie direction enregistrée
    par l'échantillon")."""
    idx_of = {"X+": 0, "X-": 1, "Y+": 2, "Y-": 3, "Z+": 4, "Z-": 5}
    warnings = []
    for key, expected in _EXPECTED_POSITION_DIRECTION.items():
        vec = _position_vector(positions, key, idx_of[key], holder)
        dev = _position_deviation_deg(vec, expected)
        if dev > _POSITION_DEVIATION_WARNING_DEG:
            warnings.append((key, dev))
    warnings.sort(key=lambda item: -item[1])
    return warnings


def _check_position_inversion(raw: Dict[str, Tuple[float, float, float]]) -> Optional[str]:
    """Equivalent de calcul.f:4386-4393 (X) et 4420-4427 (Y) : detecte une
    inversion de la position lors de l'acquisition TRM - le "+"-composant
    (NRM-soustrait) d'une position "+" doit etre positif, celui de la
    position "-" negatif (physiquement, ce sont des champs opposes).
    Fidele au Fortran, y compris son ASYMETRIE : le test X compare le
    composant BRUT de X+ a la MAGNITUDE (toujours positive) de X-,
    reduisant de facto la condition a "composant X+ negatif" ; le test Y
    compare les deux composants BRUTS (signes) comme attendu. Ce n'est
    probablement pas intentionnel cote Fortran, mais reproduit tel quel
    par fidelite - PAS de test Z equivalent dans le Fortran d'origine.
    Retourne 'X', 'Y' ou None (aucune inversion detectee)."""
    x_plus_raw = raw["X+"][0]
    x_minus_mag = polere(*raw["X-"])[0]
    if x_plus_raw < 0 and x_minus_mag > 0:
        return "X"
    y_plus_raw = raw["Y+"][1]
    y_minus_raw = raw["Y-"][1]
    if y_plus_raw < 0 and y_minus_raw > 0:
        return "Y"
    return None


def _find_zb(ech: "SelectedSample", item_rv: int) -> Optional[Measurement]:
    """Equivalent de calcul.f:4014 : mesure cod1='Z',cod2='B' a la meme
    etape que les mesures X (item_rv) - une repetition de la position Z,
    utilisee comme controle de type pTRM-check (verifie qu'aucune
    alteration/evolution n'a eu lieu entre les mesures) - PAS un
    equivalent MagIC standard, specifique a ce protocole. None si
    absente."""
    for m in ech.mesures:
        if m.etape == item_rv and m.cod1 == "Z" and m.cod2 == "B":
            return m
    return None


def _check_trm_evolution(
    positions: Dict[str, Measurement], zb: Measurement, holder: Optional[ArmHolderBackground],
) -> float:
    """Equivalent de calcul.f:4095-4108 : compare la mesure de controle ZB
    a la moyenne Z+/Z- (arnz) et a la demi-difference Z+/Z- normale
    (atrz) - `trmevol = atrzb.z / atrz.z` proche de 1.0 si ZB est
    coherente avec la paire Z+/Z- normale (pas d'alteration). Retourne
    trmevol (ratio, PAS encore convertit en %) ; 1.0 si atrz.z est nul
    (evite une division par zero, absente du Fortran d'origine)."""
    zpx, zpy, zpz = _position_vector(positions, "Z+", 4, holder)
    zmx, zmy, zmz = _position_vector(positions, "Z-", 5, holder)
    arnz_z = (zpz + zmz) / 2
    atrz_z = (zpz - zmz) / 2
    zb_x, zb_y, zb_z = zb.x, zb.y, zb.z
    if holder is not None:
        zb_z -= holder.z[4]  # meme position "Z+" pour la ligne de base (ZB est une repetition de Z+)
    atrzb_z = zb_z - arnz_z
    return atrzb_z / atrz_z if atrz_z else 1.0


@dataclass
class AnisotropyComputation:
    """Trace complete du calcul (pas seulement le resultat final) - pour
    pouvoir afficher toute la sequence comme le fait le Fortran (console
    verbeuse), demande explicite de l'utilisateur : le tenseur BRUT avant
    symetrisation (`raw`) et les diagnostics par position sont le seul
    moyen de verifier qu'un echantillon etait bien oriente lors de
    l'acquisition de la TRM."""
    tensor: AniTensor  # 'A0', symetrise - LE tenseur "officiel" (voir all_tensors)
    all_tensors: List[AniTensor]  # les 15 variantes (A0 en premier), voir _all_ani_variants
    raw: Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]
    positions: Dict[str, Measurement]
    holder_used: bool
    nrm_mean: Tuple[float, float, float]
    position_diags: List[AnisotropyPositionDiag]
    deviation_pct: float  # `deviatTRM`
    swapped_axes: List[str]  # ex. ['X'] si une inversion a ete corrigee
    trm_evolution_pct: Optional[float]  # `TRMevol`, None si pas de ZB trouvee
    zb_used: bool  # Zplus remplace par ZB (evolution > seuil + use_zb_on_evolution)
    pair_nrm: List[AnisotropyPairNrm]  # les 3 estimations "ARN X+X-/Y+Y-/Z+Z-" avant leur moyenne (nrm_mean)
    nrm_mean_diag: AnisotropyPairNrm  # nrm_mean, en moment/dec/inc - meme format que pair_nrm, pour comparaison directe
    # [(position, deviation_deg), ...] au-dela de _POSITION_DEVIATION_
    # WARNING_DEG (45 deg par defaut) - voir _find_misoriented_positions -
    # PAS deja corrige (contrairement a swapped_axes) : juste signale,
    # une deviation de plus de 45 deg ne se laisse pas "corriger" en
    # inversant simplement le signe (contrairement a X/Y ci-dessus, ou le
    # Fortran sait deja quoi faire) - demande explicite utilisateur
    # ("mettre un warning lorsque le code ne correspond pas à la vraie
    # direction enregistrée par l'échantillon").
    misoriented_positions: List[Tuple[str, float]]


def _all_ani_variants(
    id_: str,
    xpx: float, xpy: float, xpz: float, xmx: float, xmy: float, xmz: float,
    ypx: float, ypy: float, ypz: float, ymx: float, ymy: float, ymz: float,
    zpx: float, zpy: float, zpz: float, zmx: float, zmy: float, zmz: float,
) -> List[AniTensor]:
    """Port exact de calcul.f:2691-3256 (bloc d'ecriture des 15 lignes du
    .ANI, dans anisot/anisoauto) - voir le commentaire au-dessus de
    AnisotropyComputation pour le principe. Composantes nommees d'apres
    la convention Fortran (`xpx`=mes(ixp).x, etc, 0-indexe X+,X-,Y+,Y-,
    Z+,Z- dans cet ordre) - deja corrigees de la ligne de base porte-
    echantillon si applicable (memes vecteurs que ceux utilises pour A0).
    Retourne les 15 AniTensor dans l'ordre A0,A+,A-,A1,B1,A2,B2,A3,B3,A4,
    B4,A5,B5,A6,B6."""
    tensors: List[AniTensor] = []

    def add(code2: str, k11: float, k22: float, k33: float, k12: float, k23: float, k13: float) -> None:
        tensors.append(AniTensor(id=id_, code2=code2, k11=k11, k22=k22, k33=k33, k12=k12, k23=k23, k13=k13))

    # Groupe 0 : demi-difference +/- par paire d'axes (A0=symetrise,
    # A+/A- = les 2 triangles non symetrises).
    ani011, ani012, ani013 = (xpx - xmx) / 2, (xpy - xmy) / 2, (xpz - xmz) / 2
    ani021, ani022, ani023 = (ypx - ymx) / 2, (ypy - ymy) / 2, (ypz - ymz) / 2
    ani031, ani032, ani033 = (zpx - zmx) / 2, (zpy - zmy) / 2, (zpz - zmz) / 2
    add("A0", ani011, ani022, ani033, (ani012 + ani021) / 2, (ani032 + ani023) / 2, (ani013 + ani031) / 2)
    add("A+", ani011, ani022, ani033, ani012, ani023, ani013)
    add("A-", ani011, ani022, ani033, ani021, ani032, ani031)

    # Groupe 1/2 : moyenne de la paire X (bruit/remanence parasite commune
    # a X+/X-) soustraite des 3 mesures "+" (A1) et "-" (B... A2, signe
    # inverse) de chaque paire d'axes.
    anrmxx, anrmxy, anrmxz = (xpx + xmx) / 2, (xpy + xmy) / 2, (xpz + xmz) / 2
    ani111, ani112, ani113 = xpx - anrmxx, xpy - anrmxy, xpz - anrmxz
    ani121, ani122, ani123 = ypx - anrmxx, ypy - anrmxy, ypz - anrmxz
    ani131, ani132, ani133 = zpx - anrmxx, zpy - anrmxy, zpz - anrmxz
    ani211, ani212, ani213 = -(xmx - anrmxx), -(xmy - anrmxy), -(xmz - anrmxz)
    ani221, ani222, ani223 = -(ymx - anrmxx), -(ymy - anrmxy), -(ymz - anrmxz)
    ani231, ani232, ani233 = -(zmx - anrmxx), -(zmy - anrmxy), -(zmz - anrmxz)
    add("A1", ani111, ani122, ani133, ani112, ani123, ani113)
    add("B1", ani111, ani122, ani133, ani121, ani132, ani131)
    add("A2", ani211, ani222, ani233, ani212, ani223, ani213)
    add("B2", ani211, ani222, ani233, ani221, ani232, ani231)

    # Groupe 3/4 : idem, moyenne de la paire Y.
    anrmyx, anrmyy, anrmyz = (ypx + ymx) / 2, (ypy + ymy) / 2, (ypz + ymz) / 2
    ani311, ani312, ani313 = xpx - anrmyx, xpy - anrmyy, xpz - anrmyz
    ani321, ani322, ani323 = ypx - anrmyx, ypy - anrmyy, ypz - anrmyz
    ani331, ani332, ani333 = zpx - anrmyx, zpy - anrmyy, zpz - anrmyz
    ani411, ani412, ani413 = -(xmx - anrmyx), -(xmy - anrmyy), -(xmz - anrmyz)
    ani421, ani422, ani423 = -(ymx - anrmyx), -(ymy - anrmyy), -(ymz - anrmyz)
    ani431, ani432, ani433 = -(zmx - anrmyx), -(zmy - anrmyy), -(zmz - anrmyz)
    add("A3", ani311, ani322, ani333, ani312, ani323, ani313)
    add("B3", ani311, ani322, ani333, ani321, ani332, ani331)
    add("A4", ani411, ani422, ani433, ani412, ani423, ani413)
    add("B4", ani411, ani422, ani433, ani421, ani432, ani431)

    # Groupe 5/6 : idem, moyenne de la paire Z.
    anrmzx, anrmzy, anrmzz = (zpx + zmx) / 2, (zpy + zmy) / 2, (zpz + zmz) / 2
    ani511, ani512, ani513 = xpx - anrmzx, xpy - anrmzy, xpz - anrmzz
    ani521, ani522, ani523 = ypx - anrmzx, ypy - anrmzy, ypz - anrmzz
    ani531, ani532, ani533 = zpx - anrmzx, zpy - anrmzy, zpz - anrmzz
    ani611, ani612, ani613 = -(xmx - anrmzx), -(xmy - anrmzy), -(xmz - anrmzz)
    ani621, ani622, ani623 = -(ymx - anrmzx), -(ymy - anrmzy), -(ymz - anrmzz)
    ani631, ani632, ani633 = -(zmx - anrmzx), -(zmy - anrmzy), -(zmz - anrmzz)
    add("A5", ani511, ani522, ani533, ani512, ani523, ani513)
    add("B5", ani511, ani522, ani533, ani521, ani532, ani531)
    add("A6", ani611, ani622, ani633, ani612, ani623, ani613)
    add("B6", ani611, ani622, ani633, ani621, ani632, ani631)

    return tensors


def compute_anisotropy_tensor(
    ech: "SelectedSample", holder: Optional[ArmHolderBackground] = None,
    use_zb_on_evolution: bool = False,
    positions: Optional[Dict[str, Measurement]] = None,
) -> Optional[AnisotropyComputation]:
    """Tenseur 'A0' (equivalent de la branche de base d'anisot/anisoauto,
    calcul.f:4005-4595) : detecte les 6 positions (detect_six_positions) -
    ou reprend celles fournies via `positions` (saisie manuelle, "n" a
    "automated recognition of the 6 steps ?", calcul.f:2450-2459 - meme
    calcul ensuite, seule l'identification des positions differe),
    soustrait la ligne de base porte-echantillon si fournie, calcule la
    NRM residuelle moyenne et les diagnostics par position (une
    anisotropie basee sur des TRM PARTIELLES laisse souvent une part de
    NRM non remplacee - demande explicite de l'utilisateur), detecte et
    corrige une inversion X+/X- ou Y+/Y- (specimen mal oriente pendant
    l'acquisition - calcul.f:4386-4427), verifie une eventuelle
    alteration via une mesure de controle ZB si presente
    (_check_trm_evolution ; substitue Z+ par ZB si l'evolution depasse
    5% ET `use_zb_on_evolution`), puis calcule le tenseur BRUT
    (antisymetrique, `raw`) et le tenseur symetrique final. Retourne None
    si les 6 positions ne sont pas identifiees.

    `use_zb_on_evolution` correspond au choix demande UNE FOIS pour tout
    le lot dans le Fortran ("using ZB instead of R for a > + or - 5%
    evolution ? y/N"), pas par echantillon."""
    if positions is None:
        positions = detect_six_positions(ech)
        if positions is None:
            return None
    else:
        positions = dict(positions)  # copie - peut etre mutee par le swap X/Y ci-dessous

    swapped_axes: List[str] = []
    for _ in range(4):  # borne de securite - converge en 1-2 iterations en pratique
        nrm_mean, diags, deviation_pct, raw_components, pair_nrm, nrm_mean_diag = _nrm_mean_and_diag(
            positions, holder, ech.norme, ech.vol)
        axis = _check_position_inversion(raw_components)
        if axis is None:
            break
        if axis == "X":
            positions["X+"], positions["X-"] = positions["X-"], positions["X+"]
        else:
            positions["Y+"], positions["Y-"] = positions["Y-"], positions["Y+"]
        swapped_axes.append(axis)

    trm_evolution_pct: Optional[float] = None
    zb_used = False
    item_rv = positions["X+"].etape
    zb = _find_zb(ech, item_rv)
    if zb is not None:
        trmevol = _check_trm_evolution(positions, zb, holder)
        trm_evolution_pct = (trmevol - 1.0) * 100.0
        if use_zb_on_evolution and (trmevol > 1.05 or trmevol < 0.95):
            reclassified = detect_six_positions(ech, force_zplus_label=zb.cod1 + zb.cod2)
            if reclassified is not None:
                positions = reclassified
                zb_used = True
                nrm_mean, diags, deviation_pct, raw_components, pair_nrm, nrm_mean_diag = _nrm_mean_and_diag(
                    positions, holder, ech.norme, ech.vol)

    xpx, xpy, xpz = _position_vector(positions, "X+", 0, holder)
    xmx, xmy, xmz = _position_vector(positions, "X-", 1, holder)
    ypx, ypy, ypz = _position_vector(positions, "Y+", 2, holder)
    ymx, ymy, ymz = _position_vector(positions, "Y-", 3, holder)
    zpx, zpy, zpz = _position_vector(positions, "Z+", 4, holder)
    zmx, zmy, zmz = _position_vector(positions, "Z-", 5, holder)

    ani011, ani012, ani013 = (xpx - xmx) / 2, (xpy - xmy) / 2, (xpz - xmz) / 2
    ani021, ani022, ani023 = (ypx - ymx) / 2, (ypy - ymy) / 2, (ypz - ymz) / 2
    ani031, ani032, ani033 = (zpx - zmx) / 2, (zpy - zmy) / 2, (zpz - zmz) / 2

    all_tensors = _all_ani_variants(
        ech.id,
        xpx, xpy, xpz, xmx, xmy, xmz,
        ypx, ypy, ypz, ymx, ymy, ymz,
        zpx, zpy, zpz, zmx, zmy, zmz,
    )
    tensor = all_tensors[0]  # 'A0', identique a l'ancien calcul direct ci-dessus
    raw = ((ani011, ani012, ani013), (ani021, ani022, ani023), (ani031, ani032, ani033))
    misoriented_positions = _find_misoriented_positions(positions, holder)
    return AnisotropyComputation(
        tensor=tensor, all_tensors=all_tensors, raw=raw, positions=positions,
        holder_used=holder is not None,
        nrm_mean=nrm_mean, position_diags=diags, deviation_pct=deviation_pct,
        swapped_axes=swapped_axes, trm_evolution_pct=trm_evolution_pct, zb_used=zb_used,
        pair_nrm=pair_nrm, nrm_mean_diag=nrm_mean_diag,
        misoriented_positions=misoriented_positions,
    )


# Format .pmagani (tabulations, entete de colonnes explicite) - remplace
# l'ancien format .ANI list-directed Fortran (espaces, sans entete) -
# demande explicite utilisateur ("can we write the name of the ani
# extension as .pmagani", "What will be the most complete .pmagani file
# as a companion to .prmag and pmagres"). cin/caz/dip/str NE SONT PLUS
# ecrits ici - demande explicite utilisateur ("the codes cin caz dip str
# are now not needed if the pmagani is linked to the prmag file") : ce
# sont exactement les memes champs (memes noms) que ceux deja stockes
# par specimen dans .prmag (voir testlect.Pmag.cin/caz/dip/str_), donc
# une pure duplication maintenant que .pmagani est systematiquement
# associe a son .prmag (meme nom de base, jointure par specimen - voir
# ani_path_for) ; etape est conserve (diagnostic, PAS dans .prmag).
# n_positions/sigma/ftest/ftest12/ftest23 sont NOUVEAUX - "n.d" si non
# disponibles (le calcul natif STARpaleomag_Py, jackknife geometrique, ne calcule
# pas les statistiques de Hext ; seul le chemin PmagPy - voir
# anisotropy_magic.compute_aarm_pmagpy - les fournit).
_PMAGANI_HEADER = [
    "specimen", "code2", "etape",
    "k11", "k22", "k33", "k12", "k23", "k13", "s(SI*1e-5)",
    "n_positions", "sigma", "ftest", "ftest12", "ftest23", "quality", "export", "info",
]
# Largeur d'un fichier ecrit AVANT l'ajout de la colonne "export" (voir
# AniTensor.export) - MEME principe de retro-compatibilite que
# _PMAGANI_MEAN_HEADER_LEGACY_LEN plus bas ; un tel fichier reste lu
# "export=Y" par defaut (voir _read_pmagani_tensors).
_PMAGANI_HEADER_LEGACY_LEN = len(_PMAGANI_HEADER) - 1
# `quality` ('g'/'b'/n.d) : verdict PmagPy (Hext F-test, aniso_ftest_
# quality) persiste comme colonne A PART ENTIERE (pas seulement noye dans
# le texte libre `info`) - necessaire pour pouvoir le RELIRE de facon
# fiable (calcul.read_ani_tensor) plutot que de re-parser une phrase en
# texte libre - demande explicite utilisateur (prompt "inverse anisotropy
# correction anyway ?" lors de la correction paleointensite, qui a besoin
# de savoir SANS ambiguite si le tenseur relu depuis le fichier est
# satisfactory ou non).

# Deux sections dans le MEME fichier .pmagani, MEME principe que
# .pmagres (_SPECIMEN_HEADER/_MEAN_HEADER plus haut dans ce module) -
# demande explicite utilisateur ("can we have the pmagani as pmagres
# with two blocks one at a specimen level, the next one at site level
# with both blocks having their own set of parameters") : colonnes
# specimen (_PMAGANI_HEADER, k11..k13 bruts) INCHANGEES, colonnes site
# mean (_PMAGANI_MEAN_HEADER, axes propres k1>=k2>=k3 + dec/inc/ellipse
# de confiance, PAS un simple duplicata du schema specimen). Un ancien
# fichier .pmagani (avant ce changement, jamais de marqueur
# _ANI_MEAN_HEADER) reste lisible tel quel - _read_pmagani_tensors le
# traite ENTIEREMENT comme section specimen, exactement comme avant.
_ANI_SPECIMEN_HEADER = "#specimen tensor results"
_ANI_MEAN_HEADER = "#site mean tensor results"

_PMAGANI_MEAN_HEADER = [
    "site", "code2", "n",
    "k1", "k2", "k3",
    "dec1", "inc1", "dec2", "inc2", "dec3", "inc3",
    "alpha1_1", "alpha2_1", "alpha1_2", "alpha2_2", "alpha1_3", "alpha2_3",
    "P", "T", "L", "F", "Pprim",
    "tilt_correction", "export", "info",
]
# nombre de colonnes d'un ancien fichier (avant l'ajout de
# tilt_correction, colonne dediee) - demande explicite utilisateur ("oui
# ajouter une colonne avant info", suite a "in the file pmagani, there is
# no information about the IS/TC" : l'info existait mais SEULEMENT en
# texte libre dans `info` - "tilt_correction: N" - jamais une vraie
# colonne du tableau, une decision deliberee au depart cote AMS_Py
# ("encoder ceci dans info evite de re-modifier le schema deja valide
# cote STARpaleomag_Py/calcul") revenue explicitement ici). Les DEUX
# largeurs restent lisibles en lecture (voir _read_pmagani_mean_tensors) -
# un fichier ecrit AVANT ce changement n'a que 24 colonnes (info en
# dernier), un fichier ecrit APRES en a 25 (tilt_correction avant info).
# TROIS paliers desormais (meme principe que ams_selection._PMAGANI_
# MEAN_HEADER_LEGACY_LEN/_NO_EXPORT_LEN cote AMS_Py, MEME fichier
# partage) : la colonne "export" (voir AniMeanTensor.export) a ete
# ajoutee APRES tilt_correction, donc un fichier peut avoir 24
# colonnes (ni l'un ni l'autre), 25 (tilt_correction seul) ou 26 (les
# deux, format actuel) - _read_pmagani_mean_tensors verifie la plus
# longue en premier.
_PMAGANI_MEAN_HEADER_NO_EXPORT_LEN = len(_PMAGANI_MEAN_HEADER) - 1
_PMAGANI_MEAN_HEADER_LEGACY_LEN = len(_PMAGANI_MEAN_HEADER) - 2

# `s` (susceptibilite scalaire) est en SI*1e-5 (ex. 7261 = 0.07261 SI) -
# convention Bartington ("since the 80s ... measure in 1e-5 SI assuming a
# volume of 10cc ... could read between 0 and 9999"), gardee par
# l'utilisateur pour pouvoir comparer Bartington et AGICO - demande
# explicite utilisateur ("We should add in the header that s is in
# SI*1e-5"). Ecrite en commentaire d'entete par _insert_pmagani_line
# (et cote AMS_Py par import_legacy_ani/import_asc_file).
_PMAGANI_UNITS_NOTE = "# s(SI*1e-5): susceptibility in units of 1e-5 SI (Bartington/AGICO convention, e.g. 7261 = 0.07261 SI)\n"


def _fmt_pmagani_stat(v: Optional[float]) -> str:
    return "n.d" if v is None else f"{v:.6g}"


def _parse_pmagani_stat(v: str) -> Optional[float]:
    v = v.strip()
    if not v or v == "n.d":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _insert_pmagani_line(path: str, line: str, is_mean: bool) -> None:
    """Insere UNE ligne .pmagani a la bonne section - meme discipline a
    deux sections que calcul.archivres pour .pmagres (voir le commentaire
    au-dessus de _ANI_SPECIMEN_HEADER) : une ligne specimen va en fin de
    la section specimen (juste avant l'en-tete '#site mean tensor
    results' si elle existe deja, en sautant les lignes vides qui la
    precedent) ; une ligne site mean va TOUJOURS en fin de fichier
    (section mean, creee au besoin avec son propre en-tete). REMPLACE
    l'ancien `_ensure_pmagani_header` (simple `open(path, "a")` en
    aveugle) qui aurait, une fois une section mean ajoutee, laisse un
    NOUVEAU tenseur specimen s'ecrire APRES elle - cassant les deux
    sections."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        lines = [
            "# pmagani v2 - companion of .prmag/.pmagres, join key = specimen "
            "(specimen section) / site (site mean section)",
            _PMAGANI_UNITS_NOTE.rstrip("\n"),
            _ANI_SPECIMEN_HEADER,
            "\t".join(_PMAGANI_HEADER),
        ]
    else:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()

    mean_header_idx = next(
        (i for i, l in enumerate(lines) if l.strip() == _ANI_MEAN_HEADER), None)

    if is_mean:
        if mean_header_idx is None:
            if lines and lines[-1].strip():
                lines.append("")
            lines.append(_ANI_MEAN_HEADER)
            lines.append("\t".join(_PMAGANI_MEAN_HEADER))
        lines.append(line.rstrip("\n"))
    elif mean_header_idx is None:
        lines.append(line.rstrip("\n"))
    else:
        insert_idx = mean_header_idx
        while insert_idx > 0 and not lines[insert_idx - 1].strip():
            insert_idx -= 1
        lines.insert(insert_idx, line.rstrip("\n"))

    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


def _format_pmagani_line(
    specimen_id: str, tensor: AniTensor, etape: int,
    zplus_label: str, zminus_label: str, trm_evolution_pct: float, deviation_pct: float,
    info_text: Optional[str] = None,
) -> str:
    """`info_text` : override total du texte informatif habituel (etape/
    positions/evolution TRM) - utilise par write_ani_tensors_from_magic
    (tenseur DEJA calcule par une contribution MagIC importee, aucune des
    6 positions/etape STARpaleomag_Py n'est disponible ici) ; None (par defaut)
    reconstruit le texte usuel a partir de etape/zplus_label/zminus_label/
    trm_evolution_pct/deviation_pct comme avant."""
    if info_text is None:
        info_text = (
            f"TRM evo: {trm_evolution_pct:5.1f} deviation:{deviation_pct:5.1f}"
            f"  steps: X+ X- Y+ Y- {zplus_label} {zminus_label}"
        )
    if tensor.quality == "g":
        info_text += " - PmagPy Hext F-test: significant anisotropy (satisfactory)"
    elif tensor.quality == "b":
        info_text += " - PmagPy Hext F-test: NOT significant (not satisfactory)"
    info = f'"{info_text}"'
    n_positions = tensor.n_positions if tensor.n_positions is not None else 6
    fields = [
        specimen_id, tensor.code2, str(etape),
        f"{tensor.k11:.6E}", f"{tensor.k22:.6E}", f"{tensor.k33:.6E}",
        f"{tensor.k12:.6E}", f"{tensor.k23:.6E}", f"{tensor.k13:.6E}",
        "1.00000", str(n_positions),
        _fmt_pmagani_stat(tensor.sigma), _fmt_pmagani_stat(tensor.ftest),
        _fmt_pmagani_stat(tensor.ftest12), _fmt_pmagani_stat(tensor.ftest23),
        tensor.quality or "n.d",
        tensor.export or "Y",
        info,
    ]
    return "\t".join(fields) + "\n"


def _read_pmagani_tensors(path: str) -> List[AniTensor]:
    """Colonnes : specimen,code2,etape,k11,k22,k33,k12,k23,k13,s,
    n_positions,sigma,ftest,ftest12,ftest23,quality,info (voir
    _PMAGANI_HEADER) - PLUS de cin/caz/dip/str (voir commentaire au-dessus
    de _PMAGANI_HEADER) : ces champs viennent maintenant du .prmag associe
    (jointure par specimen), pas de ce fichier.

    S'ARRETE des que l'en-tete '#site mean tensor results' est rencontre
    (voir _read_pmagani_mean_tensors pour cette section) - un fichier
    ANTERIEUR a l'ajout de cette 2e section (jamais ce marqueur) est lu
    exactement comme avant, en entier."""
    out: List[AniTensor] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if stripped == _ANI_MEAN_HEADER:
                    break
                continue
            parts = line.split("\t")
            if not parts or parts[0] == "specimen":
                continue  # ligne d'entete
            if len(parts) < 9:
                continue
            try:
                k11, k22, k33, k12, k23, k13 = (float(p) for p in parts[3:9])
            except ValueError:
                continue
            n_positions = None
            if len(parts) > 10:
                try:
                    n_positions = int(float(parts[10]))
                except ValueError:
                    n_positions = None
            sigma = _parse_pmagani_stat(parts[11]) if len(parts) > 11 else None
            ftest = _parse_pmagani_stat(parts[12]) if len(parts) > 12 else None
            ftest12 = _parse_pmagani_stat(parts[13]) if len(parts) > 13 else None
            ftest23 = _parse_pmagani_stat(parts[14]) if len(parts) > 14 else None
            quality = None
            if len(parts) > 15 and parts[15].strip() in ("g", "b"):
                quality = parts[15].strip()
            # export : colonne ajoutee cote AMS_Py (voir AniTensor.export) -
            # absente (fichier plus ancien que _PMAGANI_HEADER_LEGACY_LEN)
            # => "Y" par defaut, meme retro-compatibilite que n_positions/
            # sigma/ftest* ci-dessus.
            export = "Y"
            if len(parts) > 16 and parts[16].strip() in ("Y", "N"):
                export = parts[16].strip()
            out.append(AniTensor(
                id=parts[0], code2=parts[1], k11=k11, k22=k22, k33=k33,
                k12=k12, k23=k23, k13=k13,
                n_positions=n_positions, sigma=sigma, ftest=ftest,
                ftest12=ftest12, ftest23=ftest23, quality=quality,
                export=export,
            ))
    return out


def _format_pmagani_mean_line(mean: AniMeanTensor) -> str:
    info = f'"{mean.info}"' if mean.info else '""'
    fields = [
        mean.id, mean.code2, str(mean.n),
        f"{mean.k1:.6E}", f"{mean.k2:.6E}", f"{mean.k3:.6E}",
        f"{mean.dec1:.2f}", f"{mean.inc1:.2f}", f"{mean.dec2:.2f}", f"{mean.inc2:.2f}",
        f"{mean.dec3:.2f}", f"{mean.inc3:.2f}",
        f"{mean.alpha1_1:.3f}", f"{mean.alpha2_1:.3f}",
        f"{mean.alpha1_2:.3f}", f"{mean.alpha2_2:.3f}",
        f"{mean.alpha1_3:.3f}", f"{mean.alpha2_3:.3f}",
        _fmt_pmagani_stat(mean.P), _fmt_pmagani_stat(mean.T),
        _fmt_pmagani_stat(mean.L), _fmt_pmagani_stat(mean.F), _fmt_pmagani_stat(mean.Pprim),
        mean.tilt_correction, mean.export or "Y", info,
    ]
    return "\t".join(fields) + "\n"


def _read_pmagani_mean_tensors(path: str) -> List[AniMeanTensor]:
    """Lit UNIQUEMENT la section '#site mean tensor results' - voir
    _PMAGANI_MEAN_HEADER pour les colonnes. [] si le fichier n'a jamais
    cette section (fichier ecrit avant son introduction, ou qui n'a
    encore aucune moyenne de site archivee)."""
    out: List[AniMeanTensor] = []
    in_mean_section = False
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if stripped == _ANI_MEAN_HEADER:
                    in_mean_section = True
                continue
            if not in_mean_section:
                continue
            parts = line.split("\t")
            if not parts or parts[0] == "site":
                continue  # ligne d'entete
            # Retro-compatibilite sur le NOMBRE de colonnes (meme principe
            # que _parse_mean_line pour .pmagres) : TROIS paliers, verifies
            # du plus long au plus court (voir commentaire au-dessus de
            # _PMAGANI_MEAN_HEADER_LEGACY_LEN) - le format actuel a
            # tilt_correction ET export avant info, le palier intermediaire
            # a tilt_correction mais pas export, le plus ancien n'a ni
            # l'un ni l'autre (info en dernier).
            if len(parts) >= len(_PMAGANI_MEAN_HEADER):
                tilt_correction, export, info_raw = parts[23].strip(), parts[24].strip(), parts[25]
            elif len(parts) >= _PMAGANI_MEAN_HEADER_NO_EXPORT_LEN:
                tilt_correction, export, info_raw = parts[23].strip(), "Y", parts[24]
            elif len(parts) >= _PMAGANI_MEAN_HEADER_LEGACY_LEN:
                tilt_correction, export, info_raw = "", "Y", parts[23]
            else:
                continue
            if export not in ("Y", "N"):
                export = "Y"
            try:
                out.append(AniMeanTensor(
                    id=parts[0], code2=parts[1], n=int(float(parts[2])),
                    k1=float(parts[3]), k2=float(parts[4]), k3=float(parts[5]),
                    dec1=float(parts[6]), inc1=float(parts[7]),
                    dec2=float(parts[8]), inc2=float(parts[9]),
                    dec3=float(parts[10]), inc3=float(parts[11]),
                    alpha1_1=float(parts[12]), alpha2_1=float(parts[13]),
                    alpha1_2=float(parts[14]), alpha2_2=float(parts[15]),
                    alpha1_3=float(parts[16]), alpha2_3=float(parts[17]),
                    P=_parse_pmagani_stat(parts[18]), T=_parse_pmagani_stat(parts[19]),
                    L=_parse_pmagani_stat(parts[20]), F=_parse_pmagani_stat(parts[21]),
                    Pprim=_parse_pmagani_stat(parts[22]),
                    tilt_correction=tilt_correction,
                    export=export,
                    info=info_raw.strip('"'),
                ))
            except ValueError:
                continue
    return out


def write_ani_mean_tensor(path: str, mean: AniMeanTensor) -> None:
    """Ajoute UNE moyenne de site au fichier .pmagani, dans la section
    '#site mean tensor results' (creee au besoin, toujours en fin de
    fichier - voir _insert_pmagani_line) - demande explicite utilisateur
    ("can we have the pmagani as pmagres with two blocks one at a
    specimen level, the next one at site level with both blocks having
    their own set of parameters")."""
    _insert_pmagani_line(path, _format_pmagani_mean_line(mean), is_mean=True)


def read_ani_mean_tensor(path: str, site: str, code2: str) -> Optional[AniMeanTensor]:
    """Retrouve la moyenne de site `code2` (ex. 'A0') pour `site` (ex.
    '12TU74') - meme convention que read_ani_tensor (recherche
    insensible a la casse sur l'id). None si absente/fichier introuvable."""
    if not os.path.exists(path):
        return None
    site_upper = site.strip().upper()
    for t in _read_pmagani_mean_tensors(path):
        if t.id.strip().upper() == site_upper and t.code2 == code2:
            return t
    return None


def _read_ani_tensors_legacy(path: str) -> List[AniTensor]:
    """Lit un ancien fichier .ANI (format list-directed Fortran, sans
    entete : `D id cin caz dip str etape code2 k11 k22 k33 k12 k23 k13 s
    [info]`) - lecture directe seulement (voir read_ani_tensor, dispatch
    par extension) ; la CONVERSION explicite vers .pmagani se fait
    desormais cote AMS_Py uniquement (ams_selection.import_legacy_ani) -
    demande explicite utilisateur ("in starmac_Py, you can delete the
    import legacy .ANI file, it is better to handle this in AMS_Py")."""
    out: List[AniTensor] = []
    with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 14 or parts[0] != "D":
                continue
            try:
                k11, k22, k33, k12, k23, k13 = (float(p) for p in parts[8:14])
            except ValueError:
                continue
            out.append(AniTensor(id=parts[1], code2=parts[7], k11=k11, k22=k22, k33=k33,
                                  k12=k12, k23=k23, k13=k13))
    return out


def write_ani_tensor(
    path: str, ech: "SelectedSample", tensor: AniTensor,
    positions: Optional[Dict[str, Measurement]] = None,
    trm_evolution_pct: float = 0.0, deviation_pct: float = 0.0,
) -> None:
    """Ajoute une ligne au fichier .pmagani (voir _format_pmagani_line ;
    format tabule avec entete, relu par `read_ani_tensor`/
    `_read_pmagani_tensors`). L'etape et les labels Z+/Z- (pour le texte
    informatif) viennent de `positions` (voir AnisotropyComputation.
    positions - passer ce qui a deja ete detecte evite un second appel a
    detect_six_positions) ; si non fourni, redetecte depuis `ech`.
    `trm_evolution_pct`/`deviation_pct` : voir AnisotropyComputation (0.0
    par defaut si absents, ex. pas de mesure ZB trouvee). Pour ecrire les
    15 variantes d'un coup, voir write_ani_tensors."""
    if positions is None:
        positions = detect_six_positions(ech)
    etape = positions["X+"].etape if positions else 0
    zplus_label = positions["Z+"].cod1 + positions["Z+"].cod2 if positions else ""
    zminus_label = positions["Z-"].cod1 + positions["Z-"].cod2 if positions else ""
    line = _format_pmagani_line(ech.id, tensor, etape, zplus_label, zminus_label, trm_evolution_pct, deviation_pct)
    _insert_pmagani_line(path, line, is_mean=False)


def write_ani_tensors(
    path: str, ech: "SelectedSample", tensors: List[AniTensor],
    positions: Optional[Dict[str, Measurement]] = None,
    trm_evolution_pct: float = 0.0, deviation_pct: float = 0.0,
) -> None:
    """Ecrit UNE ligne par tensor de `tensors` (typiquement les 15
    variantes de AnisotropyComputation.all_tensors - A0/A+/A-/A1/B1/A2/
    B2/A3/B3/A4/B4/A5/B5/A6/B6) - demande explicite utilisateur
    ("implement the 15 outputs for the anisotropy tensor with the A0
    code being the major one"). Equivalent du Fortran, qui ecrit les 15
    lignes d'affilee des que l'utilisateur confirme UNE FOIS "save in
    file .ANI: Y/n" (calcul.f:3038-3255) - meme etape/labels Z+/Z-/info
    pour toutes (calcules une seule fois, comme dans le Fortran), seul
    `tensor.code2` change d'une ligne a l'autre."""
    if positions is None:
        positions = detect_six_positions(ech)
    etape = positions["X+"].etape if positions else 0
    zplus_label = positions["Z+"].cod1 + positions["Z+"].cod2 if positions else ""
    zminus_label = positions["Z-"].cod1 + positions["Z-"].cod2 if positions else ""
    for tensor in tensors:
        _insert_pmagani_line(path, _format_pmagani_line(
            ech.id, tensor, etape, zplus_label, zminus_label, trm_evolution_pct, deviation_pct),
            is_mean=False)


def write_ani_tensors_from_magic_rows(path: str, rows: List[dict]) -> int:
    """Archive dans le .pmagani UN tenseur DEJA calcule par une
    contribution MagIC importee (voir extract_magic.magic_anisotropy_rows,
    qui produit `rows`) - demande explicite utilisateur ("lors de
    l'importation des fichiers Magic, archiver dans le fichier pmagani
    les donnees de tenseurs d'anisotropie avec le code A0 si c'est un
    tenseur de TRM, N0 pour l'AMS, AA pour l'ARM").

    A la difference de write_ani_tensor/write_ani_tensors (calcul natif
    STARpaleomag_Py, 6 positions X+/X-/Y+/Y-/Z+/Z- toujours disponibles), il n'y a
    ici ni etape ni positions a rapporter - `info_text` remplace le texte
    "steps: X+ X- ..." habituel par une mention explicite de la
    provenance, pour qu'une ligne "AA"/"A0"/"N0" issue d'un import MagIC
    ne soit jamais confondue, en relisant le fichier, avec un tenseur
    recalcule nativement a partir des 6 mesures. `etape` est ecrit a 0
    (aucune mesure STARpaleomag_Py associee - diagnostic seulement, voir
    _PMAGANI_HEADER). Retourne le nombre de lignes ecrites (une par
    element de `rows`, PAS un total de specimens : plusieurs protocoles
    d'anisotropie sur le meme specimen produisent plusieurs lignes)."""
    if not rows:
        return 0
    for row in rows:
        tensor = AniTensor(
            id=row["specimen"], code2=row["code2"],
            k11=row["k11"], k22=row["k22"], k33=row["k33"],
            k12=row["k12"], k23=row["k23"], k13=row["k13"],
            n_positions=row.get("n_positions"), sigma=row.get("sigma"),
            ftest=row.get("ftest"), ftest12=row.get("ftest12"), ftest23=row.get("ftest23"),
            quality=row.get("quality"),
        )
        _insert_pmagani_line(path, _format_pmagani_line(
            row["specimen"], tensor, 0, "", "", 0.0, 0.0,
            info_text="tensor imported as-is from a MagIC contribution (not recomputed by STARpaleomag_Py)",
        ), is_mean=False)
    return len(rows)


# ---------------------------------------------------------------------------
# inverseani / inverse : applique l'inverse d'un tenseur d'anisotropie
# (lu dans un fichier .ANI) a toutes les mesures d'un echantillon -
# plotpaleoint2.f:1803-1918 + calcul.f:3311-3327 (`inverse`)
# ---------------------------------------------------------------------------

def inverse_symmetric_3x3(
    a: float, b: float, c: float, d: float, e: float, f: float,
    x: float, y: float, z: float,
) -> Tuple[float, float, float]:
    """Equivalent de `inverse` : applique l'inverse de la matrice symetrique
    [[a,b,c],[b,d,e],[c,e,f]] au vecteur (x,y,z). Retourne (x,y,z) inchange
    si le determinant est nul (matrice non inversible)."""
    det = a * d * f - a * e * e - b * b * f + 2 * b * c * e - c * c * d
    if det == 0:
        return x, y, z
    a11, a12, a13 = (d * f - e * e) / det, -(b * f - c * e) / det, (b * e - c * d) / det
    a22, a23 = (a * f - c * c) / det, -(a * e - b * c) / det
    a33 = (a * d - b * b) / det
    x1 = a11 * x + a12 * y + a13 * z
    y1 = a12 * x + a22 * y + a23 * z
    z1 = a13 * x + a23 * y + a33 * z
    return x1, y1, z1


def compute_anicor_factor(
    tensor: AniTensor, direction_specimen_frame: Tuple[float, float, float],
) -> float:
    """Equivalent de `anicor` (plotpaleoint2.f:1686-1785) : facteur de
    correction d'anisotropie `fcor` a appliquer a l'intensite (H*fcor)
    d'apres le tenseur TRM 'A0' et la direction ANCREE (repere specimen,
    PAS encore tournee selon l'orientation - `xnrm_nocor` du Fortran).

    NON PARFAITEMENT VERIFIE : sur l'echantillon reel 06A (tenseur .ANI
    reel + direction ancree calculee), fcor=0.892 est obtenu contre 0.887
    dans le transcript Fortran reel (~0.5% d'ecart) - vraisemblablement
    une petite difference de precision numerique dans le calcul du
    vecteur propre principal (ACP) entre numpy et l'algorithme Fortran
    d'origine, la formule elle-meme est transcrite telle quelle depuis le
    source. A surveiller si des resultats bien plus divergents sont
    observes sur d'autres echantillons."""
    rnorm = (tensor.k11 + tensor.k22 + tensor.k33) / 3.0
    if rnorm == 0:
        return 0.0
    a, b, c = tensor.k11 / rnorm, tensor.k12 / rnorm, tensor.k13 / rnorm
    d, e, f = tensor.k22 / rnorm, tensor.k23 / rnorm, tensor.k33 / rnorm
    xx, yy, zz = direction_specimen_frame
    x1, y1, z1 = inverse_symmetric_3x3(a, b, c, d, e, f, xx, yy, zz)
    norm1 = math.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
    if norm1 == 0:
        return 0.0
    x1, y1, z1 = x1 / norm1, y1 / norm1, z1 / norm1
    rhlab = math.sqrt(c * c + e * e + f * f)
    x2 = a * x1 + b * y1 + c * z1
    y2 = d * y1 + e * z1
    z2 = f * z1
    ranc = math.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
    return rhlab / ranc if ranc else 0.0


def read_ani_tensor(path: str, sample_id: str, code2: str) -> Optional[AniTensor]:
    """Recherche un tenseur dans un fichier .pmagani (nouveau format
    tabule, voir _read_pmagani_tensors) OU un ancien .ANI (list-directed
    Fortran, voir _read_ani_tensors_legacy) - dispatch par EXTENSION, pas
    par contenu : les deux formats restent lisibles, demande explicite
    utilisateur ("can we import old style .ANI in these pmagani style").
    VERIFIE contre un vrai fichier .ANI (Miriam_2025b/SanJuan_Pmag.ANI) -
    parsing confirme correct sur des lignes reelles produites par
    anisoauto (voir compute_anisotropy_tensor).

    Comparaison d'id INSENSIBLE A LA CASSE (`.upper()` des deux cotes) -
    meme convention que PARTOUT ailleurs dans l'appli pour les id de
    specimen (selection.select_samples/select_samples_by_site/...,
    toujours `.upper()`) - bug reel confirme sur donnees reelles
    (Miriam_2026/SanJuan_Pmag) : le .prmag/.pmagres connait le specimen
    "18A" (majuscule) mais le .pmagani (issu d'un import .ANI ancien)
    l'a enregistre "18a" (minuscule) - la comparaison exacte precedente
    ne trouvait donc jamais le tenseur A0 pourtant bien present, et la
    correction d'anisotropie en paleointensite echouait silencieusement
    (aucune note, aucune erreur) - demande explicite utilisateur ("la
    correction d'anisotropie ne se fait pas... test sur sample 18A")."""
    if not os.path.exists(path):
        return None
    ext = os.path.splitext(path)[1].lower()
    tensors = _read_pmagani_tensors(path) if ext == ".pmagani" else _read_ani_tensors_legacy(path)
    sample_id_upper = sample_id.upper()
    for t in tensors:
        if t.id.upper() == sample_id_upper and t.code2 == code2:
            return t
    return None


def read_all_ani_tensors(path: str) -> List[AniTensor]:
    """TOUS les tenseurs specimen d'un .pmagani/.ANI (meme dispatch par
    extension que read_ani_tensor) - demande explicite utilisateur
    ("il peut y avoir des specimens dans pmagani qui n'ont pas de
    mesures dans prmag mais il faut aussi les exporter") : permet de
    retrouver, cote export MagIC, les specimens connus SEULEMENT via
    leur mesure d'anisotropie (jamais dans self.donnees/.prmag - AMS
    mesuree sur un sous-ensemble different, ou import .pmagani separe)
    plutot que de se limiter a `self.selection`, qui ne connait que ce
    qui vient du .prmag."""
    if not os.path.exists(path):
        return []
    ext = os.path.splitext(path)[1].lower()
    return _read_pmagani_tensors(path) if ext == ".pmagani" else _read_ani_tensors_legacy(path)


def find_orphan_ani_specimens(pmag_list: List["Pmag"], ani_path: str) -> List[str]:
    """Ids d'un .pmagani/.ANI SANS AUCUNE trace (id insensible a la
    casse) dans `pmag_list` (self.donnees) - demande explicite
    utilisateur ("en theorie pmagani est lie a prmag. Il faut interdire
    dans pmagani les specimens qui ne sont pas présents dans prmag").
    Un .pmagani decrit des specimens du .prmag COMPAGNON (meme nom de
    base, voir ani_path_for) : un id qui n'y figure meme pas signale un
    vrai desaccord entre les deux fichiers (fichiers desassortis, faute
    de frappe, .pmagani d'un autre jeu de donnees copie/renomme par
    erreur) - a controler des l'OUVERTURE du .prmag (voir app._load_
    data_file), pas seulement au moment d'un export MagIC (ou ces
    memes ids etaient jusqu'ici rejetes silencieusement, sans jamais
    alerter l'utilisateur en amont)."""
    if not os.path.exists(ani_path):
        return []
    known = {p.id.strip().upper() for p in pmag_list}
    return sorted({
        t.id.strip() for t in read_all_ani_tensors(ani_path)
        if t.id.strip() and t.id.strip().upper() not in known
    })


# ---------------------------------------------------------------------------
# Detection de doublons (meme etape/memes codes) dans .prmag/.pmagres/
# .pmagint/.pmagani - demande explicite utilisateur ("add a menu for a
# routine to evaluate the files prmag, pmagres, pmagint and pmagani to
# check for measurements with the same step and codes, just a warning to
# help the user clean the files") : QUATRE controles distincts, un par
# fichier, chacun avec sa propre notion de "meme etape et memes codes" -
# jamais correctif (juste une LISTE d'avertissements a afficher, voir
# app.ouvrir_check_duplicates_dialog), le nettoyage reste manuel.
# ---------------------------------------------------------------------------

def find_duplicate_measurements(pmag_list: List["Pmag"]) -> List[str]:
    """.prmag : pour chaque specimen, regroupe ses mesures par (etape,
    cod1, cod2) - signale tout groupe de plus d'une mesure (meme pas de
    desaimantation/acquisition mesure plusieurs fois pour le meme
    specimen, erreur d'import/de saisie frequente)."""
    warnings: List[str] = []
    for p in pmag_list:
        counts: Dict[Tuple[float, str, str], int] = {}
        for m in p.mesures:
            key = (m.etape, m.cod1, m.cod2)
            counts[key] = counts.get(key, 0) + 1
        for (etape, cod1, cod2), n in sorted(counts.items()):
            if n > 1:
                warnings.append(
                    f"{p.id}: step {_fmt1(etape)} code {cod1}{cod2} appears {n} times (.prmag)"
                )
    return warnings


def find_duplicate_results(results: List[FitResult]) -> List[str]:
    """.pmagres : pour les resultats SPECIMEN (jamais une moyenne "mean:",
    qui n'a pas de notion de step1/step2 comparable), regroupe par (id,
    step1, step2, cat1, demag, orig) - signale un meme ajustement
    archive plusieurs fois pour le meme specimen (component/numero de
    fit differents mais MEME segment/type/ancrage : probablement un
    doublon involontaire, pas deux composantes reellement distinctes)."""
    groups: Dict[Tuple[str, int, int, str, str, str], List[str]] = {}
    for r in results:
        if r.id[:5] == "mean:":
            continue
        key = (r.id, r.step_first, r.step_last, r.cat1, r.demag, r.orig)
        groups.setdefault(key, []).append(r.c or r.id)
    warnings = []
    for (rid, step1, step2, cat1, demag, orig), ids in sorted(groups.items()):
        if len(ids) > 1:
            warnings.append(
                f"{rid}: {cat1} fit {step1}-{step2} ({demag}, {orig}) archived {len(ids)} "
                f"times (.pmagres) [{', '.join(ids)}]"
            )
    return warnings


def find_duplicate_pmagint_rows(path: str) -> List[str]:
    """.pmagint : lit TOUTES les lignes brutes (contrairement a
    paleointensity.read_pmagint, qui ne garde QUE la derniere par
    specimen - une revision normale n'est PAS un doublon) - signale
    seulement un (specimen, t1, t2) identique archive plus d'une fois
    (meme intervalle de temperature reinterprete/resauvegarde sans
    aucun changement, pas une vraie revision)."""
    if not os.path.exists(path):
        return []
    from paleointensity import _PMAGINT_HEADER
    counts: Dict[Tuple[str, str, str], int] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if not parts or parts[0] == "specimen" or len(parts) < len(_PMAGINT_HEADER):
                continue
            key = (parts[0], parts[1], parts[2])  # specimen, t1, t2
            counts[key] = counts.get(key, 0) + 1
    warnings = []
    for (specimen, t1, t2), n in sorted(counts.items()):
        if n > 1:
            warnings.append(f"{specimen}: interpretation {t1}-{t2}°C archived {n} times (.pmagint)")
    return warnings


def find_duplicate_ani_tensors(path: str) -> List[str]:
    """.pmagani (section specimen) : regroupe par (id, code2) - un
    specimen ne devrait avoir qu'UN tenseur par type (A0/F0/N0...) ;
    plusieurs signalent un recalcul archive sans nettoyer le precedent."""
    if not os.path.exists(path):
        return []
    counts: Dict[Tuple[str, str], int] = {}
    for t in read_all_ani_tensors(path):
        key = (t.id, t.code2)
        counts[key] = counts.get(key, 0) + 1
    warnings = []
    for (specimen, code2), n in sorted(counts.items()):
        if n > 1:
            warnings.append(f"{specimen}: tensor '{code2}' archived {n} times (.pmagani)")
    return warnings


def apply_inverse_anisotropy(ech: SelectedSample, tensor: AniTensor) -> None:
    """Equivalent de `inverseani` : normalise le tenseur par sa trace/3,
    puis applique l'inverse a CHAQUE mesure de `ech` - MUTE ech.mesures en
    place."""
    rnorm = (tensor.k11 + tensor.k22 + tensor.k33) / 3.0
    if rnorm == 0:
        return
    a, b, c = tensor.k11 / rnorm, tensor.k12 / rnorm, tensor.k13 / rnorm
    d, e, f = tensor.k22 / rnorm, tensor.k23 / rnorm, tensor.k33 / rnorm
    for m in ech.mesures:
        m.x, m.y, m.z = inverse_symmetric_3x3(a, b, c, d, e, f, m.x, m.y, m.z)
