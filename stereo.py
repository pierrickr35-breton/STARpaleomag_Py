"""
Stereonet (canevas de Wulff/Schmidt) : port de `stermes`/`stereoplot(0)`
(plotster.f, Starmac_AWE_v22) via plotlib.PlotContext.

Porte : le cadre du reseau (cercle, graduations N/E/S/W), et le trace des
mesures individuelles (chaque point projete, relie au precedent par un arc
de grand cercle via `gdcerc`, pas une ligne droite - une projection
stereographique/equiaire deforme les trajectoires, un segment droit dans le
plan projete ne represente pas le vrai chemin sur la sphere).

`draw_stereo_net`/`draw_stereo_measurements` sont les briques reutilisees
par `build_stereo_figure` (graphique Stereo autonome, tous les echantillons
selectionnes) ET par zijderveld.py (mini-stereo insere dans le diagramme de
Zijderveld pour l'echantillon affiche, equivalent `stereozijder`).

`draw_stereo_results` porte `ichoixplot==1` (sterres) : directions moyennes
de Fisher ("F"/"m" ou "F"/"i", agregats "mean:") avec cone de confiance
(alpha95), et grands cercles/points pour les ajustements de droite/plan
(cat1 'L'/'f'/'P') - equivalent des deux blocs de `stereoplot` situes apres
le trace des mesures individuelles.
"""

import math
from typing import List, Optional

from matplotlib.figure import Figure

from selection import SelectedSample, angle, apply_orientation, polere
from calcul import FitResult
from plotlib import PlotContext


def superc(la: float, phi: float, dec: float, dip: float, iproj: int) -> tuple:
    """Port de `superc` : projette une direction (dec,dip) autour du pole
    de projection (la,phi) sur le reseau (iproj=1: equiaire/Schmidt,
    sinon stereographique/Wulff). Retourne (u, v, ifl) - ifl=5 si la
    direction tombe sur l'hemisphere "cache" (a dessiner en symbole ouvert),
    11 sinon (symbole plein)."""
    eps = 1e-6
    ifl = 11

    if abs(la) != 90.0:
        # pole de projection oblique (non exerce par les valeurs par defaut
        # de sterparam/projpole, la=-90 - transcription non testee)
        teta = phi - dec
        a = 90.0 - dip
        b = 90.0 - la
        if abs(teta) >= 180.0:
            teta = teta + 360.0 if teta < 0 else teta - 360.0
        tet = abs(teta)
        delta, ang1 = angle(tet, a, b)
        if abs(delta) < 90.0:
            ifl = 5
            delt = 180.0 - delta if delta >= 0 else -180.0 - delta
            del_, ang = angle(ang1, delt, b)
            del_ = abs(del_)
            dip = 90.0 - del_
            dec = phi + ang if teta <= 0 else phi - ang
    else:
        if dip < 0.0:
            dip = -dip
            ifl = 5

    if iproj == 1:
        la2, phi2 = -la, phi + 180.0
    else:
        la2, phi2 = la, phi

    cos0, sin0 = math.cos(math.radians(la2)), math.sin(math.radians(la2))
    sinph = math.sin(math.radians(dec - phi2))
    cosph = math.cos(math.radians(dec - phi2))
    sinla = math.sin(math.radians(dip))
    cosla = math.sqrt(max(0.0, 1.0 - sinla * sinla))
    cosa = sinla * sin0 + cosla * cos0 * cosph
    sina = math.sqrt(max(0.0, 1.0 + eps - cosa * cosa))
    sinb = cosla * sinph / sina if sina else 0.0
    cosb = (sinla * cos0 - cosla * sin0 * cosph) / sina if sina else 0.0

    if iproj != 1:
        r = (1.0 + cosa) / sina
        return r * cosb, r * sinb, ifl

    if cosa + 1.0 <= 0.0:
        return 1.0, 1.0, ifl
    r = sina / math.sqrt(0.5 + 0.5 * cosa)
    u = r * cosb / math.sqrt(2.0)
    v = -r * sinb / math.sqrt(2.0)
    return u, v, ifl


def gdcerc(
    ctx: PlotContext, la: float, phi: float, iproj: int,
    dec1: float, dip1: float, dec2: float, dip2: float, r: float,
) -> None:
    """Port de `gdcerc` : trace l'arc de grand cercle entre deux directions
    projetees, par petits pas angulaires (~3 deg)."""
    decc = dec1
    teta = dec1 - dec2
    a = 90.0 - dip2
    if a == 0.0:
        a = 0.1
    if a == 180.0:
        a = 179.9
    delti = 0.0
    icont = 1
    b = 90.0 - dip1
    if b == 180.0:
        b = 179.9
    if b == 0.0:
        b = 0.1
    if abs(teta) >= 180.0:
        teta = teta + 360.0 if teta < 0 else teta - 360.0
    tet = abs(teta)
    pas = 3.0
    delta, ang = angle(tet, a, b)
    a = 0.0
    ifll = 0

    while True:
        u, v, ifl = superc(la, phi, dec1, dip1, iproj)
        if ifl == 11 and ifll == 11:
            icont = 0
        if ifl == 11 and ifll == 5:
            icont = 1
        if ifl == 5 and ifll == 11:
            icont = 1
        ctx.plot(v * r, u * r, 2 + (icont % 2))
        if a == -1.0:
            return
        icont += 1
        if abs(delti - delta) < 3.0:
            dec1, dip1, a = dec2, dip2, -1.0
            continue
        delti += pas
        ai, pp = angle(ang, delti, b)
        dec1 = decc - pp if teta > 0 else decc + pp
        if dec1 > 360.0:
            dec1 -= 360.0
        dip1 = 90.0 - ai
        ifll = ifl


_ORIENT_LABELS = {
    1: ("Sample Coord.", 3),
    2: ("In situ", 4),
    3: ("Tilt Corrected", 5),
}


def draw_stereo_net(
    ctx: PlotContext,
    orientation: int = 1,
    la: float = -90.0,
    phi: float = 0.0,
    iproj: int = 1,
    dimster: float = 18.0,
    show_orient_label: bool = True,
) -> float:
    """Dessine le cadre du reseau (cercle, graduations, poles, N/E/S/W et
    etiquette d'orientation) a l'origine COURANTE de `ctx`. Retourne le
    rayon `r` du reseau (cm), a reutiliser pour placer les mesures."""
    r = dimster / 3.0

    ctx.thickn(1.0)
    ctx.newpen(1)
    # cadre du reseau (ifill=8 dans le Fortran : rempli de "snow" a l'ecran,
    # mais ifilpoly reste a 0 pour l'export SVG, qui ecrit donc fill="none" -
    # verifie sur un export SVG reel, on reproduit ce comportement)
    ctx.circl2(0.0, 0.0, r, 1, 0)

    # graduations (croix, symbole type 5) tous les 10 deg de pendage, sur les 4 azimuts cardinaux
    ctx.thickn(0.25)
    for teta in (0.0, 90.0, 180.0, 270.0):
        for k in range(1, 18):
            dipp = 10.0 * k
            u, v, ifl = superc(la, phi, teta, 90.0 - dipp, iproj)
            if ifl == 5 or (u * u + v * v) > 0.98:
                continue
            ctx.symbol(v * r, u * r, dimster / 72.0, 5, -1)

    u, v, ifl = superc(la, phi, 0.0, 90.0, iproj)
    if ifl != 5:
        ctx.symbol(v * r, u * r, dimster / 35.0, 5, -1)
    u, v, ifl = superc(la, phi, 0.0, -90.0, iproj)
    if ifl != 5:
        ctx.symbol(v * r, u * r, dimster / 35.0, 13, -1)

    if show_orient_label:
        text, pen = _ORIENT_LABELS.get(orientation, ("", 1))
        ctx.newpen(pen)
        ctx.plottxt(-r - 5 * dimster / 80, r + dimster / 60, dimster / 40, text)
        ctx.newpen(1)

    # petits traits vers l'interieur a N/E/S/W (memes coordonnees brutes
    # que StereoUtils_Py, oublies dans une premiere passe - pas de swap
    # `plott`, ces `plot()` du source ne passent pas par ce wrapper).
    ctx.plot(r, 0.0, 3)
    ctx.plot(r - dimster / 60.0, 0.0, 2)
    ctx.plottxt(r + dimster / 190.9, -dimster / 65.6233, dimster / 38.158692, "E")

    ctx.plot(0.0, r, 3)
    ctx.plot(0.0, r - dimster / 60.0, 2)
    ctx.plottxt(-dimster / 68.76, (r + dimster / 21.6) - dimster / 27.63, dimster / 38.158692, "N")

    ctx.plot(-r + dimster / 60.0, 0.0, 3)
    ctx.plot(-r, 0.0, 2)
    ctx.plottxt(-r - 4 * dimster / 84, -dimster / 65.6233, dimster / 38.158692, "W")

    ctx.plot(0.0, -r, 3)
    ctx.plot(0.0, -r + dimster / 60.0, 2)
    ctx.plottxt(-dimster / 68.76, (-r - dimster / 1050.0) - dimster / 27.63, dimster / 38.158692, "S")

    return r


def draw_stereo_measurements(
    ctx: PlotContext,
    selected: List[SelectedSample],
    r: float,
    orientation: int = 1,
    la: float = -90.0,
    phi: float = 0.0,
    iproj: int = 1,
    point_size: float = 0.324,
) -> None:
    """Dessine les mesures des echantillons de `selected`, reliees par des
    arcs de grand cercle, a l'origine COURANTE de `ctx` (voir draw_stereo_net
    pour `r`)."""
    ctx.thickn(0.5)
    ctx.newpen(1)
    for ech in selected:
        prev_dec = prev_dip = None
        for m in ech.mesures:
            xx, yy, zz = apply_orientation(m.x, m.y, m.z, ech, orientation)
            dec = math.degrees(math.atan2(yy, xx))
            rh = math.hypot(xx, yy)
            dip = math.degrees(math.atan2(zz, rh))

            u, v, ifl = superc(la, phi, dec, dip, iproj)
            if prev_dec is not None:
                gdcerc(ctx, la, phi, iproj, prev_dec, prev_dip, dec, dip, r)
            symtype = 8 if ifl == 5 else 14
            ctx.symbol(v * r, u * r, point_size, symtype, -1)
            prev_dec, prev_dip = dec, dip


def build_stereo_figure(
    selected: List[SelectedSample],
    orientation: int = 1,
    la: float = -90.0,
    phi: float = 0.0,
    iproj: int = 1,
    fig: Optional[Figure] = None,
) -> Figure:
    """Equivalent de `stermes` (stereoplot(0)) : cadre du reseau + mesures
    individuelles de la selection, reliees par des arcs de grand cercle."""
    dimster = 12.0 * 1.5  # diam() : diametre par defaut (cm) * 1.5
    point_size = (0.18 * dimster) / 10.0  # sizepoint() : taille par defaut

    if fig is None:
        fig = Figure(figsize=(5.5, 5.5), dpi=100)
    else:
        fig.clear()
    ax = fig.add_subplot(111)
    ctx = PlotContext(ax)
    ctx.clear()
    ctx.plot(0.0, 0.0, -3)  # origine au centre (pas de mise en page multi-panneaux ici)

    r = draw_stereo_net(ctx, orientation, la, phi, iproj, dimster)
    draw_stereo_measurements(ctx, selected, r, orientation, la, phi, iproj, point_size)

    ax.relim()
    ax.autoscale_view()
    ax.set_title(", ".join(s.id for s in selected) if len(selected) <= 3 else f"{len(selected)} samples")
    fig.tight_layout()
    return fig


def circle(al: float, ai: float, ad: float, ph: float) -> tuple:
    """Port de `circle(al,ai,ad,ph,ei,ed)` (plotster.f:652) : point
    (ei=inclinaison, ed=declinaison) sur le petit cercle de rayon angulaire
    `al` autour de la direction (ad=declinaison, ai=inclinaison), parametre
    par l'angle `ph` (0-360) - utilise pour tracer les cones de confiance
    alpha95 des directions moyennes de Fisher."""
    f = math.pi / 180.0
    sal, cal = math.sin(al * f), math.cos(al * f)
    si, ci = math.sin(ai * f), math.cos(ai * f)
    sd, cd = math.sin(ad * f), math.cos(ad * f)
    zpx = cal * cd * ci
    zpy = cal * sd * ci
    zpz = cal * si
    sp, cp = math.sin(ph * f), math.cos(ph * f)
    xp = sal * cp
    yp = sal * sp
    x = xp * (sd * sd + cd * cd * si) + yp * sd * cd * (si - 1.0) + zpx
    y = xp * (cd * sd * si - cd * sd) + yp * (cd * cd + sd * sd * si) + zpy
    z = xp * (-cd * ci) + yp * (-sd * ci) + zpz
    s = math.sqrt(x * x + y * y)
    ed = math.degrees(math.atan2(y / s, x / s))
    if ed < 0.0:
        ed += 360.0
    ei = math.degrees(math.atan2(z, s))
    return ei, ed


def _correct_dec_inc(res: FitResult, orientation: int) -> tuple:
    """Equivalent de `correct(rd,ri,cdip,caz,dip,str)` (StarUtil.f:66) :
    reapplique la correction d'orientation courante a la direction BRUTE
    stockee dans un resultat (res.dec/res.inc, repere echantillon d'origine),
    avec les cin/caz/dip/str_ PROPRES a ce resultat (pas d'un SelectedSample
    live - FitResult porte deja ces memes champs, `apply_orientation` n'a
    besoin que de leur presence, pas du type exact de l'objet).

    Bug corrige - meme correctif et memes raisons que calcul.py:
    _correct_dec_inc (duplique ici pour eviter un import circulaire) :
    un "mean:" est deja fige dans une orientation precise (par3_mean),
    pas dans le repere echantillon brut - le renvoyer tel quel plutot que
    de lui appliquer apply_orientation (qui, avec cin=caz=0 jamais
    renseignes sur une moyenne, introduirait une rotation parasite)."""
    if res.id[:5] == "mean:":
        return res.dec, res.inc
    incr, decr = math.radians(res.inc), math.radians(res.dec)
    x = math.cos(incr) * math.cos(decr)
    y = math.cos(incr) * math.sin(decr)
    z = math.sin(incr)
    xx, yy, zz = apply_orientation(x, y, z, res, orientation)
    _, dec, inc = polere(xx, yy, zz)
    return dec, inc


def draw_stereo_results(
    ctx: PlotContext,
    results: List[FitResult],
    r: float,
    orientation: int = 1,
    la: float = -90.0,
    phi: float = 0.0,
    iproj: int = 1,
    point_size: float = 0.324,
    nbech: Optional[int] = None,
) -> None:
    """Port de la partie `ichoixplot==1` de `stereoplot` (sterres,
    plotster.f:253-403) : a l'origine COURANTE de `ctx` (voir
    draw_stereo_net pour `r`) - directions moyennes de Fisher (cat1='F',
    cat2 'm' ou 'i') avec cone de confiance alpha95, points de direction
    pour les ajustements de droite/direction unique (cat1 'L'/'f'), et
    grands cercles pour les ajustements de plan (cat1='P'). `nbech` :
    nombre d'echantillons de la selection courante - le libelle "fit
    between X et Y" ne s'affiche que si nbech==1 ET un seul resultat."""
    ctx.thickn(0.5)

    # --- directions moyennes de Fisher + cone de confiance (lignes 253-306) :
    # cat2=='m' utilise dec/inc TELS QUELS (pas de correction d'orientation -
    # deja dans le bon repere, resultat "multi-echantillons") ; cat2=='i'
    # reapplique `correct()` (resultat "individuel", repere echantillon).
    means = []
    for res in results:
        if res.cat1 == "F" and res.cat2 == "m":
            means.append((res.dec, res.inc, res.mad))
        elif res.cat1 == "F" and res.cat2 == "i":
            dec, inc = _correct_dec_inc(res, orientation)
            means.append((dec, inc, res.mad))

    for rdec, rinc, alph in means:
        ctx.newpen(3)
        u, v, ifl = superc(la, phi, rdec, rinc, iproj)
        ctx.symbol(v * r, u * r, point_size * 2.0 / 3.0, 20 if ifl == 5 else 14, -1)

        ph = 0.0
        ei, ed = circle(alph, rinc, rdec, ph)
        u, v, _ifl = superc(la, phi, ed, ei, iproj)
        ctx.plot(v * r, u * r, 3)
        ctx.newpen(5)
        for _ in range(120):
            ph += 3.0
            ei, ed = circle(alph, rinc, rdec, ph)
            u, v, _ifl = superc(la, phi, ed, ei, iproj)
            ctx.plot(v * r, u * r, 2)
    if means:
        ctx.newpen(1)

    # --- points de direction pour L/f (lignes 322-335) ---
    line_results = [res for res in results if res.cat1 in ("L", "f")]
    line_dirs = [_correct_dec_inc(res, orientation) for res in line_results]
    if line_dirs:
        ctx.newpen(3)
        for (dec, inc), res in zip(line_dirs, line_results):
            u, v, ifl = superc(la, phi, dec, inc, iproj)
            ctx.symbol(v * r, u * r, point_size, 8 if ifl == 5 else 14, -1)
            # Point interactif invisible (clic-pour-info) - demande
            # explicite utilisateur ("cliquer sur des donnees d'un
            # graphique... sur un stereo results, le nom du specimen").
            ctx.pick_point(v * r, u * r, "stereo_specimen", res.id)
        ctx.newpen(1)

    # --- grands cercles pour les plans P (lignes 337-390) : les newpen()
    # bases sur numcomp dans la boucle de collecte du Fortran sont sans
    # effet (ecrases par le newpen(4) fixe juste avant la boucle de trace),
    # non reproduits ici - tous les grands cercles sortent en vert (pen 4).
    plane_poles = [
        (res, *_correct_dec_inc(res, orientation)) for res in results if res.cat1 == "P"
    ]

    if plane_poles:
        ctx.newpen(4)
        for _res, pdec, pinc in plane_poles:
            if pinc < 0.0:
                dipp = pinc + 90.0
                decc = pdec
            else:
                dipp = -pinc + 90.0
                decc = pdec + 180.0

            u, v, _ifl = superc(la, phi, decc, dipp, iproj)
            ctx.plot(v * r, u * r, 3)

            dec_cur, dip_cur = decc, dipp
            for i in range(1, 5):
                dec2 = dec_cur + 90.0
                dip2 = 0.0
                if i == 2:
                    dip2 = -dipp
                if i == 4:
                    dip2 = dipp
                gdcerc(ctx, la, phi, iproj, dec_cur, dip_cur, dec2, dip2, r)
                dec_cur, dip_cur = dec2, dip2
        ctx.newpen(1)

        if nbech == 1 and len(results) == 1:
            fit = results[0]
            titre = f" ech:{fit.id:<12s}  fit between {fit.step_first:4.0f} et {fit.step_last:4.0f}"
            dimster_equiv = r * 3.0
            u = -dimster_equiv / 3.0 - 1.0
            v = -dimster_equiv / 3.0 - 1.2
            ctx.plottxt(u, v, dimster_equiv / 60.0, titre)


def build_stereo_results_figure(
    results: List[FitResult],
    orientation: int = 1,
    la: float = -90.0,
    phi: float = 0.0,
    iproj: int = 1,
    nbech: Optional[int] = None,
    fig: Optional[Figure] = None,
) -> Figure:
    """Equivalent de `sterres` (stereoplot(1)) : cadre du reseau + directions
    moyennes de Fisher (avec cone de confiance) + grands cercles/points pour
    les ajustements de droite/plan de `results` (self.results cote app.py)."""
    dimster = 12.0 * 1.5
    point_size = (0.18 * dimster) / 10.0

    if fig is None:
        fig = Figure(figsize=(5.5, 5.5), dpi=100)
    else:
        fig.clear()
    ax = fig.add_subplot(111)
    ctx = PlotContext(ax)
    ctx.clear()
    ctx.plot(0.0, 0.0, -3)

    r = draw_stereo_net(ctx, orientation, la, phi, iproj, dimster)
    draw_stereo_results(ctx, results, r, orientation, la, phi, iproj, point_size, nbech=nbech)

    ax.relim()
    ax.autoscale_view()
    ax.set_title("Stereo Results")
    fig.tight_layout()
    return fig
