"""
Diagramme de Zijderveld : port fidele de `zijderplot` (plotorthog.f,
Starmac_AWE_v22), y compris les graduations d'echelle automatique (log10),
les numeros d'etape a cote de chaque point, et le texte recapitulatif
(echelle, id echantillon, type de desaimantation, code d'orientation).
Fonctionne avec n'importe quel objet "ctx" exposant l'API PlotContext
(plot/symbol/circl2/plottxt/number/newpen/thickn) - PlotContext (matplotlib)
ou SVGWriter (port fidele de svginit/svgplot).

Convention d'axes (logique originale de l'auteur, PAS une bibliotheque
Absoft generique) :
- Les deux projections partagent le meme axe HORIZONTAL, qui porte N (par
  defaut) ou E (bascule automatique, equivalent iorzij==0/itest1 : bascule
  quand E a plus d'amplitude que N).
- Down/Up est TOUJOURS porte par l'axe VERTICAL pour la projection
  verticale ; l'autre composante (E ou N) est portee par l'axe vertical
  pour la projection horizontale.
- Signe : screen-x = -N (ou -E si bascule), N pointe donc vers la GAUCHE
  quand N est l'axe principal - comportement d'origine du Fortran.
- Symboles fideles : type 14 (cercle plein, horizontale), type 20 (cercle
  blanc, verticale), type 5 (plus, graduations).
"""

import math
from typing import List, Optional

from matplotlib.figure import Figure

from selection import SelectedSample, apply_orientation, split_experiments, experiment_kind
from stereo import draw_stereo_net, draw_stereo_measurements, draw_stereo_results
from plotlib import PlotContext, char_width_cm

_DEMAG_TYPES = {"F": "AF", "D": "TH", "S": "TH", "C": "CH"}
_ORIENT_CODES = {1: "SC", 2: "IS", 3: "BC"}


def _scale_factor(ech: SelectedSample) -> float:
    """1e3 (masse, Am2/kg) ou 1e6 (volume, A/m) - SAUF si `ech.vol` est
    manquant (0/None), auquel cas 1.0 (aucun facteur applique) - demande
    explicite utilisateur ("when the volume or the mass of the sample is
    not given, best to show the data in total moment as done in the
    list... the previous Fortran was systematically dividing by mass and
    volume and was not expecting no vol and no mass"). Le Fortran
    d'origine divisait TOUJOURS par `ech.vol` sans jamais verifier qu'il
    etait renseigne - BUG CONFIRME (pas seulement une division par zero
    evitee, une VALEUR FAUSSE affichee : `ech.vol or 1.0` seul, sans ce
    garde-fou sur `factor`, appliquerait quand meme le facteur 1e3/1e6 a
    un moment BRUT, affichant un nombre sans rapport avec l'unite
    revendiquee "A/m"/"Am2/kg"). Sans volume/masse, le moment total brut
    (Am2) est affiche tel quel, comme le fait deja list_measurements
    (colonne "Mtot Am2", toujours calculee independamment du volume) -
    voir aussi le libelle d'unite plus bas (`unit`)."""
    if not ech.vol:
        return 1.0
    return 1.0e3 if ech.norme == "m" else 1.0e6


def _oriented_ned(x: float, y: float, z: float, ech: SelectedSample, orientation: int) -> tuple:
    """(x,y,z) brut -> (N,E,Down) corrige et mis a l'echelle (comme xzi/yzi/zzi)."""
    factor = _scale_factor(ech)
    vol = ech.vol or 1.0
    xx, yy, zz = apply_orientation(x, y, z, ech, orientation)
    return xx * factor / vol, yy * factor / vol, zz * factor / vol


def _compute_origin_uv(dimzij: float, amax: float, xn: float, ym: float, yn: float,
                        zm: float, zn: float) -> tuple:
    """Port du calcul de l'origine du NEV (zijder2, lignes ~2028-2046) :
    positionne le diagramme sur la page en fonction de l'asymetrie reelle
    des donnees (xn/ym/yn/zm/zn = extremites BRUTES xzi/yzi/zzi, PAS les
    coordonnees ecran deja negees/mises a l'echelle). Verifie exactement
    (a l'arrondi pres) contre l'origine reelle (387.638, 550.977 px) de
    zijder-11CL7801A.svg."""
    if zm < 0.0 and zn < 0.0 and ym > 0.0 and yn > 0.0:
        t = 0.0
    elif zm > -yn:
        t = zm / amax
    else:
        t = abs(yn) / amax
    tt = 2.0 * dimzij / 3.0
    v = t * tt + dimzij / 6.0
    u = (dimzij / 6.0 - tt * xn / amax) if xn < 0.0 else dimzij / 6.0
    return u, v


def _nice_scale(amax: float) -> tuple:
    """Port du calcul `ik`/`scal` (ALOG10(amax) puis arrondi a la puissance
    de 10 inferieure) : pas de graduation "rond" pour l'echelle des axes."""
    if amax <= 0:
        return 0, 1.0
    xx = math.log10(amax)
    ik = math.trunc(xx) - 1 if xx < 0 else math.trunc(xx)
    return ik, 10.0 ** ik


def _draw_ticks(ctx, cm_scale: float, scal: float, lo: float, hi: float,
                 tick_size: float, along: str, label_lo: str, label_hi: str,
                 label_size: float, label_offset: float) -> None:
    """Port des boucles de graduation de `zijderplot` (lignes ~1168-1218) :
    des croix (symbole 5) tous les `scal`, de 0 vers `lo` puis de 0 vers
    `hi`, reliees par un trait jusqu'a la valeur exacte en bout de course,
    etiquette a chaque extremite."""

    def pt(v: float) -> tuple:
        return (cm_scale * v, 0.0) if along == "x" else (0.0, cm_scale * v)

    if lo < 0:
        ctx.plot(*pt(0.0), 3)
        xt = 0.0
        while True:
            xt -= scal
            if xt < lo:
                break
            ctx.symbol(*pt(xt), tick_size, 5, -2)
            ctx.plot(*pt(xt), 3)
        ctx.plot(*pt(lo), 2)
        lx, ly = pt(lo)
        if along == "x":
            ctx.plottxt(lx - label_offset, -label_offset * 0.6, label_size, label_lo)
        else:
            ctx.plottxt(-label_offset * 0.6, ly - label_offset, label_size, label_lo)

    if hi > 0:
        ctx.plot(*pt(0.0), 3)
        xt = 0.0
        while True:
            xt += scal
            if xt > hi:
                break
            ctx.symbol(*pt(xt), tick_size, 5, -2)
            ctx.plot(*pt(xt), 3)
        ctx.plot(*pt(hi), 2)
        hx, hy = pt(hi)
        if along == "x":
            ctx.plottxt(hx + label_offset * 0.3, -label_offset * 0.6, label_size, label_hi)
        else:
            ctx.plottxt(-label_offset * 0.6, hy + label_offset * 0.3, label_size, label_hi)


def draw_zijderveld(
    ctx,
    ech: SelectedSample,
    orientation: int = 1,
    fits: Optional[List] = None,
    auto_axis: bool = True,
    show_stereo: bool = True,
    dimzij: float = 15.0,
) -> None:
    """Dessine le diagramme de Zijderveld pour l'echantillon `ech` sur `ctx`
    (deja positionne par l'appelant - un `ctx.plot(0,0,-3)` prealable
    etablit l'origine). `fits` : ajustements de droite (calcul.FitResult)
    pour cet echantillon, superposes en rouge si fournis."""
    # un echantillon peut enchainer plusieurs protocoles dans sa liste de
    # mesures - deux cas reels distincts, demandant un traitement DIFFERENT :
    # (1) AF demag, puis acquisition IRM, puis desaimantation thermique de
    #     cette IRM (ex. 19DN1607B) - les runs apres l'AF ne sont PAS de la
    #     NRM, a exclure du Zijderveld (qui ne represente que la
    #     desaimantation de la NRM).
    # (2) AF demag de la NRM insuffisante, poursuivie par une desaimantation
    #     thermique de CETTE MEME NRM (demande explicite utilisateur : "add
    #     a thermal demag when AF is not efficient... keep AF and thermal
    #     data for the same zijderveld") - ce run thermique EST de la NRM,
    #     a garder avec le run AF.
    # Distinction : experiment_kind - on exclut UNIQUEMENT 'I' (IRM
    # acquise) et 'D_IRM' (thermique d'une IRM acquise, PAS la NRM) ; tout
    # le reste (F/D mais aussi le cod1 brut d'autres protocoles NRM
    # legitimes - R/V/P/S/N/X/Y/Z/A pour la paleointensite/ATRM, affiches
    # eux aussi sur un Zijderveld, ex. dans build_paleoint_review_figure)
    # est concatene depuis le debut. Bug reel corrige : une premiere
    # version de ce filtre necessitait litteralement 'F' ou 'D', ce qui
    # vidait le Zijderveld de TOUT echantillon de paleointensite (codes
    # N/R/V/P, jamais 'F'/'D').
    runs = split_experiments(ech.mesures)
    mesures_zij = []
    for run in runs:
        if experiment_kind(run) in ("I", "D_IRM"):
            break
        mesures_zij.extend(run)

    n_vals, e_vals, d_vals, steps = [], [], [], []
    demag_code = ""
    demag_types_seen: List[str] = []  # ex. ["AF","TH"] si AF puis thermique de la meme NRM
    for m in mesures_zij:
        n, e, d = _oriented_ned(m.x, m.y, m.z, ech, orientation)
        n_vals.append(n); e_vals.append(e); d_vals.append(d)
        steps.append(m.etape)
        demag_code = m.cod1
        dt = _DEMAG_TYPES.get(m.cod1)
        if dt and dt not in demag_types_seen:
            demag_types_seen.append(dt)

    # Bascule automatique N<->E (equivalent iorzij==0 : itest1)
    xm, xn = max(n_vals + [0.0]), min(n_vals + [0.0])
    ym, yn = max(e_vals + [0.0]), min(e_vals + [0.0])
    amax = max(abs(xm), abs(xn), abs(xm - xn))
    swapped = auto_axis and (abs(ym) > amax or abs(yn) > amax or abs(ym - yn) > amax)

    def primsec(n: float, e: float) -> tuple:
        return (-e, n) if swapped else (n, e)

    primary_vals, secondary_vals = (list(v) for v in zip(*(primsec(n, e) for n, e in zip(n_vals, e_vals))))
    screen_x = [-p for p in primary_vals]
    vert_y = [-d for d in d_vals]

    # amax final (apres bascule), y compris le test de chevauchement N-E/Up
    # entre les deux projections (lignes 1113-1125 du Fortran). zm2/zn2 et
    # ym2/yn2 utilisent les valeurs BRUTES (zzi=d_vals, yzi=secondary_vals),
    # pas vert_y (=-d_vals) : c'est l'etendue des DONNEES, independante du
    # sens d'affichage a l'ecran.
    xm2, xn2 = max(screen_x + [0.0]), min(screen_x + [0.0])
    ym2, yn2 = max(secondary_vals + [0.0]), min(secondary_vals + [0.0])
    zm2, zn2 = max(d_vals + [0.0]), min(d_vals + [0.0])
    xn_raw = min(primary_vals + [0.0])  # xzi brut (non negue) - pour _compute_origin_uv
    amax = max(
        abs(xm2), abs(xn2), abs(ym2), abs(yn2), abs(zm2), abs(zn2), abs(xm2 - xn2),
        abs(zm2 + ym2), abs(zn2 + yn2), abs(zm2 - zn2), abs(ym2 - yn2),
    ) or 1.0

    secondary_label = "N" if swapped else "E"

    # cf zijderplot : le cote xm (screen-x NEGATIF, cote -primary positif)
    # porte in1 ('N' ou 'W' apres bascule) ; le cote xn (screen-x POSITIF)
    # porte in2 ('S' ou 'E' apres bascule) - verifie empiriquement contre
    # zijder-11CL7801A.svg (export original : N a gauche = screen-x negatif).
    if swapped:
        label_neg, label_pos = "W", "E"
    else:
        label_neg, label_pos = "N", "S"

    ik, scal = _nice_scale(amax)
    # dimzij (parametre, defaut 15.0) : constante de mise en page - 15.0
    # verifiee par calcul exact de l'origine du NEV (cf _compute_origin_uv)
    # contre zijder-11CL7801A.svg : (387.638,550.977) obtenu au pixel pres -
    # c'est la valeur de la vue Zijderveld PLEINE PAGE (`zijderplot`,
    # plotorthog.f, dimzij fixe a 15 dans ce contexte). La vue combinee
    # d'apercu paleointensite (visi_Paleoint.f, panneau Zijderveld en
    # quart de page) utilise elle-meme dimzij=9.0 - passe explicitement
    # par build_paleoint_review_figure, ne PAS changer ce defaut de 15.0
    # ici (casserait la fidelite verifiee de la vue pleine page). tt=2*
    # dimzij/3 est la demi-largeur physique en cm du diagramme (fixe,
    # independante de amax : c'est amax qui est mis a l'echelle pour
    # occuper cette largeur).
    cm_scale = (2.0 * dimzij / 3.0) / amax
    tick_size = 0.03 * dimzij
    label_size = dimzij / 25.0
    label_offset = dimzij / 20.0

    # Origine du diagramme sur la page : offset initial (zijder2, avant la
    # boucle d'echantillon) puis origine du NEV (adaptee a l'asymetrie des
    # donnees) - l'appelant doit avoir positionne `ctx` sur l'origine de PAGE
    # (equivalent scrhor=90,scrvor=600) avant d'appeler draw_zijderveld.
    ctx.plot(-2.0, -5.5, -3)
    u, v = _compute_origin_uv(dimzij, amax, xn_raw, ym2, yn2, zm2, zn2)
    ctx.plot(-(u - dimzij), v, -3)

    ctx.newpen(1)
    ctx.thickn(1.0)

    # graduations horizontales (axe partage N/E) : xn2 (gauche) -> label_neg, xm2 (droite) -> label_pos
    _draw_ticks(ctx, cm_scale, scal, xn2, xm2, tick_size, "x",
                label_neg, label_pos, label_size, label_offset)
    # graduations verticales : zn2/zm2 = Down (bas), yn2/ym2 = E-ou-N (haut)
    xz = max(zm2, -yn2)
    if xz != 0.0:
        _draw_ticks(ctx, cm_scale, scal, -xz, 0.0, tick_size, "y",
                    "Down", "", label_size, label_offset)
    xy = max(ym2, -zn2)
    if xy != 0.0:
        _draw_ticks(ctx, cm_scale, scal, 0.0, xy, tick_size, "y",
                    "", secondary_label, label_size, label_offset)

    marker_size = dimzij * 0.02
    ctx.thickn(0.5)

    # projection horizontale (symbole plein 14) : 1er point noir (pen 1),
    # points suivants + traits de liaison en ROUGE (pen 3) - port de
    # `zijder2` (verifie contre zijder-11CL7801A.svg : pas la routine
    # `zijderplot`, monochrome, qui n'est pas celle utilisee a l'export).
    if screen_x:
        ctx.newpen(1)
        ctx.symbol(cm_scale * screen_x[0], cm_scale * secondary_vals[0], marker_size, 14, -1)
        ctx.newpen(3)
        for x_cm, y_cm in zip(screen_x[1:], secondary_vals[1:]):
            ctx.symbol(cm_scale * x_cm, cm_scale * y_cm, marker_size, 14, -2)

    # projection verticale (Up) : traits de liaison NOIRS (pen 1), puis
    # symboles (cercle blanc, 20) en VERT (pen 4) - meme source.
    if screen_x:
        ctx.newpen(1)
        ctx.plot(cm_scale * screen_x[0], cm_scale * vert_y[0], 3)
        for x_cm, y_cm in zip(screen_x[1:], vert_y[1:]):
            ctx.plot(cm_scale * x_cm, cm_scale * y_cm, 2)
        ctx.newpen(4)
        for x_cm, y_cm in zip(screen_x, vert_y):
            ctx.symbol(cm_scale * x_cm, cm_scale * y_cm, marker_size, 20, -1)
        ctx.newpen(1)

    # Points interactifs invisibles (clic-pour-info) sur CHAQUE point trace,
    # dans les DEUX projections - demande explicite utilisateur ("cliquer
    # sur des donnees d'un graphique... par exemple sur un zijderveld, la
    # temperature"). `data` = la Measurement elle-meme (etape/cod1/cod2),
    # laisse la mise en forme au gestionnaire de clic cote app.py plutot
    # que de figer un texte ici.
    for x_cm, y_sec, y_vert, m in zip(screen_x, secondary_vals, vert_y, mesures_zij):
        ctx.pick_point(cm_scale * x_cm, cm_scale * y_sec, "zijderveld_step", m)
        ctx.pick_point(cm_scale * x_cm, cm_scale * y_vert, "zijderveld_step", m)

    # ajustements de droite : couleur par numero de composante (numcomp),
    # epaisseur doublee - port de `zijder2` (lignes ~2276-2300).
    _FIT_PENS = {1: 3, 2: 4, 3: 5, 4: 8}
    matching_fits = [f for f in (fits or []) if f.id == ech.id]
    if matching_fits:
        ctx.thickn(2.0)
        for fit in matching_fits:
            ctx.newpen(_FIT_PENS.get(fit.numcomp, 3))
            pts_n, pts_e, pts_d = [], [], []
            for tx, ty, tz in zip(fit.tx, fit.ty, fit.tz):
                n, e, d = _oriented_ned(tx, ty, tz, ech, orientation)
                pts_n.append(n); pts_e.append(e); pts_d.append(d)
            pri, sec = zip(*(primsec(n, e) for n, e in zip(pts_n, pts_e)))
            fx = [-p for p in pri]
            ctx.plot(cm_scale * fx[0], cm_scale * sec[0], 3)
            ctx.plot(cm_scale * fx[1], cm_scale * sec[1], 2)
            ctx.plot(cm_scale * fx[0], cm_scale * -pts_d[0], 3)
            ctx.plot(cm_scale * fx[1], cm_scale * -pts_d[1], 2)
        if len(matching_fits) == 1:
            fit = matching_fits[0]
            ctx.newpen(1)
            # meme niveau vertical que le texte "scale: ..." (b + dimzij/18.0
            # plus bas, ligne ~298) - demande explicite de l'utilisateur.
            ctx.plottxt(
                cm_scale * amax * 0.1, -v + dimzij / 18.0, dimzij / 40.0,
                f"fit between {fit.step_first} and {fit.step_last}",
            )
        ctx.thickn(0.5)
        ctx.newpen(1)

    # numeros d'etape a cote de chaque point (lignes 1246-1254) : attaches a
    # la trace secondaire (E/N) si son etendue domine, sinon a la trace
    # verticale (Down/Up) - seul le numero est affiche, pas le code demag.
    e_extent = max(ym2, -yn2)
    d_extent = max(zm2, -zn2)
    use_secondary = e_extent >= d_extent
    ctx.newpen(1)
    # un point sur 3, police 3x plus grande, decalee vers la droite de la
    # largeur d'un caractere - demande explicite de l'utilisateur (les
    # numeros serres/petits etaient difficiles a lire).
    step_label_height = 0.35
    step_label_dx = char_width_cm(step_label_height)
    for i, etape in enumerate(steps):
        if i % 3 != 0:
            continue
        y_cm = secondary_vals[i] if use_secondary else vert_y[i]
        ctx.number(
            cm_scale * screen_x[i] - 0.01 * dimzij + step_label_dx,
            cm_scale * y_cm - dimzij / 140.0,
            step_label_height, float(etape), 0.0, -1,
        )

    # texte recapitulatif : echelle, id, type de desaimantation, orientation
    # "Am2" (moment total brut, PAS normalise) si vol/masse absent - voir
    # _scale_factor.
    unit = "Am2" if not ech.vol else ("Am2/kg" if ech.norme == "m" else "A/m")
    scale_val = 10.0 ** ik
    scale_text = f"scale: {scale_val:g}{unit}"
    demag_text = "+".join(demag_types_seen) if demag_types_seen else _DEMAG_TYPES.get(demag_code, "")
    orient_text = _ORIENT_CODES.get(orientation, "")

    # positions du bloc de texte : "a"/"b" derives de u,v (zijder2, lignes
    # ~2224-2233), PAS une formule ad-hoc - verifie exactement (au pixel
    # pres) contre zijder-11CL7801A.svg pour "scale:", id, type et orientation.
    b = -v
    a = -u + 4.0 * dimzij / 5.0
    ctx.plottxt(-a + dimzij / 5.0, b + dimzij / 18.0, dimzij / 35.0, scale_text)
    b -= dimzij / 26.0
    a -= dimzij / 25.0
    ctx.plottxt(-a, b, dimzij / 25.0, ech.id, nchar=12)
    a -= 12.0 * dimzij / 25.0
    ctx.plottxt(-a + dimzij / 25.0, b, dimzij / 25.0, demag_text, nchar=2)
    ctx.plottxt(-a + 4.0 * dimzij / 25.0, b, dimzij / 25.0, orient_text, nchar=2)

    if show_stereo:
        # stereo (equivalent stereozijder, plotster.f) : position ABSOLUE
        # fixe sur la page (6.0, 14.0 cm depuis l'origine de PAGE), PAS
        # relative a l'origine du zijderveld - `call plot(0.,0.,-4)` revient
        # d'abord a l'origine de page, comme l'original. dimster=0.7*dimstereo
        # ; dimstereo=21.0 calibre contre zijder-11CL7801A.svg (rayon reel
        # 138.90px -> dimster=14.7cm -> dimstereo=21.0).
        dimstereo = 21.0
        dimster = 0.7 * dimstereo
        ctx.plot(0.0, 0.0, -4)
        ctx.plot(6.0, 14.0, -3)
        r_mini = draw_stereo_net(ctx, orientation, dimster=dimster, show_orient_label=True)
        draw_stereo_measurements(ctx, [ech], r_mini, orientation, point_size=0.18 * dimster / 10.0)


def build_zijderveld_figure(
    ech: SelectedSample,
    orientation: int = 1,
    fits: Optional[List] = None,
    fig: Optional[Figure] = None,
    auto_axis: bool = True,
    show_stereo: bool = True,
) -> Figure:
    """Enveloppe matplotlib de `draw_zijderveld`, pour app.py (widget
    FigureCanvasTkAgg unique). Le rendu passe par `PlotContext`, exactement
    comme pour `build_stereo_figure`."""
    if fig is None:
        fig = Figure(figsize=(6.0, 6.0), dpi=100)
    else:
        fig.clear()
    ax = fig.add_subplot(111)
    ctx = PlotContext(ax)
    ctx.clear()
    ctx.plot(0.0, 0.0, -3)

    draw_zijderveld(ctx, ech, orientation, fits, auto_axis, show_stereo)

    ax.relim()
    # marge verticale reduite (defaut matplotlib 5%) : moins d'espace vide
    # au-dessus du stereo et sous le nom de l'echantillon, demande explicite
    # de l'utilisateur - la marge horizontale reste au defaut.
    ax.margins(y=0.01)
    ax.autoscale_view()
    fig.tight_layout()
    return fig


def build_zijderveld_stereo_results_figure(
    ech: SelectedSample,
    result,
    orientation: int = 1,
    fig: Optional[Figure] = None,
) -> Figure:
    """Zijderveld + Stereo Results (avec grand cercle) cote a cote - demande
    explicite utilisateur ("in data+interpretation, when there is a plane,
    keep the plot zijderveld+stereo and plot the great circle on the
    stereo") : le Fortran d'origine (visres, plotorthog.f:107-121) n'affiche
    QUE `sterres` pour un resultat de plan (cat1='P') - PAS `zijder` - un
    choix delibere du programme d'origine, pas un bug de portage (verifie
    contre le source). Extension demandee au-dela du Fortran : voir les
    deux ensemble aide a confirmer visuellement l'ajustement de plan contre
    le trajet de desaimantation brut.

    `result` accepte soit un FitResult unique (usage historique -
    afficher_visres, un seul resultat archive a la fois), soit une LISTE de
    FitResult - demande explicite utilisateur pour l'auto-interpretation
    ("in data + interpretation, the original plot is lost; the data were
    shown on a zijderveld and stereo. I was asking to plot the plane on
    the stereo with the data when a result was a plane and keep the
    zijderveld") : ouvrir_autointerpretation_dialog superpose jusqu'a 2
    suggestions (primary+secondary) AVANT toute decision de sauvegarde, et
    utilisait jusque-la le Zijderveld seul (build_zijderveld_figure, encart
    stereo integre limite aux points bruts - JAMAIS de grand cercle, voir
    draw_zijderveld/draw_stereo_measurements) : une suggestion de PLAN y
    etait donc invisible. Memes panneaux que l'usage a un seul resultat,
    juste sur la liste complete plutot que sur `[result]`.

    Panneau gauche : Zijderveld de `ech` (`fits=results` pour l'etiquette
    "fit between X and Y" - le petit encart stereo integre est DESACTIVE
    ici, `show_stereo=False`, redondant avec le panneau droit qui montre le
    MEME reseau en plus grand et avec le(s) grand(s) cercle(s) de plan,
    absent(s) de l'encart). Panneau droit : Stereo Results
    (draw_stereo_results, meme fonction que `build_stereo_results_figure`)
    sur `results` - trace un grand cercle pour chaque resultat de plan
    (cat1=='P'), un point pour chaque ligne/direction (L/f)."""
    results = result if isinstance(result, list) else [result]
    if fig is None:
        fig = Figure(figsize=(11.0, 6.0), dpi=100)
    else:
        fig.clear()

    ax_zij = fig.add_subplot(121)
    ctx_zij = PlotContext(ax_zij)
    ctx_zij.clear()
    ctx_zij.plot(0.0, 0.0, -3)
    draw_zijderveld(ctx_zij, ech, orientation, results, True, False)
    ax_zij.relim()
    ax_zij.margins(y=0.01)
    ax_zij.autoscale_view()

    ax_ster = fig.add_subplot(122)
    ctx_ster = PlotContext(ax_ster)
    ctx_ster.clear()
    ctx_ster.plot(0.0, 0.0, -3)
    dimster = 12.0 * 1.5
    point_size = (0.18 * dimster) / 10.0
    r = draw_stereo_net(ctx_ster, orientation, dimster=dimster)
    draw_stereo_results(ctx_ster, results, r, orientation, point_size=point_size, nbech=1)
    ax_ster.relim()
    ax_ster.autoscale_view()
    ax_ster.set_title("Stereo Results")

    fig.tight_layout()
    return fig
