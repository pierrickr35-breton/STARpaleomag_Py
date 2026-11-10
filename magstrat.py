"""
Graphiques de magnetostratigraphie : declinaison et inclinaison en fonction
de la position stratigraphique (champ `stratigraphic_height` du .prmag,
metres, positif vers le haut) - demande explicite utilisateur ("un menu
specifique pour faire des graphiques de magnetostratigraphie en utilisant
la variable stratigraphic height ... a gauche Declination/Stratigraphic
height et a droite Inclination/stratigraphic height").

Deux sources de directions :
  - les resultats specimen (droites d'une composante, .pmagres) ;
  - la mesure a un palier de desaimantation donne (traitement "blanket").
Troisieme panneau au choix : latitude du VGP, ou angle de renversement
(angle a une direction normale de reference, GAD par defaut) - meme calcul
que `lismesdepth` (selection.list_measurements_depth, subroutine `angle` de
calcul.f) dans l'ancienne version Fortran.

Declinaisons ramenees dans [-90, 270[ comme dans lismesdepth (dec > 270 ->
dec - 360) : normales autour de 0, inverses autour de 180, sans coupure.
"""

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

from matplotlib.figure import Figure


@dataclass
class StratPoint:
    height: float
    dec: float
    inc: float
    specimen: str
    deviation: Optional[float] = None   # angle a la direction normale de reference (deg)
    vgp_lat: Optional[float] = None
    # directions des paliers individuels quand un intervalle est moyenne :
    # [(dec, inc, vgp_lat ou None, angle ou None)] - tracees en petits
    # cercles ouverts autour de la moyenne (demande explicite utilisateur :
    # "plotter les 3 ou 4 points de la desaimantation")
    steps: Optional[list] = None        # [dec, inc, vgp_lat, angle, libelle du palier]
    info: str = ""                      # origine du point (resultat, palier, moyenne de n paliers)


def wrap_dec(dec: float) -> float:
    dec = dec % 360.0
    return dec - 360.0 if dec >= 270.0 else dec


def angular_deviation(dec: float, inc: float, exp_dec: float, exp_inc: float) -> float:
    """Angle (deg) entre (dec, inc) et la direction attendue."""
    def u(d, i):
        d, i = math.radians(d), math.radians(i)
        return (math.cos(i) * math.cos(d), math.cos(i) * math.sin(d), math.sin(i))
    a, b = u(dec, inc), u(exp_dec, exp_inc)
    c = max(-1.0, min(1.0, sum(x * y for x, y in zip(a, b))))
    return math.degrees(math.acos(c))


def gad_inclination(lat: float) -> float:
    return math.degrees(math.atan(2.0 * math.tan(math.radians(lat))))


def polarity_of(p: StratPoint) -> Optional[str]:
    """'N' / 'R' : latitude du VGP si connue, sinon angle de renversement
    (< 90 deg = normal), sinon None."""
    if p.vgp_lat is not None:
        return "N" if p.vgp_lat >= 0 else "R"
    if p.deviation is not None:
        return "N" if p.deviation < 90.0 else "R"
    return None


def polarity_zones(points: List[StratPoint]) -> List[Tuple[float, float, str, int]]:
    """Zones (bas, haut, 'N'|'R'|'?', nb de niveaux). Polarite d'un niveau
    = majorite de ses specimens, '?' si egalite ; limites a mi-chemin entre
    niveaux de polarite opposee ; extremites aux niveaux extremes."""
    levels = {}
    for p in points:
        pol = polarity_of(p)
        if pol:
            levels.setdefault(p.height, []).append(pol)
    hs = sorted(levels)
    if not hs:
        return []
    pol = []
    for h in hs:
        n, r = levels[h].count("N"), levels[h].count("R")
        pol.append("N" if n > r else "R" if r > n else "?")
    zones = []
    start = 0
    for k in range(1, len(hs) + 1):
        if k == len(hs) or pol[k] != pol[start]:
            lo = hs[0] if start == 0 else 0.5 * (hs[start - 1] + hs[start])
            hi = hs[-1] if k == len(hs) else 0.5 * (hs[k - 1] + hs[k])
            zones.append((lo, hi, pol[start], k - start))
            start = k
    return zones


def build_magstrat_figure(
    points: List[StratPoint],
    fig: Optional[Figure] = None,
    title: str = "",
    lat: Optional[float] = None,
    third: Optional[str] = None,
    reference: Optional[Tuple[float, float]] = None,
    polarity: bool = True,
) -> Figure:
    """Dec | Inc | [latitude VGP (third='vgp') ou angle de renversement
    (third='angle')] | colonne de polarite (noir = normal, blanc = inverse,
    gris = indetermine ; zone d'un seul niveau en demi-largeur, convention
    usuelle). Axes y partages. `lat` : inclinaisons GAD +/- en tirets."""
    if fig is None:
        fig = Figure(figsize=(7.5, 8.5))
    pts = sorted(points, key=lambda p: p.height)
    h = [p.height for p in pts]
    zones = polarity_zones(pts) if polarity else []
    ratios = [3, 3] + ([3] if third else []) + ([0.8] if zones else [])
    axes = fig.subplots(1, len(ratios), sharey=True, gridspec_kw={"width_ratios": ratios})
    ax_d, ax_i = axes[0], axes[1]

    style = dict(color="0.35", lw=0.8, marker="o", ms=4.5, mfc="black", mec="black", zorder=3)
    sty_steps = dict(ls="none", marker="o", ms=5.0, mfc="none", mec="0.35", mew=0.8, zorder=2)
    sp = [(p.height, st) for p in pts for st in (p.steps or [])]

    def plot_steps(ax, k, transform=None):
        xs = [(transform(st[k]) if transform else st[k]) for _h, st in sp if st[k] is not None]
        hs = [hh for hh, st in sp if st[k] is not None]
        if xs:
            ax.plot(xs, hs, **sty_steps)

    plot_steps(ax_d, 0, wrap_dec)
    ax_d.plot([wrap_dec(p.dec) for p in pts], h, **style)
    ax_d.set_xlim(-90, 270)
    ax_d.set_xticks([-90, 0, 90, 180, 270])
    for x in (0, 180):
        ax_d.axvline(x, color="0.6", lw=0.7, ls="--", zorder=1)
    ax_d.set_xlabel("Declination (°)")
    ax_d.set_ylabel("Stratigraphic height (m)")

    plot_steps(ax_i, 1)
    ax_i.plot([p.inc for p in pts], h, **style)
    ax_i.set_xlim(-90, 90)
    ax_i.set_xticks([-90, -45, 0, 45, 90])
    ax_i.axvline(0, color="0.6", lw=0.7, zorder=1)
    if lat is not None:
        gi = gad_inclination(lat)
        for x in (gi, -gi):
            ax_i.axvline(x, color="0.6", lw=0.7, ls="--", zorder=1)
    ax_i.set_xlabel("Inclination (°)")

    col = 2
    if third == "vgp":
        ax_v = axes[col]; col += 1
        plot_steps(ax_v, 2)
        ax_v.plot([p.vgp_lat for p in pts], h, **style)
        ax_v.set_xlim(-90, 90)
        ax_v.set_xticks([-90, -45, 0, 45, 90])
        ax_v.axvline(0, color="0.6", lw=0.7, zorder=1)
        ax_v.set_xlabel("VGP latitude (°)")
    elif third == "angle":
        ax_a = axes[col]; col += 1
        plot_steps(ax_a, 3)
        ax_a.plot([p.deviation for p in pts], h, **style)
        ax_a.set_xlim(0, 180)
        ax_a.set_xticks([0, 45, 90, 135, 180])
        ax_a.axvline(90, color="0.6", lw=0.7, ls="--", zorder=1)
        ref = f" to {reference[0]:g}/{reference[1]:.0f}" if reference else ""
        ax_a.set_xlabel(f"Reversal angle{ref} (°)")

    if zones:
        ax_p = axes[col]
        ax_p.set_xlim(0, 1)
        for lo, hi, pol, nlev in zones:
            face = {"N": "black", "R": "white", "?": "0.7"}[pol]
            width = 0.5 if nlev == 1 else 1.0
            ax_p.fill_between([0, width], lo, hi, facecolor=face, edgecolor="black", lw=0.6)
        ax_p.set_xticks([])
        ax_p.set_xlabel("Polarity")
        # pas de cadre : seules les zones (bord noir) sont tracees, sinon la
        # marge blanche au-dela des niveaux extremes ressemble a une zone inverse
        for side in ("top", "bottom", "left", "right"):
            ax_p.spines[side].set_visible(False)
        ax_p.tick_params(left=False, right=False)

    # points cliquables (meme mecanisme que xygraph/irm : scatter invisible
    # portant _starmac_pick_kind/_starmac_pick_data, lu par
    # app._on_plot_pick) - demande explicite utilisateur ("peut-on cliquer
    # sur les points pour connaitre l'etape et le specimen")
    def pick_text(p, st=None):
        d, i = (st[0], st[1]) if st else (p.dec, p.inc)
        vg = st[2] if st else p.vgp_lat
        an = st[3] if st else p.deviation
        what = st[4] if st else p.info
        return (f"{p.specimen}  height {p.height:g} m  {what}  Dec {d:.1f}  Inc {i:.1f}"
                + (f"  VGP lat {vg:.1f}" if vg is not None else "")
                + (f"  angle {an:.1f}" if an is not None else ""))

    panels = [(ax_d, lambda p, st: wrap_dec(st[0] if st else p.dec)),
              (ax_i, lambda p, st: st[1] if st else p.inc)]
    if third == "vgp":
        panels.append((axes[2], lambda p, st: st[2] if st else p.vgp_lat))
    elif third == "angle":
        panels.append((axes[2], lambda p, st: st[3] if st else p.deviation))
    for ax, getx in panels:
        xs, ys, labels = [], [], []
        for p in pts:
            for st in [None] + list(p.steps or []):
                x = getx(p, st)
                if x is None:
                    continue
                xs.append(x)
                ys.append(p.height)
                labels.append(pick_text(p, st))
        if xs:
            art = ax.scatter(xs, ys, s=60, alpha=0, picker=True, zorder=1000)
            art._starmac_pick_kind = "magstrat_point"
            art._starmac_pick_data = labels

    for ax in axes:
        if ax is not (axes[col] if zones else None):
            ax.tick_params(direction="in", top=True, right=True, labelsize=8)
        ax.xaxis.label.set_size(9)
    ax_d.yaxis.label.set_size(9)
    if title:
        fig.suptitle(title, fontsize=10)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.94, bottom=0.08, wspace=0.2)
    return fig
