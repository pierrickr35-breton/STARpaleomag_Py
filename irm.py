"""
Acquisition d'aimantation remanente isotherme (IRM) : courbe d'acquisition
normalisee (M/Mmax) et spectre de coercivite (derivee dM/d(log H)), port
matplotlib idiomatique (numpy uniquement - scipy non installe dans cet
environnement) des scripts GMT de l'utilisateur
(Scripts_IRM_GMT/Script_GMT_single_IRM.txt et
Script_GMT_several_IRM.txt) :

    log10(champ) -> reechantillonnage sur grille reguliere -> lissage
    -> derivee par rapport a log10(champ) -> reconversion en champ
    lineaire (mT) pour l'affichage.

Mesures identifiees par cod1='I' (voir xygraph._IRM_CODES, deja exclues du
XYgraph - une acquisition IRM n'est pas une courbe de desaimantation).
`etape` est stocke DIRECTEMENT en mT pour ce code (PAS de conversion
Oersted->mT comme pour l'AF 'F' - confirme par l'utilisateur).

cod2='R' = acquisition en champ inverse (backfield, ex. pour determiner
Bcr) - intensite DEcroissante avec l'etape, trace separement des points
d'acquisition directe (cod2 != 'R', ex. '0'/'+'/'='/'Z' rencontres dans
des fichiers reels). Le lissage/derivee n'est PAS verifie sur donnees
reelles (aucun exemple de sequence IRM/backfield fourni pour comparaison
octet-pres, contrairement aux autres modules de ce portage) - a confirmer
visuellement sur un vrai jeu de donnees avant un usage scientifique.
"""

import math
from typing import List, Optional, Tuple

import numpy as np
from matplotlib.figure import Figure
from matplotlib import colormaps

from selection import SelectedSample

_COLOR_CYCLE = colormaps["tab10"].colors


def _sample_irm_points(
    ech: SelectedSample,
) -> Tuple[List[float], List[float], List[float], List[float]]:
    """(etapes_acquisition, valeurs_acquisition, etapes_backfield,
    valeurs_backfield) - un point par mesure cod1='I', tries par etape
    croissante ; cod2='R' = backfield (voir docstring module)."""
    forward = []
    backfield = []
    for m in ech.mesures:
        if m.cod1 != "I":
            continue
        val = math.sqrt(m.x ** 2 + m.y ** 2 + m.z ** 2)
        (backfield if m.cod2 == "R" else forward).append((m.etape, val))
    forward.sort()
    backfield.sort()
    fsteps = [s for s, _ in forward]
    fvals = [v for _, v in forward]
    bsteps = [s for s, _ in backfield]
    bvals = [v for _, v in backfield]
    return fsteps, fvals, bsteps, bvals


def _smooth_and_derive(
    steps: List[float], vals: List[float],
    log_step: float = 0.05, smooth_window: int = 5,
) -> Tuple[List[float], List[float]]:
    """Reechantillonne (steps,vals) sur une grille reguliere en
    log10(champ), lisse (moyenne glissante) puis derive par rapport a
    log10(champ) - reconvertit l'axe en champ lineaire (mT). Retourne
    (champ_mT, derivee). Equivalent numpy-only (pas de spline de lissage a
    tension comme `gmt sample1d -Fs<p>`, scipy non installe) de la meme
    intention que le script GMT d'origine : reduire le bruit point-a-point
    avant de deriver, pour un spectre de coercivite lisible."""
    if len(steps) < 3:
        return [], []
    steps_arr = np.asarray(steps, dtype=float)
    vals_arr = np.asarray(vals, dtype=float)
    mask = steps_arr > 0
    steps_arr, vals_arr = steps_arr[mask], vals_arr[mask]
    if len(steps_arr) < 3:
        return [], []

    log_field = np.log10(steps_arr)
    grid = np.arange(log_field.min(), log_field.max() + log_step, log_step)
    val_grid = np.interp(grid, log_field, vals_arr)

    w = smooth_window if smooth_window % 2 == 1 else smooth_window + 1
    if w > 1 and len(val_grid) > w:
        pad = w // 2
        padded = np.pad(val_grid, pad, mode="edge")
        val_smooth = np.convolve(padded, np.ones(w) / w, mode="valid")
    else:
        val_smooth = val_grid

    deriv = np.gradient(val_smooth, grid)
    field_mt = 10.0 ** grid
    return field_mt.tolist(), deriv.tolist()


def _pick_points(ax, steps, norm, specimen_id: str) -> None:
    """Un point interactif INVISIBLE par mesure tracee (meme mecanisme que
    plotlib.PlotContext.pick_point, mais directement sur `ax` - ce module
    reste hors PlotContext, comme xygraph.py) - clic-pour-info cote
    app.py (`_on_plot_pick`), demande explicite utilisateur ("in XY plots
    and IRM plots, can we click on a point to get the sample number")."""
    for x, y in zip(steps, norm):
        artist = ax.scatter([x], [y], s=250, alpha=0, picker=True, zorder=1000)
        artist._starmac_pick_kind = "irm_specimen"
        artist._starmac_pick_data = specimen_id


def build_irm_figure(
    selected: List[SelectedSample],
    fig: Optional[Figure] = None,
) -> Figure:
    """2 panneaux empiles (axe des champs en log10, comme
    `-JX12cl/6c` dans le script GMT d'origine) : acquisition IRM
    normalisee (M/Mmax) en bas, spectre de coercivite (derivee
    normalisee, avec le champ au pic annote) en haut - une couleur par
    echantillon, plusieurs echantillons superposes (pas empiles en pages
    separees comme le script GMT, a la demande de l'utilisateur)."""
    if fig is None:
        fig = Figure(figsize=(6.0, 8.0), dpi=100)
    else:
        fig.clear()

    ax_der = fig.add_subplot(211)
    ax_acq = fig.add_subplot(212)

    any_data = False
    for i, ech in enumerate(selected):
        fsteps, fvals, bsteps, bvals = _sample_irm_points(ech)
        if not fsteps:
            continue
        any_data = True
        color = _COLOR_CYCLE[i % len(_COLOR_CYCLE)]
        mmax = max(fvals) if max(fvals) else 1.0
        norm_f = [v / mmax for v in fvals]
        ax_acq.plot(fsteps, norm_f, "o-", color=color, markerfacecolor=color,
                    markersize=4, linewidth=1, label=ech.id)
        _pick_points(ax_acq, fsteps, norm_f, ech.id)

        if bsteps:
            norm_b = [v / mmax for v in bvals]
            ax_acq.plot(bsteps, norm_b, "s--", color=color, markerfacecolor="none",
                        markersize=4, linewidth=1)
            _pick_points(ax_acq, bsteps, norm_b, ech.id)

        field_mt, deriv = _smooth_and_derive(fsteps, norm_f)
        if deriv:
            dmax = max(deriv)
            norm_deriv = [d / dmax for d in deriv] if dmax else deriv
            ax_der.plot(field_mt, norm_deriv, "-", color=color, linewidth=1.2)
            if len(selected) <= 5:
                i_peak = int(np.argmax(deriv))
                ax_der.annotate(
                    f"{field_mt[i_peak]:.0f} mT",
                    xy=(field_mt[i_peak], norm_deriv[i_peak]),
                    xytext=(4, 4), textcoords="offset points", fontsize=8, color=color,
                )

    ax_acq.set_xscale("log")
    ax_acq.set_xlabel("Magnetizing field (mT)")
    ax_acq.set_ylabel("M / Mmax")
    ax_acq.set_title("IRM acquisition")
    ax_acq.grid(True, which="both", alpha=0.3)
    if any_data and len(selected) <= 15:
        ax_acq.legend(fontsize=8)

    ax_der.set_xscale("log")
    ax_der.set_xlabel("Magnetizing field (mT)")
    ax_der.set_ylabel("dM / d(log H), normalized")
    ax_der.set_title("Coercivity spectrum")
    ax_der.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    return fig


def has_irm_data(selected: List[SelectedSample]) -> bool:
    """True si au moins un echantillon de `selected` a une mesure IRM
    (cod1='I') - sert a app.py pour avertir si la selection n'en contient
    aucune avant d'ouvrir un graphique vide."""
    return any(any(m.cod1 == "I" for m in ech.mesures) for ech in selected)
