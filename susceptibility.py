"""
Susceptibilité vs étape : port matplotlib idiomatique (pas plotlib.PlotContext,
a la demande) de `suscep` (plotsuscep.f, Starmac_AWE_v22). Normalisation au
maximum (pas au premier point, contrairement a xygraph - c'est le
comportement du Fortran d'origine : `snorm` est initialise a sint(1) mais
toujours reecrase par le max de la serie dans la boucle qui suit).
"""

import math
from typing import List, Optional

from matplotlib.figure import Figure
from matplotlib import colormaps

from selection import SelectedSample

# R/V/P (paleointensite Thellier/IZZI) INCLUS en degC - demande explicite
# utilisateur ("when the cod1 is R, V, P, the unit is also °C") : cette
# fonction trace TOUTE la selection sans filtrage AF/thermique
# (contrairement a xygraph.py, ou un specimen R/V/P est deja ecarte en
# amont par experiment_kind - ici rien n'exclut ces mesures, donc sans
# elles le premier code rencontre pouvait rester sans unite affichee).
_UNITS = {"A": "mT", "F": "mT", "D": "°C", "S": "°C", "R": "°C", "V": "°C", "P": "°C", "C": "hr"}
_COLOR_CYCLE = colormaps["tab10"].colors


def build_susceptibility_figure(
    selected: List[SelectedSample],
    fig: Optional[Figure] = None,
) -> Figure:
    """Equivalent de `suscep` : susceptibilite (mesure.s) normalisee au
    maximum de la serie, en fonction de l'etape, une courbe par echantillon.

    Amelioration (pas dans le Fortran, voir xygraph.build_xygraph_figure
    pour le detail) : une couleur distincte par echantillon, et un axe Y
    normalise fixe a [0, ~1.05] avec un repere explicite a 1.0 (le niveau
    chi/chi_max)."""
    if fig is None:
        fig = Figure(figsize=(6.0, 4.5), dpi=100)
    else:
        fig.clear()
    ax = fig.add_subplot(111)

    demag_code = ""
    any_data = False
    ymax = 1.0
    for i, ech in enumerate(selected):
        pts = [(m.etape, m.s, m.cod1) for m in ech.mesures if m.s != 0]
        if len(pts) < 2:
            continue
        any_data = True
        steps, svals, codes = zip(*pts)
        smax = max(svals) if max(svals) else 1.0
        norm = [s / smax for s in svals]
        ymax = max(ymax, max(norm))
        color = _COLOR_CYCLE[i % len(_COLOR_CYCLE)]
        ax.plot(steps, norm, "o-", color=color, markerfacecolor=color,
                 markersize=4, linewidth=1, label=ech.id)
        demag_code = demag_code or codes[0]

    unit = _UNITS.get(demag_code, "")
    ax.set_xlabel(f"Demagnetization step" + (f" ({unit})" if unit else ""))
    ax.set_ylabel("Normalized susceptibility (χ/χmax)")
    ax.set_title("Susceptibility vs step")
    ax.set_ylim(0.0, ymax * 1.05)
    top_tick = max(10, int(math.ceil(ymax * 10)))
    ax.set_yticks([i / 10 for i in range(0, top_tick + 1, 2)])
    ax.axhline(1.0, color="0.6", linewidth=0.8, linestyle="--", zorder=0)
    if any_data:
        ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="0.7", linewidth=0.8, zorder=0)
    fig.tight_layout()
    return fig
