"""
Courbe de désaimantation (intensité vs étape) : port matplotlib idiomatique
(pas plotlib.PlotContext, a la demande) de `xygraph` (plotXY.f,
Starmac_AWE_v22), limite au cas courant (courbe de decroissance normalisee
au premier point - M/M0 - ou a chaque composante |X|/|Y|/|Z|). L'axe des
composantes n'est PAS corrige d'orientation dans le Fortran d'origine (pas
d'appel corfor/corpen dans xygraph) : on reste fidele a ce comportement.

Non porte : le mode acquisition IRM/champ retour (bloc `333` du Fortran,
calcul de Hcr, derivee du log) - fonctionnalite distincte, plus avancee.

Amelioration (pas dans le Fortran) : quand la selection melange des
echantillons desaimantes en AF et au thermique, deux sous-graphes distincts
sont traces sur la meme figure (unites d'etape incompatibles - mT vs degC -
donc un seul axe commun n'aurait pas de sens) plutot qu'un seul trace
melangeant les deux echelles.
"""

import math
from dataclasses import replace as _dc_replace
from typing import List, Optional

from matplotlib.figure import Figure
from matplotlib import colormaps

from selection import SelectedSample, split_experiments, experiment_kind

_UNITS = {"A": "mT", "F": "mT", "D": "°C", "S": "°C", "C": "hr"}
_AF_CODES = {"A", "F"}
_THERMAL_CODES = {"D", "S", "T", "K"}
# acquisition d'aimantation remanente isotherme (IRM) - pas une courbe de
# desaimantation, exclue du XYgraph (demande explicite de l'utilisateur) ;
# voir le menu dedie "Plot IRM" (irm.py) pour ces donnees.
_IRM_CODES = {"I"}
_COMPONENT_LABELS = {"total": "Mtot", "x": "|X|", "y": "|Y|", "z": "|Z|"}
_COLOR_CYCLE = colormaps["tab10"].colors


def _relevant_runs(selected: List[SelectedSample]):
    """Un echantillon peut enchainer plusieurs protocoles dans sa liste de
    mesures (ex. AF demag, puis acquisition IRM, puis desaimantation
    thermique de cette IRM acquise - cas reel signale par l'utilisateur,
    ex. 19DN1607B) : decoupe chaque echantillon en experiences distinctes
    (selection.split_experiments) et retourne (runs_normaux,
    runs_irm_demag) - 'F' (AF) et 'D' (thermique de la NRM) d'un cote,
    'D_IRM' (desaimantation thermique d'une IRM acquise, PAS de la NRM -
    ex. "0DI"/"130DI"...) de l'autre, pour affichage distinct dans le
    panneau thermique (demande explicite de l'utilisateur - visible mais
    jamais confondu avec de la vraie desaimantation thermique). 'I'
    (acquisition IRM) reste exclu (voir irm.py). Chaque run gardee devient
    un "echantillon virtuel" (meme id/orientation/volume/norme, mesures =
    ce run seul) pour que _sample_demag_code/_sample_values restent
    inchangees."""
    normal_runs, irm_demag_runs = [], []
    for ech in selected:
        for run in split_experiments(ech.mesures):
            kind = experiment_kind(run)
            if kind in ("F", "D"):
                normal_runs.append(_dc_replace(ech, mesures=run))
            elif kind == "D_IRM":
                irm_demag_runs.append(_dc_replace(ech, mesures=run))
    return normal_runs, irm_demag_runs


def _sample_demag_code(ech: SelectedSample) -> Optional[str]:
    """Code de desaimantation dominant de l'echantillon : le premier point
    est presque toujours 'N' (RMN initiale, avant tout pas de
    desaimantation), donc PAS representatif du type AF/thermique - on
    cherche le premier code de pas reel (le premier different de 'N')."""
    for m in ech.mesures:
        if m.cod1 != "N":
            return m.cod1
    return ech.mesures[0].cod1 if ech.mesures else None


def _sample_values(ech: SelectedSample, component: str):
    # m.etape est deja la valeur physique reelle (mT pour AF, degC pour
    # thermique - plus d'echelle Oersted a compenser ici, voir
    # testlect.Measurement) - meme colonne pour les deux cas desormais.
    steps = [m.etape for m in ech.mesures]
    if component == "x":
        vals = [abs(m.x) for m in ech.mesures]
    elif component == "y":
        vals = [abs(m.y) for m in ech.mesures]
    elif component == "z":
        vals = [abs(m.z) for m in ech.mesures]
    else:
        vals = [math.sqrt(m.x ** 2 + m.y ** 2 + m.z ** 2) for m in ech.mesures]
    return steps, vals


def _normalize_ref(vals: List[float]) -> float:
    """Reference de normalisation (mag@step0, comme le Fortran) - sauf si un
    pas ulterieur depasse 2x ce point, auquel cas le maximum de l'echantillon
    est utilise a la place (demande explicite de l'utilisateur : "when a
    specimen has a wrong strong data the normalisation of the plot with
    several samples go wrong... when mag@step(i)/mag@step0 > 2, better use
    maximum value of this specimen"). Un point aberrant isole normalise a
    ~1.0 au lieu de partir a 2x+ sinon ecrase (via ymax partage, cf.
    _plot_group) l'echelle Y commune de TOUS les echantillons du graphe."""
    ref = vals[0] if vals[0] else 1.0
    vmax = max(vals) if vals else ref
    return vmax if vmax > 2.0 * ref else ref


def _pick_points(ax, steps, norm, specimen_id: str) -> None:
    """Un point interactif INVISIBLE par mesure tracee (meme mecanisme que
    plotlib.PlotContext.pick_point, mais directement sur `ax` - ce module
    reste volontairement hors PlotContext, voir docstring en tete de
    fichier) - clic-pour-info cote app.py (`_on_plot_pick`), demande
    explicite utilisateur ("in XY plots and IRM plots, can we click on a
    point to get the sample number")."""
    for x, y in zip(steps, norm):
        artist = ax.scatter([x], [y], s=250, alpha=0, picker=True, zorder=1000)
        artist._starmac_pick_kind = "xygraph_specimen"
        artist._starmac_pick_data = specimen_id


def _plot_group(
    ax, samples: List[SelectedSample], component: str, title: str,
    irm_demag: Optional[List[SelectedSample]] = None,
) -> None:
    demag_code = ""
    ymax = 1.0
    for i, ech in enumerate(samples):
        if not ech.mesures:
            continue
        steps, vals = _sample_values(ech, component)
        ref = _normalize_ref(vals)
        norm = [v / ref for v in vals]
        ymax = max(ymax, max(norm))
        color = _COLOR_CYCLE[i % len(_COLOR_CYCLE)]
        ax.plot(steps, norm, "o-", color=color, markerfacecolor=color,
                 markersize=4, linewidth=1, label=ech.id)
        _pick_points(ax, steps, norm, ech.id)
        demag_code = demag_code or _sample_demag_code(ech)

    # desaimantation thermique d'une IRM deja acquise (code cod2='I' sur
    # les etapes D/S/T/K, ex. "0DI"/"130DI"...) - distincte de la
    # desaimantation thermique de la NRM (ligne pleine ci-dessus) : trace
    # en pointilles avec un marqueur different, couleurs continuant le
    # cycle pour rester distinguable, libelle suffixe "(IRM demag)" -
    # demande explicite de l'utilisateur.
    n_samples = len(samples)
    for j, ech in enumerate(irm_demag or []):
        if not ech.mesures:
            continue
        steps, vals = _sample_values(ech, component)
        ref = _normalize_ref(vals)
        norm = [v / ref for v in vals]
        ymax = max(ymax, max(norm))
        color = _COLOR_CYCLE[(n_samples + j) % len(_COLOR_CYCLE)]
        ax.plot(steps, norm, "s--", color=color, markerfacecolor="none",
                 markersize=4, linewidth=1, label=f"{ech.id} (IRM demag)")
        _pick_points(ax, steps, norm, ech.id)
        demag_code = demag_code or _sample_demag_code(ech)

    unit = _UNITS.get(demag_code, "")
    ax.set_xlabel("Demagnetization step" + (f" ({unit})" if unit else ""))
    ax.set_ylabel(f"{_COMPONENT_LABELS.get(component, component)} / "
                   f"{_COMPONENT_LABELS.get(component, component)}₀ (normalized scale)")
    ax.set_title(title)
    ax.set_ylim(0.0, ymax * 1.05)
    top_tick = max(10, int(math.ceil(ymax * 10)))
    ax.set_yticks([i / 10 for i in range(0, top_tick + 1, 2)])
    ax.axhline(1.0, color="0.6", linewidth=0.8, linestyle="--", zorder=0)
    nb_curves = len(samples) + len(irm_demag or [])
    if nb_curves and nb_curves <= 15:
        ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="0.7", linewidth=0.8, zorder=0)


def has_mixed_demag(selected: List[SelectedSample]) -> bool:
    """True si build_xygraph_figure va tracer 2 sous-graphes empiles (AF +
    thermique) plutot qu'un seul - AF present ET (thermique-de-la-NRM OU
    thermique-d'une-IRM-acquise present), par experience et non par
    echantillon entier (un seul echantillon peut a lui seul declencher
    True s'il enchaine plusieurs protocoles). Sert a l'appelant (app.py)
    pour choisir une figure plus large avant d'appeler
    build_xygraph_figure."""
    kinds = {experiment_kind(run) for ech in selected for run in split_experiments(ech.mesures)}
    return "F" in kinds and ("D" in kinds or "D_IRM" in kinds)


def build_xygraph_figure(
    selected: List[SelectedSample],
    component: str = "total",
    fig: Optional[Figure] = None,
) -> Figure:
    """Equivalent (cas courant) de `xygraph` : intensite normalisee au
    premier point (M/M0) en fonction de l'etape, une courbe par echantillon.
    `component` : 'total' (moment total, defaut), 'x', 'y' ou 'z'.

    Si la selection contient a la fois des echantillons desaimantes en AF
    (cod1 'A'/'F') et au thermique (cod1 'D'/'S'/'T'/'K'), deux sous-graphes
    cote a cote sont traces (un par type) ; sinon un seul trace comme avant.

    Amelioration (pas dans le Fortran, qui colore par COMPOSANTE via
    `newpen(ipltp1+3)` - toutes les courbes d'un meme trace partagent donc
    la meme couleur et ne se distinguent que par la legende/l'empilement) :
    une couleur distincte par ECHANTILLON, plus lisible avec plusieurs
    courbes superposees. L'axe Y normalise est fixe a [0, ~1.05] avec un
    repere explicite a 1.0 (le niveau M/M0 du tout premier point), au lieu
    de l'auto-echelle matplotlib qui pouvait legerement rogner le point de
    depart selon les echantillons affiches."""
    if fig is None:
        fig = Figure(figsize=(6.0, 4.5), dpi=100)
    else:
        fig.clear()

    normal_runs, irm_demag_runs = _relevant_runs(selected)
    af_samples = [r for r in normal_runs if _sample_demag_code(r) in _AF_CODES]
    th_samples = [r for r in normal_runs if _sample_demag_code(r) in _THERMAL_CODES]
    # une desaimantation thermique d'IRM (D_IRM) sans aucune vraie
    # thermique-de-la-NRM compte quand meme comme "thermique" ici, pour
    # qu'elle reste visible (demande explicite de l'utilisateur) plutot
    # que silencieusement omise faute de panneau pour l'accueillir.
    show_af = bool(af_samples)
    show_th = bool(th_samples or irm_demag_runs)

    if show_af and show_th:
        # Empile verticalement (pas cote a cote) : le panneau graphique de
        # l'app a une largeur FIXE (voir graph_frame dans app.py) - une
        # figure plus large que haute y serait tronquee, la moitie du trace
        # invisible hors de la zone visible sans agrandir la fenetre.
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        _plot_group(ax1, af_samples, component, "AF demagnetization")
        _plot_group(ax2, th_samples, component, "Thermal demagnetization", irm_demag=irm_demag_runs)
    elif show_th:
        ax = fig.add_subplot(111)
        _plot_group(ax, th_samples, component, "Thermal demagnetization", irm_demag=irm_demag_runs)
    else:
        ax = fig.add_subplot(111)
        _plot_group(ax, af_samples, component, "Demagnetization curve")

    fig.tight_layout()
    return fig
