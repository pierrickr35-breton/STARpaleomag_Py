"""
Controles supplementaires de "Evaluate interpretations..." tires de
l'interpretation automatique - demande explicite utilisateur (2026-11-08 :
"reexaminer la routine du menu evaluation interpretations en prenant en
compte les progres des interpretations automatiques").

Module separe d'interpretation_quality.py : auto_interpretation importe
deja interpretation_quality (import circulaire sinon).

Rien n'est modifie (ni donnees, ni resultats) : avertissements seulement.

Specimen (ajoute au rapport de chaque resultat) :
  1. palier aberrant isole DANS l'intervalle (auto_interpretation.
     spurious_steps, meme test que l'automatique - ex. 14NQ1708B 350 degC)
     -> note, au mieux "marginal" ;
  2. ajustement a cheval sur un changement de traitement ou un retour a
     un palier plus bas (auto_interpretation._treatment_breaks - ex.
     14NQ1708B 140 degC -> 120 mT) -> "poor" ;
  3. comparaison avec la proposition automatique du site (propose_
     components_for_site sur une COPIE des specimens) : angle entre le
     resultat et la composante automatique de meme lettre, note si
     > AUDIT_AUTO_DIFF_DEG (information, pas de degradation).
Site (section a part) :
  4. resultats a plus de 2 ecarts-types angulaires de la moyenne du site
     (calcul.split_outliers, meme regle que la revue des moyennes) ;
  5. droites a plus de 45 deg de la moyenne : grand cercle de
     remplacement passant pres de la moyenne, ou rejet (auto_
     interpretation._far_lines_to_planes) ;
  6. meme intervalle de desaimantation, composante differente (regle de
     l'utilisateur : "au sein d'un site tous les vecteurs dans le ~ meme
     intervalle de desaimantation devrait avoir la meme categorie de
     composante") ;
  7. specimens du site sans aucun resultat, pour lesquels l'automatique
     propose une composante ;
  8. aimantation parasite AF (auto_interpretation.af_spurious_assess) :
     niveau franc / leger ; les resultats AF qui atteignent la coupure d'un
     cas franc sont notes et ramenes au mieux a "marginal".
"""

import copy
import math
from types import SimpleNamespace
from typing import Dict, List, Optional, Tuple

import auto_interpretation as ai
from calcul import FitResult, fisher_from_results, split_outliers, _unit_vector
from interpretation_quality import QualityReport, evaluate_result, _GRADE_ORDER
from selection import zijderveld_measurements
from datatools import reduce_grm_before_fit

AUDIT_AUTO_DIFF_DEG = 15.0     # ecart a la proposition automatique signale
AUDIT_SAME_INTERVAL = 0.5      # recouvrement (fraction du plus court) = "meme intervalle"
_LINE_MAD, _PLANE_MAD = 15.0, 25.0  # seuils d'archivage (fit_line/fit_plane)


def _degrade(rep: QualityReport, grade: str) -> None:
    if _GRADE_ORDER.get(grade, 99) < _GRADE_ORDER.get(rep.grade, 99):
        rep.grade = grade


def _indices(mesures: list, res: FitResult) -> Optional[Tuple[int, int]]:
    """Indices 0-based du premier/dernier palier du resultat (premiere
    occurrence de l'etape, comme interpretation_quality._find_jdeb_jfin)."""
    i = next((k for k, m in enumerate(mesures) if m.etape == res.step_first), None)
    j = next((k for k, m in enumerate(mesures) if m.etape == res.step_last), None)
    if i is None or j is None or j < i:
        return None
    return i, j


def _specimen_checks(rep: QualityReport, res: FitResult, ech) -> None:
    mesures = zijderveld_measurements(ech)
    ij = _indices(mesures, res)
    if ij is None:
        return
    i, j = ij
    proxy = copy.copy(ech)
    proxy.mesures = list(mesures)
    inside = {id(m) for m in mesures[i:j + 1]}
    for m, why in ai.spurious_steps(proxy):
        if id(m) in inside:
            rep.notes.append(
                f"step {m.cod1} {m.etape:g} inside the interval looks like a measurement error ({why}) "
                f"- refit without it (Remove step)")
            _degrade(rep, "marginal")
    for k in sorted(ai._treatment_breaks(mesures)):
        # la NRM (palier 0) est le depart de tout traitement : N 0 -> F 5
        # n'est pas un changement de traitement
        if mesures[k].etape == 0:
            continue
        if i <= k < j:
            a, b = mesures[k], mesures[k + 1]
            rep.notes.append(
                f"the fit spans a change of treatment ({a.cod1} {a.etape:g} -> {b.cod1} {b.etape:g}) "
                f"- a result has a single demagnetization code")
            _degrade(rep, "poor")
            break


def _angle(dec1, inc1, dec2, inc2) -> float:
    c = float(_unit_vector(dec1, inc1) @ _unit_vector(dec2, inc2))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _sg_letter(sg) -> str:
    return sg.component or ("A" if sg.label == "primary" else "B")


def _cod_at(ech, step) -> str:
    m = next((m for m in zijderveld_measurements(ech) if m.etape == step), None)
    return m.cod1 if m is not None else ""


def audit_results(results: List[FitResult], donnees) -> Tuple[List[QualityReport], List[str]]:
    """Rapports specimen (interpretation_quality.evaluate_result + controles
    1-3) et lignes de la section site (controles 4-7)."""
    by_id = {s.id: s for s in donnees}
    reports: List[QualityReport] = []
    rep_of: Dict[int, QualityReport] = {}
    for res in results:
        rep = evaluate_result(res, donnees)
        if rep is None:
            continue
        ech = by_id.get(res.id)
        if ech is not None and rep.grade not in ("unrecomputable",):
            _specimen_checks(rep, res, ech)
        reports.append(rep)
        rep_of[id(res)] = rep

    # resultats specimen par site - un .pmagres accumule les reinterpretations
    # successives d'un meme specimen : pour les controles de site, seul le
    # DERNIER resultat de chaque (specimen, composante) compte
    latest: Dict[Tuple[str, str], FitResult] = {}
    for res in results:
        if res.cat1 in ("L", "P"):
            latest[(res.id, res.component)] = res
    per_site: Dict[str, List[FitResult]] = {}
    for res in latest.values():
        if res.cat1 not in ("L", "P") or id(res) not in rep_of:
            continue
        ech = by_id.get(res.id)
        site = (getattr(ech, "magic_site", "") or "").strip() if ech is not None else ""
        if site:
            per_site.setdefault(site, []).append(res)

    lines: List[str] = []
    ai.set_dataset_nrm_reference([e for e in donnees if getattr(e, "mesures", None)])
    for site in sorted(per_site):
        site_res = per_site[site]
        site_lines: List[str] = []
        echs = [s for s in donnees if (getattr(s, "magic_site", "") or "").strip() == site and s.mesures]

        # 8 : aimantation parasite AF (hauts champs mieux groupes en repere
        # echantillon qu'in situ) - resultats AF dont l'intervalle atteint
        # la coupure signales en plus
        afs = ai.af_spurious_assess(echs) if ai.AUTO_AF_SPURIOUS_CHECK else None
        if afs is not None:
            site_lines.append("  " + ai.format_af_spurious(site, afs, 0).split(None, 1)[1]
                              .replace(" left out (0 measurements)", " are suspect"))
            if afs.cut is not None:
                for res in site_res:
                    if res.id in by_id and _cod_at(by_id[res.id], res.step_last) == "F" \
                            and res.step_last >= afs.cut and id(res) in rep_of:
                        rep_of[id(res)].notes.append(
                            f"the fit reaches {res.step_last:g} mT, inside the high AF steps of this site "
                            f"suspected of spurious AF magnetization (from {afs.cut:g} mT)")
                        _degrade(rep_of[id(res)], "marginal")

        # 3 + 7 : proposition automatique du site, sur des COPIES
        suggestions = {}
        try:
            copies = [copy.deepcopy(e) for e in echs]
            # meme preparation que les menus d'interpretation automatique
            # (app._reduce_grm_before_fit : GRM methode 1), sur les copies
            reduce_grm_before_fit(copies, method=1)
            if afs is not None and afs.cut is not None:
                for c in copies:  # comme automatic_results : paliers suspects retires
                    c.mesures[:] = [m for m in c.mesures if not (m.cod1 == "F" and m.etape >= afs.cut)]
            suggestions = ai.propose_components_for_site(copies)
        except Exception as exc:  # un site atypique ne doit pas bloquer l'audit
            site_lines.append(f"  automatic interpretation failed for this site ({exc})")
        for res in site_res:
            sgs = [sg for sg in suggestions.get(res.id, []) if _sg_letter(sg) == res.component]
            if not sgs:
                continue
            sg = sgs[0]
            rep = rep_of[id(res)]
            same_kind = (sg.kind == "plane") == (res.cat1 == "P")
            if not same_kind:
                rep.notes.append(
                    f"the automatic interpretation fits a {sg.kind} here "
                    f"({sg.step_first:g}-{sg.step_last:g}), not a {'plane' if res.cat1 == 'P' else 'line'}")
                continue
            ang = _angle(res.dec, res.inc, sg.dec, sg.inc)
            if res.cat1 == "P":
                ang = min(ang, 180.0 - ang)
            if ang > AUDIT_AUTO_DIFF_DEG:
                rep.notes.append(
                    f"{ang:.0f} deg from the automatic {sg.kind} {sg.step_first:g}-{sg.step_last:g} "
                    f"(MAD {sg.mad:.1f}) - worth a second look")
        interpreted = {r.id for r in site_res}
        for e in echs:
            if e.id in interpreted:
                continue
            sgs = [sg for sg in suggestions.get(e.id, [])
                   if sg.mad <= (_PLANE_MAD if sg.kind == "plane" else _LINE_MAD)]
            if sgs:
                txt = "; ".join(f"{_sg_letter(sg)} {sg.kind} {sg.step_first:g}-{sg.step_last:g} MAD {sg.mad:.1f}"
                                for sg in sgs)
                site_lines.append(f"  {e.id}: not interpreted - automatic proposal {txt}")

        # 4 + 5 : par composante
        by_comp: Dict[str, List[FitResult]] = {}
        for res in site_res:
            by_comp.setdefault(res.component, []).append(res)
        for comp, comp_res in sorted(by_comp.items()):
            if len(comp_res) < ai.AUTO_MIN_SPECIMENS:
                continue
            ech_of = {r.id: by_id[r.id] for r in comp_res if r.id in by_id}
            sg_of = {}
            for r in comp_res:
                ech = ech_of.get(r.id)
                if ech is None:
                    continue
                ij = _indices(zijderveld_measurements(ech), r)
                if ij is not None:
                    sg_of[id(r)] = SimpleNamespace(index_first=ij[0])
            out = SimpleNamespace(specimens=list(comp_res), outliers=[], far_to_plane=[])
            kept = ai._far_lines_to_planes(list(comp_res), ech_of, sg_of, out)
            far_ids = set()
            for _s, _c, sid, d, s1, s2 in out.far_to_plane:
                far_ids.add(sid)
                site_lines.append(
                    f"  {sid} ({comp}): line {d:.0f} deg from the site mean - the great circle "
                    f"{s1:g}-{s2:g} passes near the mean: replace the line by it?")
            for _s, _c, sid, d in out.outliers:
                far_ids.add(sid)
                site_lines.append(
                    f"  {sid} ({comp}): line {d:.0f} deg from the site mean and no great circle near "
                    f"the mean - reject the specimen?")
            rest = [r for r in comp_res if r.id not in far_ids]
            if len(rest) >= ai.AUTO_MIN_SPECIMENS:
                _kept, rejected, asd = split_outliers(rest, orientation=2, combine_antipodal=True)
                for r, d in rejected:
                    site_lines.append(
                        f"  {r.id} ({comp}): {d:.1f} deg from the site mean, more than 2 sd "
                        f"({2 * asd:.1f} deg) - remove it from the mean?")
                    if id(r) in rep_of:
                        rep_of[id(r)].notes.append(f"more than 2 sd from the site {comp} mean ({d:.1f} deg)")

        # 6 : meme intervalle, composante differente (droites, meme traitement)
        flagged = set()
        lines_only = [r for r in site_res if r.cat1 == "L" and r.id in by_id]
        for r in lines_only:
            cod = _cod_at(by_id[r.id], r.step_first)
            same, other = 0, []
            for o in lines_only:
                if o is r or o.id == r.id or _cod_at(by_id[o.id], o.step_first) != cod:
                    continue
                lo, hi = max(r.step_first, o.step_first), min(r.step_last, o.step_last)
                short = min(r.step_last - r.step_first, o.step_last - o.step_first)
                if short <= 0 or hi - lo < AUDIT_SAME_INTERVAL * short:
                    continue
                if o.component == r.component:
                    same += 1
                else:
                    other.append(o)
            if other and len(other) > same and r.id not in flagged:
                flagged.add(r.id)
                comps = sorted({o.component for o in other})
                site_lines.append(
                    f"  {r.id} ({r.component} {r.step_first:g}-{r.step_last:g}): same demagnetization "
                    f"interval as {len(other)} result(s) labelled {'/'.join(comps)} in this site "
                    f"(only {same} labelled {r.component})")

        if site_lines:
            lines.append(f"Site {site}:")
            lines.extend(site_lines)
    return reports, lines


def format_site_section(lines: List[str]) -> str:
    if not lines:
        return "Site-level checks: nothing to report."
    return "Site-level checks (nothing is changed - suggestions only):\n" + "\n".join(lines)
