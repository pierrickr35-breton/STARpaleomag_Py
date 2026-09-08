"""
Menu "Anisotropy PmagPy" - demande explicite utilisateur ("I have never
been able to figure out how the anisotropy of remanent magnetization is
determined in Magic. Can you figure out from the PmagPy tools and add a
menu Anisotropy PmagPy"). 2e pipeline PARALLELE appelant le vrai pmagpy
(meme principe que paleointensity_magic.py - PAS un port, un appel direct
aux fonctions reelles de pmagpy.ipmag/pmagpy.pmag), pour repondre a la
question posee : comment MagIC determine reellement un tenseur d'AARM
(Anisotropy of Anhysteretic Remanent Magnetization) ou d'ATRM.

Trouve dans pmagpy (verifie en lisant le source live, pmagpy 4.5.2) :
- `ipmag.aarm_magic`/`ipmag.atrm_magic` : pipeline complet fichier
  measurements.txt -> specimens.txt, filtre par method_codes
  ('LP-AN-ARM'/'LP-AN-TRM'), soustrait une paire baseline/champ par
  position, puis appelle...
- `ipmag.get_matrix(n_pos)` : matrice de design A (moindres carres, Hext
  1963) pour n_pos in (6, 9, 15) - PAS un ajustement position-par-position
  comme le fait STARpaleomag_Py (A0), mais une inversion lineaire generale sur
  TOUTES les positions a la fois (surdeterminee des que n_pos>6),
  B = (AtA)^-1 At.
- `ipmag.calculate_aniso_parameters(K, n_pos)` : applique B au vecteur K
  (moments des n_pos positions, 3*n_pos valeurs), normalise par la TRACE
  (convention Jelinek : s1+s2+s3=1, CONTRAIREMENT au tenseur brut 'A0' de
  STARpaleomag_Py, en unites physiques), diagonalise (valeurs propres t1>=t2>=t3,
  vecteurs propres v1/v2/v3 = axes principaux), calcule `aniso_p` (degre
  d'anisotropie t1/t3), et - l'apport reel au-dela de ce que fait deja
  STARpaleomag_Py - le TEST F DE HEXT (`pmag.dohext`, Hext 1963) : compare le
  residu de l'ajustement moindres carres a une hypothese nulle
  d'isotropie, donnant F/F12/F23 + un indicateur qualite 'g'/'b' selon
  que F depasse le F critique - AUCUN equivalent dans le calcul natif
  STARpaleomag_Py (qui produit des variantes A0/A+/A-/A1-A6/B1-B6 a comparer a
  l'oeil, mais pas de test de significativite statistique formel).

VERIFIE : pour n_pos=6 (le seul mode que STARpaleomag_Py sait detecter via
detect_six_positions), l'ordre des positions de get_matrix(6) est EXACTEMENT
+X,+Y,+Z,-X,-Y,-Z (verifie via dir2cart sur les 6 directions codees en dur)
- correspond terme a terme a l'ordre X+,Y+,Z+,X-,Y-,Z- deja utilise par
calcul.compute_anisotropy_tensor. Sur un tenseur synthetique connu (bruit
realiste ajoute), le tenseur trace-normalise de pmagpy (`aniso_s`) est
identique (5 decimales) au tenseur 'A0' de STARpaleomag_Py normalise par sa propre
trace - confirme que les deux methodes calculent la MEME chose pour n_pos=6
(la moindre-carres se reduit a la demi-difference simple des lors qu'il y a
exactement autant de positions que d'inconnues), pmagpy ajoutant seulement
le test de significativite en plus.

BUG pmagpy CONFIRME (contourne, pas silencieusement ignore) :
`calculate_aniso_parameters` plante (`UnboundLocalError: sigma`) si le
residu de l'ajustement est EXACTEMENT nul (mesures parfaitement conformes
au tenseur ajuste, sans le moindre bruit) - `sigma` n'est assigne que dans
la branche `if S > 0`, jamais initialise avant. N'arrive normalement
jamais sur de vraies mesures (bruit instrumental toujours present), mais
intercepte ici par securite (voir compute_aarm_pmagpy)."""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pmagpy.ipmag as ipmag

from selection import Measurement

# Ordre EXACT des positions de get_matrix(6) - voir docstring module.
_POSITION_ORDER = ("X+", "Y+", "Z+", "X-", "Y-", "Z-")


@dataclass
class PmagpyAnisotropyResult:
    """Resultat de calculate_aniso_parameters, decode en champs numeriques
    (pmagpy retourne des chaines "a:b:c" pour s/v1/v2/v3 - voir docstring
    module)."""
    n_pos: int
    s: Tuple[float, float, float, float, float, float]  # s1,s2,s3,s4,s5,s6 (trace-normalises, s1+s2+s3=1)
    v1: Tuple[float, float, float]  # (valeur propre, dec, inc) - axe principal (max)
    v2: Tuple[float, float, float]  # axe intermediaire
    v3: Tuple[float, float, float]  # axe mineur (min)
    p: float  # degre d'anisotropie, t1/t3
    sigma: Optional[float]  # None si non calculable (voir bug ci-dessus)
    f_test: Optional[float]
    f12_test: Optional[float]
    f23_test: Optional[float]
    f_crit: Optional[float]
    quality: Optional[str]  # 'g'/'b' - None si non calculable


def _parse_triplet(text: str) -> Tuple[float, float, float]:
    a, b, c = text.split(":")
    return float(a), float(b), float(c)


def compute_aarm_pmagpy(
    positions: Dict[str, Measurement],
    holder: Optional[object] = None,
    n_pos: int = 6,
    nrm_mean: Optional[Tuple[float, float, float]] = None,
) -> PmagpyAnisotropyResult:
    """Construit `K` (vecteur des moments des n_pos positions, ordre
    _POSITION_ORDER pour n_pos=6) a partir de `positions` (voir
    calcul.detect_six_positions - MEME dict que celui deja utilise par
    compute_anisotropy_tensor, ligne de base porte-echantillon deja
    soustraite en amont si `holder` est fourni), puis appelle le vrai
    pmagpy (ipmag.calculate_aniso_parameters) - PAS une reimplementation.

    `nrm_mean` (x,y,z, MEME convention que calcul._nrm_mean_and_diag/
    AnisotropyComputation.nrm_mean - a passer par l'appelant, ex.
    `result.nrm_mean`) : BUG REEL corrige ici (demande explicite
    utilisateur "if the Ftest from PmagPy implemented correctly. All
    tensors are not significant and some are not that bad") - CONFIRME
    par simulation Monte-Carlo (tenseur synthetique connu + bruit
    realiste + contamination NRM constante sur les 6 positions) : la
    formule de Hext elle-meme (pmag.dohext/fcalc) est correcte et
    sensible (taux de detection proche de 100% pour un vrai signal
    anisotrope avec bruit realiste), MAIS le calcul natif 'A0'
    (compute_anisotropy_tensor) elimine la NRM residuelle (TRM/ARM
    PARTIELLE, cf. calcul._nrm_mean_and_diag) PAR CONSTRUCTION via la
    difference de paire (X+-X-)/2 (le terme NRM commun aux deux moities
    de la paire s'annule algebriquement, quelle que soit son ampleur) -
    alors que le moindres-carres de pmagpy (design lineaire sur les 6
    positions BRUTES, sans terme de decalage constant libre) n'a AUCUN
    moyen equivalent de l'ignorer : une contamination NRM constante,
    meme faible (10% du signal anisotrope suffit dans la simulation),
    fait chuter le taux de detection observe de ~100% a ~0%, sans
    rapport avec la force reelle de l'anisotropie - exactement le
    symptome signale ("tous les tenseurs 'not significant', certains pas
    si mauvais"). Corrige en soustrayant `nrm_mean` de chaque position
    AVANT de construire K, pour que le moindres-carres pmagpy voie la
    MEME quantite (deja immunisee contre la NRM residuelle) que le
    tenseur natif 'A0'. None (par defaut, retro-compatible) reproduit
    l'ancien comportement CONTAMINE - a eviter pour tout nouvel appel."""
    if n_pos != 6:
        raise NotImplementedError(
            "only n_pos=6 is wired to STARpaleomag_Py's own position detection "
            "(detect_six_positions) for now")

    # reimport local pour eviter une dependance circulaire avec calcul.py
    # (qui importe deja selection.py) - _position_vector est un utilitaire
    # generique, pas specifique a compute_anisotropy_tensor.
    from calcul import _position_vector

    idx_of = {"X+": 0, "X-": 1, "Y+": 2, "Y-": 3, "Z+": 4, "Z-": 5}
    K: List[float] = []
    for key in _POSITION_ORDER:
        x, y, z = _position_vector(positions, key, idx_of[key], holder)
        if nrm_mean is not None:
            x, y, z = x - nrm_mean[0], y - nrm_mean[1], z - nrm_mean[2]
        K.extend([x, y, z])
    K_arr = np.array(K, dtype="f")

    try:
        params = ipmag.calculate_aniso_parameters(K_arr, n_pos=n_pos)
    except UnboundLocalError:
        # bug pmagpy confirme (voir docstring module) : residu EXACTEMENT
        # nul (mesures parfaitement conformes au tenseur ajuste, sans le
        # moindre bruit - essentiellement jamais atteint sur de vraies
        # mesures instrumentales), `sigma` n'est jamais initialise avant
        # le bloc F-test dans ce cas et `calculate_aniso_parameters`
        # plante en y accedant. Reconstruit ici uniquement s/v1/v2/v3/p
        # (copie allegee de la premiere moitie, non affectee du bug) au
        # lieu de laisser planter l'appelant - pas de F-test/sigma dans
        # ce cas (aucun residu a en tirer de toute facon).
        matrices = ipmag.get_matrix(n_pos)
        s_bs = np.dot(matrices["B"], K_arr)
        trace = s_bs[0] + s_bs[1] + s_bs[2]
        s_bs = s_bs / trace
        s1, s2, s3, s4, s5, s6 = (float(v) for v in s_bs)
        from numpy.linalg import eig
        import pmagpy.pmag as pmag
        s_matrix = [[s1, s4, s6], [s4, s2, s5], [s6, s5, s3]]
        # .real : numpy.linalg.eig (le solveur GENERAL, pas eigh - meme
        # choix que pmagpy) peut retourner un dtype complexe (partie
        # imaginaire negligeable/nulle) sur une matrice parfaitement
        # conditionnee comme ici (residu exactement nul) - les valeurs
        # propres d'une matrice symetrique REELLE sont mathematiquement
        # toujours reelles, .real evite juste un artefact d'affichage
        # ("0.397+0.000j") sans rien changer numeriquement.
        t, evectors = eig(s_matrix)
        t = [float(v.real) for v in t]
        evectors = evectors.real
        t1, t3 = max(t), min(t)
        ix1, ix3 = t.index(t1), t.index(t3)
        ix2 = next(i for i in range(3) if i not in (ix1, ix3))
        t2 = t[ix2]
        v1 = pmag.cart2dir([evectors[0][ix1], evectors[1][ix1], evectors[2][ix1]])
        v2 = pmag.cart2dir([evectors[0][ix2], evectors[1][ix2], evectors[2][ix2]])
        v3 = pmag.cart2dir([evectors[0][ix3], evectors[1][ix3], evectors[2][ix3]])
        return PmagpyAnisotropyResult(
            n_pos=n_pos, s=(s1, s2, s3, s4, s5, s6),
            v1=(t1, v1[0], v1[1]), v2=(t2, v2[0], v2[1]), v3=(t3, v3[0], v3[1]),
            p=t1 / t3 if t3 else 0.0,
            sigma=None, f_test=None, f12_test=None, f23_test=None,
            f_crit=None, quality=None,
        )

    s = _parse_triplet_6(params["aniso_s"])
    v1 = _parse_triplet(params["aniso_v1"])
    v2 = _parse_triplet(params["aniso_v2"])
    v3 = _parse_triplet(params["aniso_v3"])
    return PmagpyAnisotropyResult(
        n_pos=n_pos, s=s, v1=v1, v2=v2, v3=v3,
        p=float(params["aniso_p"]),
        sigma=float(params["aniso_s_sigma"]) if "aniso_s_sigma" in params else None,
        f_test=float(params["aniso_ftest"]) if "aniso_ftest" in params else None,
        f12_test=float(params["aniso_ftest12"]) if "aniso_ftest12" in params else None,
        f23_test=float(params["aniso_ftest23"]) if "aniso_ftest23" in params else None,
        f_crit=_parse_f_crit(params.get("description")),
        quality=params.get("aniso_ftest_quality"),
    )


def _parse_triplet_6(text: str) -> Tuple[float, float, float, float, float, float]:
    a, b, c, d, e, f = (float(x) for x in text.split(":"))
    return a, b, c, d, e, f


def _parse_f_crit(description: Optional[str]) -> Optional[float]:
    if not description or "Critical F:" not in description:
        return None
    try:
        return float(description.split("Critical F:")[1].strip())
    except ValueError:
        return None


def format_aarm_pmagpy(id_: str, r: PmagpyAnisotropyResult) -> str:
    lines = [
        f"{id_}  PmagPy AARM/ATRM (n_pos={r.n_pos}, least-squares design matrix + Hext F-test)",
        f"  s (trace-normalized, s1+s2+s3=1): "
        f"s1={r.s[0]:.5f} s2={r.s[1]:.5f} s3={r.s[2]:.5f} "
        f"s4(s12)={r.s[3]:.5f} s5(s23)={r.s[4]:.5f} s6(s13)={r.s[5]:.5f}",
        f"  principal axes: "
        f"v1(max)=({r.v1[0]:.4f}, dec={r.v1[1]:.1f}, inc={r.v1[2]:.1f})  "
        f"v2(int)=({r.v2[0]:.4f}, dec={r.v2[1]:.1f}, inc={r.v2[2]:.1f})  "
        f"v3(min)=({r.v3[0]:.4f}, dec={r.v3[1]:.1f}, inc={r.v3[2]:.1f})",
        f"  anisotropy degree P (t1/t3): {r.p:.4f}",
    ]
    if r.sigma is not None and r.f_test is not None:
        lines.append(
            f"  Hext F-test: sigma={r.sigma:.3e}  F={r.f_test:.2f}  F12={r.f12_test:.2f}  "
            f"F23={r.f23_test:.2f}"
            + (f"  (critical F={r.f_crit:.4f})" if r.f_crit is not None else "")
            + (f"  -> {'significant anisotropy' if r.quality == 'g' else 'NOT significant'}"
               if r.quality else "")
        )
    else:
        lines.append(
            "  Hext F-test: not computable (pmagpy quirk - exactly-zero "
            "residual, essentially never happens on real noisy data)"
        )
    return "\n".join(lines)
