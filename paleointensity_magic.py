"""
Second traitement, PARALLELE et INDEPENDANT, des donnees de paleointensite
IZZI/Thellier : appelle directement le code de PmagPy/MagIC (deja installe,
`pmagpy` 4.5.2), plutot que de re-deriver ses formules - demande explicite
utilisateur ("I am interested in having a second parallel processing of
paleointensity. the one from Magic. Is it possible?"). Les fonctions
appelees (`pmag.sortarai`, `pmag.find_dmag_rec`, `pmag.PintPars`) sont
CELLES utilisees par le Thellier GUI officiel pour calculer les
statistiques standard (Paterson et al. 2014) : b (pente Arai), sigma, MAD,
DANG, FRAC, gap_max, SCAT, f, g, q, DRATS, MD, ZigZag - independamment du
calcul natif STARpaleomag_Py (`paleointensity.py`, porte du Fortran).

Adaptateur (`build_magic_dataframe`) : construit le DataFrame MagIC data
model 3 attendu par ces fonctions a partir d'un specimen STARpaleomag_Py. Reutilise
`magic_export._measurement_treatment` (deja valide, utilise par
convert_ren_to_r.py/l'export MagIC) pour method_codes/treat_temp(K)/
treat_dc_field(T)/phi/theta - PAS reimplemente ici.

Correction reutilisee ici (bug reel repere dans magic_export.py en
construisant ce module, PAS dans le Fortran d'origine - propre a ce
portage, voir magic_export._izzi_order_tags pour le detail complet et sa
verification sur un vrai fichier IZZI) : `_measurement_treatment` tague
STATIQUEMENT cod1='S' avec "LP-PI-TRM-IZ" et cod1='R' avec
"LP-PI-TRM-ZI", quel que soit l'ordre reel de mesure a cette
temperature - corrige, desormais a la source (magic_export.py, partagee
avec build_measurements_rows/l'export MagIC), en redecouvrant l'ordre
reel depuis la sequence.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import pandas as pd

from magic_export import _measurement_treatment, _izzi_order_tags
from selection import polere
from testlect import Pmag
from convert_ren_to_r import _parse_ifield_uT

from pmagpy import pmag

# Seuil standard MagIC (Paterson et al. 2014, criteres "Class B") pour la
# "scat box" - sans lui, PintPars ne calcule pas SCAT (retourne toujours
# "t"/pass par defaut, voir pmag.PintPars ~ligne 3666).
DEFAULT_B_BETA_THRESHOLD = 0.10


@dataclass
class MagicPintResult:
    """Resultat du traitement MagIC/PmagPy - miroir de ce que retourne
    `pmag.PintPars`, plus la paleointensite deja convertie en uT (PintPars
    ne la calcule pas lui-meme : `int_b` est un RAPPORT NRM/TRM, a
    multiplier par le champ labo)."""
    specimen: str
    n: int
    step_first: float
    step_last: float
    b: float                 # pente Arai (int_b, sans dimension, negative)
    b_sigma: float
    b_beta: float
    field_lab_uT: float
    paleointensity_uT: float
    mad: float                # int_mad_free
    dang: float
    f: float
    fvds: float
    g: float
    q: float
    n_ptrm: int
    drat: Optional[float]   # None si non fiable, voir _drat_unreliable_pmagpy_bug
    drats: Optional[float]
    md: float
    frac: float
    gap_max: float
    scat: str                 # "t" (pass) / "f" (fail)
    zigzag: float              # int_z : -1 = non teste, sinon Frat/Trat (>1 = echec)
    zigzag_method: str
    raw: Dict = field(default_factory=dict, repr=False)  # dict complet retourne par PintPars


# Pas pertinents pour la paleointensite - demande explicite utilisateur
# ("perhaps best to remove all steps except the N0 or D0, R, V, P, remove
# the X, Y, Z, Q, L perhaps") : l'ATRM (X/Y/Z) et la vitesse de
# refroidissement (L/Q) n'apportent rien au calcul Arai/PintPars et n'ont
# pas besoin d'etre exposees a PmagPy - ecartees ici plutot que de compter
# sur les filtres method_codes de sortarai/find_dmag_rec pour les ignorer
# correctement.
_PALEOINT_COD1 = {"N", "D", "R", "V", "S", "P"}


def _row_from_measurement(ech: Pmag, m, prev, ifield: float, codes_override=None, xyz_override=None) -> Dict:
    codes, temp, af_field, dc_field, phi, theta = _measurement_treatment(m, prev, ifield)
    if codes_override is not None:
        codes = codes_override
    x, y, z = xyz_override if xyz_override is not None else (m.x, m.y, m.z)
    mag, dec, inc = polere(x, y, z)
    return {
        "specimen": ech.id,
        "treat_temp": temp,
        "treat_dc_field": dc_field,
        "treat_dc_field_phi": phi,
        "treat_dc_field_theta": theta,
        "dir_dec": dec,
        "dir_inc": inc,
        "magn_moment": mag,
        "method_codes": codes,
        "treat_ac_field": af_field,
        "instrument_codes": "",
    }


def build_magic_dataframe(ech: Pmag) -> pd.DataFrame:
    """DataFrame MagIC data model 3 (colonnes : specimen, treat_temp (K),
    treat_dc_field (T), treat_dc_field_phi/theta, dir_dec, dir_inc,
    magn_moment (Am2), method_codes, treat_ac_field, instrument_codes) -
    une ligne par mesure utile (voir _PALEOINT_COD1), MEME ORDRE que dans
    `ech.mesures` (l'ordre reel importe pour le flag IZ/ZI en IZZI, et
    pour l'adjacence R/V utilisee par le controle pTRM 'P' via `prev`).

    Methode Thellier "champ en Z puis Z-" (R+V sans 'S', confirme par
    l'utilisateur sur des donnees reelles, specimen "02B" : "field in Z
    and Z-, there is no choice than (R+V)/2... I think that magic was
    able to recognize this feature") : VERIFIE que ce n'est PAS le cas -
    `pmag.sortarai` (donnees thermiques) n'a aucune notion de "champ
    oppose" (ce mecanisme, "perpendicular method"/ThetaChecks-DeltaChecks
    dans pmag.PintPars, n'existe QUE pour les donnees micro-ondes via
    `pmag.sortmwarai`, pas pour sortarai). R et V sont donc TOUS LES DEUX
    reellement mesures EN CHAMP (voir magic_export._measurement_treatment,
    theta=+90/-90) - il n'existe pas de pas zero-field distinct. Reconstruit
    ici la meme moyenne que le natif STARpaleomag_Py (paleointensity.py:464,
    `nrm_vec = (R+V)/2`) : le contenu (dec/inc/moment) du pas R est
    REMPLACE par la moyenne cartesienne (R+V)/2 et tague "LT-T-Z" (zero-
    field synthetique), V restant INCHANGE ("LT-T-I", en-champ) - PmagPy
    calcule alors pTRM = |V - R_moyenne| = |R-V|/2, la meme magnitude que
    la vraie difference physique (le signe s'inverse mais sortarai ne
    garde que la norme via cart2dir). Rien n'est modifie dans le fichier
    .prmag lui-meme (magic_export.py) : R y reste correctement decrit
    comme reellement mesure en champ - cette reconstruction est purement
    une transformation d'ANALYSE, propre a ce module."""
    ifield = next((m.treat_dc_field for m in ech.mesures if m.treat_dc_field), None)
    if ifield is None:
        ifield = _parse_ifield_uT(ech.com) or 0.0
    is_izzi = any(m.cod1 == "S" for m in ech.mesures)

    mesures = [m for m in ech.mesures if m.cod1 in _PALEOINT_COD1]
    izzi_tags = _izzi_order_tags(mesures)

    rows = []
    for idx, m in enumerate(mesures):
        prev = mesures[:idx]
        if not is_izzi and m.cod1 == "R":
            partner = next((mm for mm in mesures if mm.cod1 == "V" and mm.etape == m.etape), None)
            if partner is not None:
                xyz_avg = ((m.x + partner.x) / 2.0, (m.y + partner.y) / 2.0, (m.z + partner.z) / 2.0)
                rows.append(_row_from_measurement(
                    ech, m, prev, ifield,
                    codes_override="LT-T-Z:LP-PI-TRM:LP-PI-II:LP-PI-ALT",
                    xyz_override=xyz_avg,
                ))
                continue
        row = _row_from_measurement(ech, m, prev, ifield)
        if idx in izzi_tags:
            parts = [p for p in row["method_codes"].split(":") if p not in ("LP-PI-TRM-IZ", "LP-PI-TRM-ZI")]
            parts.append(izzi_tags[idx])
            row["method_codes"] = ":".join(parts)
        rows.append(row)
    return pd.DataFrame(rows)


def compute_magic_paleointensity(
    ech: Pmag,
    step_first: Optional[float] = None,
    step_last: Optional[float] = None,
    b_beta_threshold: float = DEFAULT_B_BETA_THRESHOLD,
) -> MagicPintResult:
    """Lance `pmag.sortarai` + `pmag.find_dmag_rec` + `pmag.PintPars` (le
    code PmagPy/MagIC reel, pas une reimplementation) sur ce specimen.
    `step_first`/`step_last` (degC, meme convention que le fit natif
    STARpaleomag_Py - None = toute la sequence disponible) selectionnent
    l'intervalle NRM utilise pour la pente/MAD/DANG/etc, convertis ici en
    indices dans `araiblock[0]` (PintPars les attend en indices, pas en
    temperature). Leve ValueError si les donnees ne ressemblent pas a une
    experience IZZI/Thellier exploitable (pas assez de pas, etc) - la
    faute a PmagPy elle-meme via son `errcode`/ses exceptions, pas
    masquee ici."""
    is_izzi = any(m.cod1 == "S" for m in ech.mesures)
    df = build_magic_dataframe(ech)
    araiblock, field_T = pmag.sortarai(df, ech.id, False, version=3)
    zijdblock, _units = pmag.find_dmag_rec(ech.id, df, version=3)

    first_Z = araiblock[0]
    if len(first_Z) < 3:
        raise ValueError(
            f"{ech.id}: not enough IZZI/Thellier steps found ({len(first_Z)}) "
            "- check method codes / cod1 R/V/S/P coverage.")

    temps = [rec[0] for rec in first_Z]
    start = 0
    if step_first is not None:
        target = step_first + 273
        start = min(range(len(temps)), key=lambda k: abs(temps[k] - target))
    end = 0
    if step_last is not None:
        target = step_last + 273
        end = min(range(len(temps)), key=lambda k: abs(temps[k] - target))

    accept = {"int_b_beta": str(b_beta_threshold)} if b_beta_threshold else {}
    try:
        pars, errcode = pmag.PintPars(df, araiblock, zijdblock, start, end, accept, version=3)
    except NameError as e:
        # Bug reel dans pmagpy 4.5.2 (derniere version publiee au moment
        # ou ceci a ete ecrit) : pmag.PintPars reference une variable
        # locale `b_key` jamais definie dans SA PROPRE portee (seulement
        # dans int_pars, une fonction imbriquee distincte) des qu'un seuil
        # int_b_beta est fourni dans `accept` - voir pmag.py ~ligne 3668
        # ("b = pars[b_key]"). PAS une erreur de cet adaptateur : on
        # retente simplement sans criteres d'acceptation (SCAT ne sera
        # alors pas calcule, cf. le commentaire pmagpy juste avant ce
        # bloc : "if threshold value for beta is not defined, then scat
        # cannot be calculated (pass)") plutot que de faire echouer tout
        # le calcul pour un bug tiers.
        if "b_key" not in str(e) or not accept:
            raise
        pars, errcode = pmag.PintPars(df, araiblock, zijdblock, start, end, {}, version=3)
        pars["_scat_skipped_pmagpy_bug"] = True
    if errcode == 1:
        raise ValueError(f"{ech.id}: PmagPy PintPars failed (errcode=1, likely too few points).")

    # Bug reel dans pmagpy 4.5.2 (verifie sur magic_contribution_19987.prmag,
    # ~moitie des 172 specimens reels touches) : le calcul de DRAT/DRATS
    # dans PintPars (`for irec in first_I: if irec[0]==step: break`) ne
    # gere pas le cas ou le pas cible d'un controle pTRM n'a JAMAIS ete
    # mesure comme pas normal (frequent : le premier controle pTRM revient
    # souvent a une temperature plus basse que le tout premier pas reel de
    # la sequence). Sans `break`, `irec` reste sur le DERNIER element de
    # first_I (une temperature elevee sans rapport), corrompant
    # silencieusement diffcum/drat_max - confirme numeriquement : DRAT
    # saute a ~95-99% (au lieu de quelques %) exactement quand un pas de
    # controle pTRM est absent de first_I. Detecte ici et signale plutot
    # que de presenter une valeur fausse comme fiable.
    first_I_temps = {rec[0] for rec in araiblock[1]}
    ptrm_check_temps = {rec[0] for rec in araiblock[2]}
    drat_reliable = ptrm_check_temps <= first_I_temps

    # 2e raison, distincte, de se mefier de DRAT/DRATS pour un specimen
    # Thellier "champ oppose" (R+V, voir build_magic_dataframe) : le
    # controle pTRM de `sortarai` soustrait TOUJOURS l'enregistrement
    # IMMEDIATEMENT PRECEDENT (`brec = datablock[step-1]`, sans egard a
    # son type), en supposant implicitement que ce voisin est LE PAS
    # ZERO-FIELD de reference. Ici, la sequence reelle (confirmee sur
    # "02B") place le controle pTRM APRES la paire R,V COMPLETE de la
    # temperature courante - le voisin immediat est donc le pas EN CHAMP
    # (V), pas le zero-field synthetique - la soustraction melange alors
    # deux etats en champ au lieu d'isoler le pTRM reacquis. Constate
    # numeriquement : magnitudes de controle ~2x superieures a
    # l'acquisition d'origine au meme palier, avec par ailleurs un
    # ajustement par ailleurs excellent (MAD/DANG/b_beta/q tous bons) -
    # signe d'un probleme localise au controle pTRM, pas au fit lui-meme.
    if not is_izzi:
        drat_reliable = False
        pars["_drat_unreliable_opposite_field_method"] = True
    if not drat_reliable:
        pars["_drat_unreliable_pmagpy_bug"] = True

    b, sigma = pars["int_b"], pars["int_b_sigma"]
    field_uT = field_T * 1.0e6
    paleointensity_uT = abs(b) * field_uT

    return MagicPintResult(
        specimen=ech.id,
        n=pars["int_n_measurements"],
        step_first=pars["meas_step_min"] - 273,
        step_last=pars["meas_step_max"] - 273,
        b=b, b_sigma=sigma, b_beta=pars["int_b_beta"],
        field_lab_uT=field_uT, paleointensity_uT=paleointensity_uT,
        mad=pars["int_mad_free"], dang=pars["int_dang"],
        f=pars["int_f"], fvds=pars["int_fvds"], g=pars["int_g"], q=pars["int_q"],
        n_ptrm=pars["int_n_ptrm"],
        drat=None if not drat_reliable else pars["int_drat"],
        drats=None if not drat_reliable else pars["int_drats"],
        md=pars["int_md"], frac=pars["int_frac"], gap_max=pars["int_gmax"],
        scat=pars["int_scat"], zigzag=pars["int_z"], zigzag_method=pars["method_codes"],
        raw=pars,
    )


def format_magic_paleointensity(r: MagicPintResult) -> str:
    """Rendu texte (meme esprit que paleointensity.py : lisible dans la
    console de l'appli, pas un widget dedie)."""
    zigzag_txt = "not tested" if r.zigzag == -1 else f"{r.zigzag:.2f} ({r.zigzag_method or 'below threshold'})"
    scat_txt = "pass" if r.scat == "t" else "FAIL"
    if r.raw.get("_scat_skipped_pmagpy_bug"):
        scat_txt = "not computed (pmagpy 4.5.2 bug, see paleointensity_magic.py)"
    return (
        f"MagIC/PmagPy paleointensity (pmag.PintPars) - {r.specimen}\n"
        f"  steps used     : {r.step_first:.0f} - {r.step_last:.0f} degC  (n={r.n}, n_ptrm={r.n_ptrm})\n"
        f"  lab field      : {r.field_lab_uT:.2f} uT\n"
        f"  paleointensity : {r.paleointensity_uT:.2f} uT   (b={r.b:.4f}, b_sigma={r.b_sigma:.4f}, b_beta={r.b_beta:.4f})\n"
        f"  MAD (free)     : {r.mad:.2f}    DANG: {r.dang:.2f}\n"
        f"  f: {r.f:.3f}   f_vds: {r.fvds:.3f}   g: {r.g:.3f}   q: {r.q:.2f}\n"
        f"  FRAC: {r.frac:.3f}   gap_max: {r.gap_max:.3f}   SCAT: {scat_txt}"
        f"  (b_beta threshold used: {DEFAULT_B_BETA_THRESHOLD})\n"
        f"  DRAT: {_pct(r.drat, r.raw)}   DRATS: {_pct(r.drats, r.raw)}   MD(tail): {r.md:.2f}%\n"
        f"  ZigZag (Shaar & Tauxe 2013): {zigzag_txt}\n"
    )


def _pct(value: Optional[float], raw: Dict) -> str:
    if value is not None:
        return f"{value:.2f}%"
    if raw.get("_drat_unreliable_opposite_field_method"):
        return "n/a (opposite-field Thellier: pTRM check baseline is the in-field step, not zero-field - see paleointensity_magic.py)"
    return "n/a (pmagpy 4.5.2 bug: pTRM check step absent from first_I, see paleointensity_magic.py)"
