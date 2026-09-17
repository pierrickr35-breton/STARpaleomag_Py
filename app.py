import io
import os
import re
import subprocess
import sys
import tkinter as tk
import webbrowser
from tkinter import ttk, filedialog, messagebox

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from testlect import read_ren_file, read_prmag_file
from convert_legacy_ren import convert_legacy_auto
from complete_sample_info import complete_site_info, complete_specimen_info, complete_sample_height
from import_new_data import parse_jr6_file, parse_legacy_new_measurements, archive_new_measurements
from convert_ren_to_r import convert_file as convert_ren_to_r_file
from convert_magic_to_r import convert_magic_file, scan_dip_sign
from convert_utrecht_to_r import convert_files as convert_utrecht_files
from convert_ipgp_to_r import convert_files as convert_ipgp_files
from convert_montpellier_to_r import convert_file as convert_montpellier_file, DEFAULT_VOLUME_CM3 as _MONTPELLIER_DEFAULT_VOLUME_CM3
from irm import build_irm_figure, has_irm_data
from selection import (
    select_samples,
    select_samples_by_site,
    delete_measurements,
    init_selection,
    list_measurements,
    list_xyz,
    list_measurements_vrm,
    list_measurements_depth,
    sample_info,
    _build_selected_sample,
)
from interpretation_quality import evaluate_result, evaluate_results, format_quality_report
from auto_interpretation import propose_components, propose_components_for_site, format_suggestions
from calcul import (
    FitResult,
    fit_line,
    fit_plane,
    fit_fisher_direction,
    fit_single_direction,
    fit_lines_auto,
    fit_from_redo_file,
    fisher_from_measurements,
    fisher_from_results,
    build_site_mean_result,
    site_lat_lon_from_donnees,
    site_of_result,
    import_published_means,
    _LINE_LIKE_CAT1,
    list_results,
    init_results,
    results_path_for,
    ani_path_for,
    pmagint_path_for,
    _iter_result_lines,
    find_duplicate_measurements,
    find_duplicate_results,
    find_duplicate_pmagint_rows,
    find_duplicate_ani_tensors,
    archivres,
    load_results,
    recompute_fit_geometry,
    available_mean_orientations,
    list_mdf,
    compute_mean_intensity,
    compute_mean_susceptibility,
    format_mean_intensity,
    compute_mean_inclination,
    format_mean_inclination,
    compute_koenigsberger,
    format_koenigsberger,
    list_diff_measurements,
    apply_viscosity_test,
    apply_subtraction,
    record_arm_holder,
    detect_cooling_rate_rows,
    compute_cooling_rate,
    format_cooling_rate,
    read_ani_tensor,
    read_ani_mean_tensor,
    read_all_ani_tensors,
    find_orphan_ani_specimens,
    apply_inverse_anisotropy,
    compute_anisotropy_tensor,
    replace_position_by_symmetry,
    detect_six_positions,
    write_ani_tensors,
    compute_anicor_factor,
    _ANI_CODE2,
    _correct_dec_inc,
    mean_components,
)
from zijderveld import build_zijderveld_figure, draw_zijderveld
from stereo import build_stereo_figure, build_stereo_results_figure
from xygraph import build_xygraph_figure, has_mixed_demag
from susceptibility import build_susceptibility_figure
from svgwriter import SVGWriter
from paleointensity import (
    compute_arai,
    fit_arai_line,
    fit_arai_direction,
    fit_arai_direction_corrected,
    angle_between_vectors,
    nrm0_vector,
    compute_rf1_rf2,
    parse_com_field,
    detect_method_and_hlab,
    compute_crm,
    arai_curvature,
    build_arai_figure,
    build_paleoint_review_figure,
    draw_arai,
    write_pmagint_line,
    read_pmagint,
)
from paleointensity_magic import compute_magic_paleointensity, format_magic_paleointensity
from datatools import (
    convert_thellier_to_nrm,
    remove_step,
    remove_bad_quality_steps,
    eliminate_grm,
    convert_z_minus,
    export_thellier_tdt,
    detect_grm,
)
from magic_export import (
    export_to_magic,
    classify_anisotropy_experiment,
    load_site_metadata_table,
    build_site_metadata_preview_rows,
    write_site_metadata_preview,
)
from anisotropy_magic import compute_aarm_pmagpy, format_aarm_pmagpy
from export_stereo import export_stereo_project
from detailed_export import export_detailed_txt
from field_notes import (
    parse_orientation_file,
    parse_complement_file,
    apply_complement,
    check_utdif_for_sites,
    build_pmag_records,
    write_prmag_from_field_notes,
    write_ged_file,
    write_diagnostics_report,
    parse_ged_file,
)

SVG_DEBUG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "svg_debug")

# Raccourcis clavier repris de starmac_OSX.f95 (AWE_addMenu(...,"CTRL+...")),
# adaptes suite a un diagnostic reel (diag_clavier.py) sur clavier AZERTY :
# Option (ALT) est le modificateur d'ACCENT sur les claviers europeens -
# Cmd+Option+Shift+E genere le keysym "Ecircumflex" (Ê), jamais "E". Tout
# raccourci combinant Option + une lettre est donc casse par construction
# sur AZERTY (confirme, pas specifique aux chiffres comme suppose au debut).
# On evite donc completement Option : CTRL+ALT+x (Fortran) est traduit ici
# en Cmd+Ctrl+<lettre> plutot qu'en Cmd+Option+Shift+<lettre>, en choisissant
# des lettres qui evitent les raccourcis systeme macOS connus (Cmd+Ctrl+Q =
# verrouiller l'ecran, Cmd+Ctrl+F = plein ecran).
#   - "CTRL+l" (lismes)          -> Cmd+L        (CTRL -> Command)
#   - "CTRL+META+z" (plotzijder) -> Cmd+Ctrl+Z    (META -> Control)
#   - "CTRL+ALT+1" (selce)       -> Cmd+Ctrl+A    (ALT -> Control aussi, pas Option)
# Note de syntaxe Tk (piege verifie empiriquement) : des que Shift est un
# modificateur explicite, Tk exige la lettre du keysym en MAJUSCULE pour que
# le binding se declenche - non applicable ici puisqu'on n'utilise plus Shift
# du tout dans les combinaisons a 3 modificateurs (sauf sur Windows, voir
# _SHORTCUTS_WIN plus bas).
# valeur : (accelerateur affiche dans le menu, sequence de bind Tk)
_SHORTCUTS_MAC = {
    "importpc":   ("Cmd+O",      "<Command-o>"),
    "starend":    ("Cmd+Ctrl+X", "<Command-Control-x>"),
    "selmes":     ("Cmd+E",      "<Command-e>"),
    "selentete":  ("Cmd+Ctrl+E", "<Command-Control-e>"),
    "effmes":     ("Cmd+D",      "<Command-d>"),
    "initmes":    ("Cmd+I",      "<Command-i>"),
    "lismes":     ("Cmd+L",      "<Command-l>"),
    "infoech":    ("Cmd+J",      "<Command-j>"),
    "selce":      ("Cmd+Ctrl+A", "<Command-Control-a>"),
    "selis":      ("Cmd+Ctrl+B", "<Command-Control-b>"),
    "selcp":      ("Cmd+Ctrl+T", "<Command-Control-t>"),
    "ajuslig":    ("Cmd+B",      "<Command-b>"),
    "fishmes":    ("Cmd+F",      "<Command-f>"),
    "fishres":    ("Cmd+Ctrl+F", "<Command-Control-f>"),
    "lisres":     ("Cmd+Ctrl+L", "<Command-Control-l>"),
    "selres":     ("Cmd+Ctrl+R", "<Command-Control-r>"),
    "initres":    ("Cmd+Ctrl+I", "<Command-Control-i>"),
    "plotzijder": ("Cmd+Ctrl+Z", "<Command-Control-z>"),
    "xygraph":    ("Cmd+Ctrl+Y", "<Command-Control-y>"),
    "stereodata": ("Cmd+Ctrl+S", "<Command-Control-s>"),
    # Cmd+Ctrl+Option+S essaye d'abord (demande explicite utilisateur) -
    # BUG CONFIRME par l'utilisateur ("the shortcuts do not work for
    # stereo results") : le gotcha Option+lettre documente en tete de ce
    # dict (caractere accentue au lieu du keysym attendu) s'applique donc
    # AUSSI avec Control comme modificateur additionnel, pas seulement
    # avec Shift comme initialement verifie. Remplace par Cmd+Ctrl+P (pas
    # d'Option) - meme lettre que _SHORTCUTS_WIN ("stereores": Ctrl+Shift+P),
    # pour rester coherent entre les deux plateformes.
    "stereores":  ("Cmd+Ctrl+P", "<Command-Control-p>"),
    "ajusfisher": ("Cmd+Ctrl+G", "<Command-Control-g>"),
    "suscep":     ("Cmd+Ctrl+K", "<Command-Control-k>"),
    "meanint":    ("Cmd+Ctrl+M", "<Command-Control-m>"),
    "anisotropy": ("Cmd+Ctrl+N", "<Command-Control-n>"),
}

# Version PC/Windows (et Linux) - demande explicite utilisateur ("for a PC
# version can you provide shortcuts... those previously defined often
# failed"). Cause racine : "Command"/"Meta" ne correspondent a AUCUNE touche
# physique fiable sous Tk sur Windows (contrairement a "Cmd" sur Aqua) - le
# bind ne leve pas d'erreur, il ne se declenche simplement JAMAIS, quel que
# soit ce qu'on appelle la touche "Meta". Reprendre le meme schema en
# substituant juste Command->Control (Ctrl+lettre seul) NE MARCHE PAS non
# plus : ca percute (a) les bindings Emacs internes de Tk sur les widgets
# Entry/Text (Ctrl+A/E/D/B/F/N/P/K/Y deja lies - meme raison que le choix de
# Command sur Mac, voir _setup_shortcuts) et (b) les raccourcis universels
# Windows/applicatifs (Ctrl+S=Enregistrer, Ctrl+F=Rechercher, Ctrl+Z=Annuler,
# Ctrl+Y=Retablir, Ctrl+A=Tout selectionner...) - exactement les exemples
# signales cassés par l'utilisateur (Zijderveld/GraphXY/Stereo data/best dir
# Fisher utilisaient F/S/Y/Z/A en 2e lettre). Ctrl+Alt+lettre est EVITE : sur
# les claviers europeens Windows (AZERTY/QWERTZ...), Ctrl+Alt EQUIVAUT a
# AltGr, le meme mecanisme de composition de caracteres accentues/speciaux
# que le gotcha Option documente ci-dessus pour Mac. Ctrl+Shift+<LETTRE> est
# donc retenu partout : aucune collision connue (Windows/apps/Tk), et Shift
# ne compose jamais de caractere special. Un seul palier de raccourcis ici
# (pas de Ctrl+lettre "simple" + Ctrl+Shift+lettre "avance" comme sur Mac) -
# chaque commande a sa propre lettre, jamais reutilisee.
_SHORTCUTS_WIN = {
    "importpc":   ("Ctrl+Shift+O", "<Control-Shift-O>"),
    "starend":    ("Ctrl+Shift+Q", "<Control-Shift-Q>"),
    "selmes":     ("Ctrl+Shift+E", "<Control-Shift-E>"),
    "selentete":  ("Ctrl+Shift+H", "<Control-Shift-H>"),
    "effmes":     ("Ctrl+Shift+D", "<Control-Shift-D>"),
    "initmes":    ("Ctrl+Shift+I", "<Control-Shift-I>"),
    "lismes":     ("Ctrl+Shift+L", "<Control-Shift-L>"),
    "infoech":    ("Ctrl+Shift+J", "<Control-Shift-J>"),
    "selce":      ("Ctrl+Shift+A", "<Control-Shift-A>"),
    "selis":      ("Ctrl+Shift+B", "<Control-Shift-B>"),
    "selcp":      ("Ctrl+Shift+T", "<Control-Shift-T>"),
    "ajuslig":    ("Ctrl+Shift+G", "<Control-Shift-G>"),
    "fishmes":    ("Ctrl+Shift+F", "<Control-Shift-F>"),
    "fishres":    ("Ctrl+Shift+Y", "<Control-Shift-Y>"),
    "lisres":     ("Ctrl+Shift+W", "<Control-Shift-W>"),
    "selres":     ("Ctrl+Shift+R", "<Control-Shift-R>"),
    "initres":    ("Ctrl+Shift+N", "<Control-Shift-N>"),
    "plotzijder": ("Ctrl+Shift+Z", "<Control-Shift-Z>"),
    "xygraph":    ("Ctrl+Shift+X", "<Control-Shift-X>"),
    "stereodata": ("Ctrl+Shift+S", "<Control-Shift-S>"),
    "stereores":  ("Ctrl+Shift+P", "<Control-Shift-P>"),
    "ajusfisher": ("Ctrl+Shift+C", "<Control-Shift-C>"),
    "suscep":     ("Ctrl+Shift+K", "<Control-Shift-K>"),
    "meanint":    ("Ctrl+Shift+M", "<Control-Shift-M>"),
    "anisotropy": ("Ctrl+Shift+V", "<Control-Shift-V>"),
}

SHORTCUTS = _SHORTCUTS_WIN if sys.platform.startswith("win") else _SHORTCUTS_MAC

# Marque une ligne de TITRE DE COLONNES ("header" au sens tableau, pas
# "en-tete de fichier") a l'interieur d'un texte de rapport, pour que
# STARpaleomagApp._afficher l'affiche en gras - demande explicite utilisateur
# ("throughout the software, is it possible to write the header in bold,
# the output text might be easier to read"). Un caractere de controle
# improbable dans un texte normal (pas de conflit avec les donnees
# reelles) plutot qu'un marquage par position/regex fragile - chaque
# fonction de rapport (selection.list_*, calcul.list_results, etc.)
# entoure sa ligne d'en-tete de colonnes avec HEADER_MARK ; _afficher
# retire les marqueurs et applique le tag "header" au texte entre eux.
HEADER_MARK = "\x01"

# Meme principe que HEADER_MARK, pour un avertissement FORT (tag "warn" -
# gras, rouge, fond surligne) - demande explicite utilisateur ("write a
# strong warning") pour signaler un signe de dip suspect a l'import MagIC
# (voir convert_magic_to_r.scan_dip_sign / ouvrir_convert_magic_to_r_dialog).
WARN_MARK = "\x02"

ORIENTATION_SHORTCUT_NAMES = {1: "selce", 2: "selis", 3: "selcp"}

ORIENTATIONS = {
    "Sample (SC)": 1,
    "In situ (IS)": 2,
    "Tilt cor. (TC)": 3,
}


def _step_token(m) -> str:
    """"<etape><cod1><cod2>" compact, ex. "250D+" ou "15.5F+" - pour les
    listes de pas disponibles montrees dans les dialogues interactifs
    (pas le rapport complet "List data", voir selection._fmt_step pour
    celui-la). `etape` est un float depuis "convert all step integer to
    float" (demande explicite utilisateur) : affiche sans decimale quand
    la valeur est entiere (le cas thermique, ex. "250" pas "250.0") et
    avec une decimale sinon (le cas AF fractionnaire, ex. "15.5")."""
    v = m.etape
    num = f"{v:.0f}" if v == int(v) else f"{v:.1f}"
    return f"{num}{m.cod1}{m.cod2}"


def _format_pmagint_row(sample_id: str, row: dict) -> str:
    """Formate une ligne DEJA archivee dans .pmagint (voir
    paleointensity.read_pmagint/_PMAGINT_HEADER pour les cles) pour
    l'affichage en lecture seule (ouvrir_view_pmagint_dialog) - reprend
    la STRUCTURE du texte affiche par ouvrir_openfilepint_dialog (memes
    lignes f/g/q/ccr/H, Hlab/b/sb, fcor/Hcorani, gamma/k, mad/dang/fvds/
    frac/gap_max/n_ptrm, rf1/rf2) mais avec les valeurs STOCKEES telles
    quelles - "n.d" affiche tel quel, jamais recalcule ni reformate en
    float (une valeur "n.d" ne doit pas planter sur un format `:.1f`)."""
    def g(key: str) -> str:
        return row.get(key, "n.d")

    out = io.StringIO()
    out.write(f"\n{sample_id} : Tmin={g('t1')} Tmax={g('t2')}  N={g('N')}\n")
    out.write(
        f" f={g('f')}  g={g('g')}  q={g('q')}  ccr={g('ccr')}  h={g('H')}  % rcrm={g('pct_crm')}\n"
    )
    out.write(f" Hlab={g('Hlab')}µT  b={g('b')}  sb={g('sb')}  sb/b={g('sb_over_b')}\n")
    if g("fcor") != "n.d":
        out.write(f" fcor={g('fcor')}  Hcorani={g('Hcorani')}µT\n")
    if g("fcorCool") != "n.d":
        out.write(f" fcorCool={g('fcorCool')}  HcorCool={g('HcorCool')}µT\n")
    out.write(f" gamma={g('gamma')}  k={g('k')}  k_sse={g('k_sse')}\n")
    out.write(
        f" mad={g('mad')}  dang={g('dang')}  fvds={g('fvds')}  frac={g('frac')}  "
        f"gap_max={g('gap_max')}  n_ptrm={g('n_ptrm')}  (PmagPy/MagIC)\n"
    )
    if g("f1") != "n.d":
        out.write(f" rf1 - rf2 : {g('f1')}  {g('f2')}\n")
    if g("tensor") != "n.d":
        out.write(f" anisotropy tensor applied: {g('tensor')}\n")
    return out.getvalue()


def _resource_path(*parts: str) -> str:
    """Chemin absolu d'une ressource livree AVEC l'appli (ex. le guide
    utilisateur HTML, voir ouvrir_user_guide) - fonctionne aussi bien
    lancee depuis les sources (repertoire de ce fichier) qu'empaquetee par
    PyInstaller (sys._MEIPASS, le dossier temporaire ou l'executable
    autonome extrait ses fichiers - voir STARpaleomag_Py.spec, datas=[
    ('help', 'help')])."""
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, *parts)


class STARpaleomagApp:
    def __init__(self, root):
        self.root = root
        self.root.title("STARpaleomag_Py - Paleomagnetism")
        self._maximize_window()

        self.donnees = []  # Stockage des données d'échantillons (List[Pmag])
        self.selection = []  # Dernière sélection (List[SelectedSample])
        self.results = []  # Ajustements de droite (List[FitResult], equivalent tr/nbres)
        self.results_path = None  # equivalent filr - fichier .r, fixe au chargement des donnees
        self.data_file_path = None  # chemin du .ren/.prmag/.txt charge - voir _load_data_file
        self._archived_ids = None  # cache des `c` deja utilises dans results_path (voir _save_result)
        self._results_path_mtime = None  # date de derniere lecture de results_path - voir ouvrir_selres_dialog
        self._arm_holder_background = None  # equivalent xholarm/yholarm/zholarm (holderarm), pour Anisotropy
        self.entete = ""  # Préfixe de sélection (equivalent selentete)
        self.orientation = tk.IntVar(value=2)  # equivalent iorient (SC/IS/TC) - defaut In situ
        self._current_graphic = None  # ("zijderveld", sample_id) / ("stereo"|"xygraph"|"susceptibility"|"arai", None)
        self._arai_state = None  # (ech, points, checks, arno, fit_ou_None, hlab) pour kind=="arai"
        self._paleoint_review_state = None  # (ech, points, checks, arno, fit) pour kind=="paleoint_review"

        self._setup_menu()
        self._setup_shortcuts()

        self.paned_window = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        self.paned_window.pack(fill=tk.BOTH, expand=True)

        # PANNEAU GAUCHE : Zone Graphique (matplotlib embarqué)
        self.graph_frame = ttk.Frame(self.paned_window, width=550)
        self.fig = Figure(figsize=(5.2, 5.2), dpi=100)
        self.canvas_fig = FigureCanvasTkAgg(self.fig, master=self.graph_frame)
        self.canvas_fig.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.paned_window.add(self.graph_frame, weight=1)
        # Clic-pour-info sur un point de graphique (Zijderveld/Stereo
        # Results) - demande explicite utilisateur ("un collegue m'a
        # demande si il etait possible de cliquer sur des donnees d'un
        # graphique pour obtenir des informations") - UN SEUL point de
        # connexion pour tous les graphiques (self.fig/self.canvas_fig
        # sont partages, voir PlotContext.pick_point pour comment chaque
        # point interactif est enregistre a la construction du graphique).
        self.canvas_fig.mpl_connect("pick_event", self._on_plot_pick)

        # PANNEAU DROIT : Console / Liste de données
        self.text_frame = ttk.Frame(self.paned_window, width=550)
        self.text_area = tk.Text(
            self.text_frame,
            bg="#ffffff",
            fg="#000000",
            insertbackground="black",
            font=("Courier", 14),
            wrap="none",
        )
        text_yscroll = ttk.Scrollbar(self.text_frame, orient=tk.VERTICAL, command=self.text_area.yview)
        text_xscroll = ttk.Scrollbar(self.text_frame, orient=tk.HORIZONTAL, command=self.text_area.xview)
        self.text_area.configure(yscrollcommand=text_yscroll.set, xscrollcommand=text_xscroll.set)
        self.text_area.tag_configure("prompt", foreground="#c0392b")
        self.text_area.tag_configure("header", font=("Courier", 14, "bold"))
        self.text_area.tag_configure(
            "warn", foreground="#c0392b", background="#fdf1d0", font=("Courier", 14, "bold"))
        self.text_area.grid(row=0, column=0, sticky="nsew")
        text_yscroll.grid(row=0, column=1, sticky="ns")
        text_xscroll.grid(row=1, column=0, sticky="ew")
        self.text_frame.rowconfigure(0, weight=1)
        self.text_frame.columnconfigure(0, weight=1)
        self.paned_window.add(self.text_frame, weight=1)

        self.load_demo_data()

        # Force l'activation de l'appli au niveau macOS. lift()/focus_force()
        # (niveau Tk) et le bascule topmost (niveau window manager) n'ont pas
        # suffi : la fenetre s'affiche et reagit a la souris mais l'appli ne
        # devient jamais le process "frontmost" pour le clavier tant qu'on
        # n'a pas clique sur la barre de menus. On demande donc directement
        # a System Events d'activer ce process par PID (osascript) - c'est
        # le mecanisme d'activation macOS reel, hors de Tk. Peut demander
        # une autorisation "Automatisation" la premiere fois (a accepter).
        self.root.after(200, self._activate_window)

    def _maximize_window(self):
        """Demarre en fenetre maximisee (occupant tout l'ecran disponible,
        SANS masquer la barre de menus - contrairement a un vrai plein
        ecran) - demande explicite utilisateur ("optimize the size of the
        window to full screen but access to the main menu. The user can
        reduce it if needed"). `state('zoomed')` est fiable sous Windows/
        Linux mais INEXISTANT sous Tk/Aqua (macOS) - repli la sur la taille
        d'ecran effective (winfo_screenwidth/height), qui laisse
        naturellement la barre de menus macOS (hors fenetre Tk) accessible.
        L'utilisateur reste libre de redimensionner/deplacer la fenetre
        ensuite (aucune contrainte imposee apres ce point de depart)."""
        try:
            if sys.platform.startswith("win"):
                self.root.state("zoomed")
                return
            if sys.platform.startswith("linux"):
                self.root.attributes("-zoomed", True)
                return
        except tk.TclError:
            pass
        self.root.update_idletasks()
        w = self.root.winfo_screenwidth()
        h = self.root.winfo_screenheight()
        self.root.geometry(f"{w}x{h}+0+0")

    def _activate_window(self):
        self.root.attributes("-topmost", True)
        self.root.after(50, lambda: self.root.attributes("-topmost", False))
        try:
            subprocess.run(
                [
                    "osascript", "-e",
                    f'tell application "System Events" to set frontmost of '
                    f'(first process whose unix id is {os.getpid()}) to true',
                ],
                check=False,
                capture_output=True,
                timeout=2,
            )
        except (subprocess.SubprocessError, OSError):
            pass

    @staticmethod
    def _labeled(text, shortcut_name):
        """Libelle avec le raccourci entre parentheses, SANS utiliser le
        parametre `accelerator=` de Tk : sur Aqua, ce dernier semble faire
        que macOS intercepte la combinaison au niveau du menu natif (la
        barre de menus s'active/se surligne) sans jamais invoquer la
        commande Tcl associee - ce qui court-circuite bind_all. Le texte
        seul, lui, est purement cosmetique."""
        return f"{text}    ({SHORTCUTS[shortcut_name][0]})"

    def _setup_menu(self):
        menubar = tk.Menu(self.root)

        # Menu PmagFile (titres repris de StarmacOSX_x.f95, bloc menu anglais
        # lignes 491-608 - AWE_addMenu(HelpMenuUnit,"PmagFile"/"Pmag data"/
        # "Graphics"/"Calcul", ...) - pas le bloc francais "Fichiers"/
        # "Selection donnees"/"Graphiques" utilise plus haut dans le fichier)
        # Reorganise a la demande explicite utilisateur ("est-ce aussi
        # possible de reorganiser le menu Pmag file") - regroupe par
        # fonction plutot qu'en une seule liste plate : ouverture, puis
        # conversions VERS le format .prmag (le format natif courant),
        # puis interoperabilite MagIC (import ET export - aucune de ces
        # commandes ne charge quoi que ce soit dans la session, meme
        # "Import..." : toutes ecrivent un nouveau fichier sur disque),
        # puis exports/rapports.
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label=self._labeled("Open Pmag file...", "importpc"),
                               command=self.ouvrir_fichier_ren)
        file_menu.add_separator()
        file_menu.add_command(label="Import Starmac legacy files...", command=self.ouvrir_convert_legacy_dialog)
        file_menu.add_separator()
        file_menu.add_command(label="Import MagIC contribution to .prmag format...",
                               command=self.ouvrir_convert_magic_to_r_dialog)
        file_menu.add_command(label="Import Utrecht/PMAG2 .col to .prmag format...",
                               command=self.ouvrir_convert_utrecht_to_r_dialog)
        file_menu.add_command(label="Import IPGP measurements/results to .prmag format...",
                               command=self.ouvrir_convert_ipgp_to_r_dialog)
        file_menu.add_command(label="Import Montpellier .montpellier to .prmag format...",
                               command=self.ouvrir_convert_montpellier_to_r_dialog)
        file_menu.add_separator()
        file_menu.add_command(label="Orientation: field notes orientation to prmag file...",
                               command=self.ouvrir_field_notes_dialog)
        file_menu.add_command(label="Create prmag from AGICO .ged file...",
                               command=self.ouvrir_ged_to_prmag_dialog)
        file_menu.add_command(label="Complete sample information...", command=self.ouvrir_complete_sample_info_dialog)
        file_menu.add_command(label="Archive new laboratory measurements", command=self.ouvrir_archive_new_data_dialog)
        file_menu.add_separator()
        file_menu.add_command(label="Check files for duplicate entries...", command=self.ouvrir_check_duplicates_dialog)
        file_menu.add_separator()
        file_menu.add_command(label="export to Magic database", command=self.ouvrir_export_magic_dialog)
        file_menu.add_separator()
        file_menu.add_command(label="export pmag content as a text file", command=self.ouvrir_export_detailed_dialog)
        file_menu.add_command(label="export results to Stereo_Py...", command=self.ouvrir_export_stereo_dialog)
        file_menu.add_separator()
        file_menu.add_command(label=self._labeled("Quit STARpaleomag_Py", "starend"), command=self.root.quit)
        menubar.add_cascade(label="PmagFile", menu=file_menu)

        # Menu Pmag data
        data_menu = tk.Menu(menubar, tearoff=0)
        data_menu.add_command(label=self._labeled("Select data...", "selmes"),
                               command=self.ouvrir_selection_dialog)
        data_menu.add_command(label="Select site...", command=self.ouvrir_selection_site_dialog)
        data_menu.add_separator()
        data_menu.add_command(label=self._labeled("Select header...", "selentete"),
                               command=self.ouvrir_entete_dialog)
        data_menu.add_separator()
        data_menu.add_command(label=self._labeled("Delete some data...", "effmes"),
                               command=self.ouvrir_effmes_dialog)
        data_menu.add_command(label=self._labeled("Init list", "initmes"),
                               command=self.reinitialiser_selection)
        data_menu.add_separator()
        data_menu.add_command(label=self._labeled("List data", "lismes"), command=self.lister_mesures)
        data_menu.add_command(label="List in XYZ", command=self.lister_xyz)
        data_menu.add_command(label="List data VRM", command=self.lister_vrm)
        data_menu.add_command(label="List and depth...", command=self.ouvrir_lismesdepth_dialog)
        data_menu.add_separator()
        data_menu.add_command(label=self._labeled("Info samples", "infoech"),
                               command=self.afficher_info_echantillons)
        data_menu.add_separator()

        # Radiobuttons directement dans le menu (pas de sous-menu imbrique) :
        # meme raisonnement que ci-dessus, en plus du probleme deja identifie
        # avec les sous-menus en cascade.
        for label, value in ORIENTATIONS.items():
            data_menu.add_radiobutton(
                label=self._labeled(label, ORIENTATION_SHORTCUT_NAMES[value]),
                variable=self.orientation, value=value,
                command=lambda v=value: self._set_orientation(v),
            )

        menubar.add_cascade(label="Pmag data", menu=data_menu)

        # Menu Results (regroupe ce qui etait eparpille entre Pmag
        # data/Calcul/Graphics - demande explicite utilisateur : "adding a
        # menu Results and another Paleointensity"). "Stereo Results"/
        # "data+interpretation" en sont ressortis vers Graphics dans une
        # reorganisation ulterieure (liste de menus/raccourcis fournie
        # explicitement par l'utilisateur).
        results_menu = tk.Menu(menubar, tearoff=0)
        results_menu.add_command(label=self._labeled("Select results...", "selres"),
                                  command=self.ouvrir_selres_dialog)
        results_menu.add_command(label=self._labeled("List results", "lisres"), command=self.lister_resultats)
        results_menu.add_command(label=self._labeled("Init results", "initres"),
                                  command=self.reinitialiser_resultats)
        results_menu.add_command(label="Delete results...", command=self.ouvrir_delete_results_dialog)
        results_menu.add_command(label="Import published site means...",
                                  command=self.ouvrir_import_published_means_dialog)
        results_menu.add_separator()
        results_menu.add_command(label=self._labeled("best lines...", "ajuslig"),
                                  command=self.ouvrir_ajuslig_dialog)
        results_menu.add_command(label="best lines auto", command=self.ouvrir_ajusligauto_dialog)
        results_menu.add_command(label="best planes...", command=self.ouvrir_ajusplans_dialog)
        results_menu.add_command(label=self._labeled("best dir Fisher...", "ajusfisher"), command=self.ouvrir_ajusfisher_dialog)
        results_menu.add_command(label="Best fit from redo file...", command=self.ouvrir_ajusligredo_dialog)
        results_menu.add_command(label="Auto-interpret (suggest components)...", command=self.ouvrir_autointerpretation_dialog)
        results_menu.add_separator()
        results_menu.add_command(label="Evaluate interpretations...", command=self.evaluer_interpretations)
        results_menu.add_command(label=self._labeled("Fisher results", "fishres"), command=self.fisher_resultats)
        menubar.add_cascade(label="Results", menu=results_menu)

        # Menu Paleointensity (demande explicite utilisateur : "view
        # paleoint va dans le menu paleointensite"). GRM va dans Calcul
        # ("Mettre GRM dans calcul"), Test Induite et ConvertZ- 2G retires
        # du menu ("supprimer test induite" / "remove the convert Z - 2G
        # menu" - fonctions conservees dans le code, simplement plus
        # exposees par un intitule de menu).
        paleoint_menu = tk.Menu(menubar, tearoff=0)
        paleoint_menu.add_command(label="Paleointensity interpretation", command=self.afficher_arai)
        paleoint_menu.add_command(label="View batch of Paleoint Results (from .pmagint, read-only)...",
                                   command=self.ouvrir_view_pmagint_dialog)
        paleoint_menu.add_command(label="Rapid view/recompute from redo file...",
                                   command=self.ouvrir_openfilepint_dialog)
        paleoint_menu.add_separator()
        paleoint_menu.add_command(label="Thellier >> NRM", command=self.ouvrir_convertthelli_dialog)
        paleoint_menu.add_command(label="Remove paleointensity step", command=self.ouvrir_removestep_dialog)
        paleoint_menu.add_command(label="Remove bad-quality (b) steps...", command=self.ouvrir_remove_bad_quality_dialog)
        paleoint_menu.add_command(label="export to ThellierTool...", command=self.ouvrir_exportthellier_dialog)
        paleoint_menu.add_separator()
        paleoint_menu.add_command(label="Cooling rate...", command=self.ouvrir_cooling_rate_dialog)
        menubar.add_cascade(label="Paleointensity", menu=paleoint_menu)

        # Menu Calcul (categorie gardee en francais dans la source elle-meme,
        # cf ligne 580 du bloc anglais - seuls les items sont traduits).
        # Items non encore implementes (au-dela d'ajuslig/fisher/resultats) :
        # stubs relies a self._not_implemented, memes intitules que la source,
        # SAUF "Orientation drill cores"/"Orient Core with LowTemp" (lies a
        # l'option Drillcore, retiree de l'appli).
        calcul_menu = tk.Menu(menubar, tearoff=0)
        calcul_menu.add_command(label=self._labeled("Fisher measures", "fishmes"), command=self.fisher_mesures)
        calcul_menu.add_command(label=self._labeled("Fisher results", "fishres"), command=self.fisher_resultats)
        calcul_menu.add_separator()
        calcul_menu.add_command(label="MdF-MdT", command=self.afficher_mdf)
        calcul_menu.add_command(label=self._labeled("Mean Intensity", "meanint"), command=self.afficher_mean_intensity)
        calcul_menu.add_command(label="Koenigsberger ratio...", command=self.ouvrir_koenigsberger_dialog)
        calcul_menu.add_separator()
        calcul_menu.add_command(label="Mean Inclination", command=self.afficher_mean_inclination)
        calcul_menu.add_separator()
        calcul_menu.add_command(label="Test viscosity", command=self.appliquer_viscosity_test)
        calcul_menu.add_command(label="Diff measurements n/n-1", command=self.afficher_diff_measurements)
        calcul_menu.add_command(label="Subtraction...", command=self.ouvrir_subtraction_dialog)
        calcul_menu.add_command(label="Autoinverse", command=lambda: self._not_implemented("Autoinverse"))
        calcul_menu.add_separator()
        calcul_menu.add_command(label=self._labeled("Anisotropy", "anisotropy"), command=self.ouvrir_anisotropy_dialog)
        calcul_menu.add_command(label="Anisotropy PmagPy...", command=self.ouvrir_anisotropy_pmagpy_dialog)
        calcul_menu.add_command(label="Holder_ARM...", command=self.ouvrir_holderarm_dialog)
        calcul_menu.add_command(label="Inverse_ANI_correction...", command=self.ouvrir_inverseani_dialog)
        calcul_menu.add_separator()
        calcul_menu.add_command(label="Detect GRM...", command=self.ouvrir_detect_grm_dialog)
        calcul_menu.add_command(label="Suppress GRM", command=self.ouvrir_elimine_grm_dialog)
        menubar.add_cascade(label="Calcul", menu=calcul_menu)

        # Menu Graphics - "Stereo Results"/"data+interpretation" y sont
        # entres depuis Results (liste de menus/raccourcis fournie
        # explicitement par l'utilisateur).
        graph_menu = tk.Menu(menubar, tearoff=0)
        graph_menu.add_command(label=self._labeled("Zijderveld", "plotzijder"), command=self.afficher_zijderveld)
        graph_menu.add_command(label=self._labeled("Stereo data", "stereodata"), command=self.afficher_stereo)
        graph_menu.add_command(label=self._labeled("XYgraph", "xygraph"), command=self.afficher_xygraph)
        graph_menu.add_command(label=self._labeled("Susceptibility", "suscep"), command=self.afficher_susceptibilite)
        graph_menu.add_command(label="Plot IRM", command=self.afficher_irm)
        graph_menu.add_separator()
        graph_menu.add_command(label=self._labeled("Stereo Results", "stereores"), command=self.afficher_stereo_results)
        graph_menu.add_command(label="data+interpretation", command=self.afficher_visres)
        graph_menu.add_separator()
        graph_menu.add_command(label="Clear Screen", command=self.clear_screen)
        graph_menu.add_separator()
        graph_menu.add_command(label="Export SVG...", command=self.exporter_svg)
        menubar.add_cascade(label="Graphics", menu=graph_menu)

        # Aide statique, locale, sans cle API ni cout recurrent (option
        # retenue par l'utilisateur parmi 3 formes possibles d'"aide en
        # ligne par Claude" - "OK pour la 3", le guide HTML deja ecrit,
        # livre avec l'appli plutot qu'un assistant en direct qui aurait
        # exige soit une cle API a la charge de chaque utilisateur, soit
        # un serveur relais a heberger pour ne jamais exposer une cle
        # partagee dans le code source public).
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="User Guide", command=self.ouvrir_user_guide)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.root.config(menu=menubar)

    def _setup_shortcuts(self):
        """Raccourcis clavier globaux (voir SHORTCUTS/_SHORTCUTS_MAC/
        _SHORTCUTS_WIN). Sur Mac, tous utilisent Command comme modificateur
        de base (jamais Control seul) ; sur Windows/Linux, tous utilisent
        Ctrl+Shift (jamais Ctrl seul) - dans les deux cas, aucun conflit avec
        les raccourcis d'edition Emacs integres a Tk pour les widgets Entry/
        Text (Control+E, Control+D... eux sont en Control seul), ni avec les
        raccourcis universels Windows (Ctrl+S/F/Z/Y/A...) : pas besoin
        d'ignorer les raccourcis selon le widget qui a le focus."""
        bindings = {
            "importpc": self.ouvrir_fichier_ren,
            "starend": self.root.quit,
            "selmes": self.ouvrir_selection_dialog,
            "selentete": self.ouvrir_entete_dialog,
            "effmes": self.ouvrir_effmes_dialog,
            "initmes": self.reinitialiser_selection,
            "lismes": self.lister_mesures,
            "infoech": self.afficher_info_echantillons,
            "selce": lambda: self._set_orientation(1),
            "selis": lambda: self._set_orientation(2),
            "selcp": lambda: self._set_orientation(3),
            "ajuslig": self.ouvrir_ajuslig_dialog,
            "fishmes": self.fisher_mesures,
            "fishres": self.fisher_resultats,
            "lisres": self.lister_resultats,
            "selres": self.ouvrir_selres_dialog,
            "initres": self.reinitialiser_resultats,
            "plotzijder": self.afficher_zijderveld,
            "xygraph": self.afficher_xygraph,
            "stereodata": self.afficher_stereo,
            "stereores": self.afficher_stereo_results,
            "ajusfisher": self.ouvrir_ajusfisher_dialog,
            "suscep": self.afficher_susceptibilite,
            "meanint": self.afficher_mean_intensity,
            "anisotropy": self.ouvrir_anisotropy_dialog,
        }
        for name, callback in bindings.items():
            _, sequence = SHORTCUTS[name]
            self.root.bind_all(sequence, lambda event, cb=callback: cb())

    def _set_orientation(self, value):
        self.orientation.set(value)
        self._refresh_current_graphic()

    def _write_svg_debug(self, ech, fits_for_sample):
        """Pour test/comparaison : ecrit le rendu SVGWriter (fidele a
        svginit/svgplot) du Zijderveld affiche, AVANT le dessin matplotlib -
        pas encore le pipeline d'affichage definitif, juste un moyen rapide
        de verifier le SVG genere contre un export Starmac_OSX original."""
        try:
            os.makedirs(SVG_DEBUG_DIR, exist_ok=True)
            # page 19x28cm, origine (90,600)px : valeurs reelles de `plots()`
            # (graphicsAWE.f95) - draw_zijderveld part de cette origine de
            # PAGE pour appliquer lui-meme le decalage initial (-2,-5.5) puis
            # l'origine du NEV (u,v), comme le fait `zijder2`.
            writer = SVGWriter(width_cm=19.0, height_cm=28.0)
            writer.set_origin_px(90.0, 600.0)
            draw_zijderveld(writer, ech, orientation=self.orientation.get(), fits=fits_for_sample)
            writer.plotnd()
            path = os.path.join(SVG_DEBUG_DIR, f"zijder-{ech.id}.svg")
            writer.save(path)
        except Exception as e:
            print(f"[svg_debug] echec ecriture SVG pour {ech.id} : {e}")

    def _refresh_current_graphic(self):
        if self._current_graphic is None:
            return
        kind, extra = self._current_graphic
        self._clear_figure()

        if kind == "zijderveld":
            ech = next((s for s in self.selection if s.id == extra), None)
            if ech is None:
                return
            fits_for_sample = [r for r in self.results if r.id == ech.id]
            self._write_svg_debug(ech, fits_for_sample)
            self.fig.set_size_inches(5.5, 8.5, forward=True)
            build_zijderveld_figure(ech, orientation=self.orientation.get(), fits=fits_for_sample, fig=self.fig)
            # liste des donnees de l'echantillon affiche dans la fenetre
            # texte (PAS sur le graphique) - demande explicite de
            # l'utilisateur, meme format que "List data" (Cmd+L).
            buffer = io.StringIO()
            list_measurements([ech], orientation=self.orientation.get(), out=buffer)
            self._afficher(buffer.getvalue())
        elif kind == "stereo":
            if not self.selection:
                return
            self.fig.set_size_inches(5.5, 5.5, forward=True)
            build_stereo_figure(self.selection, orientation=self.orientation.get(), fig=self.fig)
        elif kind == "stereo_results":
            if not self.results:
                return
            self.fig.set_size_inches(5.5, 5.5, forward=True)
            build_stereo_results_figure(
                self.results, orientation=self.orientation.get(),
                nbech=len(self.selection), fig=self.fig,
            )
        elif kind == "xygraph":
            if not self.selection:
                return
            # empile verticalement (pas cote a cote) : le panneau graphique
            # a une largeur fixe, on agrandit la HAUTEUR, pas la largeur.
            height = 8.5 if has_mixed_demag(self.selection) else 4.5
            self.fig.set_size_inches(6.0, height, forward=True)
            build_xygraph_figure(self.selection, fig=self.fig)  # appelle deja fig.tight_layout()
            self._redraw_canvas()
            return
        elif kind == "susceptibility":
            if not self.selection:
                return
            self.fig.set_size_inches(6.0, 4.5, forward=True)
            build_susceptibility_figure(self.selection, fig=self.fig)  # appelle deja fig.tight_layout()
            self._redraw_canvas()
            return
        elif kind == "irm":
            if not self.selection:
                return
            self.fig.set_size_inches(6.0, 8.0, forward=True)
            build_irm_figure(self.selection, fig=self.fig)  # appelle deja fig.tight_layout()
            self._redraw_canvas()
            return
        elif kind == "arai":
            if self._arai_state is None:
                return
            ech, points, checks, arno, fit, _hlab = self._arai_state
            self.fig.set_size_inches(6.0, 8.0, forward=True)
            build_arai_figure(ech, points, checks, arno, fit=fit, fig=self.fig)
        elif kind == "paleoint_review":
            if self._paleoint_review_state is None:
                return
            ech, points, checks, arno, fit = self._paleoint_review_state
            # page Fortran (visi_Paleoint.f: call plots(19.5,28.,fname)) -
            # remplace l'ancien (7.5,9.5) qui ecrasait le figsize deja
            # corrige a l'interieur de build_paleoint_review_figure
            # (celui-ci ne s'applique que si fig=None, or self.fig existe
            # deja ici) - demande explicite utilisateur ("each of the
            # three could be larger within the page").
            self.fig.set_size_inches(19.5 / 2.54, 28.0 / 2.54, forward=True)
            build_paleoint_review_figure(
                ech, points, checks, arno, fit=fit,
                # orientation FORCEE a 1 (echantillon), PAS
                # self.orientation.get() : le Fortran (plotpaleoint2.f et
                # visi_Paleoint.f) force "iorient=1" sans condition en
                # tete de cette routine, quel que soit le reglage
                # d'orientation par ailleurs dans l'appli - demande
                # explicite utilisateur ("the combined plots should be in
                # sample coordinates").
                orientation=1, fig=self.fig,
            )
            # PAS _fit_figure_to_data() : ecrit pour un panneau SIMPLE
            # (ax.axes[0] uniquement) - sur cette figure a 3 panneaux, elle
            # redimensionnerait toute la figure sur le seul rapport
            # largeur/hauteur du panneau Arai (le premier), ecrasant la
            # taille/mise en page deja correctes de build_paleoint_review_
            # figure (page Fortran 19.5x28cm, subplots_adjust) - meme
            # constat que xygraph/susceptibility/irm ci-dessus, qui
            # bypassent deja _fit_figure_to_data pour la meme raison -
            # demande explicite utilisateur ("each of the three could be
            # larger within the page").
            self._redraw_canvas()
            return
        else:
            return

        self._fit_figure_to_data()
        self._redraw_canvas()

    def _on_plot_pick(self, event):
        """Gestionnaire unique de clic-sur-point (pick_event), pour tous
        les graphiques - demande explicite utilisateur ("un collegue m'a
        demande si il etait possible de cliquer sur des donnees d'un
        graphique pour obtenir des informations, par exemple sur un
        zijderveld, la temperature, ou sur un stereo results, le nom du
        specimen"). Chaque point interactif est enregistre a la
        construction du graphique via PlotContext.pick_point (zijderveld.
        draw_zijderveld / stereo.draw_stereo_results), qui attache
        `_starmac_pick_kind`/`_starmac_pick_data` a un scatter invisible -
        ce gestionnaire n'a qu'a les relire, sans redecoder les
        coordonnees affichees. `event.ind[0]` : index du point le plus
        proche du clic parmi ceux de CET artiste (un artiste par point
        ici, donc toujours 0, mais reste correct si un futur point
        interactif regroupait plusieurs positions dans un seul scatter)."""
        artist = event.artist
        kind = getattr(artist, "_starmac_pick_kind", None)
        data = getattr(artist, "_starmac_pick_data", None)
        # BUG REEL corrige ici (signale par l'utilisateur, "does not show
        # anything on my macbook air") : `event.ind` est un tableau numpy
        # (ex. array([0])) - `not event.ind` teste la VALEUR de son unique
        # element (0 -> False) plutot que si le tableau est vide, contraire
        # a l'intuition d'une liste Python. Comme chaque artiste
        # interactif ne porte qu'UN SEUL point (voir PlotContext.
        # pick_point), l'index trouve est TOUJOURS 0 - `not event.ind`
        # valait donc TOUJOURS True, sortant de la fonction avant meme
        # d'afficher quoi que ce soit, sur CHAQUE clic. Diagnostique via
        # un test manuel de figure.pick()/artist.contains() (les deux
        # reussissaient) montrant que seul ce garde-fou bloquait l'affichage.
        if kind is None or len(event.ind) == 0:
            return
        idx = event.ind[0]
        item = data[idx] if isinstance(data, (list, tuple)) else data

        if kind == "zijderveld_step":
            m = item
            self._afficher(f"[Zijderveld] step {_step_token(m)}\n")
        elif kind == "stereo_specimen":
            self._afficher(f"[Stereo Results] {item}\n")
        elif kind == "xygraph_specimen":
            self._afficher(f"[XY graph] {item}\n")
        elif kind == "irm_specimen":
            self._afficher(f"[IRM] {item}\n")

    def _redraw_canvas(self):
        """Agrandit le CONTENANT (volet gauche du PanedWindow + fenetre)
        pour accueillir la Figure demandee, plafonne a l'espace ecran
        reellement disponible, force un vrai cycle de redimensionnement -
        puis MESURE la taille reellement obtenue par le widget canvas et
        cale la Figure dessus, plutot que l'inverse.

        Lecon tiree de plusieurs iterations infructueuses : matplotlib
        traite le WIDGET Tk comme source de verite pour la taille de la
        Figure (`FigureCanvasTkAgg` resynchronise la Figure dessus des
        qu'un vrai evenement <Configure> se produit) - essayer d'imposer
        une taille de Figure malgre ca revient a se battre contre le
        framework : la moindre resynchronisation ulterieure (y compris
        celle qu'on declenche nous-meme en "secouant" la fenetre pour
        forcer Tk a se mettre a jour) ecrase la valeur programmee. La bonne
        approche est donc : dimensionner le CONTENANT selon ce qu'on veut
        afficher, puis lire sa taille REELLE une fois stabilisee, et
        aligner la Figure dessus - jamais l'inverse."""
        width_in, height_in = self.fig.get_size_inches()
        dpi = self.fig.dpi

        max_width_px = max(400, self.root.winfo_screenwidth() - 100)
        max_height_px = max(300, self.root.winfo_screenheight() - 100)
        avail_fig_w_in = (max_width_px - 20) / dpi
        avail_fig_h_in = (max_height_px - 90) / dpi
        if width_in > avail_fig_w_in or height_in > avail_fig_h_in:
            scale = min(avail_fig_w_in / width_in, avail_fig_h_in / height_in, 1.0)
            width_in, height_in = width_in * scale, height_in * scale

        width_px = int(width_in * dpi) + 20
        height_px = int(height_in * dpi) + 90  # marge titre/menus/barre d'etat

        try:
            current_w = self.paned_window.sashpos(0)
            if width_px > current_w:
                self.paned_window.sashpos(0, width_px)
        except tk.TclError:
            pass

        if height_px > self.root.winfo_height():
            self.root.geometry(f"{self.root.winfo_width()}x{height_px}")

        # laisse Tk appliquer sashpos()/geometry() avant de lire la taille
        # resultante - sinon winfo_width/height() ci-dessous renverraient
        # encore l'ancienne valeur (la demande est asynchrone).
        self.root.update_idletasks()

        # "secoue" la fenetre pour forcer un vrai cycle <Configure> sur
        # tous les widgets enfants (constate : ni update_idletasks(), ni
        # sashpos()/geometry() seuls ne suffisent - seul un VRAI
        # redimensionnement le fait, d'ou "toucher la fenetre" qui
        # corrige toujours le probleme manuellement).
        w = self.root.winfo_width()
        h = self.root.winfo_height()
        self.root.geometry(f"{w + 1}x{h + 1}")
        self.root.update_idletasks()
        self.root.geometry(f"{w}x{h}")
        self.root.update_idletasks()

        # Aligne la Figure sur la taille REELLE obtenue par le widget
        # (source de verite), pas sur ce qu'on voulait au depart.
        canvas_w = self.canvas_fig.get_tk_widget().winfo_width()
        canvas_h = self.canvas_fig.get_tk_widget().winfo_height()
        if canvas_w > 1 and canvas_h > 1:
            self.fig.set_size_inches(canvas_w / dpi, canvas_h / dpi, forward=False)

        self.canvas_fig.draw()

    def _fit_figure_to_data(self, width_in: float = 5.5, max_height_in: float = 14.0):
        """Cale la taille de la Figure sur le vrai rapport largeur/hauteur des
        donnees tracees (mesure sur `ax.get_xlim()/get_ylim()` apres l'appel a
        `build_*_figure`, qui a deja fait `relim()+autoscale_view()`).

        Necessaire car `PlotContext.clear()` impose `ax.set_aspect('equal')` :
        si la Figure n'a PAS le meme rapport largeur/hauteur que les donnees,
        matplotlib retrecit la zone de tracee (`adjustable='box'`) pour
        respecter l'aspect - le dessin se retrouve reduit a une bande etroite
        au milieu d'une figure trop large/haute, avec de grandes marges vides
        (constate pour le Zijderveld+stereo, dont l'etendue verticale relle
        est bien plus grande que ce qu'une figure 5.5x8.5 fixe peut refleter)."""
        if not self.fig.axes:
            return
        ax = self.fig.axes[0]
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        w = xlim[1] - xlim[0]
        h = ylim[1] - ylim[0]
        if w <= 0 or h <= 0:
            return
        height_in = min(width_in * (h / w), max_height_in)
        self.fig.set_size_inches(width_in, height_in, forward=True)
        self.fig.tight_layout()

    # ------------------------------------------------------------------
    # Fichiers
    # ------------------------------------------------------------------

    def ouvrir_fichier_ren(self):
        """Ouvre une boîte de dialogue pour sélectionner un fichier .ren
        ou .prmag et charge ses données.

        Accepte aussi .txt : anciens formats de fichier (ex. Domeyko - une
        seule ligne 'Id:', pas de 'L:'/roche ; Corbières - 'L:' + roche
        collee sans guillemets en fin de ligne plutot que sur une ligne a
        part) que `read_ren_file` sait maintenant lire tout aussi bien.

        .prmag (nouveau format, voir convert_ren_to_r.py/testlect.
        read_prmag_file) est lu par un parseur different mais produit la
        MEME structure List[Pmag] - tout le reste de l'application
        continue de fonctionner sans changement, quel que soit le format
        ouvert - demande explicite utilisateur ("the open file still
        look for .ren")."""
        fichier_path = filedialog.askopenfilename(
            title="Select a .ren, .prmag or .txt data file",
            filetypes=[
                ("REN/PRMAG/TXT files", "*.ren *.prmag *.txt"),
                ("PRMAG files", "*.prmag"),
                ("REN/TXT files", "*.ren *.txt"),
                ("All files", "*.*"),
            ],
        )
        if fichier_path:
            self._load_data_file(fichier_path)

    def _load_data_file(self, fichier_path: str, announce: bool = True) -> bool:
        """Charge `fichier_path` (.ren/.prmag/.txt) dans self.donnees/
        self.results_path et l'affiche - factorise depuis
        ouvrir_fichier_ren pour etre reutilisable a la fin des dialogues
        de conversion (Import Starmac legacy files/MagIC/Utrecht) -
        demande explicite utilisateur ("ce serait bien d'ouvrir les
        fichiers convertis a la fin des conversions"). `announce=False`
        omet le messagebox "Success" (utilise apres une conversion, dont
        le propre messagebox "Conversion complete" a deja informe
        l'utilisateur - eviter 2 popups consecutifs pour la meme action).
        Retourne True si le chargement a reussi."""
        try:
            if fichier_path.lower().endswith(".prmag"):
                self.donnees = read_prmag_file(fichier_path)
            else:
                self.donnees = read_ren_file(fichier_path)
            self.selection = []
            # equivalent filr=fil1(1:(nlen-4))//'.r' : fichier resultats
            # derive du fichier de donnees (meme dossier, extension .r)
            self.results_path = results_path_for(fichier_path)
            self.data_file_path = fichier_path
            self._archived_ids = None  # relit le fichier .r au prochain archivage

            # Réinitialisation et affichage dans la zone de texte
            self.text_area.delete("1.0", tk.END)
            self.text_area.insert(
                tk.END,
                f"File loaded: {os.path.basename(fichier_path)}\n"
                f"Number of samples: {len(self.donnees)}\n"
                "----------------------------------------\n\n"
            )

            for sample in self.donnees:
                self.text_area.insert(
                    tk.END,
                    f"ID: {sample.id} | Measurements: {sample.nbmes} | Inc: {sample.cin} | Az: {sample.caz}\n"
                )

            # Controle de coherence .pmagani <-> .prmag DES L'OUVERTURE
            # (pas seulement au moment d'un export MagIC) - demande
            # explicite utilisateur ("important pas seulement pour
            # l'export dans Magic mais il faut un controle lors de
            # l'ouverture du fichier pour verifier la compatibilite
            # pmagani <-> prmag") : un .pmagani compagnon (meme nom de
            # base, voir ani_path_for) portant un specimen absent du
            # .prmag qu'on vient de charger signale un vrai desaccord
            # entre les deux fichiers - averti tout de suite plutot que
            # decouvert seulement en export.
            ani_path = ani_path_for(fichier_path)
            if os.path.exists(ani_path):
                orphans = find_orphan_ani_specimens(self.donnees, ani_path)
                if orphans:
                    shown = ", ".join(orphans[:20])
                    more = f", ... ({len(orphans) - 20} more)" if len(orphans) > 20 else ""
                    self.text_area.insert(
                        tk.END,
                        f"\nWARNING: {len(orphans)} specimen(s) in {os.path.basename(ani_path)} "
                        f"have NO counterpart in this .prmag - mismatched files? "
                        f"{shown}{more}\n"
                    )

            if announce:
                self._showinfo(
                    "Success", f"{len(self.donnees)} sample(s) loaded!"
                )
            return True
        except Exception as e:
            self._showerror("Error", f"Could not read the file:\n{e}")
            return False

    def ouvrir_convert_legacy_dialog(self):
        """Equivalent GUI de Convert2newpmagformat/convert_oldpmag_2ren.f95 :
        met a niveau un vieux fichier de donnees (1, 2 ou 3 lignes d'entete
        par echantillon, eventuellement MELANGEES dans le meme fichier)
        vers le format actuel a 3 lignes, en completant Site/Sample/Fm/
        Age/GC/SMT/Li/Loc/Obs depuis un fichier "complement" externe
        OPTIONNEL. Un seul point d'entree, `convert_legacy_auto` - demande
        explicite utilisateur ("single import legacy files (able to
        import 1 ligne, 2 lignes and 3 lines header) - the app will
        automatically recognize if there is one, two or three lines") :
        remplace l'ancien choix manuel "Case 1/2/3" (les 3 sous-routines
        dediees convert_case1_specimen_info/convert_case2_site_info/
        convert_case3_new_format restent disponibles directement dans
        convert_legacy_ren.py pour un complement SANS en-tete reconnu,
        seul cas que l'auto-detection ne sait pas lever - voir sa
        docstring). Seul l'ancien "cas 3" (fichier deja a 2 lignes) a pu
        etre verifie octet pres contre un exemple reel fourni par
        l'utilisateur (Convert2newpmagformat/old_pmag.txt + complement.txt
        -> old_pmag.ren).

        Le complement est OPTIONNEL - demande explicite utilisateur
        ("import these files as they are and later we can complete the
        missing information... most of the tasks within the program just
        need the core correction, bedding and volume or mass... most of
        the missing information in the oldest file are for an export to
        Magic") : la ligne `Id:` du vieux fichier porte deja in:/az:/dip:/
        str:/v: (correction de carotte, pendage, volume/masse), jamais
        touchee par le complement - celui-ci ne sert qu'a renseigner
        Site/Fm/Age/GC/SMT/Li/Loc/Obs pour un futur export MagIC. Sans
        complement (ou pour un specimen sans correspondance), la
        conversion procede quand meme, avec "Not Specified" pour les
        champs geologie - a completer plus tard par une routine dediee
        plutot que de bloquer l'import maintenant.

        Enchaine AUSSI la conversion .ren -> .prmag (convert_ren_to_r.
        convert_file) - demande explicite utilisateur ("convert all types
        of legacy files to the .prmag format") : jusqu'ici il fallait
        relancer manuellement "Convert .ren to new format .r..." sur le
        fichier "..._converted.ren" tout juste produit. Le .ren
        intermediaire est CONSERVE sur disque (pas juste un fichier de
        travail jete) : c'est le point d'entree naturel pour completer a
        la main les champs "Not Specified"/les codes non reconnus avant
        de reconvertir en .prmag, sans avoir a repartir du fichier legacy
        d'origine."""
        self.text_area.insert(tk.END, "\n--- Import legacy files (Escape to cancel) ---\n", "prompt")

        old_path = filedialog.askopenfilename(
            title="Select the old data file (1, 2 or 3 header lines)",
            filetypes=[("Text/REN", "*.txt *.ren"), ("All files", "*.*")],
        )
        if not old_path:
            return

        complement_path = None
        if messagebox.askyesno(
            "Complement file",
            "Provide a complement file (Site/Formation/Age/... for a future MagIC export)?\n\n"
            "Not required: core correction, bedding and volume/mass are already in the data "
            "file and will be imported either way.\n\n"
            "The complement file must have a recognized column header (e.g. a first line "
            "'specimen\\tsite\\tsample\\tlat\\tlon\\t...' or 'site\\tformation\\tage\\t...').",
        ):
            complement_path = filedialog.askopenfilename(
                title="Select the complement file",
                filetypes=[("Text", "*.txt"), ("All files", "*.*")],
            )

        base, _ext = os.path.splitext(old_path)
        ren_path = base + "_converted.ren"
        prmag_path = base + "_converted.prmag"

        try:
            nb, unmatched, code_anomalies, skipped_multi_suffix = convert_legacy_auto(
                old_path, complement_path, ren_path)
        except Exception as e:
            self._showerror("Error", f"Conversion failed:\n{e}")
            return

        legacy_results_path = base + ".r"
        try:
            nb_prmag, nb_results, results_warnings = convert_ren_to_r_file(
                ren_path, prmag_path, legacy_results_path=legacy_results_path)
        except Exception as e:
            self._showerror("Error", f"Conversion to .prmag failed:\n{e}")
            return

        msg = f"Converted: {nb} sample(s) -> {ren_path}\n"
        msg += f"Converted: {nb_prmag} sample(s) -> {prmag_path}\n"
        if nb_results:
            msg += (
                f"Converted: {nb_results} result(s) -> {results_path_for(prmag_path)} "
                f"(anti-collision id renumbered as <specimen>_<letter>)\n"
            )
        if results_warnings:
            msg += f"WARNING: {len(results_warnings)} site mean(s) with an incomplete codes list:\n"
            for w in results_warnings[:20]:
                msg += f"  - {w}\n"
            if len(results_warnings) > 20:
                msg += f"  ... and {len(results_warnings) - 20} more (see console).\n"
                print("\n".join(results_warnings))
        missing_info_path = None
        if unmatched:
            # Liste COMPLETE des echantillons (sample, pas specimen - un
            # "95VN1901A"+"95VN1901B" ne comptent qu'une fois) sans
            # correspondance dans le complement, ecrite dans un fichier
            # texte a part - demande explicite utilisateur ("provide in a
            # text file the list of sample or site... to let the user
            # have a complete list of samples with missing information") :
            # la console tronque a 20 exemples, insuffisant pour aller
            # completer un fichier complement ou une table
            # "Complete sample information" sur un import a des centaines
            # de specimens sans correspondance.
            samples = sorted({s[:-1] if len(s) > 1 else s for s in unmatched})
            missing_info_path = base + "_missing_info.txt"
            with open(missing_info_path, "w", encoding="utf-8") as f:
                f.write("\n".join(samples) + "\n")
            shown = ", ".join(unmatched[:20])
            more = f", ... ({len(unmatched) - 20} more)" if len(unmatched) > 20 else ""
            msg += (f"No complement match ({len(unmatched)} specimen(s), {len(samples)} sample(s)), "
                    f"geology fields set to \"Not Specified\" (to complete later).\n"
                    f"Full sample list written to: {missing_info_path}\n"
                    f"  e.g. {shown}{more}\n")
        if code_anomalies:
            by_pair = {}
            for specimen, etape, cod1, cod2 in code_anomalies:
                by_pair.setdefault((cod1, cod2), []).append(f"{specimen} (step {etape})")
            msg += (f"\nUnrecognized cod1/cod2 ({len(code_anomalies)} measurement(s), "
                    f"{len(by_pair)} distinct code(s)) - please check and replace these in "
                    f"the source file:\n")
            for (cod1, cod2), examples in sorted(by_pair.items(), key=lambda kv: -len(kv[1])):
                shown = ", ".join(examples[:5])
                more = f", ... ({len(examples) - 5} more)" if len(examples) > 5 else ""
                msg += f"  '{cod1}{cod2}' x{len(examples)}: {shown}{more}\n"
        skipped_path = None
        if skipped_multi_suffix:
            # Specimens EXCLUS du fichier de sortie (pas juste "a
            # completer plus tard" comme unmatched/code_anomalies ci-
            # dessus) - demande explicite utilisateur ("do not import
            # samples with multi-char suffices, just provide the list of
            # non imported samples") : liste COMPLETE ecrite a part, meme
            # raison que missing_info_path (console tronquee a 20).
            skipped_path = base + "_not_imported.txt"
            with open(skipped_path, "w", encoding="utf-8") as f:
                f.write("\n".join(sorted(skipped_multi_suffix)) + "\n")
            shown = ", ".join(skipped_multi_suffix[:20])
            more = f", ... ({len(skipped_multi_suffix) - 20} more)" if len(skipped_multi_suffix) > 20 else ""
            msg += (f"\nNOT imported ({len(skipped_multi_suffix)} specimen(s) with a multi-character "
                    f"suffix, e.g. \"10CL0601AP\" - sample name can't be reliably derived):\n"
                    f"Full list written to: {skipped_path}\n"
                    f"  e.g. {shown}{more}\n")
        # Ouvre directement le .prmag converti - demande explicite
        # utilisateur ("ce serait bien d'ouvrir les fichiers convertis a
        # la fin des conversions"). AVANT self._afficher(msg) : le
        # chargement efface la zone de texte (voir _load_data_file), le
        # rapport de conversion (unmatched/code_anomalies/skipped) doit
        # etre affiche APRES, pas remplace par la simple liste de
        # specimens charges.
        self._load_data_file(prmag_path, announce=False)
        self._afficher(msg)
        info = f"{nb} sample(s) converted.\nOutput: {ren_path}\n{prmag_path}"
        if unmatched:
            info += (f"\n{len(unmatched)} specimen(s) had no complement match - geology fields "
                      f"left \"Not Specified\".\nFull sample list: {missing_info_path}")
        if code_anomalies:
            info += (f"\n{len(code_anomalies)} measurement(s) with an unrecognized cod1/cod2 "
                      f"- please check and fix them in the source file (see console).")
        if skipped_multi_suffix:
            info += (f"\n{len(skipped_multi_suffix)} specimen(s) NOT imported (multi-character suffix)."
                      f"\nFull list: {skipped_path}")
        self._showinfo("Conversion complete", info)

    def ouvrir_complete_sample_info_dialog(self):
        """Complete les metadonnees Site/Formation/Age/GC/SMT/Li/Loc/Obs
        (et lat/lon, et en mode specimen sample/site) d'un .prmag DEJA
        converti, depuis une table externe - demande explicite
        utilisateur ("can we put a menu like complete sample information
        with data in a table"), suite du "we can build a routine to
        complete the sample information later on" evoque lors du passage
        a un import legacy "as-is" (voir complete_sample_info.py pour le
        detail exact des trois modes et le format de table attendu).
        Mode 'h' (sample/height) ajoute a la demande explicite utilisateur
        ("add a specific case for stratigraphic_height filled from sample
        and height as all specimens have the same height"). Mode 's'
        accepte aussi une table separee par des virgules et une seconde
        colonne 'site' pour RENOMMER le site (ancien nom -> nouveau nom)
        - demande explicite utilisateur ("i cannot select site and
        replace it by an other site"), voir complete_sample_info.py.
        Patche le .prmag EN PLACE (une sauvegarde .bak est ecrite avant
        toute modification, une seule fois). Rapporte aussi, en plus des
        specimens sans correspondance dans la table, ceux qui ont encore
        au moins un champ a "n.d" APRES completion (table incomplete ou
        specimen non trouve) - demande explicite utilisateur ("it will be
        good to know what are the samples with missing information, when
        we complete the file"). Les deux listes sont ecrites en entier
        dans un fichier compagnon "<base>_completion_report.txt" (pas
        tassees sur une ligne de console tronquee a 20 - demande
        explicite utilisateur, cas reel a 251 specimens sans
        correspondance : "can you write all missing in a file instead of
        a long line")."""
        self.text_area.insert(tk.END, "\n--- Complete sample information (Escape to cancel) ---\n", "prompt")
        mode = self._console_input(
            "Table indexed by site, repeat 'site' header twice to rename (s) "
            "/ by specimen (also fills sample/site) (p) "
            "/ by sample, stratigraphic_height only (h): ", "s")
        if mode is None:
            return
        mode = (mode.strip().lower() or "s")[:1]
        if mode not in ("s", "p", "h"):
            self._showerror("Error", "Must be 's' (site), 'p' (specimen) or 'h' (sample/height).")
            return

        prmag_path = filedialog.askopenfilename(
            title="Select the .prmag file to complete",
            filetypes=[("prmag", "*.prmag"), ("All files", "*.*")],
        )
        if not prmag_path:
            return
        table_path = filedialog.askopenfilename(
            title="Select the information table (tab-separated)",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")],
        )
        if not table_path:
            return

        func = {"s": complete_site_info, "p": complete_specimen_info, "h": complete_sample_height}[mode]
        try:
            n_updated, unmatched, still_missing = func(prmag_path, table_path)
        except Exception as e:
            self._showerror("Error", f"Completion failed:\n{e}")
            return

        # Les listes (potentiellement des centaines de specimens - cas reel
        # signale par l'utilisateur, 251 "no table match") sont ecrites
        # dans un fichier compagnon plutot que tassees sur une ligne de
        # console tronquee a 20 - demande explicite utilisateur ("can you
        # write all missing in a file instead of a long line"), meme
        # convention de nommage que write_diagnostics_report
        # (<base>_diagnostics.txt) plus haut dans ce module.
        report_path = None
        if unmatched or still_missing:
            base, _ext = os.path.splitext(prmag_path)
            report_path = base + "_completion_report.txt"
            report_lines = [f"Completion report for {prmag_path}"]
            if unmatched:
                report_lines.append(f"\nNo table match ({len(unmatched)}):")
                report_lines.extend(unmatched)
            if still_missing:
                report_lines.append(f"\nStill missing information ({len(still_missing)}):")
                report_lines.extend(f"{sid}\t{', '.join(fields)}" for sid, fields in still_missing)
            with open(report_path, "w", encoding="utf-8") as f:
                f.write("\n".join(report_lines) + "\n")

        msg = f"{os.path.basename(prmag_path)}: {n_updated} specimen(s) updated.\n"
        if unmatched:
            msg += f"No table match ({len(unmatched)}) - see {report_path}\n"
        if still_missing:
            msg += f"Still missing information ({len(still_missing)}) - see {report_path}\n"
        self._afficher(msg)
        self._showinfo(
            "Completion done",
            f"{n_updated} specimen(s) updated in {prmag_path}\n"
            f"(backup: {prmag_path}.bak)"
            + (f"\n{len(unmatched)} specimen(s) had no match in the table." if unmatched else "")
            + (f"\n{len(still_missing)} specimen(s) still have missing information (n.d) "
               "after this completion." if still_missing else "")
            + (f"\nFull list(s) written to {report_path}" if report_path else "")
            + "\n\nRe-open this file to see the changes (the currently loaded data is not refreshed automatically).",
        )

    def ouvrir_check_duplicates_dialog(self):
        """Evalue le .prmag charge et ses compagnons .pmagres/.pmagint/
        .pmagani (memes noms de base, voir results_path_for/pmagint_
        path_for/_find_ani_path) pour des entrees en double (meme
        etape/mêmes codes - voir calcul.find_duplicate_measurements/
        find_duplicate_results/find_duplicate_pmagint_rows/
        find_duplicate_ani_tensors pour la notion exacte de "doublon"
        propre a chaque fichier) - demande explicite utilisateur ("add a
        menu for a routine to evaluate the files prmag, pmagres, pmagint
        and pmagani to check for measurements with the same step and
        codes, just a warning to help the user clean the files"). Chaque
        fichier compagnon absent est simplement saute (pas une erreur) ;
        SIMPLE LISTE D'AVERTISSEMENTS, aucune correction automatique -
        le nettoyage reste a la charge de l'utilisateur."""
        if not self.donnees:
            self._showwarning("No data", "Load a .ren/.prmag file first.")
            return

        self.text_area.insert(tk.END, "\n--- Check files for duplicate entries ---\n", "prompt")
        report = []

        dup_mesures = find_duplicate_measurements(self.donnees)
        report.append(f"[.prmag] {len(dup_mesures)} duplicate measurement(s):" if dup_mesures
                       else "[.prmag] no duplicate measurement found.")
        report.extend(f"  {w}" for w in dup_mesures)

        if self.results_path and os.path.exists(self.results_path):
            file_results = list(_iter_result_lines(self.results_path))
            dup_results = find_duplicate_results(file_results)
            report.append(
                f"[.pmagres] {len(dup_results)} duplicate result(s) in "
                f"{os.path.basename(self.results_path)}:" if dup_results
                else f"[.pmagres] no duplicate result found in {os.path.basename(self.results_path)}."
            )
            report.extend(f"  {w}" for w in dup_results)
        else:
            report.append("[.pmagres] no results file found - skipped.")

        pmagint_path = pmagint_path_for(self.results_path) if self.results_path else None
        if pmagint_path and os.path.exists(pmagint_path):
            dup_pmagint = find_duplicate_pmagint_rows(pmagint_path)
            report.append(
                f"[.pmagint] {len(dup_pmagint)} duplicate interpretation(s):" if dup_pmagint
                else "[.pmagint] no duplicate interpretation found."
            )
            report.extend(f"  {w}" for w in dup_pmagint)
        else:
            report.append("[.pmagint] no paleointensity file found - skipped.")

        ani_path = self._find_ani_path()
        if ani_path:
            dup_ani = find_duplicate_ani_tensors(ani_path)
            report.append(
                f"[.pmagani] {len(dup_ani)} duplicate tensor(s):" if dup_ani
                else "[.pmagani] no duplicate tensor found."
            )
            report.extend(f"  {w}" for w in dup_ani)
        else:
            report.append("[.pmagani] no anisotropy file found - skipped.")

        self._afficher("\n".join(report) + "\n")

    def ouvrir_archive_new_data_dialog(self):
        """Archive de nouvelles mesures (acquises APRES la creation du
        .prmag, sur le meme instrument de labo) dans un .prmag DEJA
        existant - demande explicite utilisateur ("After the creation of
        the .prmag file, it is possible that some new data will be
        acquire on the magnetometer in the lab... I need the possibility
        to archive new data acquired in the legacy files. I need also to
        upload those acquire with the JR6 magnetometer") - voir
        import_new_data.py pour le detail exact du format JR6 et de la
        logique de fusion/deduplication (port de
        reference/ImportJR6/ImportJR6data.f95 et importinpmagren.f).

        Un specimen de la nouvelle acquisition SANS correspondance dans
        le .prmag est ECARTE et journalise, PAS cree - demande explicite
        utilisateur ("we assume that data will be archived only for
        specimens already defined in the .prmag file"), a la difference
        du Fortran d'origine. Patche le .prmag EN PLACE (sauvegarde .bak
        avant toute modification, une seule fois) - accepte plusieurs
        fichiers source a la fois (fusionnes avant l'archivage)."""
        self.text_area.insert(tk.END, "\n--- Archive new lab data (Escape to cancel) ---\n", "prompt")
        source = self._console_input(
            "Source format: legacy Rennes file(s) (r) / JR6 file(s) (j): ", "r")
        if source is None:
            return
        source = (source.strip().lower() or "r")[:1]
        if source not in ("r", "j"):
            self._showerror("Error", "Must be 'r' (legacy Rennes) or 'j' (JR6).")
            return

        prmag_path = filedialog.askopenfilename(
            title="Select the .prmag file to archive new data into",
            filetypes=[("prmag", "*.prmag"), ("All files", "*.*")],
        )
        if not prmag_path:
            return

        if source == "r":
            source_paths = filedialog.askopenfilenames(
                title="Select one or more legacy Rennes file(s) with the new data",
                filetypes=[("REN/Text", "*.ren *.txt"), ("All files", "*.*")],
            )
        else:
            source_paths = filedialog.askopenfilenames(
                title="Select one or more JR6 file(s) with the new data",
                filetypes=[("JR6", "*.jr6 *.txt"), ("All files", "*.*")],
            )
        if not source_paths:
            return

        parser = parse_legacy_new_measurements if source == "r" else parse_jr6_file
        new_by_specimen = {}
        try:
            for path in source_paths:
                for specimen, mesures in parser(path).items():
                    new_by_specimen.setdefault(specimen, []).extend(mesures)
        except Exception as e:
            self._showerror("Error", f"Could not read the source file(s):\n{e}")
            return
        if not new_by_specimen:
            self._showwarning("No data", "No usable measurement found in the selected file(s).")
            return

        try:
            n_specimens, n_added, n_dup, unmatched = archive_new_measurements(prmag_path, new_by_specimen)
        except Exception as e:
            self._showerror("Error", f"Archiving failed:\n{e}")
            return

        msg = (
            f"{os.path.basename(prmag_path)}: {n_added} measurement(s) archived "
            f"across {n_specimens} specimen(s) ({n_dup} already present, skipped).\n"
        )
        if unmatched:
            shown = ", ".join(unmatched[:20])
            more = f", ... ({len(unmatched) - 20} more)" if len(unmatched) > 20 else ""
            msg += (f"NOT archived - specimen not found in this .prmag ({len(unmatched)}): "
                    f"{shown}{more}\n")
        self._afficher(msg)
        self._showinfo(
            "Archiving done",
            f"{n_added} measurement(s) archived in {prmag_path}\n"
            f"(backup: {prmag_path}.bak)"
            + (f"\n{len(unmatched)} specimen(s) had no match - not archived." if unmatched else "")
            + "\n\nRe-open this file to see the changes (the currently loaded data is not refreshed automatically).",
        )

    def ouvrir_convert_magic_to_r_dialog(self):
        """Convertit un fichier de contribution MagIC COMBINE (comme pour
        Import MagIC contribution file..., toutes les tables locations/
        sites/samples/specimens/measurements dans un seul .txt) vers le
        nouveau format .r (voir convert_magic_to_r.py) - reconstruit un
        cod1/cod2 a la Rennes a partir de method_codes/treat_*, detecte
        les protocoles hors perimetre STARpaleomag_Py (AMS, hysteresis, MPMS,
        susceptibilite vs champ/frequence/temperature) et les ecarte sans
        faire echouer la conversion, plutot que les demagnetiser/mesurer
        de paleointensite. Le rapport de fin de conversion (protocoles
        ecartes, method_codes non reconnus) est affiche dans la console
        et dans la boite de dialogue.

        Genere AUSSI le .pmagres compagnon (specimens.txt deja interprete
        + sites.txt deja calcule) - demande explicite utilisateur ("is
        the pmagres file also generated during the magic import") - voir
        convert_magic_to_r._convert_magic_results pour le detail.

        Archive AUSSI, dans le .pmagani compagnon, les tenseurs
        d'anisotropie DEJA calcules par la contribution (colonne
        aniso_s de specimens.txt - A0/TRM, N0/AMS, AA/ARM) - demande
        explicite utilisateur ("lors de l'importation des fichiers
        Magic, archiver dans le fichier pmagani les donnees de tenseurs
        d'anisotropie") - voir extract_magic.magic_anisotropy_rows.

        Utilise DESORMAIS l'orientation (azimuth/dip) au niveau SPECIMEN
        quand la contribution en fournit une (specimens.txt), au lieu de
        toujours retomber sur le niveau sample - demande explicite
        utilisateur ("in some Magic contributions, the orientation of the
        core azimuth and dip are given at the specimen level, when this
        is the case, please take it into account"), voir
        convert_magic_to_r._orientation_source.

        AVANT toute conversion (demande explicite utilisateur "can you
        check the file before the import?"), scan_dip_sign verifie le
        signe des dip effectivement utilises (specimen-level en
        priorite) sur TOUS les specimens du fichier : un signe uniforme
        (tous positifs OU tous negatifs) est affiche en AVERTISSEMENT
        FORT (tag "warn") - "In MagIC, the dip of the core is the angle
        of the X axis from the horizontal" (negation standard vers la
        convention STARpaleomag_Py, cin) - un jeu de donnees reel montre
        normalement un melange de signes ; un signe uniforme peut
        indiquer que cette contribution a deja enregistre le dip dans la
        convention STARpaleomag_Py plutot que la convention MagIC brute.
        Demande alors explicitement s'il faut appliquer la negation
        standard - PAS assume dans un sens ou l'autre."""
        path = filedialog.askopenfilename(
            title="Select the MagIC contribution file (.txt)",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            dip_scan = scan_dip_sign(path)
        except Exception as e:
            self._showerror("Error", f"Could not read {path}:\n{e}")
            return

        negate_dip = True
        dip_warning_note = ""
        if dip_scan["all_positive"] or dip_scan["all_negative"]:
            sign_word = "positive" if dip_scan["all_positive"] else "negative"
            dip_warning_note = (
                f"WARNING: all {dip_scan['n_total']} dip value(s) found in this file are "
                f"{sign_word} ({dip_scan['n_from_specimen']} from specimen-level "
                f"orientation, the rest from sample-level).\n"
                "In MagIC, the dip of the core is the angle of the X axis from the "
                "horizontal - the standard import negates it to match the STARpaleomag_Py "
                "convention (cin). A real drilled data set normally shows a mix of "
                "positive and negative dip values; a uniform sign like this one may mean "
                "this contribution already recorded dip in the STARpaleomag_Py convention "
                "instead of the raw MagIC one.\n"
            )
            # Affiche IMMEDIATEMENT (tag "warn", voir WARN_MARK) pour donner
            # le contexte de la question qui suit - MAIS ce console sera
            # entierement efface par _load_data_file plus bas
            # (self.text_area.delete("1.0", tk.END), le meme mecanisme qui
            # affiche "File loaded: ..."), donc `dip_warning_note` (sans
            # marqueurs, pour rester lisible dans la messagebox finale) est
            # aussi repris dans `msg` ci-dessous - sinon un avertissement
            # "fort" disparaitrait silencieusement de l'ecran une fois la
            # conversion terminee, perdant precisement l'effet recherche.
            self._afficher(f"{WARN_MARK}{dip_warning_note}{WARN_MARK}")
            ans = self._console_input(
                "Import by taking the negative of the dip value (standard MagIC "
                "convention)? Y/n: ", "Y")
            if ans is None:
                return
            negate_dip = ans.strip().lower() != "n"

        base, _ext = os.path.splitext(path)
        output_path = base + ".prmag"
        try:
            nb, report, nb_results, nb_means, nb_pint, redo_pint_path, nb_aniso = convert_magic_file(
                path, output_path, negate_dip=negate_dip)
        except Exception as e:
            self._showerror("Error", f"Conversion failed:\n{e}")
            return
        results_msg = (
            f"Converted: {nb_results} result(s) and {nb_means} site mean(s) -> "
            f"{results_path_for(output_path)}\n"
            if (nb_results or nb_means) else ""
        )
        pint_msg = (
            f"Converted: {nb_pint} paleointensity determination(s) -> {redo_pint_path}\n"
            "(use View Paleoint Results... to replay them)\n"
            if redo_pint_path else ""
        )
        aniso_msg = (
            f"Converted: {nb_aniso} anisotropy tensor(s) -> {ani_path_for(output_path)}\n"
            if nb_aniso else ""
        )
        negate_note = (
            "" if not dip_warning_note
            else f"Dip {'negated (standard MagIC convention)' if negate_dip else 'kept as-is (NOT negated)'} on import.\n"
        )
        msg = (
            f"{dip_warning_note}{negate_note}{report}\n"
            f"Converted: {nb} specimen(s) -> {output_path}\n{results_msg}{pint_msg}{aniso_msg}"
        )
        # Ouvre directement le .prmag converti - demande explicite
        # utilisateur ("ce serait bien d'ouvrir les fichiers convertis a
        # la fin des conversions"). AVANT self._afficher(msg), meme
        # raison que ouvrir_convert_legacy_dialog.
        self._load_data_file(output_path, announce=False)
        self._afficher(msg)
        self._showinfo("Conversion complete", msg)

    def ouvrir_convert_utrecht_to_r_dialog(self):
        """Convertit un ou plusieurs fichiers .col d'Utrecht
        (paleomagnetism.org/PMAG2, JSON malgre l'extension) vers un SEUL
        fichier .prmag combine, PLUS un seul .pmagres compagnon pour les
        interpretations deja calculees par leur logiciel (voir
        convert_utrecht_to_r.py : "import the interpretation in
        Pmagres", "import all the col files in a single .prmag file").
        Site = nom du fichier .col source (sans extension) - ces fichiers
        n'ont pas de notion de site propre, ni extractible du nom de
        specimen ; approximation assumee ("assume that the site is the
        name of the col file"), a raffiner manuellement au besoin.
        Orientation de carotte (azimuth/dip) convertie via la regle
        donnee par l'utilisateur (caz=coreAzimuth+90, cin=90-coreDip) et
        verifiee numeriquement sur les 156 interpretations reelles
        fournies (voir docstring de convert_utrecht_to_r.py) - les vues
        in-situ/apres pendage sont fiables pour ces echantillons."""
        paths = filedialog.askopenfilenames(
            title="Select one or more Utrecht/PMAG2 .col files (JSON)",
            filetypes=[("Utrecht collection", "*.col"), ("All files", "*.*")],
        )
        if not paths:
            return
        default_name = (
            os.path.splitext(os.path.basename(paths[0]))[0] + ".prmag" if len(paths) == 1
            else "utrecht_combined.prmag"
        )
        output_path = filedialog.asksaveasfilename(
            title="Save combined .prmag as",
            defaultextension=".prmag",
            initialfile=default_name,
            initialdir=os.path.dirname(paths[0]),
            filetypes=[("STARpaleomag_Py prmag", "*.prmag"), ("All files", "*.*")],
        )
        if not output_path:
            return
        try:
            nb, nb_results = convert_utrecht_files(list(paths), output_path)
        except Exception as e:
            self._showerror("Error", f"Conversion failed:\n{e}")
            return
        sites = ", ".join(sorted({os.path.splitext(os.path.basename(p))[0] for p in paths}))
        msg = (
            f"Converted: {nb} specimen(s) from {len(paths)} file(s) -> {output_path}\n"
            f"Converted: {nb_results} interpretation(s) -> {results_path_for(output_path)}\n"
            f"\nSite = source .col file name: {sites}\n"
            "Refine manually in STARpaleomag_Py if a file bundles more than one real site.\n"
        )
        # Ouvre directement le .prmag converti - demande explicite
        # utilisateur ("ce serait bien d'ouvrir les fichiers convertis a
        # la fin des conversions"). AVANT self._afficher(msg), meme
        # raison que ouvrir_convert_legacy_dialog.
        self._load_data_file(output_path, announce=False)
        self._afficher(msg)
        self._showinfo("Conversion complete", msg)

    def ouvrir_convert_ipgp_to_r_dialog(self):
        """Convertit les fichiers de mesures brutes IPGP/CryoMag(.pmd)
        (developpes par Johan Guyodo/Dragomir Dragomirov a l'IPGP, et le
        format CryoMag/PaleoMac "pour les fichiers .pmd") vers
        .prmag/.pmagres - demande explicite utilisateur. Deux fichiers
        INDEPENDANTS, chacun optionnel a l'appel (au moins un requis) :
        le fichier de mesures brutes -> .prmag (2 variantes reconnues,
        voir convert_ipgp_to_r.py : paliers nus avec ligne de site, ou
        paliers prefixes T.../A... type CryoMag/.pmd sans ligne de
        site), et le fichier "Results" CSV IPGP (interpretations deja
        calculees) -> .pmagres. Transformation d'orientation
        (caz=a+90/cin=b/bed_dip_strike=s/bed_dip=d) VERIFIEE
        numeriquement sur les DEUX variantes (sources differentes,
        Dg/Ig/Ds/Is reproduits, voir docstring module)."""
        measurements_path = filedialog.askopenfilename(
            title="Select IPGP raw measurement file (Cancel to skip)",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
        )
        results_path = filedialog.askopenfilename(
            title="Select IPGP Results CSV file (Cancel to skip)",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not measurements_path and not results_path:
            return
        default_base = measurements_path or results_path
        output_path = filedialog.asksaveasfilename(
            title="Save converted .prmag as",
            defaultextension=".prmag",
            initialfile=os.path.splitext(os.path.basename(default_base))[0] + ".prmag",
            initialdir=os.path.dirname(default_base),
            filetypes=[("STARpaleomag_Py prmag", "*.prmag"), ("All files", "*.*")],
        )
        if not output_path:
            return
        try:
            nb_sp, nb_res, nb_dup = convert_ipgp_files(
                measurements_path or None, results_path or None, output_path
            )
        except Exception as e:
            self._showerror("Error", f"Conversion failed:\n{e}")
            return
        msg = (
            f"Converted: {nb_sp} specimen(s) -> {output_path}\n"
            f"Converted: {nb_res} interpretation(s) -> {results_path_for(output_path)}"
            f"{f' ({nb_dup} exact duplicate row(s) skipped)' if nb_dup else ''}\n"
            "\nAssumptions used (format examples were flagged incomplete - verify and "
            "report if wrong):\n"
            "- site = the single-token line just before a specimen header, if any "
            "(none found -> \"n.d\")\n"
            "- step type: read directly from a letter prefix when present "
            "(T...=thermal, A...=AF); for a bare numeric step (no prefix), AF is "
            "assumed if most non-zero steps are < 150\n"
        )
        self._load_data_file(output_path, announce=False)
        self._afficher(msg)
        self._showinfo("Conversion complete", msg)

    def ouvrir_convert_montpellier_to_r_dialog(self):
        """Convertit un fichier ".montpellier" (table plate Sample/Step/
        Intensity/Dec/Inc/Azimuth/Plunge/Strike/Dip) vers .prmag - demande
        explicite utilisateur ("convert also the format from Montpellier").
        Voir convert_montpellier_to_r.py pour le detail : colonnes
        VERIFIEES (pas supposees) contre BN.col (meme collection BN4.*,
        deja validee pour convert_utrecht_to_r.py) - Dec/Inc = repere
        specimen brut, Intensity = magnetisation A/m deja, Azimuth/Plunge
        suivent la convention PMAG2/Utrecht (caz=Azimuth+90,
        cin=90-Plunge), PAS celle de convert_ipgp_to_r.py. Ce format ne
        porte AUCUN volume : demande a l'utilisateur (defaut 10.80 cm3,
        meme defaut que ams_prmag.create_prmag_from_legacy_ani) - la
        direction/les ajustements n'en dependent pas, seul MAG(A/m)
        affiche serait faux si le volume reel differe."""
        path_in = filedialog.askopenfilename(
            title="Select a Montpellier .montpellier file",
            filetypes=[("Montpellier", "*.montpellier"), ("All files", "*.*")],
        )
        if not path_in:
            return
        vol_s = self._console_input(
            "Assumed specimen volume in cm3 (this format does not carry one): ",
            f"{_MONTPELLIER_DEFAULT_VOLUME_CM3:.2f}",
        )
        if vol_s is None:
            return
        try:
            volume_cm3 = float(vol_s)
        except ValueError:
            self._showerror("Error", f"Not a number: {vol_s}")
            return
        output_path = filedialog.asksaveasfilename(
            title="Save converted .prmag as",
            defaultextension=".prmag",
            initialfile=os.path.splitext(os.path.basename(path_in))[0] + ".prmag",
            initialdir=os.path.dirname(path_in),
            filetypes=[("STARpaleomag_Py prmag", "*.prmag"), ("All files", "*.*")],
        )
        if not output_path:
            return
        try:
            nb_sp = convert_montpellier_file(path_in, output_path, volume_cm3)
        except Exception as e:
            self._showerror("Error", f"Conversion failed:\n{e}")
            return
        msg = (
            f"Converted: {nb_sp} specimen(s) -> {output_path}\n"
            f"\nAssumed specimen volume: {volume_cm3:.2f} cm3 (not carried by this "
            "format - correct manually per specimen if the real volume differs; "
            "direction/fits are unaffected, only the displayed MAG(A/m)).\n"
            "site: \"n.d\" (this format has no site column either).\n"
        )
        self._load_data_file(output_path, announce=False)
        self._afficher(msg)
        self._showinfo("Conversion complete", msg)

    def ouvrir_field_notes_dialog(self):
        """Equivalent GUI de `subroutine orientation` (orient_paleomag.f -
        demande explicite utilisateur "ajouter le menu field notes to
        Pmag files. From the field orientations, it prepares the prmag
        files (check the Fortran file). It also prepares the .ged file
        for the AGICO instruments"). Voir field_notes.py pour le detail du
        format d'entree (fichier "simple", REVU a la demande explicite de
        l'utilisateur pour etre plus rapide a remplir sur le terrain - "the
        field notes are more like the original file. as we now have the
        possibility to complete the different fields later, we can have
        two files (a simple file written following the field trip) and a
        complement file that can be written later") et des sorties
        produites :

        - un .prmag (format MODERNE STARpaleomag_Py, mesures vides - un
          instrument n'a pas encore mesure ces specimens a ce stade) ;
        - UN SEUL .ged combinant tous les sites pour les instruments AGICO
          (meme comportement que le Fortran, `unfichierged` code en dur a
          'Y' - demande explicite utilisateur) ;
        - un rapport diagnostic (azimut geo/declinaison locale/IGRF par
          specimen, "err>3/5/10" si l'ecart depasse un seuil - complete le
          calcul manifestement voulu mais jamais assigne dans le Fortran,
          `d_err`, variable non initialisee dans le source).

        `tool`/`ellipsoid` ne sont plus dans le fichier "simple" lui-meme
        (ne variaient pas d'un site a l'autre dans les fichiers reels vus
        jusqu'ici) - demandes UNE FOIS ici, appliques a tous les sites."""
        self.text_area.insert(tk.END, "\n--- Field notes to Pmag files (Escape to cancel) ---\n", "prompt")
        ori_path = filedialog.askopenfilename(
            title="Select the field notes file",
            filetypes=[("Field notes", "*.txt *.ori"), ("All files", "*.*")],
        )
        if not ori_path:
            return

        tool = self._console_input(
            "Orientation tool - ASC (vertical needle, rotating stage) or "
            "other (pivoting square, fixed stage): ", "ASC")
        if tool is None:
            return
        ellipsoid = self._console_input(
            "Ellipsoid/datum (WGS84 or other, e.g. SA56): ", "WGS84")
        if ellipsoid is None:
            return
        strati = self._console_input(
            "Bedding convention - Strikedip (last column already is the "
            "strike) or Dipdip (last column is the dip direction): ", "Strikedip")
        if strati is None:
            return
        ispec_s = self._console_input("Specimens per core (1, 2 or 3): ", "2")
        if ispec_s is None:
            return
        try:
            ispec = int(ispec_s)
        except ValueError:
            self._showerror("Error", "Must be an integer (1, 2 or 3).")
            return
        if ispec not in (1, 2, 3):
            self._showerror("Error", "Must be 1, 2 or 3.")
            return

        try:
            with open(ori_path, "r", encoding="utf-8", errors="replace") as f:
                lines = f.readlines()
        except Exception as e:
            self._showerror("Error", f"Could not read the file:\n{e}")
            return

        sites, warnings = parse_orientation_file(
            lines, tool=tool or "ASC", ellipsoid=ellipsoid or "WGS84", strati=strati or "Strikedip")
        if not sites:
            self._showwarning(
                "No site found",
                "No usable site/core data found in this file (see the console "
                "for details if warnings were printed).",
            )
            for w in warnings:
                print(f"WARNING: {w}")
            return

        # Fichier complement OPTIONNEL (geologie + carottes ajoutees apres
        # coup) - meme principe que le complement de
        # ouvrir_convert_legacy_dialog - demande explicite utilisateur.
        if messagebox.askyesno(
            "Complement file",
            "Provide a complement file (site geology - Formation/Age/... - "
            "and/or additional cores added after the field trip)?",
        ):
            complement_path = filedialog.askopenfilename(
                title="Select the complement file",
                filetypes=[("Text", "*.txt"), ("All files", "*.*")],
            )
            if complement_path:
                try:
                    with open(complement_path, "r", encoding="utf-8", errors="replace") as f:
                        clines = f.readlines()
                except Exception as e:
                    self._showerror("Error", f"Could not read the complement file:\n{e}")
                    return
                complement, cwarnings = parse_complement_file(clines)
                warnings = warnings + cwarnings + apply_complement(sites, complement)

        # Verifie que le utdif enregistre est bien celui qui minimise
        # l'ecart declinaison locale/IGRF - demande explicite utilisateur
        # ("one of the most common error is with the UTM dif...").
        warnings = warnings + check_utdif_for_sites(sites)

        base, _ext = os.path.splitext(ori_path)
        prmag_path = base + "_pmag.prmag"
        ged_path = base + ".ged"
        res_path = base + "_diagnostics.txt"

        records = build_pmag_records(sites, ispec)
        n_prmag = write_prmag_from_field_notes(records, prmag_path)
        n_ged = write_ged_file(sites, ispec, ged_path)
        n_diag = write_diagnostics_report(sites, res_path)

        msg = (
            f"{len(sites)} site(s), {n_prmag} specimen(s) -> {prmag_path}\n"
            f"{n_ged} specimen(s) -> {ged_path} (AGICO, all sites combined)\n"
            f"Diagnostics ({n_diag} core(s), az geo / local declination / IGRF) -> {res_path}\n"
        )
        if warnings:
            msg += f"\n{len(warnings)} warning(s):\n"
            for w in warnings[:20]:
                msg += f"  - {w}\n"
            if len(warnings) > 20:
                msg += f"  ... and {len(warnings) - 20} more (see console).\n"
                for w in warnings:
                    print(f"WARNING: {w}")

        self._load_data_file(prmag_path, announce=False)
        self._afficher(msg)
        self._showinfo("Field notes converted", msg)

    def ouvrir_ged_to_prmag_dialog(self):
        """Cree un .prmag (mesures VIDES) a partir d'un fichier .ged AGICO
        seul - demande explicite utilisateur ("pourrait on mettre un menu
        de creation de prmag file a partir d'un fichier agico .ged. et
        ensuite archive data") : couvre le cas ou seules les field notes
        BRUTES (le fichier "simple" attendu par `Orientation: field notes
        orientation to prmag file...`) ne sont pas disponibles, mais un
        .ged deja produit (par ce meme menu, ou recu d'un collegue, ou
        exporte par le logiciel de l'instrument AGICO) l'est. Voir
        field_notes.parse_ged_file pour le detail de ce que le .ged porte
        (specimen + geometrie SEULEMENT - pas le site/la date/la
        geologie). Le .prmag cree n'a pas de mesures : l'etape suivante
        naturelle est `Archive new laboratory measurements` pour y
        attacher les mesures reelles de l'instrument (qui exige des
        specimens DEJA presents dans le fichier cible - exactement ce que
        ce menu prepare)."""
        self.text_area.insert(tk.END, "\n--- Create prmag from AGICO .ged file (Escape to cancel) ---\n", "prompt")
        ged_path = filedialog.askopenfilename(
            title="Select the AGICO .ged file",
            filetypes=[("AGICO .ged", "*.ged"), ("All files", "*.*")],
        )
        if not ged_path:
            return

        try:
            records = parse_ged_file(ged_path)
        except Exception as e:
            self._showerror("Error", f"Could not read the .ged file:\n{e}")
            return
        if not records:
            self._showwarning("No specimen found", "No usable specimen line found in this .ged file.")
            return

        base, _ext = os.path.splitext(ged_path)
        prmag_path = base + "_pmag.prmag"
        n_prmag = write_prmag_from_field_notes(records, prmag_path)

        msg = (
            f"{n_prmag} specimen(s) -> {prmag_path}\n"
            "Only specimen id, core azimuth/dip and bedding strike/dip come "
            "from the .ged - site, date, geology, volume/mass are left as "
            "'n.d'/defaults (fill in later with Complete sample "
            "information...).\n"
            "Next step: PmagFile -> Archive new laboratory measurements, to "
            "attach the instrument's actual measurements to these specimens.\n"
        )
        self._load_data_file(prmag_path, announce=False)
        self._afficher(msg)
        self._showinfo("prmag created from .ged", msg)

    def ouvrir_export_magic_dialog(self):
        """Equivalent GUI de `export2magic` ("export Rennes to Magic",
        fichiers_magic.f) - mode classique (sites/samples/specimens/
        measurements.txt + locations.txt) PLUS, si un .pmagint existe,
        les resultats de paleointensite deja archives dedans (voir
        magic_export.py pour le detail des ecarts assumes par rapport au
        Fortran : tri explicite par site, les 2 bugs Fortran reperes sont
        corriges, formatage numerique simple ; magnetostratigraphie
        toujours HORS PERIMETRE). Les champs Site/Sample/Fm/Age/GC/SMT/
        Li/Loc utilises viennent de la ligne « roche » decodee a
        l'ouverture du fichier (testlect.decode_roche) - un echantillon
        sans site MagIC decode n'aura simplement pas de ligne dans
        sites.txt.

        Deux ameliorations demandees explicitement par l'utilisateur
        ("il faudrait guider un peu l'utilisateur en l'invitant a
        selectionner les donnees resultats, paleointensite, anisotropie
        [des le debut]... avertir les utilisateurs sur le contenu des
        caracteristiques des echantillons... exporter d'abord un tableau
        avec site, lat lon Geologic classes, types lithology... plutot
        que de poser cette question a la fin") :

        1) Un bilan de disponibilite (nb de resultats/paleointensite/
           anisotropie trouves pour la selection) est affiche EN
           PREMIER, suivi de TROIS choix consolides (inclure
           resultats/paleointensite/anisotropie - Y/n chacun) - remplace
           les anciens prompts Y/n eparpilles plus loin dans le
           dialogue (un par categorie, au moment ou chacune etait
           decouverte), qui ne donnaient a l'utilisateur aucune vue
           d'ensemble avant de decider.
        2) Une table d'apercu des metadonnees de site (site/lat/lon/
           geologic_classes/geologic_types/lithologies/formation, voir
           magic_export.build_site_metadata_preview_rows) peut etre
           exportee pour revue - suivie d'une confirmation "Are the
           specimen metadata OK" et, seulement si la reponse est non,
           d'un prompt pour completer via un fichier complement. Ce bloc
           entier a d'abord ete place juste apres le bilan de
           disponibilite (1), PUIS deplace en TOUTE DERNIERE etape
           interactive - demande explicite utilisateur, une fois le
           premier comportement teste en pratique ("can you ask these
           questions at the end just before writing the files") : il
           s'execute desormais juste avant l'appel a export_to_magic,
           apres le choix du dossier de sortie - plus aucune autre
           decision a prendre entre cette question et l'ecriture reelle
           des fichiers. Le fichier ecrit (write_site_metadata_preview)
           utilise EXACTEMENT les memes noms de colonnes que
           load_site_metadata_table : on peut le completer a la main
           puis le redonner tel quel au prompt du fichier complement."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        self.text_area.insert(tk.END, "\n--- export Rennes to Magic (Escape to cancel) ---\n", "prompt")

        # --- Bilan de disponibilite + choix consolides (voir docstring,
        # point 1) : calcule ce qui EXISTE pour la selection courante
        # avant de poser la moindre question, pour que l'utilisateur
        # sache d'emblee ce qu'il peut inclure. `n_results`/
        # `pmagint_available` portent sur TOUT self.results/.pmagint
        # (pas seulement l'intersection avec self.selection) UNIQUEMENT
        # pour que ce bilan annonce le total reel avant filtrage - un
        # resultat/une paleointensite dont le specimen n'a PAS de mesure
        # de desaimantation/paleointensite dans la selection courante
        # n'est PAS inclus dans l'export (voir plus bas : signale par
        # une liste, jamais archive) - demande explicite utilisateur,
        # apres une premiere tentative erronee d'auto-inclusion (voir git
        # history) : "on ne devrait pas avoir de resultats de direction
        # ou d'intensite si il n'y a eu de desaimantation ou de manip de
        # paleointensite. C'est different pour l'AMS car on peut faire
        # une mesure d'AMS sur des specimens non traites en remanence."
        n_results = len({r.id for r in self.results if r.cat1 in ("L", "P", "f", "F")})
        pmagint_path = pmagint_path_for(self.results_path) if self.results_path else None
        pmagint_available = (
            dict(read_pmagint(pmagint_path)) if pmagint_path and os.path.exists(pmagint_path) else {}
        )
        ani_path = self._find_ani_path()
        n_aniso_raw = sum(1 for ech in self.selection if any(m.cod1 in ("X", "Y", "Z") for m in ech.mesures))

        self._afficher(
            f"{len(self.selection)} specimen(s) selected for export. Available data:\n"
            f"  - directional results (.pmagres): {n_results} specimen(s)\n"
            f"  - paleointensity (.pmagint): {len(pmagint_available)} specimen(s)\n"
            f"  - anisotropy: {n_aniso_raw} specimen(s) with raw X/Y/Z measurements"
            + (f", archived tensor(s) found in {os.path.basename(ani_path)}" if ani_path else "")
            + "\n"
            "(a result/paleointensity for a specimen with no demagnetization "
            "measurement in the current selection is reported below, never "
            "included in the export.)\n"
        )

        want_results = False
        if n_results:
            ans = self._console_input(f"Include directional results ({n_results} specimen(s))? Y/n: ", "Y")
            if ans is None:
                return
            want_results = ans.strip().lower() != "n"

        want_pmagint = False
        if pmagint_available:
            ans = self._console_input(
                f"Include paleointensity results for {len(pmagint_available)} specimen(s) "
                f"(from {os.path.basename(pmagint_path)})? Y/n: ", "Y")
            if ans is None:
                return
            want_pmagint = ans.strip().lower() != "n"

        want_aniso = False
        if n_aniso_raw or ani_path:
            ans = self._console_input(
                "Include anisotropy data (raw X/Y/Z measurements and/or archived .pmagani "
                "tensor)? Y/n: ", "Y")
            if ans is None:
                return
            want_aniso = ans.strip().lower() != "n"
        # --- fin bilan/choix ---

        lab_analysts = self._console_input(
            "Analysts' names (in quotes, separated by ':'): ", "")
        if lab_analysts is None:
            return

        continent = self._console_input("Continent/ocean (optional): ", "")
        if continent is None:
            return
        country = self._console_input("Country (optional): ", "")
        if country is None:
            return
        region = self._console_input("Region (optional): ", "")
        if region is None:
            return

        # Anisotropie X/Y/Z (LP-AN-TRM vs LP-AN-IRM) - demande explicite
        # utilisateur ("anisotropy with code X+,Y+,Z+ X-,Y-,Z- is usually
        # done by TRM acquisition ... but there are some done with high
        # field IRM ... ask the user to confirm the kind of anisotropy
        # experiment when it is dubious? and whether it wants to archive
        # these data"). classify_anisotropy_experiment tranche seule dans
        # les cas non-ambigus (companion paleointensite -> trm ; etape >
        # 700 -> irm) ; sinon on demande specimen par specimen. Rien de
        # tout ca si `want_aniso` est deja Non (choix consolide ci-dessus).
        anisotropy_kind_by_specimen = {}
        anisotropy_skip = set()
        if want_aniso:
            for ech in self.selection:
                if not any(m.cod1 in ("X", "Y", "Z") for m in ech.mesures):
                    continue
                kind = classify_anisotropy_experiment(ech.mesures)
                if kind is not None:
                    anisotropy_kind_by_specimen[ech.id] = kind
                    continue
                archive = self._console_input(
                    f"{ech.id}: anisotropy experiment kind is unclear. "
                    f"Archive this anisotropy data (y/N): ", "N")
                if archive is None:
                    return
                if archive.strip().lower() != "y":
                    anisotropy_skip.add(ech.id)
                    continue
                kind_s = self._console_input(
                    f"{ech.id}: TRM acquisition or high-field IRM anisotropy "
                    f"(trm/irm): ", "trm")
                if kind_s is None:
                    return
                anisotropy_kind_by_specimen[ech.id] = (
                    "irm" if kind_s.strip().lower().startswith("i") else "trm")

        # Paleointensite (.pmagint, deja archivee via "Paleointensity
        # interpretation"/"View batch of Paleoint Results..."/"Rapid
        # view/recompute from redo file...") - demande explicite
        # utilisateur ("in export to Magic, it is not asking for
        # Paleointensity") : ce dialogue ne proposait jusqu'ici que les
        # donnees directionnelles/anisotropie, jamais les resultats de
        # paleointensite deja calcules et sauvegardes. Inclusion decidee
        # par `want_pmagint` ci-dessus (bilan/choix consolide), plus de
        # second prompt ici.
        pmagint_rows = pmagint_available if want_pmagint else {}

        # Anisotropie (.pmagani, tenseur natif deja calcule) - demande
        # explicite utilisateur ("il manque aussi les données
        # d'anisotropie au niveau du fichier specimens (prendre les A0)")
        # : specimens.txt n'avait jusqu'ici AUCUNE colonne aniso_* (les
        # prompts anisotropy_kind_by_specimen/anisotropy_skip ci-dessus ne
        # servent qu'a etiqueter les mesures BRUTES X/Y/Z de
        # measurements.txt, pas le tenseur RESULTAT). Cherche, PAR ORDRE
        # DE PRIORITE, un tenseur A0 (ATRM) puis F0 (AARM) puis N0 (AMS
        # basse frequence) - PAS seulement 'A0' comme le premier essai :
        # un .pmagani importe depuis AMS_Py/un ASC Agico (ex.
        # rn01_15_converted.pmagani reel) est presque toujours code2='N0'
        # (susceptibilite), jamais 'A0' - demande explicite utilisateur
        # en verifiant "how to join the two exports from prmag & pmagres
        # with the Anisotropy" sur un vrai fichier ("oui dans le dossier
        # Concon_creixell"). Tout ce bloc (tenseurs specimen/site ET
        # l'extension export_samples pour les specimens AMS-only) est
        # desormais gate par `want_aniso` (choix consolide ci-dessus) -
        # plus de prompts Y/n separes par sous-categorie.
        _ANISO_CODE_PRIORITY = ("A0", "F0", "N0")
        aniso_tensors = {}
        aniso_mean_tensors = {}
        export_samples = list(self.selection)

        # Specimen avec un resultat/une paleointensite archives mais SANS
        # mesure de desaimantation/paleointensite dans la selection
        # courante - demande explicite utilisateur, REVENANT sur la
        # tentative precedente (voir git history) qui auto-ajoutait ces
        # specimens comme pour l'AMS : "l'approche precedente etait
        # correcte, on ne devrait pas avoir de resultats de direction ou
        # d'intensite si il n'y a eu de desaimantation ou de manip de
        # paleointensite. C'est different pour l'AMS car on peut faire
        # une mesure d'AMS sur des specimens non traites en remanence."
        # Un resultat/une paleointensite pour un tel specimen signale
        # donc une VRAIE incoherence (donnees desassorties, specimen
        # supprime du .prmag apres coup, faute de frappe...), pas une
        # simple selection trop etroite a completer automatiquement :
        # SIGNALE (liste des id concernes) mais jamais archive - PAS
        # d'ajout a `export_samples` ici (contrairement au bloc AMS-only
        # plus bas, legitime lui : une mesure AMS n'exige aucun
        # traitement de remanence prealable).
        known_ids = {e.id.strip().upper() for e in export_samples}
        if want_results:
            inconsistent_ids = sorted({
                r.id.strip() for r in self.results
                if r.id.strip() and r.id.strip().upper() not in known_ids
                and r.cat1 in ("L", "P", "f", "F")
            })
            if inconsistent_ids:
                shown = ", ".join(inconsistent_ids[:20])
                more = f", ... ({len(inconsistent_ids) - 20} more)" if len(inconsistent_ids) > 20 else ""
                self._afficher(
                    f"WARNING: {len(inconsistent_ids)} specimen(s) have an archived "
                    f"directional result but no demagnetization measurement in the "
                    f"current selection - NOT included in this export (check the "
                    f".pmagres file / the sample selection for an inconsistency): "
                    f"{shown}{more}\n"
                )

        if want_pmagint:
            inconsistent_ids = sorted({
                sid.strip() for sid in pmagint_available
                if sid.strip() and sid.strip().upper() not in known_ids
            })
            if inconsistent_ids:
                shown = ", ".join(inconsistent_ids[:20])
                more = f", ... ({len(inconsistent_ids) - 20} more)" if len(inconsistent_ids) > 20 else ""
                self._afficher(
                    f"WARNING: {len(inconsistent_ids)} specimen(s) have archived "
                    f"paleointensity but no measurement in the current selection - "
                    f"NOT included in this export (check the .pmagint file / the "
                    f"sample selection for an inconsistency): {shown}{more}\n"
                )

        if want_aniso and ani_path:
            # Specimens connus SEULEMENT via leur tenseur AMS (.pmagani),
            # pas dans self.selection (jamais dans le .prmag, ou juste
            # filtres par la selection courante) - demande explicite
            # utilisateur ("attention: l'exportation concerne seulement
            # les spécimens/samples avec des données prmag par contre il
            # peut y avoir des specimens dans pmagani qui n'ont pas de
            # mesures dans prmag ... mais il faut aussi les exporter").
            # `export_samples` (self.selection + ces extras) remplace
            # self.selection pour TOUT le reste de cette methode
            # (recherche des tenseurs, des sites, et l'appel final a
            # export_to_magic) - self.selection lui-meme reste inchange.
            # Cherche DIRECTEMENT dans self.donnees par id exact (PAS via
            # select_samples : verifie sur donnees reelles - des
            # specimens AMS-only d'un vrai .pmagani sont bien dans le
            # .prmag mais SANS aucune mesure de demag, un cas que
            # select_samples exclut TOUJOURS sauf id vide, voir sa
            # docstring - "AMS sans jamais de vraie mesure NRM" est le
            # meme genre de cas que "site sans specimen", ici au niveau
            # specimen) - vraies metadonnees site/geologie alors
            # recuperees via _build_selected_sample.
            #
            # Un specimen du .pmagani INTROUVABLE dans self.donnees
            # (aucune trace, meme sans mesure) est REJETE, pas fabrique -
            # demande explicite utilisateur ("en theorie pmagani est lie
            # a prmag. Il faut interdire dans pmagani les specimens qui
            # ne sont pas présents dans prmag") : un .pmagani est cense
            # decrire des specimens du .prmag compagnon ; un id absent du
            # .prmag signale un vrai probleme (fichiers desassortis,
            # faute de frappe, .pmagani d'un autre jeu de donnees) que
            # fabriquer une ligne specimens.txt minimaliste masquerait
            # plutot que resoudre - SIGNALE et EXCLU de l'export, jamais
            # silencieusement invente.
            known_ids = {e.id.strip().upper() for e in export_samples}
            # `t.export == "N"` (voir calcul.AniTensor.export, colonne
            # ajoutee cote AMS_Py par "Mark selection for MagIC
            # export...") exclut ce tenseur - un specimen connu
            # SEULEMENT via un tel tenseur n'a donc plus aucune raison
            # d'etre ajoute ici (demande explicite utilisateur "only
            # these selected data will be taken into account in the
            # main export from Starpaleomag").
            extra_ids = sorted({
                t.id.strip() for t in read_all_ani_tensors(ani_path)
                if t.id.strip() and t.id.strip().upper() not in known_ids
                and t.export != "N"
            })
            donnees_by_id = {p.id.strip().upper(): p for p in self.donnees}
            n_from_prmag = 0
            orphan_ids = []
            for extra_id in extra_ids:
                p = donnees_by_id.get(extra_id.upper())
                if p is not None:
                    export_samples.append(_build_selected_sample(p, p.mesures))
                    n_from_prmag += 1
                else:
                    orphan_ids.append(extra_id)
            if n_from_prmag:
                self._afficher(
                    f"[Export to MagIC] {n_from_prmag} specimen(s) with AMS data only "
                    f"(no measurement in the current selection) added for export.\n"
                )
            if orphan_ids:
                shown = ", ".join(orphan_ids[:20])
                more = f", ... ({len(orphan_ids) - 20} more)" if len(orphan_ids) > 20 else ""
                self._afficher(
                    f"WARNING: {len(orphan_ids)} specimen(s) in {os.path.basename(ani_path)} "
                    f"have NO counterpart in the loaded .prmag - excluded from export "
                    f"(check for mismatched files or a typo): {shown}{more}\n"
                )

            for ech in export_samples:
                for code2 in _ANISO_CODE_PRIORITY:
                    tensor = read_ani_tensor(ani_path, ech.id, code2)
                    # export == "N" : specimen explicitement ECARTE de
                    # l'export MagIC via "Mark selection for MagIC
                    # export..." cote AMS_Py - continue vers le code2
                    # suivant plutot que de considerer qu'il n'y a
                    # simplement aucun tenseur A0/F0/N0 (un autre code2
                    # du meme specimen peut, lui, etre marque export=Y).
                    if tensor is not None and tensor.export == "N":
                        continue
                    if tensor is not None:
                        aniso_tensors[ech.id] = tensor
                        break

            # Anisotropie de SITE (.pmagani, section "#site mean tensor
            # results", moyenne deja calculee - voir calcul.AniMeanTensor/
            # read_ani_mean_tensor) - demande explicite utilisateur ("je
            # voudrais ajouter l'exportation de l'AMS (niveau specimen et
            # sites). Les données d'AMS vont aussi dans les fichiers
            # sites.txt et specimens.txt") : meme principe que le bloc
            # specimen ci-dessus (meme ordre de priorite A0/F0/N0, meme
            # raison reelle - le .pmagani reel teste est en 'N0'), une
            # ligne PAR SITE distinct de la selection plutot que par
            # specimen.
            for site in {ech.magic_site.strip() for ech in export_samples if ech.magic_site.strip()}:
                for code2 in _ANISO_CODE_PRIORITY:
                    mean_tensor = read_ani_mean_tensor(ani_path, site, code2)
                    # Meme filtre export=="N" que pour les tenseurs
                    # specimen ci-dessus (voir calcul.AniMeanTensor.export).
                    if mean_tensor is not None and mean_tensor.export == "N":
                        continue
                    if mean_tensor is not None:
                        aniso_mean_tensors[site] = mean_tensor
                        break

        out_dir = filedialog.askdirectory(title="MagIC output folder (sites/samples/specimens/measurements.txt)")
        if not out_dir:
            return

        # Table d'apercu des metadonnees de site (voir docstring, point
        # 2) - demande explicite utilisateur ("can you ask these
        # questions at the end just before writing the files") : ce bloc
        # entier (bilan, export de la table, confirmation, fichier
        # complement) est desormais la DERNIERE etape interactive, juste
        # avant l'appel a export_to_magic ci-dessous - plus rien d'autre
        # a decider entre cette question et l'ecriture reelle des
        # fichiers.
        preview_rows = build_site_metadata_preview_rows(export_samples)
        missing_sites = [
            r for r in preview_rows
            if not r["lat"] or not r["geologic_classes"] or not r["geologic_types"] or not r["lithologies"]
        ]
        preview_path = None
        if preview_rows:
            self._afficher(
                "\nMagIC requires site-level metadata (coordinates + geologic "
                "classification) for a valid contribution.\n"
                f"Site metadata check: {len(preview_rows)} site(s) in this export, "
                f"{len(missing_sites)} with at least one missing field\n"
                "(lat/lon, geologic classes, geologic types, or lithology).\n"
            )
            export_preview = self._console_input(
                "Export a site metadata table (site/lat/lon/geologic_classes/"
                "geologic_types/lithologies/formation)\n"
                "to review/complete before continuing? Y/n: ", "Y" if missing_sites else "n")
            if export_preview is None:
                return
            if export_preview.strip().lower() != "n":
                preview_path = filedialog.asksaveasfilename(
                    title="Save site metadata table as",
                    defaultextension=".txt",
                    initialfile="site_metadata_preview.txt",
                    # Dossier du .prmag charge par defaut - demande
                    # explicite utilisateur ("where is it saved?") : sans
                    # initialdir, la boite de dialogue s'ouvrait sur le
                    # dernier dossier utilise par macOS (souvent sans
                    # rapport), ne laissant aucun indice visible d'ou le
                    # fichier finirait par etre ecrit.
                    initialdir=os.path.dirname(self.results_path) if self.results_path else None,
                    filetypes=[("Text", "*.txt"), ("All files", "*.*")],
                )
                if preview_path:
                    write_site_metadata_preview(preview_path, preview_rows)
                    self._afficher(
                        f"Site metadata table written to {preview_path}.\n"
                        "Edit it now if needed, then use it as the complement table "
                        "below - the column names already match.\n"
                    )

        # Confirmation explicite AVANT de proposer un fichier complement -
        # demande explicite utilisateur : plutot que de toujours poser la
        # question "Fill in... from a complement table", demande d'abord
        # si les metadonnees sont bonnes telles quelles (l'utilisateur a
        # pu les completer a la main entre-temps, dans le fichier
        # exporte ci-dessus ou directement dans le .prmag) - la question
        # du fichier complement n'est plus posee du tout si la reponse
        # est oui. Defaut "n" (pas OK) si des champs manquaient au bilan
        # ci-dessus, "Y" sinon.
        metadata_ok = True
        if preview_rows:
            ok_answer = self._console_input(
                "Are the specimen metadata OK (Y/n): ", "n" if missing_sites else "Y")
            if ok_answer is None:
                return
            metadata_ok = ok_answer.strip().lower() != "n"

        # Metadonnees de site (formation/lithologies/geologic_classes/
        # geologic_types/age) depuis une table externe - demande
        # explicite utilisateur ("can we also let the site-only path
        # pull from a complement table... so those fields aren't just
        # blank") : ne COMBLE que les champs encore vides (voir
        # magic_export._apply_site_metadata) - utile en particulier pour
        # un site archive uniquement via sa moyenne (aucun specimen
        # charge, voir build_sites_rows), mais s'applique aussi aux
        # sites avec specimens dont un champ n'a jamais ete renseigne.
        # Pose desormais QUE si `metadata_ok` est faux ci-dessus. Defaut
        # "y" (et fichier pre-selectionne) si la table d'apercu ci-dessus
        # vient d'etre exportee - cas d'usage attendu : l'utilisateur
        # vient de la completer a la main.
        site_metadata = None
        add_meta = "n"
        if not metadata_ok:
            add_meta = self._console_input(
                "Fill in missing site formation/lithology/age\n"
                "from an updated complement table (e.g. the table "
                "just exported above)? y/N: ", "y" if preview_path else "n")
            if add_meta is None:
                return
        if add_meta.strip().lower() == "y":
            meta_path = filedialog.askopenfilename(
                title="Select the site metadata table",
                initialdir=os.path.dirname(preview_path) if preview_path else None,
                initialfile=os.path.basename(preview_path) if preview_path else "",
                filetypes=[("Text", "*.txt"), ("All files", "*.*")],
            )
            if meta_path:
                try:
                    site_metadata = load_site_metadata_table(meta_path)
                except Exception as e:
                    self._showerror("Error", f"Could not load site metadata table:\n{e}")
                    return

        try:
            result = export_to_magic(
                export_samples, self.results if want_results else [], out_dir,
                lab_analysts=lab_analysts,
                continent_ocean=continent, country=country, region=region,
                anisotropy_kind_by_specimen=anisotropy_kind_by_specimen,
                anisotropy_skip=anisotropy_skip,
                pmagint_rows=pmagint_rows,
                aniso_tensors=aniso_tensors,
                aniso_mean_tensors=aniso_mean_tensors,
                site_metadata=site_metadata,
            )
        except OSError as e:
            self._showerror("Error", f"MagIC export failed:\n{e}")
            return

        summary = "\n".join(
            f"  {name}.txt: {result.counts[name]} line(s)" for name in result.paths)
        self._afficher(f"MagIC export finished in {out_dir}:\n{summary}\n")

    def _prompt_magstrat_heights(self):
        """Equivalent du prompt "with magnetostratigraphic data (y/N)"
        commun a exportpmagren et exporttolatex. Utilise en priorite le
        champ `stratigraphic_height` du .prmag (voir testlect.Pmag.
        stratigraphic_height) quand au moins un specimen de la selection
        l'a renseigne - meme logique que ouvrir_lismesdepth_dialog,
        demande explicite utilisateur ("check for all the functions where
        depth is used and invite the users to fill this parameter using
        complement") : plus besoin de (re)fournir un fichier "id depth"
        externe a chaque export des que l'information vit deja dans le
        fichier de donnees. Ne pose la question "with magnetostratigraphic
        data (y/N)" (ancien comportement, fichier externe) que si AUCUN
        specimen de la selection n'a de position renseignee - et invite
        alors explicitement a la renseigner une fois pour toutes via
        `Complete sample information...` plutot que de refournir un
        fichier a chaque export.
        Retourne None si annule (Échap), sinon un dict {id: depth} (vide
        si aucune position n'est disponible et que l'utilisateur repond
        non)."""
        heights = {
            ech.id: ech.stratigraphic_height
            for ech in self.selection if ech.stratigraphic_height is not None
        }
        if heights:
            n_missing = len(self.selection) - len(heights)
            self._afficher(
                f"Using stratigraphic_height from the .prmag file "
                f"({len(heights)} specimen(s) with a known position"
                + (f", {n_missing} without)" if n_missing else ")")
                + ".\n"
            )
            if n_missing:
                self._afficher(
                    "Tip: fill in the missing position(s) with PmagFile -> "
                    "Complete sample information... (specimen mode, "
                    "'stratigraphic_height' column) so future exports include "
                    "them automatically.\n"
                )
            return heights

        cmagstrat = self._console_input("With magnetostratigraphic data (y/N): ", "N")
        if cmagstrat is None:
            return None
        if cmagstrat.strip().lower() != "y":
            return {}
        self._afficher(
            "Tip: instead of a one-off height file, you can fill in "
            "'stratigraphic_height' directly in the .prmag with PmagFile -> "
            "Complete sample information... - it will then be picked up "
            "automatically here (and in List and depth...) without asking "
            "for a file again.\n"
        )
        path = filedialog.askopenfilename(
            title="Sample/height file (positive upward)",
            filetypes=[("Text", "*.txt *.dat"), ("All files", "*.*")],
        )
        if not path:
            return {}
        heights = {}
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            heights[parts[0].strip()] = float(parts[1])
                        except ValueError:
                            continue
        except OSError as e:
            self._showwarning("Error", f"Could not read {path}:\n{e}")
        return heights

    def ouvrir_export_detailed_dialog(self):
        """Equivalent GUI de `exportpmagren` ("export detailed Rennes",
        dataselect.f) : fichier texte unique avec, pour chaque échantillon
        sélectionné, un bloc de paramètres complet puis le tableau de
        mesures (Dsc/Isc, Dis/Iis, Dtc/Itc = dec/inc en repère échantillon/
        in-situ/après pendage).
        La déclinaison IGRF n'est pas calculée (IGRF non porté) - affichée
        "n.d" comme le fait déjà le Fortran pour ses propres cas de
        données manquantes.

        Genere aussi automatiquement une copie PDF du meme texte (police
        MONOSPACE, meme alignement garanti quel que soit l'ordinateur/
        editeur de texte qui l'ouvre - voir export_detailed_txt) - demande
        explicite utilisateur ("a PDF copy of the text export: do it
        automatically with the export text"), qui remplace l'export Latex
        desormais retire (necessitait une chaine LaTeX installee juste
        pour obtenir un PDF, alors que STARpaleomag_Py lui-meme reste le
        bon outil pour consulter les donnees - demande explicite
        utilisateur : "a paleomagnetist who want to see the data can do
        it through STARpaleomag_Py")."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        self.text_area.insert(tk.END, "\n--- export detailed Rennes (Escape to cancel) ---\n", "prompt")
        location = self._console_input("Main study location (country, region): ", "")
        if location is None:
            return
        heights = self._prompt_magstrat_heights()
        if heights is None:
            return
        # Nom/dossier par defaut derives du fichier de donnees charge
        # (.prmag/.ren -> "..._listing.txt", et donc "..._listing.pdf"
        # pour la copie PDF automatique - meme base, voir
        # export_detailed_txt) plutot qu'un dialogue vide a chaque fois -
        # demande explicite utilisateur ("can you save the text file by
        # replacing .prmag by _listing.txt. and _listing.pdf"). Reste
        # modifiable dans le dialogue, ce n'est qu'une valeur initiale.
        initial_file, initial_dir = None, None
        if self.data_file_path:
            initial_file = os.path.splitext(os.path.basename(self.data_file_path))[0] + "_listing.txt"
            initial_dir = os.path.dirname(self.data_file_path)
        out_path = filedialog.asksaveasfilename(
            title="Detailed .txt file", defaultextension=".txt",
            initialfile=initial_file, initialdir=initial_dir,
            filetypes=[("Text", "*.txt"), ("All files", "*.*")],
        )
        if not out_path:
            return
        source_file = os.path.basename(self.data_file_path) if self.data_file_path else None
        try:
            export_detailed_txt(
                self.selection, location, out_path, heights=heights or None, source_file=source_file)
        except OSError as e:
            self._showerror("Error", f"Export failed:\n{e}")
            return
        pdf_path = os.path.splitext(out_path)[0] + ".pdf"
        self._afficher(f"Detailed export written: {out_path}\n  (and PDF copy: {pdf_path})\n")

    def ouvrir_export_stereo_dialog(self):
        """Equivalent GUI de `exportres` (fichiers.f:1203), adapte pour
        etre DIRECTEMENT compatible avec StereoUtils_Py plutot que de
        reprendre tel quel le format texte du Fortran d'origine - demande
        explicite utilisateur ("is it possible to add the export results
        to Stereo_Py as in the original Fortran with in addition a
        specific file for the poles"), affinee au fil de plusieurs
        demandes explicites ulterieures jusqu'a l'actuel fichier UNIQUE a
        3 blocs (voir export_stereo.export_stereo_project pour le detail
        complet du format) : "individual directions" (lignes/plans de
        specimen), "mean directions" (moyennes de site), "VGP" (poles,
        quand la moyenne en porte un) - chacun avec son propre en-tete de
        colonnes, une colonne `tilt_correction` (0/100, meme convention
        MagIC que partout ailleurs dans ce projet) et une colonne `info`
        VOLONTAIREMENT vide, reservee a l'utilisateur ("can you manage a
        stereo project with three blocks... a column for the tilt
        correction 0 to 100 and an information column to be filled
        manually by the user")."""
        if not self.results:
            self._showwarning("No results", "Run one or more fits first.")
            return
        orient_s = self._console_input(
            "Export directions in-situ (i) or after tilt correction (t) : ", "i")
        if orient_s is None:
            return
        orientation = 3 if orient_s.strip().lower().startswith("t") else 2

        out_path = filedialog.asksaveasfilename(
            title="Export results to Stereo_Py - Project file",
            defaultextension=".txt",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")],
        )
        if not out_path:
            return

        counts = export_stereo_project(self.results, out_path, orientation=orientation)
        if not any(counts.values()):
            self._showinfo("Nothing to export", "No exportable line/plane/mean result found.")
            return
        self._afficher(
            f"Stereo project written to {out_path} "
            f"({counts['directions']} direction(s), {counts['means']} mean(s), "
            f"{counts['vgp']} VGP(s)) - StereoUtils_Py: Project > Load Project...\n"
        )

    # ------------------------------------------------------------------
    # Graphiques : Zijderveld (equivalent de `plotzijder`/`zijderplot`)
    # ------------------------------------------------------------------

    def load_demo_data(self):
        self.text_area.insert(tk.END, "Ready to import data...\n")

    def _clear_figure(self):
        """Efface juste le CONTENU de la Figure (axes), SANS appeler
        _redraw_canvas - car celle-ci ajuste la geometrie de la fenetre/du
        panneau en fonction de `self.fig.get_size_inches()`, qui a ce stade
        vaut encore la taille de l'ANCIEN graphique affiche, pas celle du
        nouveau qui va etre construit juste apres (`set_size_inches` suit
        toujours cet appel). Appeler `_redraw_canvas()` ici ferait un
        premier ajustement avec la MAUVAISE taille, en course avec le bon
        ajustement suivant - c'etait la cause du glitch d'affichage (taille
        de Figure corrompue, ex. [6.2, 7.7] au lieu de [6.0, 8.5]) qui
        persistait meme apres avoir corrige `_redraw_canvas` elle-meme."""
        self.fig.clear()

    def clear_screen(self):
        self._clear_figure()
        self._redraw_canvas()
        self._current_graphic = None
        self._arai_state = None
        self._paleoint_review_state = None

    def exporter_svg(self):
        """Equivalent fonctionnel de laser1/laser2/... (menu Graphics_SVG) :
        contrairement au Fortran, pas besoin d'un code d'export separe -
        matplotlib ecrit un SVG directement depuis la Figure affichee."""
        if not self.fig.axes:
            self._showwarning("No graphic", "Display a graphic first (e.g. Zijderveld).")
            return

        sample_id = self._current_graphic[1] if self._current_graphic and self._current_graphic[1] else None
        default_name = f"{sample_id}.svg" if sample_id else "graphic.svg"
        path = filedialog.asksaveasfilename(
            title="Export as SVG",
            defaultextension=".svg",
            initialfile=default_name,
            filetypes=[("SVG", "*.svg"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            self.fig.savefig(path, format="svg")
        except Exception as e:
            self._showerror("Error", f"SVG export failed:\n{e}")
            return
        self._showinfo("Export successful", f"Graphic exported:\n{path}")

    def ouvrir_user_guide(self):
        """Ouvre le guide utilisateur (STARpaleomag_Py Guide) dans le
        navigateur systeme - demande explicite utilisateur ("quel serait
        la possibilite de developper une aide en ligne par Claude pour
        des futurs utilisateurs ?"). Option retenue apres discussion des
        alternatives (une cle API partagee dans le code source serait
        rapidement retrouvee et abusee ; une cle par utilisateur ajoute
        de la friction pour un public qui n'est pas familier des API) :
        un fichier HTML STATIQUE livre EN LOCAL avec l'appli (voir
        _resource_path, help/STARpaleomag_Py_Guide.html - copie du guide
        deja redige comme artifact Claude, sans ses scripts de plateforme)
        plutot qu'une aide en direct necessitant reseau/cle/cout recurrent.
        A resynchroniser manuellement avec l'artifact si celui-ci est mis
        a jour par la suite (pas de lien automatique entre les deux)."""
        guide_path = _resource_path("help", "STARpaleomag_Py_Guide.html")
        if not os.path.exists(guide_path):
            self._showerror("Error", f"User guide not found:\n{guide_path}")
            return
        webbrowser.open(f"file://{guide_path}")

    def afficher_zijderveld(self):
        """Equivalent GUI de la boucle `do i=1,nbech ... call zijder2(...)`
        (linesplans.f) : quand plusieurs echantillons sont selectionnes, on
        les parcourt SEQUENTIELLEMENT, un `read(*,*)` (ici `_console_input`)
        servant de PAUSE entre chaque - pas pour lire une valeur utile, juste
        pour laisser le temps de regarder le diagramme avant de passer au
        suivant (Echap pour arreter la sequence en cours de route)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        samples = [s for s in self.selection if s.mesures]
        if not samples:
            self._showwarning("No measurements", "No selected sample has any measurement.")
            return

        for i, ech in enumerate(samples):
            self._current_graphic = ("zijderveld", ech.id)
            try:
                self._refresh_current_graphic()
            except Exception as e:
                # Un echantillon ne devrait normalement plus faire planter
                # draw_zijderveld (voir le garde-fou "no NRM data" ajoute
                # dans zijderveld.py), mais si un cas imprevu se presente
                # quand meme, on ne veut pas interrompre la sequence pour
                # les autres echantillons selectionnes - demande explicite
                # utilisateur : "it seems that the zijderveld cannot be
                # plotted and the app jump out... can you just skip the
                # sample". On log et on passe au suivant plutot que de
                # laisser l'exception remonter jusqu'au callback Tkinter.
                self._afficher(f"[Zijderveld] {ech.id}: skipped ({e})\n")
                self._current_graphic = None
                if i == len(samples) - 1:
                    break
                continue
            if i == len(samples) - 1:
                break
            reste = len(samples) - i - 1
            answer = self._console_input(
                f"[{ech.id}] Enter for next sample ({reste} remaining, Escape to stop): "
            )
            if answer is None:
                break

    def afficher_stereo(self):
        """Equivalent de `stermes` (stereoplot(0))."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        self._current_graphic = ("stereo", None)
        self._refresh_current_graphic()

    def afficher_stereo_results(self):
        """Equivalent de `sterres` (stereoplot(1)) : directions moyennes de
        Fisher avec cone de confiance alpha95 (cat1='F'), et grands cercles/
        points pour les ajustements de droite/plan (cat1 'L'/'f'/'P') de
        self.results (a charger au prealable via Ajustement... ou
        Select results...)."""
        if not self.results:
            self._showwarning(
                "No results",
                "No results in memory - run a fit or "
                "load some via Pmag data > Select results...")
            return
        self._current_graphic = ("stereo_results", None)
        self._refresh_current_graphic()

    def afficher_xygraph(self):
        """Equivalent de `xygraph` (courbe de désaimantation)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        self._current_graphic = ("xygraph", None)
        self._refresh_current_graphic()

    def afficher_susceptibilite(self):
        """Equivalent de `suscep`."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        self._current_graphic = ("susceptibility", None)
        self._refresh_current_graphic()

    def afficher_irm(self):
        """Acquisition d'aimantation remanente isotherme (IRM, mesures
        cod1='I') - pas dans le Fortran, port des scripts GMT de
        l'utilisateur (Scripts_IRM_GMT), demande explicite."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        if not has_irm_data(self.selection):
            self._showwarning(
                "No IRM data",
                "No selected sample has an IRM measurement (step code 'I').",
            )
            return
        self._current_graphic = ("irm", None)
        self._refresh_current_graphic()

    def afficher_arai(self):
        """Equivalent GUI de `paleoin` (diagramme d'Arai/Thellier) - reprend
        la sequence interactive du Fortran (plotpaleoint2.f) presque telle
        quelle : champ com:, tableau des pas NRM/TRM, verifications pTRM,
        courbure de Paterson (Taubin+LMA, adjustcircle.f95, deux appels
        comme l'original), directions ancree/libre + DANG, tableau CRM
        (%rcrm), statistiques de Coe (b/sigma/ccr/f/g/q/h), correction
        d'anisotropie (tenseur .ANI 'A0') et gamma. Verifie octet-pres
        contre un transcript Fortran reel (echantillon 06A,
        SanJuan_Pmag.ren, ni,nj=4,10) - y compris la coquille reelle du
        source (constante '3.14152927' au lieu de pi dans `compute_crm`),
        preservee pour fidelite.

        Volontairement HORS PERIMETRE (non portes) : ARN a 2/3 composantes
        (icomposante 2/3 - traites comme 1), rf1/rf2 et le bloc Prevot Ql
        (code mort cote Fortran - la lecture de susceptibilite y est
        commentee, `sus`/`ql` valent toujours 0)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        self.text_area.insert(tk.END, "\n--- Arai diagram / paleointensity (Escape to cancel) ---\n", "prompt")

        ani_path = self._find_ani_path()
        pmagint_path = pmagint_path_for(self.results_path) if self.results_path else None

        traite = 0
        for ech in self.selection:
            if len(ech.mesures) < 2:
                continue
            # orientation=1 (echantillon) FORCE, pas self.orientation.get()
            # - Fortran (plotpaleoint2.f/visi_Paleoint.f) : "iorient=1"
            # sans condition en tete de routine, ignore le reglage
            # d'orientation de l'appli - demande explicite utilisateur
            # ("the combined plots should be in sample coordinates").
            points, checks, arno = compute_arai(ech, orientation=1)
            if len(points) < 2:
                continue  # pas de codes N/R/V/P exploitables pour cet echantillon
            traite += 1

            # champ/methode auto-detectes depuis les mesures elles-memes
            # (treat_dc_field pour .prmag, cod1 S/R/V) plutot que du seul
            # commentaire - demande explicite utilisateur ("extract the
            # field value from treat_dc_field... the method is IZZI if
            # there are S codes and Thellier if there are R and V
            # codes"). Repli sur l'ancien parsing du commentaire si aucune
            # mesure ne porte de treat_dc_field (.ren).
            method_auto, ichamp_auto = detect_method_and_hlab(ech.mesures, ech.com)
            confirm = self._console_input(
                f" method={method_auto}, lab field detected = {ichamp_auto:g} uT "
                f"(comment=[{ech.com}]) - is this OK Y/n : ", "Y")
            if confirm is None:
                return
            ncom = ichamp_auto != 0 and confirm.strip().lower() != "n"
            ichamp = ichamp_auto if ncom else 0

            icomp_s = self._console_input(
                " ARN one (1) two (2) or several components (3) : ", "1")
            if icomp_s is None:
                return
            try:
                icomposante = int(icomp_s.strip() or "1")
            except ValueError:
                icomposante = 1
            if icomposante not in (2, 3):
                icomposante = 1
            if icomposante != 1:
                self._afficher(
                    "(ARN 2/3 components not yet ported - treated as 1 component.)\n")

            plot_inv = self._console_input(" Plot TRM in - and + : y/N", "N")
            if plot_inv is None:
                return

            if ichamp != 0:
                hlab = float(ichamp)
                parsed = parse_com_field(ech.com)  # imois/iannee : purement informatif
                self._afficher(
                    f" ech:{ech.id:<12}  field of {ichamp}t "
                    f"month:{parsed['imois']}year:{parsed['iannee']}\n"
                )
            else:
                hlab_s = self._console_input(
                    "  intensity of the laboratory field (in microteslas): ", "0")
                if hlab_s is None:
                    return
                try:
                    hlab = float(hlab_s or 0.0)
                except ValueError:
                    self._showerror("Error", "Hlab must be a number.")
                    continue

            volumasse = 1.0e-3 if ech.norme == "m" else 1.0

            out = io.StringIO()
            out.write(" paleointensity diagram NRM-TRM\n")
            out.write(f"  numero: {ech.id}\n\n")
            out.write("    n  temp.    arn     dec    inc      atr      mom.atr    dec    inc\n\n")
            for i, p in enumerate(points, start=1):
                mom = p.xp * arno * volumasse
                out.write(
                    f"  {i:3d}  {p.temp:4.0f}   {p.yp:9.3E}  {p.decl:7.2f} {p.aincl:7.2f}   "
                    f"{p.xp:9.3E}   {mom:9.3E}   {p.dec:7.1f} {p.winc:7.1f}\n"
                )
            if checks:
                pt_by_k = {p.k: p for p in points}
                out.write("\n   results of the ptrm checks\n")
                # 2 colonnes de temperature distinctes - demande explicite
                # utilisateur ("you also miss the step at which the ptrm
                # is done") : `done_at` = pas COURANT ou le controle a ete
                # physiquement effectue (c.k, voir compute_arai), `target`
                # = pas revisite par le controle (c.temp) - auparavant
                # affichees identiques par erreur (les deux calculees a
                # partir du meme "target").
                out.write("   done_at  target     mom       %atr   ecart/atrt  ecart/atrp\n")
                for c in checks:
                    done_at = pt_by_k[c.k].temp if c.k in pt_by_k else c.temp
                    mom = c.xt * arno
                    ecart = c.xt - c.xtptrm
                    ratio = ecart / c.xtptrm if c.xtptrm else 0.0
                    out.write(
                        f"   {done_at:5.0f}  {c.temp:5.0f}   {mom:9.3E}   {c.xt:6.2f}   "
                        f"{ecart:6.2f}      {ratio:6.2f}\n"
                    )
            self._afficher(out.getvalue())

            self._arai_state = (ech, points, checks, arno, None, hlab)
            self._current_graphic = ("arai", ech.id)
            self._refresh_current_graphic()

            while True:
                range_s = self._console_input(
                    f" calculation of the slope between steps ni,nj (1-{len(points)}, 0 = skip): ",
                    f"1 {len(points)}")
                if range_s is None:
                    return
                parts = range_s.split()
                if len(parts) == 1 and parts[0].strip() == "0":
                    break
                try:
                    if len(parts) >= 2:
                        n1, n2 = int(parts[0]), int(parts[1])
                    else:
                        n1, n2 = 1, len(points)
                except ValueError:
                    self._showerror("Error", "Enter two integers (ni nj).")
                    continue
                if not (1 <= n1 <= n2 <= len(points)) or n2 - n1 < 1:
                    self._showerror("Error", "Invalid point range (at least 2 points).")
                    continue

                fit = fit_arai_line(points, n1, n2, hlab=hlab)
                direction = fit_arai_direction(points, n1, n2, ech, orientation=1)  # SC forcee, voir compute_arai ci-dessus
                gamma = 90.0 - points[n2 - 1].winc
                rf1, rf2 = compute_rf1_rf2(points, n1, n2, nrm0_vector(ech, orientation=1))

                fcor, hcorani, aniso_note, aniso_tensor = self._apply_anisotropy_correction(
                    ech, ani_path, direction, fit.h)

                # courbure de Paterson (adjustcircle.f95, AraiCurvature) :
                # appel 1 sur (0,1)+points[1..n2] (utilise dans la ligne de
                # resultats finale), appel 2 sur points[n1..n2] uniquement
                # (affiche mais pas repris dans la ligne finale - transcrit
                # tel quel, meme choix que le Fortran).
                xpat1 = [0.0] + [p.xp for p in points[:n2]]
                ypat1 = [1.0] + [p.yp for p in points[:n2]]
                curv0 = arai_curvature(xpat1, ypat1)
                xpat2 = [p.xp for p in points[n1 - 1:n2]]
                ypat2 = [p.yp for p in points[n1 - 1:n2]]
                curv1 = arai_curvature(xpat2, ypat2)

                crm = compute_crm(points, n1, n2, arno) if n2 - n1 >= 1 else None

                self._arai_state = (ech, points, checks, arno, fit, hlab)
                self._refresh_current_graphic()

                res = io.StringIO()
                res.write(f"\n calculation of the slope between steps ni,nj :{n1} {n2}\n\n")
                res.write(" -----------------------------------\n")
                res.write(
                    f" circle Parameter for initialization after Taubin: {curv0.taubin_a:9.6f} "
                    f"{curv0.taubin_b:9.6f} {curv0.taubin_r:8.5f}\n"
                )
                res.write(
                    f" circle Parameter (a,b,r) after LMA: {curv0.lma_a:9.6f} "
                    f"{curv0.lma_b:9.6f} {curv0.lma_r:8.5f}\n"
                )
                res.write(" curvature calculated from point (0,1) to n2\n")
                res.write(f" parameter k (1/r) : {curv0.k:7.4f}  error SSE :{curv0.sse:7.5f}\n")
                res.write("\n -----------------------------------\n")
                res.write(
                    f" circle Parameter for initialization after Taubin: {curv1.taubin_a:9.6f} "
                    f"{curv1.taubin_b:9.6f} {curv1.taubin_r:8.5f}\n"
                )
                res.write(
                    f" circle Parameter (a,b,r) after LMA: {curv1.lma_a:9.6f} "
                    f"{curv1.lma_b:9.6f} {curv1.lma_r:8.5f}\n"
                )
                res.write(" curvature calculated from point n1 to n2\n")
                res.write(f" parameter k (1/r) : {curv1.k:7.4f}  error SSE :{curv1.sse:7.5f}\n")
                res.write("\n -----------------------------------\n\n")
                if direction.anchored_dec is not None:
                    res.write(
                        f" anchored direction: dec={direction.anchored_dec:6.1f}  "
                        f"inc={direction.anchored_inc:6.1f}  mad={direction.anchored_mad:5.1f}  "
                        f"nb points: {direction.nb}\n"
                    )
                if direction.free_dec is not None:
                    res.write(
                        f"free direction:  dec={direction.free_dec:6.1f}  "
                        f"inc={direction.free_inc:6.1f}  mad={direction.free_mad:5.1f}  "
                        f"nb points: {direction.nb}\n"
                    )
                if direction.dang is not None:
                    res.write(f"\n Lisa Tauxe DANG {direction.dang:6.1f}\n")

                # Directions AVANT/APRES correction d'anisotropie - demande
                # explicite utilisateur ("ce serait bien d'afficher les
                # directions (anchored not anchored) avant et apres
                # correction d'anisotropie ; ça donne une idée de la
                # déviation des directions du champ par l'anisotropie de
                # la roche"). Seulement si un tenseur A0 a reellement ete
                # utilise (aniso_tensor non None - donc PAS si le
                # specimen est deja flaganiso, ni si l'utilisateur a
                # decline "not satisfactory... anyway?").
                if aniso_tensor is not None:
                    corrected_direction = fit_arai_direction_corrected(
                        points, n1, n2, ech, aniso_tensor, orientation=1)
                    res.write("\n --- after anisotropy correction (tensor 'A0') ---\n")
                    if corrected_direction.anchored_dec is not None:
                        res.write(
                            f" anchored direction: dec={corrected_direction.anchored_dec:6.1f}  "
                            f"inc={corrected_direction.anchored_inc:6.1f}  "
                            f"mad={corrected_direction.anchored_mad:5.1f}\n"
                        )
                    if corrected_direction.free_dec is not None:
                        res.write(
                            f"free direction:  dec={corrected_direction.free_dec:6.1f}  "
                            f"inc={corrected_direction.free_inc:6.1f}  "
                            f"mad={corrected_direction.free_mad:5.1f}\n"
                        )
                    dev_anchored = angle_between_vectors(
                        direction.anchored_specimen_frame, corrected_direction.anchored_specimen_frame)
                    dev_free = angle_between_vectors(
                        direction.free_specimen_frame, corrected_direction.free_specimen_frame)
                    if dev_anchored is not None:
                        res.write(f" deviation from anisotropy (anchored): {dev_anchored:5.1f} deg\n")
                    if dev_free is not None:
                        res.write(f" deviation from anisotropy (free):     {dev_free:5.1f} deg\n")

                if rf1 is not None:
                    res.write(f" rf1 - rf2 : {rf1:5.2f}  {rf2:5.2f}\n")
                if crm is not None:
                    res.write("\n temp   mom.crm       temp   mom.2,6       temp   mom.crm\n")
                    ks = list(range(n1, n2 + 1))
                    k = n1
                    while k <= n2:
                        kkj = min(k + 2, n2)
                        cols = []
                        for j in range(k, kkj + 1):
                            cols.append(f"{points[j-1].temp:4.0f}   {crm.values[j]:10.4E}")
                        res.write("  " + "     ".join(cols) + "\n")
                        k += 2
                res.write(f"\n slope b= {fit.b:8.4f}\n")
                if fit.sigma:
                    res.write(f"\n sigma = {fit.sigma:7.4f}\n")
                    res.write(f"\n linear correlation coefficient: {fit.ccr:9.5f}\n")
                if crm is not None:
                    res.write(f"\n crm max ={crm.crmmax:8.4f}    deltatrm={crm.dtrm:8.4f}\n")
                res.write(
                    f"\n{ech.id:<12}   f={fit.f:7.4f}    g={fit.g:7.4f}    q={fit.qq:7.4f}"
                    f"    ccr={fit.ccr:8.5f}  h={fit.h:8.3f}"
                    + (f"     % rcrm={crm.rcrm:7.3f}\n" if crm is not None else "\n")
                )
                if aniso_note:
                    res.write(f"\n{aniso_note}")
                res.write(
                    "\nNum          t1   t2   N    f       g       q     mad   dang   Hlab"
                    "     b       sb     sb/b     ccr      H      fcor  Hcorani  gamma     k     k_sse\n"
                )
                res.write(
                    f"{ech.id:<12} {points[n1-1].temp:4.0f} {points[n2-1].temp:4.0f}  "
                    f"{n2 - n1 + 1:2d}  {fit.f:6.3f}  {fit.g:6.3f}  {fit.qq:6.3f}  "
                    f"{(direction.free_mad or 0.0):5.1f}  {(direction.dang or 0.0):5.1f}  "
                    f"{hlab:5.1f}  {fit.b:7.4f}  {fit.sigma:6.4f}  "
                    f"{(fit.sigma / fit.b if fit.b else 0.0):7.4f}  {fit.ccr:8.5f}  "
                    f"{fit.h:6.2f}  "
                    + (f"{fcor:6.3f}  {hcorani:6.2f}" if fcor is not None else "   -       -  ")
                    + f"  {gamma:5.1f}  {curv0.k:7.4f}  {curv0.sse:7.5f}\n"
                )
                self._afficher(res.getvalue())

                # Meme traitement PmagPy/MagIC parallele que le batch
                # review (compute_magic_paleointensity) - AJOUTE ici (pas
                # calcule auparavant dans ce dialogue interactif) pour que
                # mad/dang/fvds/frac/gap_max/n_ptrm du .pmagint viennent
                # de la MEME source PmagPy dans les deux dialogues -
                # demande explicite utilisateur (voir write_pmagint_line).
                magic_result = None
                try:
                    magic_result = compute_magic_paleointensity(
                        ech, step_first=points[n1 - 1].temp, step_last=points[n2 - 1].temp)
                    self._afficher(format_magic_paleointensity(magic_result))
                except Exception as e:
                    self._afficher(f"MagIC/PmagPy paleointensity: not computed for {ech.id} ({e})\n")

                if pmagint_path:
                    write_pmagint_line(
                        pmagint_path, ech.id, points[n1 - 1].temp, points[n2 - 1].temp,
                        n2 - n1 + 1, fit.f, fit.g, fit.qq, hlab, fit.b, fit.sigma, fit.ccr, fit.h,
                        mad=(magic_result.mad if magic_result else None),
                        dang=(magic_result.dang if magic_result else None),
                        pct_crm=(crm.rcrm if crm is not None else None),
                        fcor=fcor, hcorani=hcorani,
                        gamma=gamma, k=curv0.k, k_sse=curv0.sse,
                        fvds=(magic_result.fvds if magic_result else None),
                        frac=(magic_result.frac if magic_result else None),
                        gap_max=(magic_result.gap_max if magic_result else None),
                        n_ptrm=(magic_result.n_ptrm if magic_result else None),
                        tensor=aniso_tensor, f1=rf1, f2=rf2,
                    )

                decision = self._console_input(
                    "Redo (r) / next sample (Enter): ", "")
                if decision is None:
                    return
                if decision.strip().lower() == "r":
                    continue
                break

        if traite == 0:
            self._showinfo(
                "No diagram",
                "No selected sample has usable paleointensity measurements (N/R/V/P codes).",
            )

    def ouvrir_view_pmagint_dialog(self):
        """Vue en lot READ-ONLY des interpretations DEJA archivees dans
        .pmagint - demande explicite utilisateur ("differencier la
        visualisation des donnees deja dans pmagint, d'un redo a partir
        du fichier simplifie... dans la premiere option, l'utilisateur
        n'est pas sense modifier l'interpretation deja dans pmagint").
        Sibling de `ouvrir_openfilepint_dialog` ("Rapid view/recompute
        from redo file...", ex-"View batch of Paleoint Results...", nom
        trompeur - CETTE routine recalcule un traitement complet depuis
        un fichier redo EXTERNE et REECRIT .pmagint a chaque revue,
        contrairement a ce que son ancien libelle laissait penser) :
        meme UI (menu numerote, Entree = suivant, q/Echap = quitter),
        mais deux differences essentielles :

        1) La SOURCE de la liste est `.pmagint` lui-meme (`read_pmagint`),
           pas un fichier redo separe a fournir par l'utilisateur - Tmin/
           Tmax de CHAQUE ligne viennent directement du fichier deja
           archive.
        2) AUCUNE ecriture : `write_pmagint_line` n'est jamais appele ici.
           Le texte affiche (f/g/q/ccr/H/mad/dang/fcor/Hcorani/gamma/k/
           fvds/frac/gap_max/n_ptrm/tensor...) est la valeur STOCKEE,
           relue telle quelle - PAS recalculee (evite tout risque de
           deviation silencieuse si les donnees brutes ou un tenseur
           .pmagani ont change depuis l'archivage). Ni prompt Hlab, ni
           second calcul PmagPy/MagIC, ni correction d'anisotropie
           recalculee ici - tout ca appartient au traitement (Arai
           interactif ou "Rapid view/recompute from redo file..."), pas
           a une simple relecture.

        Le diagramme (Arai/Zijderveld/stereo, `build_paleoint_review_
        figure`) DOIT neanmoins etre redessine depuis les mesures brutes
        (.pmagint ne stocke pas la serie de points, seulement les
        statistiques resumees) : `fit_arai_line(points, n1, n2, hlab=
        Hlab_stocke)` est appele UNIQUEMENT pour reconstruire la
        geometrie du trace, fonction PURE et DETERMINISTE de (points,
        n1, n2, Hlab) - donc reproduit exactement le meme resultat que
        l'archivage d'origine tant que les mesures/le fichier .pmagani
        n'ont pas change, sans compter comme une "modification" de
        l'interpretation au sens ou l'entend l'utilisateur (aucune
        ecriture, aucun nouveau choix Tmin/Tmax/Hlab)."""
        if not self.donnees:
            self._showwarning("No data", "Load a .ren file first.")
            return
        pmagint_path = pmagint_path_for(self.results_path) if self.results_path else None
        if not pmagint_path or not os.path.exists(pmagint_path):
            self._showwarning(
                "No .pmagint file", "No paleointensity results file (.pmagint) found for this data set.")
            return
        stored = read_pmagint(pmagint_path)
        if not stored:
            self._showwarning("Empty", f"No interpretation found in {pmagint_path}.")
            return
        entries = list(stored.items())  # (specimen_id, row dict)

        self.text_area.insert(
            tk.END,
            "\n--- View batch of Paleoint Results (from .pmagint, read-only - Escape to quit) ---\n",
            "prompt",
        )

        last_iligne = None
        while True:
            menu = io.StringIO()
            for i in range(0, len(entries), 3):
                row = entries[i:i + 3]
                menu.write("   ".join(
                    f"{i + j + 1:3d}  {e[0]:<12}" for j, e in enumerate(row)) + "\n")
            self._afficher(menu.getvalue())

            choice = self._console_input(
                "Type the line number to select the sample (Enter = next, q to quit) : ", "")
            if choice is None:
                return
            choice = choice.strip()
            if choice.lower().startswith("q"):
                break
            if not choice:
                if last_iligne is None:
                    break
                iligne = last_iligne + 1
                if iligne > len(entries):
                    self._afficher("(no more entries)\n")
                    break
            else:
                try:
                    iligne = int(choice)
                except ValueError:
                    continue
                if not (1 <= iligne <= len(entries)):
                    continue
            last_iligne = iligne

            sample_id, row = entries[iligne - 1]
            matches = select_samples(
                self.donnees, sample_id, step_min=0, step_max=2000,
                demag1="*", demag2="*", verbose=False)
            ech = matches[0] if matches else None
            if ech is None:
                self._afficher(f"{sample_id}: sample not found in the loaded data.\n")
                continue
            self.selection = [ech]

            try:
                tmin, tmax = float(row.get("t1", "n.d")), float(row.get("t2", "n.d"))
            except ValueError:
                self._afficher(f"{sample_id}: invalid Tmin/Tmax stored in .pmagint\n")
                continue
            try:
                hlab_stored = float(row.get("Hlab", "n.d"))
            except ValueError:
                hlab_stored = 0.0

            points, checks, arno = compute_arai(ech, orientation=1)
            if len(points) < 2:
                self._afficher(f"{sample_id}: no usable paleointensity measurements (N/R/V/P).\n")
                continue
            n1 = next((i + 1 for i, p in enumerate(points) if p.temp == tmin), None)
            if n1 is None and points and tmin <= points[0].temp:
                n1 = 1
            n2 = next((i + 1 for i, p in enumerate(points) if p.temp == tmax), None)
            if n1 is None or n2 is None or n2 - n1 < 1:
                self._afficher(
                    f"{tmin:g} {tmax:g} please check the temperature interval for: {sample_id}\n"
                    "(raw data may have changed since this specimen was archived)\n")
                continue

            # fit UNIQUEMENT pour redessiner la droite du diagramme - voir
            # docstring ci-dessus (fonction pure, ne "modifie" rien).
            fit = fit_arai_line(points, n1, n2, hlab=hlab_stored)

            self._paleoint_review_state = (ech, points, checks, arno, fit)
            self._current_graphic = ("paleoint_review", sample_id)
            self._refresh_current_graphic()

            self._afficher(_format_pmagint_row(sample_id, row))

            pause = self._console_input(
                "Press Return to continue (Escape to quit): ", "")
            if pause is None:
                return

    def ouvrir_openfilepint_dialog(self):
        """Equivalent GUI de `openfilepint`/`visi_paleoin` ("View Paleoint
        Results", visi_Paleoint.f) : PAS un traitement en lot - l'objectif
        (comme le Fortran) est de revisiter RAPIDEMENT un traitement deja
        effectue, en laissant l'utilisateur choisir l'echantillon dans une
        liste numerotee (fichier "echantillon Tmin Tmax [taux refroid.]",
        une ligne par determination deja faite), rejouee en boucle jusqu'a
        'q' - exactement le motif `openfilepint` (menu numerote, prompt
        "Type the line number...(q to quit)") plutot qu'un traitement de
        toutes les lignes d'un coup.

        Le graphique associe est la mise en page combinee de
        `visi_paleoin` (`boite(1)/(2)/(3)`) : diagramme d'Arai en haut,
        Zijderveld en bas a gauche, stereo NRM/TRM en bas a droite
        (`build_paleoint_review_figure`) - PAS le panneau Arai seul.

        Comme le Fortran (qui cherche l'echantillon dans `pmag(:)` - la
        TOTALITE des donnees chargees, pas une selection prealable),
        chaque echantillon est retrouve dans self.donnees (meme raison que
        afficher_visres, cf. son docstring). Le champ `com:` fournit Hlab
        silencieusement (pas de prompt de confirmation, ni de choix ARN
        composante/signe TRM - tous fixes cote Fortran dans cette
        routine, contrairement a `afficher_arai`/`paleoin`) ; le taux de
        refroidissement optionnel (4e colonne) n'est qu'un multiplicateur
        final sur H (`rHfinal=H*corcool`), PAS la correction complete de
        `vitref`/Cooling rate.

        5e colonne optionnelle (aniso_corr, `int_corr_aniso` MagIC - voir
        extract_magic.magic_pint_results_to_redo_lines) : facteur
        d'anisotropie DEJA calcule par le logiciel source de la
        contribution MagIC importee - demande explicite utilisateur
        ("while writing the redo file for paleointensity, add the
        int_corr_cooling_rate, and int_corr_aniso in the redo file").
        Applique UNIQUEMENT EN REPLI, si STARpaleomag_Py n'a PAS de tenseur
        .pmagani natif pour ce specimen (voir _apply_anisotropy_correction
        - `aniso_tensor is None` ET l'echantillon n'a pas deja ete
        corrige vectoriellement) - demande explicite utilisateur, pour ne
        jamais corriger deux fois le meme effet physique (tenseur natif
        recalcule ET facteur MagIC deja applique par le logiciel
        source)."""
        if not self.donnees:
            self._showwarning("No data", "Load a .ren file first.")
            return
        list_path = filedialog.askopenfilename(
            title="List file (sample Tmin Tmax [cooling rate] [aniso factor])",
            filetypes=[("Text", "*.txt *.lst *.dat"), ("All files", "*.*")],
        )
        if not list_path:
            return
        try:
            with open(list_path, "r", encoding="iso-8859-1", errors="replace") as f:
                lines = [l for l in f.read().splitlines() if l.strip() and not l.strip().startswith("!")]
        except OSError as e:
            self._showerror("Error", f"Could not read {list_path}:\n{e}")
            return

        entries = []  # (sample_id, tmin, tmax, cooling, aniso_corr)
        for line in lines:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                sample_id = parts[0]
                tmin, tmax = int(float(parts[1])), int(float(parts[2]))
                cooling = float(parts[3]) if len(parts) > 3 else 0.0
                aniso_corr = float(parts[4]) if len(parts) > 4 else 0.0
            except ValueError:
                continue
            entries.append((sample_id, tmin, tmax, cooling, aniso_corr))
        if not entries:
            self._showwarning("Empty list", f"No usable entry found in {list_path}.")
            return

        ani_path = self._find_ani_path()
        pmagint_path = pmagint_path_for(self.results_path) if self.results_path else None

        self.text_area.insert(
            tk.END, "\n--- Rapid view of previous Paleointensity determinations (Escape to quit) ---\n", "prompt")

        last_iligne = None  # dernier numero de ligne affiche - permet d'avancer par simple Retour
        while True:
            menu = io.StringIO()
            for i in range(0, len(entries), 3):
                row = entries[i:i + 3]
                menu.write("   ".join(
                    f"{i + j + 1:3d}  {e[0]:<12}" for j, e in enumerate(row)) + "\n")
            self._afficher(menu.getvalue())

            choice = self._console_input(
                "Type the line number to select the sample (Enter = next, q to quit) : ", "")
            if choice is None:
                return
            choice = choice.strip()
            if choice.lower().startswith("q"):
                break
            if not choice:
                # Retour simple : avance sur la ligne suivante (demande
                # explicite utilisateur - "a simple return just
                # incremented the selection... You could also start to 11
                # and the next one was 12") - seulement s'il y a deja eu
                # un choix explicite ; sans historique, un Retour a vide
                # quitte comme avant (rien a partir de quoi avancer).
                if last_iligne is None:
                    break
                iligne = last_iligne + 1
                if iligne > len(entries):
                    self._afficher("(no more entries)\n")
                    break
            else:
                try:
                    iligne = int(choice)
                except ValueError:
                    continue
                if not (1 <= iligne <= len(entries)):
                    continue
            last_iligne = iligne

            sample_id, tmin, tmax, cooling, aniso_corr = entries[iligne - 1]
            matches = select_samples(
                self.donnees, sample_id, step_min=0, step_max=2000,
                demag1="*", demag2="*", verbose=False)
            ech = matches[0] if matches else None
            if ech is None:
                self._afficher(f"{sample_id} sample not found\n")
                continue
            # rend l'echantillon revisite ACTIF (self.selection), comme
            # toutes les autres routines qui affichent un echantillon
            # (afficher_zijderveld/afficher_arai...) - demande explicite
            # utilisateur ("does not select the sample") : cette routine
            # ne le faisait pas, seul self._current_graphic/
            # _paleoint_review_state (pilotant uniquement le graphique)
            # etaient mis a jour, donc les autres menus (Zijderveld,
            # Delete some data, export SVG...) ne voyaient pas
            # l'echantillon comme selectionne apres l'avoir revisite ici.
            self.selection = [ech]

            # orientation=1 (echantillon) FORCE, pas self.orientation.get()
            # - Fortran (plotpaleoint2.f/visi_Paleoint.f) : "iorient=1"
            # sans condition en tete de routine, ignore le reglage
            # d'orientation de l'appli - demande explicite utilisateur
            # ("the combined plots should be in sample coordinates").
            points, checks, arno = compute_arai(ech, orientation=1)
            if len(points) < 2:
                self._afficher(f"{sample_id}: no usable paleointensity measurements (N/R/V/P).\n")
                continue

            n1 = next((i + 1 for i, p in enumerate(points) if p.temp == tmin), None)
            # tmin==0 (ou tout tmin <= au premier point disponible) designe
            # generalement le pas NRM (etape 0) comme borne de depart - ce
            # pas ne figure jamais dans `points` (compute_arai ne l'utilise
            # que pour `arno`, jamais comme point du diagramme d'Arai) et
            # ne matchera donc jamais par egalite exacte - bug reel signale
            # par l'utilisateur ("on ne peut pas lire le dernier échantillon
            # de la liste" - en realite le SEUL echantillon de la liste
            # avec tmin=0, confirme sur kr24_09b6/magic_contribution_19987 :
            # le premier point disponible demarre a 200degC, jamais 0).
            # Sans borne anterieure existante, demarrer simplement au
            # premier point disponible.
            if n1 is None and points and tmin <= points[0].temp:
                n1 = 1
            n2 = next((i + 1 for i, p in enumerate(points) if p.temp == tmax), None)
            if n1 is None or n2 is None or n2 - n1 < 1:
                self._afficher(
                    f"{tmin} {tmax} please check the temperature interval for: {ech.id}\n")
                continue

            _method_auto, ichamp_auto = detect_method_and_hlab(ech.mesures, ech.com)
            if ichamp_auto:
                hlab = float(ichamp_auto)
            else:
                hlab_s = self._console_input(
                    "  intensity of the laboratory field (in microteslas): ", "0")
                if hlab_s is None:
                    return
                try:
                    hlab = float(hlab_s or 0.0)
                except ValueError:
                    hlab = 0.0

            fit = fit_arai_line(points, n1, n2, hlab=hlab)
            direction = fit_arai_direction(points, n1, n2, ech, orientation=1)  # SC forcee, voir compute_arai ci-dessus
            gamma = 90.0 - points[n2 - 1].winc
            rf1, rf2 = compute_rf1_rf2(points, n1, n2, nrm0_vector(ech, orientation=1))

            fcor, hcorani, aniso_note, aniso_tensor = self._apply_anisotropy_correction(
                ech, ani_path, direction, fit.h)

            # Repli sur le facteur d'anisotropie DEJA calcule par la
            # contribution MagIC (int_corr_aniso, 5e colonne du fichier
            # redo - voir extract_magic.magic_pint_results_to_redo_lines)
            # UNIQUEMENT si STARpaleomag_Py n'a pas de tenseur .pmagani natif pour
            # ce specimen (aniso_tensor is None) ET que l'echantillon n'a
            # pas deja ete corrige vectoriellement (flaganiso) - demande
            # explicite utilisateur, pour ne jamais corriger deux fois le
            # meme effet physique (tenseur natif recalcule ET facteur
            # MagIC deja applique par le logiciel source).
            if aniso_tensor is None and aniso_corr and not getattr(ech, "flaganiso", False):
                fcor = aniso_corr
                hcorani = fit.h * aniso_corr
                aniso_note = (
                    f"Anisotropy: no native .pmagani tensor for {ech.id} - falling back to "
                    f"the MagIC contribution's own int_corr_aniso: fcor={fcor:.3f}  "
                    f"H (raw)={fit.h:.2f}µT -> Hcorani (corrected)={hcorani:.2f}µT\n"
                )

            h_final = hcorani if hcorani is not None else fit.h
            if cooling:
                h_final = h_final * cooling

            crm = compute_crm(points, n1, n2, arno)
            xpat1 = [0.0] + [p.xp for p in points[:n2]]
            ypat1 = [1.0] + [p.yp for p in points[:n2]]
            curv0 = arai_curvature(xpat1, ypat1)

            self._paleoint_review_state = (ech, points, checks, arno, fit)
            self._current_graphic = ("paleoint_review", ech.id)
            self._refresh_current_graphic()

            res = io.StringIO()
            res.write(f"\n{ech.id} : Tmin={tmin} Tmax={tmax}\n")
            if direction.anchored_dec is not None:
                res.write(
                    f" anchored direction: dec={direction.anchored_dec:6.1f}  "
                    f"inc={direction.anchored_inc:6.1f}  mad={direction.anchored_mad:5.1f}\n"
                )
            if direction.free_dec is not None:
                res.write(
                    f" free direction:     dec={direction.free_dec:6.1f}  "
                    f"inc={direction.free_inc:6.1f}  mad={direction.free_mad:5.1f}\n"
                )
            res.write(
                f" f={fit.f:.3f}  g={fit.g:.3f}  q={fit.qq:.3f}  ccr={fit.ccr:.3f}  "
                f"h={fit.h:.2f}  % rcrm={crm.rcrm:.1f}\n"
            )
            res.write(
                f" Hlab={hlab:.1f}µT"
                + (f"  fcor={fcor:.3f}  Hcorani={hcorani:.1f}µT" if fcor is not None else "")
                + (f"  fcorCool={cooling:.3f}  HcorCool={h_final:.1f}µT" if cooling else "")
                + f"  gamma={gamma:.1f}  k={curv0.k:.4f}\n"
            )
            if rf1 is not None:
                res.write(f" rf1 - rf2 : {rf1:5.2f}  {rf2:5.2f}\n")
            if aniso_note:
                res.write(aniso_note)

            # Directions AVANT/APRES correction d'anisotropie - meme
            # demande explicite utilisateur que dans le dialogue Arai
            # interactif ("afficher les directions... avant et apres
            # correction d'anisotropie").
            if aniso_tensor is not None:
                corrected_direction = fit_arai_direction_corrected(
                    points, n1, n2, ech, aniso_tensor, orientation=1)
                res.write(" --- after anisotropy correction (tensor 'A0') ---\n")
                if corrected_direction.anchored_dec is not None:
                    res.write(
                        f" anchored direction: dec={corrected_direction.anchored_dec:6.1f}  "
                        f"inc={corrected_direction.anchored_inc:6.1f}  "
                        f"mad={corrected_direction.anchored_mad:5.1f}\n"
                    )
                if corrected_direction.free_dec is not None:
                    res.write(
                        f" free direction:     dec={corrected_direction.free_dec:6.1f}  "
                        f"inc={corrected_direction.free_inc:6.1f}  "
                        f"mad={corrected_direction.free_mad:5.1f}\n"
                    )
                dev_anchored = angle_between_vectors(
                    direction.anchored_specimen_frame, corrected_direction.anchored_specimen_frame)
                dev_free = angle_between_vectors(
                    direction.free_specimen_frame, corrected_direction.free_specimen_frame)
                if dev_anchored is not None:
                    res.write(f" deviation from anisotropy (anchored): {dev_anchored:5.1f} deg\n")
                if dev_free is not None:
                    res.write(f" deviation from anisotropy (free):     {dev_free:5.1f} deg\n")

            self._afficher(res.getvalue())

            # Second traitement, PARALLELE et INDEPENDANT du natif STARpaleomag_Py
            # ci-dessus : appelle le code PmagPy/MagIC reel (pmag.PintPars)
            # sur le MEME intervalle Tmin/Tmax - demande explicite
            # utilisateur ("a second parallel processing of paleointensity.
            # the one from Magic"). Echec attendu et non bloquant pour une
            # partie reelle des specimens (protocoles hors IZZI/Thellier
            # standard, voir paleointensity_magic.py) - la boucle de revue
            # continue, avec juste une ligne d'explication. Calcule AVANT
            # write_pmagint_line : mad/dang/fvds/frac/gap_max/n_ptrm de ce
            # fichier viennent de CE resultat PmagPy, pas du calcul natif
            # (demande explicite utilisateur, voir write_pmagint_line).
            magic_result = None
            try:
                magic_result = compute_magic_paleointensity(ech, step_first=tmin, step_last=tmax)
                self._afficher(format_magic_paleointensity(magic_result))
            except Exception as e:
                self._afficher(f"MagIC/PmagPy paleointensity: not computed for {ech.id} ({e})\n")

            if pmagint_path:
                write_pmagint_line(
                    pmagint_path, ech.id, points[n1 - 1].temp, points[n2 - 1].temp,
                    n2 - n1 + 1, fit.f, fit.g, fit.qq, hlab, fit.b, fit.sigma, fit.ccr, fit.h,
                    mad=(magic_result.mad if magic_result else None),
                    dang=(magic_result.dang if magic_result else None),
                    pct_crm=crm.rcrm,
                    fcor=fcor, hcorani=hcorani,
                    fcor_cool=(cooling if cooling else None), hcor_cool=(h_final if cooling else None),
                    gamma=gamma, k=curv0.k, k_sse=curv0.sse,
                    fvds=(magic_result.fvds if magic_result else None),
                    frac=(magic_result.frac if magic_result else None),
                    gap_max=(magic_result.gap_max if magic_result else None),
                    n_ptrm=(magic_result.n_ptrm if magic_result else None),
                    tensor=aniso_tensor, f1=rf1, f2=rf2,
                )

            # Pause avant de reafficher la liste numerotee - demande
            # explicite utilisateur ("could you put a pause before to
            # print the list again in order to be able to read the
            # listing") : sans elle, la liste + le prompt suivant
            # s'affichaient IMMEDIATEMENT apres le detail MagIC/PmagPy
            # (`format_magic_paleointensity`), poussant ce detail hors de
            # la vue au prochain auto-scroll (`_afficher`/`see(tk.END)`)
            # avant meme d'avoir pu le lire. Meme convention Echap que le
            # reste de cette boucle (quitte entierement, voir `choice is
            # None` plus haut).
            pause = self._console_input(
                "Press Return to continue (Escape to quit): ", "")
            if pause is None:
                return

    def afficher_visres(self):
        """Equivalent GUI de `visres` ("data+interpretation",
        plotorthog.f:8-149) : boucle sur CHAQUE résultat de self.results
        INDIVIDUELLEMENT (pas regroupé par échantillon) - pour cat1 in
        ('L','P') affiche le Zijderveld de CET échantillon (les autres
        ajustements du même échantillon ne sont PAS superposés, comme le
        Fortran qui réduit `tr` à un seul élément avant d'appeler
        `zijder`), panneau UNIQUE avec son mini-stereo intégré (voir
        build_zijderveld_figure/draw_zijderveld) qui montre maintenant
        AUSSI le résultat par-dessus les données brutes : un point de
        direction pour une ligne, le grand cercle pour un plan - MÊME
        plot "classique" pour les deux désormais - demande explicite
        utilisateur ("put the great circle on the stereo within the same
        classic plot that has already zijderveld and stereo; when the
        result is a line, the result goes to the zijderveld, when it is
        a plane it is plotted on the stereo", puis "I do not understand
        why the stereo is still on the right. The routine with the
        zijder and stereo above the zijder already exist" - le panneau
        séparé "Stereo Results" à côté forçait une figure large à 2
        panneaux, donc une fenêtre plus large que partout ailleurs dans
        l'appli). Pour cat1=='f' (direction
        individuelle sans Zijderveld associé) affiche Stereo Results
        seul. Pour cat1=='F' (moyenne de site) - IGNORÉ par le
        Fortran d'origine, désormais traité - affiche un Stereo Results
        combinant la moyenne (cône de confiance) ET ses composants
        individuels (voir calcul.mean_components, cross-référencés via
        `r.liste`) - demande explicite utilisateur ("when the results is
        a mean, can you plot the stereo with the mean and individual
        results"). Même motif pause/décision (Entrée = suivant, Échap =
        arrêter) que les autres boucles (ajuslig, Arai...).

        Comme le Fortran (qui reconstruit l'échantillon DIRECTEMENT depuis
        `pmag(:)`, la totalité des données chargées, avec etapmin=0/
        etapmax=2000/demag1='*'/demag2='*' figés dans le code - PAS depuis
        une sélection préalablement filtrée), l'échantillon de chaque
        résultat est retrouvé dans self.donnees et non dans self.selection :
        sinon, un résultat chargé via « Select results... » (selres) sans
        avoir aussi fait « Select samples... » sur le même échantillon ne
        trouvait jamais ses données et data+interpretation semblait ne
        rien voir malgré des résultats bien sélectionnés.

        Affiche aussi l'évaluation de qualité de CE résultat (voir
        interpretation_quality.evaluate_result) - demande explicite
        utilisateur ("dans l'option data+interpretation, est-ce possible
        de lister l'évaluation de l'échantillon?")."""
        if not self.results:
            self._showwarning(
                "No results", "No results saved (see Calcul > Ajustement...).")
            return
        if not self.donnees:
            self._showwarning(
                "No data", "Load a .ren file first.")
            return

        self.text_area.insert(tk.END, "\n--- Data + interpretation (Escape to stop) ---\n", "prompt")
        traite = 0
        for r in self.results:
            if r.cat1 not in ("L", "P", "f", "F"):
                continue

            if r.cat1 == "F":
                # N'affiche QUE les moyennes archivees dans l'orientation
                # COURANTE - demande explicite utilisateur ("If you are
                # in in situ, you select only the mean that are in in
                # situ... only the means in IS are plotted if you are in
                # IS") : self.results peut contenir un melange de
                # moyennes IS/TC (accumulees via plusieurs "Select
                # results..." a des orientations differentes) - une
                # moyenne qui ne correspond pas est ignoree ici (pas
                # meme tentee, plutot que de tracer un stereo vide -
                # voir aussi le meme filtre dans draw_stereo_results).
                if r.par3_mean != float(self.orientation.get()):
                    continue
                # Moyenne de site : plus de specimen unique a chercher dans
                # self.donnees (r.id = "mean: <site>", pas un id specimen) -
                # trace la moyenne (cone de confiance) ET ses composants
                # individuels ensemble sur un seul stereo, plutot que de
                # l'ignorer comme le faisait le Fortran d'origine - demande
                # explicite utilisateur ("when the results is a mean, can
                # you plot the stereo with the mean and individual
                # results"). build_stereo_results_figure accepte deja une
                # liste melant moyenne(s) et resultats individuels (voir
                # "Stereo Results") ; mean_components() retrouve les
                # composants presents dans self.results via r.liste (les
                # `c` combines) - vide si self.results ne contient que des
                # moyennes (ex. charge avec carselect='m'), auquel cas la
                # moyenne seule est quand meme tracee.
                components = mean_components(r, self.results)
                traite += 1
                self._clear_figure()
                self.fig.set_size_inches(5.5, 5.5, forward=True)
                build_stereo_results_figure(
                    [r] + components, orientation=self.orientation.get(),
                    nbech=len(components) or None, fig=self.fig)
                self._fit_figure_to_data()
                self._redraw_canvas()

                dec, inc = _correct_dec_inc(r, self.orientation.get())
                n_txt = f", {len(components)} individual result(s)" if components else ""
                self._afficher(f"{r.id}: dec={dec:.1f}  inc={inc:.1f}{n_txt}\n")
                decision = self._console_input("Next (Enter) / stop (Escape): ", "")
                if decision is None:
                    return
                continue

            matches = select_samples(
                self.donnees, r.id, step_min=0, step_max=2000, demag1="*", demag2="*", verbose=False)
            ech = matches[0] if matches else None
            if ech is None:
                continue
            traite += 1
            self._clear_figure()

            if r.cat1 in ("L", "P"):
                # Zijderveld SEUL (panneau unique), avec son mini-stereo
                # INTEGRE au-dessus qui montre desormais aussi le resultat
                # (point pour une ligne, grand cercle pour un plan - voir
                # draw_zijderveld/show_stereo) - demande explicite
                # utilisateur ("I do not understand why the stereo is
                # still on the right. The routine with the zijder and
                # stereo above the zijder already exist") : le panneau
                # separe "Stereo Results" a cote (build_zijderveld_
                # stereo_results_figure) forcait une figure large a 2
                # panneaux et, avec elle, une fenetre plus large que les
                # autres vues - remplace par EXACTEMENT le meme layout a
                # panneau unique que la vue "Zijderveld" normale
                # (_refresh_current_graphic, kind=='zijderveld'), donc le
                # meme comportement de redimensionnement de fenetre
                # partout dans l'appli.
                self.fig.set_size_inches(5.5, 8.5, forward=True)
                build_zijderveld_figure(ech, orientation=self.orientation.get(), fits=[r], fig=self.fig)
                self._fit_figure_to_data()
                self._redraw_canvas()
            else:
                self.fig.set_size_inches(5.5, 5.5, forward=True)
                build_stereo_results_figure(
                    [r], orientation=self.orientation.get(), nbech=1, fig=self.fig)
                self._fit_figure_to_data()
                self._redraw_canvas()

            dec, inc = _correct_dec_inc(r, self.orientation.get())
            self._afficher(f"{r.id} ({r.cat1}): dec={dec:.1f}  inc={inc:.1f}\n")
            # Evaluation de l'interpretation affichee - demande explicite
            # utilisateur ("dans l'option data+interpretation, est-ce
            # possible de lister l'evaluation de l'echantillon?"). None
            # pour un resultat 'f' (moyenne Fisher, pas un fit ligne/plan -
            # voir interpretation_quality.evaluate_result), silencieusement
            # pas de ligne d'evaluation dans ce cas, comme pour "mean:".
            report = evaluate_result(r, self.donnees)
            if report is not None:
                self._afficher("\n" + format_quality_report([report]) + "\n")
            decision = self._console_input("Next (Enter) / stop (Escape): ", "")
            if decision is None:
                return

        if traite == 0:
            self._showinfo(
                "No displayable result",
                "No line/plane/mean type result found (or the referenced "
                "specimens are not loaded in self.donnees).",
            )

    def ouvrir_convertthelli_dialog(self):
        """Equivalent GUI de `convertthelli` ("Thellier >> NRM") : convertit
        la séquence Thellier de chaque échantillon sélectionné en séquence
        NRM/demag simple (cod1='D'), utilisable comme un Zijderveld classique."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        if not messagebox.askyesno(
            "Thellier >> NRM",
            f"Convert the Thellier sequence of {len(self.selection)} sample(s) "
            "into a simple NRM/demag sequence - irreversible for this session. Continue?",
        ):
            return
        for ech in self.selection:
            before = len(ech.mesures)
            convert_thellier_to_nrm(ech)
            self._afficher(f"{ech.id}: {before} -> {len(ech.mesures)} measurement(s).\n")

    def ouvrir_removestep_dialog(self):
        """Equivalent GUI de `removestep` ("Remove step") : supprime toutes
        les lignes d'un palier de démagnétisation pour un échantillon."""
        if len(self.selection) != 1:
            self._showwarning("Invalid selection", "Select a single sample.")
            return
        ech = self.selection[0]
        steps_text = "  ".join(
            f"{i + 1}:{_step_token(m)}" for i, m in enumerate(ech.mesures))
        self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")
        etape_s = self._console_input("Step (temperature/field) to remove: ", "")
        if etape_s is None:
            return
        try:
            etape = int(etape_s)
        except ValueError:
            self._showerror("Error", "Must be an integer.")
            return
        if not messagebox.askyesno(
            "Remove step",
            f"Remove all measurements of {ech.id} at step {etape} - irreversible "
            "for this session. Continue?",
        ):
            return
        removed = remove_step(ech, etape)
        self._afficher(f"{ech.id}: {removed} measurement(s) removed ({len(ech.mesures)} remaining).\n")

    def ouvrir_remove_bad_quality_dialog(self):
        """Detecte automatiquement, pour TOUS les echantillons de la
        selection courante, les etapes marquees qualite 'b' (colonne
        "quality" du .prmag - voir testlect.Measurement.quality) et les
        supprime via datatools.remove_bad_quality_steps (meme
        renumerotation des cod2 R/V/P que "Remove paleointensity step",
        meme suppression en cascade d'une eventuelle verification pTRM
        qui se referait a l'etape rejetee) - demande explicite
        utilisateur ("je voudrais utiliser le critere de qualite b/g
        dans les donnees pour ne pas prendre en compte cette etape...
        un des problemes en paleointensite est eventuellement qu'une
        serie d'echantillons ne soient pas mis correctement dans le
        four... il faut alors renumeroter les etapes automatiquement et
        aussi supprimer la PTRM check si elle se refere a l'etape
        rejetee"). Traite TOUTE la selection en un coup (le cas motivant
        est justement une SERIE d'echantillons affectee par le meme
        probleme de four), pas un seul echantillon a la fois."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        preview = []
        for ech in self.selection:
            bad_etapes = sorted({
                m.etape for m in ech.mesures if (m.quality or "g").strip().lower() == "b"
            })
            if bad_etapes:
                preview.append((ech, bad_etapes))
        if not preview:
            self._showinfo(
                "No bad-quality steps",
                "No measurement marked quality 'b' found in the current selection.")
            return
        lines = "\n".join(
            f"  {ech.id}: step(s) {', '.join(f'{e:g}' for e in bad_etapes)}"
            for ech, bad_etapes in preview)
        if not messagebox.askyesno(
            "Remove bad-quality steps",
            "The following step(s) are marked quality 'b' and will be removed "
            "(steps renumbered, any pTRM check referencing them also removed) - "
            f"irreversible for this session:\n\n{lines}\n\nContinue?",
        ):
            return
        total_removed = 0
        for ech, _bad_etapes in preview:
            removed, bad_etapes = remove_bad_quality_steps(ech)
            total_removed += removed
            self._afficher(
                f"{ech.id}: {removed} measurement(s) removed (bad-quality step(s) "
                f"{', '.join(f'{e:g}' for e in bad_etapes)}), {len(ech.mesures)} remaining.\n"
            )
        self._afficher(f"Total: {total_removed} measurement(s) removed across {len(preview)} sample(s).\n")

    def ouvrir_detect_grm_dialog(self):
        """Diagnostic (PAS de correction - voir "Suppress GRM" juste en
        dessous, qui ne traite qu'un cas different, voir datatools.py) de
        contamination par GRM (gyroremanent magnetization) sur un
        degausser 3 axes EN LIGNE avec le magnetometre - demande explicite
        utilisateur ("Most laboratories using 3 axis degausser online with
        the Cryogenic magnetometer... It is often not recognized by users
        or student... Can you write a test to detect such behavior
        especially for the Magic database contribution", puis "not during
        export but as an independant menu, perhaps above suppress GRM").
        N'empeche rien (pas de blocage d'export) - un signal a inspecter
        visuellement au Zijderveld, voir datatools.detect_grm pour la
        methode."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        results = [detect_grm(ech) for ech in self.selection]
        usable = [r for r in results if r is not None]
        if not usable:
            self._afficher(
                "Detect GRM: no specimen in the selection has enough AF steps on a "
                "3-axis in-line degausser (instrument code starting with 'C') to run this test.\n"
            )
            return
        flagged = [r for r in usable if r.suspected]
        lines = "\n".join(
            f"  {r.specimen}: perp. component reaches {r.max_ratio * 100:.0f}% of NRM at "
            f"highest AF step (correlation with field: {r.correlation:.2f}, "
            f"{r.mean_perp_fraction * 100:.0f}% of the deviation is orthogonal to the coil axis)"
            for r in sorted(flagged, key=lambda r: -r.max_ratio)
        )
        if flagged:
            self._afficher(
                f"Detect GRM: possible GRM contamination on {len(flagged)}/{len(usable)} "
                f"specimen(s) checked (AF demagnetization on a 3-axis in-line degausser, "
                f"component orthogonal to the last coil axis growing with AF field) - "
                f"inspect the Zijderveld plot before trusting these directions:\n{lines}\n"
            )
        else:
            self._afficher(f"Detect GRM: no contamination suspected on the {len(usable)} specimen(s) checked.\n")

    def ouvrir_elimine_grm_dialog(self):
        """Equivalent GUI de `elimineGRM` ("Suppress GRM") : réduit les
        triplets de mesures GRM (cod1='F', cod2 'X','Y','Z' consécutifs)
        en un seul point corrigé, méthode au choix (1 = substitution
        axe-par-axe, x du point X / y du point Y / z du point Z ; 2 =
        moyenne simple des 3 points)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        method_s = self._console_input("Method 1 (axis-by-axis substitution) or 2 (average): ", "1")
        if method_s is None:
            return
        try:
            method = int(method_s)
        except ValueError:
            method = 1
        if method not in (1, 2):
            method = 1
        if not messagebox.askyesno(
            "Suppress GRM",
            f"Reduce the GRM triplets (X/Y/Z) of {len(self.selection)} sample(s) "
            "with method " + ("axis-by-axis substitution" if method == 1 else "average")
            + " - irreversible for this session. Continue?",
        ):
            return
        for ech in self.selection:
            n = eliminate_grm(ech, method=method)
            self._afficher(f"{ech.id}: {n} triplet(s) reduced ({len(ech.mesures)} measurements remaining).\n")

    def ouvrir_convzmoins_dialog(self):
        """Equivalent GUI de `convzmoins` ("ConvertZ- 2G") : inverse y/z et
        recode en 'R' les mesures cod1='Z' cod2='-'."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        if not messagebox.askyesno(
            "ConvertZ- 2G",
            f"Convert Z- measurements to R for {len(self.selection)} sample(s) "
            "- irreversible for this session. Continue?",
        ):
            return
        total = 0
        for ech in self.selection:
            n = convert_z_minus(ech)
            total += n
            if n:
                self._afficher(f"{ech.id}: {n} Z- measurement(s) converted.\n")
        if not total:
            self._afficher("No Z- measurement found in the selection.\n")

    def ouvrir_exportthellier_dialog(self):
        """Equivalent GUI de `exportthellier` ("export to ThellierTool") :
        écrit un fichier .tdt par échantillon sélectionné, dans un dossier
        choisi par l'utilisateur."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        out_dir = os.path.dirname(self.results_path) if self.results_path else os.getcwd()
        chosen = filedialog.askdirectory(initialdir=out_dir, title="ThellierTool export folder (.tdt)")
        if not chosen:
            return
        written = []
        for ech in self.selection:
            try:
                path = export_thellier_tdt(ech, chosen)
                written.append(path)
            except OSError as e:
                self._showerror("Error", f"Export failed for {ech.id}:\n{e}")
                return
        self._afficher("Exported .tdt files:\n" + "\n".join(written) + "\n")

    # ------------------------------------------------------------------
    # Aides communes
    # ------------------------------------------------------------------

    def _restore_focus(self):
        """Apres un popup natif (messagebox), le focus clavier ne revient
        PAS automatiquement sur la fenetre principale (constate sur macOS/
        Aqua) - les raccourcis lies via bind_all restent inertes tant qu'on
        n'a pas clique manuellement dans une fenetre de l'appli - demande
        explicite utilisateur ("we need to manually click in one of the
        window before to activate the shortcuts. Is there a way to change
        this behavior?"). Force le focus sur la fenetre principale, et sur
        la console texte (cible naturelle de la frappe/des raccourcis)."""
        self.root.lift()
        self.root.focus_force()
        self.text_area.focus_set()

    def _apply_anisotropy_correction(self, ech, ani_path, direction, h_raw):
        """Equivalent de plotpaleoint2.f:1253-1263 : `anicor` (correction
        scalaire, voir calcul.compute_anicor_factor) n'est appele QUE si le
        specimen n'a PAS deja ete corrige via "Inverse_ANI_correction..."
        (`ech.flaganiso` - qui corrige VECTORIELLEMENT chaque mesure brute,
        voir calcul.apply_inverse_anisotropy) - `if(ech(inh).flaganiso .eq.
        .false.)` cote Fortran. Le port precedent ignorait ce garde-fou et
        appliquait donc systematiquement la correction scalaire, y compris
        sur des mesures DEJA corrigees vectoriellement - une double
        correction silencieuse - demande explicite utilisateur ("it is
        unclear for the user that a correction for anisotropy is done...
        the whole list of measurement were corrected by the inverse
        anisotropy tensor and the anisotropy corrected paleointensity was
        automatically calculated from the same temperature interval. This
        is not obvious in the text output").

        Retourne (fcor, hcorani, note, tensor) : `note` est une ligne de
        texte EXPLICITE (absente du Fortran, qui se contentait de remettre
        fcor/tensor a 0/"none" sans rien signaler de plus clair) indiquant
        sans ambiguite SI et COMMENT H a ete corrige pour l'anisotropie -
        a afficher juste apres le resultat principal. `tensor` (calcul.
        AniTensor, le tenseur A0 BRUT utilise, ou None si aucune
        correction) est renvoye pour l'archivage .pmagint (voir
        paleointensity.write_pmagint_line) - evite une seconde lecture du
        fichier .pmagani pour la meme information."""
        if getattr(ech, "flaganiso", False):
            return None, None, (
                "Anisotropy: measurements already corrected (Inverse_ANI_correction "
                "applied) - H below is the anisotropy-corrected intensity.\n"
            ), None
        if not ani_path or direction.anchored_specimen_frame is None:
            return None, None, "", None
        tensor = read_ani_tensor(ani_path, ech.id, "A0")
        if tensor is None:
            return None, None, "", None

        if tensor.quality == "b":
            # demande explicite utilisateur ("in the paleointensity
            # correction: ask : [id] anisotropy NOT significant (not
            # satisfactory) - inverse anisotropy correction anyway ?") :
            # le verdict PmagPy (Hext F-test, persiste dans .pmagani -
            # voir calcul._PMAGANI_HEADER "quality") indique que ce
            # tenseur n'est pas distinguable du bruit de mesure - corriger
            # quand meme reste possible, mais ne doit plus se faire
            # silencieusement.
            confirm = self._console_input(
                f"[{ech.id}] anisotropy NOT significant (not satisfactory).\n"
                f"Inverse anisotropy correction anyway? y/N: ", "N")
            if confirm is None or confirm.strip().lower() != "y":
                return None, None, (
                    "Anisotropy: NOT applied (Hext F-test not satisfactory, declined) - "
                    "H below is the raw (uncorrected) intensity.\n"
                ), None

        fcor = compute_anicor_factor(tensor, direction.anchored_specimen_frame)
        hcorani = h_raw * fcor
        note = (
            "Anisotropy correction applied (tensor 'A0'"
            + (", NOT significant - applied anyway" if tensor.quality == "b" else "")
            + f"): fcor={fcor:.3f}  H (raw)={h_raw:.2f}µT -> Hcorani (corrected)={hcorani:.2f}µT\n"
        )
        return fcor, hcorani, note, tensor

    def _showinfo(self, *args, **kwargs):
        result = messagebox.showinfo(*args, **kwargs)
        self._restore_focus()
        return result

    def _showwarning(self, *args, **kwargs):
        result = messagebox.showwarning(*args, **kwargs)
        self._restore_focus()
        return result

    def _showerror(self, *args, **kwargs):
        result = messagebox.showerror(*args, **kwargs)
        self._restore_focus()
        return result

    def _not_implemented(self, name):
        """Entree de menu presente dans le Fortran d'origine (StarmacOSX_x.f95,
        menu Calcul) mais pas encore portee en Python."""
        self._showinfo("Not yet implemented", f"« {name} » has not been ported yet.")

    def _archive_only(self, fit):
        """Equivalent de `call archivres` seul (sans toucher self.results) -
        pour les flux qui ajoutent `fit` a self.results en amont a titre
        d'apercu (ex. superposition sur le Zijderveld avant confirmation),
        et n'archivent qu'une fois la sauvegarde confirmee par l'utilisateur."""
        if not self.results_path:
            return
        try:
            _, self._archived_ids = archivres(fit, self.results_path, self._archived_ids)
            # Met a jour la date connue APRES notre propre ecriture - sinon
            # le prochain "Select results..." croirait le fichier modifie
            # PAR L'UTILISATEUR (voir ouvrir_selres_dialog) et viderait
            # self.results a tort, alors que ce changement vient de
            # l'appli elle-meme.
            self._results_path_mtime = os.path.getmtime(self.results_path)
        except OSError as e:
            self._showwarning(
                "Archiving failed",
                f"The result is in memory but could not be written to "
                f"{self.results_path}:\n{e}",
            )

    def _save_result(self, fit):
        """Equivalent de `tr(nbres+1)=res; nbres=nbres+1; call archivres` :
        ajoute `fit` a self.results ET l'archive immediatement dans le
        fichier .r (equivalent filr) - le Fortran appelle archivres des
        qu'un resultat est sauvegarde, pas seulement en fin de session."""
        self.results.append(fit)
        self._archive_only(fit)

    def _save_results(self, fits):
        for fit in fits:
            self._save_result(fit)

    def _find_ani_path(self):
        """Retrouve le fichier .pmagani associe (voir calcul.ani_path_for) -
        essaie d'abord directement depuis self.results_path (le cas normal :
        .prmag et .pmagres/.pmagani partagent le meme nom de base), puis via
        results_path_for(self.results_path) au cas ou self.results_path
        porterait encore le nom du fichier de DONNEES plutot que du fichier
        de resultats (ancien comportement, garde par securite). Si aucun
        .pmagani n'existe, retombe sur un ANCIEN .ANI de meme nom de base
        (retro-compatibilite - "can we import old style .ANI in these
        pmagani style" : les fichiers .ANI pas encore convertis restent
        utilisables directement, voir calcul.read_ani_tensor qui dispatche
        deja par extension). None si rien n'existe."""
        if not self.results_path:
            return None
        for base in (self.results_path, results_path_for(self.results_path)):
            candidate = ani_path_for(base)
            if os.path.exists(candidate):
                return candidate
            legacy = os.path.splitext(base)[0] + ".ANI"
            if os.path.exists(legacy):
                return legacy
        return None

    def _afficher(self, text):
        """Ajoute `text` a la suite du contenu existant (n'efface plus la
        console a chaque appel de menu), et fait defiler jusqu'au nouveau
        contenu. Les segments entoures de HEADER_MARK (lignes de titres de
        colonnes des fonctions de rapport - list_measurements, list_results,
        etc.) sont inseres avec le tag "header" (gras) plutot qu'en texte
        normal, marqueurs retires - demande explicite utilisateur
        ("throughout the software, is it possible to write the header in
        bold"). Meme principe pour WARN_MARK/tag "warn" (gras, rouge, fond
        surligne) - demande explicite utilisateur ("write a strong
        warning", voir convert_magic_to_r.scan_dip_sign). Chaque type de
        marqueur bascule INDEPENDAMMENT son propre tag (les deux ne sont
        jamais imbriques dans les appelants actuels, mais rien ne
        l'empeche)."""
        if self.text_area.get("1.0", "end-1c").strip():
            self.text_area.insert(tk.END, "\n" + "-" * 60 + "\n")
        mark_tags = {HEADER_MARK: "header", WARN_MARK: "warn"}
        active_tag = {HEADER_MARK: False, WARN_MARK: False}
        for part in re.split(f"([{HEADER_MARK}{WARN_MARK}])", text):
            if part in mark_tags:
                active_tag[part] = not active_tag[part]
                continue
            if not part:
                continue
            tag = next((mark_tags[m] for m, on in active_tag.items() if on), None)
            if tag:
                self.text_area.insert(tk.END, part, tag)
            else:
                self.text_area.insert(tk.END, part)
        self.text_area.see(tk.END)

    @staticmethod
    def _parse_demag(code):
        code = (code or "").strip()
        if not code or code == "*":
            return "*", "*"
        if len(code) == 1:
            return code[0], "*"
        return code[0], code[1]

    def _read_step_range(self, demag1, default_min="0", default_max="9999"):
        """Demande Step min/Step max, avec un libelle adapte a `demag1`
        quand il est deja connu : mT pour un pas AF (A/F), degC pour un
        pas thermique (D/S/T/K, et R/V/P - paleointensite Thellier/IZZI,
        demande explicite utilisateur "when the cod1 is R, V, P, the unit
        is also °C"), step brut sinon (code inconnu, '*', ou `demag1=None`
        - l'appelant demande alors le code demag APRES le step range,
        comme les anciennes applis, voir ouvrir_selection_dialog :
        demande explicite utilisateur "can you ask for the steps before
        the demag code. (as done in the old apps; I get confused
        otherwise)"). `Measurement.etape` est desormais directement la
        valeur physique reelle (plus d'echelle Oersted-equivalente a
        compenser ici - demande explicite utilisateur "convert all step
        integer to float") : ce que l'utilisateur tape est compare tel
        quel a `etape`, le libelle sert seulement a clarifier l'unite
        attendue. Retourne (step_min, step_max) ou None si annule/invalide
        (message d'erreur deja affiche dans ce cas)."""
        is_af = demag1 in ("A", "F")
        is_thermal = demag1 in ("D", "S", "T", "K", "R", "V", "P")
        unit_label = "mT" if is_af else ("degC" if is_thermal else "step")
        step_min_s = self._console_input(f"Step min ({unit_label}): ", default_min)
        if step_min_s is None:
            return None
        step_max_s = self._console_input(f"Step max ({unit_label}): ", default_max)
        if step_max_s is None:
            return None
        try:
            step_min_val = float(step_min_s or default_min)
            step_max_val = float(step_max_s or default_max)
        except ValueError:
            self._showerror("Error", "Step min / Step max must be numbers.")
            return None
        return step_min_val, step_max_val

    def _console_input(self, prompt, default=""):
        """Lit une ligne tapee au clavier DANS la fenetre texte, comme une
        vraie console (equivalent du `write(*,...)` + `read(*,*)` du Fortran
        d'origine) : affiche `prompt`, laisse l'utilisateur taper, puis rend
        la main des qu'il appuie sur Entree. Renvoie None si annule (Echap).

        `default` n'est PAS pre-insere dans le texte (essaye au debut, mais
        source de confusion : si l'utilisateur tapait sans d'abord effacer
        la valeur pre-remplie, elle s'affichait concatenee avec sa saisie,
        ex. "dm" pour un defaut "d" suivi d'un "m" tape) - applique
        uniquement si la ligne validee est vide (Entree seule).

        Bloque l'appelant sans geler l'interface : `wait_variable` continue
        de faire tourner la boucle d'evenements Tk pendant l'attente (pas de
        thread ni de callback a gerer cote appelant).

        Le prompt commence TOUJOURS sur une ligne neuve - demande explicite
        utilisateur ("there are still long lines text output and expecting
        an answer at the end the line... display the question (in red) in
        a new line") : sans ce saut de ligne, un prompt arrivant juste
        apres un texte deja affiche SANS retour a la ligne final (ex. le
        dernier chiffre d'un rapport) se retrouvait concatene a sa suite
        sur la meme ligne, illisible. La couleur rouge existe deja (tag
        "prompt", voir tag_configure) - seul le retour a la ligne manquait
        ici, point d'entree UNIQUE de tous les prompts interactifs de
        l'appli (corrige une fois pour toutes plutot qu'a chaque site
        d'appel individuellement)."""
        last_char = self.text_area.get("end-2c", "end-1c")
        if last_char not in ("", "\n"):
            self.text_area.insert(tk.END, "\n")
        self.text_area.insert(tk.END, prompt, "prompt")
        start_index = self.text_area.index("end-1c")
        self.text_area.mark_set(tk.INSERT, tk.END)
        self.text_area.see(tk.END)
        self.text_area.focus_set()

        outcome = {"value": None}
        done = tk.BooleanVar(value=False)

        def on_return(event):
            typed = self.text_area.get(start_index, "end-1c")
            outcome["value"] = typed if typed.strip() else default
            self.text_area.insert(tk.END, "\n")
            done.set(True)
            return "break"

        def on_escape(event):
            self.text_area.delete(start_index, "end-1c")
            outcome["value"] = None
            self.text_area.insert(tk.END, "\n")
            done.set(True)
            return "break"

        ret_id = self.text_area.bind("<Return>", on_return)
        esc_id = self.text_area.bind("<Escape>", on_escape)
        self.text_area.wait_variable(done)
        self.text_area.unbind("<Return>", ret_id)
        self.text_area.unbind("<Escape>", esc_id)
        self.text_area.see(tk.END)
        return outcome["value"]

    def _lister(self, func):
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        buffer = io.StringIO()
        func(self.selection, orientation=self.orientation.get(), out=buffer)
        self._afficher(buffer.getvalue())

    # ------------------------------------------------------------------
    # Selection donnees : selmes / selentete / effmes / initmes
    # ------------------------------------------------------------------

    def ouvrir_selection_dialog(self):
        """Equivalent GUI de `selmes` - saisie au clavier DANS la fenetre
        texte (console), pas de fenetre popup - questions posees les unes
        apres les autres comme le faisait le Fortran d'origine."""
        if not self.donnees:
            self._showwarning(
                "No data",
                "Load a .ren file or import a MagIC folder first.",
            )
            return

        self.text_area.insert(tk.END, "\n--- Data selection (Escape to cancel) ---\n", "prompt")
        pattern = self._read_prefixed_pattern("Sample number (* = all, ? = wildcard): ")
        if pattern is None:
            return
        # Step avant Demag code (ordre des anciennes applis - demande
        # explicite utilisateur, "I get confused otherwise") : le libelle
        # ne peut alors plus s'adapter a l'unite (mT/degC), demag1 pas
        # encore connu a ce point - `_read_step_range(None)` retombe sur
        # le libelle generique "step" (voir sa docstring), comme le fait
        # deja `select_samples_interactive` (selection.py, equivalent CLI
        # de selmes) avec son " Step min value:" generique.
        step_range = self._read_step_range(None)
        if step_range is None:
            return
        step_min, step_max = step_range
        demag_s = self._console_input("Demag code (e.g. N0, T, AF, * = all): ", "*")
        if demag_s is None:
            return
        demag1, demag2 = self._parse_demag(demag_s)
        new_matches = select_samples(
            self.donnees,
            pattern=pattern or "*",
            step_min=step_min,
            step_max=step_max,
            demag1=demag1,
            demag2=demag2,
            verbose=False,
        )
        # equivalent de selmes (dataselect.f) : accumule sur la selection
        # existante (nbech n'est jamais remis a zero dans selmes lui-meme,
        # seul initmes/"Init list" le fait) - ne remplace pas self.selection.
        self.selection = self.selection + new_matches
        nb_mesures = sum(len(s.mesures) for s in self.selection)
        self.text_area.insert(
            tk.END, f"Selection: +{len(new_matches)} sample(s) - "
                    f"total {len(self.selection)} sample(s), {nb_mesures} measurement(s)\n")
        self.text_area.see(tk.END)

    def _pick_from_list(self, title, items, allow_all=True):
        """Fenetre modale avec une Listbox pour choisir UN element a la
        souris (double-clic ou bouton Select) - PAS dans le Fortran
        (console texte uniquement), demande explicite de l'utilisateur
        pour remplacer la saisie au clavier du nom de site. Retourne
        l'element choisi, '*' si `allow_all` et l'entree "* (all)" est
        choisie, ou None si annule (bouton Cancel / fermeture)."""
        top = tk.Toplevel(self.root)
        top.title(title)
        top.transient(self.root)
        top.grab_set()

        listbox = tk.Listbox(top, activestyle="dotbox", exportselection=False)
        if allow_all:
            listbox.insert(tk.END, "* (all)")
        for item in items:
            listbox.insert(tk.END, item)
        listbox.selection_set(0)
        listbox.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        result = {"value": None}

        def confirm(_event=None):
            sel = listbox.curselection()
            if sel:
                text = listbox.get(sel[0])
                result["value"] = "*" if allow_all and text == "* (all)" else text
            top.destroy()

        def cancel():
            top.destroy()

        listbox.bind("<Double-Button-1>", confirm)
        listbox.bind("<Return>", confirm)

        btn_frame = ttk.Frame(top)
        btn_frame.pack(fill=tk.X, padx=8, pady=(0, 8))
        ttk.Button(btn_frame, text="Select", command=confirm).pack(side=tk.RIGHT, padx=(4, 0))
        ttk.Button(btn_frame, text="Cancel", command=cancel).pack(side=tk.RIGHT)

        top.protocol("WM_DELETE_WINDOW", cancel)
        listbox.focus_set()
        self.root.wait_window(top)
        return result["value"]

    def ouvrir_selection_site_dialog(self):
        """Pas dans le Fortran (aucune selection par site n'existe dans le
        menu d'origine) : selectionne tous les echantillons d'un site
        MagIC (champ `magic_site`, decode depuis la ligne roche a
        l'ouverture du fichier - voir testlect.decode_roche). Accumule sur
        self.selection comme `selmes` (voir ouvrir_selection_dialog),
        plutot que de la remplacer."""
        if not self.donnees:
            self._showwarning(
                "No data",
                "Load a .ren file or import a MagIC file first.",
            )
            return

        sites = sorted({p.magic_site.strip() for p in self.donnees if p.magic_site.strip()})
        if not sites:
            self._showwarning(
                "No site", "No sample has a decoded MagIC site (see the « roche » line).")
            return

        self.text_area.insert(tk.END, "\n--- Selection by site (Escape to cancel) ---\n", "prompt")
        site = self._pick_from_list("Select site", sites)
        if site is None:
            return
        # Step avant Demag code - voir ouvrir_selection_dialog.
        step_range = self._read_step_range(None)
        if step_range is None:
            return
        step_min, step_max = step_range
        demag_s = self._console_input("Demag code (e.g. N0, T, AF, * = all): ", "*")
        if demag_s is None:
            return
        demag1, demag2 = self._parse_demag(demag_s)
        new_matches = select_samples_by_site(
            self.donnees, site=site or "*", step_min=step_min, step_max=step_max,
            demag1=demag1, demag2=demag2, verbose=False,
        )
        self.selection = self.selection + new_matches
        nb_mesures = sum(len(s.mesures) for s in self.selection)
        self._afficher(
            f"Site selection « {site} »: +{len(new_matches)} sample(s) - "
            f"total {len(self.selection)} sample(s), {nb_mesures} measurement(s)"
        )

    def ouvrir_entete_dialog(self):
        """Equivalent GUI de `selentete` : une saisie VIDE efface le
        header (lentete=0 dans le Fortran - dataselect.f:1295, blanc ==
        pas de prefixe actif), plutot que de re-appliquer l'ancienne
        valeur. Bug reel corrige : `self.entete` etait utilise comme
        valeur par defaut, rendant le header impossible a reinitialiser
        (Entree seul renvoyait toujours l'ancienne valeur)."""
        current = f" (current: '{self.entete}')" if self.entete else ""
        entete = self._console_input(f"Header{current}, blank to clear: ", "")
        if entete is None:
            return
        self.entete = entete.upper()[:12]
        self._afficher(
            f"Header set: '{self.entete}'\n" if self.entete else "Header cleared.\n"
        )

    def _read_prefixed_pattern(self, prompt):
        """Lit un motif echantillon en appliquant le header actif comme
        VRAI prefixe, fidele a `numero=entete(1:lentete)//chaine(1:len(chaine))`
        (dataselect.f, `selmes`/`selres` branche 'd') : le header est
        concatene devant TOUT ce qui est tape, y compris une saisie vide -
        dans ce cas le motif final est le header seul, pas '*'. Avant ce
        correctif, `self.entete` n'etait utilise que comme valeur par
        defaut (uniquement si la saisie etait vide), donc taper un motif
        n'appliquait jamais le prefixe : bug rapporte par l'utilisateur
        (header 'TRO' sans effet, obligeant a taper 'TRO14' en entier)."""
        label = f"{prompt}{self.entete}" if self.entete else prompt
        raw = self._console_input(label, "")
        if raw is None:
            return None
        if self.entete:
            return self.entete + raw
        return raw if raw.strip() else "*"

    def ouvrir_effmes_dialog(self):
        """Equivalent GUI de `effmes`."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        self.text_area.insert(tk.END, "\n--- Delete measurements (Escape to cancel) ---\n", "prompt")
        pattern = self._console_input("Sample to erase (* = all): ", "*")
        if pattern is None:
            return
        # Step avant Demag code - voir ouvrir_selection_dialog.
        step_range = self._read_step_range(None, default_max="9000")
        if step_range is None:
            return
        step_min, step_max = step_range
        demag_s = self._console_input("Demag code (* = all): ", "*")
        if demag_s is None:
            return
        demag1, demag2 = self._parse_demag(demag_s)
        occurrence = self._console_input(
            "Occurrence to delete (* = all matches, 1 = first only, 2 = second only...): ", "*")
        if occurrence is None:
            return

        self.selection = delete_measurements(
            self.selection,
            pattern=pattern or "*",
            step_min=step_min,
            step_max=step_max,
            demag1=demag1,
            demag2=demag2,
            occurrence=occurrence or "*",
            verbose=False,
        )
        nb_mesures = sum(len(s.mesures) for s in self.selection)
        self.text_area.insert(
            tk.END, f"Selection after deletion: {len(self.selection)} sample(s), {nb_mesures} measurement(s)\n")
        self.text_area.see(tk.END)

    def reinitialiser_selection(self):
        """Equivalent GUI de `initmes`."""
        self.selection = init_selection()
        self._afficher("Selection reset - no sample selected.\n")

    # ------------------------------------------------------------------
    # Selection donnees : lismes / listeXYZ / lismesVRM / lismesdepth / infoech
    # ------------------------------------------------------------------

    def lister_mesures(self):
        self._lister(list_measurements)

    def lister_xyz(self):
        self._lister(list_xyz)

    def lister_vrm(self):
        self._lister(list_measurements_vrm)

    def ouvrir_lismesdepth_dialog(self):
        """Equivalent GUI de `lismesdepth`. Utilise en priorite le champ
        `stratigraphic_height` du .prmag (voir testlect.Pmag.
        stratigraphic_height) quand au moins un specimen de la selection
        l'a renseigne - demande explicite utilisateur ("add an additional
        variable: stratigraphic_position (or the Magic equivalent)... this
        field will replace the need to load a file for magnetostratigraphic
        studies") : plus besoin de fichier externe specimen/depth des que
        cette information est deja dans le fichier de donnees lui-meme.
        Ne demande un fichier externe (ancien comportement) que si AUCUN
        specimen de la selection n'a de position renseignee - garde une
        voie de secours pour les fichiers plus anciens qui n'ont pas encore
        ce champ."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        depths = {
            ech.id: ech.stratigraphic_height
            for ech in self.selection if ech.stratigraphic_height is not None
        }
        if depths:
            n_missing = len(self.selection) - len(depths)
            self._afficher(
                f"Using stratigraphic_height from the .prmag file "
                f"({len(depths)} specimen(s) with a known position"
                + (f", {n_missing} without - shown as -999.0" if n_missing else "")
                + ").\n"
            )
            if n_missing:
                self._afficher(
                    "Tip: fill in the missing position(s) with PmagFile -> "
                    "Complete sample information... (specimen mode, "
                    "'stratigraphic_height' column) so they show up here too.\n"
                )
        else:
            self._afficher(
                "Tip: instead of a one-off depth file, you can fill in "
                "'stratigraphic_height' directly in the .prmag with PmagFile -> "
                "Complete sample information... - it will then be picked up "
                "automatically here (and in export detailed/Latex) without "
                "asking for a file again.\n"
            )
            depth_path = filedialog.askopenfilename(
                title="No stratigraphic_height in this file - depth file (2 columns: sample_name  depth)",
                filetypes=[("Text", "*.txt"), ("All files", "*.*")],
            )
            if not depth_path:
                return
            try:
                with open(depth_path, "r", encoding="utf-8") as f:
                    for line in f:
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                depths[parts[0]] = float(parts[1])
                            except ValueError:
                                continue
            except Exception as e:
                self._showerror("Error", f"Could not read the depth file:\n{e}")
                return

        self.text_area.insert(tk.END, "\n--- Expected site direction (Escape to cancel) ---\n", "prompt")
        dec_s = self._console_input("Expected declination (D): ", "0.0")
        if dec_s is None:
            return
        inc_s = self._console_input("Expected inclination (I): ", "0.0")
        if inc_s is None:
            return
        try:
            expected_dec = float(dec_s or 0.0)
            expected_inc = float(inc_s or 0.0)
        except ValueError:
            self._showerror("Error", "Declination / Inclination must be numbers.")
            return

        buffer = io.StringIO()
        list_measurements_depth(
            self.selection, depths,
            expected_dec=expected_dec, expected_inc=expected_inc,
            orientation=self.orientation.get(), out=buffer,
        )
        self._afficher(buffer.getvalue())

    def afficher_info_echantillons(self):
        """Equivalent GUI de `infoech`."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        self._afficher(sample_info(self.selection))

    # ------------------------------------------------------------------
    # Calcul : ajuslig / fishmes / fishres / lisres / initres
    # ------------------------------------------------------------------

    def ouvrir_ajuslig_dialog(self):
        """Equivalent GUI de `ajuslig` (ajustement de droite par ACP,
        linesplans.f) : boucle sur TOUS les echantillons de la selection,
        comme le `do i=1,nbech` d'origine - zijderveld affiche pour chaque
        echantillon, premier step a 0 pour passer au suivant sans ajuster,
        puis apres chaque ajustement : sauver/refaire/suivant (Y/r/n), au
        lieu de traiter un seul echantillon choisi a l'avance."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        samples = [s for s in self.selection if len(s.mesures) >= 2]
        if not samples:
            self._showwarning(
                "Not enough measurements", "No selected sample has at least 2 measurements.")
            return

        for ech in samples:
            self._current_graphic = ("zijderveld", ech.id)
            self._refresh_current_graphic()

            steps_text = "  ".join(
                f"{i + 1}:{_step_token(m)}" for i, m in enumerate(ech.mesures))
            self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")

            while True:  # boucle "refaire" (r)
                jdeb_s = self._console_input(
                    f"[{ech.id}] First step (0 = skip this sample): ", "1")
                if jdeb_s is None:
                    return
                try:
                    jdeb = int(jdeb_s)
                except ValueError:
                    self._showerror("Error", "First step must be an integer.")
                    continue
                if jdeb == 0:
                    break  # echantillon suivant, sans ajustement

                jfin_s = self._console_input("Last step: ", str(len(ech.mesures)))
                if jfin_s is None:
                    return
                ancr_s = self._console_input("Anchored to origin (Y/n): ", "Y")
                if ancr_s is None:
                    return
                component = self._console_input(
                    "Magnetization component (A/B/C...): ", "A")
                if component is None:
                    return
                try:
                    jfin = int(jfin_s)
                except ValueError:
                    self._showerror("Error", "Last step must be an integer.")
                    continue
                if not (1 <= jdeb <= jfin <= len(ech.mesures)):
                    self._showerror("Error", "Invalid step range.")
                    continue
                anchored = ancr_s.strip().lower() != "n"

                fit = fit_line(ech, jdeb, jfin, anchored=anchored)
                if fit is None:
                    self._showwarning(
                        "Fit rejected",
                        "MAD > 15° or non-linear trend (linearity test failed).",
                    )
                    continue
                fit.component = component.strip().upper() or "A"

                dec, inc = _correct_dec_inc(fit, self.orientation.get())
                self._afficher(
                    f"{ech.id}: dec={dec:.1f}  inc={inc:.1f}  mad={fit.mad:.1f}  "
                    f"nb points={fit.nb}  ({'anchored' if anchored else 'not anchored'})\n"
                )
                # affiche le zijderveld AVEC cet ajustement superpose (pas
                # encore sauvegarde - ajoute temporairement a self.results,
                # que _refresh_current_graphic utilise deja pour filtrer les
                # ajustements a superposer par echantillon)
                self.results.append(fit)
                self._current_graphic = ("zijderveld", ech.id)
                self._refresh_current_graphic()

                decision = self._console_input(
                    "Save (Y) / redo (r) / next sample (n): ", "Y")
                if decision is None:
                    self.results.remove(fit)
                    return
                decision = decision.strip().lower()
                if decision == "r":
                    self.results.remove(fit)
                    continue  # redo pour le meme echantillon
                if decision != "n":
                    self._archive_only(fit)
                    self._afficher(f"Result saved: {len(self.results)}\n")
                else:
                    self.results.remove(fit)
                break  # echantillon suivant

    def ouvrir_autointerpretation_dialog(self):
        """Pas dans le Fortran (aucun equivalent) - demande explicite
        utilisateur ("une routine intelligente qui ferait des
        interpretations des diagrammes de desaimantation"), deuxieme des
        deux routines demandees (voir auto_interpretation.py pour
        l'algorithme - PCA en extension gloutonne, PAS un modele IA/ML).

        Boucle sur la selection, propose jusqu'a 2 composantes par
        echantillon (primary/secondary). Les DEUX sont converties en
        FitResult et superposees ENSEMBLE sur le Zijderveld AVANT toute
        question de sauvegarde - demande explicite utilisateur ("is it
        possible to visualize the interpretation on the zijderveld before
        to accept it") - puis seule(s) celle(s) choisie(s) sont
        conservees/archivees, les autres retirees de self.results.
        AUCUNE suggestion n'est archivee automatiquement (limite assumee
        et documentee dans auto_interpretation.py : fiable sur une
        decroissance a composante unique, faillible sur des composantes
        qui se chevauchent - chaque suggestion reste a valider).

        Regroupe d'abord par site MagIC (`magic_site`) avant d'appeler
        propose_components_for_site - demande explicite utilisateur
        ("you miss the point at a site level. if some samples have a
        well defined secondary magnetization in a wide temperature
        range... it is likely that this behavior is the same for all
        samples"). N'echange JAMAIS automatiquement primary/secondary
        (voir auto_interpretation.py, docstring de
        propose_components_for_site, pour deux contre-exemples reels -
        14NQ04 puis 14NQ0403B - qui ont ecarte cette idee) : se contente
        d'ajouter une note des deux cotes quand le "primary" isole d'un
        specimen s'ecarte du palier typique des autres specimens du site,
        pour que l'utilisateur compare lui-meme les deux directions sur
        le Zijderveld avant de choisir.
        Les specimens sans site decode (magic_site vide) restent traites
        individuellement - un site vide melangerait des specimens sans
        rapport dans un meme "consensus"."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        viable = [ech for ech in self.selection if len(ech.mesures) >= 2]

        by_site = {}
        no_site = []
        for ech in viable:
            site = (ech.magic_site or "").strip()
            if site:
                by_site.setdefault(site, []).append(ech)
            else:
                no_site.append(ech)

        suggestions_by_id = {}
        for echs in by_site.values():
            suggestions_by_id.update(propose_components_for_site(echs))
        for ech in no_site:
            suggestions_by_id[ech.id] = propose_components(ech)

        self.text_area.insert(tk.END, "\n--- Auto-interpretation (suggestions, Escape to stop) ---\n", "prompt")
        for ech in viable:
            suggestions = suggestions_by_id[ech.id]
            if not suggestions:
                self._afficher(format_suggestions(ech.id, suggestions))
                continue

            # Construit un FitResult PROVISOIRE par suggestion et les
            # superpose TOUS ensemble sur le Zijderveld avant de demander
            # quoi que ce soit (numcomp=1/2 pour que les deux se
            # distinguent visuellement, meme convention que ajuslig).
            candidates = {}  # label -> FitResult
            for s in suggestions:
                jdeb = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == s.step_first), None)
                jfin = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == s.step_last), None)
                if jdeb is None or jfin is None:
                    continue
                numcomp = 1 if s.label == "primary" else 2
                if s.kind == "plane":
                    fit = fit_plane(ech, jdeb, jfin, numcomp=numcomp)
                else:
                    fit = fit_line(ech, jdeb, jfin, anchored=s.anchored, numcomp=numcomp)
                if fit is not None:
                    candidates[s.label] = fit

            self._afficher(format_suggestions(ech.id, suggestions))
            if not candidates:
                self._afficher(f"{ech.id}: suggested interval(s) no longer fit (data changed?).\n")
                continue

            self.results.extend(candidates.values())
            # Zijderveld SEUL (panneau unique) - le mini-stereo INTEGRE
            # (au-dessus, voir draw_zijderveld/show_stereo) trace
            # desormais lui-meme un grand cercle pour toute suggestion de
            # PLAN, en plus des points bruts - demande explicite
            # utilisateur ("in data + interpretation, the original plot
            # is lost... I was asking to plot the plane on the stereo
            # with the data" puis "I do not understand why the stereo is
            # still on the right. The routine with the zijder and stereo
            # above the zijder already exist") : le panneau separe
            # "Stereo Results" a cote (build_zijderveld_stereo_results_
            # figure) n'est plus necessaire pour ca, et forcait une
            # figure large a 2 panneaux (donc une fenetre plus large que
            # partout ailleurs dans l'appli). Contourne _current_graphic/
            # _refresh_current_graphic (meme motif qu'avant, ce mode
            # combine n'est pas un des "kind" qu'il connait) mais utilise
            # EXACTEMENT le meme layout/redimensionnement a panneau
            # unique que la vue "Zijderveld" normale.
            self._current_graphic = None
            self._clear_figure()
            self.fig.set_size_inches(5.5, 8.5, forward=True)
            build_zijderveld_figure(
                ech, orientation=self.orientation.get(), fits=list(candidates.values()), fig=self.fig)
            self._fit_figure_to_data()
            self._redraw_canvas()

            choice = self._console_input(
                "Save which (p=primary, s=secondary, b=both, Enter=discard all): ", "")
            if choice is None:
                for fit in candidates.values():
                    self.results.remove(fit)
                return
            choice = choice.strip().lower()
            keep_labels = set()
            if choice == "b":
                keep_labels = set(candidates.keys())
            elif choice == "p" and "primary" in candidates:
                keep_labels = {"primary"}
            elif choice == "s" and "secondary" in candidates:
                keep_labels = {"secondary"}

            for label, fit in candidates.items():
                if label in keep_labels:
                    self._archive_only(fit)
                else:
                    self.results.remove(fit)
            if keep_labels:
                self._afficher(f"Saved: {', '.join(sorted(keep_labels))}  (total results: {len(self.results)})\n")

    def ouvrir_ajusplans_dialog(self):
        """Equivalent GUI de `ajusplans` (linesplans.f) : ajustement de plan
        (grand cercle) par ACP, meme boucle pause/decision que ajuslig. Pas
        de prompt d'ancrage (le Fortran ancre toujours a l'origine pour un
        plan, aucun choix propose a l'utilisateur)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        samples = [s for s in self.selection if len(s.mesures) >= 2]
        if not samples:
            self._showwarning(
                "Not enough measurements", "No selected sample has at least 2 measurements.")
            return

        for ech in samples:
            self._current_graphic = ("zijderveld", ech.id)
            self._refresh_current_graphic()

            steps_text = "  ".join(
                f"{i + 1}:{_step_token(m)}" for i, m in enumerate(ech.mesures))
            self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")

            while True:  # boucle "refaire"
                jdeb_s = self._console_input(
                    f"[{ech.id}] First step (0 = skip this sample): ", "1")
                if jdeb_s is None:
                    return
                try:
                    jdeb = int(jdeb_s)
                except ValueError:
                    self._showerror("Error", "Must be an integer.")
                    continue
                if jdeb == 0:
                    break

                jfin_s = self._console_input("Last step: ", str(len(ech.mesures)))
                if jfin_s is None:
                    return
                try:
                    jfin = int(jfin_s)
                except ValueError:
                    self._showerror("Error", "Must be an integer.")
                    continue
                if not (1 <= jdeb <= jfin <= len(ech.mesures)):
                    self._showerror("Error", "Invalid step range.")
                    continue

                norm_s = self._console_input("Normalize (Y/n): ", "Y")
                if norm_s is None:
                    return
                normalize = norm_s.strip().lower() != "n"

                component = self._console_input(
                    "Magnetization component (A/B/C...): ", "A")
                if component is None:
                    return

                fit = fit_plane(ech, jdeb, jfin, normalize=normalize)
                if fit is None:
                    self._showwarning("Fit rejected", "MAD > 25° (poorly defined plane).")
                    continue
                fit.component = component.strip().upper() or "A"

                dec, inc = _correct_dec_inc(fit, self.orientation.get())
                self._afficher(
                    f"{ech.id} : pole dec={dec:.1f}  inc={inc:.1f}  mad={fit.mad:.1f}  "
                    f"nb points={fit.nb}\n"
                )

                decision = self._console_input(
                    "Save (Y) / redo (r) / next sample (n): ", "Y")
                if decision is None:
                    return
                decision = decision.strip().lower()
                if decision == "r":
                    continue
                if decision != "n":
                    self._save_result(fit)
                    self._afficher(f"Result saved: {len(self.results)}\n")
                break

    def ouvrir_ajusfisher_dialog(self):
        """Equivalent GUI de `ajusfisher` (linesplans.f) : moyenne de Fisher
        sur une plage d'etapes. Repli automatique sur une direction unique
        (branche 735 du Fortran) pour les echantillons a moins de 3 mesures,
        ou si l'utilisateur choisit un seul step (premier == dernier)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        for ech in self.selection:
            if not ech.mesures:
                continue

            if len(ech.mesures) < 3:
                fit = fit_single_direction(ech)
                dec, inc = _correct_dec_inc(fit, self.orientation.get())
                self._afficher(
                    f"{ech.id}: single direction (fewer than 3 measurements) dec={dec:.1f}  "
                    f"inc={inc:.1f}\n"
                )
                decision = self._console_input(
                    f"[{ech.id}] Save the single direction (Y/n): ", "Y")
                if decision is None:
                    return
                if decision.strip().lower() != "n":
                    self._save_result(fit)
                    self._afficher(f"Result saved: {len(self.results)}\n")
                continue

            steps_text = "  ".join(
                f"{i + 1}:{_step_token(m)}" for i, m in enumerate(ech.mesures))
            self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")

            while True:
                jdeb_s = self._console_input(
                    f"[{ech.id}] First step (0 = skip this sample): ", "1")
                if jdeb_s is None:
                    return
                try:
                    jdeb = int(jdeb_s)
                except ValueError:
                    self._showerror("Error", "Must be an integer.")
                    continue
                if jdeb == 0:
                    break

                jfin_s = self._console_input("Last step: ", str(len(ech.mesures)))
                if jfin_s is None:
                    return
                try:
                    jfin = int(jfin_s)
                except ValueError:
                    self._showerror("Error", "Must be an integer.")
                    continue
                if jfin == 0:
                    break
                if not (1 <= jdeb <= jfin <= len(ech.mesures)):
                    self._showerror("Error", "Invalid step range.")
                    continue

                component = self._console_input(
                    "Magnetization component (A/B/C...): ", "A")
                if component is None:
                    return

                if jfin == jdeb:
                    fit = fit_single_direction(ech, index=jdeb - 1)
                    dec, inc = _correct_dec_inc(fit, self.orientation.get())
                    self._afficher(f"{ech.id}: single direction dec={dec:.1f}  inc={inc:.1f}\n")
                else:
                    fit = fit_fisher_direction(ech, jdeb, jfin)
                    dec, inc = _correct_dec_inc(fit, self.orientation.get())
                    self._afficher(
                        f"{ech.id} : dec={dec:.1f}  inc={inc:.1f}  a95={fit.mad:.1f}  "
                        f"k={fit.tx[0]:.1f}  nb={fit.nb}\n"
                    )
                fit.component = component.strip().upper() or "A"

                decision = self._console_input(
                    "Save (Y) / redo (r) / next sample (n): ", "Y")
                if decision is None:
                    return
                decision = decision.strip().lower()
                if decision == "r":
                    continue
                if decision != "n":
                    self._save_result(fit)
                    self._afficher(f"Result saved: {len(self.results)}\n")
                break

    def ouvrir_ajusligauto_dialog(self):
        """Equivalent GUI de `ajusligauto` : ajustement de droite automatique
        sur la totalite des mesures de chaque echantillon selectionne, meme
        ancrage pour tous, sauvegarde sans confirmation."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        ancr_s = self._console_input(
            "Automatic fit anchored to origin (Y/n): ", "Y")
        if ancr_s is None:
            return
        anchored = ancr_s.strip().lower() != "n"

        fits = fit_lines_auto(self.selection, anchored=anchored)
        self._save_results(fits)
        self._afficher(
            f"Automatic fit: {len(fits)} line(s) fitted and saved "
            f"({'anchored' if anchored else 'not anchored'}).\n"
        )

    def ouvrir_ajusligredo_dialog(self):
        """Equivalent GUI de `ajusligredo` : rejoue des ajustements (droite
        ou plan) a partir d'un fichier texte 'redo' (voir docstring de
        calcul.fit_from_redo_file pour le format des lignes). Comme le
        Fortran (`call seloneech(samnum)` a chaque ligne), chaque
        specimen designe dans le fichier redo est recherche et charge
        directement depuis self.donnees - PAS depuis self.selection, qui
        n'a donc pas besoin d'etre preparee au prealable."""
        if not self.donnees:
            self._showwarning("No data", "Load a .ren file first.")
            return

        path = filedialog.askopenfilename(
            title="Redo file (sample L|P o|n tempmin tempmax numcomp)",
            filetypes=[("Text files", "*.txt *.redo"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                lines = f.readlines()
        except Exception as e:
            self._showerror("Error", f"Could not read the file:\n{e}")
            return

        fits = fit_from_redo_file(self.donnees, lines)
        self._save_results(fits)
        self._afficher(
            f"Redo file « {os.path.basename(path)} »: {len(fits)} fit(s) replayed "
            f"and saved.\n"
        )

    def afficher_mdf(self):
        """Equivalent GUI de `mdf` (calcul.f:519-678)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        self._afficher(list_mdf(self.selection))

    def afficher_mean_intensity(self):
        """Equivalent GUI de `mdi`/`mds` (calcul.f:683-865)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        mi = compute_mean_intensity(self.selection)
        if mi is None:
            self._showerror(
                "Error",
                "Could not compute (mass/volume mix in the selection, "
                "or not enough measurements).",
            )
            return
        ms = compute_mean_susceptibility(self.selection)
        ech0 = self.selection[0]
        self._afficher(format_mean_intensity(mi, ms, lat=ech0.lat, rlong=ech0.rlong))

    def ouvrir_koenigsberger_dialog(self):
        """Equivalent GUI de `Koenigs` (calcul.f:3869-3926)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        valk_s = self._console_input("Reference field for the Koenigsberger ratio (µT): ", "40")
        if valk_s is None:
            return
        try:
            valk = float(valk_s)
        except ValueError:
            self._showerror("Error", "Must be a number.")
            return
        if valk == 0.0:
            self._showerror("Error", "The reference field cannot be zero.")
            return
        rows = compute_koenigsberger(self.selection, valk)
        self._afficher(
            format_koenigsberger(rows) if rows
            else "(no measurement with non-zero susceptibility in the selection)"
        )

    def afficher_mean_inclination(self):
        """Equivalent GUI de `meaninc` (calcul.f:1011-1129, estimateur de
        McFadden & Reid 1982)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        r = compute_mean_inclination(self.selection, orientation=self.orientation.get())
        if r is None:
            self._showerror("Error", "Not enough measurements, or the iteration does not converge.")
            return
        self._afficher(format_mean_inclination(r))

    def afficher_diff_measurements(self):
        """Equivalent GUI de `diffmes` (dataselect.f:1156-1236)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        self._afficher(list_diff_measurements(self.selection, orientation=self.orientation.get()))

    def appliquer_viscosity_test(self):
        """Equivalent GUI de `viscos` (viscos.f) : MUTE les 2 premières
        mesures de chaque échantillon (N+/N-) en (vecteur moyen, vecteur
        différence) - operation irréversible sur la sélection en mémoire,
        confirmation demandée."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        if not messagebox.askyesno(
            "Viscosity test",
            "This replaces the first 2 measurements (N+/N-) of each selected "
            "sample with a mean vector and a difference vector - irreversible "
            "for this session. Continue?",
        ):
            return
        warnings = apply_viscosity_test(self.selection)
        msg = "Viscosity test applied."
        if warnings:
            msg += "\n" + "\n".join(warnings)
        self._afficher(msg + "\n")

    def ouvrir_subtraction_dialog(self):
        """Equivalent GUI de `soustra` (dataselect.f:1239-1283) : soustrait
        le vecteur d'une mesure de toutes les autres du même échantillon,
        puis supprime cette ligne - un seul échantillon a la fois, comme le
        Fortran."""
        if len(self.selection) != 1:
            self._showwarning(
                "Invalid selection", "Select a single sample for the subtraction.")
            return
        ech = self.selection[0]
        steps_text = "  ".join(
            f"{i + 1}:{_step_token(m)}" for i, m in enumerate(ech.mesures))
        self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")
        row_s = self._console_input("Line number to subtract: ", "")
        if row_s is None:
            return
        try:
            row = int(row_s)
        except ValueError:
            self._showerror("Error", "Must be an integer.")
            return
        if not messagebox.askyesno(
            "Subtraction",
            f"Subtract line {row} of {ech.id} from all other measurements, "
            f"and remove it - irreversible for this session. Continue?",
        ):
            return
        err = apply_subtraction(ech, row)
        if err:
            self._showerror("Error", err)
            return
        self._afficher(f"{ech.id}: line {row} subtracted and removed ({len(ech.mesures)} measurements remaining).\n")

    def ouvrir_holderarm_dialog(self):
        """Equivalent GUI de `holderarm` (calcul.f:2337-2414) : enregistre
        les 6 mesures ARM d'un porte-échantillon vide - utile pour une
        future correction dans Anisotropy (non encore porté)."""
        if len(self.selection) != 1:
            self._showwarning(
                "Invalid selection",
                "Select the single sample corresponding to the empty holder.")
            return
        ech = self.selection[0]
        steps_text = "  ".join(
            f"{i + 1}:{_step_token(m)}" for i, m in enumerate(ech.mesures))
        self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")
        idx_s = self._console_input("Line numbers x+ x- y+ y- z+ z-: ", "")
        if idx_s is None:
            return
        try:
            ixp, ixm, iyp, iym, izp, izm = (int(v) for v in idx_s.split())
        except ValueError:
            self._showerror("Error", "Exactly 6 space-separated integers are required.")
            return
        bg = record_arm_holder(ech, ixp, ixm, iyp, iym, izp, izm)
        if bg is None:
            self._showerror("Error", "One of the line numbers is invalid.")
            return
        self._arm_holder_background = bg
        self._afficher(
            "ARM holder recorded:\n"
            + "\n".join(f"{i + 1}: {bg.x[i]:.3e}  {bg.y[i]:.3e}  {bg.z[i]:.3e}" for i in range(6))
            + "\n"
        )

    def ouvrir_anisotropy_dialog(self):
        """Equivalent GUI de `anisoauto` (calcul.f:3951+), tenseur 'A0'
        uniquement (voir calcul.compute_anisotropy_tensor) : detection
        automatique des 6 positions X+/X-/Y+/Y-/Z+/Z- (Z+/Z- pouvant etre
        substituees par des mesures R/V a la meme etape que X, ex. 'RH'/
        'VH'), correction optionnelle par la ligne de base porte-
        echantillon (Holder_ARM...), ecriture dans le fichier .ANI (meme
        emplacement/format que celui deja lu par Inverse_ANI_correction...).
        Les 14 variantes jackknife (A+/A-/A1-A6/B1-B6) ne sont PAS portees
        (etape ulterieure).

        Calcule aussi PmagPy (anisotropy_magic.compute_aarm_pmagpy) sur la
        MEME detection de 6 positions, et attache son test F de Hext
        (sigma/F/F12/F23 + verdict satisfactory/not satisfactory) au
        tenseur 'A0' avant l'ecriture .pmagani (colonnes sigma/ftest/
        ftest12/ftest23, jusque-la "n.d" pour ce chemin natif, ET une
        phrase dans la colonne info) - demande explicite utilisateur
        ("ajouter le calcul de PmagPy et les erreurs dans pmagani, ainsi
        que l'estimation dans le comment satisfactory or not
        satisfactory").

        Affiche TOUTE la sequence de calcul dans la fenetre texte, comme
        le fait la console du Fortran d'origine (pas seulement le resultat
        final) - demande explicite de l'utilisateur : les donnees
        d'anisotropie sont souvent basees sur des TRM PARTIELLES (une part
        de la NRM n'a pas ete remplacee), et le tenseur BRUT avant
        symetrisation + les diagnostics par position sont le seul moyen de
        verifier qu'un echantillon etait bien oriente lors de l'acquisition
        de la TRM."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        if not self.results_path:
            self._showwarning(
                "No file", "Load a .ren file first (the .ANI file is derived from it).")
            return
        ani_path = ani_path_for(self.results_path)

        self.text_area.insert(tk.END, "\n--- Anisotropy (Escape to cancel) ---\n", "prompt")
        auto_s = self._console_input("Automated recognition of the 6 steps? y/N: ", "N")
        if auto_s is None:
            return
        use_auto = auto_s.strip().lower() == "y"

        use_zb = False
        if use_auto:
            use_zb_s = self._console_input(
                "Using ZB instead of R for a > +/- 5% evolution? y/N: ", "N")
            if use_zb_s is None:
                return
            use_zb = use_zb_s.strip().lower() == "y"

        done, skipped, declined = [], [], []
        for ech in self.selection:
            lines = [f"\n--- {ech.id}: {'automated ' if use_auto else ''}calcul of TRM or ARM anisotropy with 6 positions ---"]
            lines.append(
                "Remanent magnetization list, normalized by mass (Am2/kg)"
                if ech.norme == "m" else
                "Remanent magnetization list (A/m)"
            )
            # orientation=1 (echantillon) FORCE, pas self.orientation.get()
            # - meme raison que ouvrir_openfilepint_dialog (calcul.f/
            # anisot-anisoauto n'appelle jamais corfor/corpen sur les
            # mesures : le tenseur et les diagnostics par position sont
            # TOUJOURS calcules en coordonnees echantillon, quel que soit
            # le reglage d'orientation courant de l'appli - demande
            # explicite utilisateur ("force the orient to sample
            # coordinates at the beginning of the anisotropy routine").
            # Sans ce forçage, cette liste de mesures pouvait s'afficher
            # dans une orientation differente (In situ/Tilt cor.) de
            # celle, TOUJOURS echantillon, des diagnostics affiches juste
            # apres - une incoherence au sein du meme dialogue.
            buffer = io.StringIO()
            list_measurements([ech], orientation=1, out=buffer)
            lines.append(buffer.getvalue().rstrip("\n"))
            self._afficher("\n".join(lines) + "\n")
            lines = []

            manual_positions = None
            if not use_auto:
                # equivalent du repli manuel de `anisot` (carreco != "y") :
                # pas de substitution R/V->Z+/Z- dans cette branche du
                # Fortran (zplus/zminus n'y sont jamais lus), seulement les
                # codes litteraux X+/X-/Y+/Y-/Z+/Z-.
                # Return (pas de numeros saisis) = passer CE specimen et
                # continuer la boucle ; Escape = sortir de TOUTE la boucle
                # batch (les specimens deja traites restent acquis,
                # resume affiche quand meme) - demande explicite
                # utilisateur ("can you add return to skip, escape to go
                # out of the loop") : avant ce changement, Escape ne
                # sautait qu'UN specimen (comme Return, qui en plus
                # affichait une erreur "6 valid line numbers are
                # required" avant de sauter) - aucun moyen d'interrompre
                # le lot entier depuis ce prompt.
                idx_s = self._console_input(
                    f"[{ech.id}] Line numbers x+ x- y+ y- z+ z- "
                    f"(Return to skip, Escape to stop): ", "")
                if idx_s is None:
                    break
                if not idx_s.strip():
                    skipped.append(ech.id)
                    continue
                try:
                    idx = [int(v) for v in idx_s.split()]
                    if len(idx) != 6 or any(not (1 <= i <= len(ech.mesures)) for i in idx):
                        raise ValueError
                except ValueError:
                    self._showerror("Error", "6 valid line numbers are required.")
                    skipped.append(ech.id)
                    continue
                manual_positions = dict(zip(
                    ("X+", "X-", "Y+", "Y-", "Z+", "Z-"),
                    (ech.mesures[i - 1] for i in idx),
                ))

            result = compute_anisotropy_tensor(
                ech, holder=self._arm_holder_background, use_zb_on_evolution=use_zb,
                positions=manual_positions,
            )
            if result is None:
                skipped.append(ech.id)
                lines.append("Decoding incomplete: could not identify all 6 positions (X+/X-/Y+/Y-/Z+/Z-).")
                self._afficher("\n".join(lines) + "\n")
                continue

            unit = "Am2/kg" if ech.norme == "m" else "A/m"

            def _position_diag_lines(res):
                out = []
                if res.holder_used:
                    bg = self._arm_holder_background
                    out.append("Holder ARM background subtracted (X+,X-,Y+,Y-,Z+,Z-):")
                    out.extend(
                        f"  {i + 1}: {bg.x[i]:.3E}  {bg.y[i]:.3E}  {bg.z[i]:.3E}" for i in range(6)
                    )
                if res.swapped_axes:
                    out.extend(
                        f"!! inversion {ax}+ and {ax}- for sample: {ech.id}" for ax in res.swapped_axes
                    )
                if res.misoriented_positions:
                    # Deviation angulaire par rapport a l'axe theorique de
                    # la position (voir calcul._find_misoriented_
                    # positions) - PAS deja corrige (contrairement a
                    # swapped_axes ci-dessus, un simple echange de signe
                    # ne "repare" pas une deviation > 45 deg) - demande
                    # explicite utilisateur ("mettre un warning lorsque le
                    # code ne correspond pas à la vraie direction
                    # enregistrée par l'échantillon... jamais à + de 45°").
                    out.extend(
                        f"!! WARNING {key} direction is {dev:.0f}° from the expected axis "
                        f"for sample {ech.id} - likely mis-oriented in the oven (anisotropy alone "
                        f"should never deviate a TRM by more than 45°)"
                        for key, dev in res.misoriented_positions
                    )
                out.append("NRM residual estimated from each pair (before averaging the 3):")
                out.extend(
                    f"  {p.key}: {p.intensity:10.3E} {unit}  dec={p.dec:6.1f}  inc={p.inc:6.1f}"
                    for p in res.pair_nrm
                )
                m = res.nrm_mean_diag
                out.append(
                    f"Mean NRM (residual, subtracted below): "
                    f"{m.intensity:10.3E} {unit}  dec={m.dec:6.1f}  inc={m.inc:6.1f}"
                )
                out.append("Mean NRM subtracted - TRM values:")
                out.extend(
                    f"  {d.key}: {_step_token(d.measurement)}  "
                    f"{d.intensity:10.3E} {unit}  dec={d.dec:6.1f}  inc={d.inc:6.1f}"
                    for d in res.position_diags
                )
                out.append(f"Deviation (pair asymmetry / mean TRM intensity): {res.deviation_pct:5.1f} %")
                if res.trm_evolution_pct is not None:
                    out.append(f"ZB check found - TRM evolution: {res.trm_evolution_pct:5.1f} %")
                    if res.zb_used:
                        out.append("  -> evolution exceeds 5%: Z+ replaced by the ZB control measurement.")
                else:
                    out.append("No ZB control measurement found (no pTRM-check available).")
                return out

            self._afficher("\n".join(lines + _position_diag_lines(result)) + "\n")
            lines = []

            # Equivalent GUI de la question DESACTIVEE dans le Fortran
            # d'origine ("voulez-vous supprimer une des 6 mesures ?
            # laquelle:", calcul.f:4481-4552 - commentee "c", jamais
            # active dans la version livree, MAIS le bloc de
            # reconstruction qui la suit est bien present, complet, en
            # dur) - demande explicite utilisateur : "in the original
            # Fortran code... after [listing the NRM/TRM values] I was
            # asking if I wanted to replace a step" puis, precision
            # ulterieure avec un exemple reel du Fortran a l'appui,
            # "when we remove one line, it is replaced by the
            # antiparallel... also easier to write the number of the
            # line from 1 to 6" : PAS une nouvelle mesure a choisir
            # (contrairement a un essai precedent de ce port) - le
            # numero 1-6 correspond exactement a l'ordre deja affiche
            # dans "Mean NRM subtracted - TRM values" ci-dessus (memes
            # positions que _ANI_POSITION_KEYS). Voir
            # calcul.replace_position_by_symmetry pour la formule exacte
            # portee depuis ce bloc Fortran. Le tenseur est RECALCULE en
            # entier (positions/NRM moyenne/inversion X-Y/ZB re-evaluees)
            # apres chaque remplacement, pas juste la position changee,
            # puisque tout en depend.
            position_order = ("X+", "X-", "Y+", "Y-", "Z+", "Z-")
            while True:
                replace_s = self._console_input(
                    f"[{ech.id}] Do you want to remove one of the 6 measurements?\n"
                    "Which one (1-6, N to skip): ", "N")
                if replace_s is None:
                    # Escape == "N" ici (le prompt le dit deja : "N to
                    # skip") - PAS un `return` (bug signale par
                    # l'utilisateur, "now it skip the process here" :
                    # `return` sortait de TOUTE la methode en silence,
                    # sans le resume final done/skipped/declined ni les
                    # specimens suivants du lot batch - une simple
                    # decision "ne pas retirer de mesure" ne doit annuler
                    # ni le specimen courant ni le reste du lot).
                    break
                replace_s = replace_s.strip()
                if replace_s.lower() in ("", "n"):
                    break
                try:
                    line_no = int(replace_s)
                    if not (1 <= line_no <= 6):
                        raise ValueError
                except ValueError:
                    self._showerror("Error", "Must be a number from 1 to 6 (or N).")
                    continue
                pos_key = position_order[line_no - 1]
                new_positions = replace_position_by_symmetry(
                    result.positions, self._arm_holder_background, pos_key)
                result = compute_anisotropy_tensor(
                    ech, holder=self._arm_holder_background, use_zb_on_evolution=use_zb,
                    positions=new_positions,
                )
                self._afficher(
                    f"{pos_key} reconstructed from its pair partner and the other two "
                    "pairs' mean NRM (not a new measurement line).\n"
                    + "\n".join(_position_diag_lines(result)) + "\n"
                )

            lines.append("-------- Tensor (raw, before symmetrization) --------")
            lines.extend(f"  {r[0]:10.3E}  {r[1]:10.3E}  {r[2]:10.3E}" for r in result.raw)
            lines.append("(forcing symmetric)")

            t = result.tensor
            lines.append("-------- Tensor 'A0' (symmetric) --------")
            lines.append(f"  {t.k11:10.3E}  {t.k12:10.3E}  {t.k13:10.3E}")
            lines.append(f"  {t.k12:10.3E}  {t.k22:10.3E}  {t.k23:10.3E}")
            lines.append(f"  {t.k13:10.3E}  {t.k23:10.3E}  {t.k33:10.3E}")

            # PmagPy (Hext 1963 F-test, voir anisotropy_magic.
            # compute_aarm_pmagpy) sur la MEME detection de 6 positions -
            # attache sigma/F/F12/F23 + verdict satisfactory/not
            # satisfactory au tenseur 'A0' (all_tensors[0]) avant
            # l'ecriture .pmagani - demande explicite utilisateur
            # ("ajouter le calcul de PmagPy et les erreurs dans pmagani,
            # ainsi que l'estimation dans le comment satisfactory or not
            # satisfactory"). `nrm_mean=result.nrm_mean` : BUG REEL corrige
            # (voir docstring compute_aarm_pmagpy) - sans elle, la moindre
            # NRM residuelle (TRM/ARM partielle) contamine le moindres-
            # carres pmagpy et fait chuter le F-test vers "not significant"
            # meme pour une vraie anisotropie moderee a forte.
            pmagpy_result = compute_aarm_pmagpy(
                result.positions, holder=self._arm_holder_background, n_pos=6,
                nrm_mean=result.nrm_mean)
            a0_tensor = result.all_tensors[0]
            a0_tensor.n_positions = 6
            a0_tensor.sigma = pmagpy_result.sigma
            a0_tensor.ftest = pmagpy_result.f_test
            a0_tensor.ftest12 = pmagpy_result.f12_test
            a0_tensor.ftest23 = pmagpy_result.f23_test
            a0_tensor.f_crit = pmagpy_result.f_crit
            a0_tensor.quality = pmagpy_result.quality
            # Detail complet (sigma/F/F12/F23/axes propres/s trace-
            # normalise), pas juste F/critical F - demande explicite
            # utilisateur ("can you list the PmagPy details") : MEME
            # fonction que le menu separe "Anisotropy PmagPy" (voir
            # format_aarm_pmagpy), pour ne pas avoir deux formats de
            # resume differents pour le meme calcul.
            lines.append(format_aarm_pmagpy(ech.id, pmagpy_result))

            self._afficher("\n".join(lines) + "\n")
            lines = []

            # Pause TOUJOURS ici, quelle que soit `quality` - demande
            # explicite utilisateur ("is it possible to still pause at the
            # save to pmagani waiting for the validation by a return. The
            # user should pay attention at the output even if the
            # anisotropy is satisfactory") : avant ce changement, un
            # tenseur 'g' (satisfactory) s'ecrivait sans pause des que le
            # F-test passait, laissant filer le detail PmagPy (sigma/F/F12/
            # F23/axes) sans que l'utilisateur ait eu l'occasion de le
            # regarder - desormais CHAQUE specimen s'arrete sur Entree,
            # 'g' comme 'b', seul le texte du prompt change. Defaut Y
            # (Entree = sauvegarder) dans les deux cas - demande explicite
            # anterieure ("reverse the y/N to Y/n by default") : la reponse
            # est presque toujours "garder", la pause sert a FAIRE LIRE le
            # resultat, pas a decourager de sauvegarder.
            if pmagpy_result.quality == "b":
                prompt = (
                    f"[{ech.id}] anisotropy NOT significant (not satisfactory).\n"
                    f"Save to .pmagani anyway? Y/n: "
                )
            elif pmagpy_result.quality == "g":
                prompt = f"[{ech.id}] anisotropy significant (satisfactory). Save to .pmagani? Y/n: "
            else:
                prompt = f"[{ech.id}] Hext F-test not computable. Save to .pmagani? Y/n: "
            confirm = self._console_input(prompt, "Y")
            if confirm is None or confirm.strip().lower() == "n":
                declined.append(ech.id)
                self._afficher(f"{ech.id}: not saved (declined).\n")
                continue

            write_ani_tensors(
                ani_path, ech, result.all_tensors, positions=result.positions,
                trm_evolution_pct=result.trm_evolution_pct or 0.0,
                deviation_pct=result.deviation_pct,
            )
            done.append(ech.id)

        summary = (
            f"15 tensor variants (A0/A+/A-/A1/B1/A2/B2/A3/B3/A4/B4/A5/B5/A6/B6, "
            f"'A0' being the main one) written for {len(done)} sample(s) -> {ani_path}\n"
        )
        if skipped:
            summary += f"Skipped (6 positions not identified): {', '.join(skipped)}\n"
        if declined:
            summary += f"Not saved (anisotropy not satisfactory, declined): {', '.join(declined)}\n"
        self._afficher(summary)

    def ouvrir_anisotropy_pmagpy_dialog(self):
        """Menu "Anisotropy PmagPy..." - PAS dans le Fortran, demande
        explicite utilisateur ("I have never been able to figure out how
        the anisotropy of remanent magnetization is determined in Magic.
        Can you figure out from the PmagPy tools and add a menu
        Anisotropy PmagPy"). 2e pipeline PARALLELE (comme le panneau
        paleointensite PmagPy) appelant le vrai pmagpy (voir
        anisotropy_magic.py pour le detail de ce qui a ete trouve dans le
        source pmagpy - `ipmag.get_matrix`/`ipmag.calculate_aniso_
        parameters`, moindres carres a la Hext 1963 + test de
        significativite F, absent du calcul natif STARpaleomag_Py).

        Detection des 6 positions IDENTIQUE a "Anisotropy" (detect_six_
        positions, meme substitution R/V->Z+/Z- si necessaire), avec la
        MEME ligne de base porte-echantillon optionnelle - affiche le
        tenseur trace-normalise de pmagpy (`aniso_s`), ses axes
        principaux (v1/v2/v3), le degre d'anisotropie P et le test F de
        Hext, cote a cote avec le tenseur natif 'A0' de STARpaleomag_Py pour
        montrer qu'ils calculent la MEME chose (verifie sur donnees
        synthetiques - voir anisotropy_magic.py) - pmagpy apportant en
        plus le test de significativite statistique."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return

        self.text_area.insert(tk.END, "\n--- Anisotropy PmagPy (Escape to cancel) ---\n", "prompt")
        done, skipped = [], []
        for ech in self.selection:
            positions = detect_six_positions(ech)
            if positions is None:
                skipped.append(ech.id)
                self._afficher(
                    f"{ech.id}: decoding incomplete - could not identify all 6 "
                    f"positions (X+/X-/Y+/Y-/Z+/Z-).\n"
                )
                continue

            starmac = compute_anisotropy_tensor(
                ech, holder=self._arm_holder_background, positions=positions)
            # nrm_mean=starmac.nrm_mean : meme correction que dans
            # ouvrir_anisotropy_dialog (voir docstring compute_aarm_pmagpy)
            # - sans elle, cette comparaison cote-a-cote afficherait un F-
            # test pmagpy contamine par la NRM residuelle alors que le
            # tenseur natif 'A0' juste au-dessus en est deja immunise.
            pmagpy_result = compute_aarm_pmagpy(
                positions, holder=self._arm_holder_background, n_pos=6,
                nrm_mean=starmac.nrm_mean if starmac is not None else None)

            lines = [f"\n--- {ech.id} ---"]
            if starmac is not None:
                t = starmac.tensor
                trace = t.k11 + t.k22 + t.k33
                lines.append(
                    f"STARpaleomag_Py 'A0' (raw, physical units): k11={t.k11:.3E}  k22={t.k22:.3E}  "
                    f"k33={t.k33:.3E}  k12={t.k12:.3E}  k23={t.k23:.3E}  k13={t.k13:.3E}"
                )
                if trace:
                    lines.append(
                        f"STARpaleomag_Py 'A0' (trace-normalized, for comparison with pmagpy's s "
                        f"below): s1={t.k11/trace:.5f}  s2={t.k22/trace:.5f}  "
                        f"s3={t.k33/trace:.5f}  s4(s12)={t.k12/trace:.5f}  "
                        f"s5(s23)={t.k23/trace:.5f}  s6(s13)={t.k13/trace:.5f}"
                    )
            lines.append(format_aarm_pmagpy(ech.id, pmagpy_result))
            done.append(ech.id)
            self._afficher("\n".join(lines) + "\n")

        summary = f"Anisotropy PmagPy: {len(done)} sample(s) processed.\n"
        if skipped:
            summary += f"Skipped (6 positions not identified): {', '.join(skipped)}\n"
        self._afficher(summary)

    def ouvrir_inverseani_dialog(self):
        """Equivalent GUI de `inverseani` (plotpaleoint2.f:1803-1918) :
        applique l'inverse d'un tenseur d'anisotropie (lu dans un fichier
        .ANI, dérivé du fichier de données comme le fichier .r) à toutes
        les mesures de l'échantillon sélectionné."""
        if len(self.selection) != 1:
            self._showwarning(
                "Invalid selection", "Select a single sample.")
            return
        ech = self.selection[0]
        if getattr(ech, "flaganiso", False):
            self._showwarning("Already corrected", f"{ech.id} has already been corrected for anisotropy.")
            return
        if not self.results_path:
            self._showwarning(
                "No file", "Load a .ren file first (the .ANI file is derived from it).")
            return
        ani_path = self._find_ani_path() or ani_path_for(self.results_path)
        choice = self._console_input(
            "Correction with 1: TRM tensor (A0)  2: ARM tensor (F0)  3: susceptibility (N0)  "
            "4: ARM tensor imported from MagIC (AA): ", "1")
        if choice is None:
            return
        try:
            ichoice = int(choice)
        except ValueError:
            ichoice = 1
        if ichoice not in _ANI_CODE2:
            ichoice = 1
        code2 = _ANI_CODE2[ichoice]

        tensor = read_ani_tensor(ani_path, ech.id, code2)
        if tensor is None:
            self._showerror(
                "Error", f"No '{code2}' tensor for {ech.id} in {ani_path}.")
            return
        apply_inverse_anisotropy(ech, tensor)
        ech.flaganiso = True
        self._afficher(f"{ech.id}: measurements corrected (inverse of the {code2} tensor).\n")

        # Equivalent de `call lismes` en fin de `inverseani` (plotpaleoint2.f:
        # 1897) : liste les mesures corrigees (dec/inc recalcules depuis les
        # x,y,z desormais corriges) - demande explicite utilisateur ("can you
        # list also the corrected directions after Inverse anisotropy"),
        # meme format que "List data" (selection.list_measurements).
        buffer = io.StringIO()
        list_measurements([ech], orientation=self.orientation.get(), out=buffer)
        self._afficher("\nCorrected directions:\n" + buffer.getvalue())

    def ouvrir_cooling_rate_dialog(self):
        """Equivalent GUI de `vitref` (calcul.f:3461-3625, branche "live") :
        boucle sur la sélection, détection automatique du motif L/Q + R/V
        (sinon saisie manuelle des 5 numéros de ligne)."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        auto_s = self._console_input("Automated procedure (Y/n): ", "Y")
        if auto_s is None:
            return
        auto = auto_s.strip().lower() != "n"

        for ech in self.selection:
            rows = detect_cooling_rate_rows(ech) if auto else None
            if rows is None:
                steps_text = "  ".join(
                    f"{i + 1}:{_step_token(m)}" for i, m in enumerate(ech.mesures))
                self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")
                if auto:
                    self._afficher(f"{ech.id}: L/Q + R/V pattern not found automatically.\n")
                rows_s = self._console_input(
                    f"[{ech.id}] Line numbers ATR+fast ATR-fast ATR+slow before/after loop "
                    "(empty = skip): ", "")
                if rows_s is None:
                    return
                if not rows_s.strip():
                    continue
                try:
                    rows = tuple(int(v) for v in rows_s.split())
                    if len(rows) != 5:
                        raise ValueError
                except ValueError:
                    self._showerror("Error", "Exactly 5 integers are required.")
                    continue

            result = compute_cooling_rate(ech, *rows)
            self._afficher(format_cooling_rate(result) + "\n")

    @staticmethod
    def _format_fisher(stats, title):
        flip_note = (
            f" ({stats.n_flipped} direction(s) inverted - antipodal group combine)\n"
            if stats.n_flipped else "\n"
        )
        return (
            f"{title}\n"
            "----------------------------------------\n"
            f"n = {stats.n}{flip_note}"
            f"Dec = {stats.dec:.1f}   Inc = {stats.inc:.1f}\n"
            f"R = {stats.r:.3f}   k = {stats.k:.1f}\n"
            f"a95 = {stats.a95:.1f}   csd = {stats.csd:.1f}\n"
        )

    def _ask_combine_antipodal(self):
        """Prompt commun a fisher_mesures/fisher_resultats - voir
        calcul.combine_antipodal_groups. Defaut 'n' : ne change pas le
        comportement existant tant que l'utilisateur ne le demande pas
        explicitement (typiquement une courte section magnetostratigraphique
        avec un intervalle de polarite minoritaire)."""
        ans = self._console_input(
            "Combine antipodal N/R groups before averaging "
            "(invert the smaller group)? (y/N): ", "n")
        if ans is None:
            return None
        return (ans.strip().lower() or "n").startswith("y")

    def fisher_mesures(self):
        """Equivalent GUI de `fishmes`."""
        if not self.selection:
            self._showwarning("No selection", "Select some samples first.")
            return
        combine_antipodal = self._ask_combine_antipodal()
        if combine_antipodal is None:
            return
        try:
            stats = fisher_from_measurements(
                self.selection, orientation=self.orientation.get(),
                combine_antipodal=combine_antipodal)
        except ValueError as e:
            self._showerror("Error", str(e))
            return
        self._afficher(self._format_fisher(stats, "Fisher on the selected measurements"))

    def fisher_resultats(self):
        """Equivalent GUI de `fishres` (limité aux résultats de type ligne).
        Propose ensuite d'archiver la moyenne calculee comme un resultat
        "mean:" dans le .pmagres - demande explicite utilisateur ("in
        fisher results, a mean is calculated for the site, if the user
        want to archive this mean, this mean needs to be recorded in the
        .pmagres file") : jusqu'ici cette fonction se contentait
        d'afficher les statistiques, sans jamais les archiver (seul
        l'import MagIC produisait des lignes "mean:" - voir
        calcul.build_site_mean_result)."""
        if not self.results:
            self._showwarning("No results", "Run one or more line fits first.")
            return
        combine_antipodal = self._ask_combine_antipodal()
        if combine_antipodal is None:
            return
        try:
            stats = fisher_from_results(
                self.results, orientation=self.orientation.get(),
                combine_antipodal=combine_antipodal)
        except ValueError as e:
            self._showerror("Error", str(e))
            return
        self._afficher(self._format_fisher(stats, "Fisher on the results (lines)"))
        self._archive_fisher_mean(stats, combine_antipodal)

    def _archive_fisher_mean(self, stats, combine_antipodal=False):
        """Propose d'archiver `stats` (calcule par fisher_resultats) comme
        un resultat "mean:" - demande le site, la composante de
        magnetisation (A/B/C..., voir calcul.FitResult.component -
        JAMAIS deduite de numcomp des resultats individuels, "I think it
        is best not to use the numcomp of individual samples for the
        mean"), et cherche lat/lon du site dans self.donnees pour le
        calcul du VGP (dir_to_vgp).

        Si la moyenne archivee n'est pas deja en coordonnees apres pendage
        (TC, orientation 3), archive AUSSI automatiquement la moyenne de
        Fisher recalculee apres correction de pendage, a partir des memes
        resultats individuels - demande explicite utilisateur ("for site
        when the bedding_dip is not equal to 0, it is also best to
        automatically calculate the site-mean result after bedding
        correction. best to recalculate the fisher from individual
        results after tilt correction; the user will not need to do
        it"), etendue ensuite a TOUS les sites meme a pendage nul ("in
        fact better to do it in all cases even when the dip is zero. It
        will be easier to have all sites in the same coordinates") : plus
        de condition sur le pendage, seule l'orientation de la moyenne
        deja archivee (orientation==3, deja TC) evite un second
        archivage redondant."""
        if not self.results_path:
            self._afficher(
                "(no .r file for the loaded data - the mean cannot be "
                "archived; run a fit first, or load a .ren file)\n")
            return
        decision = self._console_input(
            "Archive this mean to the .pmagres file (y/N): ", "n")
        if decision is None or decision.strip().lower() != "y":
            return

        # meme regroupement que fisher_from_results (L/f/s + P, PAS
        # uniquement 'L') - demande explicite utilisateur ("when lines and
        # planes are selected for a fisher result, use the combined L & P").
        contributing = [r for r in self.results if r.cat1 in _LINE_LIKE_CAT1 or r.cat1 == "P"]
        # site_of_result (magic_site du specimen correspondant, PAS les 6
        # premiers caracteres de l'id) - demande explicite utilisateur,
        # apres renommage de site ("in the calculation of a mean for a
        # site, it does not use the information from the site") : un site
        # renomme via "Complete sample information" (mode site, colonne
        # 'site' repetee) garde son id de specimen INCHANGE (voir
        # complete_sample_info.py) - grouper par id[:6] retombait donc sur
        # l'ancien nom de campagne, jamais sur le nom renomme.
        sites = sorted({site_of_result(r, self.donnees) for r in contributing})
        if len(sites) == 1:
            # un seul site dans les resultats combines -> pas d'ambiguite,
            # nom attribue automatiquement sans demander confirmation -
            # demande explicite utilisateur ("by default when there are
            # results from a single site, it is best to automatically [use]
            # the name of the site to the mean value").
            site = sites[0]
            self._afficher(
                f"Site name automatically set to « {site} » (single site in the results).\n")
        else:
            site = self._console_input("Site name: ", "")
            if site is None:
                return
            site = site.strip()
            if not site:
                self._showerror("Error", "A site name is required.")
                return

        latlon = site_lat_lon_from_donnees(site, self.donnees)
        if latlon is None:
            self._afficher(
                f"(no specimen of site « {site} » found in the loaded data - "
                "site lat/lon set to 0.0/0.0, the archived VGP will not be "
                "meaningful)\n"
            )
            site_lat, site_lon = 0.0, 0.0
        else:
            site_lat, site_lon = latlon

        component = self._console_input(
            "Magnetization component (A/B/C...): ", "A")
        if component is None:
            return
        component = component or "A"

        orientation = self.orientation.get()
        mean_fit = build_site_mean_result(
            stats, contributing, site, orientation,
            site_lat=site_lat, site_lon=site_lon, component=component,
        )
        self._save_result(mean_fit)
        self._afficher(
            f"Mean archived: mean: {site}  [component {mean_fit.component}]  "
            f"n={mean_fit.nb}  dec={mean_fit.dec:.1f}  inc={mean_fit.inc:.1f}  "
            f"VGP=({mean_fit.par4:.1f},{mean_fit.par5:.1f})\n"
        )

        if orientation == 3:
            return
        try:
            tc_stats = fisher_from_results(
                contributing, orientation=3, combine_antipodal=combine_antipodal)
        except ValueError:
            return
        tc_mean = build_site_mean_result(
            tc_stats, contributing, site, 3,
            site_lat=site_lat, site_lon=site_lon, component=component,
        )
        self._save_result(tc_mean)
        self._afficher(
            f"Tilt-corrected mean also automatically archived: mean: {site}  "
            f"[component {tc_mean.component}]  (TC)  n={tc_mean.nb}  "
            f"dec={tc_mean.dec:.1f}  inc={tc_mean.inc:.1f}  "
            f"VGP=({tc_mean.par4:.1f},{tc_mean.par5:.1f})\n"
        )

    def lister_resultats(self):
        """Equivalent GUI de `lisres` (limité aux résultats de type ligne).
        Passe l'orientation courante : les dec/inc affichés sont recalculés
        dans ce repère à chaque appel (voir list_results), pas figés dans
        le repère échantillon d'origine."""
        if not self.results:
            self._showwarning("No results", "Run one or more line fits first.")
            return
        self._afficher(list_results(self.results, orientation=self.orientation.get(),
                                     donnees=self.donnees))

    def evaluer_interpretations(self):
        """Pas dans le Fortran (aucun equivalent) - demande explicite
        utilisateur ("une routine qui permettrait d'evaluer les
        interpretations"). Audite self.results (memes resultats deja
        charges/produits que List results) : recalcule chaque fit ligne/
        plan depuis les mesures vivantes (self.donnees) et note la
        qualite (voir interpretation_quality.py pour le detail des
        criteres - MAD recalcule, angle ancre/libre pour les droites,
        ratio de linearite). Les moyennes de site ("mean:") sont ignorees
        (pas un fit ligne/plan)."""
        if not self.results:
            self._showwarning("No results", "Run one or more line/plane fits first, or load some via Select results...")
            return
        reports = evaluate_results(self.results, self.donnees)
        self._afficher(format_quality_report(reports) + "\n")

    def reinitialiser_resultats(self):
        """Equivalent GUI de `initres`."""
        self.results = init_results()
        self._afficher("Results reset - no result saved.\n")

    def ouvrir_delete_results_dialog(self):
        """Menu "Delete results..." - PAS dans le Fortran (`initres` n'efface
        que la totalite) - demande explicite utilisateur ("ajouter un menu
        delete results from the list in memory"). Retire des resultats de
        self.results (MEMOIRE uniquement - le fichier .pmagres deja archive
        n'est jamais modifie/reecrit ici), selectionnes par leur numero
        affiche par "List results" (list_results), pour cibler precisement
        un ou plusieurs resultats plutot que de tout reinitialiser."""
        if not self.results:
            self._showwarning("No results", "No result currently loaded.")
            return
        self._afficher(list_results(self.results, orientation=self.orientation.get(),
                                     donnees=self.donnees))
        indices_s = self._console_input(
            "Indices to delete (e.g. 1,3,5-7 - see the numbers above; "
            "Escape to cancel): ", "")
        if indices_s is None:
            return
        indices_s = indices_s.strip()
        if not indices_s:
            return

        to_delete = set()
        try:
            for token in indices_s.split(","):
                token = token.strip()
                if not token:
                    continue
                if "-" in token:
                    a, b = token.split("-", 1)
                    to_delete.update(range(int(a), int(b) + 1))
                else:
                    to_delete.add(int(token))
        except ValueError:
            self._showerror("Error", "Invalid index list (use e.g. 1,3,5-7).")
            return

        n = len(self.results)
        invalid = sorted(i for i in to_delete if i < 1 or i > n)
        if invalid:
            self._showerror("Error", f"Index out of range: {invalid} (1-{n}).")
            return

        kept = [r for i, r in enumerate(self.results, start=1) if i not in to_delete]
        n_removed = n - len(kept)
        self.results = kept
        self._afficher(
            f"{n_removed} result(s) removed from memory (the .pmagres file, if "
            f"already archived, is unaffected) - {len(self.results)} remaining.\n"
        )

    def ouvrir_import_published_means_dialog(self):
        """Archive des moyennes de site PUBLIEES (table externe) comme
        resultats "mean:", SANS aucun specimen/mesure requis - demande
        explicite utilisateur ("how do we fill the prmag in that case",
        suite a "can we have a site in prmag without data?" pour un
        fichier legacy dont les mesures brutes sont perdues) : le SEUL
        moyen d'obtenir un resultat "mean:" exploitable jusqu'ici etait de
        calculer une moyenne depuis des specimens reellement charges
        ("best dir Fisher..." + archivage), un import MagIC deja publie,
        ou un script ad hoc - voir calcul.import_published_means pour le
        format de table exact et le detail complet."""
        table_path = filedialog.askopenfilename(
            title="Select the published site means table "
                  "(site, orientation [IS/TC], dec, inc, a95, k, n, ...)",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")],
        )
        if not table_path:
            return

        pmagres_path = self.results_path
        if not pmagres_path:
            self._showwarning(
                "No results file",
                "No .r/.pmagres file for the loaded data (load a .ren/.prmag file first, "
                "so there is somewhere to archive these means to).",
            )
            return

        try:
            n_imported, errors = import_published_means(pmagres_path, table_path)
        except Exception as e:
            self._showerror("Error", f"Import failed:\n{e}")
            return

        self._afficher(
            f"{n_imported} published site mean(s) archived to {pmagres_path}.\n"
        )
        if errors:
            self._afficher(f"{len(errors)} row(s) skipped:\n" + "\n".join(errors) + "\n")
        # self.results n'est PAS mis a jour automatiquement (les moyennes
        # sont ecrites sur disque, pas chargees en memoire) - meme
        # convention que l'archivage normal (_save_result les ajoute a
        # self.results EN PLUS d'ecrire sur disque, mais ici il n'y a pas
        # d'equivalent "resultat courant" a ajouter un par un) : utiliser
        # "Select results..." (mode m/s) pour les recharger si besoin de
        # les voir/tracer immediatement dans cette session.
        self._afficher(
            "(use \"Select results...\" mode 'm' or 's' to load them into this session)\n"
        )

    def ouvrir_selres_dialog(self):
        """Equivalent GUI de `selres` (dataselect.f) : charge des resultats
        DEPUIS le fichier .r (equivalent filr) dans self.results. Trois
        modes, comme le Fortran (`carselect` : "Data (default), mean (m) ou
        site [data+mean] (s)") :
        - Data : resultats normaux (L/P/f/s), filtres par echantillon/type/
          composante - les moyennes "mean:" sont exclues.
        - Mean : uniquement les moyennes de site, filtrees par orientation
          courante (equivalent `res.par3==float(iorient)`) - pas de filtre
          type/composante numerique (non demandes par le Fortran dans ce
          mode).
        - Site : la/les moyenne(s) matchee(s) PLUS les resultats individuels
          qui la composent (equivalent `decodelisteres`).

        Filtre supplementaire "component" (A/B/C..., PAS dans le Fortran -
        voir calcul.FitResult.component) disponible dans les 3 modes :
        demande explicite utilisateur pour pouvoir selectionner
        independamment UNE seule moyenne de site (et ses resultats
        associes en mode Site) quand un meme site porte plusieurs moyennes
        de composantes de magnetisation differentes ("when there is
        different components of magnetizations within the same site, we
        may have two or three means by site... during the select results
        with m or s, we can also differentiate by component")."""
        if not self.results_path or not os.path.exists(self.results_path):
            self._showwarning(
                "No results file",
                "No .r file found for the loaded data "
                "(run a fit first, or load a .ren file)."
                if not self.results_path else
                f"{self.results_path} does not exist yet (no result archived).",
            )
            return

        # Si results_path a change sur le disque depuis la derniere lecture
        # (mtime) - typiquement une edition manuelle/externe du .pmagres,
        # ex. un dedoublonnage ou une correction faite hors de l'appli -
        # self.results et self._archived_ids sont vides/recalcules avant
        # de charger, plutot que d'accumuler par-dessus un etat perime
        # (des resultats deja retires du fichier resteraient sinon
        # visibles en memoire) - demande explicite utilisateur ("is it
        # possible to reload the file before any selection of results:
        # this file is usually edited by the user"). Comportement
        # d'accumulation NORMAL (selres, "the selection of results should
        # not initialize the previous list") inchange tant que le fichier
        # n'a pas change entre deux "Select results..." consecutifs.
        current_mtime = os.path.getmtime(self.results_path)
        if self._results_path_mtime is not None and current_mtime != self._results_path_mtime:
            if self.results:
                self._afficher(
                    f"({self.results_path} was modified on disk since it was last read - "
                    f"reloading from disk, {len(self.results)} previously loaded result(s) "
                    "discarded)\n"
                )
            self.results = []
            self._archived_ids = None
        self._results_path_mtime = current_mtime

        self.text_area.insert(tk.END, "\n--- Select results (Escape to cancel) ---\n", "prompt")
        carselect = self._console_input(
            "Data (default), mean (m) or site [data+mean] (s): ", "d")
        if carselect is None:
            return
        carselect = (carselect.strip().lower() or "d")[:1]

        if carselect in ("m", "s"):
            site = self._console_input("Site (name, without « mean: », * = all): ", "*")
            if site is None:
                return
            component = self._console_input(
                "Magnetization component (A/B/C..., blank = all): ", "")
            if component is None:
                return
            loaded = load_results(
                self.results_path, pattern=site or "*", carselect=carselect,
                iorient=self.orientation.get(), component=component or "*",
            )
            if not loaded:
                available = available_mean_orientations(self.results_path, site or "*")
                if available and self.orientation.get() not in available:
                    noms = {1: "Sample (SC)", 2: "In situ (IS)", 3: "Tilt cor. (TC)"}
                    self._afficher(
                        "No mean found for the current orientation "
                        f"« {noms.get(self.orientation.get(), self.orientation.get())} » - "
                        "but some exist for: "
                        + ", ".join(noms.get(o, str(o)) for o in available)
                        + ". Change the orientation (Pmag data menu) then try again.\n"
                    )
        else:
            pattern = self._read_prefixed_pattern("Sample (* = all): ")
            if pattern is None:
                return
            cat1 = self._console_input("Type (L/P/f/s, * = all): ", "*")
            if cat1 is None:
                return
            component = self._console_input(
                "Magnetization component (A/B/C..., blank = all): ", "")
            if component is None:
                return
            loaded = load_results(self.results_path, pattern=pattern or "*",
                                   carselect="d", cat1=cat1 or "*", numcomp=None,
                                   component=component or "*")

        # tx/ty/tz (segment ajuste, pour le trace sur un Zijderveld) ne sont
        # plus stockes dans le fichier .r (voir calcul.recompute_fit_geometry)
        # - recalcules ici a partir des mesures brutes (self.donnees), une
        # seule fois au chargement plutot qu'a chaque affichage - demande
        # explicite utilisateur ("to draw the line on the zijderveld plot we
        # can redo the calculation").
        new_matches = [recompute_fit_geometry(r, self.donnees) for r in loaded]
        # equivalent de selres (dataselect.f) : accumule sur la selection
        # existante, comme selmes (voir ouvrir_selection_dialog) - seul
        # initres/"Init results" remet self.results a zero. Demande
        # explicite utilisateur ("the selection of results should not
        # initialize the previous list").
        self.results = self.results + new_matches
        self._afficher(
            f"Selection: +{len(new_matches)} result(s) - "
            f"total {len(self.results)} result(s)\n")


if __name__ == "__main__":
    root = tk.Tk()
    # Force l'encodage systeme de Tcl a utf-8 - meme correctif que
    # AMS_Py/StereoUtils_Py (demande explicite utilisateur "l'appli
    # installee ne fonctionne pas tres bien, par exemple probleme de
    # texte Latin lors de l'importation. Pas de pb depuis le terminal") :
    # Tcl/Tk devine son "system encoding" depuis LANG/LC_ALL au
    # demarrage - absent pour un .app lance depuis le Finder/Dock
    # (contrairement a un Terminal, qui herite la locale du shell).
    root.tk.call("encoding", "system", "utf-8")
    app = STARpaleomagApp(root)
    root.mainloop()
