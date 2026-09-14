"""
Export de resultats STARpaleomag_Py vers StereoUtils_Py, au format "Stereo
Project" a TROIS BLOCS - demande explicite utilisateur ("can you manage a
stereo project with three blocks (individual directions; mean directions,
and VGP, each block with its own labels. a column for the tilt correction
0 to 100 and an information column to be filled manually by the user").

Historique (voir StereoUtils_Py/stereo_project.py pour le format "Project"
D'ORIGINE, toujours lu/trace tel quel par StereoUtils_Py) : une premiere
version ecrivait deux fichiers separes (un "Project" plat sans en-tete pour
directions/moyennes, un fichier VGP a part pour "Plot VGPs on Map") -
demande explicite utilisateur initiale ("is it possible to add the export
results to Stereo_Py as in the original Fortran with in addition a specific
file for the poles", precisee : "with a format compatible to the project in
Stereo_Py"), puis affinee ("I would prefer the export to Stereo in the form
of the project where user can choose manually the symbol and color... It
will be especially usefull for the results in the same export file, we can
have lines, plane, mean direction and mean VGP"). Remplacee ICI par UN SEUL
fichier, 3 blocs marques par un en-tete "#" auto-descriptif (meme
convention que .pmagani - #specimen/#site mean tensor results - et que
`stereo_selection.split_header`/`detect_header`, deja utilises ailleurs
dans ce port) :

    # individual directions
    layer  id  type  dec  inc  alpha95  tilt_correction  symbol  color  size  info
    # mean directions
    layer  id  dec  inc  alpha95  n  tilt_correction  symbol  color  size  info
    # VGP
    site  id  paleolon  paleolat  dp  dm  n  a95  site_lat  site_lon  tilt_correction  symbol  color  size  info

- Bloc "individual directions" : UNE ligne par resultat de specimen
  (cat1 in L/P/f/s), `type` d(irection)/g(rand cercle) selon cat1=='P' ou
  non (voir _line_or_plane_row) - meme raisonnement que l'ancien
  export_results_to_stereo : un PLAN doit rester un vrai grand cercle
  ('g'), jamais seulement le point de son pole. `layer` = site (6 premiers
  caracteres du specimen, meme convention que magic_export._site_mean_row).
- Bloc "mean directions" : UNE ligne par moyenne de site ("mean:",
  cat1=='F').
- Bloc "VGP" : UNE ligne par moyenne de site QUI PORTE un VGP (par4/par5
  renseignes) - paleolon/paleolat, pas dec/inc (point GEOGRAPHIQUE, jamais
  une direction de stereonet - voir StereoUtils_Py/stereo_pmagpy.
  plot_vgp_project). `site_lat`/`site_lon` (r.lat/r.rlong - coordonnees du
  SITE, PAS du pole) AJOUTES - demande explicite utilisateur ("is it
  possible to plot the VGP with their dp,dm ellipse") : la vraie ellipse
  dp/dm (asymetrique, orientee le long du meridien site->pole - voir
  Butler 1992 fig. A.2) a besoin de la position du SITE pour calculer
  cette orientation (angle de rotation via la loi des cosinus spherique
  sur le triangle site/pole nord/paleopole) - PAS calculable depuis le
  seul VGP. Sans site connu (0.0/0.0, sentinelle documentee par
  calcul.build_site_mean_result quand le site n'a pas de lat/lon), voir
  StereoUtils_Py/stereo_pmagpy.plot_vgp_project pour le repli sur un
  simple cercle de rayon p95, deja en place avant cet ajout.

`tilt_correction` : 0/100 (in-situ/apres pendage complet), MEME convention
MagIC dir_tilt_correction que calcul._ORIENT_TO_FILE_CODE (colonne "IS/TC"
de .pmagres) et AMS_Py/ams_selection._ORIENT_TO_FILE_CODE (.pmagani, "#site
mean tensor results") - demande explicite utilisateur ("for the mean, can
we use 0 and 100 for IS and full TC like in magic", reprise ici pour rester
sur UNE SEULE convention 0-100 a travers tout le projet plutot que d'en
inventer une 4e). "1" (coordonnees echantillon) reste gere en ecriture pour
un `orientation=1` explicite, meme s'il n'est normalement pas utilise pour
une moyenne.

`info` : TOUJOURS vide a l'export - demande explicite utilisateur ("an
information column to be filled manually by the user") - colonne reservee
a l'utilisateur, jamais pre-remplie automatiquement (le fichier reste un
fichier texte ordinaire, directement editable).

Couleurs par defaut (retouchables a la main ensuite, comme pour
`symbol`/`color`/`size`) : "black" pour les directions/grands cercles de
specimen (nom de couleur, PAS un triplet "0_0_0" - StereoUtils_Py accepte
directement les noms de NAMED_COLORS, plus lisible/plus facile a retoucher
a la main - demande explicite utilisateur : "replace 0_0_0 by black as a
default"), "red" pour les moyennes de site, "blue" pour les VGP - memes
couleurs que les versions precedentes de cet export (colonne `rgb`
renommee `color` ici - demande explicite utilisateur : "can you change rgb
by color" - les valeurs elles-memes restent utilisables telles quelles,
nom ou triplet "R_G_B", voir stereo_project.decode_color), gardees pour
rester visuellement distinctes d'emblee."""

from typing import List

from calcul import FitResult, _correct_dec_inc, _ORIENT_TO_FILE_CODE, dp_dm_from_a95

_BLOCK_HEADERS = {
    "directions": ["layer", "id", "type", "dec", "inc", "alpha95",
                   "tilt_correction", "symbol", "color", "size", "info"],
    "means": ["layer", "id", "dec", "inc", "alpha95", "n",
              "tilt_correction", "symbol", "color", "size", "info"],
    "vgp": ["site", "id", "paleolon", "paleolat", "dp", "dm", "n", "a95",
            "site_lat", "site_lon", "tilt_correction", "symbol", "color", "size", "info"],
}


def _tilt_code(orientation: int) -> str:
    return _ORIENT_TO_FILE_CODE.get(float(orientation), "0")


def export_stereo_project(results: List[FitResult], out_path: str, orientation: int = 2) -> dict:
    """Ecrit UN fichier "Stereo Project" a 3 blocs (voir docstring module).
    Retourne {"directions": n, "means": n, "vgp": n} - le nombre de lignes
    ecrites dans chaque bloc (0 si aucun resultat exploitable pour ce
    bloc ; les 3 compteurs peuvent etre 0 en meme temps, le fichier n'est
    alors pas ecrit)."""
    tilt = _tilt_code(orientation)

    direction_rows = []
    for r in results:
        if r.id[:5] == "mean:" or r.cat1 not in ("L", "P", "f", "s"):
            continue
        site = r.id[:6]
        dec, inc = _correct_dec_inc(r, orientation)
        etype = "g" if r.cat1 == "P" else "d"
        size = "0.20" if r.cat1 == "P" else "0.30"
        # alpha95 : jamais 0.0 pour un plan ('g') - stereo_project.
        # load_project (StereoUtils_Py) retombe silencieusement sur
        # type='d' des que alpha95==0.0, quel que soit le type ecrit dans
        # le fichier (quirk du Fortran d'origine, reproduit tel quel) ;
        # r.mad (jamais exactement nul pour un ajustement reel) sert de
        # repli, inutilise par le trace du grand cercle lui-meme.
        alph = f"{r.mad:.1f}" if etype == "g" else "0.0"
        direction_rows.append([
            site, r.id, etype, f"{dec:.1f}", f"{inc:.1f}", alph,
            tilt, "c", "black", size, "",
        ])

    mean_rows = []
    vgp_rows = []
    for r in results:
        if r.id[:5] != "mean:":
            continue
        site = r.id[6:].strip()
        mean_rows.append([
            "Site_Means", site, f"{r.dec:.1f}", f"{r.inc:.1f}", f"{r.mad:.1f}", str(r.nb),
            tilt, "e", "red", "0.55", "",
        ])
        if r.par4 or r.par5:
            vgp_dp, vgp_dm = (
                (r.vgp_dp, r.vgp_dm) if (r.vgp_dp or r.vgp_dm)
                else dp_dm_from_a95(r.mad, r.inc)
            )
            vgp_rows.append([
                site, site, f"{r.par5:.1f}", f"{r.par4:.1f}",  # paleolon=par5, paleolat=par4
                f"{vgp_dp:.1f}", f"{vgp_dm:.1f}", str(r.nb), f"{r.mad:.1f}",
                f"{r.lat:.4f}", f"{r.rlong:.4f}",  # site_lat/site_lon - voir docstring module
                tilt, "c", "blue", "0.50", "",
            ])

    counts = {"directions": len(direction_rows), "means": len(mean_rows), "vgp": len(vgp_rows)}
    if not any(counts.values()):
        return counts

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("# individual directions\n")
        f.write("#" + "\t".join(_BLOCK_HEADERS["directions"]) + "\n")
        for row in direction_rows:
            f.write("\t".join(row) + "\n")
        f.write("\n# mean directions\n")
        f.write("#" + "\t".join(_BLOCK_HEADERS["means"]) + "\n")
        for row in mean_rows:
            f.write("\t".join(row) + "\n")
        f.write("\n# VGP\n")
        f.write("#" + "\t".join(_BLOCK_HEADERS["vgp"]) + "\n")
        for row in vgp_rows:
            f.write("\t".join(row) + "\n")
    return counts
