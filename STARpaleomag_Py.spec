# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_data_files

# pmagpy ships non-.py data files (IGRF coefficients in field_models/,
# the MagIC data model in data_model/) that orient_sample.py (IGRF
# declination) and paleointensity_magic.py depend on at runtime - PyInstaller
# does not bundle these automatically, only .py modules. Both call sites
# were added after the previous build (2026-08-22), so this was never
# exercised in a packaged app before.
datas = collect_data_files('pmagpy')
# Guide utilisateur statique (Help > User Guide, voir app._resource_path/
# ouvrir_user_guide) - demande explicite utilisateur ("aide en ligne...
# OK pour la 3", l'option sans cle API ni cout recurrent). Doit etre
# EXTRAIT sous le meme nom de dossier ('help/') pour que _resource_path
# (sys._MEIPASS + 'help' + nom de fichier) le retrouve une fois empaquete.
datas += [('help', 'help')]

# matplotlib charge les backends SVG/PDF dynamiquement au moment de
# fig.savefig(path, format="svg") (Export SVG...) / pdf.savefig(fig)
# (detailed_export.py), PAS via un `import` statique visible dans le
# source - PyInstaller ne les detecte donc pas tout seuls et l'app
# packagee echoue a l'export avec "No module named matplotlib.backends.
# backend_svg" (repere sur un vrai .app construit, fonctionnait depuis
# les sources).
hiddenimports = ['matplotlib.backends.backend_svg', 'matplotlib.backends.backend_pdf']

a = Analysis(
    ['app.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='STARpaleomag_Py',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='STARpaleomag_Py',
)
app = BUNDLE(
    coll,
    name='STARpaleomag_Py.app',
    icon='resources/AppIcon.icns',
    bundle_identifier=None,
)
