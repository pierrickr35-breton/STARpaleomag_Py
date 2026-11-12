# STARpaleomag_Py

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24170591.svg)](https://doi.org/10.5281/zenodo.24170591)

Paleomagnetism and rock-magnetism data analysis (demagnetization diagrams, component fitting, site means, paleointensity, anisotropy of remanence, MagIC import/export) — Python port of the Starmac software.

STARpaleomag_Py is a plain **Python 3** program (Tkinter + Matplotlib). It runs on **Windows, Linux and macOS**. For Macs with an Apple Silicon processor, a ready-to-use application is also provided.

## Running on any platform (Windows, Linux, macOS)

1. Install [Python 3](https://www.python.org/downloads/) if it is not already on your computer. On Linux, Tkinter may be a separate package (for example `sudo apt install python3-tk`).
2. Get the program. **Recommended: with git**, which makes updates a double-click (see below). Open the Terminal, go to the folder where you want the program (for example `Pmag_PR_Python` in your home folder) and clone it:

   ```bash
   mkdir -p ~/Pmag_PR_Python
   cd ~/Pmag_PR_Python
   git clone https://github.com/pierrickr35-breton/STARpaleomag_Py.git
   ```

   This creates `~/Pmag_PR_Python/STARpaleomag_Py`. On a Mac, the first `git` command may offer to install the Apple command line tools: accept, then run it again.
   Without git, you can also **[download the source code (zip)](https://github.com/pierrickr35-breton/STARpaleomag_Py/archive/refs/heads/main.zip)** and unzip it, but it then has to be downloaded again for every update.
3. Install the required libraries (once):

   ```bash
   pip install numpy matplotlib pandas scipy pmagpy
   ```

4. Start the program from its folder:

   ```bash
   python3 app.py
   ```

   (on Windows: `python app.py`)

To update later: `git pull` in the program folder (on a Mac, double-click `Update_STARpaleomag_Py.command`, see below), or download the zip again. The libraries do not need to be reinstalled.

## macOS: Dock shortcut and updates

Once the program is installed with git (step 2), two scripts in its folder are all you need. Double-click them in the Finder (if macOS refuses the first time: right-click > *Open*, then *Open*).

- **`Run_STARpaleomag_Py.command`** starts the program. Drag it to the **right part of the Dock** (next to the Trash) to start the program from there. It runs the files of this folder, so the Dock always starts the latest version. The program icon is in `resources/AppIcon.icns`: select `Run_STARpaleomag_Py.command`, *File > Get Info*, and drag `AppIcon.icns` onto the small icon at the top left of the window.
- **`Update_STARpaleomag_Py.command`** downloads the latest version from GitHub. If you modified program files, it offers to set your changes aside (`git stash`) first; your data files in the folder are not touched. It then lists what changed and warns if a new library is needed.

The same update in the Terminal is a single command:

```bash
git -C ~/Pmag_PR_Python/STARpaleomag_Py pull
```

and, to update all the programs installed in that folder at once:

```bash
for d in ~/Pmag_PR_Python/*/; do git -C "$d" pull; done
```

## macOS: ready-to-use application

For Macs with an **Apple Silicon** processor (M1 and later), a standalone application is available: Python and all libraries are included, there is nothing else to install.

1. **[Go to the latest release](https://github.com/pierrickr35-breton/STARpaleomag_Py/releases/latest)** and, under *Assets*, click `STARpaleomag_Py.zip`.
2. Double-click the zip to unpack it, and move `STARpaleomag_Py.app` to your *Applications* folder.
3. The application is not signed with an Apple developer certificate, so macOS blocks it the first time. Try to open it once, then go to *System Settings > Privacy & Security*, scroll down and click **Open Anyway**. This is needed only once.

   If macOS says the application "is damaged and can't be opened", open the Terminal and type:

   ```bash
   xattr -dr com.apple.quarantine /Applications/STARpaleomag_Py.app
   ```

On an Intel Mac, use the Python instructions above.

## User Guide

In the program: *Help > User Guide*. The same guide is in this repository: `help/STARpaleomag_Py_Guide.html`.

## License

MIT — see `LICENSE`.

## How to cite

Roperch, P. (2026). *STARpaleomag_Py* [Software]. Zenodo. https://doi.org/10.5281/zenodo.24170591

This DOI covers all versions; please also give the version you used (each release has its own DOI on Zenodo). GitHub's "Cite this repository" button gives the same reference in APA and BibTeX.
