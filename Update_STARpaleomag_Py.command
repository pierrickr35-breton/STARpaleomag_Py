#!/bin/bash
# Double-click to update STARpaleomag_Py from GitHub (git pull).
NAME="STARpaleomag_Py"
cd "$(dirname "$0")" || exit 1
echo "Updating $NAME in $(pwd)"
echo
finish() { echo; read -r -p "Press Return to close."; exit "$1"; }
if [ ! -d .git ]; then
    echo "This folder was not installed with git (zip download?): it cannot be updated"
    echo "automatically. Install it once with git (see README.md), or download the"
    echo "zip again from https://github.com/pierrickr35-breton/$NAME"
    finish 1
fi
# on a Mac without the developer tools, the first git call offers to install them
if ! git --version >/dev/null 2>&1; then
    echo "git is not available. Install the Apple command line tools (a window may have"
    echo "opened), or type in the Terminal:  xcode-select --install   then try again."
    finish 1
fi
git fetch --quiet origin || { echo "Could not reach GitHub (network?)."; finish 1; }
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    echo "These program files were modified on this computer:"
    git status --short --untracked-files=no
    echo
    read -r -p "Set these changes aside (git stash, recoverable) and update? [y/N] " a
    case "$a" in
        [yYoO]*) git stash push --quiet -m "before update $(date '+%Y-%m-%d %H:%M')" || finish 1
                 echo "Changes saved: 'git stash list' / 'git stash pop' to get them back." ;;
        *) echo "Nothing changed."; finish 1 ;;
    esac
fi
BEFORE=$(git rev-parse HEAD)
git merge --ff-only --quiet origin/main || { echo "Update failed (see message above)."; finish 1; }
if [ "$BEFORE" = "$(git rev-parse HEAD)" ]; then
    echo "Already up to date."
else
    echo "Changes:"
    git log --format="  %cd  %s" --date=short "$BEFORE..HEAD"
fi
echo
git log -1 --format="Version now: %h of %cd" --date=short
# python3 that has the libraries (Finder does not use the Terminal PATH,
# so the full path is looked for here once)
PY=""
for P in "$(command -v python3)" /opt/homebrew/bin/python3 /usr/local/bin/python3 \
         /Library/Frameworks/Python.framework/Versions/Current/bin/python3 /usr/bin/python3; do
    if [ -n "$P" ] && [ -x "$P" ] && "$P" -c "import tkinter, numpy, matplotlib, pandas, scipy, pmagpy" >/dev/null 2>&1; then
        PY="$P"; break
    fi
done
if [ -z "$PY" ]; then
    echo
    echo "Some libraries are missing. Install them with:"
    echo "    python3 -m pip install numpy matplotlib pandas scipy pmagpy"
fi
finish 0
