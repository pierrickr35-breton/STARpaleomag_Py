#!/bin/bash
# Double-click to start STARpaleomag_Py from this folder.
# (if macOS refuses the first time: right-click > Open)
cd "$(dirname "$0")" || exit 1
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
    echo "No python3 with the required libraries was found. Install them with:"
    echo "    python3 -m pip install numpy matplotlib pandas scipy pmagpy"
    read -r -p "Press Return to close."
    exit 1
fi
exec "$PY" app.py
