python -m nuitka `
    --onefile `
    --enable-plugin=pyqt6 `
    --include-package=jolymer `
    --output-dir=build `
    --output-filename=jolymer-secplot.exe `
    src/jolymer/gui/app.py

