python -m nuitka `
    --onefile `
    --enable-plugin=pyqt6 `
    --include-package=jolymer `
    --include-package-data=jolymer `
    --output-dir=build `
    --output-filename=jolymer-secplot.exe `
    src/jolymer/gui/sec_plot/app.py
