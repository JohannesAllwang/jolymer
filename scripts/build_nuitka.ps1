python -m nuitka `
    --standalone `
    --enable-plugin=pyqt6 `
    --output-dir=build `
    --include-package=jolymer `
    src/jolymer/gui/app.py
