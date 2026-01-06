import sys
from pathlib import Path
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QFormLayout,
    QPushButton, QLineEdit, QFileDialog, QSplitter,
    QMessageBox
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction

# --- matplotlib ---
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

# --- IPython console (Spyder-style) ---
from qtconsole.rich_jupyter_widget import RichJupyterWidget
from qtconsole.manager import QtKernelManager

from jolymer.gui.regals.views.components import ComponentsDock, ComponentRow, ComponentsEditorWindow
from jolymer.gui.regals.views.uv import UVWindow
from jolymer.gui.regals.views.saxs import SAXSWindow
from jolymer.gui.regals.views.alignment import AlignmentWindow
from jolymer.gui.regals.state import BioREGALSState
from jolymer.gui.regals.console.ipython_widget import IPythonConsole
from jolymer.sas.bioREGALS import bioREGALS


# -----------------------------
# Matplotlib plot widget
# -----------------------------
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        self.fig, self.axes = plt.subplots(nrows=2, ncols=2, figsize=(5, 4))
        super().__init__(self.fig)


# -----------------------------
# Main Window
# -----------------------------
class RegalsMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.state = BioREGALSState()
        self.setWindowTitle("bioREGALS (demo GUI)")
        self.resize(1200, 800)
        # Central widget
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        # === LEFT: controls ===
        controls = QWidget()
        controls_layout = QFormLayout(controls)
        self.saxs_folder = QLineEdit()
        self.sample_name = QLineEdit("test_sample")
        btn_saxs = QPushButton("Browse…")
        btn_uv = QPushButton("Browse…")
        btn_run = QPushButton("Run (dummy)")
        btn_saxs.clicked.connect(self.select_saxs_folder)
        # btn_uv.clicked.connect(self.select_uv_file)
        btn_run.clicked.connect(self.run_dummy)
        controls_layout.addRow("X-ray folder:", self._with_button(self.saxs_folder, btn_saxs))
        # controls_layout.addRow("UV file:", self._with_button(self.uv_file, btn_uv))
        controls_layout.addRow("Sample name:", self.sample_name)
        controls_layout.addRow(btn_run)
        # === CENTER: plot ===
        self.canvas = MplCanvas()
        # === BOTTOM: IPython console ===
        self.console = IPythonConsole(
            namespace={
                "state": self.state,
            },
            # startup_file="~/.ipython/profile_default/startup/00-jolymer.py"
        )
        # === SPLITTERS ===
        vertical_splitter = QSplitter(Qt.Orientation.Vertical)
        vertical_splitter.addWidget(self.canvas)
        vertical_splitter.addWidget(self.console)
        vertical_splitter.setSizes([500, 300])
        horizontal_splitter = QSplitter(Qt.Orientation.Horizontal)
        horizontal_splitter.addWidget(controls)
        horizontal_splitter.addWidget(vertical_splitter)
        horizontal_splitter.setSizes([300, 900])
        main_layout.addWidget(horizontal_splitter)  # ← THIS WAS MISSING
        self.setup_ui()

    def setup_ui(self):
        self.components_dock = ComponentsDock(self.state)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea,
                           self.components_dock)
        self.components_dock.mixture_changed.connect(
            lambda: print("Mixture updated")
        )
        self._setup_menus()


    # -------------------------
    # Helpers
    # -------------------------

    def _setup_menus(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu("&File")
        # Save
        save_action = QAction("&Save state…", self)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.save_state)
        file_menu.addAction(save_action)
        # Load
        load_action = QAction("&Load state…", self)
        load_action.setShortcut("Ctrl+O")
        load_action.triggered.connect(self.load_state)
        file_menu.addAction(load_action)
        file_menu.addSeparator()
        # Quit (standard)
        quit_action = QAction("&Quit", self)
        quit_action.setShortcut("Ctrl+Q")
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)
        view_menu = self.menuBar().addMenu("View")
        for dock in [
            self.components_dock,
        ]:
            view_menu.addAction(dock.toggleViewAction())
        menu = self.menuBar().addMenu("Data")
        toolbar = self.addToolBar("Data")
        open_uv = QAction("UV", self)
        open_uv.triggered.connect(self.open_uv_window)
        menu.addAction(open_uv)
        toolbar.addAction(open_uv)
        open_saxs = QAction("SAXS", self)
        open_saxs.triggered.connect(self.open_saxs_window)
        menu.addAction(open_saxs)
        toolbar.addAction(open_saxs)
        open_align = QAction("ALIGN", self)
        open_align.triggered.connect(self.open_align_window)
        menu.addAction(open_align)
        toolbar.addAction(open_align)
        open_components_editor = QAction("COMPONENTS", self)
        open_components_editor.triggered.connect(self.open_components_editor_window)
        menu.addAction(open_components_editor)
        toolbar.addAction(open_components_editor)

    def _with_button(self, lineedit, button):
        w = QWidget()
        l = QHBoxLayout(w)
        l.setContentsMargins(0, 0, 0, 0)
        l.addWidget(lineedit)
        l.addWidget(button)
        return w


    # -------------------------
    # Slots
    # -------------------------
    def select_saxs_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select X-ray data folder")
        if folder:
            self.saxs_folder.setText(folder)

    def select_uv_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Select UV file")
        if fname:
            self.uv_file.setText(fname)

    def run_dummy(self):
        """Dummy action to show interaction"""
        from scipy.optimize import least_squares
        q = self.state.to_regals['q']
        x = self.state.to_regals['x']
        I = self.state.to_regals['I']
        sigma = self.state.to_regals['sigma']
        uv_meas = self.state.to_regals['uv_meas']
        def update_uv_scale(mix, Abs):
            import numpy as np
            from scipy.optimize import least_squares
            def resid(u):
                # linear combination of concentrations should match measured Abs
                return Abs - mix.concentrations @ u
            u0 = np.array([c.uv_scale for c in mix.components])
            sol = least_squares(resid, u0)
            # update components
            for k, comp in enumerate(mix.components):
                comp.uv_scale = sol.x[k]
            return mix
        stop_fun = lambda num_iter, params: [num_iter >= 3000, 'max_iter']
        def update_fun(num_iter, new_mix, params, resid):
            if num_iter % 500 == 0:
                print(num_iter)
            update_uv_scale(new_mix, uv_meas)
            return True  # continue as normal
        MM = self.state.mixture
        MM.lambda_profile = [1 * lam for lam in MM.estimate_profile_lambda(sigma)]
        MM.lambda_concentration = [1 * lam for lam in MM.estimate_concentration_lambda(sigma)]
        RR = bioREGALS(I,sigma)
        self.state.bioREGALS = RR
        MM = RR.fit_concentrations(MM)
        [MM,params,resid,exit_cond] = RR.run(MM, stop_fun, update_fun);
        self.state.mixture = MM
        self.state.regals_result = [params, resid, exit_cond]
        self.update_plot()

    def update_plot(self):
        # self.canvas.ax.clear()
        try:
            import numpy as np
            import pandas as pd
            q = self.state.to_regals['q']
            x = self.state.to_regals['x']
            I = self.state.to_regals['I']
            sigma = self.state.to_regals['sigma']
            fig, axs = self.canvas.fig, self.canvas.axes
            MM =  self.state.mixture
            params, resid, exit_cond = self.state.regals_result
            uv_scales = [comp.uv_scale for comp in MM.components]
            axs[0, 0].plot(x, MM.uv_meas, label="UV")
            axs[0, 0].plot(x, MM.concentrations*np.array(uv_scales))
            axs[0,0].legend()
            #chi2 vs x
            axs[1, 0].plot(x, np.mean(resid ** 2, 0))
            axs[1, 0].set_xlabel('$x$')
            axs[1, 0].set_ylabel('$\chi^2$')
            axs[1, 0].set_ylim(0,10)
            #profiles vs q
            axs[0, 1].semilogx(q, MM.profiles)
            axs[0, 1].set_xlabel('$q (Å^{-1})$')
            axs[0, 1].set_ylabel('Profiles')
            #extracted profiles
            for i in range(len(MM.components)):
                [Ii,sigmai] = MM.extract_profile(I,sigma,i);
                df = pd.DataFrame({'q':q,
                                   'I': Ii,
                                   'err_I': sigmai})
                # m.save_data(f'regals_output/{NAME}-P{i}.dat', df=df)
                m = self.state.saxs.ms[0]
                Rgdict = m.get_rg(df=df, qmin=0.02, qmax=0.1)
                Rg = Rgdict['Rg']
                print(f'P{i}: I_max = ', Ii.max())
                axs[1, 1].errorbar(q, Ii/Ii.max(), sigmai,
                                   label=f'P{i}; rg={Rg:.2f} A')
            axs[1, 1].set_xscale('log')
            axs[1, 1].set_yscale('log')
            axs[1, 1].set_xlabel('$q (Å^{-1})$')
            axs[1, 1].set_ylabel('Extracted profiles')
            axs[1,1].legend()
            self.canvas.draw()
            # Print to IPython console
            self.console._append_plain_text(
                f"\nRan dummy fit for {self.sample_name.text()}\n"
            )
        except Exception as e:
            QMessageBox.critical(self, "Plot failed", str(e))


    def open_uv_window(self):
        if not hasattr(self, "_uv_window"):
            self._uv_window = UVWindow(self.state, self)
        self._uv_window.show()
        self._uv_window.raise_()

    def open_saxs_window(self):
        if not hasattr(self, "_saxs_window"):
            self._saxs_window = SAXSWindow(self.state, self)
        self._saxs_window.show()
        self._saxs_window.raise_()

    def open_align_window(self):
        if not hasattr(self, "_align_window"):
            self._align_window = AlignmentWindow(self.state, self)
        self._align_window.show()
        self._align_window.raise_()

    def open_components_editor_window(self):
        if not hasattr(self, "_components_editor_window"):
            self._components_editor_window = ComponentsEditorWindow(self.state, self)
        self._components_editor_window.show()
        self._components_editor_window.raise_()

    def save_state(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save REGALS state", "", "REGALS state (*.json)"
        )
        if not path:
            return
        try:
            self.state.to_json(Path(path))
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))

    def load_state(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load REGALS state", "", "REGALS state (*.json)"
        )
        if not path:
            return

        try:
            print('test 3')
            self.state.load_from_json(Path(path))
            print('test 3')
            # self.state.rebuild_from_state()
        except Exception as e:
            QMessageBox.critical(self, "Load failed", str(e))


# -----------------------------
# Entry point
# -----------------------------
def main():
    app = QApplication(sys.argv)
    win = RegalsMainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

