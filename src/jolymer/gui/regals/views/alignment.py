from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QPushButton,
    QFileDialog, QFormLayout, QDoubleSpinBox
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from jolymer.uv.onlineUV import onlineUV
from jolymer.gui.regals.state import BioREGALSState
from jolymer.gui.regals.console.ipython_widget import IPythonConsole

from pathlib import Path

class AlignmentWindow(QMainWindow):
    def __init__(self, state: BioREGALSState, parent=None):
        super().__init__(parent)
        self.state = state
        self.setWindowTitle("Alignment")
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        # ---------------------------
        # Controls
        # ---------------------------
        form = QFormLayout()
        self.shift0 = QDoubleSpinBox()
        self.shift0.setRange(-2000, 2000)
        self.shift0.setValue(-300)
        self.scale0 = QDoubleSpinBox()
        self.scale0.setRange(0, 10)
        self.scale0.setValue(0.946)
        form.addRow("shift0", self.shift0)
        form.addRow("scale0", self.scale0)
        layout.addLayout(form)
        self.align_btn = QPushButton("Align")
        self.align_btn.clicked.connect(self.align)
        layout.addWidget(self.align_btn)
        # ---------------------------
        # Plot
        # ---------------------------
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasQTAgg(self.fig)
        layout.addWidget(self.canvas)
        # ---------------------------
        # Console
        # ---------------------------
        self.console = IPythonConsole(
            namespace={
                "state": self.state,
                "cm_saxs": self.state.cm_saxs,
            },
            # startup_file="~/.ipython/profile_default/startup/00-jolymer.py"
        )
        layout.addWidget(self.console)

    # ---------------------------
    def align(self):
        try:
            scale0 = self.scale0.value()
            shift0 = self.shift0.value()
            cm_saxs = self.state.cm_saxs
            cm_saxs.align_uv_to_saxs(
                    qmin=cm_saxs._alignment['qmin'],
                    qmax=cm_saxs._alignment['qmax'],
                    qstar=cm_saxs._alignment['qstar'],
                    saxs_kind=cm_saxs._alignment['saxs_kind'],
                    load=True, scale0shift0=[scale0,shift0])
            self.scale0.setValue(cm_saxs._alignment['scale'])
            self.shift0.setValue(cm_saxs._alignment['shift'])
            cm_saxs.interpolate_uv()
            self.update_plot()
            # self.console.bind_uv(uv, self.fig, self.ax)
        except Exception as e:
            # Print to console if you have one
            if hasattr(self, "console") and self.console is not None:
                self.console._append_plain_text(f"Error loading UV: {e}")
            else:
                print(f"Error loading UV: {e}")

    def update_plot(self):
        try:
            if self.state.cm_saxs is None:
                return

            self.ax.clear()
            self.state.cm_saxs.plot_alignment(ax=self.ax)
            # self.ax.set_xlabel("time [s]")
            # self.ax.set_ylabel("Abs")
            self.canvas.draw_idle()
            self.state.to_regals = self.state.cm_saxs.to_regals()
            self.console._append_plain_text(str(self.state.to_regals))
        except Exception as e:
            if hasattr(self, "console") and self.console is not None:
                self.console._append_plain_text(f"Error updating plot: {e}")
            else:
                print(f"Error updating plot: {e}")


