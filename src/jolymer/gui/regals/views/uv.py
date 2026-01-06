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

class UVWindow(QMainWindow):
    def __init__(self, state: BioREGALSState, parent=None):
        super().__init__(parent)
        self.state = state
        self.setWindowTitle("UV")
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        # ---------------------------
        # Controls
        # ---------------------------
        form = QFormLayout()
        self.refwl = QDoubleSpinBox()
        self.refwl.setRange(0, 2000)
        self.outwl = QDoubleSpinBox()
        self.outwl.setRange(0, 2000)
        self.alignment_time = QDoubleSpinBox()
        self.alignment_time.setMaximum(1e6)
        self.mintime = QDoubleSpinBox()
        self.mintime.setMaximum(1e9)
        self.maxtime = QDoubleSpinBox()
        self.maxtime.setMaximum(1e9)
        form.addRow("refwl", self.refwl)
        form.addRow("outwl", self.outwl)
        form.addRow("alignment_time", self.alignment_time)
        form.addRow("mintime", self.mintime)
        form.addRow("maxtime", self.maxtime)
        layout.addLayout(form)
        self.load_btn = QPushButton("Load UV (.spc)")
        self.load_btn.clicked.connect(self.load_uv)
        layout.addWidget(self.load_btn)
        if self.state.uv is None:
            # first-time defaults
            self.refwl.setValue(280)
            self.outwl.setValue(280)
            self.maxtime.setValue(1e6)
        else:
            uv = self.state.uv
            # hydrate from state
            if hasattr(uv, "refwl") and uv.refwl is not None:
                self.refwl.setValue(uv.refwl)
            if hasattr(uv, "outwl") and uv.outwl is not None:
                self.outwl.setValue(uv.outwl)
            if hasattr(uv, "alignment_time") and uv.alignment_time is not None:
                self.alignment_time.setValue(uv.alignment_time)
            if hasattr(uv, "mintime") and uv.mintime is not None:
                self.mintime.setValue(uv.mintime)
            if hasattr(uv, "maxtime") and uv.maxtime is not None:
                self.maxtime.setValue(uv.maxtime)
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
    def load_uv(self):
        fn, _ = QFileDialog.getOpenFileName(
            self, "Load UV file", "", "SPC (*.spc)"
        )
        if not fn:
            return
        try:
            uv = onlineUV(
                spec_filename=fn.split("/")[-1],
                path=str(Path(fn).parent),
                refwl=self.refwl.value(),
                outwl=self.outwl.value(),
                alignment_time=int(self.alignment_time.value()),
                mintime=self.mintime.value(),
                maxtime=self.maxtime.value(),
            )
            uv.get_data()
            self.state.uv = uv
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
            if self.state.uv is None:
                return

            self.ax.clear()
            df = self.state.uv.get_scaled_Abs(show=False)
            self.ax.plot(df.time, df.Abs)
            self.ax.set_xlabel("time [s]")
            self.ax.set_ylabel("Abs")
            self.canvas.draw_idle()
        except Exception as e:
            if hasattr(self, "console") and self.console is not None:
                self.console._append_plain_text(f"Error updating plot: {e}")
            else:
                print(f"Error updating plot: {e}")

