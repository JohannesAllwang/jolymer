from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QPushButton,
    QFileDialog, QFormLayout, QDoubleSpinBox,
    QSpinBox, QHBoxLayout, QLineEdit,
    QStyle
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import numpy as np

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
        path_layout = QHBoxLayout()
        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText("Select data file or directory…")
        self.path_btn = QPushButton()
        self.path_btn.setIcon(self.style().standardIcon(QStyle.SP_DirOpenIcon))
        self.path_btn.setText("")
        self.path_btn.setFixedWidth(30)
        self.path_btn.clicked.connect(self.browse_path)
        path_layout.addWidget(self.path_edit)
        path_layout.addWidget(self.path_btn)
        self.refwl = QDoubleSpinBox()
        self.refwl.setRange(0, 2000)
        self.outwl_from = QDoubleSpinBox()
        self.outwl_from.setRange(0, 2000)
        self.outwl_to = QDoubleSpinBox()
        self.outwl_to.setRange(0, 2000)
        self.alignment_time = QDoubleSpinBox()
        self.alignment_time.setMaximum(1e6)
        self.mintime = QDoubleSpinBox()
        self.mintime.setMaximum(1e9)
        self.maxtime = QDoubleSpinBox()
        self.maxtime.setMaximum(1e9)
        form.addRow("File path", path_layout)
        form.addRow("refwl", self.refwl)
        form.addRow("outwl", self.outwl_from)
        form.addRow("to", self.outwl_to)
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
            self.outwl_from.setValue(280)
            self.outwl_to.setValue(280)
            self.maxtime.setValue(1e6)
        else:
            uv = self.state.uv
            # hydrate from state
            if hasattr(uv, "refwl") and uv.refwl is not None:
                self.refwl.setValue(uv.refwl)
            if hasattr(uv, "outwl") and uv.outwl is not None:
                if not isinstance(uv.outwl, (list, tuple, np.ndarray)):
                    self.outwl_from.setValue(uv.outwl)
                    self.outwl_to.setValue(uv.outwl)
                else:
                    self.outwl_from.setValue(uv.outwl[0])
                    self.outwl_to.setValue(uv.outwl[-1])
            if hasattr(uv, "alignment_time") and uv.alignment_time is not None:
                self.alignment_time.setValue(uv.alignment_time)
            if hasattr(uv, "mintime") and uv.mintime is not None:
                self.mintime.setValue(uv.mintime)
            if hasattr(uv, "maxtime") and uv.maxtime is not None:
                self.maxtime.setValue(uv.maxtime)
            if hasattr(uv, "filename") and uv.filename is not None:
                self.path_edit.setText(str(Path(uv.path, uv.spec_filename)))
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
            },
            # startup_file="~/.ipython/profile_default/startup/00-jolymer.py"
        )
        layout.addWidget(self.console)

    def browse_path(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select file",
            "",
            "All files (*)"
        )
        if path:
            self.path_edit.setText(path)

    # ---------------------------
    def load_uv(self):
        try:
            path = Path(self.path_edit.text())
            if not path.exists():
                raise ValueError("Invalid file path")
            wl_from = self.outwl_from.value()
            wl_to = self.outwl_to.value()
            if np.isclose(wl_from, wl_to):
                outwl = wl_from
            else:
                outwl = list(range(
                    int(round(wl_from)),
                    int(round(wl_to)) + 1
                ))
            uv = onlineUV(
                spec_filename=path.name,
                path=str(path.parent),
                refwl=self.refwl.value(),
                outwl=outwl,
                alignment_time=float(self.alignment_time.value()),
                mintime=self.mintime.value(),
                maxtime=self.maxtime.value(),
            )
            uv.get_data()
            self.state.uv = uv
            self.update_plot()
        except Exception as e:
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

