from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QPushButton,
    QFileDialog, QFormLayout, QDoubleSpinBox,
    QLineEdit, QLabel, QSpinBox, QHBoxLayout,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from jolymer.uv.onlineUV import onlineUV
from jolymer.sas.CoupledMeasurement import CoupledMeasurement, ms_from_folder
from jolymer.gui.regals.state import BioREGALSState
from jolymer.gui.regals.console.ipython_widget import IPythonConsole

from pathlib import Path

class SAXSWindow(QMainWindow):
    def __init__(self, state: BioREGALSState, parent=None):
        super().__init__(parent)
        self.state = state
        self.setWindowTitle("SAXS")
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        form = QFormLayout()
        # --- Folder and load button ---
        self.path_edit = QLineEdit()
        self.name_edit = QLineEdit()
        self.load_btn = QPushButton("Load SAXS frames")
        self.load_btn.clicked.connect(self.load_saxs)
        layout.addWidget(QLabel("SAXS folder:"))
        prepath = "/path/to/datafolder/"
        if not self.state.saxs is None:
            prepath = str(self.state.saxs.ms[0].path)
        self.path_edit.setText(prepath)
        layout.addWidget(self.path_edit)
        prename = "filename_prefix"
        if not self.state.saxs is None:
            prename = str(self.state.saxs.name)
        layout.addWidget(QLabel("File prefix:"))
        self.name_edit.setText(prename)
        layout.addWidget(self.name_edit)
        layout.addWidget(self.load_btn)

        # --- Frame selection ---
        self.start_idx = QSpinBox()
        self.start_idx.setRange(0, 20000)
        self.end_idx = QSpinBox()
        self.end_idx.setRange(0, 20000)
        frame_layout = QHBoxLayout()
        frame_layout.addWidget(QLabel("Frames from"))
        frame_layout.addWidget(self.start_idx)
        frame_layout.addWidget(QLabel("to"))
        frame_layout.addWidget(self.end_idx)
        layout.addLayout(frame_layout)

        # --- SAXS kind and q parameters ---
        self.saxs_kind_edit = QLineEdit("I0")
        self.q_beamstop_edit = QDoubleSpinBox()
        self.q_beamstop_edit.setDecimals(3)
        self.q_beamstop_edit.setSingleStep(0.0001)
        self.qstar_edit = QDoubleSpinBox()
        self.qstar_edit.setDecimals(3)
        self.qstar_edit.setSingleStep(0.01)
        self.qmin_edit = QDoubleSpinBox()
        self.qmin_edit.setDecimals(3)
        self.qmin_edit.setSingleStep(0.01)
        self.qmax_edit = QDoubleSpinBox()
        self.qmax_edit.setDecimals(3)
        self.qmax_edit.setSingleStep(0.01)

        if self.state.cm_saxs is None:
            self.q_beamstop_edit.setValue(0.006)
            self.qmax_edit.setValue(0.4)
            self.qmin_edit.setValue(0.04)
            self.qstar_edit.setValue(0.1)
        else:
            cm_saxs = self.state.cm_saxs
            self.q_beamstop_edit.setValue(cm_saxs.saxs_list.ms[0].qmin)
            self.qmax_edit.setValue(cm_saxs.qmax)
            self.qmin_edit.setValue(cm_saxs.qmin)
            self.qstar_edit.setValue(cm_saxs.qstar)
        self.start_idx.setValue(0)
        self.end_idx.setValue(500)


        form_layout = QFormLayout()
        form_layout.addRow("q beamstop:", self.q_beamstop_edit)
        form_layout.addRow("SAXS kind:", self.saxs_kind_edit)
        form_layout.addRow("q*:", self.qstar_edit)
        form_layout.addRow("q min:", self.qmin_edit)
        form_layout.addRow("q max:", self.qmax_edit)
        layout.addLayout(form_layout)

        # --- Plot canvas ---
        self.canvas = FigureCanvasQTAgg(Figure(figsize=(5,3)))
        self.ax = self.canvas.figure.add_subplot(111)
        layout.addWidget(self.canvas)

        # --- Prepare CoupledMeasurement ---
        self.export_btn = QPushButton("Prepare CoupledMeasurement")
        self.export_btn.clicked.connect(self.prepare_coupled)
        layout.addWidget(self.export_btn)
        self.console = IPythonConsole(
            namespace={
                "state": self.state,
            },
        )
        layout.addWidget(self.console)



    # -------------------------
    # Slots
    # -------------------------
    def load_saxs(self):
        from jolymer.sas.CoupledMeasurement import ms_from_folder
        path = self.path_edit.text()
        name = self.name_edit.text()
        q_beamstop = self.q_beamstop_edit.value()
        start, end = self.start_idx.value(), self.end_idx.value()
        OUV = self.state.uv
        sample = self.state.sample
        try:
            ms_saxs = ms_from_folder(path, file_prefix=name, max_seqi=end,
                                     min_seqi=start,
                                     q_beamstop=q_beamstop)
            cm_saxs = CoupledMeasurement(ms_saxs, OUV, sample)
            self.state.saxs = ms_saxs
            self.state.cm_saxs = cm_saxs
            self.update_plot()
        except Exception as e:
            if hasattr(self.state, "console"):
                self.state.console._append_plain_text(f"Error loading SAXS: {e}")
            else:
                print(f"Error loading SAXS: {e}")

    def update_plot(self):
        if self.state.cm_saxs is None:
            return
        try:
            saxs_kind = self.saxs_kind_edit.text()
            qstar = self.qstar_edit.value()
            qmin = self.qmin_edit.value()
            qmax = self.qmax_edit.value()
            df = self.state.cm_saxs.get_saxs_scalar(kind=saxs_kind,
                                                    qstar=qstar,
                                                    qmin=qmin,
                                                    qmax=qmax)
            self.ax.clear()
            self.ax.plot(df.time, df.I0, label=f"SAXS {saxs_kind}")
            self.ax.set_xlabel("Time [s]")
            self.ax.set_ylabel("I(t)")
            self.ax.legend()
            self.canvas.draw()
        except Exception as e:
            if hasattr(self.state, "console"):
                self.state.console._append_plain_text(f"Error plotting SAXS: {e}")
            else:
                print(f"Error plotting SAXS: {e}")

    def prepare_coupled(self):
        self.load_saxs()
