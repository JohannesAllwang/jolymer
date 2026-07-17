from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QPushButton,
    QFileDialog, QFormLayout, QDoubleSpinBox,
    QLineEdit, QLabel, QSpinBox, QHBoxLayout,
    QBoxLayout
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from jolymer.uv.onlineUV import onlineUV
from jolymer.sas.CoupledMeasurement import CoupledMeasurement, ms_from_folder
from jolymer.gui.regals.state import BioREGALSState
from jolymer.gui.regals.console.ipython_widget import IPythonConsole

from pathlib import Path

class WAXSWindow(QMainWindow):
    def __init__(self, state: BioREGALSState, parent=None):
        super().__init__(parent)
        self.state = state
        self.setWindowTitle("WAXS")
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        form = QFormLayout()
        # --- Folder and load button ---
        self.path_edit = QLineEdit()
        self.name_edit = QLineEdit()
        self.load_btn = QPushButton("Load WAXS frames")
        self.load_btn.clicked.connect(self.load_waxs)
        layout.addWidget(QLabel("WAXS folder:"))
        prepath = "/path/to/datafolder/"
        if not self.state.waxs is None:
            prepath = str(self.state.waxs.ms[0].path)
        self.path_edit.setText(prepath)
        layout.addWidget(self.path_edit)
        prename = "filename_prefix"
        if not self.state.waxs is None:
            prename = str(self.state.waxs.name)
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
        self.waxs_kind_edit = QLineEdit("Iq")
        self.qstar_edit = QDoubleSpinBox()
        self.qstar_edit.setDecimals(3)
        self.qstar_edit.setSingleStep(0.01)
        self.qmin_edit = QDoubleSpinBox()
        self.qmin_edit.setDecimals(3)
        self.qmin_edit.setSingleStep(0.01)
        self.qmax_edit = QDoubleSpinBox()
        self.qmax_edit.setDecimals(3)
        self.qmax_edit.setSingleStep(0.01)

        if self.state.cm_waxs is None:
            self.qmax_edit.setValue(0.4)
            self.qmin_edit.setValue(0.04)
            self.qstar_edit.setValue(0.9)
        else:
            cm_waxs = self.state.cm_waxs
            self.qmax_edit.setValue(cm_waxs.qmax)
            self.qmin_edit.setValue(cm_waxs.qmin)
            self.qstar_edit.setValue(cm_waxs.qstar)
        self.start_idx.setValue(state.saxs.min_seqi)
        self.end_idx.setValue(state.saxs.max_seqi)


        form_layout = QFormLayout()
        form_layout.addRow("WAXS kind:", self.waxs_kind_edit)
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
    def load_waxs(self):
        from jolymer.sas.CoupledMeasurement import ms_from_folder
        path = self.path_edit.text()
        name = self.name_edit.text()
        start, end = self.start_idx.value(), self.end_idx.value()
        OUV = self.state.uv
        sample = self.state.sample
        cm_saxs = self.state.cm_saxs
        if cm_saxs is None:
            self.state.console._append_plain_text("Load SAXS first")
            return
        try:
            ms_waxs = ms_from_folder(path, file_prefix=name, max_seqi=end,
                                     min_seqi=start)
            cm_waxs = CoupledMeasurement(ms_waxs, OUV, sample)
            cm_waxs._alignment = cm_saxs._alignment.copy()
            cm_waxs._alignment.update({
                "saxs_kind": self.waxs_kind_edit.text(),
                "qstar": self.qstar_edit.value(),
                "qmin": self.qmin_edit.value(),
                "qmax": self.qmax_edit.value(),
            })
            cm_waxs._uv_on_saxs = cm_saxs._uv_on_saxs
            to_regals_waxs = cm_waxs.to_regals()
            q_waxs = to_regals_waxs["q"]
            I_waxs = to_regals_waxs["I"]
            sigma_waxs = to_regals_waxs["sigma"]
            x_waxs = to_regals_waxs["x"]
            uv_meas_waxs = to_regals_waxs["uv_meas"]

            self.state.waxs = ms_waxs
            self.state.to_regals_waxs = to_regals_waxs
            self.state.cm_waxs = cm_waxs
            self.update_plot()
        except Exception as e:
            if hasattr(self.state, "console"):
                self.state.console._append_plain_text(f"Error loading WAXS: {e}")
            else:
                print(f"Error loading WAXS: {e}")

    def update_plot(self):
        if self.state.cm_waxs is None:
            return
        try:
            self.ax.clear()
            self.state.cm_waxs.plot_alignment(ax=self.ax)
            self.ax.legend()
            self.canvas.draw()
        except Exception as e:
            if hasattr(self.state, "console"):
                self.state.console._append_plain_text(f"Error plotting WAXS: {e}")
            else:
                print(f"Error plotting WAXS: {e}")

    def prepare_coupled(self):
        self.load_waxs()
