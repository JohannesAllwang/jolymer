from PyQt6.QtWidgets import (
    QDockWidget, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QSpinBox, QDoubleSpinBox,
    QFileDialog, QGroupBox
)
from PyQt6.QtCore import pyqtSignal

from jolymer.samples.bioMOLECULE import bioMOLECULE
from jolymer.sas.bioREGALS import bioComponent, concentration_class, profile_class, bioMIXTURE

from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QSpinBox, QDoubleSpinBox, QCheckBox, QWidget, QScrollArea
)

class ComponentRow(QWidget):
    """One row of widgets for a single bioComponent."""
    def __init__(self, parent=None,
                 xmin0=0, xmax0=400):
        super().__init__(parent)
        layout = QHBoxLayout(self)

        self.xmin = QDoubleSpinBox()
        self.xmin.setRange(0, 10000)
        self.xmin.setValue(xmin0)
        self.xmax = QDoubleSpinBox()
        self.xmax.setRange(0, 10000)
        self.xmax.setValue(xmax0)
        self.zero_min = QCheckBox("Zero @ min")
        self.zero_max = QCheckBox("Zero @ max")
        self.Nw = QSpinBox()
        self.Nw.setRange(1, 101)
        self.Nw.setValue(31)
        self.dmax = QDoubleSpinBox()
        self.dmax.setRange(0, 1000)
        self.dmax.setValue(5)
        self.remove_btn = QPushButton("Remove")

        for w in [self.xmin, self.xmax, self.zero_min, self.zero_max, self.Nw, self.dmax, self.remove_btn]:
            layout.addWidget(w)

class ComponentsEditorWindow(QDialog):
    def __init__(self, state, parent=None):
        super().__init__(parent)
        self.state = state
        to_regals = self.state.to_regals
        if not to_regals is None:
            self.x = to_regals['x']
            self.q = to_regals['q']
            self.component_widgets = []
            main_layout = QVBoxLayout(self)
            # Scrollable area for components
            self.scroll = QScrollArea()
            self.scroll.setWidgetResizable(True)
            self.scroll_content = QWidget()
            self.scroll_layout = QVBoxLayout(self.scroll_content)
            self.scroll.setWidget(self.scroll_content)
            main_layout.addWidget(self.scroll)
            # Buttons
            btn_layout = QHBoxLayout()
            self.add_btn = QPushButton("Add Component")
            self.generate_btn = QPushButton("Generate bioMIXTURE")
            btn_layout.addWidget(self.add_btn)
            btn_layout.addWidget(self.generate_btn)
            main_layout.addLayout(btn_layout)
            # Signals
            self.add_btn.clicked.connect(self.add_component)
            self.generate_btn.clicked.connect(self.generate_mixture)
            # Add initial component
            self.add_component()

    def add_component(self):
        row = ComponentRow(xmin0=self.x.min(), xmax0=self.x.max())
        row.remove_btn.clicked.connect(lambda: self.remove_component(row))
        self.scroll_layout.addWidget(row)
        self.component_widgets.append(row)

    def remove_component(self, row):
        self.scroll_layout.removeWidget(row)
        row.setParent(None)
        self.component_widgets.remove(row)

    def generate_mixture(self):
        uv_meas = self.state.to_regals['uv_meas']
        peaks = []
        for row in self.component_widgets:
            conc = concentration_class(
                'smooth', self.x,
                xmin=row.xmin.value(),
                xmax=row.xmax.value(),
                is_zero_at_xmin=row.zero_min.isChecked(),
                is_zero_at_xmax=row.zero_max.isChecked(),
                Nw=row.Nw.value()
            )
            prof = profile_class('realspace', self.q, dmax=row.dmax.value(), Nw=101)
            peaks.append(bioComponent(conc, prof))
        # Store in state
        self.state.mixture = bioMIXTURE(peaks, uv_meas=uv_meas)
        print(f"Generated bioMIXTURE with {len(peaks)} components")

class ComponentsDock(QDockWidget):
    mixture_changed = pyqtSignal()
    def __init__(self, state, parent=None):
        super().__init__("REGALS Components", parent)
        self.state = state
        self.container = QWidget()
        self.layout = QVBoxLayout(self.container)
        self.setWidget(self.container)
        self.refresh()
        add_btn = QPushButton("+ Add component")
        # add_btn.clicked.connect(self.add_component)
        self.layout.addWidget(add_btn)

    # -------------------------
    def refresh(self):
        # Clear existing widgets
        while self.layout.count():
            item = self.layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        if self.state.mixture is None:
            return
        for i, comp in enumerate(self.state.mixture.peaks):
            self.layout.addWidget(
                self._component_editor(i, comp)
            )

    # -------------------------
    def _component_editor(self, idx, comp):
        box = QGroupBox(f"Component {idx+1}")
        layout = QVBoxLayout(box)
        # ---- PDB / bioMOLECULE ----
        pdb_row = QHBoxLayout()
        pdb_label = QLabel(
            getattr(comp, "biomolecule", None).name
            if hasattr(comp, "biomolecule")
            else "No PDB loaded"
        )
        load_btn = QPushButton("Load PDB")
        load_btn.clicked.connect(
            lambda: self.load_pdb(comp, pdb_label)
        )
        pdb_row.addWidget(pdb_label)
        pdb_row.addWidget(load_btn)
        layout.addLayout(pdb_row)
        # ---- dmax ----
        dmax = QDoubleSpinBox()
        dmax.setRange(1, 500)
        dmax.setValue(comp.profile.dmax)
        dmax.valueChanged.connect(
            lambda v: setattr(comp.profile, "dmax", v)
        )
        layout.addWidget(QLabel("dmax"))
        layout.addWidget(dmax)
        # ---- x-range ----
        xmin = QSpinBox()
        xmax = QSpinBox()
        xmin.setRange(0, 10000)
        xmax.setRange(0, 10000)
        xmin.setValue(comp.concentration.xmin)
        xmax.setValue(comp.concentration.xmax)
        xmin.valueChanged.connect(
            lambda v: setattr(comp.concentration, "xmin", v)
        )
        xmax.valueChanged.connect(
            lambda v: setattr(comp.concentration, "xmax", v)
        )
        layout.addWidget(QLabel("x-range"))
        layout.addWidget(xmin)
        layout.addWidget(xmax)
        return box

    # -------------------------
    def load_pdb(self, comp, label):
        fname, _ = QFileDialog.getOpenFileName(
            self, "Load PDB", filter="PDB files (*.pdb)"
        )
        if not fname:
            return
        bm = bioMOLECULE.from_pdb(fname)
        comp.biomolecule = bm
        label.setText(bm.name)
        # ---- autosuggest ----
        if hasattr(bm, "dmax"):
            comp.profile.dmax = bm.dmax
        if hasattr(bm, "uv_scale"):
            comp.uv_scale = bm.uv_scale
        self.mixture_changed.emit()

