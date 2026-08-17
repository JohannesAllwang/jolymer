import sys
from pathlib import Path

import numpy as np

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QFormLayout,
    QPushButton, QLineEdit, QFileDialog, QSplitter,
    QDoubleSpinBox, QCheckBox, QStyle, QSizePolicy, QLabel,
    QToolButton, QScrollArea, QComboBox
)

# --- matplotlib ---
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.ticker

# --- IPython console (Spyder-style) ---
from jolymer.gui.regals.console.ipython_widget import IPythonConsole

from jolymer.sas.CoupledMeasurement import CoupledMeasurement, ms_from_folder
from jolymer.uv.onlineUV import onlineUV
from jolymer.sas.SEC_SWAXS import SEC_SWAXS

import traceback

# -----------------------------
# Small reusable "path + browse button" widget
# -----------------------------
class PathRow(QWidget):
    def __init__(self, mode="file", parent=None):
        """mode: 'file' or 'dir'"""
        super().__init__(parent)
        self.mode = mode
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.edit = QLineEdit()
        self.edit.setPlaceholderText("Select path…")
        btn = QPushButton()
        btn.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon))
        btn.setFixedWidth(30)
        btn.clicked.connect(self.browse)
        layout.addWidget(self.edit)
        layout.addWidget(btn)

    def browse(self):
        if self.mode == "file":
            path, _ = QFileDialog.getOpenFileName(self, "Select file", "", "All files (*)")
        else:
            path = QFileDialog.getExistingDirectory(self, "Select directory")
        if path:
            self.edit.setText(path)

    def text(self):
        return self.edit.text()

    def setText(self, text):
        self.edit.setText(text)


# -----------------------------
# Collapsible panel (click the header to fold/unfold)
# -----------------------------
class CollapsibleBox(QWidget):
    def __init__(self, title="", expanded=True, parent=None):
        super().__init__(parent)
        self.toggle_button = QToolButton()
        self.toggle_button.setText(title)
        self.toggle_button.setCheckable(True)
        self.toggle_button.setChecked(expanded)
        self.toggle_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.toggle_button.setArrowType(
            Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow
        )
        self.toggle_button.setStyleSheet(
            "QToolButton { border: none; font-weight: 600; padding: 4px; }"
        )
        self.toggle_button.clicked.connect(self._on_toggle)

        self.content = QWidget()
        self.form = QFormLayout(self.content)
        self.content.setVisible(expanded)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addWidget(self.toggle_button)
        outer.addWidget(self.content)

        line = QWidget()
        line.setFixedHeight(1)
        line.setStyleSheet("background-color: palette(mid);")
        outer.addWidget(line)

    def _on_toggle(self):
        expanded = self.toggle_button.isChecked()
        self.toggle_button.setArrowType(
            Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow
        )
        self.content.setVisible(expanded)


# -----------------------------
# Matplotlib canvas
# -----------------------------
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(7, 4))
        super().__init__(self.fig)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)


class SecMainWindow(QMainWindow):
    def __init__(self, state: SEC_SWAXS, parent=None):
        super().__init__(parent)
        self.state = state
        self.setWindowTitle("SEC-SAXS / UV overview")
        self.resize(1400, 900)

        central = QWidget()
        self.setCentralWidget(central)
        root_layout = QHBoxLayout(central)
        root_layout.setContentsMargins(0, 0, 0, 0)

        # ---------------------------------------------------------
        # LEFT: collapsible parameter panels (UV, SAXS, plot options),
        # wrapped in a scroll area so it never forces the window taller.
        # ---------------------------------------------------------
        left_content = QWidget()
        left_layout = QVBoxLayout(left_content)
        left_layout.setSpacing(0)
        left_layout.addWidget(self._build_uv_group())
        left_layout.addWidget(self._build_saxs_group())
        left_layout.addWidget(self._build_alignment_group())
        left_layout.addWidget(self._build_plot_group())
        left_layout.addStretch()

        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setWidget(left_content)
        left_scroll.setMaximumWidth(360)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        # ---------------------------------------------------------
        # RIGHT: plot on top, console below (resizable split), console kept
        # ---------------------------------------------------------
        right_splitter = QSplitter(Qt.Orientation.Vertical)
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        self.canvas = MplCanvas()
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        right_splitter.addWidget(self.toolbar)
        right_splitter.addWidget(self.canvas)

        self.console = IPythonConsole(namespace={
            "state": self.state,
            "cm": self.state.CM,
            "np": np})
        right_splitter.addWidget(self.console)
        right_splitter.setStretchFactor(0, 3)
        right_splitter.setStretchFactor(1, 1)

        main_splitter = QSplitter(Qt.Orientation.Horizontal)
        main_splitter.addWidget(left_scroll)
        main_splitter.addWidget(right_splitter)
        main_splitter.setStretchFactor(0, 0)
        main_splitter.setStretchFactor(1, 1)
        root_layout.addWidget(main_splitter)

        self._hydrate_from_state()

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
            self.state.load_from_json(Path(path))
            # self.state.rebuild_from_state()
        except Exception as e:
            QMessageBox.critical(self, "Load failed", str(e))

    # ---------------------------------------------------------------
    # Panel builders
    # ---------------------------------------------------------------
    def _build_uv_group(self):
        box = CollapsibleBox("UV", expanded=True)
        form = box.form

        self.uv_path = PathRow(mode="file")
        self.uv_refwl = QDoubleSpinBox()
        self.uv_refwl.setRange(0, 2000)
        self.uv_refwl.setValue(280)
        self.uv_outwl_from = QDoubleSpinBox()
        self.uv_outwl_from.setRange(0, 2000)
        self.uv_outwl_from.setValue(280)
        self.uv_outwl_to = QDoubleSpinBox()
        self.uv_outwl_to.setRange(0, 2000)
        self.uv_outwl_to.setValue(280)
        self.uv_alignment_time = QDoubleSpinBox()
        self.uv_alignment_time.setMaximum(1e6)
        self.uv_alignment_time.setValue(1600)
        self.uv_mintime = QDoubleSpinBox()
        self.uv_mintime.setMaximum(1e9)
        self.uv_maxtime = QDoubleSpinBox()
        self.uv_maxtime.setMaximum(1e9)
        self.uv_maxtime.setValue(1e6)

        form.addRow("File", self.uv_path)
        form.addRow("refwl", self.uv_refwl)
        form.addRow("outwl from", self.uv_outwl_from)
        form.addRow("outwl to", self.uv_outwl_to)
        form.addRow("alignment_time", self.uv_alignment_time)
        form.addRow("mintime", self.uv_mintime)
        form.addRow("maxtime", self.uv_maxtime)

        load_btn = QPushButton("Load UV (.spc)")
        load_btn.clicked.connect(self.load_uv)
        form.addRow(load_btn)
        return box

    def _build_saxs_group(self):
        box = CollapsibleBox("SAXS", expanded=False)
        form = box.form

        self.saxs_dir = PathRow(mode="dir")
        self.saxs_pattern = QLineEdit("B3")
        self.autorg_name = QLineEdit("autorg.dat")

        form.addRow("Folder", self.saxs_dir)
        form.addRow("File pattern", self.saxs_pattern)
        form.addRow("autorg output", self.autorg_name)

        run_btn_saxs = QPushButton("Load SAXS")
        run_btn_saxs.clicked.connect(self.run_load_saxs)
        form.addRow(run_btn_saxs)
        run_btn = QPushButton("Run autorg")
        run_btn.clicked.connect(self.run_autorg)
        form.addRow(run_btn)
        return box

    def _build_alignment_group(self):
        box = CollapsibleBox("Alignment (UV → SAXS)", expanded=False)
        form = box.form
        # ------------------------------------------------------------
        # Alignment method
        # ------------------------------------------------------------
        self.align_method = QComboBox()
        self.align_method.addItems([
            "From Para",
            "From H5",
            "Stochastic fit",
        ])
        self.align_method.setCurrentText("From Para")
        self.align_method.currentTextChanged.connect(
            self._update_alignment_options
        )
        form.addRow("Method", self.align_method)
        # ------------------------------------------------------------
        # Stochastic-fit options
        # ------------------------------------------------------------
        self.align_para_widget = QWidget()
        para_form = QFormLayout(self.align_para_widget)
        self.refine_para = QCheckBox(
            "Refine Para by fitting"
        )
        para_form.addRow(self.refine_para)
        form.addRow(self.align_para_widget)
        self.align_options_widget = QWidget()
        options_form = QFormLayout(self.align_options_widget)
        self.align_min_t = QDoubleSpinBox()
        self.align_min_t.setMaximum(1e9)
        self.align_min_t.setValue(400)
        self.align_max_t = QDoubleSpinBox()
        self.align_max_t.setMaximum(1e9)
        self.align_max_t.setValue(1e6)
        self.align_i0_min = QDoubleSpinBox()
        self.align_i0_min.setDecimals(6)
        self.align_i0_min.setMaximum(1e9)
        self.align_scale0 = QDoubleSpinBox()
        self.align_scale0.setRange(-100, 100)
        self.align_scale0.setDecimals(4)
        self.align_scale0.setValue(1.0)
        self.align_shift0 = QDoubleSpinBox()
        self.align_shift0.setRange(-1e6, 1e6)
        self.align_shift0.setValue(0.0)
        self.auto_apply_shift = QCheckBox(
            "apply fitted shift to UV alignment_time"
        )
        self.auto_apply_shift.setChecked(True)
        options_form.addRow("min t (SAXS, s)", self.align_min_t)
        options_form.addRow("max t (SAXS, s)", self.align_max_t)
        options_form.addRow("I0 min", self.align_i0_min)
        options_form.addRow("scale0", self.align_scale0)
        options_form.addRow("shift0", self.align_shift0)
        options_form.addRow(self.auto_apply_shift)
        form.addRow(self.align_options_widget)
        # ------------------------------------------------------------
        # Result
        # ------------------------------------------------------------
        self.align_result_label = QLabel("not aligned yet")
        self.align_result_label.setWordWrap(True)
        form.addRow(self.align_result_label)
        self.align_btn = QPushButton("Align UV → SAXS")
        self.align_btn.clicked.connect(self.run_alignment)
        form.addRow(self.align_btn)
        self._update_alignment_options(self.align_method.currentText())
        return box

    def _update_alignment_options(self, method):
        is_stochastic = method == "Stochastic fit"
        self.align_options_widget.setVisible(is_stochastic)
        is_para = method == "From Para"
        self.align_para_widget.setVisible(is_para)

    def _build_plot_group(self):
        box = CollapsibleBox("Overview plot", expanded=False)
        form = box.form

        self.wl_min = QDoubleSpinBox()
        self.wl_min.setRange(0, 2000)
        self.wl_min.setValue(260)
        self.wl_max = QDoubleSpinBox()
        self.wl_max.setRange(0, 2000)
        self.wl_max.setValue(290)
        self.plot_refwl = QDoubleSpinBox()
        self.plot_refwl.setRange(0, 2000)
        self.plot_refwl.setValue(280)
        self.xmin = QDoubleSpinBox()
        self.xmin.setMaximum(500)
        self.xmax = QDoubleSpinBox()
        self.xmax.setMaximum(1e9)
        self.xmax.setValue(3000)
        self.minutes_check = QCheckBox("x-axis in minutes")
        self.sample_name = QLineEdit()

        form.addRow("UV wl min", self.wl_min)
        form.addRow("UV wl max", self.wl_max)
        form.addRow("UV refwl", self.plot_refwl)
        form.addRow("x min (s)", self.xmin)
        form.addRow("x max (s)", self.xmax)
        form.addRow(self.minutes_check)
        form.addRow("label", self.sample_name)

        plot_btn = QPushButton("Update plot")
        plot_btn.clicked.connect(self.update_plot)
        form.addRow(plot_btn)
        return box

    def _hydrate_from_state(self):
        """Pre-fill fields from an already-populated CoupledMeasurement (if any)."""
        uv = self.state.CM.uv
        if uv is not None:
            if getattr(uv, "refwl", None) is not None:
                self.uv_refwl.setValue(uv.refwl)
            outwl = getattr(uv, "outwl", None)
            if outwl is not None:
                if isinstance(outwl, (list, tuple, np.ndarray)):
                    self.uv_outwl_from.setValue(outwl[0])
                    self.uv_outwl_to.setValue(outwl[-1])
                else:
                    self.uv_outwl_from.setValue(outwl)
                    self.uv_outwl_to.setValue(outwl)
            if getattr(uv, "alignment_time", None) is not None:
                self.uv_alignment_time.setValue(uv.alignment_time)
            if getattr(uv, "mintime", None) is not None:
                self.uv_mintime.setValue(uv.mintime)
            if getattr(uv, "maxtime", None) is not None:
                self.uv_maxtime.setValue(uv.maxtime)
            if getattr(uv, "path", None) and getattr(uv, "spec_filename", None):
                self.uv_path.setText(str(Path(uv.path, uv.spec_filename)))

        try:
            saxs_dir = self.state.CM.get_saxs_path()
            if saxs_dir:
                self.saxs_dir.setText(str(saxs_dir))
        except Exception:
            pass

    # ---------------------------------------------------------------
    # Logging helper (writes into the console instead of popping dialogs)
    # ---------------------------------------------------------------
    def _log(self, msg):
        if getattr(self, "console", None) is not None:
            self.console._append_plain_text(str(msg) + "\n")
        else:
            print(msg)

    def _log_exception(self, prefix, exc):
        """Log an error message plus the full traceback (to the console
        and stderr), instead of just str(exc), so real bugs aren't
        reduced to an uninformative one-liner."""
        tb = traceback.format_exc()
        self._log(f"{prefix}: {exc}\n{tb}")
        print(tb, file=sys.stderr)


    # ---------------------------------------------------------------
    # Actions
    # ---------------------------------------------------------------
    def load_uv(self):
        try:
            path = Path(self.uv_path.text())
            if not path.exists():
                raise ValueError("Invalid UV file path")
            wl_from, wl_to = self.uv_outwl_from.value(), self.uv_outwl_to.value()
            if np.isclose(wl_from, wl_to):
                outwl = wl_from
            else:
                outwl = list(range(int(round(wl_from)), int(round(wl_to)) + 1))
            uv = onlineUV(
                spec_filename=path.name,
                path=str(path.parent),
                refwl=self.uv_refwl.value(),
                outwl=outwl,
                alignment_time=float(self.uv_alignment_time.value()),
                mintime=self.uv_mintime.value(),
                maxtime=self.uv_maxtime.value(),
            )
            uv.get_data()
            self.state.CM.uv = uv
            self._log(f"Loaded UV: {path}")
        except Exception as e:
            self._log_exception(f"Error loading UV:", e)

    def run_autorg(self):
        try:
            directory = Path(self.saxs_dir.text())
            if not directory.exists():
                raise ValueError("Invalid SAXS directory")
            out = self.state.CM.run_autorg(
                autorg_name=self.autorg_name.text(),
                pattern=f"{self.saxs_pattern.text()}_*.dat",
                directory=directory,
            )
            self._log(f"autorg written to {out}")
        except Exception as e:
            self._log_exception(f"Error running autorg:", e)

    def run_load_saxs(self):
        try:
            directory = Path(self.saxs_dir.text())
            if not directory.exists():
                raise ValueError("Invalid SAXS directory")
            ms = ms_from_folder(path=directory, file_prefix=self.saxs_pattern.text(),
                   min_seqi=0, max_seqi=100000, q_beamstop=0.006,
                   exclude_seqi=[], frame_time=2.1,
                   angular_unit='A')
            self.state.CM.saxs_list = ms
            self._log(f"saxs loaded to")
        except Exception as e:
            self._log_exception(f"Error running autorg:", e)

    def run_alignment(self):
        method = self.align_method.currentText()
        if method == "From H5":
            self._align_from_h5()
        elif method == "From Para":
            self._align_from_para()
        elif method == "Stochastic fit":
            self._align_stochastic()
        self.update_plot()

    def _align_from_para(self):
        try:
            if self.state.CM.uv is None:
                raise ValueError("Load UV data first")
            if self.state.CM.saxs_list is None:
                raise ValueError("Load SAXS data first")
            self.state.align_from_Para()
            if self.refine_para.isChecked():
                self.state.CM.refine_uv_alignment(
                    shift_bounds=(-100, 100)
                )
                for m in self.state.CM.saxs_list:
                    m.time = m.time - self.state.CM._alignment['shift']
        except Exception as e:
            self._log_exception(f"Error aligning UV/SAXS:", e)

    def _align_from_h5(self):
        try:
            if self.state.CM.uv is None:
                raise ValueError("Load UV data first")
            if self.state.CM.saxs_list is None:
                raise ValueError("Load SAXS data first")
            self.state.align_from_h5()
        except Exception as e:
            self._log_exception(f"Error aligning UV/SAXS:", e)

    def _align_stochastic(self):
        try:
            if self.state.CM.uv is None:
                raise ValueError("Load UV data first")
            directory = Path(self.saxs_dir.text())
            autorg_path = directory / self.autorg_name.text()
            if not autorg_path.exists():
                self.run_autorg()
            # make sure self.state.CM._autorg_df reflects the chosen directory/name
            self.state.CM.load_autorg(
                autorg_filename=self.autorg_name.text(),
                directory=directory,
            )
            alignment = self.state.CM.align_uv_to_saxs(
                saxs_kind="autorg",
                min_t_saxs=self.align_min_t.value(),
                max_t_saxs=self.align_max_t.value(),
                I0_min=self.align_i0_min.value(),
                scale0shift0=[self.align_scale0.value(), self.align_shift0.value()],
            )
            self.align_result_label.setText(
                f"scale={alignment['scale']:.4f}  shift={alignment['shift']:.1f}  "
                f"score={alignment['score']:.3f}  success={alignment['success']}"
            )
            self._log(f"Alignment: {alignment}")
            if self.auto_apply_shift.isChecked():
                for m in self.state.CM.saxs_list:
                    m.time = m.time / alignment['scale'] - alignment['shift'] / alignment['scale']
        except Exception as e:
            self._log_exception(f"Error aligning UV/SAXS:", e)

    def update_plot(self):
        try:
            if self.state.CM.uv is None:
                raise ValueError("Load UV data first")
            directory = Path(self.saxs_dir.text())
            autorg_path = directory / self.autorg_name.text()
            if not autorg_path.exists():
                self.run_autorg()
            autorg_df = self.state.CM.load_autorg(
                autorg_filename=self.autorg_name.text(),
                directory=directory,
            )
            fig = self.canvas.fig
            fig.clear()
            ax = fig.add_subplot(111)
            axU = ax.twinx()
            self.state.CM.uv.plot_full_wavelength_map(
                wl_min=self.wl_min.value(),
                wl_max=self.wl_max.value(),
                refwl=self.plot_refwl.value(),
                ax=axU,
                alignment_time=self.uv_alignment_time.value(),
            )
            ax.plot(autorg_df.time, autorg_df.Rg, color="green", label="Rg")
            print(autorg_df)
            ax.set_ylabel(r"$R_g$ (nm)", color="green")
            ax.tick_params(axis="y", colors="green")
            axI = ax.twinx()
            axI.spines["right"].set_position(("outward", 50))
            axI.plot(autorg_df.time, autorg_df.I0, color="black", label="I(0)")
            axI.set_ylabel(r"$I(0)$", color="black")
            axI.tick_params(axis="y", colors="black")
            ax.set_xlabel("Time (s)")
            xmin, xmax = self.xmin.value(), self.xmax.value()
            if xmax > xmin:
                ax.set_xlim(xmin, xmax)
            saxs = self.state.CM.saxs_list
            saxs_time = np.array([m.time for m in saxs], dtype=float)
            saxs_frame = np.array([m.seqi for m in saxs], dtype=float)
            ax_frame = ax.twiny()
            def update_frame_axis(ax):
                xmin, xmax = ax.get_xlim()
                mask = (
                    (saxs_time >= xmin) &
                    (saxs_time <= xmax)
                )
                times = saxs_time[mask]
                frames = saxs_frame[mask]
                if len(times) == 0:
                    ax_frame.set_xticks([])
                    return
                target_labels = 8
                step = max(1, int(np.ceil(len(frames) / target_labels)))
                ax_frame.set_xticks(times[::step])
                ax_frame.set_xticklabels(
                    [str(int(f)) for f in frames[::step]]
                )
            def on_xlim_changed(ax):
                ax_frame.set_xlim(ax.get_xlim())
                update_frame_axis(ax)
            ax_frame.set_xlim(ax.get_xlim())
            update_frame_axis(ax)
            ax.callbacks.connect("xlim_changed", on_xlim_changed)
            ax_frame.set_xlabel("X-ray frame")
            ax.set_ylim(0,150)
            name = self.sample_name.text()
            if name:
                ax.annotate(name, xy=(0.5, 0.9), xycoords="axes fraction")
            fig.tight_layout()
            self.canvas.draw_idle()
            self._log("Plot updated.")
        except Exception as e:
            self._log_exception(f"Error updating plot:", e)


def main():
    app = QApplication(sys.argv)
    CM = CoupledMeasurement([], None, None)
    state = SEC_SWAXS(CM=CM)
    win = SecMainWindow(state=state)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
