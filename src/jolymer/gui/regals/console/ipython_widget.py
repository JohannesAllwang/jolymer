# gui/regals/console/ipython_widget.py

from qtconsole.rich_jupyter_widget import RichJupyterWidget
from qtconsole.inprocess import QtInProcessKernelManager


class IPythonConsole(RichJupyterWidget):
    def __init__(self, namespace=None, startup_file=None, parent=None):
        super().__init__(parent)

        # In-process kernel (Qt-safe)
        self.kernel_manager = QtInProcessKernelManager()
        self.kernel_manager.start_kernel(show_banner=False)

        self.kernel = self.kernel_manager.kernel
        self.kernel.gui = "qt"

        self.kernel_client = self.kernel_manager.client()
        self.kernel_client.start_channels()

        # Inject namespace
        if namespace:
            self.kernel.shell.push(namespace)

        # Run startup file (like Spyder)
        if startup_file:
            try:
                with open(startup_file) as f:
                    code = f.read()
                self.kernel.shell.run_cell(code)
            except Exception as e:
                print(f"Failed to load startup file: {e}")

