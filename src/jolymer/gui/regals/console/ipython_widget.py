# gui/regals/console/ipython_widget.py

# from qtconsole.rich_jupyter_widget import RichJupyterWidget
# from qtconsole.inprocess import QtInProcessKernelManager


# class IPythonConsole(RichJupyterWidget):
#     def __init__(self, namespace=None, startup_file=None, parent=None):
#         super().__init__(parent)
#         self.kernel_manager = QtInProcessKernelManager()
#         self.kernel_manager.start_kernel(show_banner=False)
#         self.kernel = self.kernel_manager.kernel
#         self.kernel.gui = "qt"
#         self.kernel_client = self.kernel_manager.client()
#         self.kernel_client.start_channels()
#         if namespace:
#             self.kernel.shell.push(namespace)
#         if startup_file:
#             try:
#                 with open(startup_file) as f:
#                     code = f.read()
#                 self.kernel.shell.run_cell(code)
#             except Exception as e:
#                 print(f"Failed to load startup file: {e}")

from pyqtgraph.console import ConsoleWidget

class SimpleConsole(ConsoleWidget):
    def __init__(self, namespace=None, startup_file=None, parent=None):
        super().__init__(
            parent=parent,
            namespace=namespace or {},
        )
        if startup_file:
            try:
                with open(startup_file) as f:
                    self.execCommand(f.read())
            except Exception as e:
                print(f"Failed to load startup file: {e}")

    def _append_plain_text(self, text):
        """Compatibility with the old IPythonConsole."""
        self.output.insertPlainText(str(text))
        self.output.ensureCursorVisible()
