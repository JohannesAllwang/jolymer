import numpy as np
import subprocess
import getpass
import io
import os
from os.path import join
import pandas as pd

sasview_bin = f"C:\\Users\\{getpass.getuser()}\\sasview"


def run_program(name, *args):
    cmd = ['powershell', f'./{name}.exe', *args]
    return subprocess.run(cmd, cwd=sasview_bin)

def run_sasview(*args):
    out = run_program('sasview', *args)
    return out
