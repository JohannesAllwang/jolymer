import originpro as op
import pandas as pd
import numpy as np

import sys

def origin_shutdown_exception_hook(exctype, value, traceback):
    op.exit()
    sys.__exceptionhook__(exctype, value, traceback)

if op and op.next:
    sys.excthook = origin_shutdown_exception_hook

if op.oext:
    op.set_show(True)


