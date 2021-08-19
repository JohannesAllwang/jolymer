from .desy import Desy
from .. import database_operations as dbo
from .. import os_utility as osu
from . import sas_plotlib as sp

import matplotlib.pyplot as plt
import os
import pandas as pd

def get_m(tid, T):
    query = f"""
    SELECT * FROM desy_measurements
    WHERE sample = 'tresy_{tid}'
    AND comment = '{T} deg'
    """
    did = dbo.execute(query)[0][0]
    m = Desy(did)
    m.targetT = T
    return m

def from_query(query, T):
    with dbo.dbopen() as c:
        conn = c.connection
        df = pd.read_sql(query, conn)
    dids = list(df.id)
    ms = [get_m(did, T) for did in dids]
    return ms

# model = ?

def plot_tresyfits(tresy_numbers, model, p0, bounds, fixed_pars={}):
    meas_nums= tresy_numbers
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (12,9), 
                             sharex=True, sharey=True)
    for i, ax in zip(tresy_numbers, axes.flatten()):
        index = i-1
        # iqmin=iqmins[i-1]
        # color = cm[i-1]
        m = get_m(i, 20)
        color = 'blue'
        iqmin=0
        fit_dict, fit_df = model.fit(m, bounds=bounds, 
                            p0=p0, iqmin=iqmin, fixed_parameters=fixed_pars)
        label = f'{m.sample.get_NPname()} pH {m.sample.buffer.pH} {m.sample.buffer.salt_concentration} salt'
        marker = '.'
        df = m.get_data(cout=False)[iqmin::]
        ax.errorbar(df.q, df.I, df.err_I, marker = marker, color=color,
                        linestyle='', label = label, elinewidth=0.2)
        model.plot_fit(fit_df, (fig, ax), color='salmon')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        #$\\pm$ {3:.2f}
        text=model.get_text(fit_dict)
        ax.annotate(text, xy=(0.0, 0.0), 
                    xycoords='axes fraction')
        ax.set_ylim(1e-7, 1)
        ax.grid()

    axes[1][0].set_xlabel('$q$ [1/nm]')
    axes[1][1].set_xlabel('$q$ [1/nm]')
    axes[1][0].set_ylabel('$I$ [1/cm]')
    axes[0][0].set_ylabel('$I$ [1/cm]')
