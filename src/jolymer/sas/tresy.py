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
    return Desy(did)

def from_query(query, T):
    with dbo.dbopen() as c:
        conn = c.connection
        df = pd.read_sql(query, conn)
    dids = list(df.id)
    ms = [get_m(did, T) for did in dids]
    return ms
