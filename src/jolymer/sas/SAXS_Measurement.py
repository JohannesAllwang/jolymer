"""
"""

from scipy import optimize
import numpy as np
from dataclasses import dataclass

from .. import database_operations as dbo
from ..Measurement import Measurement

@dataclass
class SAXS_Measurement(Measurement):

    instrument = 'unspecified'
    sample: str = None
    mid: str = None
    datestring: str = None
    given_name: str = None
    samplestring: str = None
    comment: str = None
    timestring: str = None
    rawdatapath: str = None

    def __post_init__(self):
        if self.samplestring == None:
            self.sample = None
        else:
            sample_type, sample_id = self.samplestring.split('_')
            self.sample = Sample.get_sample(sample_type, sample_id)
        self.path = os.path.join(self.rawdatapath,
                                 '{self.instrument}{0:03}/'.format(self.mid))
        self.processed_subtracted_file = os.path.join(self.path, 'processed_subtracted.dat')
        self.processed_file = os.path.join(self.path, 'processed.dat')
        self.frames_path = os.path.join(self.path, 'frames')
        self.absolute_path = os.path.join(self.path, 'absolute')
        self.averaged_path = os.path.join(self.path, 'averaged')
        self.buffer_frames_path = os.path.join(self.path,'buffer', 'frames')
        self.buffer_absolute_path = os.path.join(self.path,'buffer', 'absolute')
        self.origpath = os.path.join(self.rawdatapath,
                                     f'{self.instrument}{self.datestring}',
                                     'datacollection', 'data', 'absolute')

    def get_filename():
        pass

    def get_data(self, cout=True, altername='No', nrows=2653):
        "get data from data.csv and apply some filter."
        if altername=='No':
            path = self.get_filename()
        else:
            path = os.path.join(self.path, f'{altername}.dat')
        df = pd.read_csv(path, sep='\s+', header=None,
                         skiprows=3, nrows=nrows, names=['q', 'I', 'err_I'])
        # log.log(f'{len_before-len_after} negative I values!')
        # log.log(f'{len_before-len_after} excluded I values!')
        return df
        pass

    @staticmethod
    def get_distribution(df):
        pass
