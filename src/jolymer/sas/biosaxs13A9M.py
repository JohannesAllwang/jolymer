import h5py
import datetime
from os.path import join

from .SAXS_Measurement import SAXS_Measurement


class biosaxs13A9M(SAXS_Measurement):

    def get_dat_header(self, filepath=None):
        if filepath is None:
            filepath = self.get_filename()
        outdict = {}
        with open(filepath) as f:
            for line in f:
                value = line.strip().split('    ')[0]
                key = line[len(value)+1::].strip()
                outdict[key] = value
                if line.replace(' ', '')=='QIERROR\n':
                    break
            return outdict

    def get_data_collection_date(self, filepath=None, frame=0):
        if filepath is None:
            filepath = self.get_dat_header()['Sample filename']
            filepath = filepath.split('_00')[0]+'.h5'
            filepath = join(self.path, '..', filepath)
            print(filepath)
        with h5py.File(filepath, 'r') as hdf:
            date_bytes = hdf['entry']['instrument']['detector']['detectorSpecific']['data_collection_date'][()]
        date_str = date_bytes.decode('utf-8')
        date = datetime.datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S.%f')
        if frame > 0:
            frame_time = self.get_frame_time(filepath=filepath)
            print('timeshift', frame, frame_time)
            frameshift = datetime.timedelta(seconds=frame*frame_time)
            date = date + frameshift
        return date

    def get_count_time(self, filepath=None):
        if filepath is None:
            filepath = self.get_dat_header()['Sample filename']
            filepath = filepath.split('_00')[0]+'.h5'
            filepath = join(self.path, '..', filepath)
            print(filepath)
        with h5py.File(filepath, 'r') as hdf:
            length = hdf['entry']['instrument']['detector']['count_time'][()]
        return length

    def get_frame_time(self, filepath=None):
        if filepath is None:
            filepath = self.get_dat_header()['Sample filename']
            filepath = filepath.split('_00')[0]+'.h5'
            filepath = join(self.path, '..', filepath)
            print(filepath)
        with h5py.File(filepath, 'r') as hdf:
            length = hdf['entry']['instrument']['detector']['frame_time'][()]
        return length

