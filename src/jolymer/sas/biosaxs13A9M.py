import h5py
import datetime
from os.path import join
from pathlib import Path
import numpy as np
import pandas as pd

from .SAXS_Measurement import SAXS_Measurement


class biosaxs13A9M(SAXS_Measurement):

    def get_poni_filename(self, filepath=None):
        parent_dir = Path(self.path).parent
        outpath = Path(parent_dir, '1M.poni')
        return outpath

    def get_rigi_filename(self, filepath=None, buffer=False):
        header_dict = self.get_dat_header(filepath=filepath)
        rigi_name = header_dict['Sample information (beam flux)']
        if buffer:
            rigi_name = header_dict['Empty Cell information (beam flux)']
        parent_dir = Path(self.path).parent
        outpath = Path(parent_dir, rigi_name)
        return outpath

    def get_9M_filename(self, filepath=None, buffer=False):
        header_dict = self.get_dat_header(filepath=filepath)
        # print(header_dict)
        filename_9M = header_dict['Sample filename']
        if buffer:
            filename_9M = header_dict['Empty Cell filename']
        filename_9M = f'{filename_9M[0:3]}_master.h5'
        parent_dir = Path(self.path).parent
        outpath = Path(parent_dir, filename_9M)
        return outpath

    def get_1M_filename(self, filepath=None, buffer=False):
        header_dict = self.get_dat_header(filepath=filepath)
        # print(header_dict)
        filename_9M = header_dict['Sample filename']
        if buffer:
            filename_9M = header_dict['Empty Cell filename']
        filename_1M = f'{filename_9M[1:3]}{filename_9M[0]}_master.h5'
        parent_dir = Path(self.path).parent
        outpath = Path(parent_dir, filename_1M)
        return outpath

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

    def get_len_frames(self, filepath=None, buffer=False):
        header_dict = self.get_dat_header(filepath=filepath)
        out = int(header_dict['Total Frame Numbers'])
        if buffer:
            with open(self.get_rigi_filename(filepath=filepath,
                                             buffer=buffer)) as f:
                out = sum(1 for _ in f)
                out = out / 7
        return out

    def get_rigi(self, filepath=None, buffer=False):
        header_dict = self.get_dat_header(filepath=filepath)
        len_frames = self.get_len_frames(filepath=filepath, buffer=buffer)
        civiSMPs, rigiSMPs, expSMPs = [], [], []
        with open(self.get_rigi_filename(filepath=filepath, buffer=buffer)) as rigi_file:
            for i, line in enumerate(rigi_file):
                if i < len_frames:
                    pass
                elif i < 2*len_frames:
                    civiSMPs.append(float(line))
                elif i < 3*len_frames:
                    rigiSMPs.append(float(line))
                elif i < 4*len_frames:
                    pass
                elif i < 5*len_frames:
                    expSMPs.append(float(line))
                elif i < 6*len_frames:
                    pass
        out = pd.DataFrame({'frame_number' : np.linspace(1,len_frames,
                                                         num=len(civiSMPs)),
                         'civiSMS': civiSMPs,
                         'rigiSMP': rigiSMPs,
                         'expSMP': expSMPs})
        return out

    def get_data_collection_date(self, filepath=None, frame=0):
        if filepath is None:
            filepath = self.get_dat_header()['Sample filename']
            filepath = filepath.split('_00')[0]+'.h5'
            filepath = join(self.path, '..', filepath)
            # print(filepath)
        with h5py.File(filepath, 'r') as hdf:
            date_bytes = hdf['entry']['instrument']['detector']['detectorSpecific']['data_collection_date'][()]
        date_str = date_bytes.decode('utf-8')
        date = datetime.datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S.%f')
        if frame > 0:
            frame_time = self.get_frame_time(filepath=filepath)
            # print('timeshift', frame, frame_time)
            frameshift = datetime.timedelta(seconds=float(frame*frame_time))
            date = date + frameshift
        return date

    def get_count_time(self, filepath=None):
        """
        This is the exposure time
        """
        if filepath is None:
            filepath = self.get_dat_header()['Sample filename']
            filepath = filepath.split('_00')[0]+'.h5'
            filepath = join(self.path, '..', filepath)
            # print(filepath)
        with h5py.File(filepath, 'r') as hdf:
            length = hdf['entry']['instrument']['detector']['count_time'][()]
        return length

    def get_frame_time(self, filepath=None):
        """
        this is the exposure time plus the wait after it
        """
        if filepath is None:
            filepath = self.get_dat_header()['Sample filename']
            filepath = filepath.split('_00')[0]+'.h5'
            filepath = join(self.path, '..', filepath)
            # print(filepath)
        with h5py.File(filepath, 'r') as hdf:
            length = hdf['entry']['instrument']['detector']['frame_time'][()]
        return length

    def integrate1d(self, filename=None, waxs=False, frame=None,
                    buffer=False, npt=200):
        import pyFAI
        from pyFAI.azimuthalIntegrator import AzimuthalIntegrator

        if not frame is None:
            frame += 1 # because it usually starts to count at 1
        ai = AzimuthalIntegrator()
        ai.load(self.get_poni_filename())

        path = self.get_9M_filename(buffer=buffer)
        sasImage = self.get_sasImage(path, frame=frame)
        if waxs:
            path = self.get_1M_filename(buffer=buffer)
            sasImage = self.get_sasImage(path, frame=frame)
        data = sasImage.data
        med = np.median(data)
        mad = np.median(np.abs(data - med))
# robust threshold
        mask = np.abs(data - med) > (8 * mad)
        data = np.ma.array(data, mask=mask)
        # data = np.ma.masked_greater(sasImage.data, 100000)
        # print('data is masked')
        q, I, errI = ai.integrate1d(data, npt=npt,
                                  mask=data.mask,
                                  error_model='poisson',
                                  correctSolidAngle=True)
        dat_header = self.get_dat_header()
        T = float(dat_header['Sample Transmission coefficient'])
        if buffer:
            T = float(dat_header['Empty Cell Transmission coefficient'])
        # count_time = self.get_count_time(filepath=path)
        count_time = 1 # wile using rigi to normalize
        I = I / T / count_time
        errI = errI/T/count_time
        df = pd.DataFrame({'q': q, 'I_sample': I, 'err_I_sample': errI})
        return df

    def integrade_subtract(self, filename=None, filename_buffer=None, waxs=False, frame=None,
                           npt=200, adjustTMbuffer=1.0):
        df_sample = self.integrate1d(filename=filename, waxs=waxs, frame=frame,
                                     npt=npt)
        df_rigi_sample = self.get_rigi()
        rigi_sample = df_rigi_sample.loc[
            df_rigi_sample.frame_number == frame + 1,  # frame numbering convention
            'rigiSMP'
        ].iloc[0]
        df_buffer = self.integrate1d(filename=filename, waxs=waxs,
                                     npt=npt, buffer=True)
        df_rigi_buffer = self.get_rigi(buffer=True)
        # print(df_rigi_buffer)
        rigi_buffer = df_rigi_buffer.loc[
            df_rigi_buffer.frame_number == 1,  # frame numbering convention
            'rigiSMP'
        ].iloc[0]
        # print('rigi sample:', rigi_sample)
        # print('rigi buffer:', rigi_buffer)
        df = df_sample.copy()
        df['I_sample'] = df.I_sample / rigi_sample
        df['err_I_sample'] = df.err_I_sample / rigi_sample
        df['q_buffer'] = df_buffer.q
        df['I_buffer'] = df_buffer.I_sample / rigi_buffer
        df['err_I_buffer'] = df_buffer.err_I_sample / rigi_buffer
        df['I'] = df.I_sample - df.I_buffer/adjustTMbuffer
        df['err_I'] = np.sqrt(df.err_I_sample**2 + (df.err_I_buffer/adjustTMbuffer)**2)
        # df['err_I'] = df.err_I_sample + df.err_I_buffer
        return df



