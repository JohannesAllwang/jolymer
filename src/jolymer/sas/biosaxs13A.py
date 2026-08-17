import h5py
import datetime
from os.path import join
from pathlib import Path
import numpy as np
import pandas as pd

from .SAXS_Measurement import SAXS_Measurement


class biosaxs13A(SAXS_Measurement):
    """
    SAXS_Measurement subclass with additional methods specific to integrating detector immages at the 13A biosaxs beamline of the NSRRC in Taiwan
    """

    def get_poni_filename(self, filepath=None, M=9):
        parent_dir = Path(self.path).parent
        if M==9:
            outpath = Path(parent_dir, '9Mjolymer.poni')
        elif M==1:
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

    def get_Para(self, filepath=None):
        import re
        if filepath is None:
            filepath = self.get_filename()
        filepath = Path(filepath)
        header_dict = self.get_dat_header(filepath=filepath)
        filename_9M = header_dict["Sample filename"]
        # O64_master_0001.tif -> O64
        prefix = filename_9M.split("_")[0]
        match = re.fullmatch(r"([A-Z])(\d+)", prefix)
        if match is None:
            raise ValueError(f"Could not parse sample prefix: {prefix}")
        letter, number = match.groups()
        number = int(number)
        for i in range(1, 100):
            prev_number = number - i
            prev_letter = letter
            # O01 -> N99 -> N98 ...
            if prev_number < 1:
                prev_letter = chr(ord(letter) - 1)
                prev_number = 100 + prev_number
            candidate_prefix = f"{prev_letter}{prev_number:02d}"
            matches = list(
                filepath.parent.parent.glob(f"{candidate_prefix}_Para_*result.csv")
            )
            if matches:
                if i >= 5:
                    print(
                        f"WARNING: Para file is {i} samples backwards: "
                        f"{matches[0].name}"
                    )
                return pd.read_csv(matches[0])
        raise FileNotFoundError(
            f"No preceding Para result found for {prefix}"
        )

    def get_Para_time0(self):
        para = self.get_Para()
        filename_9M = self.get_dat_header(filepath=self.filepath)["Sample filename"]
        prefix = filename_9M.split("_")[0]
        row = para.loc[para["File name"].astype(str).str.strip() == prefix]
        if row.empty:
            raise ValueError(
                f"Could not find {prefix} in Para 'File name' column"
            )
        return float(row.iloc[0]["IntervalTime(s)"])

    def get_Para_times(self):
        "Time from the sample injection acording to Para*csv file"
        out_times = []
        Para_time0 = self.get_Para_time0()
        for frame_number in range(self.get_len_frames()):
            time = Para_time0 + frame_number * self.get_count_time()
            out_times.append(time)
        return out_times

    def get_len_frames(self, filepath=None, buffer=False):
        header_dict = self.get_dat_header(filepath=filepath)
        out = int(header_dict['Total Frame Numbers'])
        if buffer:
            with open(self.get_rigi_filename(filepath=filepath,
                                             buffer=buffer)) as f:
                out = sum(1 for _ in f)
                out = out / 7
        return out

    def get_rigi(self/, filepath=None, buffer=False):
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
                         'civiSMP': civiSMPs,
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

    def get_data_collection_dates(self, filepath=None, frame=0, time0=None):
        if filepath is None:
            filepath = self.get_dat_header()['Sample filename']
            filepath = filepath.split('_00')[0]+'.h5'
            filepath = join(self.path, '..', filepath)
            # print(filepath)
        with h5py.File(filepath, 'r') as hdf:
            date_bytes = hdf['entry']['instrument']['detector']['detectorSpecific']['data_collection_date'][()]
            data_group = hdf['entry']['data']
            nframes = sum(
                data_group[key].shape[0]
                for key in data_group.keys()
            )
            # print(nframes.keys())
            # print(np.shape(nframes.attr))
        date_str = date_bytes.decode('utf-8')
        if time0 is None:
            start_date = datetime.datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S.%f')
        else:
            start_date = time0
        frame_time = self.get_frame_time()
        frames = np.arange(nframes)
        times = frames * frame_time
        datetimes = [
            start_date + datetime.timedelta(seconds=float(t))
            for t in times
        ]
        return pd.DataFrame({
            'frame': frames,
            'time': times,
            'datetime': datetimes
        })

    def save_data_collection_dates(self, **kwargs):
        df = self.get_data_collection_dates(**kwargs)
        outpath = Path(self.path) / 'data_collection_dates.dat'
        df.to_csv(outpath, sep='\t', index=False)
        return df

    def save_time_lookup(self, filepath=None):
        get_data_collection_date(self, filepath=None, frame=0)

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
                    buffer=False, npt=200, mask=None,
                    normalize_by='civi'):
        import pyFAI
        from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
        if frame is not None:
            frame += 1
        ai = AzimuthalIntegrator()
        M = 1 if waxs else 9
        ai.load(self.get_poni_filename(M=M))
        path = self.get_9M_filename(buffer=buffer)
        sasImage = self.get_sasImage(path, frame=frame)
        if waxs:
            path = self.get_1M_filename(buffer=buffer)
            sasImage = self.get_sasImage(path, frame=frame)
        data = np.asarray(sasImage.data)
        if mask is None:
            # mask = np.zeros(data.shape, dtype=bool)
            mask = np.zeros(data.shape, dtype=np.int8)
        else:
            mask = mask.copy()
        med = np.median(data)
        mad = np.median(np.abs(data - med))
        # if mad > 0:
        if mad > 0:
            mask |= np.abs(data - med) > (300 * mad)
        q, I, errI = ai.integrate1d(data, npt=npt,
                                    mask=mask,
                                    error_model='poisson',
                                    method="csr",
                                    correctSolidAngle=True)
        dat_header = self.get_dat_header()
        T = float(dat_header['Sample Transmission coefficient'])
        if buffer:
            T = float(dat_header['Empty Cell Transmission coefficient'])
        count_time = 1
        if normalize_by == 'time':
            count_time = self.get_count_time(filepath=path)
        I = I / T / count_time
        errI = errI / T / count_time
        df = pd.DataFrame({'q': q, 'I_sample': I, 'err_I_sample': errI})
        return df

    def integrade_subtract(self, filename=None, filename_buffer=None, waxs=False, frame=None,
                           npt=200, adjustTMbuffer=1.0, normalize_by='civi', mask=None):
        df_sample = self.integrate1d(filename=filename, waxs=waxs, frame=frame,
                                     npt=npt, normalize_by=normalize_by, mask=mask)
        df_buffer = self.integrate1d(filename=filename, waxs=waxs,
                                     npt=npt, buffer=True, normalize_by=normalize_by,
                                     mask=mask)
        if normalize_by == 'time':
            rigi_sample = 100000
            rigi_buffer = 100000
        elif normalize_by == 'rigi':
            df_rigi_sample = self.get_rigi()
            rigi_sample = df_rigi_sample.loc[
                df_rigi_sample.frame_number == frame + 1,  # frame numbering convention
                'rigiSMP'
            ].iloc[0]
            df_rigi_buffer = self.get_rigi(buffer=True)
            # print(df_rigi_buffer)
            rigi_buffer = df_rigi_buffer.loc[
                df_rigi_buffer.frame_number == 1,  # frame numbering convention
                'rigiSMP'
            ].iloc[0]
        elif normalize_by == 'civi':
            df_rigi_sample = self.get_rigi()
            rigi_sample = df_rigi_sample.loc[
                df_rigi_sample.frame_number == frame + 1,  # frame numbering convention
                'civiSMP'
            ].iloc[0]
            df_rigi_buffer = self.get_rigi(buffer=True)
            # print(df_rigi_buffer)
            rigi_buffer = df_rigi_buffer.loc[
                df_rigi_buffer.frame_number == 1,  # frame numbering convention
                'civiSMP'
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

    def psaxs_to_poni(self, psaxs_file=None, poni_filename='9Mjolymer.poni',
                  detector="Eiger9M",
                  pixel_size=75e-6):
        if psaxs_file is None:
            psaxs_file = self.get_psaxs_file()
        poni_file = Path(self.path).parent / poni_filename
        with open(psaxs_file, "r") as f:
            lines = [l.strip() for l in f if l.strip()]
        energy_keV = float(lines[3].split()[0])
        beam_x, beam_y = map(float, lines[8].split('Beam Center')[0].split())
        distance_mm = float(lines[9].split("Sample-to")[0].split()[0])
        distance = distance_mm * 1e-3
        poni1 = beam_y * pixel_size
        poni2 = beam_x * pixel_size
        wavelength = 12.39842 / energy_keV * 1e-10
        with open(poni_file, "w") as f:
            f.write(
f"""# Converted from TPS pSAXS calibration
poni_version: 2
Detector: {detector}
Detector_config: {{}}
Distance: {distance}
Poni1: {poni1}
Poni2: {poni2}
Rot1: 0.0
Rot2: 0.0
Rot3: 0.0
Wavelength: {wavelength}
""")

    def get_psaxs_file(self):
        return Path(self.path).parent / 'pSAXSpar.txt'

    def get_reject_mask(self, filename='REJECT.dat', shape=None):
        """
        Convert TPS pSAXS REJECT.dat to pyFAI boolean mask.
        Parameters
        ----------
        filename : str
            REJECT.dat file
        shape : tuple
            Detector image shape (ny, nx)
        Returns
        -------
        mask : ndarray(bool)
            True = masked pixel
        """
        if shape is None:
            path = self.get_9M_filename(buffer=False)
            shape = self.get_sasImage(path).shape
            print('shape', shape)
        mask = np.zeros(shape, dtype=bool)
        reject_path = Path(self.path).parent / filename
        with open(reject_path) as f:
            for line in f:
                if not line.strip():
                    continue
                x, y, value = line.split()
                x = int(float(x))
                y = int(float(y))
                if 0 <= y < shape[0] and 0 <= x < shape[1]:
                    mask[y, x] = True
        return mask.astype(np.int8)
