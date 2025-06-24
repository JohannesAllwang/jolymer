from jolymer.sas.GROMACS_SWAXS import GROMACS_SWAXS

samplename = 'ac6'
mdpath = f'/home/johannes/gromacs_xps/odna_2025/{samplename}/'
gs = GROMACS_SWAXS(
            path='/home/johannes/jophd/reports/odna/fastsaxs/',
            filename=f'{samplename}.dat',
            mdpath=mdpath,
        )


