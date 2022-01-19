from .lsi import lsi, lsi3d


class PhiDLS(lsi):

    def __init__(self, phid):
        self.path = ''
        self.solvent_m = ''
        self.toluene_m = ''
