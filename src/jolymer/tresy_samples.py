from .Sample import Sample

class TresySample(Sample):

    def __init__(self, tresy_id):

        self.type = 'tresy_sample'
        with dbo.dbopen() as c:
            c.execute("SELECT * FROM tresy_sample WHERE id=?", (tresy_id,))
            values = c.fetchone()
        self.id, buffer_id, polysacharide_id, polysaccharide_gpl, self.tt, loaded, self.comment = values
        self.loaded = loaded
        if loaded == 1:
            self.loaded_dict = {
                    'material' : 'curcumin',
                    'gpl' : 0.1 }
        self.datestring = "20210527"
        self.protein = Protein(1)
        self.protein_gpl = 1
        self.PS = Polysaccharide(polysaccharide_id)
        self.PS_gpl = polycaccharide_gpl
        self.buffer = TresyBuffer(buffer_id)

    def get_NPname(self):
        out = f'{self.PS.short_name}/{self.protein.short_name}'
        return out
        
class TresyBuffer(Sample):

    def __init__(self, buffer_id):
        self.type = 'tresy_buffer'
        with dbo.dbopen() as c:
            c.execute("SELECT * FROM tresy_buffers WHERE id=?", (buffer_id,))
            values = c.fetchone()
        self.id, self.pH, self.salt_concentration = values

    def get_viscosity(self, TK):
        print("using the viscosity of water")
        if TK == 293.15:
            return 0.0010016 

    def get_n(self, TK):
        print("using the refractive index of water")
        if TK == 293.15:
            return 1.333


