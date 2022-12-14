from .beaucage import beaucage, beaucage_bg, fw_beaucage_bg
from .gauss_lorentz_gel import gauss, lorentz, gauss_lorentz, gauss_lorentz_bg, fw_gauss_lorentz, fw_gauss_lorentz_bg
from .guinierPorod import gupo, gupo_forward
from .ornsteinZernike import oz, oz_porod, oz_bg
from .double_beaucage import dbeaucage
from .triple_beaucage import tbeaucage
from .SAXS_Model import porod, bg, fw, AddSaxsModels
from .gaussian_peak import gel_model
from .dab import dab
from .sphere import sphere


beaucage_gauss_lorentz_bg = beaucage_bg.plusmodel(gauss_lorentz, name = 'beaucage_gauss_lorentz_bg', longname = 'Beaucage and gauss lorentz gel')
# beaucage_gauss_lorentz_bg.get_text = get_text_beaucage_gauss_lorentz_bg
