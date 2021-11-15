from . import SAXS_Model as sasmodel
import numpy as np
from scipy import special


def _fa_sphere(qr):
    """
    scattering amplitude sphere with catching the zero
    qr is array dim 1
    """
    fa=np.ones(qr.shape)
    qr0 = (qr!=0)
    fa[qr0] = 3 / qr[qr0] ** 3 * (np.sin(qr[qr0]) - qr[qr0] * np.cos(qr[qr0]))
    return fa


def sphere(q, radius, contrast=1):
    r"""
    Scattering of a single homogeneous sphere.

    Parameters
    ----------
    q : float
        Wavevector  in units of 1/nm
    radius : float
        Radius in units nm
    contrast : float, default=1
        Difference in scattering length to the solvent = contrast

    Returns
    -------
    dataArray
        Columns [q, Iq, fa]
        Iq    scattering intensity
        - fa formfactor amplitude
        - .I0   forward scattering


    Notes
    -----
    .. math:: I(q)=  4\pi\rho^2V^2\left[\frac{3(sin(qR) - qR cos(qR))}{(qR)^3}\right]^2

    with contrast :math:`\rho` and sphere volume :math:`V=\frac{4\pi}{3}R^3`

    The first minimum of the form factor is at qR=4.493

    Examples
    --------
    ::

     import jscatter as js
     import numpy as np
     q=js.loglist(0.1,5,300)
     p=js.grace()
     R=3
     sp=js.ff.sphere(q, R)
     p.plot(sp.X*R,sp.Y,li=1)
     p.yaxis(label='I(q)',scale='l',min=1e-4,max=1e5)
     p.xaxis(label='qR',scale='l',min=0.1*R,max=5*R)
     p.legend(x=0.15,y=0.1)
     #p.save(js.examples.imagepath+'/sphere.jpg')

    .. image:: ../../examples/images/sphere.jpg
     :align: center
     :width: 50 %
     :alt: sphere


    References
    ----------
    .. [1] Guinier, A. and G. Fournet, "Small-Angle Scattering of X-Rays", John Wiley and Sons, New York, (1955).

    """
    R = radius
    qr = np.atleast_1d(q) * R
    fa0 = (4 / 3. * np.pi * R ** 3 * contrast)  # forward scattering amplitude q=0
    faQR = fa0 * _fa_sphere(qr)
    result = dA(np.c_[q, faQR** 2, faQR].T)
    result.columnname = 'q; Iq; fa'
    result.setColumnIndex(iey=None)
    result.radius = radius
    result.I0 = fa0**2
    result.fa0 = fa0
    result.contrast = contrast
    return result

def sphere_fitfunc(q, scale, R, V=1, Drho=1):
    """
    This is the equation from the sasview documentation.
    The actual sasview implementation uses C and I don't know how to do that myself.
    """
    qR = q * R
    Iq = (np.sin(qR) - qR * np.cos(qR)) / (qR ** 3)
    Iq = 3 * V * Drho * Iq
    Iq = scale * Iq * Iq / V
    return Iq

class Sphere(sasmodel.SAXS_Model):

    def __init__(self):
        self.name = 'sphere'
        self.longname = 'Sphere'
        self.parameters = ['sphere_scale', 'sphere_radius']
        self.fitfunc = sphere_fitfunc

    def get_text(self, fit_dict):
        text = """
        $R = {:.3f} \\pm {:.3f}\\mathrm{{\\,nm}}
        $A_{{Sphere}} = {:.3e}$
        """.format(fit_dict['sphere_radius'], fit_dict['std_sphere_radius'],
                fit_dict['sphere_scale'])
        return text

sphere = Sphere()
