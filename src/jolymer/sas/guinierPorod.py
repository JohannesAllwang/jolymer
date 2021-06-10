# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:19:12 2021

@author: xcill
"""

from . import SAXS_Model as sasmodel
import numpy as np

def genGuinier(q, Rg, A, alpha):
    r"""
    Generalized Guinier approximation for low wavevector q scattering q*Rg< 1-1.3

    For absolute scattering see introduction :ref:`formfactor (ff)`.

    Parameters
    ----------
    q : array of float
        Wavevector
    Rg : float
        Radius of gyration in units=1/q
    alpha : float
        Shape [α = 0] spheroid,    [α = 1] rod-like    [α = 2] plane
    A : float
        Amplitudes

    Returns
    -------
    dataArray
        Columns [q,Fq]

    Notes
    -----
    Quantitative analysis of particle size and shape starts with the Guinier approximations.
     - For three-dimensional objects the Guinier approximation is given by
       :math:`I(q) = A e^{-Rg^2q^2/3}`
     - This approximation can be extended also to rod-like and plane objects by
       :math:`I(q) =(\alpha \pi q^{-\alpha})  A e^{-Rg^2q^2/(3-\alpha) }`

    If the particle has one dimension of length L that is much larger than
    the others (i.e., elongated, rod-like, or worm-like), then there is a q
    range such that qR_c < 1 <<  qL, where α = 1.

    Examples
    --------
    ::

     import jscatter as js
     import numpy as np
     q=js.loglist(0.01,5,300)
     spheroid=js.ff.genGuinier(q, Rg=2, A=1, alpha=0)
     rod=js.ff.genGuinier(q, Rg=2, A=1, alpha=1)
     plane=js.ff.genGuinier(q, Rg=2, A=1, alpha=2)
     p=js.grace()
     p.plot(spheroid,le='sphere')
     p.plot(rod,le='rod')
     p.plot(plane,le='plane')
     p.yaxis(scale='l',min=1e-4,max=1e4)
     p.xaxis(scale='l')
     p.legend(x=0.03,y=0.1)
     #p.save(js.examples.imagepath+'/genGuinier.jpg')

    .. image:: ../../examples/images/genGuinier.jpg
     :align: center
     :width: 50 %
     :alt: genGuinier


    References
    ----------
    .. [1] Form and structure of self-assembling particles in monoolein-bile salt mixtures
           Rex P. Hjelm, Claudio Schteingart, Alan F. Hofmann, and Devinderjit S. Sivia
           J. Phys. Chem., 99:16395--16406, 1995

    """
    q = np.atleast_1d(q)
    if alpha == 0:
        pre = 1
    elif alpha == 1 or alpha == 2:
        pre = alpha * np.pi * q ** -alpha
    else:
        raise TypeError('alpha needs to be in 0,1,2')
    I = pre * A * np.exp(-Rg ** 2 * q ** 2 / (3 - alpha))
    return I





def guinierPorod3d(q, Rg1, s1, Rg2, s2, G2, dd):
    r"""
    Generalized Guinier-Porod Model with high Q power law with 3 length scales.

    The model represents the most general case containing three Guinier regions [1]_.

    Parameters
    ----------
    q : float
        Wavevector  in units of 1/nm
    Rg1 : float
        Radii of gyration for the short size of scattering object in units nm.
    Rg2 : float
        Radii of gyration for the overall size of scattering object in units nm.
    s1 : float
        Dimensionality parameter for the short size of scattering object (s1=1 for a cylinder)
    s2 : float
        dimensionality parameter for the overall size of scattering object (s2=0 for a cylinder)
    G2 : float
        Intensity for q=0.
    dd : float
        Porod exponent

    Returns
    -------
    dataArray
        Columns [q,Iq]
         Iq scattering intensity

    Notes
    -----
    For a cylinder with length L and radius R (see [1]_)
    :math:`R_{g2} = (L^2/12+R^2/2)^{\frac{1}{2}}`  and :math:`R_{g1}=R/\sqrt{2}`


    Examples
    --------
    ::

     import jscatter as js
     q=js.loglist(0.01,5,300)
     I=js.ff.guinierPorod3d(q,Rg1=1,s1=1,Rg2=10,s2=0,G2=1,dd=4)
     p=js.grace()
     p.plot(I)
     p.xaxis(scale='l',label='q / nm\S-1')
     p.yaxis(scale='l',label='I(q) / a.u.')
     #p.save(js.examples.imagepath+'/guinierPorod3d.jpg')

    .. image:: ../../examples/images/guinierPorod3d.jpg
     :align: center
     :width: 50 %
     :alt: guinierPorod3d

    References
    ----------
    .. [1]  A new Guinier/Porod Model
            B. Hammouda J. Appl. Cryst. (2010) 43, 716-719

    Author M. Kruteva JCNS 2019

    """
    q = np.atleast_1d(q)

    # define parameters for smooth transitions
    Q1 = (1 / Rg1) * ((dd - s1) * (3 - s1) / 2) ** 0.5
    Q2 = ((s1 - s2) / (2 / (3 - s2) * Rg2 ** 2 - 2 / (3 - s1) * Rg1 ** 2)) ** 0.5
    G1 = G2 / (np.exp(-Q2 ** 2 * (Rg1 ** 2 / (3 - s1) - Rg2 ** 2 / (3 - s2))) * Q2 ** (s2 - s1))
    D = G1 * np.exp(-Q1 ** 2 * Rg1 ** 2 / (3 - s1)) * Q1 ** (dd - s1)

    # define functions in different regions
    def _I1_3regions(q):
        res = G2 / q ** s2 * np.exp(-q ** 2 * Rg2 ** 2 / (3 - s2))
        return res

    def _I2_3regions(q):
        res = G1 / q ** s1 * np.exp(-q ** 2 * Rg1 ** 2 / (3 - s1))
        return res

    def _I3_3regions(q):
        res = D / q ** dd
        return res

    I = np.piecewise(q, [q < Q2, (Q2 <= q) & (q < Q1), q >= Q1], [_I1_3regions, _I2_3regions, _I3_3regions])
    return I


def guinierPorod(q, Rg, s, G, dd):
    r"""
    Generalized Guinier-Porod Model with high Q power law.

    Parameters
    ----------
    q : float
        Wavevector  in units of 1/nm
    Rg : float
        Radii of gyration in units nm.
    s : float
        Dimensionality parameter describing the low Q region.
    dd : float
        Porod exponent describing the high Q slope.
    G : float
        intensity

    Returns
    -------
    dataArray
        Columns [q,Iq]
        Iq    scattering intensity

    Examples
    --------
    ::

     import jscatter as js
     q=js.loglist(0.01,5,300)
     I=js.ff.guinierPorod(q,s=0,Rg=5,G=1,dd=4)
     p=js.grace()
     p.plot(I)
     p.xaxis(scale='l',label='q / nm\S-1')
     p.yaxis(scale='l',label='I(q) / a.u.')
     #p.save(js.examples.imagepath+'/guinierPorod.jpg')

    .. image:: ../../examples/images/guinierPorod.jpg
     :align: center
     :width: 50 %
     :alt: guinierPorod

    References
    ----------
    .. [1]  A new Guinier/Porod Model
            B. Hammouda J. Appl. Cryst. (2010) 43, 716-719


    Author M. Kruteva JCNS 2019
    """
    q = np.atleast_1d(q)

    # define parameters for smooth transitions
    Q1 = (1 / Rg) * ((dd - s) * (3 - s) / 2) ** 0.5
    D = G * np.exp(-Q1 ** 2 * Rg ** 2 / (3 - s)) * Q1 ** (dd - s)

    # define functions in different regions
    def _I1_2regions(q):
        res = G / q ** s * np.exp(-q ** 2 * Rg ** 2 / (3 - s))
        return res

    def _I2_2regions(q):
        res = D / q ** dd
        return res

    I = np.piecewise(q, [q < Q1, q >= Q1], [_I1_2regions, _I2_2regions])
    return I


class GuinierPorod(sasmodel.SAXS_Model):
    
    def __init__(self):
        self.name = 'gupo'
        self.longname = 'Guinier Porod'
        self.parameters = ['Rg', 's', 'scale', 
                           'porod_exp']
        self.pdict = {'xi' : ['$\\xi$', 'nm']}
        self.fitfunc = guinierPorod
        
    def get_clustering_strength(self, fitdict, q):
        A = fitdict['porod_scale']
        n = fitdict['porod_exp']
        out = A/q**n
        return out
    

background = sasmodel.Background()
gupo = GuinierPorod().plusmodel(background)
def get_text_gupo(fit_dict):
        text = """
        $R_g =$ {0:.2f} $\\pm$ {1:.2f} nm
        $n =$ {2:.2f} $\\pm$ {3:.2f}
        $C =$ {4:.2E}
        $\\chi^2 = $ {5:.4}
        """.format(fit_dict['Rg'], fit_dict['std_Rg'],
                   fit_dict['porod_exp'], fit_dict['std_porod_exp'], 
                   fit_dict['scale'],
                  fit_dict['chi2'])
        return text
gupo.get_text = get_text_gupo

forward = sasmodel.Porod()
forward.parameters = ['fw_scale', 'fw_exp']

gupo_forward = gupo.plusmodel(forward, name='gupo_fw', longname = 'Guinier Porod')
def get_text_gupof(fit_dict):
        text = """
        $R_g =$ {0:.2f} $\\pm$ {1:.2f} nm
        $n =$ {2:.2f} $\\pm$ {3:.2f}
        $m_F = $ {4:.2f} $\\pm$ {5:.2f}
        $C =$ {6:.2E}
        $A_F = $ {7:.2E}
        $\\chi^2 = $ {8:.4}
        """.format(fit_dict['Rg'], fit_dict['std_Rg'],
                   fit_dict['porod_exp'], fit_dict['std_porod_exp'], 
                   fit_dict['fw_exp'], fit_dict['std_fw_exp'],
                   fit_dict['scale'], fit_dict['fw_scale'],
                  fit_dict['chi2'])
        return text
gupo_forward.get_text = get_text_gupof