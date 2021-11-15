import numpy as np
import .SAXS_Model as sasmodel


def ellipsoid(q, Ra, Rb, SLD=1, solventSLD=0, alpha=None, tol=1e-6):
    r"""
    Form factor for a simple ellipsoid (ellipsoid of revolution).

    Parameters
    ----------
    q : float
        Scattering vector unit e.g.  1/A or 1/nm  1/Ra
    Ra : float
        Radius rotation axis   units in 1/unit(q)
    Rb : float
        Radius rotated axis    units in 1/unit(q)
    SLD : float, default =1
        Scattering length density of unit nm^-2
        e.g. SiO2 = 4.186*1e-6 A^-2 = 4.186*1e-4 nm^-2 for neutrons
    solventSLD : float, default =0
        Scattering length density of solvent. unit nm^-2
        e.g. D2O = 6.335*1e-6 A^-2 = 6.335*1e-4 nm^-2 for neutrons
    alpha : [float,float] , default [0,90]
        Angle between rotation axis Ra and scattering vector q in unit grad
        Between these angles orientation is averaged
        alpha=0 axis and q are parallel, other orientation is averaged
    tol : float
        relative tolerance for integration between alpha

    Returns
    -------
    dataArray
        Columns [q; Iq; beta ]
         - .RotationAxisRadius
         - .RotatedAxisRadius
         - .EllipsoidVolume
         - .I0         forward scattering q=0
         - beta is asymmetry factor according to [3]_.
           :math:`\beta = |<F(Q)>|^2/<|F(Q)|^2>` with scattering amplitude :math:`F(Q)` and
           form factor :math:`P(Q)=<|F(Q)|^2>`

    Examples
    --------
    Simple ellipsoid in vacuum::

     import jscatter as js
     import numpy as np
     x=np.r_[0.1:10:0.01]
     Rp=6.
     Re=8.
     ashell=js.ff.multiShellEllipsoid(x,Rp,Re,1)
     #plot it
     p=js.grace()
     p.new_graph(xmin=0.24,xmax=0.5,ymin=0.2,ymax=0.5)
     p[1].subtitle('contrastprofile')
     p[0].plot(ashell)
     p[0].yaxis(scale='l',label='I(q)',min=0.01,max=100)
     p[0].xaxis(scale='l',label='q / nm\S-1',min=0.1,max=10)
     p[0].title('ellipsoid')
     p[1].plot(ashell.contrastprofile,li=1) # a contour of the SLDs
     #p.save(js.examples.imagepath+'/ellipsoid.jpg')

    .. image:: ../../examples/images/ellipsoid.jpg
     :width: 50 %
     :align: center
     :alt: ellipsoid


    References
    ----------
    .. [1] Structure Analysis by Small-Angle X-Ray and Neutron Scattering
           Feigin, L. A, and D. I. Svergun, Plenum Press, New York, (1987).
    .. [2] http://www.ncnr.nist.gov/resources/sansmodels/Ellipsoid.html
    .. [3] M. Kotlarchyk and S.-H. Chen, J. Chem. Phys. 79, 2461 (1983).

    """
    if alpha is None:
        alpha = [0, 90]
    result = multiShellEllipsoid(q, Ra, Rb, shellSLD=SLD, solventSLD=solventSLD, alpha=alpha, tol=tol)
    attr = result.attr
    result.EllipsoidVolume = result.outerVolume
    result.RotationAxisRadius = Ra
    result.RotatedAxisRadius = Rb
    result.contrast = result.shellcontrast
    result.angles = alpha
    attr.remove('columnname')
    attr.remove('I0')
    for at in attr:
        delattr(result, at)
    result.modelname = inspect.currentframe().f_code.co_name
    return result
