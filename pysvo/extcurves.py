##################################################################################
#
#   pysvo.extcurves: rudimentary extinction curve module
#
#   Main Functions:
#
#       - Cardelli1989:    Cardelli et al. (1989) extinction curve
#       - Schlafly2016:    Schlafly et al. (2016) extinction curve
#
#   For more options, use the extinction module available on pip.
#
###################################################################################
import numpy as np
from pysvo import fitting_tools

"""
Defining the transformation constants between a set of filters
relative to to the V-band extinction.
"""

def Cardelli1989( wavelength, RV=3.1 ):
    """
    Get A(lambda) / AV using recipe of Cardelli et al (1989).

    Input:
        wavelength - in Angstroms
    Optional:
        RV         - the preferred Rv value
    Output:
        Extinction relative to V-band: A(lambda) / AV
        
    History:
        2016       - Written             - F. Anders (AIP)
        2019-06-10 - Renamed for starsed - F. Anders (ICCUB)
    """
    x = 10000. / wavelength     # x is in (mum)^-1
    # Different coefficients for UV, optical, and infrared waveleghts:
    # opt: 0.25 - 0.9 micron, IR: 0.9 - 3.3 micron
    a = np.piecewise( x, [x<1.1, x>=1.1],
                      [lambda x: 0.574 * x**1.61,
                       lambda x: 1 + 0.17699*(x-1.82)  - 0.50447*(x-1.82)**2
                                 - 0.02427*(x-1.82)**3 + 0.72085*(x-1.82)**4
                                 + 0.01979*(x-1.82)**5 - 0.77530*(x-1.82)**6
                                 + 0.32999*(x-1.82)**7  ] )
    b = np.piecewise( x, [x<1.1, x>=1.1],
                      [lambda x: -0.527 * x**1.61,
                       lambda x:   1.41338*(x-1.82)    + 2.28305*(x-1.82)**2
                                 + 1.07233*(x-1.82)**3 - 5.38434*(x-1.82)**4
                                 - 0.62251*(x-1.82)**5 + 5.30260*(x-1.82)**6
                                 - 2.09002*(x-1.82)**7  ] )
    return a + b / RV

def Schlafly2016( wavelength, RV=3.32 ):
    """
    Based on E. Schlafly's code: http://e.schlaf.ly/apored/extcurve.html
    
    Returns the extinction curve, A(lambda)/A(5420 A), according to
    Schlafly+2016. RV is transformed to the parameter "x," which controls the overall shape of
    the extinction curve in an R(V)-like way. The extinction curve is based on broad band photometry between the PS1 g
    band and the WISE W2 band, which have effective wavelengths between 5000
    and 45000 A, and is blindly extrapolated outside that
    range.  The gray component of the extinction curve is fixed by enforcing
    A(H)/A(K) = 1.55 (Indebetouw+2005). 

    Args:
        wavelength - in Angstrom
    Optional:
        RV         - default is 3.32 (Schlafly+2016, stdev 0.18)

    Returns:
        A(lambda)/A(5420 A)
        
    History:
        2017       - adapted for starhorse - F. Anders (AIP)
        2019-06-10 - ported to starsed     - F. Anders (ICCUB)
    """
    
    # Transform RV' to Schlafly's x parameter
    x = (RV - 3.3) / 9.1
    #   ra: extinction vector at anchor wavelengths
    #   dra: derivative of extinction vector at anchor wavelengths
    #   lam: anchor wavelengths (angstroms)
    ra = np.array([ 0.65373283,  0.39063843,  0.20197893,  0.07871701, -0.00476316,
                   -0.14213929, -0.23660605, -0.28522577, -0.321301  , -0.33503192])
    dra = np.array([-0.54278669,  0.03404903,  0.36841725,  0.42265873,  0.38247769,
                     0.14148814, -0.04020524, -0.13457319, -0.26883343, -0.36269229])
    lam = np.array([  5032.36441067,   6280.53335141,   7571.85928312,   8690.89321059,
                          9635.52560909,  12377.04268274,  16381.78146718,  21510.20523237,
                         32949.54009328,  44809.4919175 ])
    # Indebetouw (2005) A_H / A_K
    rhk = 1.55 

    anchors = ra + x*dra
    # fix gray component so that A(H)/A(K) = 1.55
    anchors += (-anchors[6] + rhk*anchors[7])/(1 - rhk)
    cs0 = fitting_tools.CubicSpline(lam, anchors, yp='3d=0')
    # normalize at 5420 angstroms
    ec  = fitting_tools.CubicSpline(lam, anchors/cs0(5420.), yp='3d=0')
    return ec( wavelength )
