##################################################################################
#
#   pysvo.fitting_tools: some simple fitting tools
#
#   Main Functions:
#
#       - polyfit2d:    2D polynomial fitter, using numpy.polyvander2d
#       - splint:       Simple spline interpolator
#
#   Classes:
#       - CubicSpline:  Wrapper class for slpint
#
#   For more options, use the extinction module available on pip.
#
###################################################################################
import numpy as np
from scipy.linalg import solve_banded
from numpy.polynomial import polynomial

def polyfit2d(x, y, f, deg):
    """
    Performs a simple 2D polynomial fit using numpy.polyvander2d

    Input:
        x, y    - coordinate vectors
        f       - function value at respective coordinate tuple
        deg     - degrees of polynomial to be fitted, e.g. [2,2]
    Output:
        Array of polynomial coefficients
        
    History:
       2017     - Written - Anders (AIP)
    """
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = polynomial.polyvander2d(x, y, deg)
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f, rcond=None)[0]
    return c.reshape(deg+1)

def splint(spl, x):
    """
    Simple spline interpolator inspired by Numerical Recipes
    
    Credits unclear (E. Schlafly?)
    """
    npts = len(spl.x)
    lo = np.searchsorted(spl.x, x)-1
    lo = np.clip(lo, 0, npts-2)
    hi = lo + 1
    dx = spl.x[hi] - spl.x[lo]
    a = (spl.x[hi] - x)/dx
    b = (x-spl.x[lo])/dx
    y = (a*spl.y[lo]+b*spl.y[hi]+
         ((a**3-a)*spl.y2[lo]+(b**3-b)*spl.y2[hi])*dx**2./6.)
    return y

class CubicSpline:
    """
    Simple spline interpolator inspired by Numerical Recipes
    
    Credits unclear. (E. Schlafly?)
    """
    def __init__(self, x, y, yp=None):
        npts = len(x)
        mat = np.zeros((3, npts))
        # enforce continuity of 1st derivatives
        mat[1,1:-1] = (x[2:  ]-x[0:-2])/3.
        mat[2,0:-2] = (x[1:-1]-x[0:-2])/6.
        mat[0,2:  ] = (x[2:  ]-x[1:-1])/6.
        bb = np.zeros(npts)
        bb[1:-1] = ((y[2:  ]-y[1:-1])/(x[2:  ]-x[1:-1]) -
                    (y[1:-1]-y[0:-2])/(x[1:-1]-x[0:-2]))
        if yp is None: # natural cubic spline
            mat[1,0] = 1.
            mat[1,-1] = 1.
            bb[0] = 0.
            bb[-1] = 0.
        elif yp == '3d=0':
            mat[1, 0] = -1./(x[1]-x[0])
            mat[0, 1] =  1./(x[1]-x[0])
            mat[1,-1] =  1./(x[-2]-x[-1])
            mat[2,-2] = -1./(x[-2]-x[-1])
            bb[ 0] = 0.
            bb[-1] = 0.
        else:
            mat[1, 0] = -1./3.*(x[1]-x[0])
            mat[0, 1] = -1./6.*(x[1]-x[0])
            mat[2,-2] =  1./6.*(x[-1]-x[-2])
            mat[1,-1] =  1./3.*(x[-1]-x[-2])
            bb[ 0] = yp[0]-1.*(y[ 1]-y[ 0])/(x[ 1]-x[ 0])
            bb[-1] = yp[1]-1.*(y[-1]-y[-2])/(x[-1]-x[-2])
        y2 = solve_banded((1,1), mat, bb)
        self.x, self.y, self.y2 = (x, y, y2)
    def __call__(self, x):
        return splint(self, x)
    
    
    
    
