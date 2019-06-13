##################################################################################
#
#   pysvo.path: return the path of various files related to pysvo
#
##################################################################################
import numpy as np
import os.path

def photofolderPath():
    """
    Input:
       (none)
    Output:
       path string
    History:
       2017-03-08 - Written             - F. Anders (AIP)
       2019-06-10 - Adapted for pysvo   - F. Anders (ICCUB)
    """
    return './pysvo/photometry/'

def filtervotPath(telescope, photofilter, zeropoint):
    """
    Input:
       (none)
    Output:
       path string
    History:
       2015-05-05 - Written - Anders (AIP)
       2019-06-10 - Adapted for pysvo   - F. Anders (ICCUB)
    """
    return os.path.join(photofolderPath(), telescope + '_' + photofilter +
                               '_' + zeropoint + '.xml')

def girarditablePath():
    """
    Input:
       (none)
    Output:
       path string
    History:
       2015-05-05 - Written - Anders (AIP)
       2019-06-10 - Adapted for pysvo   - F. Anders (ICCUB)
    """
    return os.path.join(photofolderPath(), 'girardi_table_new.csv')

def spectrallib_path(lib="Kurucz"):
    """
    Path of the synthetic spectra.
    
    Optional arg:
        lib:               Spectral model name. Default: "Kurucz"
    Output:
        Path name
    History:
       2018-07-30 - Written                - F. Anders (AIP)
       2019-06-10 - Ported to pysvo.path - F. Anders (ICCUB)
    """
    return "./pysvo/spectrallib/" + lib + "/"

def spectrum_path(lib="Kurucz", mh=0.0, teff=6000., logg=5.0):
    """
    Path of a synthetic spectrum with given {[Fe/H], Teff, logg}
        Optional args:
        lib:               Spectral model name.   Default: "Kurucz"
        mh:                Metallicity.           Default: 0.0
        teff:              Effective temperature. Default: 6000.0
        logg:              Surface gravity.       Default: 5.0
    Output:
        Path name
    History:
       2018-07-30 - Written                  - F. Anders (AIP)
       2019-06-10 - Adapted for pysvo.path - F. Anders (ICCUB)
    """
    if lib=="Kurucz":
        if mh < 0:
            metstring = "m" + str(int(np.ceil(mh))) + str(int(10.0*mh % 10))
        else:
            metstring = "p" + str(int(np.floor(mh))) + str(int(10.0*mh % 10))
        return os.path.join( spectrallib_path(lib=lib), "f" + metstring + 
                             "k2odfnew.pck.teff=" + str(int(teff)) + "..logg="
                             + str(np.round(logg,1)) + "0000.dat.txt" ) 
    if lib=="BT-Settl":
        if teff<3000:
            pass
        else:
            raise ValueError("Out of temperature range.")
        return os.path.join( spectrallib_path(lib=lib), "lte0" + 
                             str(int(teff))[:2] + "-" + str(np.round(logg,1)) 
                             + "-" + str(np.round(mh,1)) 
                             + ".BT-Settl.7.dat.txt" )
    else:
        raise ValueError("No spectral library of name " + lib + " found")

