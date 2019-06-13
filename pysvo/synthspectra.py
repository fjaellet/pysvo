##################################################################################
#
#   pysvo.synthspectra: deal with SVO synthetic spectra
#
#   Main Functions:
#
#       - read_one_spec:      Open specified synthetic spectrum
#
#   Main classes:
#          - SynSpectrum:      Class for snthetic spectra
#               Funcs:  PlotSpectrum: Plots the spectrum
#                       Interpolate:  Interpolate the spectrum 
#                       GetFlux:      Get the flux value at a given wavelength 
#
##################################################################################
from matplotlib import pyplot as plt
import numpy as np
import os.path
from scipy.interpolate import interp1d
from pysvo import path, fitting_tools

def read_one_spectrum(filename):
    """
    Opens the specified SVO-like spectrum XML file, returns data recarray.
    
    Input:
        filename:          Filename of the ascii spectrum table.
    Output:
        tab:               Table without meta-data
    
    History:
        2018-07-30 - Written                - F. Anders (AIP)
        2019-06-10 - Ported to pysvo        - F. Anders (ICCUB)
    """
    tab = np.genfromtxt(filename, names=["lambda", "flux"])
    return tab

class SynSpectrum(object):
    """
    Class for synthetic spectra.
    Based on SVO synthetic spectra database:
    --- http://svo2.cab.inta-csic.es/theory/newov2/index.php ---
    
    History:
        2018-07-30 - Written                - F. Anders (AIP)
        2019-06-10 - Ported to pysvo        - F. Anders (ICCUB)
    """
    def __init__(self, lib="Kurucz", mh=0.0, teff=5750., logg=4.5):
        """
        Find out everything about a specified synthetic spectrum.
        
        Parameters
        ----------
        lib : string
            Which spectral models will be used. Default: "Kurucz"
        mh  : float
            Metallicity [dex]
        teff: float
            Effective temperature [K]
        logg: float
            Surface gravity [dex]

        """
        self.lib       = lib
        self.teff      = teff
        self.mh        = mh
        self.logg      = logg

        # Path where the information should be found:
        self.path = path.spectrum_path(lib=lib, mh=mh, teff=teff, logg=logg)
        
        # Check if the spectrum file already exists:
        if not os.path.isfile(self.path):
            # If not, download it from the SVO pages
            raise ValueError("""Please download the synthetic spectra by hand from
                                http://svo2.cab.inta-csic.es/theory/newov2/index.php 
                                before going ahead.""")
        # Read the Spectrum Table
        self.spectrum   = read_one_spectrum(self.path)
        self.wavelength = self.spectrum["lambda"]
        self.flux       = self.spectrum["flux"]

    def PlotSpectrum(self):
        """
        Plots the spectrum.
        """
        fig, ax = plt.subplots()
        plt.plot(self.wavelength, self.flux)
        ax.set_xlabel(r"$\lambda\quad [{\rm \AA}]$")
        ax.set_ylabel(r"Flux")

        
    def Interpolate(self, kind='cubic'):
        """
        Interpolate the spectrum either using scipy.interpolate (linear or cubic (slow!)), 
        or with the cubic splines module.
        """
        if kind in ['linear', 'cubic']:
            return interp1d(self.wavelength, self.flux, kind=kind)
        elif kind=='cubicspline':
            return fitting_tools.CubicSpline(self.wavelength, self.flux, yp='3d=0')
        else:
            raise ValueError("Interpolation type not found")
        
    def GetFlux(self, wavelengtharray, **kwargs):
        """
        Gets the flux for a given set of wavelengths
        """
        interpolation = self.Interpolate(**kwargs)
        return interpolation(wavelengtharray)
        

