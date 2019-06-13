##################################################################################
#
#   pysvo.photometry: tools to deal with various sources of photometry
#
#   Main Functions:
#
#       - read_one_filter:    Open filter curve + metadata of one specified
#                             photometric filter
#       - download_svo_filter:Downloads transmission curves + metadata
#       - read_girarditable:  Read the 'Available photometric systems' table
#                             for the Padova/PARSEC tracks
#       - PhotoDict:          For a girardi-style filter list, crossmatch to
#                             the SVO photometry database to get more details
#                             on the photometric filters...
#
#   Main classes:
#          - Photometry:      Class collecting photometric measurements.
#                             Initialise it by typing, e.g.,
#                             >>> photometry.Photometry('2mass+sloan')
#               Funcs:  get_filterinfo: Opens VO table file containing
#                                       information on the photometric
#                                       filter(s) specified.
#                                       If the file is not present, it will
#                                       be downloaded from the SVO website.
#          - PhotoFilter:     Class for photometric filters/measurements.
#                             Based on SVO-like photometry files
#               Funcs:  Transmissioncurve: Opens the transmission curve of the
#                                          filter.
#                       PlotTransmissioncurve: Plots the transmissioncurve of
#                                              the filter.
#                       getAlambda: Calculates the extinction in one filter
#                                   relative to the V-band, using a specified
#                                   extinction law.
#
###################################################################################

import astropy.io.votable as vot
from collections import OrderedDict
import os.path
import numpy as np
from scipy import interpolate
from numpy.polynomial import polynomial

from pysvo import path, extcurves, synthspectra, fitting_tools

def read_one_filter(filename, transmissioncurve=False):
    """
    Opens the specified SVO-like photometry file, returns data recarray.
    
    Input:
        filename:Filename of the VOTable.
    Optional:
        transmissioncurve: Boolean value. if True, download and load also
                           the transmission curve for the specified filters.
        zeropoint:         String. Either 'AB', 'Vega', or 'ST'.
    Output:
        votab:             Table with all the relevant (meta-)data
    """
    votab = vot.parse(filename)
    if transmissioncurve == True:
        # Return only the transmission curve as a np array
        transmission = votab.get_first_table()
        return transmission
    else:
        # Return everything as an astropy.votable Object
        return votab

def download_svo_filter(telescope, photofilter, zeropoint='AB'):
    """
    Downloads the transmission curve + metadata for the specified filter.

    Requires arguments:
        telescope:   SVO-like name of Telescope/Source of photometric system.
        photofilter: SVO-like name of photometric filter.
    Optional:
        zeropoint:   String. Either 'AB', 'Vega', or 'ST'.
        destination: Destination folder.
    Output:
        url:   URL of the relevant file.
    """
    import urllib.request, urllib.parse, urllib.error
    # Specify source and target 
    url = svo_filter_url(telescope, photofilter, zeropoint)
    target_file = path.filtervotPath(telescope, photofilter, zeropoint)
    # then do it.
    urllib.request.urlretrieve(url, target_file)
    return

def svo_filter_url(telescope, photofilter, zeropoint='AB'):  
    """
    Returns the URL where the filter transmission curve is hiding.

    Requires arguments:
        telescope:   SVO-like name of Telescope/Source of photometric system.
        photofilter: SVO-like name of photometric filter.
    Optional:
        zeropoint:   String. Either 'AB', 'Vega', or 'ST'.
    Output:
        url:   URL of the relevant file.
    """
    url = 'http://svo2.cab.inta-csic.es/theory/fps3/fps.php?' + \
          'PhotCalID=' + telescope + '/' + photofilter + '/' + zeropoint
    return url

def read_girarditable():
    """
    Opens the modified table of 'Available photometric systems' in the
    Padova/PARSEC tracks, originally downloaded from
    --- http://stev.oapd.inaf.it/~lgirardi/cmd_2.7/photsys.html ---,
    returns numpy recarray.

    Optional:
        path:    Path string, if not default directory.
    Output:
        table:   Table with all the relevant data
    """
    table = np.genfromtxt( path.girarditablePath(), comments='#', names=True,
                           delimiter=',', dtype=None ,encoding='utf8' )
    return  table 

class Photometry(object):
    """
    Class collecting photometric measurements.
    Initialise it by typing, e.g.,
        >>> photometry.Photometry('2mass+sloan')
    """
    def __init__(self, nicknames):
        """
        Collects and delivers information on the specified filter systems.
        
        Parameters
        ----------
        nicknames : string
            should contain the filter system names in the 'Girardi' style
            convention, separated by '+'
        """
        # Split the nicknames into an array of strings (quick and dirty)
        self.Nickname = np.array(nicknames.replace('_','+').split('+'))
        # Read the photometry table from the CMD website
        gir = read_girarditable()
        # Check if photometry nicknames exists in this table.
        for system in self.Nickname:
            if system not in gir['Nickname']:
                raise IOError("""The wanted Photometry nickname {} does not
                          exist in the listed nicknames.
                          You should either correct the spelling or
                          add your info in photometry/girardi_table_new.csv
                          """.format(system))

        # Crossmatch with girarditable:
        #  gi = np.array([g for g in gir if any( g['Nickname']==ph for ph in
        #                                          self.Nickname )] )
        gi=[]
        for name in self.Nickname:
            gi.extend(gir[ gir['Nickname'] == name ])
        gi = np.array( gi )
        
        # Create properties from this table:
        self.PhotSys      = gi['PhotSys']
        self.SVO_Nickname = gi['SVO_Nickname']
        self.MagSys       = gi['Zeropoint']
        self.Reference    = gi['Reference']
        self.Comments     = gi['Comments']
        
       
        # The property SVO_Names returns an array.
        self.SVO_Names = np.array([ i.tolist().split() for i in gi['SVO_Name']]).flatten()
        # The property FilterList gives you a plain list instead.
        self.FilterList = [ i.split() for i in gi['SVO_Name']]
        
    def get_filterinfo(zeropoints='girardi'):
        """
        Each of the photometric filters inherits the PhotoFilter properties. 

        Requires arguments:
            --
                       
        Optional args:
            zeropoints:        String. Currently restricted to 'girardi' -
                               for comparison with PARSEC/Padova tracks 
        Output:
            -- ( enables the properties of PhotoFilter for each filter )
        """
        if zeropoints == 'girardi':
            for Filter in self.FilterList:
                f = PhotoFilter(self.SVO_Nickname, Filter, self.MagSys)
                print(f.ZeroPoint)


    
class PhotoFilter(object):
    """
    Class for photometric filters/measurements.
    Based on SVO-like photometry files: see
    --- http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=voservice ---
    """
    def __init__(self, telescope, photofilter, zeropoint):
        """
        Find out everything about a specified photometric filter.
        
        Parameters
        ----------
        telescope : string
            Which telescope?
        photofilter : string
            Which filter?
        zeropoint : string
            Which Magnitude system? (Vega/AB/ST)
        """
        self.name      = photofilter
        self.telescope = telescope
        self.zeropoint = zeropoint

        # Path where the information should be found:
        self.path = path.filtervotPath(telescope, photofilter, zeropoint)
        
        # Check if the VOTable file already exists:
        if not os.path.isfile(self.path):
            # If not, download it from the SVO pages
            download_svo_filter(telescope, photofilter, zeropoint=zeropoint)
        # Read the VOTable
        votab = read_one_filter(self.path)
        self.votab = votab

        # Extract the metadata from the VOTable PARAMs:
        self.filterID = votab.get_field_by_id('filterID').value
        self.Description = votab.get_field_by_id('Description').value
        self.PhotSystem = votab.get_field_by_id('PhotSystem').value
        self.WavelengthMean = votab.get_field_by_id('WavelengthMean').value
        self.WavelengthEff = votab.get_field_by_id('WavelengthEff').value
        self.WavelengthMin = votab.get_field_by_id('WavelengthMin').value
        self.WavelengthMax = votab.get_field_by_id('WavelengthMax').value
        self.WidthEff = votab.get_field_by_id('WidthEff').value
        self.WavelengthCen = votab.get_field_by_id('WavelengthCen').value
        self.WavelengthPivot = votab.get_field_by_id('WavelengthPivot').value
        self.WavelengthPeak = votab.get_field_by_id('WavelengthPeak').value
        self.WavelengthPhot = votab.get_field_by_id('WavelengthPhot').value
        self.FWHM = votab.get_field_by_id('FWHM').value
        self.MagSys = votab.get_field_by_id('MagSys').value
        self.ZeroPoint = votab.get_field_by_id('ZeroPoint').value
        self.ZeroPointUnit = votab.get_field_by_id('ZeroPointUnit').value
        self.ZeroPointType = votab.get_field_by_id('ZeroPointType').value

    def References(self):
        """
        Opens the transmission curve of the filter.
        Input:
            self
        Output:
            transmission:  lambda T(lambda) numpy recarray
        """
        self.ProfileReference = self.votab.get_field_by_id('ProfileReference').value
        self.CalibrationReference = self.votab.get_field_by_id(
            'CalibrationReference').value
        return self.ProfileReference, self.CalibrationReference

    def Transmissioncurve(self):
        """
        Opens the transmission curve of the filter.
        Input:
            self
        Output:
            transmission:  lambda T(lambda) numpy recarray
        """
        return self.votab.get_first_table()

    def PlotTransmissioncurve(self):
        """
        Plots the transmissioncurve of the filter.
        """
        from matplotlib import pyplot as plt
        t = self.Transmissioncurve()
        fig, ax = plt.subplots()
        ax.set_xlabel(r"$\lambda\quad [{\rm \AA}]$")
        ax.set_ylabel(r"$T_{\lambda}$")
        plt.plot(t.array['Wavelength'], t.array['Transmission'])
        
    def getAlambda_rough(self, extlaw='schlafly', RV=3.32):
        """
        Calculates the approximate extinction in our filter relative to the V-band,
        using a specified extinction law.
        Input:
            self
        Output:
            Alambda:  A(filter) / AV = Integral ( A_lambda t(lambda) dlambda )
        """
        if hasattr(self, 'Afilter_rough'):
            pass
        else:
            # Get filter transmission curve (wavelengths in AA)
            t      = self.Transmissioncurve()
            lb, tm = t.array['Wavelength'], t.array['Transmission']
            # Get extinction law
            if extlaw.lower() == 'cardelli':
                extfunc = extcurves.Cardelli
            elif extlaw.lower() == 'schlafly':
                extfunc = extcurves.Schlafly2016
            else:
                raise ValueError("Did not find the requested extinction law.")

            # Calculate the filter extinction (numerical integration)
            ti      = 0.5 * ( tm[1:] + tm[:-1] ) 
            Ai      = extfunc( 0.5 * ( lb[1:] + lb[:-1] ), RV=RV )
            dli     = lb[1:] - lb[:-1] 
            self.Afilter_rough = np.sum( ti * Ai * dli ) / np.sum( ti * dli )

        return self.Afilter_rough

    def getAlambda_1spec(self, extlaw='schlafly', RV=3.32, lib="Kurucz", 
                        teff=5750., logg=4.5, mh=0.0, AV=1., 
                        kind="cubicspline"):
        """
        Calculates the extinction in our filter relative to the V-band,
        using a specified extinction law, a specified source spectrum, and a
        specified V-band extinction.
        Input:
            self
        Optional:
            extlaw:  String. Either "schlafly" (default) or "cardelli"
            RV:      R(V) parameter of the extinction curve
            lib:     synthetic spectral library. Default: "Kurucz"
            teff:    Teff [K]
            logg:    logg [dex]
            mh:      [M/H] [dex]
            AV:      A_V [mag]
            kind:    source spectrum interpolation type. Either "cubicspline", 
                      "linear", or "cubic" (very slow)

        Output:
                                        2.5                 Integral ( flux * transmission dlambda )
                      A(filter) / AV = ----- * log -------------------------------------------------------------
                                         AV         Integral ( flux * 10^(-0.4*Alambda) * transmission dlambda )

        """
        # Get filter transmission curve (wavelengths in AA)
        t      = self.Transmissioncurve()
        lb, tm = t.array['Wavelength'], t.array['Transmission']
        # Get extinction law
        if extlaw.lower() == 'cardelli':
            extfunc = extcurves.Cardelli
        elif extlaw.lower() == 'schlafly':
            extfunc = extcurves.Schlafly2016
        else:
            raise ValueError("Did not find the requested extinction law.")
        # Get source spectrum flux at the same wavelengths as the transmission curve
        source = synthspectra.SynSpectrum(lib=lib, teff=teff, logg=logg, mh=mh)
        sm = source.GetFlux(lb, kind=kind)
        # Calculate the filter extinction (numerical integration)
        ti      = 0.5 * ( tm[1:] + tm[:-1] ) 
        Si      = 0.5 * ( sm[1:] + sm[:-1] ) 
        Ai      = 10.**(-0.4 * AV * extfunc( 0.5 * ( lb[1:] + lb[:-1] ), RV=RV ) ) 
        #li      = 0.5 * ( lb[1:] + lb[:-1] )
        dli     = lb[1:] - lb[:-1] 
        Afilter = 2.5 * np.log10( np.sum( ti * Si * dli ) / \
                                  np.sum( ti * Si * Ai *dli )) / AV
        return Afilter

    def getAlambda_allspec(self, fit='poly2', **kwargs):
        """
        Calculates the extinction in our filter relative to the V-band,
        using a specified extinction law, a specified combination of stellar
        parameters, and a specified V-band extinction.

        This interpolates the results obtained by getAlambda_1spec for different
        stellar parameters.
        Input:
            self
        Optional:
            extlaw:  String. Either "schlafly" (default) or "cardelli"
            RV:      R(V) parameter of the extinction curve
            lib:     synthetic spectral library. Default: "Kurucz"
            kind:     source spectrum interpolation type. Either "cubicspline", 
                      "linear", or "cubic" (very slow)

        Output: Alambda / AV interpolator. You can now evaluate self.Alambda_2D(Teff, AV)
        
        History:
            2018-07-31 - Written                            - F. Anders (AIP)
            2019-06-10 - Adapted for pysvo                - F. Anders (ICCUB)
            2019-06-10 - Deleted deprecated fitting options - F. Anders (ICCUB)
        """
        print("Calculating extinction coefficient grid for filter ", self.name)
        # Define stellar parameter grid (Kurucz case):
        teff_vals = [ 3500., 3750., 4000., 4250., 4500., 5000., 5500., 6000., 6500., 
                      7000., 8000., 9000.,10000.,12000.,14000.,16000.,18000.,20000.,
                     25000.,30000.,35000.,40000.,45000.]
        #logg_vals = [ 5., 4., 3., 2., 1., 0.]
        #met_vals  = [-2.5, -1.5, -0.5, 0., 0.5]
        AV_vals   = [0.01, 0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.5, 15.]
        #teff_arr  = np.array(teff_vals, np.newaxis, np.newaxis)
        #logg_arr  = np.array(np.newaxis, logg_vals, np.newaxis)
        #met_arr   = np.array(np.newaxis, np.newaxis,  met_vals)
        Teff_arr, AV_arr    = np.zeros((len(teff_vals), len(AV_vals))), np.zeros((len(teff_vals), len(AV_vals)))
        Teff_arr[:,:] = np.array(teff_vals)[:,np.newaxis]
        AV_arr[:,:]   = np.array(AV_vals)[np.newaxis,:]
        Alambda_arr = -99. * np.ones((len(teff_vals), len(AV_vals)))
        for ii in np.arange(len(teff_vals)):
            for jj in np.arange(len(AV_vals)):
                # For the moment, forget about logg and M/H:
                Alambda_arr[ii, jj] = self.getAlambda_1spec(teff=teff_vals[ii],
                                                            logg=4.5, mh=0.0, AV=AV_vals[jj])
        if fit=='spline':
            self.Alambda_2D = interpolate.bisplrep(Teff_arr.flatten(), AV_arr.flatten(), Alambda_arr.flatten())
        elif 'poly' in fit:
            print("Fitting polynomial of order ", str(fit[-1]))
            self.Alambda_2D = fitting_tools.polyfit2d(Teff_arr.flatten(), AV_arr.flatten(), 
                                                      Alambda_arr.flatten(), 
                                                      [int(fit[-1]), int(fit[-1])])
        elif fit=='rough':
            self.Alambda_2D = self.getAlambda_rough()
        else:
            raise ValueError("Fitting option unknown.")
        return

    def evalAlambda_2D(self, Teff, AV, fit='poly2'):
        """
        Evaluate the extinction coefficient Alambda / AV for a combination
        of effective teperature and V-band extinction.
        
        Input:
            Teff, AV
        Output:
            Alambda/AV
        """
        if hasattr(self, 'Alambda_2D'):
            pass
        else:
            self.getAlambda_allspec(fit=fit)
        if fit=='spline':
            return interpolate.bisplev(Teff, AV, self.Alambda_2D)
        elif 'poly' in fit:
            return polynomial.polyval2d(Teff, AV, self.Alambda_2D)
        #elif fit=='poly2':
        #    return fitting_tools.polynomial2((Teff, AV), *self.Alambda_2D)
        #elif fit=='poly1':
        #    return fitting_tools.polynomial1((Teff, AV), *self.Alambda_2D)
        elif fit=='rough':
            return self.Alambda_2D 

def PhotoDict( filterlist ):
    """
    The magic function.
    
    For a girardi-style filter list, crossmatch to the SVO photometry
    database to get more details on the photometric filters...
    Input:
        filterlist: e.g., [ 'J_2mass', 'H_2mass', 'Ks_2mass' ]
    Output:
        dictionary connecting filterlist items to PhotoFilter objects
    """
    PhotoDict    = dict()
    filters = "+".join(list(OrderedDict.fromkeys(
                   [ii.split('_')[-1] for ii in filterlist] )))
    photolist    = Photometry(filters)
    kk = 0
    
    for ii in np.arange(len(photolist.Nickname)):
        sys_ii   = photolist.SVO_Nickname[ii]
        zero_ii  = photolist.MagSys[ii]
        print(photolist.SVO_Nickname, photolist.FilterList)
        for jj in np.arange(len(photolist.FilterList[ii])):
            filt_jj = photolist.FilterList[ii][jj]
            print(filterlist[kk].split('_')[-1], photolist.Nickname[ii])
            assert filterlist[kk].split('_')[-1] == photolist.Nickname[ii]
            if filt_jj != "-":
                PhotoDict[filterlist[kk]] = \
                            PhotoFilter(sys_ii, filt_jj,zero_ii)
            kk +=1
    return PhotoDict
    
