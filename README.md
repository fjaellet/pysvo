# pysvo - Tools for dealing with stellar photometry, transmission curves, extinction coefficients, etc.

## Description

pysvo is a toolset that was developed for use with StarHorse ([Queiroz et al. 2019](https://ui.adsabs.harvard.edu/abs/2018MNRAS.476.2556Q/abstract), [Anders et al. 2019](https://ui.adsabs.harvard.edu/abs/2019arXiv190411302A/abstract) ).

Currently, it can do the following things:
- Download and handle transmission curves supplied by the [Filter Profile Service](http://svo2.cab.inta-csic.es/theory/fps3/) of the Spanish Virtual Observatory (SVO)
- Handle some of the stellar model spectra available from the [Theoretical Spectra Webserver](http://svo2.cab.inta-csic.es/theory/newov2/index.php) of the SVO
- Handle some of the most widely used extinction curves
- Calculate extinction coefficients for arbitrary photometric bands on the fly

## Dependencies

You need to have the following python packages installed in order to use pysvo:
- numpy, scipy, astropy, matplotlib, collections

## Tutorial 

Check out the tutorial [here](TBD)

## Credit
If you make use of pysvo as part of your work please cite 
[Anders et al. 2019](https://ui.adsabs.harvard.edu/abs/2019arXiv190411302A/abstract), 
and link to this repository.
