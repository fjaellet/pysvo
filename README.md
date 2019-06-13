# pysvo - tools for dealing with the Theoretical Model Services of the Spanish Virtual Observatory

## Description

pysvo is a toolset that allows you to handle data coming from the Theoretical Model Services maintained by the Spanish Virtual Observatory. It was originally developed for use with StarHorse ([Queiroz et al. 2019](https://ui.adsabs.harvard.edu/abs/2018MNRAS.476.2556Q/abstract), [Anders et al. 2019](https://ui.adsabs.harvard.edu/abs/2019arXiv190411302A/abstract) ).

Currently, it can do the following things:
- Download and handle transmission curves supplied by the [Filter Profile Service](http://svo2.cab.inta-csic.es/theory/fps3/) of the Spanish Virtual Observatory (SVO)
- Handle some of the stellar model spectra available from the [Theoretical Spectra Webserver](http://svo2.cab.inta-csic.es/theory/newov2/index.php) of the SVO
- Handle some of the most widely used extinction curves
- Calculate extinction coefficients for arbitrary photometric bands on the fly

## Setup

You need to have the following python packages installed in order to use pysvo:
- numpy, scipy, astropy, matplotlib, collections

Optional:
- extinction

In order to work properly, currently you have to download the synthetic spectra [by hand](http://svo2.cab.inta-csic.es/theory//newov2/index.php?models=Kurucz) to the folder  ./pysvo/spectrallib/SUBFOLDERNAME/SPECTRUMFILE.txt
where SUBFOLDERNAME is either "Kurucz" or "BT-Settl", and the SPECTRUMFILE names depend on the library used...

This will be extended and automated soon.

## Tutorial 

Check out the tutorial [in this directory](https://github.com/fjaellet/pysvo/blob/master/pysvo_tutorial.ipynb).

## Credit
The code would be useless without the huge effort by the SVO group at CAB to create and maintain the Spanish VO services. Thus, if you make use of 'pysvo' as part of your work, please link to this repository and include the following acknowledgement in your publication:

    This research has made use of the SVO Filter Profile Service (http://svo2.cab.inta-csic.es/theory/fps/) supported from the Spanish MINECO through grant AYA2017-84089. 

If you use the Filter Profile Service, then the SVO team also appreciates if you could include the following references in your publication:

    The SVO Filter Profile Service. Rodrigo, C., Solano, E., Bayo, A. http://ivoa.net/documents/Notes/SVOFPS/index.html
    The Filter Profile Service Access Protocol. Rodrigo, C., Solano, E. http://ivoa.net/documents/Notes/SVOFPSDAL/index.html

Initially developed by Friedrich Anders (while at AIP). A few code snippets dealing with extinction curve of Schlafly et al. (2016) was written by Eddie Schlafly (MPIA).

