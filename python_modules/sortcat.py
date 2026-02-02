#!/usr/bin/python3
"""
 sortcat.py catalog.fits
"""

import os,sys
from astropy.io.fits import open
from numpy import arange,where

def usage():
    print (__doc__)
    sys.exit()

def sortcat(catalog,idfile=''):
    """ sort a catalog by flux """

    hdu = open(catalog)

    tofix = hdu[1].data['VECTOR_ASSOC']==0
    wtofix = where(tofix)[0]

    flx = hdu[1].data['FLUX_APER'][tofix,0]
    s = flx.argsort()[::-1]

    hdu[1].data['VECTOR_ASSOC'][wtofix[s]] = -1-arange(len(s),dtype='int32')
   
    hdu.writeto(catalog,overwrite=True)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<2): usage()

    catalog=sys.argv[1]
    if (os.path.exists(catalog)==False): usage()

    sortcat(catalog)
