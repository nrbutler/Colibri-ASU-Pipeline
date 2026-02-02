#!/usr/bin/python3
"""
 sortcat.py catalog.fits psf_stars.fits
"""

import os,sys
from astropy.io.fits import open,getdata
from numpy import arange,where

def usage():
    print (__doc__)
    sys.exit()

def catpsfonly(catalog,psffile):
    """ get only psf star entries from a catalog """

    ii = getdata(psffile)
    hdu = open(catalog)
    hdu[1].data = hdu[1].data[ii]

    dir=os.path.dirname(catalog)
    file=os.path.basename(catalog)
    hdu.writeto(dir+'/psf_'+file,overwrite=True)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<3): usage()

    catalog=sys.argv[1]
    if (os.path.exists(catalog)==False): usage()

    psffile=sys.argv[2]
    if (os.path.exists(psffile)==False): usage()

    catpsfonly(catalog,psffile)
