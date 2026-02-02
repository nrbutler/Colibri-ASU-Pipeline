#!/usr/bin/python3
"""
 quick_mode.py fits_table [field]
"""

import os,sys
from numpy import median
from astropy.io.fits import getdata

def usage():
    print (__doc__)
    sys.exit()

def quick_mode(x,frac=0.3,return_range=False):
    """
    """
    y = x.flatten()
    n0 = len(y)
    n=int(len(y)*frac)

    if (n0==0):
        return 0.
    elif (n==0):
        return median(y)
    else:
        y.sort()
        i0 = (y[n:] - y[:-n]).argmin()
        x0 = median( y[i0:i0+n+1] )
        if (return_range):
            return x0,y[i0],y[i0+n+1]
        else:
            return x0

if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<2): usage()

    file=sys.argv[1]
    if (os.path.exists(file)==False): usage()

    field='FWHM_IMAGE'
    if (len(sys.argv)>2): field=sys.argv[2]

    dat = getdata(file)[field]
    dat0=quick_mode(dat)
    h = (dat>dat0/2.)*(dat<dat0*4)
    print ("""FWHM=%.2f NSTARS=%d\n""" % (dat0,h.sum()))
