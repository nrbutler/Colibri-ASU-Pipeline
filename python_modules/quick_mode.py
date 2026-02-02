#!/usr/bin/python3
"""
 quick_mode.py image
"""

import os,sys
from astropy.io.fits import getdata
from numpy import median

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

    img=sys.argv[1]
    if (os.path.exists(img)==False): usage()
    i0=1100
    if (len(sys.argv)>2): i0=int(sys.argv[2])

    dat = getdata(img)
    x0 = quick_mode(dat[:i0])
    x1 = 1.
    if (i0<len(dat)): x1 = quick_mode(dat[i0:])

    print (x0,x1)
