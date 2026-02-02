#!/usr/bin/python3
"""
 ps1wt2wt.py weight_image maskfile new_weight_image
"""
import os, sys
from astropy.io.fits import getdata,writeto
from scipy.ndimage import uniform_filter
from numpy import median

def usage():
    print (__doc__)
    sys.exit()

def ps1wt2wt(weightfile,maskfile,new_weightfile):
    """
      remove the source contribution from a ps1 weightmap
    """
    wt,hdr = getdata(weightfile,header=True)
    msk = getdata(maskfile)==0

    wt0 = 1.*wt
    for i in range(5):
        wt0[msk] = uniform_filter(wt0,25)[msk]

    if ('GAIN' not in hdr): hdr['GAIN'] = hdr['EXPTIME']

    a,b = wt.shape
    i0,i1 = a//2-100,a//2+101
    j0,j1 = b//2-100,b//2+101

    var0 = median( wt0[i0:i1,j0:j1] )
    if (var0>0): var0=1./var0
    else: var0=1.

    hdr['VAR0'] = var0
    writeto(new_weightfile,wt0,hdr,overwrite=True)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<4): usage()

    weightfile = sys.argv[1]
    if (os.path.exists(weightfile)==False): usage()

    maskfile = sys.argv[2]
    if (os.path.exists(maskfile)==False): usage()

    new_weightfile = sys.argv[3]
    ps1wt2wt(weightfile,maskfile,new_weightfile)
