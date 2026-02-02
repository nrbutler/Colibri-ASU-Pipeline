#!/usr/bin/python3
"""
imxy2expos.py fitsfile_wtfile xy.txt
"""

import sys,os
from astropy.io.fits import getheader,getdata
from numpy import floor,loadtxt

def usage():
    print (__doc__)
    sys.exit()

def imxy2expos(x,y,idx,wtfile):
    """
    just work the exposure out using the weight map
    """

    dfile=wtfile.replace('.wt.fits','.fits')
    hdr=getheader(dfile)
    try: expos = hdr['EXPTIME']
    except: expos=1.0

    try: sx=hdr['NAXIS1']
    except: sx=1024.
    try: sy=hdr['NAXIS2']
    except: sy=1024.

    wt=getdata(wtfile)

    i,j = floor(y).astype('int16'),floor(x).astype('int16')
    expos_ar = wt[i,j]*expos/wt.max()

    x/=sx; y/=sy
    x2=0.5*(3*x**2-1)
    y2=0.5*(3*y**2-1)

    for k in range(len(x)):
        print ("""%.6f %.6f %.6f %.6f %.6f %d""" % (x[k],y[k],x2[k],y2[k],expos_ar[k],idx[k]))


def main():

    if (len(sys.argv)<3): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==0): usage()

    xyfile=sys.argv[2]
    if (os.path.exists(xyfile)==0): usage()

    x,y,idx = loadtxt(xyfile,unpack=True,usecols=(0,1,2),ndmin=2)
    idx=idx.astype('int32')

    imxy2expos(x,y,idx,infile)


if __name__ == "__main__":
    main()
