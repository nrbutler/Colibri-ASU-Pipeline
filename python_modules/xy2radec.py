#!/usr/bin/python3
"""
xy2radec.py fits_file xy_file
"""

import sys,os
from fit_wcs import xy2ad
from astropy.io.fits import getheader,getdata
from numpy import loadtxt

def usage():
    print (__doc__)
    sys.exit()

def xy2radec(fits_file,xyfile):
    """
    """
    hdr=getheader(fits_file)
    distort = 'PV1_1' in hdr

    filename, file_extension = os.path.splitext(xyfile)
    if (file_extension=='.fits'): data = getdata(xyfile)
    else: data = loadtxt(xyfile,ndmin=2).T

    x,y = data[:2]
    r,d = xy2ad(x,y,hdr,distort=distort)

    for i in range(len(x)): print ("""%f %f""" % (r[i],d[i]))


def main():

    if (len(sys.argv)<3): usage()

    fits_file=sys.argv[1]
    xyfile=sys.argv[2]

    xy2radec(fits_file,xyfile)

if __name__ == "__main__":
    main()
