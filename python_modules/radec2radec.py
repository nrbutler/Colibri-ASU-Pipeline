#!/usr/bin/python3
"""
radec2radec.py fits_file radec_file
"""

import sys,os
from fit_wcs import ad2xy
from astropy.io.fits import getheader,getdata,writeto

def usage():
    print (__doc__)
    sys.exit()

def radec2radec(fits_file,radecfile,verify=True,fac=0.2):
    """
    """
    hdr=getheader(fits_file)
    distort = 'PV1_1' in hdr

    data = getdata(radecfile)
    ra,dec = data[:2]
    x,y = ad2xy(ra,dec,hdr,distort=distort)

    nx=hdr['NAXIS1']; ny=hdr['NAXIS2']
    h = (x>=1-fac*nx)*(x<=nx*(1+fac))*(y>=1-fac*ny)*(y<=ny*(1+fac))
    writeto(radecfile,data[:,h],overwrite=True)


def main():

    if (len(sys.argv)<3): usage()

    fits_file=sys.argv[1]
    radecfile=sys.argv[2]

    radec2radec(fits_file,radecfile)

if __name__ == "__main__":
    main()
