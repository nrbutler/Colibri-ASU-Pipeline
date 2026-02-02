#!/usr/bin/python3
"""
radec2xy.py fits_file radec_file [noverify]
"""

import sys,os
from fit_wcs import ad2xy
from astropy.io.fits import getheader,getdata
from numpy import loadtxt,arange

def usage():
    print (__doc__)
    sys.exit()

def radec2xy(fits_file,radecfile,verify=True):
    """
    """
    hdr=getheader(fits_file)
    distort = 'PV1_1' in hdr

    filename, file_extension = os.path.splitext(radecfile)
    if (file_extension=='.fits'): data = getdata(radecfile)
    else: data = loadtxt(radecfile,ndmin=2).T

    ra,dec = data[:2]
    if (len(data)>=4): mag,dmag = data[2:4]
    else: mag,dmag = 0*ra+19,0*dec+0.1

    if (len(data)>=5): color=data[4]
    else: color = 0.*ra

    x,y = ad2xy(ra,dec,hdr,distort=distort)
    ii=1+arange(len(x),dtype='int32')
    if (verify):
        nx=hdr['NAXIS1']; ny=hdr['NAXIS2']
        h = (x>=1)*(x<=nx)*(y>=1)*(y<=ny)
        x = x[h]; y = y[h]; mag = mag[h]; dmag=dmag[h]; ra=ra[h]; dec=dec[h]; color=color[h]; ii=ii[h]

    for i in range(len(x)): print ("""%f %f %d %f %f %f %f %f""" % (x[i],y[i],ii[i],mag[i],dmag[i],ra[i],dec[i],color[i]))


def main():

    if (len(sys.argv)<3): usage()

    fits_file=sys.argv[1]
    radecfile=sys.argv[2]

    verify=True
    if (len(sys.argv)>3): 
        if (sys.argv[3]=="noverify"): verify=False

    radec2xy(fits_file,radecfile,verify=verify)

if __name__ == "__main__":
    main()
