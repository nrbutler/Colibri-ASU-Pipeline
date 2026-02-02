#!/usr/bin/python3
"""
grab_gaia_local.py fits_file outfile ra dra dec ddec [mjd]
"""

import sys,os
from astropy.io.fits import getdata,writeto
from numpy import array,vstack,sin,cos,arctan2,arcsin,zeros,pi,sqrt,searchsorted

def usage():
    print (__doc__)
    sys.exit()


def grab_gaia_local(fitsfile,outfile,ra1,ra2,dec1,dec2,mjd=57388.,ct=0.24):
    """
    """
    data = getdata(fitsfile).astype('float32')
    i0,i1 = searchsorted(data[0],array([ra1,ra2]))
    data = data[:,i0:i1]
    data = data[:,(data[1]>=dec1)*(data[1]<=dec2)]

    # proper motion
    fac = 180/pi
    xi = data[2]*((mjd-57388.)/365.25)/3.6e6/fac;
    eta = data[3]*((mjd-57388.)/365.25)/3.6e6/fac;

    seta0, ceta0 = sin(data[1]/fac), cos(data[1]/fac)
    cosc = 1./sqrt(1.+xi**2.+eta**2.)

    data[0] += fac*arctan2( xi, ceta0 - eta*seta0 )
    data[1] = fac*( arcsin(cosc*(seta0 + eta*ceta0)) )

    # ab corrections:
    data[4] += 0.1136
    data[6] -= 0.8406

    writeto(outfile, data[[0,1,4,5,6]])


def main():
    """
    """
    if (len(sys.argv)<7): usage()

    fits_file=sys.argv[1]
    if (os.path.exists(fits_file)==0): usage()

    outfile=sys.argv[2]

    file0 = os.path.basename(fits_file)
    fs = file0.split('_')
    ra0,dec0 = float(fs[2]), float(fs[3].strip('.fits.gz'))

    ra=float(sys.argv[3])
    dra=float(sys.argv[4])
    dec=float(sys.argv[5])
    ddec=float(sys.argv[6])

    mjd=57388.
    if (len(sys.argv)>7): mjd=float(sys.argv[7])

    if (ra-ra0>180): ra-=360
    elif (ra0-ra>180): ra+=360

    ra1,ra2 = ra-dra/2.,ra+dra/2.
    dec1,dec2 = dec-ddec/2.,dec+ddec/2.

    if (ra1<0): ra1,ra2 = 0.0,ra2-ra1
    if (ra2>360): ra1,ra2 = ra1-(ra2-360.0),360.0
  
    grab_gaia_local(fits_file,outfile,ra1,ra2,dec1,dec2,mjd=mjd)


if __name__ == "__main__":
    main()
