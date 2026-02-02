#!/usr/bin/python3
"""
grab_usno_local.py fits_file outfile ra dra dec ddec [year]
"""

import sys,os
from astropy.io.fits import getdata,writeto
from numpy import array,vstack,cos,zeros,searchsorted

def usage():
    print (__doc__)
    sys.exit()


def grab_usno_local(fitsfile,outfile,ra1,ra2,dec1,dec2,year=2018.):
    """
    """
    data = getdata(fitsfile).astype('float32')
    i0,i1 = searchsorted(data[0],array([ra1,ra2]))
    data = data[:,i0:i1]

    r,d,m = data[:3]

    prop_motion=False
    fitsfile1 = fitsfile.replace('.fits.gz','.pm.fits.gz')
    if (os.path.exists(fitsfile1)):
        prop_motion=True
        dr,dd,pm = 0.*r,0.*d,zeros(len(r),dtype='int16')
        sys.stderr.write("""Making proper motion corrections for %s.\n""" % fitsfile)
        ii,dx,dy = getdata(fitsfile1).astype('int32')
        j0,j1 = searchsorted(ii,array([i0,i1]))
        ii,dx,dy = ii[j0:j1]-i0,dx[j0:j1],dy[j0:j1]
        pm[ii]=1
        cd=cos(0.5*(dec1+dec2)/57.296)
        dr[ii] = dx*(year-2000.)/3.6e6
        dd[ii] = dy*(year-2000.)/3.6e6
        r[ii] += dr[ii]/cd
        d[ii] += dd[ii]

    h = (d>=dec1)*(d<=dec2)*(r>=ra1)*(r<=ra2)
    r,d,m = r[h],d[h],m[h]
    if (prop_motion): dr,dd,pm = dr[h],dd[h],pm[h]
    dm = 0*m+999
    h = m<=21
    dm[h] = 0.05 * 10**(0.11*(m[h]-10).clip(0))

    writeto(outfile, vstack((r,d,m,dm)) )
    if (prop_motion):
        outfile1 = outfile.replace('.fits','.pm.fits')
        h = pm==1
        if (h.sum()>0): writeto(outfile1, vstack((r[h],d[h],m[h],dm[h],dr[h],dd[h])) )


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

    year=2018.
    if (len(sys.argv)>7): year=float(sys.argv[7])

    if (ra-ra0>180): ra-=360
    elif (ra0-ra>180): ra+=360

    ra1,ra2 = ra-dra/2.,ra+dra/2.
    dec1,dec2 = dec-ddec/2.,dec+ddec/2.

    if (ra1<0): ra1,ra2 = 0.0,ra2-ra1
    if (ra2>360): ra1,ra2 = ra1-(ra2-360.0),360.0

    grab_usno_local(fits_file,outfile,ra1,ra2,dec1,dec2,year=year)


if __name__ == "__main__":
    main()
