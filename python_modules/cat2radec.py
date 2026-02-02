#!/usr/bin/python3
"""
 cat2radec.py catalog.fits radecfile.txt
"""

import os,sys
from astropy.io.fits import getdata,getheader
from numpy import log10,log,savetxt,vstack

def usage():
    print (__doc__)
    sys.exit()

def cat2radec(catalog,radecfile):
    """ convert binary sextractor catalog file to txt radec format """    

    dat = getdata(catalog)
    hdr = getheader(catalog,0)

    try:
        dt,T0,T1,gain,cam = hdr['EXPTIME'], hdr['T0'], hdr['T1'], hdr['GAIN'], hdr['CCD']
    except:
        dt,T0,T1,gain,cam = 0,0,0,0,0

    hdr = f" RA DEC mag dmag mag_big dmag_big FWHM x y xa ya x2a y2a expos num (ccd {cam} gain {gain} exposure {dt} 0.0 0.0 sex_zero 25.0 t0 {T0} t1 {T1} )"
  
    xcol, ycol = dat.columns.names[0], dat.columns.names[1]
    x,y = dat[xcol], dat[ycol]
    Nx,Ny = x.max()-x.min(),y.max()-y.min()

    ra,dec = dat['ALPHA_J2000'], dat['DELTA_J2000']
    fwhm = dat['FWHM_IMAGE']
    f,f0 = dat['FLUX_APER'].T
    df,df0 = dat['FLUXERR_APER'].T
    if ( 'NUMBER' in dat.columns.names ): nn = dat['NUMBER']
    else: nn = dat['VECTOR_ASSOC']

    apc=[]
    apcfile=catalog+'.apweight'
    if (os.path.exists(apcfile)):
        wdata = getdata(apcfile)
        apc,apcb = wdata[:,2],wdata[:,3]
        good = apc>0; apcg=apc[good]
        f[good] /= apcg
        df[good] /= apcg
        f[~good] = 0
        df[~good] = 0
        good = apcb>0; apcg=apcb[good]
        f0[good] /= apcg
        df0[good] /= apcg
        f0[~good] = 0
        df0[~good] = 0

    m,m0,dm,dm0 = 0.*f,0.*f,0.*f,0.*f
    good = (f>0)*(df>0)*(f>0.2*df)
    goodb = (f0>0)*(df0>0)

    mcatalog = 'med_'+catalog
    if (os.path.exists(mcatalog)):
        fm = getdata(mcatalog)['FLUX_APER'][:,0]
        keep = fm>f/10
        gs = good.sum()
        good *= keep
        print (f" Dumping {gs-good.sum()} sources not in median stack")

    m[good],dm[good]= 25-2.5*log10(f[good]), 2.5/log(10)*df[good]/f[good]
    #m[good] -= dm[good]
    m0[goodb],dm0[goodb]= 25-2.5*log10(f0[goodb]), 2.5/log(10)*df0[goodb]/f0[goodb]
    #m0[good] -= dm0[good]

    xa,ya = x/Nx, y/Ny
    x2a,y2a = 0.5*(3*xa*xa-1),0.5*(3*ya*ya-1)
    fmt='%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %d'
    savetxt(radecfile, vstack((ra,dec,m,dm,m0,dm0,fwhm,x,y,xa,ya,x2a,y2a,0*x+dt,nn)).T,fmt=fmt,header=hdr)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<3): usage()

    catalog=sys.argv[1]
    if (os.path.exists(catalog)==False): usage()

    radecfile=sys.argv[2]

    cat2radec(catalog,radecfile)
