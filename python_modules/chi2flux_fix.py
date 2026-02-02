#!/usr/bin/python3
"""
 chi2flux_fix.py infile
"""
import os, sys
from astropy.io.fits import getdata,writeto
from numpy import sqrt

def usage():
    print (__doc__)
    sys.exit()


def chi2flux_fix(in_file,im_min=1.e-3):
    """
      scale chi2 fluxes near threshold
    """
    dat, hdr = getdata(in_file,header=True)

    if ("SCALED" not in hdr):

        flx, dflx = 1.*dat['FLUX_APER'][:,0], 1.*dat['FLUXERR_APER'][:,0]

        h = (flx>0)*(dflx>0)
        snr = (flx[h]/dflx[h]).clip(3)
        corr_fac = 10/snr*sqrt( (1+0.1*snr)**2-1 )
        dat['FLUX_APER'][h,0] *= corr_fac
        dat['FLUXERR_APER'][h,0] *= corr_fac

        hdr["SCALED"]=1.0
        writeto(in_file,dat,hdr,overwrite=True)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<2): usage()

    in_file= sys.argv[1]
    if (os.path.exists(in_file)==False): usage()

    chi2flux_fix(in_file)
