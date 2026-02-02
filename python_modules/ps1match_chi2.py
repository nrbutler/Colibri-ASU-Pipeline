#!/usr/bin/python3
"""
 ps1match_chi2.py catfile ps1_catfile <ofile> 
"""
import os, sys
from astropy.io.fits import getdata,getheader,writeto
from numpy import log10,median

def usage():
    print (__doc__)
    sys.exit()


def ps1match_chi2(catfile,ps1_catfile,ofile='catalog_matches.fits',mag_diff_max=1.,snr_min=3.):
    """
      determine if sources in catfile are detected in ps1_catfile
    """
    cat = getdata(catfile)
    cat0 = getdata(ps1_catfile)

    c2,dc2 = cat['FLUX_APER'][:,0], cat['FLUXERR_APER'][:,0]
    c20,dc20 = cat0['FLUX_APER'][:,0],cat0['FLUXERR_APER'][:,0]

    # get an estimate of the calibration scatter
    j = (c20>100*dc20)*(c2+c20>0)
    delta_c2 = 2.5*log10( (c2[j]+c20[j])/c20[j] )
    calscatt = 1.48*median(abs(delta_c2-median(delta_c2)))

    fac = 10**(0.4*mag_diff_max)
    matched = ( c2 < (fac-1)*c20.clip(0) )*(c20>snr_min*dc20)*(dc2>0)

    detlim = 25-2.5*log10(snr_min*median(dc20))

    cat_id=1
    n0 = len(cat)
    if (os.path.exists(ofile)):
        os.system("""sethead BSCALE=1 %s""" % ofile)
        res = getdata(ofile).astype('int16')
        ohdr = getheader(ofile)
        cat_id += ohdr["NCAT"]
        cat_name="ps1_chi2"
        matched = matched * (res[1]==0)

    mag1 = 0.*c2
    mag1[matched] = 25 - 2.5*log10(c20[matched].clip(1))

    dis = 0.*c2
    res0 = (dis*10).round().astype('int16')
    res1 = (mag1*10).round().astype('int16')

    h = res1!=0
    res[0,h] = res0[h]
    res[1,h] = res1[h]
    res[2,h] = 10*cat_id

    if (cat_id>1): writeto(ofile, res,ohdr,overwrite=True)
    else: writeto(ofile, res,overwrite=True)

    os.system("""sethead BSCALE=0.1 NCAT=%d CATID%d=%s DETLIM=%.2f CALSCAT=%.4f %s""" % (cat_id,cat_id,cat_name,detlim,calscatt,ofile))


if __name__ == '__main__':
    """
    """
    if (len(sys.argv)<3): usage()

    catfile= sys.argv[1]
    if (os.path.exists(catfile)==False): usage()

    ps1_catfile= sys.argv[2]
    if (os.path.exists(ps1_catfile)==False): usage()

    ofile='catalog_matches.fits'
    if (len(sys.argv)>3): ofile=sys.argv[3]

    ps1match_chi2(catfile,ps1_catfile,ofile=ofile)
