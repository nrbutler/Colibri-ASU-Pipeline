#!/usr/bin/python3
"""
 ps1match.py catfile ps1_catfile <ofile> <norm>
"""
import os, sys
from astropy.io.fits import getdata,getheader,writeto
from numpy import median, logical_or, savetxt, log10, zeros, sqrt
from quick_mode import quick_mode

def usage():
    print (__doc__)
    sys.exit()


def ps1match(catfile,ps1_catfile,ofile='catalog_matches.fits',mag_diff_max=1.,snr_min=2.,norm=-1):
    """
      determine if sources in catfile are detected in ps1_catfile
    """
    cat = getdata(catfile)
    flx, dflx = cat['FLUX_APER'][:,0], cat['FLUXERR_APER'][:,0]
    cat0 = getdata(ps1_catfile)
    flx0, dflx0 = cat0['FLUX_APER'][:,0], cat0['FLUXERR_APER'][:,0]

    if (norm<0):
        fmin=10**(-0.4*(20-25))
        j=(flx0>0)*(dflx0>0)*(flx>0)*(dflx0>0)*(flx0>fmin)
        norm = median(flx[j]/flx0[j])

    almost_matched = (flx>0)*(flx0>snr_min*dflx0)
    fac = 10**(0.4*mag_diff_max)

    cat_id=1
    n0 = len(cat)
    mag1 = 0.*flx
    maglim = 0.
    if (os.path.exists(ofile)):
        os.system("""sethead BSCALE=1 %s""" % ofile)
        res = getdata(ofile).astype('int16')
        ohdr = getheader(ofile)
        cat_id += ohdr["NCAT"]
        cat_name="ps1_stack"
        matched = almost_matched * ( flx<norm*flx0*fac ) * ( cat['VECTOR_ASSOC']<=0 )
        savetxt(catfile+'.ps1_matched.txt',matched,fmt='%d',header=' index')
    else:
        res = zeros((4,n0),dtype='int16')
        cat_name="ps1_subtraction"
        matched = almost_matched*( flx<norm*flx0*fac )
        # get 3-sigma limit
        h = (flx>0)*(flx<100*dflx)*(dflx0>0)
        maglim = quick_mode(25-2.5*log10(10*dflx0[h]))

    mag1[matched] = 25 - 2.5*log10(flx0[matched].clip(1)) - 2.5*log10(norm)

    dis = 0.38*sqrt( (cat['X_IMAGE']-cat0['X_IMAGE'])**2 + (cat['Y_IMAGE']-cat0['Y_IMAGE'])**2 )
    res0 = (dis*10).round().astype('int16')
    res1 = (mag1*10).round().astype('int16')

    h = res1!=0
    res[0,h] = res0[h]
    res[1,h] = res1[h]
    res[2,h] = 10*cat_id

    if (cat_id>1): writeto(ofile, res,ohdr,overwrite=True)
    else: writeto(ofile, res,overwrite=True)

    os.system("""sethead BSCALE=0.1 NCAT=%d CATID%d=%s %s""" % (cat_id,cat_id,cat_name,ofile))
    if (maglim>0): os.system(f"sethead MAGLIM={maglim:.2f} {ofile}")


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

    norm=-1
    if (len(sys.argv)>4): norm=float(sys.argv[4])

    ps1match(catfile,ps1_catfile,ofile=ofile,norm=norm)
