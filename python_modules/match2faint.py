#!/usr/bin/python3
"""
match2faint.py fitsfile catalog_file faintfile [snr_min] [nmatch]
"""
import sys,os
from numpy import vstack,log10,arange,zeros,sqrt,where,abs
from scipy.spatial import cKDTree
from astropy.io.fits import getdata,getheader,writeto
from fit_wcs import ad2ij

def usage():
    print (__doc__)
    sys.exit()


def match2faint(fitsfile,catalog_file,faintfile,snr_min=2.,mag_diff_max=2.,nmatch=1,ofile='catalog_matches.fits'):
    """
       allow a faint source to match to as many as nmatch
       results saved to ofile
    """
    hdr=getheader(fitsfile)
    ps=sqrt(hdr['CD2_1']**2+hdr['CD2_2']**2)
    fwhm0 = hdr['FWHM']*ps

    cdata = getdata(catalog_file)
    n0 = len(cdata)
    id = cdata['VECTOR_ASSOC']

    cat_id=1
    if (os.path.exists(ofile)):
        os.system("""sethead BSCALE=1 %s""" % ofile)
        res = getdata(ofile).astype('int16')
        ohdr = getheader(ofile)
        cat_id += ohdr["NCAT"]
    else: res = zeros((4,n0),dtype='int16')

    flx = cdata['FLUX_APER'][:,0]
    dflx = cdata['FLUXERR_APER'][:,0]
    good = (flx>snr_min*dflx)*(id<=0)*(res[1]<=0)

    #first time through we should check mags of id>0

    if (good.sum()==0):
        sys.stderr.write("""No sources to match to %s\n""" % faintfile)
    else:

        cdata = cdata[good]
        mag = ( 25-2.5*log10(flx[good]) ).astype('float64')

        fwhm = cdata['FWHM_IMAGE'].clip(fwhm0/ps)*ps
        i,j = ad2ij(cdata['ALPHA_J2000'],cdata['DELTA_J2000'],hdr)

        # read in faint data
        ra,dec,m0,dm = getdata(faintfile).astype('float64')
        i0,j0 = ad2ij(ra,dec,hdr)

        # want to do the same thing, but building the tree on i,j
        # find the nmatch closest ddoti source to each faint source
        tre = cKDTree( vstack((i,j)).T )
        dis,ii = tre.query( vstack((i0,j0)).T,k=nmatch )

        # matched sources
        sel = dis<=fwhm[ii]
        if (nmatch==1): m = 1.*m0
        else:
            m = zeros((len(i0),nmatch),dtype='float64')
            for i in xrange(nmatch): m[:,i] = m0

        # associate each ddoti sources with the faint source of nearest mag
        mag1,dis1 = 0.*mag,0.*mag
        ii = ii[sel]; dis = dis[sel]; m = m[sel]
        iis = abs(mag[ii]-m).argsort()[::-1]
        ii = ii[iis]
        mag1[ii] = m[iis]
        dis1[ii] = dis[iis]

        # flag matches with out of range magnitudes for later (used in ddoti_lc_plot.py)
        # f_fac takes into account psf-mag error wrt fwhm
        f_fac = 1.25*log10((fwhm/(2*ps)).clip(1.))
        kind_of_bad = mag1>=mag+mag_diff_max+f_fac
        #mag1[kind_of_bad] *= -1
        mag1[kind_of_bad] = 0

        sys.stderr.write("""Matched %d (%d possible) previously-uncatalogued sources to %s\n""" % (len(ii),kind_of_bad.sum(),faintfile))

        res0 = (dis1*10*3600).round().astype('int16')
        res1 = (mag1*10).round().astype('int16')
        wg=where(good)[0]

        h = res1!=0
        res[0,wg[h]] = res0[h]
        res[1,wg[h]] = res1[h]
        res[2,wg[h]] = 10*cat_id

        cat_name = faintfile.replace('_radec.fits','')

        if (cat_id>1): writeto(ofile, res,ohdr,overwrite=True)
        else: writeto(ofile, res,overwrite=True)

        os.system("""sethead BSCALE=0.1 NCAT=%d CATID%d=%s %s""" % (cat_id,cat_id,cat_name,ofile))


if __name__ == "__main__":

    if (len(sys.argv)<4): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==0): usage()
    catfile=sys.argv[2]
    if (os.path.exists(catfile)==0): usage()
    faintfile=sys.argv[3]
    if (os.path.exists(faintfile)==0): usage()

    snr_min=2.
    if (len(sys.argv)>4): snr_min=float(sys.argv[4])

    mag_diff_max=2.
    if (len(sys.argv)>5): mag_diff_max=float(sys.argv[5])

    nmatch=1
    if (len(sys.argv)>6): nmatch=int(sys.argv[6])

    match2faint(infile,catfile,faintfile,mag_diff_max=mag_diff_max,nmatch=nmatch,snr_min=snr_min)
