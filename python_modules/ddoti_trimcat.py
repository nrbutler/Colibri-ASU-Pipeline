#!/usr/bin/python3
"""
ddoti_trimcat.py stack_file phot_file [snr_min] [snr_ref]
"""

import sys,os
from numpy import ones,where,arange,vstack,unique,logical_or,log10,sqrt,savetxt
from astropy.io.fits import getdata,getheader,writeto
from fit_wcs import ad2ij
from scipy.spatial import cKDTree

def usage():
    print (__doc__)
    sys.exit()


def ddoti_trimcat(stack_file,phot_file,snr_min=2.,snr_ref=500.,edge=150,cluster=True):
    """
      only retain catalog sources with >snr_min, and
        dump sources not present in the median stack
        dump clustered detections near bright stars
      note: dumping == id=0
    """
    hdr = getheader(stack_file)
    nx,ny = hdr['NAXIS1'],hdr['NAXIS2']

    # platescale
    ps=sqrt(hdr['CD2_1']**2+hdr['CD2_2']**2)*3600.

    phot_dat = getdata(phot_file)
    if (len(phot_dat)==0): sys.exit()

    id=phot_dat['VECTOR_ASSOC'].astype('int32')

    flx,dflx = phot_dat['FLUX_APER'][:,0],phot_dat['FLUXERR_APER'][:,0]

    # reference catalogued,bright stars to prunce around
    ref = (flx>snr_ref*dflx)*(id>0)
    source_rad = (sqrt(phot_dat['ISOAREAF_IMAGE'][ref])*ps).clip(1.)
    r0,d0 = phot_dat['ALPHA_J2000'][ref],phot_dat['DELTA_J2000'][ref]
    i0,j0 = ad2ij(r0,d0,hdr)

    # new sources
    #
    # ignore low snr new sources
    bad0 = logical_or(flx<=snr_min*dflx,phot_dat['FLAGS']==1)
    id[bad0] = 0
  
    good0 = (~bad0)*(id<=0)
    wg0 = where(good0)[0]

    # give new sources negative indices
    ii = flx[good0].argsort()[::-1]
    id[wg0[ii]] = -1-arange(len(wg0))

    try:
        faint_data = getdata('catalog_matches.fits')
        fhdr = getheader('catalog_matches.fits')
        good = good0*(faint_data[1]<=0)
        have_faint = True
    except:
        print ("unable to open faint source list")
        good = good0
        have_faint = False

    sys.stderr.write("""Starting with %d\n""" % good.sum())

    ng = good.sum()
    if (ng>0):
        wg = where(good)[0]

        #
        # potentially drop sources not in median stack
        stack_dir=stack_file.replace('.fits','_dir')
        stack_cat=stack_dir+'/catalog.fits'
        mstack_dir='med_'+stack_dir
        mstack_cat=mstack_dir+'/catalog.fits'
        if (os.path.exists(mstack_cat)):
            flx0 = getdata(stack_cat)['FLUX_APER'][good,0].clip(1.)
            dflx0 = getdata(stack_cat)['FLUXERR_APER'][good,0].clip(1.)
            # snr_fac is 2 for snr>=0, rising to low snr, and ~10 for snr=5
            # snr_fac = ( flx0/(flx0-5.*dflx0).clip(0.1*flx0) ).clip(2.)
            # more agressive, but safe for 5<snr<10:
            snr_fac = 2.
            flx1 = getdata(mstack_cat)['FLUX_APER'][good,0]
            not_in_median = flx1*snr_fac<flx0
            id[wg[not_in_median]]=0
            sys.stderr.write("""Dumping %d sources not in the median stack\n""" % not_in_median.sum())
            good[good] *= ~not_in_median

    #
    # now, potentially drop new sources around bright stars
    ng = good.sum()

    if (ng>1 and cluster):
        wg = where(good)[0]

        snr_good = flx[good]/dflx[good]
        r1,d1 = phot_dat['ALPHA_J2000'][good],phot_dat['DELTA_J2000'][good]
        i1,j1 = ad2ij(r1,d1,hdr)

        tre = cKDTree( vstack((i1,j1)).T )

        # find the 100 closest, new sources to each catalogued source
        dis,ii = tre.query( vstack((i0,j0)).T,k=min(100,ng) )
        dis *= 3600

        # 
        maybe_bad = ( dis.T < source_rad/(1+(snr_good/500.)[ii].T) ).T
        mb1 = maybe_bad.sum(axis=1)>2
        dis=dis[mb1]; ii=ii[mb1]; source_rad=source_rad[mb1]

        bad = ( dis.T < source_rad/(1+(snr_good/500.)[ii].T) ).T
        iib0 = unique(ii[bad])
        iib = wg[iib0]
        sys.stderr.write("""Dumping %d clustered sources near bright stars\n""" % len(iib))

        id[iib] = 0

        #savetxt('trimmed.txt',vstack((r1[iib0],d1[iib0])).T)

    # finally remove sources near the edge
    n0 = (id==0).sum()
    x,y = phot_dat['X_IMAGE'],phot_dat['Y_IMAGE']
    id[x<=edge]=0; id[y<=edge]=0
    id[x>=nx-edge]=0; id[y>=ny-edge]=0
    n1 = (id==0).sum()
    sys.stderr.write("""Dumping %d sources near the image edge (%d %d %d)\n""" % (n1-n0,edge,nx,ny))

    phot_dat['VECTOR_ASSOC']=id
    writeto(phot_file,phot_dat,overwrite=True)


if __name__ == "__main__":
    """
    """

    if (len(sys.argv)<3): usage()

    stack_file=sys.argv[1]
    if (os.path.exists(stack_file)==0): usage()

    phot_file=sys.argv[2]
    if (os.path.exists(phot_file)==0): usage()

    snr_min=2.
    if (len(sys.argv)>3): snr_min = float(sys.argv[3])

    edge=150
    if (len(sys.argv)>4): edge = int(sys.argv[4])

    cluster=True
    if (len(sys.argv)>5):
        if (sys.argv[5]=="AUTO"): cluster=False

    ddoti_trimcat(stack_file,phot_file,snr_min=snr_min,edge=edge,cluster=cluster)
