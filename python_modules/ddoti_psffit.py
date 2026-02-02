#!/usr/bin/python3
"""
ddoti_psffit.py fitsfile radec.txt [prior_file] [find_psf_stars]
   psf photometry for ddoti
"""

import sys,os
from astropy.io.fits import getheader,getdata,ImageHDU,HDUList,writeto
from numpy import abs,ones,zeros,median,arange,vstack,dstack,where,unique,hstack,sqrt
from scipy.spatial import cKDTree

def usage():
    print (__doc__)
    sys.exit()


def getcomp_frac(flx,flxb,frac=0.9):
    """
       find the flux level where flxb>flx for >=frac of sources
    """
    s=flx.argsort()[::-1]
    h = flxb>flx

    i0 = where( h[s].cumsum() < frac*arange(1,1+len(s)) )[0]
    if (len(i0)>0): i0=i0[0]
    else: i0=-1

    return min(1.e4,max(0.,flx[s[i0]]))


def twodbin(x,y,r,nbinx=20,nbiny=20,nmin=3):
    """
        create a 2d, median histogram
    """
    z = zeros((nbinx,nbiny),dtype='float64')-999

    dx=2./nbinx
    ix = ((x+1)/dx).astype('int16').clip(0,nbinx-1)
    dy=2./nbiny
    iy = ((y+1)/dy).astype('int16').clip(0,nbiny-1)

    for i in range(nbinx):
        kx = abs(ix-i)<=1
        if (kx.sum()>0):
            yk=iy[kx]; rk=r[kx]
            for j in range(nbiny):
                ky = abs(yk-j)<=1
                if (ky.sum()>=nmin): z[i,j] = median(rk[ky])

    if (nbinx>1):
        h = z[0]==-999; z[0,h] = z[1,h]
        h = z[-1]==-999; z[-1,h] = z[-2,h]
    if (nbiny>1):
        h = z[:,0]==-999; z[h,0] = z[h,1]
        h = z[:,-1]==-999; z[h,-1] = z[h,-2]
    z[z==-999]=0

    return z,ix,iy


def median_bin(x,y,z,nbinx=30,nbiny=30,mdl=[],ix=[],iy=[]):
    """
      pass the data through twodbin and linearly interpolate
    """
    n = len(x)
    mdl1=zeros(n,dtype='float64'); mdl2=zeros(n,dtype='float64')
    output=zeros(n,dtype='float64')

    # median bin the data
    dx,dy=2./nbinx,2./nbiny
    if (len(mdl)==0): mdl,ix,iy = twodbin(x,y,z,nbinx=nbinx,nbiny=nbiny)

    ix1 = (2*( (x+1)/dx-ix ).round()-1+ix).astype('int16').clip(0,nbinx-1)
    iy1 = (2*( (y+1)/dy-iy ).round()-1+iy).astype('int16').clip(0,nbiny-1)

    # interpolation of binned data in x
    h = ix1!=ix; ih=~h
    xh=(x[h]+1)/dx-0.5; ixh=ix[h]; iyh=iy[h]; ixh1=ix1[h];  iyh1=iy1[h]
    norm = ixh1-ixh
    mdl1[h] = ( (ixh1-xh)*mdl[ixh,iyh ] + (xh-ixh)*mdl[ixh1,iyh ] )/norm
    mdl2[h] = ( (ixh1-xh)*mdl[ixh,iyh1] + (xh-ixh)*mdl[ixh1,iyh1] )/norm
    ixth=ix[ih]
    mdl1[ih] = mdl[ixth,iy[ih]]; mdl2[ih] = mdl[ixth,iy1[ih]] 

    # interpolation of binned data in y
    h = iy1!=iy; ih=~h
    yh=(y[h]+1)/dy-0.5; iyh=iy[h]; iyh1 = iy1[h]
    output[h] = ( (iyh1-yh)*mdl1[h] + (yh-iyh)*mdl2[h] )/(iyh1-iyh)
    output[ih] = mdl1[ih]

    return output


def ddoti_psffit(fitsfile,radecfile,prior_file='',bin_fac=12/4096.,nmode=4,snr_min_psf=100.,snr_max_psf=10000.,find_psf_stars=False):
    """
      PSF-like fitting routine, where we consider a central aperture and an annulus.
      The ratio of flux in big aperture (flxb) to the central aperture (flx) is fit
        using a 2d polynomial.  result is saved to file.

      Define: ( flx + r1*(flxb-flx) )*(1+r1)/(1+r1**2) to be the PSF estimate of flxb

      If a prior file is supplied (e.g., from the stack), then an aperture correction (ap_corr)
        is calculated to related flxb to ap_corr*flxb0, with flxb0 from the prior file.

      The fit results are then used (by calibrate.py) to determine the big aperture flux using the central aperture only.
    """
    nmode = max(2,2*(nmode//2)) # must be even, >=2

    # read in the sextractor output
    phot_dat=getdata(radecfile)

    hdr = getheader(fitsfile)
    a,b = hdr['NAXIS1'],hdr['NAXIS2']
    sat_level=hdr['SATURATE']

    nbinx,nbiny = int(a*bin_fac), int(b*bin_fac)

    rdir=os.path.dirname(radecfile)
    fits_file_large=rdir+'/'+fitsfile
    x0,y0 = 0.,0.
    if (os.path.exists(fits_file_large)):
        hdr0=getheader(fits_file_large)
        x0 = hdr0['CRPIX1']-hdr['CRPIX1']
        y0 = hdr0['CRPIX2']-hdr['CRPIX2']

    flx,flxb=phot_dat['FLUX_APER'][:,0],phot_dat['FLUX_APER'][:,1]
    dflx,dflxb=phot_dat['FLUXERR_APER'][:,0],phot_dat['FLUXERR_APER'][:,1]
    x=(phot_dat['X_IMAGE']-x0-a//2)*2/a; y=(phot_dat['Y_IMAGE']-y0-b//2)*2/b

    good=(flx<10*sat_level)*(flx>0)*(flxb>0)
    fmin = median(flx[good][flxb[good]>10*dflxb[good]])
    good *= (flx>fmin)*(flx<10*sat_level)
    gs = good.sum()

    if (gs==0):
       sys.stderr.write("""No good data found, exiting...\n""")
       sys.exit()

    # use a tree classifier to get nmed=3*nmode nearest (good) neighbors to each source
    xy=vstack((x,y)).T
    tre = cKDTree(xy[good])
    if (3*nmode>gs): nmode=gs//3
    ii = tre.query(xy,k=3*nmode)[1]

    n=len(x)
    ii0 = arange(n)
    ap_corr=ones(n,dtype='float64')
    ap_corrb=ones(n,dtype='float64')
    mdl_ap_corr = ones((nbinx,nbiny),dtype='float32')
    mdl_ap_corrb = ones((nbinx,nbiny),dtype='float32')

    if (gs<nbinx*nbiny/2):
        sys.stderr.write(""" Not enough sources (%d), psf fitting skipped for %s\n""" % (gs,fitsfile))
    else:

        if (find_psf_stars):
            # estimate the mode of r1 at each point
            r10 = flxb[good]/flx[good]-1
            r1 = r10[ii]
            ir1 = r1.argsort(axis=1)
            r1.sort(axis=1)
            i0=(r1[:,nmode:]-r1[:,:-nmode]).argmin(axis=1) + nmode//2

            if (nbinx>1 and nbiny>1):
                mdl_r1,ix,iy = twodbin(x,y,r1[ii0,i0],nbinx=nbinx,nbiny=nbiny)
                r1 = median_bin(x,y,[],nbinx=nbinx,nbiny=nbiny,mdl=mdl_r1,ix=ix,iy=iy)
            else:
                r1 = r1[ii0,i0]

            # what are the matches for the good stars?
            wg = where(good)[0]
            dr1 = r1[good] - r10[ii[wg,ir1[wg,i0[wg]]]]
            dr10 = median(dr1); dr1 = abs(dr1-dr10)
            h = dr1>0
            if (h.sum()>0): delta_dr1 = 1.48*median(dr1[h])
            else: delta_dr1 = 1.
            good_psf = good.copy()
            h = (dr1<2*delta_dr1)*(flx[good]>snr_min_psf*dflx[good])*(flx[good]<snr_max_psf*dflx[good])
            if (h.sum()>1): good_psf[wg] = h

            sel_indx = wg[unique(ii[ii0,ir1[ii0,i0]])]
            sel_indx = sel_indx[good_psf[sel_indx]]
            sfile = fitsfile.replace('.fits','.psf_stars.fits')
            writeto(sfile,sel_indx.astype('int32'),overwrite=True)
            os.system("""sethead NX0=%d NY0=%d NBINX=%d NBINY=%d %s""" % (a,b,nbinx,nbiny,sfile))

        elif (len(prior_file)>0):
            if (os.path.exists(prior_file)):
                sys.stderr.write("""Using prior file %s\n""" % prior_file)
                pdat = getdata(prior_file,0)
                flx0 = pdat[:,0].clip(0.1)
                flxb0 = pdat[:,1].clip(0.1)
                ap0 = pdat[:,2]
                apb0 = pdat[:,2]
                #try: mdl_ap0 = getdata(prior_file,1)[1]
                #except:
                #   sys.stderr.write(" No prior aperture correction data\n")
                #   mdl_ap0 = 1.
                good0 = good*(flxb0>1)
                gs0 = good0.sum()
                if (gs0<gs):
                    good = good0
                    tre = cKDTree(xy[good])
                    if (3*nmode>gs0): nmode=gs0//3
                    ii = tre.query(xy,k=3*nmode)[1]

                # estimate the mode of ap_corr at each point
                ap_corr = (flx[good]/flx0[good])[ii]
                ap_corr.sort(axis=1)
                i0=(ap_corr[:,nmode:]-ap_corr[:,:-nmode]).argmin(axis=1) + nmode//2

                if (nbinx>1 and nbiny>1):
                    mdl_ap_corr,ix,iy = twodbin(x,y,ap_corr[ii0,i0],nbinx=nbinx,nbiny=nbiny)
                    ap_corr = apb0*median_bin(x,y,[],nbinx=nbinx,nbiny=nbiny,mdl=mdl_ap_corr,ix=ix,iy=iy).clip(0)
                else: 
                    ap_corr = apb0*ap_corr[ii0,i0]

                ap_corr /= median(ap_corr)
                mdl_ap_corr = twodbin(x,y,ap_corr,nbinx=nbinx,nbiny=nbiny)[0]

                ap_corrb = (flxb[good]/flxb0[good])[ii]
                ap_corrb.sort(axis=1)
                i0=(ap_corrb[:,nmode:]-ap_corrb[:,:-nmode]).argmin(axis=1) + nmode//2

                if (nbinx>1 and nbiny>1):
                    mdl_ap_corrb = twodbin(x,y,ap_corrb[ii0,i0],nbinx=nbinx,nbiny=nbiny)[0]
                    ap_corrb = ap0*median_bin(x,y,[],nbinx=nbinx,nbiny=nbiny,mdl=mdl_ap_corrb,ix=ix,iy=iy).clip(0)
                else:
                    ap_corrb = ap0*ap_corrb[ii0,i0]

                ap_corrb /= median(ap_corrb)
                mdl_ap_corrb = twodbin(x,y,ap_corrb,nbinx=nbinx,nbiny=nbiny)[0]

    outfile=radecfile+'.apweight'
    hdu0 = ImageHDU(); hdu0.data = vstack((flx,flxb,ap_corr,ap_corrb)).T
    hdu1 = ImageHDU(); hdu1.data = dstack((mdl_ap_corr,mdl_ap_corrb)).T
    hdu = HDUList(); hdu.append(hdu0); hdu.append(hdu1)
    hdu.writeto(outfile,overwrite=True)


if __name__ == "__main__":

    if (len(sys.argv)<3): usage()

    fitsfile=sys.argv[1]
    if (os.path.exists(fitsfile)==0): usage()
    radecfile=sys.argv[2]
    if (os.path.exists(radecfile)==0): usage()

    prior_file=''
    if (len(sys.argv)>3): prior_file=sys.argv[3]

    hdr = getheader(fitsfile) 
    ps=sqrt(hdr['CD1_1']**2+hdr['CD1_2']**2)*3600.
    bin_fac=12/4096.
    #bin_fac=0.5*ps*25/6144

    find_psf_stars=False
    if (len(sys.argv)>4): 
        if (sys.argv[4]=="find_psf_stars"):
            find_psf_stars=True

    ddoti_psffit(fitsfile,radecfile,prior_file=prior_file,bin_fac=bin_fac,find_psf_stars=find_psf_stars)
