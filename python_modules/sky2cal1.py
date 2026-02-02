#!/usr/bin/python3
"""
sky2cal.py stack.fits catalog.fits sky.list calfile
"""

import sys,os
from fit_wcs import ad2ij,ij2ad
from numpy import loadtxt,ones,vstack,cos,median,sqrt,where,dstack,arange,meshgrid,array,dot,log10,zeros
from numpy.linalg import lstsq

from astropy.io.fits import open,getdata,getheader,writeto
from scipy.spatial import cKDTree

from ddoti_psffit import median_bin,twodbin

import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import subplots


def usage():
    print (__doc__)
    sys.exit()


def dist_corr(dx_img,dy_img,hdr):
    """
    """
    a,b = hdr['NAXIS1'],hdr['NAXIS2']
    nx,ny = dx_img.shape

    xx = a*( arange(nx)*1./nx + 1./nx/2. )
    yy = b*( arange(ny)*1./ny + 1./ny/2. )
    x0,dx0 = hdr['CRPIX1'], hdr['CD1_1']
    y0,dy0 = hdr['CRPIX2'], hdr['CD2_2']
    i,j = (xx-x0)*dx0, (yy-y0)*dy0

    xx,yy = meshgrid(i, j, copy=False, indexing='ij')
    A = array([xx, yy, xx**2, xx*yy, yy**2])

    cxy = lstsq(A.reshape((5,nx*ny)).T, array([dx_img,dy_img]).reshape(2,nx*ny).T,rcond=None)[0]

    return dot(A.T,cxy).T,cxy


def sky2cal(stack_file,fits_file,sky_file,cal_file,mfile='catalog_matches.fits',nstars=10):
    """
    """
    hdr=getheader(stack_file)
    cfilter=hdr['CALFILT']

    ps=sqrt(hdr['CD2_1']**2+hdr['CD2_2']**2)
    a,b = hdr['NAXIS1'],hdr['NAXIS2']

    hdu = open(fits_file)
    id = hdu[1].data['VECTOR_ASSOC'].astype('int32')

    x=(hdu[1].data['X_IMAGE']-a//2)*2/a; y=(hdu[1].data['Y_IMAGE']-b//2)*2/b

    # we will track matches in a separate file (mfile), with array res
    n0 = len(x)
    res = zeros((4,n0),dtype='int16')

    sdata = getdata(sky_file)
    ra,dec,mag,dmag = sdata[:4]
    if (len(sdata)>4): color = sdata[4]
    else: color = 0.*ra
    i,j = ad2ij(ra,dec,hdr)

    #
    # first we tune up the ra0,dec0 positions
    #
    ra0,dec0 = hdu[1].data['ALPHA_J2000'], hdu[1].data['DELTA_J2000']
    i0,j0 = ad2ij(ra0,dec0,hdr)

    h = id>0
    idh = id[h]-1
    xy=vstack((i[idh],j[idh])).T
    tre=cKDTree(xy)
    xy0=vstack((i0,j0)).T
    ii = tre.query(xy0,k=21)[1]

    dx, dy = i[idh]-i0[h], j[idh]-j0[h]
    dx0 = median(dx[ii],axis=1)
    dy0 = median(dy[ii],axis=1)
    dx0_img,ix,iy = twodbin(x,y,dx0,nbinx=20,nbiny=20)
    dy0_img = twodbin(x,y,dy0,nbinx=20,nbiny=20)[0]
    dis_img=sqrt(dx0_img**2+dy0_img**2)
    dx0 = median_bin(x,y,[],mdl=dx0_img,ix=ix,iy=iy,nbinx=20,nbiny=20)
    dy0 = median_bin(x,y,[],mdl=dy0_img,ix=ix,iy=iy,nbinx=20,nbiny=20)

    print (""" Median distortion: %.2f arcsec""" % ( median(sqrt(dx**2+dy**2))*3600. ))
    print (""" Median residual distortion: %.2f arcsec""" % ( median(sqrt((dx-dx0[h])**2+(dy-dy0[h])**2))*3600. ))

    mdl_dx_dy, cxy = dist_corr(dx0_img,dy0_img,hdr)

    st1="""PV1_1=%.8f PV1_2=%.8f PV1_4=%.8f PV1_5=%.8f PV1_6=%.8f""" % (1+cxy[0,0],cxy[1,0],cxy[2,0],cxy[3,0],cxy[4,0])
    st2="""PV2_2=%.8f PV2_1=%.8f PV2_6=%.8f PV2_5=%.8f PV2_4=%.8f""" % (cxy[0,1],1+cxy[1,1],cxy[2,1],cxy[3,1],cxy[4,1])
    sys.stderr.write(""" sethead %s %s %s\n""" % (st1,st2,stack_file))
    os.system("""sethead %s %s %s""" % (st1,st2,stack_file))

    i0 += dx0; j0 += dy0

    ra0,dec0 = ij2ad(i0,j0,hdr)
    hdu[1].data['ALPHA_J2000'] = ra0
    hdu[1].data['DELTA_J2000'] = dec0

    # make a couple of plots
    fig,ax = subplots()
    axs = ax.imshow(dis_img.T*3600.,cmap='gray_r',origin='lower',extent=[0,a,0,b])
    ax.set_title("Source Position Offsets")
    ax.set_xlabel("Image X Position")
    ax.set_ylabel("Image Y Position")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb=fig.colorbar(axs,cax=cax,orientation='vertical')
    cb.set_label('arcsec')
    fig.savefig("distortion.jpg")

    fwhm00 = hdr['FWHM']*ps
    fwhm = hdu[1].data['FWHM_IMAGE'][h].clip(fwhm00/ps)*ps
    fwhm0 = median(fwhm[ii],axis=1)
    fwhm_img = twodbin(x,y,fwhm0,nbinx=20,nbiny=20)[0]
    fwhm0 = median_bin(x,y,[],nbinx=20,nbiny=20,mdl=fwhm_img,ix=ix,iy=iy).clip(1./3600)
    
    fig,ax = subplots()
    ax.set_title("Source FWHM")
    axs = ax.imshow(fwhm_img.T*3600,cmap='gray_r',origin='lower',extent=[0,a,0,b])
    ax.set_xlabel("Image X Position")
    ax.set_ylabel("Image Y Position")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb=fig.colorbar(axs,cax=cax,orientation='vertical')
    cb.set_label('arcsec')
    fig.savefig("fwhm.jpg")

    #
    # check for bad matches, dis>fwhm
    #
    if (os.path.exists(mfile)):
        sys.stderr.write(""" keeping existing mfile: %s\n""" % mfile)
    else:
        h = id>0
        idh = id[h]-1
        dis = sqrt((i[idh]-i0[h])**2+(j[idh]-j0[h])**2)

        good_match = dis<=fwhm0[h]
        sys.stderr.write("""Ignoring %d bad matches (large offset) to %s stars\n""" % (h.sum()-good_match.sum(),cfilter))
        h[h] = good_match
        id[~h] = 0

        # record for the match file
        res[0,h] = (dis[good_match]*10*3600).round().astype('int16')
        res[1,h] = (mag[idh[good_match]]*10).round().astype('int16')
        res[2,h] = 10

        #
        # check for new matches, image i0,j0 compared to catalog i,j
        #
        matched = zeros(len(i),dtype='int16')
        matched[id[h]-1]=1
        unmatched = matched==0

        h1 = id==0
        xy0=vstack((i0[h1],j0[h1])).T
        tre=cKDTree(xy0)
        xy=vstack((i[unmatched],j[unmatched])).T
        dis,ii = tre.query(xy)  # nearest detected source to each catalogued source
        wh1 = where(h1)[0]
        new = dis<=fwhm0[wh1[ii]]
        if (new.sum()>0):
            sys.stderr.write("""Considering %d new matches to %s stars\n""" % (new.sum(),cfilter) )
            wh2 = where(unmatched)[0]
            wh1 = wh1[ii[new]]
            id[wh1] = wh2[new] + 1
            # record for the match file
            res[0,wh1] = (dis[new]*10*3600).round().astype('int16')
            res[1,wh1] = (mag[wh2[new]]*10).round().astype('int16')
            res[2,wh1] = 10
            h = id>0

        hdu[1].data['VECTOR_ASSOC'][h]=id[h]
        hdu[1].data['VECTOR_ASSOC'][~h]=0

        #
        # save out a summary of the matches to mfile
        #
        cat_name = sky_file.replace('_radec.fits','')
        writeto(mfile, res,overwrite=True)
        os.system("""sethead BSCALE=0.1 NCAT=1 CATID1=%s %s""" % (cat_name,mfile))

    hdu.writeto(fits_file,overwrite=True)

    #
    # now create catalog for calibration
    #
    adata=[]
    if (os.path.exists(cal_file)):
        #ra1,dec1,mag0,dmag0 = loadtxt(cal_file,unpack=True)
        adata = loadtxt(cal_file)
        n18_0 = (mag<=18).sum()
        n18 = (adata[:,2]<=18).sum()
        if (n18<n18_0/4):
            sys.stderr.write("""The catalog does not seem complete for this field, using %s.\n""" % cfilter)
            adata=[]

    if (len(adata)>0):
        ra1,dec1,mag0,dmag0 = adata.T
        gs = len(ra1)

        if (gs>=nstars):
            i1,j1 = ad2ij(ra1,dec1,hdr)
            xy1=vstack((i1,j1)).T
            tre=cKDTree(xy1)
            xy=vstack((i0,j0)).T
            dis,ii1=tre.query(xy)
     
            h1 = dis<10/3600.; ii1 = ii1[h1]
            gs = h1.sum()

        if (gs<nstars):
            sys.stderr.write("""The catalog does not seem complete for this field (nstars=%d<%d), using %s.\n""" % (gs,nstars,cfilter))
            adata=[]

    # now write out the results

    # first record usno or gaia star mags
    idh=id[h]-1
    midh,dmidh=mag[idh],dmag[idh]
    h[h] = (midh<999)*(dmidh<999)
    idh=id[h]-1

    # now record results
    if (len(adata)>0):
        hdu[1].data['ALPHA_J2000'][~h1]=0; hdu[1].data['DELTA_J2000'][~h1]=0
        hdu[1].data['ALPHA_J2000'][h1] = ra1[ii1]; hdu[1].data['DELTA_J2000'][h1] = dec1[ii1]
        hdu[1].data['FLUX_APER'][~h1]=0; hdu[1].data['FLUXERR_APER'][~h1]=0
        hdu[1].data['FLUX_APER'][h1,0] = 10**(-0.4*(mag0[ii1]-25))
        hdu[1].data['FLUX_APER'][h1,1] = 10**(-0.4*(mag0[ii1]-25))
        hdu[1].data['FLUXERR_APER'][h1,0] = 0.9210*dmag0[ii1]*hdu[1].data['FLUX_APER'][h1,0]
        hdu[1].data['FLUXERR_APER'][h1,1] = 0.9210*dmag0[ii1]*hdu[1].data['FLUX_APER'][h1,1]
    else:
        hdu[1].data['ALPHA_J2000'][~h]=0; hdu[1].data['DELTA_J2000'][~h]=0
        hdu[1].data['ALPHA_J2000'][h]=ra[idh]; hdu[1].data['DELTA_J2000'][h]=dec[idh]
        hdu[1].data['FLUX_APER'][~h]=0; hdu[1].data['FLUXERR_APER'][~h]=0
        hdu[1].data['FLUX_APER'][h,0] = 10**(-0.4*(mag[idh]-25))
        hdu[1].data['FLUX_APER'][h,1] = 10**(-0.4*(mag[idh]-25))
        hdu[1].data['FLUXERR_APER'][h,0] = 0.9210*dmag[idh]*hdu[1].data['FLUX_APER'][h,0]
        hdu[1].data['FLUXERR_APER'][h,1] = 0.9210*dmag[idh]*hdu[1].data['FLUX_APER'][h,1]

    dir=os.path.dirname(fits_file) or "."
    hdu.writeto(dir+'/calibration.fits',overwrite=True)
    os.system("""sethead CFILTER='%s' %s""" % (cfilter,dir+'/calibration.fits'))

    o=ones(len(hdu[1].data))
    writeto(dir+'/calibration.fits.apweight',vstack((hdu[1].data['FLUX_APER'][:,0],hdu[1].data['FLUX_APER'][:,0],o,o)).T,overwrite=True)


def main():
    """
    """

    if (len(sys.argv)<4): usage()

    stack_file=sys.argv[1]
    if (os.path.exists(stack_file)==False): usage()

    fits_file=sys.argv[2]
    if (os.path.exists(fits_file)==False): usage()

    sky_file=sys.argv[3]
    if (os.path.exists(sky_file)==False): usage()

    cal_file='null.txt'
    if (len(sys.argv)>4): cal_file=sys.argv[4]

    mfile='catalog_matches.fits'
    if (len(sys.argv)>5): mfile=sys.argv[5]

    sky2cal(stack_file,fits_file,sky_file,cal_file=cal_file,mfile=mfile)

if __name__ == "__main__":
    main()
