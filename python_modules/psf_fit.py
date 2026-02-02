#!/usr/bin/python3
"""
 psf_fit.py fits_file xyfile [sat_level] [varypos] [aperture] [nospatialvary]

   if psf.fits exists, use that as the psf model

"""

import argparse
import sys,os
from astropy.io.fits import PrimaryHDU,getdata
from quick_mode import quick_mode
from math import ceil
from numpy import zeros,linspace,exp,sqrt,median,abs,log,log10,loadtxt,ones,array,real,where,arange,log2,floor,dot,newaxis,vstack,inf,isinf,isnan,diag,unique,pi,cos,s_,arctan2,sin,delete,sign,savetxt
from scipy.interpolate import RectBivariateSpline
from scipy.linalg import eigh
from scipy.optimize import fmin
from refine_pos import refine_pos
from fit_wcs import ad2xy
from weight2rms import make_rms


def get90range(x):
    """ get a robust 90conf interval for x """
    x0 = median(x)
    dx0 = 1.48*median(abs(x-x0))
    j = abs(x-x0) < 2*dx0
    nj = j.sum()
    ii = x[j].argsort()
    x0 = x[j][ii[int(0.05*nj)]]
    x1 = x[j][ii[int(0.95*nj)]]
    return x0,x1


def gauss2dfitpsf(dat):
    """
     fit a 2d gaussian + const. to a psf core
    """

    n=len(dat)
    sz = (n-1)/2.

    #
    # determine starting guesses for gaussian widths and rotation angle
    #
    i,j = where(dat>dat.max()/2.)
    rad2 = (i-sz)**2+(j-sz)**2
    i0 = rad2.argmax()
    fwhm = sqrt(rad2[i0])*2.*0.9

    theta = arctan2(i[i0]-sz,j[i0]-sz)
    i1 = cos(theta)*(i-sz) - sin(theta)*(j-sz)
    fwhm1 = abs(i1).max()*2*0.9

    #
    # now fit
    #
    xx = linspace(-sz,sz,n)
    h = ( (xx**2)[:,newaxis] + (xx**2)[newaxis,:] ) <= (1+sqrt(rad2[i0]))**2
    dd = dat[h] - dat[h].mean()

    def fitfun(par):
        s,s1,theta,x0,y0 = par
        x = ( ((xx-x0)*cos(theta)/s)[:,newaxis] + ((xx-y0)*sin(theta)/s)[newaxis,:] )[h]
        y = ( (-(xx-x0)*sin(theta)/s1)[:,newaxis] + ((xx-y0)*cos(theta)/s1)[newaxis,:] )[h]
        gh = exp(-0.5*(x**2+y**2)); gh1 = gh-gh.mean()
        r = dot(dd,gh1)/dot(gh1,gh1)
        return ( (dd-r*gh1)**2 ).mean()

    s,s1,theta,x0,y0 = fmin(fitfun,[fwhm/2.35,fwhm1/2.35,theta,0.,0.],disp=False)
    fwhm = sqrt(s*s1)*2.35

    x = ((xx-x0)*cos(theta)/s)[:,newaxis] + ((xx-y0)*sin(theta)/s)[newaxis,:]
    y = (-(xx-x0)*sin(theta)/s1)[:,newaxis] + ((xx-y0)*cos(theta)/s1)[newaxis,:]
    gauss = exp(-0.5*(x**2+y**2))

    gh = gauss[h]; g0 = gh.mean(); gh1 = gh-g0
    r = dot(dd,gh1)/dot(gh1,gh1)

    return r*gauss,fwhm


def dewing(rad,psf,fwhm=10.):
    """
     fit a powerlaw to the psf wings
    """
    h = (psf<=0)*(rad>fwhm)
    resw = []
    if (h.sum()>0):
        rmin = rad[h].min()
        r1 = floor(rmin*4/fwhm)*fwhm/4
        sys.stderr.write("""Smoothing PSF for radius>%.1f fwhm\n""" % (r1/fwhm))
        h = rad>=r1
        def fitfun(par):
            if (par[1]>-2.1): par[1]=-2.1
            m = par[0]*(rad[h]/fwhm)**par[1]
            return ( (m-psf[h])**2 ).mean()

        fn = (rad[h]/fwhm)**(-4.)
        p0 = dot(psf[h],fn)/dot(fn,fn)

        resw = fmin(fitfun,[p0,-4.],disp=False)
        psf[h] = resw[0]*(rad[h]/fwhm)**resw[1]
    else:
        sys.stderr.write("""Not smoothing psf wings.\n""")

    psf[psf<=0] = 0.

    return resw


def get_psf(dat,rdat,dmag,x,y,xa,ya,x2a,y2a,idx,mask=[],fwhm=10.,outfile='psf.fits',max_flux=5.e4,oversamp=4,nsig_clip=5.,nfwhm=3,nmax=100,mag_min=16.,spatial_vary=True):
    """
      Determines the mean psf for an image.
        infile: input fits image
        mask: optional input fits mask image (0 for star, 1 otherwise)
        xyfile: sextractor file with star locations, run_sex1.sh format
        fwhm: stellar fwhm
        nfwhm: determine psf out to a radius of nfwhm*fwhm
    """

    psf=[]
    psf_extra=[]
    rese=[]
    resw=[]
    sz0 = nfwhm*int(round(fwhm))
    sz = oversamp*sz0
    sz00 = 8*int(round(fwhm))
    sx,sy = dat.shape

    xx = linspace(-sz0,sz0,2*sz+1)
    rad = sqrt( (xx**2)[:,newaxis] + (xx**2)[newaxis,:] )
    hr = rad<=fwhm

    # select only unsaturated sources away from edges and other sources
    n=len(x)
    jj = zeros(n,dtype='bool')
    bg = zeros(n,dtype='float64')
    for i in range(n):
        j,k = int(round(y[i]))-1,int(round(x[i]))-1
        r = rdat[j-sz00:j+sz00+1,k-sz00:k+sz00+1]
        d = dat[j-sz0:j+sz0+1,k-sz0:k+sz0+1]
        dist = sqrt( (x[i]-x)**2+(y[i]-y)**2 )
        dist[i] = 999.
        if ( (r==0).sum() == 0 and r.shape == (1+2*sz00,1+2*sz00) and d.max()<=max_flux and dist.min()>fwhm*2 ):
            jj[i] = True
            if (len(mask)>0):
                d0 = dat[j-sz00:j+sz00+1,k-sz00:k+sz00+1]
                m = mask[j-sz00:j+sz00+1,k-sz00:k+sz00+1]
                if (m.sum()>0):
                    dd0 = d0[m]
                    bg[i] = quick_mode(dd0)
                    #bg[i] = median(dd0)
                    #dbg = 1.48*median(abs(dd0-bg[i]))
                    #h = abs(dd0-bg[i])<2*dbg
                    #if (h.sum()>0): bg[i] = 2.5*bg[i] - 1.5*dd0[h].mean()

    if (jj.sum()>0):
        x = x[jj]; y = y[jj]; dmag = dmag[jj]
        xa = xa[jj]; ya = ya[jj]
        x2a = x2a[jj]; y2a = y2a[jj]; bg = bg[jj]; idx = idx[jj]
        n = jj.sum()

        # only need nmax sources at most
        if (n>nmax):
            x = x[:nmax]; y = y[:nmax]; dmag = dmag[:nmax]
            xa = xa[:nmax]; ya = ya[:nmax]
            x2a = x2a[:nmax]; y2a = y2a[:nmax]; idx = idx[:nmax]; bg = bg[:nmax]
            n = nmax

        # now, create and fill postage stamps around each source
        dat1 = zeros((n,1+2*sz0,1+2*sz0),dtype='float64')
        vdat1 = zeros((n,1+2*sz0,1+2*sz0),dtype='float64')

        xx0 = linspace(-sz0,sz0,2*sz0+1)
        rad0 = sqrt( (xx0**2)[:,newaxis] + (xx0**2)[newaxis,:] )

        # be sure to unmask a suitable region around each star to allow for fitting
        #  (the maskfile blocks all stars by default)
        if (len(mask)>0):
            nf = int(floor((rad0.max()-fwhm/4)*4/fwhm))-1
            for i in range(n):
                j,k = int(round(y[i]))-1,int(round(x[i]))-1
                dat1[i] = dat[j-sz0:j+sz0+1,k-sz0:k+sz0+1]
                msk = mask[j-sz0:j+sz0+1,k-sz0:k+sz0+1].copy()

                r1 = fwhm*nfwhm
                for l in range(nf):
                    h = (rad0>fwhm+l*fwhm/4.)*(rad0<=fwhm+(l+1)*fwhm/4.)
                    if (msk[h].sum()>0.5*h.sum()):
                        r1 = fwhm+(l+1)*fwhm/4.
                        break

                msk[rad0<=r1] = True
                r = rdat[j-sz0:j+sz0+1,k-sz0:k+sz0+1]
                msk[r<=0] = False
                vdat1[i,msk] = r[msk]**2
        else:
            for i in range(n):
                j,k = int(round(y[i]))-1,int(round(x[i]))-1
                dat1[i] = dat[j-sz0:j+sz0+1,k-sz0:k+sz0+1]
                r = rdat[j-sz0:j+sz0+1,k-sz0:k+sz0+1]
                h=r>0
                vdat1[i,h] = r[h]**2

        h = dat1>=max_flux
        dat1[h] = 0.
        vdat1[h] = 0.

        # correctly (sub-pixel) center the stamps and weights at the sextractor positions
        # we will oversample
        dat = zeros((n,2*sz+1,2*sz+1),dtype='float64')
        vdat = zeros((n,2*sz+1,2*sz+1),dtype='float64')
        dx,dy = y-y.round(),x-x.round()
        for i in range(n):
            dat[i] = RectBivariateSpline(xx0,xx0,dat1[i])(xx+dx[i],xx+dy[i]) - bg[i]
            vdat[i] = RectBivariateSpline(xx0,xx0,vdat1[i])(xx+dx[i],xx+dy[i]).clip(0)

        vdat1=0.; dat1=0.

        # make sure each source is shaped like a majority of the others
        #  sources i and j with matr[i,j] > 0.9 or so have similar shapes (see, "score" below)
        dhr = dat[:,hr]; shr = sqrt(vdat[:,hr])
        matr = dot(dhr,dhr.T) - dot(shr,shr.T); matrt = matr.T
        dhr,shr=0.,0.
        dd = diag(matr)
        h = dd>0; norm = sqrt(dd[h])
        matr[:,h] /= norm; matrt[:,h] /= norm; matr[~h,:] = 0; matr[:,~h] = 0

        f = median(matr,axis=0); matr=0.
        g = f>0.9; ng=g.sum()
        n0 = n
        if (ng>0.5*n):
            if (g.sum()>5):
                dat = dat[g]; vdat = vdat[g]; n = ng;
                x = x[g]; y = y[g]; dmag = dmag[g]
                xa = xa[g]; ya = ya[g]
                x2a = x2a[g]; y2a = y2a[g]; idx = idx[g]

        # iteratively determine the psf:
        #   guess psf, determine flux f, refine psf, determine flux,...
        psf = exp(-0.5*(rad*2.35/fwhm)**2)
        f = zeros(n,dtype='float64')
        for i in range(3):
            norm = dot(vdat[:,hr],psf[hr]**2)
            g = norm>0
            f[g] = dot(dat[g][:,hr]*vdat[g][:,hr],psf[hr]) / norm[g]
            g *= f>0

            datf = dat[g]/f[g,newaxis,newaxis]; datf1 = datf.copy()
            h = vdat[g]==0
            datf[h] = inf; datf1[h] = -inf
            psf = median(vstack((datf,datf1)),axis=0)
            psf[isinf(psf)] = 0
            psf[isnan(psf)] = 0

        #
        # smooth out the psf wings
        #
        resw = dewing(rad,psf,fwhm=fwhm)

        # now with decent psf and flux estimates, sigma-clip bad pixels 
        # don't use deviant or low snr points to estimate the psf
        h = (dat - f[:,newaxis,newaxis]*psf[newaxis,:,:])**2 > nsig_clip**2*( vdat + (0.1*dat.clip(0))**2 )
        vdat[h] = 0.
        h = vdat > ((f/3.)**2)[:,newaxis,newaxis]*(psf**2)[newaxis,:,:]
        vdat[h] = 0.

        vnorm = dot(vdat[g].T,f[g]**2).T
        h = vnorm>0
        if (h.sum()>0): psf[h] = dot((dat[g]*vdat[g]).T,f[g]).T[h] / vnorm[h]

        #
        # now use the psf core (gaussian fit) to determine possible (quadratic) psf spatial variation
        #  also refine the fwhm estimate
        #
        ng = g.sum()
        if (spatial_vary and ng>5):
            psfg0,fwhm1 = gauss2dfitpsf(psf)
            if (fwhm1<2*fwhm*oversamp): fwhm = fwhm1/oversamp
            else: psfg0 = psf.max()*exp(-0.5*(2.35*rad/fwhm)**2)

            hr2 = rad<=fwhm*2.
            psfhr = psf[hr2]; ps = psf.sum(); ps2 = dot(psfhr,psfhr)
            psfghr = psfg0[hr2] - psfhr*dot(psfg0[hr2],psfhr)/ps2
            psg = psfghr.sum()
            area = len(psf)**2

            fac,fac1 = zeros(ng,dtype='float64'), zeros(ng,dtype='float64')
            k=0
            for i in range(n):
                if (g[i]):
                    d,v = dat[i,hr2], vdat[i,hr2]
                    nm1,nm2,nm3 = dot(psfhr,psfhr*v),dot(psfghr,psfghr*v),v.sum()
                    f0=0.
                    if (nm1>0 and nm2>0 and nm3>0):
                        f0 = dot(d,psfhr*v)/nm1
                        r = d-f0*psfhr; fg = dot(r,psfghr*v)/nm2
                        r -= fg*psfghr; fk = dot(r,v)/nm3
                        f0 += fg*psg/ps + fk*area/ps
                    if (f0>0):
                        fac[k],fac1[k] = fg/f0,fk/f0
                        k+=1
                    else: g[i],fac,fac1 = False,delete(fac,k),delete(fac1,k)

            if (len(fac)>0):
                g1 = abs(fac)<1
                if (g1.sum()==0): g1 = ones(n,dtype='bool')
                for i in range(5):
                    fc0 = median(fac[g1])
                    dfc0 = 1.48*median(abs(fac[g1]-fc0))
                    g2 = abs(fac-fc0)<3*dfc0
                    if (g2.sum()>5): g1=g2

                fac = fac[g1]; fac1 = fac1[g1]
                g[g] = g1; ng = g.sum()

                psf_extra = zeros((2*sz+1,2*sz+1),dtype='float64')
                psf_extra[hr2] = psfghr

                xm,ym,x2m,y2m = xa[g].mean(),ya[g].mean(),x2a[g].mean(),y2a[g].mean()
                xa=xa[g]-xm; ya=ya[g]-ym; x2a=x2a[g]-x2m; y2a=y2a[g]-y2m
                xya = vstack((xa,ya,x2a,xa*ya,y2a))
                def param_mdl_fit(par,ff):
                    # t-dist with 5 dof
                    return 6.*log( 1 + (par[0]+dot(par[1:],xya) - ff)**2/5. ).mean()

                sys.stderr.write("""Determining PSF spatial variation.\n""")
                par0,rese = zeros(6,dtype='float64'), zeros((2,6),dtype='float64')
                rese[0] = fmin(param_mdl_fit,par0,args=(fac,),disp=False); rese[0] = fmin(param_mdl_fit,rese[0],args=(fac,),disp=False)
                rese[1] = fmin(param_mdl_fit,par0,args=(fac1,),disp=False); rese[1] = fmin(param_mdl_fit,rese[1],args=(fac1,),disp=False)
                rese[:,0] -= rese[:,1]*xm + rese[:,2]*ym + rese[:,3]*x2m - rese[:,4]*xm*ym + rese[:,5]*y2m
                rese[:,1] -= rese[:,4]*ym; rese[:,2] -= rese[:,4]*xm

    # almost done
    if (len(psf)==0):
        sys.stderr.write("""PSF determination failed, using Gaussian\n""")
        psf = exp(-0.5*(rad*2.35/fwhm)**2)

    norm = psf.sum()/oversamp**2
    psf /= norm
    if (len(resw)>0): resw[0]/=norm
    if (len(rese)>0): rese[1]/=norm
    if (len(psf_extra)>0): psf_extra/=norm
    mdpsf = 1./sqrt( (1./dmag**2).sum() )

    # finally save out some potentially useful information
    sys.stderr.write("PSF Precision: %.2e\n""" % mdpsf )
    a,b = psf.shape
    psfa = zeros((2,a,b),dtype='float64')
    psfa[0] = psf
    if (len(psf_extra)>0):psfa[1] = psf_extra
    hdu = PrimaryHDU(psfa)
    hdu.header['DPSF'] = mdpsf
    hdu.header['FWHM'] = fwhm
    hdu.header['OVERSAM'] = oversamp
    if (len(resw)>0):
        hdu.header['WINGN'] = resw[0]
        hdu.header['WINGI'] = resw[1]
    if (len(psf_extra)>0):
        hdu.header['RES0'] = str(rese[0].tolist())
        hdu.header['RES1'] = str(rese[1].tolist())
        hdu.header['PVARY'] = 1
    else: hdu.header['PVARY'] = 0
    hdu.writeto(outfile,overwrite=True)

    if (len(rese)>0):
        k=0
        fl = open('scorr_file.txt','w')
        for i in range(n):
            if (g[i]>0):
                xya = array([1,xa[k]+xm,ya[k]+ym,x2a[k]+x2m,(xa[k]+xm)*(ya[k]+ym),y2a[k]+y2m])
                ff = dot(rese,xya)
                fl.write("""%f %f %f %f %f %f %f %f %f %d\n""" % (xya[1],xya[2],xya[3],xya[5],fac[k],fac1[k]/norm,ff[0],ff[1],dmag[i],idx[i]))
                k+=1

        fl.close()

    return 1


def fit_psf(dat,vdat,mag,dmag,x,y,xa,ya,x2a,y2a,psffile,fit_rad=1.,ap_rad=1.,max_flux=5.e4,nsigma=3.,add_sys=True,aperture=False,remove_contam=True,spatial_vary=True,nsig_clip=3.):
    go=True
    """
      Fits the mean psf for an image to each source to determine flux.
        dat: input fits image
        rdat: input fits variance image
        aperture: True for aperture photometry with radius fit_rad
        remove_contam: True to remove psf contamination during aperture photometry
        fit_rad: factor times fwhm in psffile for fitting

      Note: returns the nsigma limit if flux is <= nsigma * dflux
    """
    if (os.path.exists(psffile)):
        psf,hdra = getdata(psffile,header=True)
        psf,epsf = psf[0], psf[1]
        fwhm = hdra['FWHM']
        fit_rad *= fwhm
        ap_rad *= fwhm
        dpsf = hdra['DPSF']
        oversamp = hdra['OVERSAM']
        if (hdra['PVARY']==1 and spatial_vary):
            sys.stderr.write("""Fitting PSF with spatial variations.\n""")
            res = zeros((2,6),dtype='float64')
            res[0] = array(hdra['RES0'].strip('[]').split(',')).astype('float64')
            res[1] = array(hdra['RES1'].strip('[]').split(',')).astype('float64')
            eps = epsf.sum()/oversamp**2
        else: res=[]
    else:
        sys.stderr.write("""Warning: psf file %s not found!\n""" % psffile)

    sx,sy = dat.shape
    sz = len(psf)
    sz0  = (sz-1)//(2*oversamp)
    if (fit_rad>sz0): fit_rad=1.*sz0
    area = (2*sz0+1)**2
    score,chi2 = 0.*mag, 0.*mag

    vdat[dat>max_flux] = 0.

    # create functions to shift psf's for use below
    xx0 = linspace(-sz0,sz0,sz)
    xx = linspace(-sz0,sz0,2*sz0+1)
    psf_shift0 = RectBivariateSpline(xx0,xx0,psf)
    def psf_shift(xx1,yy1,xa0=0.,ya0=0.,x2a0=0.,y2a0=0.):
        return psf_shift0(xx1,yy1)

    # potential psf spatial variation
    if (len(res)>0):
        psfe_shift = RectBivariateSpline(xx0,xx0,epsf)
        def psf_shift(xx1,yy1,xa0=0.,ya0=0.,x2a0=0.,y2a0=0.):
            ff = dot(res, array([1,xa0,ya0,x2a0,xa0*ya0,y2a0]) )
            return psf_shift0(xx1,yy1)*(1-ff[0]*eps-ff[1]*area) + psfe_shift(xx1,yy1)*ff[0] + ff[1]

    ap_corr0=1.
    if (aperture):
        sys.stderr.write("""Performing aperture photometry with radius %.2f\n""" % ap_rad)
        rad0 = sqrt( (xx0**2)[:,newaxis] + (xx0**2)[newaxis,:] )
        psfa = 0.*rad0
        psfa[rad0<=ap_rad-0.5] = 0.5
        psfa[rad0<=ap_rad+0.5] += 0.5
        psfa_shift = RectBivariateSpline(xx0,xx0,psfa,kx=1,ky=1)

        # estimate aperture correction
        p = psf_shift(xx,xx,xa0=xa.mean(),ya0=ya.mean(),x2a0=x2a.mean(),y2a0=y2a.mean())
        pa = psfa_shift(xx,xx)
        ap_corr0 = (p*pa).sum()

    # holders for final flux and error
    n=len(x)
    flx = 10**(-0.4*(mag-25.))
    dflx = zeros(n,dtype='float64')

    #
    # now compute flx and dflx for each source
    #
    i,j = y.round().astype('int16')-1, x.round().astype('int16')-1
    dx,dy = y-i-1,x-j-1

    # watch out for image boundaries
    i1,i2 = (i-sz0).clip(0),(i+sz0+1).clip(0,sx)
    j1,j2 = (j-sz0).clip(0),(j+sz0+1).clip(0,sy)
    i10,i20 = i1-i+sz0,i2-i+sz0
    j10,j20 = j1-j+sz0,j2-j+sz0

    # assign parents to children (overlapping psf's)
    parent = zeros(n,dtype='int16')-1
    ii = mag.argsort()
    sz02 = sz0**2
    for i in ii:
        dis2 = (x-x[i])**2 + (y-y[i])**2
        j = (dis2<sz02)*(parent==-1)
        p0 = parent[i]
        if (p0==-1): parent[j] = i
        else: parent[j] = p0
        parent[i] = p0

    # figure out which ones we just want to punt on
    ii = parent.argsort()
    h = (i2[ii]>i1[ii])*(j2[ii]>j1[ii])
    mag[ii[~h]],dmag[ii[~h]] = 999.,999.
    ii = ii[h]; parent = parent[ii]
    n1 = len(ii)

    # mapping of that input psf array to the data:
    hr = ( (xx**2)[:,newaxis] + (xx**2)[newaxis,:] ) <= fit_rad**2

    #
    # fit the psf, fitting multiple psf's to grouped sources
    #
    nparent=[]
    ib = 0
    while (ib<n1):
        ie = ib
        while (parent[ie]==parent[ib]):
            ie+=1
            if (ie==n1 or parent[ib]==-1): break

        neb = ie-ib
        ii1 = ii[ib:ie]
        nparent.append(neb)
        #
        # do work with neb overlapping sources numbered ii[ib:ie]
        #

        matr = zeros((neb,neb),dtype='float64')
        vmatr = zeros((neb,neb),dtype='float64')
        vec = zeros(neb,dtype='float64')
        vec0 = zeros(neb,dtype='float64')

        if (aperture):
            if (remove_contam):
                matra = zeros((neb,neb),dtype='float64')
                vmatra = zeros((neb,neb),dtype='float64')
            ap_corr = zeros(neb,dtype='float64')
            veca = zeros(neb,dtype='float64')
            dveca = zeros(neb,dtype='float64')

        i1a,i2a = i1[ii1].min(),i2[ii1].max()
        j1a,j2a = j1[ii1].min(),j2[ii1].max()

        # define the extraction regions and mask, wt
        p = []
        wt = zeros((i2a-i1a,j2a-j1a),dtype='bool')
        for k in range(neb):
            i = ii[ib+k]
            pstamp1 = s_[i10[i]:i20[i],j10[i]:j20[i]]
            dstamp1 = s_[i1[i]-i1a:i2[i]-i1a,j1[i]-j1a:j2[i]-j1a]
            wt[dstamp1] += hr[pstamp1]

            # p: psf to be extracted on selection region
            p.append( psf_shift(xx-dx[i],xx-dy[i],xa[i],ya[i],x2a[i],y2a[i])[pstamp1] )

        wtb = ~wt
        wt0 = wt.copy()
        h = vdat[i1a:i2a,j1a:j2a]<=0
        wt[h] = False; wtb[h] = False

        # now calculate the psf or aperture photometry background
        if (wtb.sum()>0):
            d = dat[i1a:i2a,j1a:j2a][wtb]
            #bg = quick_mode(d)
            bg = 0.
            #bg = median(d)
            #dbg = 1.48*median(abs(d-bg))
            #h = abs(d-bg)<2*dbg
            #if (h.sum()>0): bg = 2.5*bg-1.5*d[h].mean()
        else: bg=0.

        # now calculate the psf or aperture photometry
        for k in range(neb):
            i = ii[ib+k]
            pstamp = s_[i10[i]:i20[i],j10[i]:j20[i]]
            dstamp = s_[i1[i]:i2[i],j1[i]:j2[i]]
            dstamp1 = s_[i1[i]-i1a:i2[i]-i1a,j1[i]-j1a:j2[i]-j1a]

            d,vd = dat[dstamp]-bg,vdat[dstamp]
            hr0 = wt[dstamp1]
            phr,phr1 = p[k][hr0], p[k][wt0[dstamp1]]

            # we will solve flux = matr^(-1) vec
            vec[k] = dot(d[hr0],phr)
            matr[k,k], vmatr[k,k] = dot(phr,phr), dot(phr,phr*vd[hr0])
            # systematic error estimate for saturated sources
            if (matr[k,k]>0): vec0[k] = 0.01*sqrt(abs(dot(phr1,phr1)/matr[k,k] - 1.))

            if (aperture):
                # pa, aperture photometry window function
                pa = psfa_shift(xx-dx[i],xx-dy[i])[pstamp]*(vd>0)
                phra = pa.ravel()
                veca[k] = dot(d.ravel(),phra)
                dveca[k] = sqrt( dot(vd.ravel(),phra**2) )
                ap_corr[k] = dot(p[k].ravel(),phra)/ap_corr0

            # calculate the psf overlap between sources (covariance matrix)
            for l in range(k):
                j = ii[ib+l]
                ist,isp = max(i1[i],i1[j]), min(i2[i],i2[j])
                jst,jsp = max(j1[i],j1[j]), min(j2[i],j2[j])
                if (isp>ist and jsp>jst):
                    slca = s_[ist-i1[i]:ist-i1[i]+isp-ist,jst-j1[i]:jst-j1[i]+jsp-jst]
                    slcb = s_[ist-i1[j]:ist-i1[j]+isp-ist,jst-j1[j]:jst-j1[j]+jsp-jst]
                    hr1 = hr0[slca]; pk,pl = p[k][slca][hr1],p[l][slcb][hr1]
                    matr[k,l] = matr[l,k] = dot(pk,pl)
                    vmatr[k,l] = vmatr[l,k] = dot( vd[slca][hr1]*pk,pl )

            if (aperture and remove_contam and neb>1):
                # determine fraction of psf l present in aperture k
                for l in range(neb):
                    if (l==k): continue
                    j = ii[ib+l]
                    ist,isp = max(i1[i],i1[j]), min(i2[i],i2[j])
                    jst,jsp = max(j1[i],j1[j]), min(j2[i],j2[j])
                    if (isp>ist and jsp>jst):
                        slca = s_[ist-i1[i]:ist-i1[i]+isp-ist,jst-j1[i]:jst-j1[i]+jsp-jst]
                        slcb = s_[ist-i1[j]:ist-i1[j]+isp-ist,jst-j1[j]:jst-j1[j]+jsp-jst]
                        matra[k,l] = dot(pa[slca].ravel(),p[l][slcb].ravel())
                        vmatra[k,l] = dot((pa[slca]*vd[slca]).ravel(),p[l][slcb].ravel())

        if (neb>1):
            # replace matr with it's inverse
            eig,mm = eigh(matr); eig = abs(eig)
            h=eig==0
            if ((~h).sum()>0):
                eig[h] += eig[~h].min()*1.e-6
                matr = dot(mm/eig,mm.T)
            else:
                matr = zeros((neb,neb),dtype='float64')
        elif(matr[0,0]!=0): matr[0,0] = 1./matr[0,0]

        # the psf flux
        flx0 = dot(matr,vec)

        # determine the chi^2 for each source and the score
        for k in range(neb):
            i = ii[ib+k]
            pstamp = s_[i10[i]:i20[i],j10[i]:j20[i]]
            dstamp = s_[i1[i]:i2[i],j1[i]:j2[i]]

            d,vd = dat[dstamp],vdat[dstamp]

            mdl = flx0[k]*p[k]
            for l in range(neb):
                if (l==k): continue
                j = ii[ib+l]
                ist,isp = max(i1[i],i1[j]), min(i2[i],i2[j])
                jst,jsp = max(j1[i],j1[j]), min(j2[i],j2[j])
                if (isp>ist and jsp>jst):
                    slca = s_[ist-i1[i]:ist-i1[i]+isp-ist,jst-j1[i]:jst-j1[i]+jsp-jst]
                    slcb = s_[ist-i1[j]:ist-i1[j]+isp-ist,jst-j1[j]:jst-j1[j]+jsp-jst]
                    mdl[slca] += flx0[l]*p[l][slcb]

            hr0 = hr[pstamp]*(vd>0)
            hs = hr0.sum()
            if (hs>0):
                chi2[i] = dot( (d[hr0] - mdl[hr0])**2,1./(vd[hr0]+(0.1*mdl[hr0])**2) )/hs

                pkhr0 = p[k][hr0]
                src = d[hr0]-mdl[hr0]+flx0[k]*pkhr0
                dd = dot(src,src) - vd[hr0].sum()
                dd1 = dot(pkhr0,pkhr0)
                if (dd>0 and dd1>0): score[i] = dot(src,pkhr0) / sqrt(dd*dd1)

        vmatr1 = dot(matr,dot(vmatr,matr))
        if (aperture):
            # do not include psf estimates of fluxes from other sources in the aperture
            flx[ii1],dflx[ii1] = veca,dveca
            if (remove_contam and neb>1):
                flx[ii1] -= dot(matra,flx0)
                dflx[ii1] = sqrt( ( dflx[ii1]**2 + ((dot(matra,vmatr1)-2*dot(vmatra,matr))*matra).sum(axis=1) ).clip(0) )
            h = ap_corr>0
            flx[ii1[h]] /= ap_corr[h]; dflx[ii1[h]] /= ap_corr[h]
        else:
            flx[ii1],dflx[ii1] = flx0, sqrt( (diag(vmatr1)).clip(0) )

        if (add_sys): dflx[ii1] = sqrt( dflx[ii1]**2 + (vec0*flx[ii1])**2 )

        ib = ie

    # summarize the PSF grouping
    nparent = array(nparent)
    np = unique(nparent)
    sys.stderr.write("""PSF Group(Frequency): """)
    for p in np:
        j=nparent==p
        sys.stderr.write("""%d(%d) """ % (p,j.sum()))

    sys.stderr.write("""\n""")

    # translate fluxes to mags
    for i0 in range(n1):
        i = ii[i0]
        if (dflx[i]<=0): mag[i],dmag[i] = 999,999
        elif (flx[i]<nsigma*dflx[i]): mag[i],dmag[i] = 25-2.5*log10(flx[i].clip(0)+nsigma*dflx[i]),999
        else: mag[i],dmag[i] = 25-2.5*log10(flx[i]),2.5/log(10.)*dflx[i]/flx[i]

    return flx,dflx,score,chi2,2.5*log10(ap_corr0)


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("fitsfile", type=str, help="input fits file")
    parser.add_argument("xyfile", type=str, help="input photometry file")
    parser.add_argument("--dmag_max", type=float, help="maximum magnitude error for psf estimate (default 0.01)", default=0.01)
    parser.add_argument("--mag_min", type=float, help="minimum magnitude error for psf estimate (default 16.)", default=16.)
    parser.add_argument("--satlevel", type=float, help="pixel saturation level (default 6.e4)", default=6.e4)
    parser.add_argument("--fitrad", type=float, help="fitting radius (units of fhwm, default 1)", default=1.)
    parser.add_argument("--aprad", type=float, help="aperture radius (units of fhwm, default 1)", default=1.)
    parser.add_argument("--novarypos", help="do not refine centroid", action="store_true")
    parser.add_argument("--useradec", help="use radec instead of xy source positions", action="store_true")
    parser.add_argument("--aperture", help="do aperture photometry", action="store_true")
    parser.add_argument("--noremove_contam", help="do not remove psf contamination in apertures", action="store_true")
    parser.add_argument("--spatial_vary", help="fit psf with spatial variations", action="store_true")
    args = parser.parse_args()

    dmag_max = args.dmag_max
    mag_min = args.mag_min
    sat_level = args.satlevel
    fit_rad = args.fitrad
    ap_rad = args.aprad
    varypos=True
    if (args.novarypos): varypos=False
    useradec=False
    if (args.useradec): useradec=True
    aperture = args.aperture
    remove_contam=True
    if (args.noremove_contam): remove_contam=False
    spatial_vary=False
    if (args.spatial_vary): spatial_vary=True

    fitsfile = args.fitsfile
    if (os.path.exists(fitsfile)==False):
        sys.stderr.write("""fitsfile %s not found, exiting...\n""" % fitsfile)
        sys.exit()
    xyfile = args.xyfile
    if (os.path.exists(fitsfile)==False):
        sys.stderr.write("""xyfile %s not found, exiting...\n""" % xyfile)
        sys.exit()

    dat,hdr = getdata(fitsfile,header=True)

    # the sextractor fwhm will be used for determining oversampling and the psf extraction
    fwhm0 = hdr['FWHM']
    oversamp = int(ceil(2.35*3*2/fwhm0))

    dirname=os.path.dirname(fitsfile) or '.'
    base=os.path.basename(fitsfile).replace('.fits','')
    wtfile=dirname+'/'+base+'.wt.fits'
    maskfile=dirname+'/'+base+'.wmap_mask.fits'

    # load in the sextractor photometry as a starting point
    (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx)=loadtxt(xyfile,unpack=True,ndmin=2)
    if (useradec): x,y = ad2xy(ra,dec,hdr)

    # recenter all but idx<0 (corresponding to forced photometry)
    n = len(mag)
    h=idx>0
    if (varypos):
        sys.stderr.write("""Refining source centroids.\n""")
        ra1 = ra[h]; dec1 = dec[h]; x1 = x[h]; y1 = y[h]
        if (os.path.exists('dual.fits')): dath=getdata('dual.fits')
        else: dath=dat
        stat = refine_pos(dath,hdr,ra1,dec1,x1,y1,fit_rad=fwhm0)
        ra[h] = ra1; dec[h] = dec1; x[h] = x1; y[h] = y1

    g = ones(n,dtype='bool')
    if ((~h).sum()>0):
        ii = where(~h)[0]
        wh = where(h)[0]
        cd = cos(median(dec)*pi/180.)
        for i in ii:
            d2 = sqrt( ((ra[i]-ra[h])*cd)**2 + (dec[i]-dec[h])**2 )*3600.
            g[wh[d2<2]] = False

        ra= ra[g]; dec= dec[g]; mag= mag[g]; dmag= dmag[g]; mag_big= mag_big[g]; dmag_big= dmag_big[g]
        fwhm= fwhm[g]; x= x[g]; y= y[g]; xa= xa[g]; ya= ya[g]; x2a= x2a[g]; y2a= y2a[g]
        expos= expos[g]; idx= idx[g]
        n1 = len(ra)
        if (n1<n): sys.stderr.write("""Dropping %d source(s) overlapping with forced positions\n""" % (n-n1))

    # load in the rms and source mask images
    wt = getdata(wtfile)
    rdat = make_rms(dat,wt,gain=hdr['GAIN'])

    mask=[]
    if (os.path.exists(maskfile)):
        # use a sextrator file generated with -CHECKIMAGE_NAME mask.fits in sextractor 
        #   to mask out sources (see run_sex1.sh and sources2mask.py)
        mask,hdrm=getdata(maskfile,header=True)
        x1,y1 = hdr['CRPIX1'],hdr['CRPIX2']
        x0,y0 = hdrm['CRPIX1'],hdrm['CRPIX2']
        dy,dx = int(x0-x1),int(y0-y1)
        sx,sy = dat.shape
        mask = mask[dx:dx+sx,dy:dy+sy].astype('bool')

    # calculate the PSF or read it in if already present
    if (os.path.exists('psf.fits')==True):
        sys.stderr.write("""Using psf file %s\n""" % 'psf.fits')
    else:

        mag1 = 1.*mag; mag_big1 = 1.*mag_big; dmag1 = 1.*dmag; x1 = 1.*x; y1 = 1.*y
        x1a = 1.*xa; y1a = 1.*ya; x12a = 1.*x2a; y12a = 1.*y2a; idxa = 1*idx
        h = idx>=0   # catalogued sources only
        if (h.sum()>0):
            mag1=mag1[h]; mag_big1=mag_big1[h]; dmag1=dmag1[h]; x1=x1[h]; y1=y1[h]; x1a=x1a[h]; y1a=y1a[h]; x12a=x12a[h]; y12a=y12a[h]; idxa=idxa[h]

        ii = mag1.argsort()
        mag1 = mag1[ii]; mag_big1 = mag_big1[ii]; dmag1 = dmag1[ii]
        x1 = x1[ii]; y1 = y1[ii]
        x1a = x1a[ii]; y1a = y1a[ii]; x12a = x12a[ii]; y12a = y12a[ii]; idxa=idxa[ii]

        #
        # select good values using ddoti routine
        #
        dis = sqrt( (x1-x1.mean())**2+(y1-y1.mean())**2 ); dmax = dis.max()*0.75
        good = (dis<=dmax)*(dmag1<=dmag_max)*(mag1>=mag_min)
        ma,mb = get90range((mag1-mag_big1)[good])
        ii = (mag1-mag_big1>=ma)*(mag1-mag_big1<=mb)*good
        #sys.stderr.write(f" 90% range for mag-mag_big: {ma:.3f} - {mb:.3f} ({ii.mean():.2f}, distance<{dmax:.1f}, mag>{mag_min}, dmag<{dmag_max})\n")
        mag1 = mag1[ii]; mag_big1 = mag_big1[ii]; dmag1 = dmag1[ii]
        x1 = x1[ii]; y1 = y1[ii]
        x1a = x1a[ii]; y1a = y1a[ii]; x12a = x12a[ii]; y12a = y12a[ii]; idxa=idxa[ii]
        #savetxt('out.txt',vstack((x1,y1)).T)

        stat = get_psf(dat,rdat,dmag1,x1,y1,x1a,y1a,x12a,y12a,idxa,mask=mask,fwhm=fwhm0,max_flux=sat_level,oversamp=oversamp,spatial_vary=spatial_vary)

    # fit the psf and report the results
    flx,dflx,score,chi2,ap_corr = fit_psf(dat,rdat**2,mag,dmag,x,y,xa,ya,x2a,y2a,fit_rad=fit_rad,ap_rad=ap_rad,psffile='psf.fits',max_flux=sat_level,aperture=aperture,remove_contam=remove_contam,spatial_vary=spatial_vary)

    f = open(xyfile,'r')
    str0 = f.readline()
    f.close()

    psf = getdata('psf.fits')[0]
    FWHMeff = 1./sqrt( (psf**2).sum() )*oversamp/1.505

    print (str0.strip())
    print ("""# 10-sigma limiting mag: 22.5 - 2.5*log10(FWHMeff*1.505*rms), APCorr= %.6f FWHMeff= %.2f""" % (ap_corr,FWHMeff))
    n = len(mag)
    for i in range(n):
        if (mag[i]<999):
            print ("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d""" % (ra[i],dec[i],mag[i],dmag[i],mag_big[i],dmag_big[i],fwhm[i],x[i],y[i],xa[i],ya[i],x2a[i],y2a[i],expos[i],idx[i]))

    f = open('psf_stats.txt','w')
    f.write("""# score chi2 Flux dFlux id \n""")
    for i in range(n):
        if (mag[i]<999):
            f.write("""%10.2f %20.2f %20.6e %20.6e %d\n""" % (score[i],chi2[i],flx[i],dflx[i],idx[i]))
    f.close()
