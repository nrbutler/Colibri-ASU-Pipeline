#!/usr/bin/python3
"""
 refine_pos.py fits_file [fwhm]
"""

import sys,os
from math import ceil
from astropy.io.fits import getdata
from numpy import linspace,exp,sqrt,loadtxt,dot,newaxis,where,pi,arange,zeros
from fit_wcs import xy2ad
from scipy.special import erf


def usage():
    print (__doc__)
    sys.exit()


def do_refine_pos(d,dx,dy,fit_rad=10.):
    """
    """
    s = fit_rad/2.35
    s2 = 2*s*sqrt(2./pi)
    s2s = s*sqrt(2.)

    sz = len(d)
    sz0 = (sz-1)//2
    xx = linspace(-sz0-0.5,sz0+0.5,2*sz0+2)

    x1,y1 = xx - dx, xx - dy
    rad2 = ((x1[1:]-0.5)**2)[:,newaxis] + ((y1[1:]-0.5)**2)[newaxis,:]
    i1 = (rad2<=(fit_rad-0.5)**2).astype('float64')
    i2 = (rad2<=(fit_rad+0.5)**2).astype('float64')
    dhr = d*(i1+i2)/2.

    # sub pixel weighting
    wx,wy = -exp(-0.5*(x1/s)**2),-exp(-0.5*(y1/s)**2)
    wx,wy = wx[1:]-wx[:-1], wy[1:]-wy[:-1]

    # bias it back to zero a bit in case there is a nearby source
    h = x1[1:]-0.5<-fit_rad/2.
    if (h.sum()>0): wx[h]-=2*wx[h].mean()
    h = x1[1:]-0.5>fit_rad/2.
    if (h.sum()>0): wx[h]-=2*wx[h].mean()
    h = y1[1:]-0.5<-fit_rad/2.
    if (h.sum()>0): wy[h]-=2*wy[h].mean()
    h = y1[1:]-0.5>fit_rad/2.
    if (h.sum()>0): wy[h]-=2*wy[h].mean()

    cwx,cwy = erf(x1/s2s),erf(y1/s2s)
    cwx,cwy = cwx[1:]-cwx[:-1], cwy[1:]-cwy[:-1]

    dwx,dwy = dot(cwx,dhr),dot(dhr,cwy)
    norm = dot(cwx,dwy)

    deltax,deltay = 0.,0.
    if (norm>0):
        dxa = dot(wx,dwy) / norm
        dya = dot(dwx,wy) / norm
        deltax = dxa*s2; deltay = dya*s2

    return dx + deltax, dy + deltay


def refine_pos(dat,hdr,ra,dec,x,y,fit_rad=10.,niter=5):
    """
        dat: input image
        hdr: input image header
        refines x,y,ra,dec using dat
    """
    sx,sy = dat.shape
    sz0 = int(ceil(2*fit_rad))

    #
    i,j = y.round().astype('int16')-1, x.round().astype('int16')-1
    dx,dy = y-i-1,x-j-1

    # figure out which ones we just want to punt on
    h = (i-sz0>=0)*(i+sz0<sx)*(j-sz0>=0)*(j+sz0<sy)
    ii = where(h)[0]
    n1 = len(ii)

    for k0 in range(n1):
        k = ii[k0]

        # postage stamp data d
        d = dat[i[k]-sz0:i[k]+sz0+1,j[k]-sz0:j[k]+sz0+1].copy()

        # remove any linear trend in the psf wings/background
        a,b = d.shape
        xa = arange(a,dtype='float64'); xa-=xa.mean()
        xb = arange(b,dtype='float64'); xb-=xb.mean()
        xx = xa[:,newaxis] + zeros(b)[newaxis,:]
        yy = xb[newaxis,:] + zeros(a)[:,newaxis]
        rad2 = (xa**2)[:,newaxis] + (xb**2)[newaxis,:]
        hr = rad2>fit_rad**2
        if (h.sum()>0):
            d0 = d[hr].mean()
            d -= xb*dot(d[hr]-d0,yy[hr])/dot(yy[hr],yy[hr])
            dT = d.T
            dT -= xa*dot(d[hr]-d0,xx[hr])/dot(xx[hr],xx[hr])

        for l in range(niter):
            dx[k],dy[k] = do_refine_pos(d,dx[k],dy[k],fit_rad=fit_rad)

    x1,y1 = dy+j+1,dx+i+1
    h1 = ( abs(x1-x)<1 )*( abs(y1-y)<1 )*h

    ra0,dec0 = xy2ad(x[h1],y[h1],hdr)
    x[h1],y[h1] = x1[h1],y1[h1]
    ra1,dec1 = xy2ad(x[h1],y[h1],hdr)
    ra[h1] += ra1-ra0; dec[h1] += dec1-dec0

    return 1


if __name__ == '__main__':
    """
    """
    if (len(sys.argv)<2): usage()

    fitsfile=sys.argv[1]

    fit_rad=10.
    if (len(sys.argv)>2): fit_rad=float(sys.argv[2])

    dirname=os.path.dirname(fitsfile) or '.'
    base=os.path.basename(fitsfile).replace('.fits','')
    xyfile=base+'_dir/'+base+'_radec.txt'
    if (not os.path.exists(xyfile)): xyfile=base+'_radec.txt'

    (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx)=loadtxt(xyfile,unpack=True,ndmin=2)

    dat,hdr = getdata(fitsfile,header=True)
    stat = refine_pos(dat,hdr,ra,dec,x,y,fit_rad=fit_rad)

    f = open(xyfile,'r')
    str0 = f.readline().strip()
    f.close()

    print (str0)
    for i in range(len(ra)):
        print ("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d""" % (ra[i],dec[i],mag[i],dmag[i],mag_big[i],dmag_big[i],fwhm[i],x[i],y[i],xa[i],ya[i],x2a[i],y2a[i],expos[i],idx[i]))
