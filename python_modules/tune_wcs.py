#!/usr/bin/python3
"""
tune_wcs.py fitsfile radecfile.txt calfile.txt
"""

import sys,os

from numpy import sqrt,array,dot,cos,sin,arange,median,log10,vstack,loadtxt
from fit_wcs import ad2ij,xy2ij,xy2ad,ij2ad
from astropy.io.fits import getheader,getdata,writeto
from scipy.optimize import fmin
from nsigma_clip import nsigma_clip
from scipy.spatial import cKDTree


def usage():
     print (__doc__)
     sys.exit()


def tune_wcs(fits_file,radecfile,catfile):
    """
    tune up a wcs solution in one file to match the catfile
    """

    hdr = getheader(fits_file)

    # image positions
    if (radecfile.split('.')[-1]=='txt'):
        x,y = loadtxt(radecfile,unpack=True,usecols=(7,8))
    else:
        radec = getdata(radecfile)
        x,y = radec['X_IMAGE'], radec['Y_IMAGE']

    i,j = xy2ij(x,y,hdr,distort=True)

    # catalog positions
    if (catfile.split('.')[-1]=='txt'):
        ra0,dec0 = loadtxt(catfile,unpack=True,usecols=(0,1))
    else:
        ra0,dec0 = getdata(catfile)[:2,:]

    i0,j0 = ad2ij(ra0,dec0,hdr)

    # match
    xy = vstack((i,j)).T; xy0 = vstack((i0,j0)).T
    tre=cKDTree(xy0)
    dis,ii = tre.query(xy,k=1)
    i0,j0 = i0[ii],j0[ii]

    delta_i,delta_j = i-i0,j-j0
    k = (abs(delta_i)<1./60)*(abs(delta_j)<1./60)
    k[k] = nsigma_clip(delta_i[k],2) * nsigma_clip(delta_j[k],2)
        
    delta_i,delta_j = delta_i[k].mean(),delta_j[k].mean()
    i -= delta_i; j -= delta_j

    par=array([0.,0.,0.,0.,0.])

    def myfun(par):
        delti,deltj,a,b,th = par
        s,c = sin(th),cos(th)
        di = i0 - (1+a)*(c*i+s*j) - delti
        dj = j0 - (1+b)*(c*j-s*i) - deltj
        k = (abs(di)<1./60)*(abs(dj)<1./60)
        k[k] = nsigma_clip(di[k],3) * nsigma_clip(dj[k],3)
        return ( (di[k]*3600.)**2+(dj[k]*3600.)**2 ).mean()

    res = fmin(myfun,par,disp=0)
    res = fmin(myfun,res,disp=0)
    delti,deltj,a,b,th = res
    s,c = sin(th),cos(th)

    di = i0 - (1+a)*(c*i+s*j) - delti
    dj = j0 - (1+b)*(c*j-s*i) - deltj
    err = sqrt(di**2+dj**2)*3600.
    writeto(fits_file+'.matches',err,overwrite=True)

    crval1,crval2 = ij2ad(delti-delta_i,deltj-delta_j,hdr)

    matr = array([ [(1+a)*c,(1+a)*s],[-(1+b)*s,(1+b)*c] ] )
    c11,c12,c21,c22=hdr['CD1_1'],hdr['CD1_2'],hdr['CD2_1'],hdr['CD2_2']
    cd = array([ [c11,c12],[c21,c22] ])
    cd1 = dot(matr,cd)
    hdr['CD1_1'],hdr['CD1_2'],hdr['CD2_1'],hdr['CD2_2'],hdr['CRVAL1'],hdr['CRVAL2'] = cd1[0,0],cd1[0,1],cd1[1,0],cd1[1,1],crval1,crval2
    str="""sethead CD1_1=%.8f CD1_2=%.8f CD2_1=%.8f CD2_2=%.8f CRVAL1=%.8f CRVAL2=%.8f %s""" % (cd1[0,0],cd1[0,1],cd1[1,0],cd1[1,1],crval1,crval2,fits_file)
    os.system(str)


if __name__ == "__main__":
    """
    """
    if (len(sys.argv)<4): usage()

    fits_file=sys.argv[1]
    if (os.path.exists(fits_file)==False): usage()

    radecfile=sys.argv[2]
    if (os.path.exists(radecfile)==False): usage()

    catfile=sys.argv[3]
    if (os.path.exists(catfile)==False): usage()

    tune_wcs(fits_file,radecfile,catfile)
