#!/usr/bin/python3
"""
   fits_resize.py fitsfile [factor]
"""

import sys,os
from astropy.io.fits import getdata,writeto
from numpy import zeros

def usage():
    print (__doc__)
    sys.exit()

def fits_resize(x,wx,hdr,filename,fac=15):
    """
       create thumbnails from fits files and weights
    """
    fac = int(round(fac))

    a,b = x.shape
    a1,b1 = a//fac,b//fac
    dx,dy = a-a1*fac,b-b1*fac
    dx1,dy1 = dx//2,dy//2
    dx2,dy2 = dx-dx1,dy-dy1

    cx = zeros((a+1,b+1),dtype='float64')

    # rebin the stack*weight
    cx[1:,1:] = (x*wx).astype('float64').cumsum(axis=0).cumsum(axis=1)
    dcx = 0.25*( cx[dx1::fac,dy1::fac] + cx[dx1::fac,dy2::fac] + cx[dx2::fac,dy1::fac] + cx[dx2::fac,dy2::fac] )
    dcx = dcx[1:] - dcx[:-1]
    x1 = ( (dcx[:,1:]-dcx[:,:-1])/(fac*fac) ).astype('float32')

    # rebin the weight
    cx[1:,1:] = wx.astype('float64').cumsum(axis=0).cumsum(axis=1)
    dcx = 0.25*( cx[dx1::fac,dy1::fac] + cx[dx1::fac,dy2::fac] + cx[dx2::fac,dy1::fac] + cx[dx2::fac,dy2::fac] )
    dcx = dcx[1:] - dcx[:-1]
    wx1 = ( (dcx[:,1:]-dcx[:,:-1])/(fac*fac) ).astype('float32')
    cx=0

    h = wx1>1.e-5
    x1[h] /= wx1[h]
    x1[~h] = 0

    base=filename.replace('.fits','')
    file1=base+'_thumb.fits'
    file2=base+'_thumb.wt.fits'
    writeto(file1,x1,hdr,overwrite=True)
    writeto(file2,wx1,hdr,overwrite=True)

    cd11,cd12,cd21,cd22 = fac*hdr['CD1_1'],fac*hdr['CD1_2'],fac*hdr['CD2_1'],fac*hdr['CD2_2']
    x0,y0 = hdr['CRPIX1'],hdr['CRPIX2']
    a1,b1 = x1.shape
    x1 = (1+b1)/2. + (x0-(1+b)/2.)/fac
    y1 = (1+a1)/2. + (y0-(1+a)/2.)/fac

    str="""sethead CD1_1=%.8f CD1_2=%.8f CD2_1=%.8f CD2_2=%.8f CRPIX1=%.4f CRPIX2=%.4f %s %s""" % (cd11,cd12,cd21,cd22,x1,y1,file1,file2)
    os.system(str)


if __name__ == "__main__":
    """
     create thumbnails from fits files and weights
    """
    if (len(sys.argv)<2): usage()
    infile=sys.argv[1]

    fac=30
    if (len(sys.argv)>2): fac=int(sys.argv[2])

    x,hdr = getdata(infile,header=True)
    winfile = infile.replace('.fits','.wt.fits')
    wx = getdata(winfile)

    fits_resize(x,wx,hdr,infile,fac=fac)
