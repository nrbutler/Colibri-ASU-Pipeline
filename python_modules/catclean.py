#!/usr/bin/python3
"""
 catclean.py catalog.fits seg.fits stack.fits
"""

import os,sys
from astropy.io.fits import open,getdata,getheader,writeto
from numpy import arange,zeros,sqrt

def usage():
    print (__doc__)
    sys.exit()

#def catclean(catalog,segfile,stackfile,snr_min=4.5,snr_min0=3.0,err_max=13.1579):
def catclean(catalog,segfile,stackfile,snr_min=4.5,snr_min0=3.0,err_max=52.6315):
    """ clean a catalog file , rejecting sources with snr<snr_min
        no rejection is applied above snr_min0 within error region
        if error region (radius) is smaller than err_max (pixels)
    """

    hdu = open(catalog)

    flx = hdu[1].data['FLUX_APER'][:,0]
    dflx = hdu[1].data['FLUXERR_APER'][:,0]

    good = flx>snr_min*dflx

    hdr = getheader(stackfile)
    if ('SEXERR' in hdr):
        x0,y0,err = hdr['SEXX0'],hdr['SEXY0'],hdr['SEXERR']
        if (err<=err_max):
            dis = sqrt( (hdu[1].data['X_IMAGE']-x0)**2+(hdu[1].data['Y_IMAGE']-y0)**2 )
            ng0 = good.sum()
            good[ (dis<err_max)*(flx>snr_min0*dflx) ] = True
            ng = good.sum()
            print (f"Retaining {ng-ng0} sources with {snr_min0}<snr<={snr_min} within {err_max} pixels of x,y=({x0},{y0})")
    else:
        print ("No error region found in stack header")
    
    hdu[1].data = hdu[1].data[good]
    hdu.writeto(catalog,overwrite=True)

    seg,hdr = getdata(segfile,header=True)
    Ng,Ng0 = good.sum(),len(good)
    #old_ids = arange(Ng0)
    new_ids = zeros(Ng0,dtype='int64')
    new_ids[good] = 1 + arange(Ng)
 
    h = seg>0
    seg[h] = new_ids[seg[h]-1]

    writeto(segfile,seg,hdr,overwrite=True)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<4): usage()

    catalog=sys.argv[1]
    if (os.path.exists(catalog)==False): usage()

    segfile=sys.argv[2]
    if (os.path.exists(segfile)==False): usage()

    stackfile=sys.argv[3]
    if (os.path.exists(stackfile)==False): usage()

    snr_min=5.0
    if (len(sys.argv)>4): snr_min = float(sys.argv[4])

    catclean(catalog,segfile,stackfile,snr_min=snr_min)
