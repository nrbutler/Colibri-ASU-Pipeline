#!/usr/bin/python3
"""
get_ps1imagelist.py ra dec <sz> <filter>
"""

def usage():
    print (__doc__)
    sys.exit()

import sys,os
from astropy.table import Table
from numpy import argsort,unique,sqrt
from astropy.io.fits import getheader
from fit_wcs import ij2ad

from findskycell import findskycell

def geturl(ra,dec,filter='r'):
    tbl = findskycell(ra,dec)
    pcell = tbl['projcell'][0]
    scell = tbl['subcell'][0]
    file = f"rings.v3.skycell.{pcell:04d}.{scell:03d}.stk.{filter}.unconv.fits"
    return file

def getfile_urls(fits_file,ps=0.38/3600):
    """
      pull info from a fits header
    """
    if (fits_file.replace('.fits.fz','.fits')==fits_file):
        hdr = getheader(fits_file)
    else:
        hdr = getheader(fits_file,1)

    try:
        ra0, dec0 = hdr['CRVAL1'],hdr['CRVAL2']
    except:
        ra0, dec0 = hdr['ETROBRA'],hdr['ETROBDE']
        hdr['CRVAL1'],hdr['CRVAL2'] = ra0, dec0
        hdr['CD1_1'] = -ps
        hdr['CD1_2'] = 0
        hdr['CD2_1'] = 0
        hdr['CD2_2'] = ps

    pos_err = max(hdr['NAXIS1'],hdr['NAXIS2'])*ps/2.

    try:
        filter = hdr['FILTER']
    except:
        filter = 'r'

    if (filter=='zy'): filter='z'
    if (filter=='gri'): filter='r'
    if (filter=='B'): filter='g'

    i = [-1,0,1]
    files = []

    for i0 in i:
        for j0 in i:
            ra,dec = ij2ad(i0*pos_err,j0*pos_err,hdr)
            files.append(geturl(ra,dec,filter))

    return " ".join(list(unique(files)))

            
if __name__ == "__main__":

    if (len(sys.argv)<3): usage()
 
    ra = float(sys.argv[1])
    dec = float(sys.argv[2])

    sz=1024
    if (len(sys.argv)>3): sz=int(sys.argv[3])
    filter="r"
    if (len(sys.argv)>4): filter=sys.argv[4]

    infile="null_nofile.txt"
    if (len(sys.argv)>5): infile=sys.argv[5]
 
    if (os.path.exists(infile)):
        url = getfile_urls(infile) 
    else:
        url = geturl(ra,dec,filter)

    #https://ps1images.stsci.edu/rings.v3.skycell/1495/075/rings.v3.skycell.1495.075.stk.r.unconv.fits
    #https://ps1images.stsci.edu/rings.v3.skycell/1495/076/rings.v3.skycell.1495.076.stk.r.unconv.fits

    print (url)
