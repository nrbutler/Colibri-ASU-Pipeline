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

def get_image_table(ra,dec,filters="grizy"):
    """
    Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    filters = string with filters to include. includes all by default
    Returns a table with the results
    """
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    # The final URL appends our query to the PS1 image service
    url = f"{service}?ra={ra}&dec={dec}&filters={filters}"
    # Read the ASCII table returned by the url
    table = Table.read(url, format='ascii')
    return table


def geturl(ra, dec, table, size=1024, output_size=None, filters="grizy", format="fits", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           f"ra={ra}&dec={dec}&size={size}&format={format}")
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


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

    i = [-1,0,1]
    files = []

    for i0 in i:
        for j0 in i:
            ra,dec = ij2ad(i0*pos_err,j0*pos_err,hdr)
            tbl = get_image_table(ra,dec,filter)
            files.append(geturl(ra,dec,tbl)[0].split('/')[-1])

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
        tbl = get_image_table(ra,dec,filter)
        url = geturl(ra,dec,tbl,size=sz)[0]

    #https://ps1images.stsci.edu/rings.v3.skycell/1495/075/rings.v3.skycell.1495.075.stk.r.unconv.fits
    #https://ps1images.stsci.edu/rings.v3.skycell/1495/076/rings.v3.skycell.1495.076.stk.r.unconv.fits

    print (url)
