#!/usr/bin/python3
"""
 gaussfilt.py image <fwhm>
"""
import os, sys
from astropy.io.fits import getdata,writeto
from scipy.ndimage import gaussian_filter

def usage():
    print (__doc__)
    sys.exit()


def gaussfilt(data_file,outfile='gauss.fits',fwhm=2.):
    """
    convolve an image with a gaussian filter
    """
    dat,hdr = getdata(data_file,header=True)

    dat1 = gaussian_filter(dat,fwhm/2.3548)
    writeto(outfile, dat1, hdr, overwrite=True)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<2): usage()

    data_file = sys.argv[1]
    if (os.path.exists(data_file)==False): usage()

    fwhm=2.0
    if (len(sys.argv)>2): fwhm=float(sys.argv[2])

    dir = os.path.dirname(data_file) or "."
    file0 = os.path.basename(data_file)
    outfile = dir + '/gauss_' + file0

    gaussfilt(data_file,fwhm=fwhm,outfile=outfile)
