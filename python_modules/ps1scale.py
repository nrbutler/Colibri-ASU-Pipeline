#!/usr/bin/python3
"""
 ps1scale.py infile
"""
import os, sys
from astropy.io.fits import getdata,getheader,writeto
from numpy import sinh,log,isnan

def usage():
    print (__doc__)
    sys.exit()


def ps1scale(in_file):
    """
      apply the asinh scaling for ps1 files
    """
    out_file = in_file.replace('.fits.fz','.fits')
    if (out_file!=in_file):
        os.system(f"funpack {in_file}")
        in_file=out_file

    dat,hdr = getdata(in_file,header=True)
    nn = isnan(dat)
    dat[nn] = 0.

    gain0 = 1.
    gain_word = 'HIERARCH CELL.GAIN'
    if ( gain_word in hdr ): gain0 = hdr[gain_word]

    if ('BSOFTEN' in hdr):
        a = 2.5/log(10)
        dat[~nn] = hdr['BOFFSET'] + 2*hdr['BSOFTEN']*sinh(dat[~nn]/a)
        hdr['FLXSCALE'] = 1/hdr['EXPTIME']
        hdr['GAIN'] = gain0*hdr['EXPTIME']
        del hdr['BOFFSET']
        del hdr['BSOFTEN']
        writeto(out_file,dat,hdr,overwrite=True)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<2): usage()

    in_file= sys.argv[1]
    if (os.path.exists(in_file)==False): usage()

    ps1scale(in_file)
