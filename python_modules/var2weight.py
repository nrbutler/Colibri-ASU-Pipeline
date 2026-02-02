#!/usr/bin/python3
"""
 var2weight.py variance_image weight_image
"""
import os, sys
from astropy.io.fits import getdata,getheader,writeto
from numpy import sinh,log,isnan

def usage():
    print (__doc__)
    sys.exit()


def var2weight(var_file,weight_file):
    """
    Take a sextractor exposure variance file and transform to weight (invert)
    """

    print ("""Making weight image: %s from variance image %s""" % (weight_file,var_file))

    dat,hdr = getdata(var_file,header=True)
    nn = isnan(dat)
    dat[nn] = 0.

    if ('BSOFTEN' in hdr):
        a = 2.5/log(10)
        dat[~nn] = hdr['BOFFSET'] + 2*hdr['BSOFTEN']*sinh(dat[~nn]/a)
        del hdr['BOFFSET']
        del hdr['BSOFTEN']

    wt = 0.*dat
    h = dat>0
    wt[h] = 1/dat[h]

    writeto(weight_file,wt,hdr,overwrite=True)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<3): usage()

    var_file= sys.argv[1]
    if (os.path.exists(var_file)==False): usage()
    wt_file= sys.argv[2]

    var2weight(var_file,wt_file)
