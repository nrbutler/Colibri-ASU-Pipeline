#!/usr/bin/python3
"""
 radec2sex.py ra dec
"""
import os, sys
from fit_wcs import dec2sex,ra2sex

def usage():
    print (__doc__)
    sys.exit()


def radec2sex(ra,dec):
    """
      convert decimal to sexagesimal
    """
    ras = ra2sex(ra)
    decs = dec2sex(dec)

    return ras,decs


if __name__ == '__main__':
    """
    """
    if (len(sys.argv)<3): usage()

    ra = float(sys.argv[1])
    dec = float(sys.argv[2])

    ras, decs = radec2sex(ra,dec)
    print (ras,decs)
