#!/usr/bin/python3
"""
 ut2gps.py ut_time
"""

from numpy import array,empty,atleast_1d
from time import mktime

import sys

def usage():
    print (__doc__)
    sys.exit()

def ut2gps(ut):
    """
    crude transformation between UT and GPS times
     eg., ut='20130504T033516' or '20130504_033516'
     gps = ut2gps(ut)
    """
    uta = atleast_1d(array(ut))
    gps = empty(len(uta),dtype='float64')

    t0 = mktime( (1980, 1, 5, 23, 59, 47,1,1,-1) )

    for i in range(len(uta)):
        ut0 = uta[i]
        yr,mo,da =int(ut0[0:4]), int(ut0[4:6]), int(ut0[6:8])
        hr,mi,se =int(ut0[9:11]), int(ut0[11:13]), float(ut0[13:])
        dse = se - int(se)
        gps[i] = mktime( (yr,mo,da,hr,mi,int(se),1,1,-1) ) - t0 + dse + 5

    return gps


if __name__ == "__main__":
    """
    just run
    """
    if (len(sys.argv)<2): usage()
    ut = sys.argv[1]
    gps = ut2gps(ut)[0]
    print (f"{gps:.3f}")



