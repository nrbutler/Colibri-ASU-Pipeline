#!/usr/bin/python3
"""
 ast_err_plot.py
"""
import sys,os
from numpy import loadtxt,sqrt,log10

import matplotlib as mpl
mpl.use('Agg')

from glob import glob
from astropy.io.fits import getdata
from matplotlib.pyplot import boxplot,xlabel,ylabel,savefig

def usage():
    print (__doc__)
    sys.exit()


def ast_err_plot(data):
    """
       make a box plot of the positional offsets
    """

    boxplot(data)
    xlabel("Image Number",fontsize=14)
    ylabel("Star Offsets [arcsec]",fontsize=14)

    savefig("boxplot.jpg")


if __name__ == "__main__":
    """
    make a boxplot for the astrometry
    """
    files=glob('20*fits.matches')
    files.sort()

    data=[]
    for file in files:
        data.append( getdata(file) )

    if (len(data)>0): ast_err_plot(data)
