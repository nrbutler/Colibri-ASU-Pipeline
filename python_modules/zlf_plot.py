#!/usr/bin/python3
"""
zlf_plot.py datafile
"""

import sys,os
import matplotlib as mpl
mpl.use('Agg')
from astropy.io.fits import getdata
from matplotlib.pyplot import figure,subplot,errorbar,savefig,subplots_adjust,scatter

from ut2gps import ut2gps

def zlf_plot(infile):
    """
    """
    # read in the summary of images reduced (ignoring stack)
    data=getdata(infile,1)
    t1,t2 = data['T0'],data['T1']

    i0 = t1.argmin()
    t1,t2=t1-t1[i0],t2-t1[i0]
    t10 = data['DATE-OBS'][i0]

    zp,lm,f=data['magzero'],data['maglim'],data['fwhm']

    t = 0.5*(t1+t2)/3600.
    dt = 0.5*(t2-t1)/3600.

    ff=figure(figsize=(8,6))
    ax0 = subplot(311); ax1 = subplot(312); ax2 = subplot(313)
    ax0.set_title(infile)

    dtm=dt.max()
    ax0.errorbar (t,zp,yerr=0*zp+0.01,xerr=dt,fmt='bo',capsize=0,linestyle='None',markersize=0,mew=1,alpha=0.2)
    ax0.scatter(t,zp,s=(dt/dtm)*200)
    ax1.errorbar (t,lm,yerr=0*lm+0.01,xerr=dt,fmt='ro',capsize=0,linestyle='None',markersize=0,mew=1,alpha=0.2)
    ax1.scatter(t,lm,s=(dt/dtm)*200)
    ax2.errorbar (t,f,yerr=0*f+0.01,xerr=dt,fmt='go',capsize=0,linestyle='None',markersize=0,mew=1,alpha=0.2)
    ax2.scatter(t,f,s=(dt/dtm)*200)

    ax0.set_ylabel("Zero Point",fontsize=14)
    ax1.set_ylabel("10-sigma LM",fontsize=14)
    ax2.set_xlabel("Hours Since %s" % t10,fontsize=16); ax2.set_ylabel("FWHM [pixels]",fontsize=14)
    tmin,tmax=t.min(),t.max()
    ax0.set_xlim((tmin,tmax))
    ax1.set_xlim((tmin,tmax))
    ax2.set_xlim((tmin,tmax))

    ax0.set_yticks(ax0.get_yticks()[1:-1])
    ax1.set_yticks(ax1.get_yticks()[1:-1])
    ax2.set_yticks(ax2.get_yticks()[1:-1])
    ax0.set_xticklabels([])
    ax1.set_xticklabels([])

    subplots_adjust(hspace=0,top=0.95,bottom=0.1)
    savefig(infile+'.jpg')


if __name__ == "__main__":
    """
    """
    file=sys.argv[1]
    zlf_plot(file)
