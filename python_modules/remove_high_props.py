#!/usr/bin/python3
"""
remove_high_props.py radecfile.txt high_prop_file.fits

   remove high proper motion stars
"""
import sys,os
from numpy import loadtxt,ones,sqrt,abs,pi,cos
from astropy.io.fits import getdata

def usage():
    print (__doc__)
    sys.exit()


def remove_high_props(infile='radecfile.txt',propfile='high_prop_file.fits'):
    """
    """

    phot_dat=loadtxt(infile,ndmin=2).T

    if (len(phot_dat[0])>0):

        r,d = phot_dat[:2]
        idx = phot_dat[-1]

        nr = len(r)
        good = ones(nr,dtype='bool')

        r0,d0,m,dm,dr0,dd0 = getdata(propfile)
        dis0=sqrt(dd0**2+dr0**2)*3600.

        for i in range(nr):
            if (idx[i]>0): continue
            cd=cos(d[i]*pi/180.) 
            dis = sqrt( ((r[i]-r0)*cd)**2+(d[i]-d0)**2 )*3600.
            dis1= abs( dd0*(r[i]-r0)*cd-dr0*(d[i]-d0) )/sqrt( dd0**2+dr0**2 )*3600.
            h = (dis<10*dis0)*(dis1<dis0)  # require a match for a possible bad pm estimate
            if (h.sum()>0): good[i] = False

        nr1 = good.sum()
        sys.stderr.write("""Removing %d detections near high-proper-motion stars...\n""" % (nr-nr1))

        of=open(infile,'r')
        l1=of.readline()
        l2=of.readline()
        if (l2[0]=='#'): l1+=l2
        of.close()

        of=open(infile,'w')
        of1=open('pm_'+infile,'w')

        phot_dat_pm = phot_dat[:,~good]
        phot_dat = phot_dat[:,good]
        if (l1[0]=='#'): of.write(l1)
        for i in range(len(phot_dat[0])):
            (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx) = phot_dat[:,i]
            of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx))

        for i in range(len(phot_dat_pm[0])):
            (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx) = phot_dat_pm[:,i]
            of1.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx))

        of.close()
        of1.close()


if __name__ == "__main__":

    if (len(sys.argv)<3): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==0): usage()
    propfile=sys.argv[2]
    if (os.path.exists(propfile)==0): usage()

    remove_high_props(infile,propfile)
