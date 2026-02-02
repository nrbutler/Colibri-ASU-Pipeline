#!/usr/bin/python3
"""
 sources2mask.py wt_file sources_file <outfile> <mask_fwhm>
  [requires mask.fits from sextractor]
"""
import os, sys
from astropy.io.fits import getdata,writeto
from scipy.ndimage import gaussian_filter
from scipy.signal import medfilt2d
from numpy import median,loadtxt,atleast_1d,sqrt

def usage():
    print (__doc__)
    sys.exit()


def sources2mask(wt_file,source_file,outfile='mask.fits',fwhm=7.):
    """
      Turn a sextractor mask.fits file into a 0,1 mask (0==masked),
    """

    print ("""Making object mask for: %s""" % source_file)

    if (fwhm<1): fwhm = 1.
    sigma = fwhm/2.35

    # read in the weightmap
    #wt=getdata(wt_file)
    wt,hdr=getdata(wt_file,header=True)
    wt_type = wt_file.split('.')[-2]
    if (wt_type=='rms'):
        dat_file = wt_file.replace('.rms.fits','.fits')
        rms = wt
    else:
        dat_file = wt_file.replace('.wt.fits','.fits')
        rms = 0.*wt
        h = wt>0
        rms[h] = 1./sqrt(wt[h])

    # make the source mask, from the sextractor source_file
    #mask0,hdr = getdata(source_file,header=True)
    mask0 = getdata(source_file)
    if (sigma>0):
        mask = gaussian_filter(mask0,sigma)
        mask[(mask0==0)*(mask<rms)]=0.
        mask[mask0>0]=1
        mask0 = 0.
    else:
        mask = mask0

    mask+=1; mask[mask!=1]=0

    sex_dir=wt_file.replace('.wt.fits','_dir/')
    sex_file=os.path.basename(wt_file).replace('.wt.fits','_radec.txt')
    sex_file=sex_dir+sex_file
    if (os.path.exists(sex_file)==False):
        sex_dir=wt_file.replace('.rms.fits','_dir/')
        sex_file=os.path.basename(wt_file).replace('.rms.fits','_radec.txt')
        sex_file=sex_dir+sex_file

    if (os.path.exists(sex_file) and fwhm>=1):
        print ("""Using sex_file %s""" % sex_file)

        try: dm,x,y = loadtxt(sex_file,unpack=True,usecols=(3,7,8))
        except: dm,x,y=0.1,0.,0.

        dm,x,y = atleast_1d(dm).clip(1.e-4), atleast_1d(x), atleast_1d(y)
        sx,sy = mask.shape
        for i in range(len(x)):
            y0,x0 = int(round(x[i])),int(round(y[i]))
            rad = int(fwhm**0.75*dm[i]**(-0.25))
            for j in range(2*rad):
                j0 = x0-rad+j
                if (j0>=0 and j0<sx):
                    for k in range(2*rad):
                        k0 = y0-rad+k
                        rr = (k-rad)**2+(j-rad)**2
                        if (k0>=0 and k0<sy and rr<=rad*rad): mask[j0,k0]=0

    if (sigma>0):
        mask = gaussian_filter(mask,sigma)
        mask[mask>0.9]=1
        mask[mask<=0.9]=0

    mask[wt==0]=0

    if (os.path.exists(outfile)): os.remove(outfile)
    writeto(outfile,mask,hdr)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<3): usage()

    wtfile= sys.argv[1]
    if (os.path.exists(wtfile)==False): usage()
    sourcesfile= sys.argv[2]
    if (os.path.exists(sourcesfile)==False): usage()

    maskfile='mask.fits'
    if (len(sys.argv)>=3): maskfile=sys.argv[3]

    fwhm=3.0
    if (len(sys.argv)>4): fwhm=float(sys.argv[4])

    sources2mask(wtfile,sourcesfile,outfile=maskfile,fwhm=fwhm)
