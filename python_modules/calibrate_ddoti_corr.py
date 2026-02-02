#!/usr/bin/python3
"""
calibrate.py catalog.fits calibration.fits [doplot]
   for comparing a sextractor radec list to catalogs
"""

import sys,os
from astropy.io.fits import getheader,getdata,writeto
from quick_mode import quick_mode
from robust_mean import robust_mean
from numpy import sqrt,log10,median,zeros,sign
from linfit import linfit

import matplotlib as mpl
mpl.use('Agg')
from matplotlib.pyplot import subplots,ioff,close
ioff()

def usage():
    print (__doc__)
    sys.exit()

def calibrate(infile,calfile,sys_err=0.01,cal_mag_min=14.0,cal_mag_max=19.0,do_plot=False,quick_cal=False,mag0=0.):
    """
    sys_err used only to down-weight very bright calibration sources
    """
    infile0=os.path.basename(infile)
    dirname=os.path.dirname(infile) or "."

    phot_dat=getdata(infile)

    wtfile=infile+'.apweight'
    if (os.path.exists(wtfile)):
        wdata = getdata(wtfile)
        apc,apcb = wdata[:,2],wdata[:,3]
        good = apc>0; apcg=apc[good]
        phot_dat['FLUX_APER'][good,0] /= apcg
        phot_dat['FLUXERR_APER'][good,0] /= apcg
        phot_dat['FLUX_APER'][~good,0] = 0
        phot_dat['FLUXERR_APER'][~good,0] = 0
        good = apcb>0; apcg=apcb[good]
        phot_dat['FLUX_APER'][good,1] /= apcg
        phot_dat['FLUXERR_APER'][good,1] /= apcg
        phot_dat['FLUX_APER'][~good,1] = 0
        phot_dat['FLUXERR_APER'][~good,1] = 0

    if (quick_cal):
        fac = 10**(-0.4*mag0)
        phot_dat['FLUX_APER'][:,0] *= fac
        phot_dat['FLUX_APER'][:,1] *= fac
        phot_dat['FLUXERR_APER'][:,0] *= fac
        phot_dat['FLUXERR_APER'][:,1] *= fac
        ofile=dirname+'/calibrated_'+infile0
        writeto(ofile,phot_dat,overwrite=True)
        return 0

    flx=phot_dat['FLUX_APER'][:,0]
    flx_big=phot_dat['FLUX_APER'][:,1]
    dflx=phot_dat['FLUXERR_APER'][:,0]

    # calibration data
    cdata=getdata(calfile)
    #flx0 = cdata['FLUX_APER'][:,1]
    #dflx0 = cdata['FLUXERR_APER'][:,1]
    flx0 = cdata['FLUX_APER'][:,0]
    dflx0 = cdata['FLUXERR_APER'][:,0]

    hdr = getheader(infile,1)
    if ('FWHM' not in hdr): hdr = getheader(infile)
    fwhm=hdr.get('FWHM',1.0)

    hdr = getheader(infile,1)
    if ('GAIN' not in hdr): hdr = getheader(infile)
    gain=hdr.get('GAIN',132.)
    dt=hdr.get('EXPTIME',60.)
    am=hdr.get('AMEXT',0.)
    sex_zero=hdr.get('SEXZERO',25.)
    t0=hdr.get('T0',0.)
    t1=hdr.get('T1',0.)

    hdrc = getheader(calfile)
    if ('CFILTER' not in hdrc): hdrc = getheader(calfile,1)
    cfilter=hdrc.get('CFILTER','unknown system')

    flx0_min,flx0_max = 10**(-0.4*(cal_mag_max-25)),10**(-0.4*(cal_mag_min-25))
    gphot = (flx0>flx0_min)*(flx0<=flx0_max)*(dflx0>0)*(flx>100*dflx)*(dflx>0)*(phot_dat['FWHM_IMAGE']<10*fwhm)
    if (gphot.sum()<10):
        print ("Not enough calibration sources for snr>100, falling back to snr>10")
        gphot = (flx0>flx0_min)*(flx0<=flx0_max)*(dflx0>0)*(flx>10*dflx)*(dflx>0)*(phot_dat['FWHM_IMAGE']<10*fwhm)

    mag_all = 25-2.5*log10(flx.clip(1))
    dmag_all0 = 1.0857*dflx/flx.clip(1)
    dmag_all = sqrt( dmag_all0**2 + sys_err**2 )

    mag = mag_all[gphot]
    mag_big = 25 - 2.5*log10(flx_big[gphot].clip(1))
    dmag = dmag_all[gphot]
    mag0 = 25-2.5*log10(flx0[gphot])
    dmag0 = sqrt( (1.0857*dflx0[gphot]/flx0[gphot])**2 + sys_err**2 )

    nn=gphot.sum()
    if (nn>0):

        print (" Matches %d" % nn)
        print (""" Time Range %s - %s""" % (t0,t1))

        print (""" Median Sextractor FWHM %.2f""" % fwhm)

        gcal_file='good_cal.fits'
        if (os.path.exists(gcal_file)):
            gcal_indx = getdata(gcal_file).astype('int64')
            gphot0 = zeros(len(gphot),dtype='bool')
            gphot0[gcal_indx] = True
            good_cal = gphot0[gphot]
            if (good_cal.sum()<10): good_cal = gphot[gphot]
        else: good_cal = gphot[gphot]

        mag_offset = mag0[good_cal]-mag[good_cal]-dmag[good_cal]/2
        dmag_offset = sqrt(dmag0[good_cal]**2+dmag[good_cal]**2)
        mag_offset, dmag_offset = robust_mean(mag_offset,dmag_offset)

        print (""" Magnitude Offset %.4f +/- %0.4f""" % (mag_offset,dmag_offset))

        fac = 10**(-0.4*mag_offset)

        f = phot_dat['FLUX_APER'][:,0]; sf = sign(f)
        phot_dat['FLUX_APER'][:,0] = sf*(abs(f) - phot_dat['FLUXERR_APER'][:,0]/2).clip(0)*fac
        phot_dat['FLUXERR_APER'][:,0] *= fac

        f = phot_dat['FLUX_APER'][:,1]; sf = sign(f)
        phot_dat['FLUX_APER'][:,1] = sf*(abs(f) - phot_dat['FLUXERR_APER'][:,1]/2).clip(0)*fac
        phot_dat['FLUXERR_APER'][:,1] *= fac

        flx=phot_dat['FLUX_APER'][:,0]
        mag_all = 25-2.5*log10(flx.clip(1)) - mag_offset
        dmag_all0 = 1.0857*dflx/flx.clip(1)
        dmag_all = sqrt( dmag_all0**2 + sys_err**2 )
        mag = mag_all[gphot]
        dmag = dmag_all[gphot]

        # calculate the aperture correction
        #ap_corr = 0.
        ap_corr = quick_mode(mag_big-mag)
        zero_pt = sex_zero + mag_offset + 2.5*log10(gain/dt)+am - ap_corr
        print (""" Median Zero Point %.3f [gain=%.2f, dt=%.1f, am_corr=%.3f]""" % (zero_pt,gain,dt,am))
        print (""" Sextractor Aperture Correction (mag) %.3f""" % ap_corr)

        gl = (mag_all<25)*(dmag_all0<0.1)
        mag_all += mag_offset
        mag10s = quick_mode( mag_all[gl]-2.5*log10(9.21*dmag_all0[gl]) )
        print (" 10-sigma limiting magnitude %.2f" % mag10s)

        ofile=dirname+'/calibrated_'+infile0
        writeto(ofile,phot_dat,overwrite=True)

        if (do_plot):
            fig,ax = subplots()
            x0=max(cal_mag_min,mag0.min())-0.1
            x1=min(cal_mag_max,mag0.max())+0.1
            ax.errorbar (mag0,mag0-mag_all[gphot],yerr=sqrt(dmag**2+dmag0**2),xerr=dmag0,fmt='bo',capsize=0,linestyle='None',markersize=5,mew=1,alpha=0.1)
            ax.plot (mag0,mag0-mag_all[gphot],'ro',alpha=0.1)
            ax.plot ([x0,x1],[0.,0.],':')
            ax.set_xlim((x0,x1))
            #ax.set_ylim((-0.5,0.5))
            ax.set_ylim((-0.2,0.2))
            ax.set_xlabel("""Catalog Mag [%s]""" % cfilter,fontsize=14);
            ax.set_ylabel("Catalog Mag - Mag",fontsize=14);
            ax.set_title("""10-sigma: %.2f ; FWHM=%.1f arcsec""" % (mag10s,fwhm*2.01))
            fig.savefig(ofile.replace('.fits','.jpg'))
            close(fig)

            fig,ax = subplots()
            ra = phot_dat['ALPHA_J2000'][gphot]
            x0,x1 = ra.min(),ra.max()
            ax.errorbar (ra,mag0-mag_all[gphot],yerr=sqrt(dmag**2+dmag0**2),xerr=0*dmag0,fmt='bo',capsize=0,linestyle='None',markersize=5,mew=1,alpha=0.1)
            ax.plot (ra,mag0-mag_all[gphot],'ro',alpha=0.1)
            ax.plot ([x0,x1],[0.,0.],':')
            ax.set_xlim((x0,x1))
            #ax.set_ylim((-0.5,0.5))
            ax.set_ylim((-0.2,0.2))
            ax.set_xlabel("Catalog Right Ascension",fontsize=14)
            ax.set_ylabel("Catalog Mag - Mag",fontsize=14);
            ax.set_title("""10-sigma: %.2f ; FWHM=%.1f arcsec""" % (mag10s,fwhm*2.01))
            fig.savefig(ofile.replace('.fits','_ra.jpg'))
            close(fig)

            fig,ax = subplots()
            dec = phot_dat['DELTA_J2000'][gphot]
            x0,x1 = dec.min(),dec.max()
            ax.errorbar (dec,mag0-mag_all[gphot],yerr=sqrt(dmag**2+dmag0**2),xerr=0*dmag0,fmt='bo',capsize=0,linestyle='None',markersize=5,mew=1,alpha=0.1)
            ax.plot (dec,mag0-mag_all[gphot],'ro',alpha=0.1)
            ax.plot ([x0,x1],[0.,0.],':')
            ax.set_xlim((x0,x1))
            #ax.set_ylim((-0.5,0.5))
            ax.set_ylim((-0.2,0.2))
            ax.set_xlabel("Catalog Declination",fontsize=14)
            ax.set_ylabel("Catalog Mag - Mag",fontsize=14);
            ax.set_title("""10-sigma: %.2f ; FWHM=%.1f arcsec""" % (mag10s,fwhm*2.01))
            fig.savefig(ofile.replace('.fits','_dec.jpg'))
            close(fig)

    else:
        sys.stderr.write("""No matches for %s!\n""" % infile)

    return 1


if __name__ == "__main__":

    if (len(sys.argv)<3): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==0): usage()

    mag0=0.
    quick_cal=False
    calfile=sys.argv[2]
    if (os.path.exists(calfile)==0):
        if ("mag0" in calfile): 
            mag0 = float( calfile.split('=')[-1] )
            quick_cal=True
        else: usage()

    do_plot=False
    if (len(sys.argv)>3):
        if (sys.argv[3]=="doplot"): do_plot=True

    stat = calibrate(infile,calfile,do_plot=do_plot,quick_cal=quick_cal,mag0=mag0)
