#!/usr/bin/python3
"""
calibrate.py radec.txt catalog_radec.txt <phot_outfile.txt> <cat_outfile.txt> <cal_mag_max> <ab_offset> <do_plot>
   for comparing a sextractor radec list to catalogs
"""
import sys,os
from quick_mode import quick_mode
from robust_mean import robust_mean
from numpy import loadtxt,sqrt,abs,log10,arange,where,median,zeros,ones

def usage():
    print (__doc__)
    sys.exit()


def calibrate(infile='radec.txt',calfile='cat_radec.txt',outfile='photometry.txt',matchfile='match.txt',sys_err=0.0001,cal_mag_max=18.0,ab_offset=0.,do_plot=False):
    """
    """
    if (do_plot):
        import matplotlib as mpl
        mpl.use('Agg')
        from matplotlib.pyplot import plot,xlabel,ylabel,savefig,semilogy,title,legend,clf,scatter,semilogx,xlim,ylim,colorbar,hist,errorbar

    phot_dat=loadtxt(infile,ndmin=2).T
    idx = phot_dat[-1].astype('int32')
    gphot = (phot_dat[2]<999)*(phot_dat[3]<999)*(idx>0)*(phot_dat[2]>0)*(phot_dat[3]>0)
    idx = idx[gphot]
    mag,dmag = phot_dat[2:4,gphot]; mag_big,dmag_big = phot_dat[4:6,gphot]; fwhm = phot_dat[6,gphot]

    lr1='#\n'
    f=open(infile,'r')
    try:
        lr1=f.readline(); dat=lr1.split();
        gain, dt, am, sex_zero, t0, t1 = float(dat[-12]), float(dat[-10]), float(dat[-8]), float(dat[-6]), dat[-4], dat[-2]
    except:
        gain,dt,am=1.,1.,1.
        sex_zero = 25.0
        t0, t1 = "0.0","0.0"
    f.close()

    # insert x,y dependent systematic error
    h = dmag>0
    dmag[h] = sqrt( dmag[h]**2 + sys_err**2 )

    cdat=loadtxt(calfile,usecols=(2,3,4),ndmin=2).T
    cdat = cdat[:,(cdat[1]<999)*(cdat[2]<999)*(cdat[1]>0)*(cdat[2]>0)*(cdat[0]>0)]
    idx0 = cdat[0].astype('int32')

    # w matchlist refers to index in full star list
    # idx and idx0 don't necessarily have all those
    ii0 = zeros(max(idx0.max(),idx.max()),dtype='int32')-1
    ii0[idx0-1] = arange(len(idx0))
    w = ii0[idx-1]

    w0=w>=0
    nn = w0.sum()
    if (nn>0):

        print (" Matches %d" % nn)
        print (""" Time Range %s - %s""" % (t0,t1))

        w = w[w0]
        idx0 = idx0[w]
        mag0,dmag0 = cdat[1:,w]

        ii = where(gphot)[0]
        gphot[ii[~w0]] = False

        mag,dmag,mag_big,dmag_big,fwhm = mag[w0],dmag[w0],mag_big[w0],dmag_big[w0],fwhm[w0]

        fwhm0 = quick_mode(fwhm)
        print (""" Median Sextractor FWHM %.2f""" % fwhm0)

        wh = (mag0<=cal_mag_max)*(fwhm<fwhm0*2.)
        if (wh.sum()<5): wh = mag0<=cal_mag_max
        if (wh.sum()<5): wh = ones(len(mag0),dtype='bool')

        mag0 += ab_offset
        mag_offset = mag0-mag
        dmag_offset = sqrt(dmag0**2+dmag**2)
        mag_offset0, dmag_offset0 = robust_mean(mag_offset[wh],dmag_offset[wh])

        wh1 = abs(mag_offset-mag_offset0)<0.25
        if (wh1.sum()>5):
            wh *= wh1
            mag_offset, dmag_offset = robust_mean(mag_offset[wh],dmag_offset[wh])
        else:
            mag_offset,dmag_offset = mag_offset0,dmag_offset0

        mag += mag_offset
        mag_big += mag_offset
        print (""" Magnitude Offset %.6f +/- %0.6f""" % (mag_offset,dmag_offset))

        h = (phot_dat[3]<999)*(phot_dat[3]>0)*(phot_dat[2]>0)
        phot_dat[2,h] += mag_offset
        phot_dat[4,h] += mag_offset
        phot_dat[3,h] = sqrt( phot_dat[3,h]**2 + dmag_offset**2 + sys_err**2 )

        # calculate the aperture correction
        ap_corr = median(mag_big[wh]-mag[wh])
        zero_pt = sex_zero + mag_offset + 2.5*log10(gain/dt)+am - ap_corr
        print (""" Median Zero Point %.3f [gain=%.2f, dt=%.1f, am_corr=%.3f]""" % (zero_pt,gain,dt,am))
        print (""" Sextractor Aperture Correction (mag) %.3f""" % ap_corr)

        mag10s = quick_mode( mag-2.5*log10(dmag/0.10857) )
        print (" 10-sigma limiting magnitude %.2f" % mag10s)

        if (do_plot):
            x0=mag0.min()-0.3
            x1=mag0.max()+0.3
            h = dmag<999
            errorbar (mag0[h],mag0[h]-mag[h],yerr=sqrt(dmag[h]**2+dmag0[h]**2),xerr=dmag[h],fmt='bo',capsize=0,linestyle='None',markersize=5,mew=1)
            plot(mag0[wh],mag0[wh]-mag[wh],'ro',markersize=5) 
            if( (~h).sum()>0 ): plot (mag0[~h],mag0[~h]-mag[~h],'gv',ms=6,alpha=0.6)
            plot ([x0,x1],[0.,0.],':')
            xlim((x0,x1))
            xlabel("Catalog Mag");
            ylabel("Catalog Mag - Mag");
            title("""10-sigma: %.2f ; FWHM=%.1f pixels""" % (mag10s,median(fwhm)))
            ylim((-1.,1.))
            savefig(outfile+"_dm.jpg")

        of = open(outfile,'w')
        alr1 = lr1.split()
        alr1[-6] = """%f""" % (sex_zero+mag_offset)
        of.write(' '.join(alr1)+'\n')
        for i in range(len(phot_dat[0])):
            (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx) = phot_dat[:,i]
            of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx))

        of.close()

        f=open(matchfile,'w')

        f.write("#RA Dec mag dmag RA_cat Dec_cat cat_mag cat_dmag x y xa ya x2a y2a expos num\n")
        f.write("## Time Range %s - %s\n""" % (t0,t1))
        f.write("## Median Sextractor FWHM %.2f\n" % fwhm0)
        f.write("## Sextractor Mag Zero Point = %.4f +/- %.4f\n" % (sex_zero+mag_offset,dmag_offset))
        f.write("## Median Zero Point %.2f\n" % zero_pt)
        f.write("## 10-sigma limiting magnitude %.2f\n" % mag10s)
        f.write("## Magnitude Offset %.6f +/- %0.6f\n" % (mag_offset,dmag_offset))

        phot_dat = phot_dat[:,gphot]
        for j in range(len(mag0)):
            (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx) = phot_dat[:,j]
            f.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,ra,dec,mag0[j],dmag0[j],x,y,xa,ya,x2a,y2a,expos,idx0[j]))

        f.close()

    else:
        print ("No matches!")


if __name__ == "__main__":

    if (len(sys.argv)<3): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==0): usage()
    calfile=sys.argv[2]
    if (os.path.exists(calfile)==0): usage()

    outfile=infile+'.photometry.txt'
    matchfile=infile+'.match.txt'

    cal_mag_max=18.0
    if (len(sys.argv)>3): cal_mag_max=float(sys.argv[3])

    ab_offset=0.
    if (len(sys.argv)>4): ab_offset=float(sys.argv[4])

    do_plot=False
    if (len(sys.argv)>5): 
        if (sys.argv[5]=="do_plot"): do_plot=True

    calibrate(infile,calfile,outfile,matchfile,cal_mag_max=cal_mag_max,ab_offset=ab_offset,do_plot=do_plot)
