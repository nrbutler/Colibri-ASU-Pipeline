#!/usr/bin/python3
"""
do_match1.py radec.txt catalog_radec.txt <cal_mag_max> <ab_offset> <do_plot>
   for comparing a sextractor radec list to catalogs
"""
import sys,os
from coord import match
from quick_mode import quick_mode
from robust_mean import robust_mean
from numpy import loadtxt,sqrt,where,log10,arange,isnan,array,median,cos,pi,ones,zeros,savetxt

def usage():
    print (__doc__)
    sys.exit()


def do_match(infile,calfile,outfile,matchfile,outfile_matched,outfile_notmatched,aperture=10.0,sys_err=0.0001,mag_min=14.,mag_max=18.,ab_offset=0.,do_plot=False):
    """
    """
    if (do_plot):
        import matplotlib as mpl
        mpl.use('Agg')
        from matplotlib.pyplot import plot,xlabel,ylabel,savefig,semilogy,title,legend,clf,scatter,semilogx,xlim,ylim,colorbar,hist,errorbar

    (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx)=loadtxt(infile,unpack=True)

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

    ap_corr=999

    # insert x,y dependent systematic error
    h = dmag>0
    dmag[h] = sqrt( dmag[h]**2 + sys_err**2 )
    phot_err = sqrt( dmag**2 + 0.01*(xa**2+ya**2) )

    (ra0,dec0,mag0,dmag0)=loadtxt(calfile,unpack=True,usecols=(0,1,2,3))
    j = dmag0<999
    ra0=ra0[j]; dec0=dec0[j]; mag0=mag0[j]; dmag0=dmag0[j]
    try:
        idx0=loadtxt(calfile,usecols=(14,))[j]
    except:
        idx0=arange(len(ra0))

    mdec0 = median(dec0); cmdec0=cos(mdec0*pi/180); mra0 = median(ra0)
    dra = (ra-mra0)*cmdec0; dra0 = (ra0-mra0)*cmdec0;
    ddec = dec-mdec0; ddec0 = dec0-mdec0
    ii=match(dra0,ddec0,dra,ddec,3.0)

    w=where( (ii['sep']>=0) )[0]

    matched=zeros(len(mag),dtype='bool')

    nn = len(w)
    if (nn>0):
        w1 = ii['ind'][w]
        kk=w1.argsort()

        mag0 = mag0[w]; dmag0 = dmag0[w];
        ra0 = ra0[w]; dec0=dec0[w]
        idx0 = idx0[w]

        matched[w1] = True

        wh = (mag0>=mag_min)*(mag0<=mag_max)*(fwhm[w1]<median(fwhm[w1])*2.)*(mag[w1]>0)
        if (wh.sum()<5): wh = mag0<=mag_max
        if (wh.sum()<5): wh = ones(len(mag0),dtype='bool')
        wa,w1a = where(wh)[0],w1[wh]

        # calculate the aperture correction
        if (ap_corr>990): ap_corr,dap_corr = robust_mean(mag_big[w1a]-mag[w1a],phot_err[w1a])

        h = dmag[w1]<0.01
        if (h.sum()>10): fwhm0 = quick_mode(fwhm[w1[h]])
        else: fwhm0,dfwhm0 = robust_mean(fwhm[w1],phot_err[w1])

        print (" Matches %d" % nn)
        print (""" Time Range %s - %s""" % (t0,t1))
        print (""" Median Sextractor FWHM %.2f""" % fwhm0)

        mag0 += ab_offset
        mag_offset = mag0[wa] - mag[w1a]
        dmag_offset = sqrt(dmag0[wa]**2+dmag[w1a]**2)

        mag_offset0, dmag_offset0 = robust_mean(mag_offset,dmag_offset)
        wh1 = abs(mag_offset-mag_offset0)<0.25
        if (wh1.sum()>5):
            mag_offset, dmag_offset = robust_mean(mag_offset[wh1],dmag_offset[wh1])
            wa,w1a = wa[wh1],w1a[wh1]
        else:
            mag_offset,dmag_offset = mag_offset0,dmag_offset0

        if (isnan(dmag_offset)): dmag_offset=0.

        print (""" Magnitude Offset %.6f +/- %0.6f""" % (mag_offset,dmag_offset))
        #dmag[dmag<=0] = 1.e-6

        # last term is the aperture correction from sextractor
        zero_pt = sex_zero + mag_offset + 2.5*log10(gain/dt)+am - ap_corr

        print (""" Median Zero Point %.3f [gain=%.2f, dt=%.1f, am_corr=%.3f]""" % (zero_pt,gain,dt,am))
        print (""" Sextractor Aperture Correction (mag) %.3f""" % ap_corr)

        h = mag>0
        mag[h] += mag_offset
        mag_big[h] += mag_offset
        dmag[h] = sqrt( dmag[h]**2+dmag_offset**2 )

        hm = h * (dmag<0.2)
        if (hm.sum()<10): hm = h
        mag10s = quick_mode( mag[hm]-2.5*log10(dmag[hm]/0.10857) )
        print (" 10-sigma limiting magnitude %.2f" % mag10s)

        if (do_plot):
            x0=mag0.min()-0.3
            x1=mag0.max()+0.3
            errorbar (mag0,mag0-mag[w1],yerr=sqrt(dmag[w1]**2+dmag0**2),xerr=dmag[w1],fmt='bo',capsize=0,linestyle='None',markersize=5,mew=1,alpha=0.2)
            plot(mag0[wa],mag0[wa]-mag[w1a],'ro',markersize=5)
            plot ([x0,x1],[0.,0.],':')
            #xlim((x0,x1))
            xlim((14,21.5))
            xlabel("Catalog Mag");
            ylabel("Catalog Mag - Mag");
            title("""10-sigma: %.2f ; FWHM=%.1f pixels""" % (mag10s,median(fwhm)))
            ylim((-1,1))
            #ylim((-0.2,0.2))
            savefig(outfile+"_dm.jpg")

        of = open(outfile,'w')
        #of1 = open(outfile_matched,'w')
        #of2 = open(outfile_notmatched,'w')
        alr1 = lr1.split()
        alr1[-6] = """%f""" % (sex_zero+mag_offset)
        of.write(' '.join(alr1)+'\n')
        for i in range(len(mag)):
            of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra[i],dec[i],mag[i],dmag[i],mag_big[i],dmag_big[i],fwhm[i],x[i],y[i],xa[i],ya[i],x2a[i],y2a[i],expos[i],idx[i]))
            #if (matched[i]): of1.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra[i],dec[i],mag[i],dmag[i],mag_big[i],dmag_big[i],fwhm[i],x[i],y[i],xa[i],ya[i],x2a[i],y2a[i],expos[i],idx[i]))
            #else: of2.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra[i],dec[i],mag[i],dmag[i],mag_big[i],dmag_big[i],fwhm[i],x[i],y[i],xa[i],ya[i],x2a[i],y2a[i],expos[i],idx[i]))

        of.close()
        #of1.close()
        #of2.close()
        savetxt(outfile_matched,matched,fmt='%d',header=' index')

        f=open(matchfile,'w')

        f.write("#RA Dec mag dmag RA_cat Dec_cat cat_mag cat_dmag x y xa ya x2a y2a expos num\n")
        f.write("## Time Range %s - %s\n""" % (t0,t1))
        f.write("## Median Sextractor FWHM %.2f\n" % fwhm0)
        f.write("## Sextractor Mag Zero Point = %.4f +/- %.4f\n" % (sex_zero+mag_offset,dmag_offset))
        f.write("## Median Zero Point %.2f\n" % zero_pt)
        f.write("## 10-sigma limiting magnitude %.2f\n" % mag10s)
        f.write("## Magnitude Offset %.6f +/- %0.6f\n" % (mag_offset,dmag_offset))

        for i in range(len(mag0)):
            j=w1[i]
            f.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra[j],dec[j],mag[j],dmag[j],ra0[i],dec0[i],mag0[i],dmag0[i],x[j],y[j],xa[j],ya[j],x2a[j],y2a[j],expos[j],idx0[i]))

        f.close()

    else:
        print ("No matches!")


def main():

    if (len(sys.argv)<3): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==0): usage()
    calfile=sys.argv[2]
    if (os.path.exists(calfile)==0): usage()

    outfile=infile+'.photometry.txt'
    outfile_matched=infile+'.matched.txt'
    outfile_notmatched=infile+'.notmatched.txt'
    matchfile=infile+'.match.txt'

    cal_mag_max=18.0 
    if (len(sys.argv)>3): cal_mag_max=float(sys.argv[3])

    ab_offset=0.
    if (len(sys.argv)>4): ab_offset=float(sys.argv[4])

    do_plot=False
    if (len(sys.argv)>5):
        if (sys.argv[5]=="do_plot"): do_plot=True

    do_match(infile,calfile,outfile,matchfile,outfile_matched,outfile_notmatched,mag_max=cal_mag_max,ab_offset=ab_offset,do_plot=do_plot)

if __name__ == "__main__":
    main()

