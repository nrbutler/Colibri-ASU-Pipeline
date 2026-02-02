#!/usr/bin/python3
"""
 ddoti_lc_plots.py photfile <snr_cut>
"""
import sys,os
from astropy.io.fits import getheader,getdata
from numpy import loadtxt,sqrt,median,log10,where,ones,zeros,savetxt,vstack,arange,array
from scipy.signal import medfilt
from numpy.ma import median as masked_median
from numpy.ma import masked_array

from fit_wcs import ra2sex,dec2sex

from multiprocessing import Pool,Lock

import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import subplots,ioff,close
ioff()

from linfit import linfit

def usage():
    print (__doc__)
    sys.exit()


def init(lock0):
    """
    to avoid drawing more than one lc plot at once
    """
    global lock
    lock = lock0


def make_lc_plot(lc_file):
    """
    """
    lcf = open(lc_file,'r')
    ln0 = lcf.readline()
    lcf.close()

    pdat = ln0.split()[10:]
    tt,dtt,m,dm,phot = loadtxt(lc_file,unpack=True,usecols=(2,3,5,6,7))
    bad_phot = phot<1

    j2 = dm>=0.5; j = ~j2

    fig,ax = subplots()
    ax.set_ylim(((m+2*dm).max(),(m-2*dm).min()))
    ax.set_xlim(((tt-dtt/2.).min(),(tt+dtt/2.).max()))
    if (bad_phot.sum()>0): ax.plot (tt[bad_phot],m[bad_phot],'ro',alpha=0.5,label='non-photometric')
    ax.errorbar(tt,m,xerr=dtt/2,yerr=dm,marker='o',capsize=0,linestyle='None',markersize=3,mew=1)
    ax.plot (tt[j2],m[j2]+0.5,'bv')
    ax.plot (tt[j2],m[j2]-0.5,'bo')
    ii = tt[j].argsort()
    ax.plot (tt[j][ii],float(pdat[6])-float(pdat[7])*2.5*log10(tt[j][ii]),label="""%s +/- %s""" % (pdat[7],pdat[8]))
    ax.legend()

    ax.set_xlabel("""Time Since %s - %s [hours]""" % (pdat[4],pdat[5]),fontsize=14)
    ax.set_ylabel(pdat[0],fontsize=14)
    ax.set_title("""Light Curve for Source %s (RA=%s, Dec=%s)""" % (pdat[1],pdat[2],pdat[3]))

    lock.acquire()
    fig.savefig("""lc_%s.jpg""" % pdat[1])
    lock.release()

    close(fig)


def lc_plots_ddoti(photfile,outfile,cfilter='USNO-R',snr_cut=10.,snr_min=2.,nmax=1000,do_lc_plots=True,do_variability_plot=False):
    """
       do some fitting and plotting
    """
    hdr=getheader(photfile)
    # read in the summary of images reduced (ignoring stack)
    imdata=getdata(photfile,1)
    # read in the summary of sources detected
    srcdata=getdata(photfile,2)

    srt = srcdata['SRCID'].argsort()[::-1]
    srcdata = srcdata[srt]
    id = srcdata['SRCID']

    # start by comparing first and second half fluxes
    # read in the lightcurve data again
    data = getdata(photfile,3)[:,srt]
    flx0, dflx0 = 1.*data[0,:,0], 1.*data[0,:,1]
    m0 = 25-2.5*log10(flx0.clip(1))
    dm0 = 1.0857*dflx0/flx0.clip(1)
    flx = 1.*data[1:,:,0]
    dflx = 1.*data[1:,:,1]
    wt = dflx>0

    m,n=flx.shape; n2=m//2
    f0 = zeros(n,dtype='float64')
    df0 = zeros(n,dtype='float64')
    f1 = zeros(n,dtype='float64')
    df1 = zeros(n,dtype='float64')

    norm = wt[:n2].sum(axis=0)
    h = norm>0
    f0[h] = ( (flx[:n2,h]*wt[:n2,h]).sum(axis=0) / norm[h] ).clip(1.)
    df0[h] = sqrt( ( dflx[:n2,h]**2*wt[:n2,h] ).sum(axis=0) )/ norm[h]
    norm = wt[n2:].sum(axis=0)
    h = norm>0
    f1[h] = ( (flx[n2:,h]*wt[n2:,h]).sum(axis=0) / norm[h] ).clip(1.)
    df1[h] = sqrt( ( dflx[n2:,h]**2*wt[n2:,h] ).sum(axis=0) )/ norm[h]

    good = (df0>0)*(df1>0)
    wgood = where(good)[0]
    delta_mag = -2.5*log10(f1[good]/f0[good])
    delta_mag_err = 1.0857*sqrt( (df0[good]/f0[good])**2 + (df1[good]/f1[good])**2 )

    darker = delta_mag>(2*delta_mag_err).clip(0.1)
    brighter = delta_mag<-(2*delta_mag_err).clip(0.1)
    fading = zeros(n,dtype='int16')
    fading[wgood[darker]] = 1
    fading[wgood[brighter]] = -1
    if (do_variability_plot):
        snr_fade = zeros(n,dtype='float32')
        snr_fade[wgood] = delta_mag/delta_mag_err

    ng0 = len(id)

    hpm = srcdata['HIGHPM']
    catdis = srcdata['CATDIS']
    catmag = srcdata['CATMAG']
    catid = srcdata['CATID']
    submag = srcdata['SUBMAG']

    ncat=hdr.get('NCAT') or 0
    cnames = ["nocat"]
    for i in range(ncat):
       cnames.append( hdr["""CATID%d""" % (i+1)] )

    cnames = array(cnames,dtype='S20')
    catname = cnames[catid].astype('str')

    #
    # determine the new and calibration source lists        
    #
    good = (id!=0)*(srcdata['FLAGS']!=1)*(flx0>snr_min*dflx0)
    if (len(dflx)>0): good *= ((dflx>0).sum(axis=0)>1)
    deblended = (srcdata['FLAGS']==2).astype('int16')
    cal = (id>0)*good
    new = (id<0)*good

    new_faint = new*(dm0>=1.0857/snr_cut)
    new_bright = new*(dm0<1.0857/snr_cut)

    # faint, new wrt ps1,USNO,2mass,etc?
    cat_detected = catmag>0
    catmag[catmag<0] *= -1
    new_faint *= ~cat_detected
    new_bright *= ~cat_detected
    sub_detected = ( (submag>0)*(~cat_detected) ).astype('int16')

    x,y = srcdata['x'],srcdata['y']
    try:
        sex_x0,sex_y0,ps = hdr['SEXX0'],hdr['SEXY0'],hdr['PS']
        # in pixels, multiply by arcsec/pixel (ps)
        swift_dis = sqrt( (x-sex_x0)**2+(y-sex_y0)**2 )*ps
    except:
        swift_dis = zeros(n)

    ra,dec,fwhm = srcdata['RA'],srcdata['DEC'],srcdata['fwhm']
    id1 = id.astype('str')

    ras,decs=[],[]
    if (len(ra)>0):
        ras = ra2sex(ra)
        decs = dec2sex(dec)

    hn = hpm<0
    if (hn.sum()>0):
        #mhdr=' id   ra           dec    x,y    mag          dmag         fwhm high_pm fade subdet deblend grb_offset catmag catoffset catname'
        mhdr=' id   RA           DEC                                 mag          dmag         fwhm high_pm fade subdet deblend grb_offset catmag catoffset catname'
        mdat=vstack((id1[hn],ras[hn],decs[hn],ra[hn],dec[hn],x[hn],y[hn],m0[hn],dm0[hn],fwhm[hn],hpm[hn],fading[hn],submag[hn],deblended[hn],swift_dis[hn],catmag[hn],catdis[hn],catname[hn])).T
        savetxt('mp_phot.txt',mdat,fmt='%s',header=mhdr)

    # new (faint or bright, only in cat)
    cat_only = new*cat_detected
    new *= ~cat_detected

    # these sources can be uncatalogued, cat detected, usno detected (id>0)
    bhdr=' x y mag dmag id fwhm uncatalogued lcgenerated catmag catoffset catname'
    uncatalogued=zeros(len(m0),dtype='int16')
    uncatalogued[new] = 1
    lcgenerated=zeros(len(m0),dtype='int16')
    if (len(imdata)>1): lcgenerated[new_bright]=1
    savetxt('all_sources.txt',vstack((x[good],y[good],m0[good],dm0[good],id1[good],fwhm[good],uncatalogued[good],lcgenerated[good],catmag[good],catdis[good],catname[good])).T,fmt='%s',header=bhdr)

    # can make a plot about fading
    #  have new_faint, new_bright, cal, and cat_only using snr_fade versus m0
    if (do_variability_plot):
        fig,ax = subplots()
        ax.plot(m0[cal],snr_fade[cal],'o',color='k',alpha=0.01,label='Cal Source')
        ax.plot(m0[cat_only],snr_fade[cat_only],'o',alpha=0.2,label='PS1 Source')
        ax.plot(m0[new_faint],snr_fade[new_faint],'o',alpha=0.5,label='New Faint')
        ax.plot(m0[new_bright],snr_fade[new_bright],'o',alpha=0.5,label='New Bright')
        ax.set_xlabel("Source Magnitude",fontsize=14)
        ax.set_ylabel("Magnitude Change (S/N)",fontsize=14)
        ax.set_ylim(( -0.5+min(snr_fade[new_bright].min(),snr_fade[new_faint].min()), 0.5+max(snr_fade[new_bright].max(),snr_fade[new_faint].max()) ))
        ax.set_xlim((max(m0[new_bright].max(),m0[new_faint].max()),13.))
        ax.legend(loc='right')
        fig.savefig('variability.jpg')
        close(fig)

    #new_hdr=' id   ra           dec    x,y    mag          dmag         fwhm high_pm fade subdet deblend grb_offset catmag catoffset catname'
    new_hdr=' id   RA           DEC                                 mag          dmag         fwhm high_pm fade subdet deblend grb_offset catmag catoffset catname'
    if (new_faint.sum()>0):
        new_dat=vstack((id1[new_faint],ras[new_faint],decs[new_faint],ra[new_faint],dec[new_faint],x[new_faint],y[new_faint],m0[new_faint],dm0[new_faint],fwhm[new_faint],hpm[new_faint],fading[new_faint],submag[new_faint],deblended[new_faint],swift_dis[new_faint],catmag[new_faint],catdis[new_faint],catname[new_faint])).T
        savetxt('new_phot.txt',new_dat,fmt='%s',header=new_hdr)
    else: sys.stderr.write("No faint sources\n")
    if (cat_only.sum()>0):
        cat_dat=vstack((id1[cat_only],ras[cat_only],decs[cat_only],ra[cat_only],dec[cat_only],x[cat_only],y[cat_only],m0[cat_only],dm0[cat_only],fwhm[cat_only],hpm[cat_only],fading[cat_only],submag[cat_only],deblended[cat_only],swift_dis[cat_only],catmag[cat_only],catdis[cat_only],catname[cat_only])).T
        savetxt('cat_only_phot.txt',cat_dat,fmt='%s',header=new_hdr)
    else: sys.stderr.write("No cat sources\n")
    if (cal.sum()>0):
        cal_dat=vstack((id1[cal],ras[cal],decs[cal],ra[cal],dec[cal],x[cal],y[cal],m0[cal],dm0[cal],fwhm[cal],hpm[cal],fading[cal],submag[cal],deblended[cal],swift_dis[cal],catmag[cal],catdis[cal],catname[cal])).T
        savetxt('cal_phot.txt',cal_dat,fmt='%s',header=new_hdr)
    else: sys.stderr.write("No calibration sources\n")

    t0=hdr['TRIGTM']
    cfilter=hdr['CFILTER']

    # read in the summary of images reduced (ignoring stack)
    t10,t20,dt0 = imdata['DATE-OBS'][0], imdata['DATE-OBE'][0],imdata['exptime'][0]
    fwhm_stack = imdata['fwhm'][0]
    imdata=imdata[1:]
    if (len(imdata)==0):

        j0 = new_bright
        idx_ar,flx0_ar,dflx0_ar,fader = srcdata['SRCID'][j0],flx0[j0],dflx0[j0],fading[j0]
        x_ar,y_ar = srcdata['x'][j0],srcdata['y'][j0]
        ra0_ar,dec0_ar,fwhm_ar,highpm_ar = srcdata['RA'][j0],srcdata['DEC'][j0],srcdata['FWHM'][j0],srcdata['HIGHPM'][j0]
        ra0_ars,dec0_ars = ras[j0],decs[j0]
        swift_dis_ar = swift_dis[j0]
        deblended_ar = deblended[j0]
        sub_detected_ar = sub_detected[j0]
        submag_ar = submag[j0]
        catmag_ar = catmag[j0]
        catdis_ar = catdis[j0]
        catname_ar = catname[j0]

        mag0_ar = 25-2.5*log10(flx0_ar)
        dmag0_ar = 1.0857*dflx0_ar/flx0_ar

        ofile=open(outfile,'w')
        ofile.write("#   id    RA           DEC   x,y     mag      dmag     fwhm     slope    dslope       chi2 n_detect high_pm fade subdet deblend grb_offset catmag catoffset catname\n")
        for i in range(j0.sum()):
            ofile.write("""%6d %11s %11s %9.5f %9.5f %8.2f %8.2f %8.4f %8.4f %8.4f %8.4f %8.4f %10.2f %4d/%-4d %5d %5d %5.2f %5d %8.1f %8.1f %8.1f %s\n""" % (idx_ar[i],ra0_ars[i],dec0_ars[i],ra0_ar[i],dec0_ar[i],x_ar[i],y_ar[i],mag0_ar[i],dmag0_ar[i],fwhm_ar[i],0,0,0,0,0,highpm_ar[i],0,submag_ar[i],deblended_ar[i],swift_dis_ar[i],catmag_ar[i],catdis_ar[i],catname_ar[i]))

        ofile.close()
        return 0

        mag0_ar = 25-2.5*log10(flx0_ar)
        dmag0_ar = 1.0857*dflx0_ar/flx0_ar

        #number of sources
        M=len(idx_ar)
        if (M>nmax): M=nmax

        # read in the lightcurve data

    photometric = (imdata['maglim']>18+1.25*log10(imdata['exptime']/1.e3)-2.5*log10(imdata['fwhm']/3.5))*(imdata['magzero']>21.25)
    print (f"Fraction of photometric images: {photometric.mean():.2f}")

    times,times1,expos = imdata['DATE-OBS'], imdata['DATE-OBE'],imdata['exptime']/3600.
    gps1,gps2 = imdata['T0']/3600.,imdata['T1']/3600.
    gps = 0.5*(gps1+gps2)

    gps0=0.
    if (gps.min()<0):
        gps0 -= gps.min()
        gps1 += gps0
        gps2 += gps0

    # only consider new sources for plotting
    j0 = new_bright

    lc_files=[]
    if (j0.sum()>0):

        idx_ar,flx0_ar,dflx0_ar,fader = srcdata['SRCID'][j0],flx0[j0],dflx0[j0],fading[j0]
        x_ar,y_ar = srcdata['x'][j0],srcdata['y'][j0]
        ra0_ar,dec0_ar,fwhm_ar,highpm_ar = srcdata['RA'][j0],srcdata['DEC'][j0],srcdata['FWHM'][j0],srcdata['HIGHPM'][j0]
        ra0_ars,dec0_ars = ras[j0],decs[j0]
        swift_dis_ar = swift_dis[j0]
        deblended_ar = deblended[j0]
        sub_detected_ar = sub_detected[j0]
        submag_ar = submag[j0]
        catmag_ar = catmag[j0]
        catdis_ar = catdis[j0]
        catname_ar = catname[j0]

        mag0_ar = 25-2.5*log10(flx0_ar)
        dmag0_ar = 1.0857*dflx0_ar/flx0_ar

        #number of sources
        M=len(idx_ar)
        if (M>nmax): M=nmax

        # read in the lightcurve data
        data = getdata(photfile,3)[:,srt[j0],:]

        #
        # first plot any lightcurves
        #

        ofile=open(outfile,'w')
        ofile.write("#   id    RA           DEC   x,y     mag      dmag     fwhm     slope    dslope       chi2 n_detect high_pm fade subdet deblend grb_offset catmag catoffset catname\n")

        for i in range(M):
            j = data[1:,i,1]>0
            n1=j.sum()

            id0,mag0,dmag0,ra0s,dec0s,ra0,dec0,x,y,fwhm0,highpm,sdis,c1mag,c1dis,c1name=idx_ar[i],mag0_ar[i],dmag0_ar[i],ra0_ars[i],dec0_ars[i],ra0_ar[i],dec0_ar[i],x_ar[i],y_ar[i],fwhm_ar[i],highpm_ar[i],swift_dis_ar[i],catmag_ar[i],catdis_ar[i],catname_ar[i]

            dt,t1,t2,t1g,t2g,phot = 1.*expos[j],times[j].copy(),times1[j].copy(),1.*gps1[j],1.*gps2[j],1.*photometric[j]
            n = where(j)[0]

            flx = 1.*data[1:,i,0][j]
            dflx = 1.*data[1:,i,1][j]
            wt = 1/imdata['var0']

            # divide the n1 data points into groups of ngroup
            snr_group = 5.
            if (dmag0>0.10857): snr_group=3.
            ngroup = max(1,int( n1*(snr_group*dmag0/1.0857)**2 ))
            #  will have nbins bins:
            nbins = n1//ngroup
            if (ngroup>1 and nbins>1):

                ia = arange(n1//ngroup)*ngroup
                ib = arange(1,n1//ngroup+1)*ngroup; ib[-1] = n1
                flx_a=zeros(nbins,dtype=flx.dtype)
                dflx_a=zeros(nbins,dtype=dflx.dtype)
                dt_a=zeros(nbins,dtype=dt.dtype)
                phot_a=zeros(nbins,dtype=dt.dtype)
                t1g_a=zeros(nbins,dtype=t1g.dtype)
                t2g_a=zeros(nbins,dtype=t2g.dtype)
                t1_a=zeros(nbins,dtype=t1.dtype)
                t2_a=zeros(nbins,dtype=t2.dtype)
                for k in range(nbins):
                    w0 = wt[ia[k]:ib[k]].sum()
                    flx_a[k] = (flx[ia[k]:ib[k]]*wt[ia[k]:ib[k]]).sum()/w0
                    dflx_a[k] = sqrt( ((dflx[ia[k]:ib[k]]*wt[ia[k]:ib[k]])**2).sum() )/w0
                    dt_a[k] = dt[ia[k]:ib[k]].sum()
                    phot_a[k] = phot[ia[k]:ib[k]].mean()
                    t1g_a[k],t2g_a[k] = t1g[ia[k]],t2g[ib[k]-1]
                    t1_a[k],t2_a[k] = t1[ia[k]],t2[ib[k]-1]

                flx,dflx = flx_a,dflx_a
                t1,t2,t1g,t2g,dt,phot = t1_a,t2_a,t1g_a,t2g_a,dt_a,phot_a

            j=flx>2.1715*dflx; j2=~j
            nj = j.sum()

            tt,dtt = 0.5*(t1g+t2g),t2g-t1g
            m = 0*flx; dm = 0*flx
            m[j] = 25-2.5*log10(flx[j])
            dm[j] = 1.0857*dflx[j]/flx[j]
            dm[j2] = 0.5; m[j2] = 25-2.5*log10(flx[j2].clip(0)+3*dflx[j2]) + 0.5

            if (nj>=2):
                # only fit sources with 2 or more >2 sigma detections
                xx=-2.5*log10(tt)
                res = linfit(xx[j],m[j],dy=dm[j],slope_prior_err=1.)
                sig = (m[j]-res[0]-res[1]*xx[j]).std()
                if (sig<0.01): sig=0.01
                var = dm[j]**2 + sig**2

                xm = xx[j].mean(); ym = m[j].mean()
                #res = linfit(xx[j]-xm,m[j]-ym,dy=sqrt(var),slope_prior_err=1.)
                res = linfit(xx[j]-xm,m[j]-ym,dy=sqrt(var))
                slp=res[1]
                dslp=sqrt(res[2][1][1])

                dt1 = dt[j].max()
                dt2 = dt[j].sum()-dt1
                if(dt1>10*dt2): slp=0.

                # handle over-fitting
                if (j.sum()<=2):
                    chi2=0.
                    if (slp<0): slp=min(slp+dslp,0.)
                    else: slp=max(slp-dslp,0.)
                else:
                    chi2 = ( ((m[j]-ym-res[0]-res[1]*(xx[j]-xm))/dm[j])**2 ).sum()/(j.sum()-2)
                #   slp /= 1.+dslp**2*chi2
                #   dslp /= 1.+dslp**2*chi2

                slp0 = res[0]-slp*xm+ym
                if (slp > -dslp and fader[i]==1): fader[i] = 0
                if (slp < dslp and fader[i]==-1): fader[i] = 0
            else:
                chi2,slp,slp0,dslp=999.,999.,999.,999.

            of = """lc_%d.txt""" % id0
            lc_files.append(of)
            of_hdr="""time1 time2 t_av dt exposure mag dmag photometric : %s %d %11s %11s %s %.2f %.4f %.4f %.4f""" % (cfilter.replace(' ','_'),id0,ra0s,dec0s,t0,gps0/3600.,slp0,slp,dslp)
            yy = vstack((t1,t2,tt,dtt,dt*3600.,m,dm,phot)).T
            savetxt(of,yy,header=of_hdr,fmt='%s')

            ofile.write("""%6d %11s %11s %9.5f %9.5f %8.2f %8.2f %8.4f %8.4f %8.4f %8.4f %8.4f %10.2f %4d/%-4d %5d %5d %5.2f %5d %8.1f %8.1f %8.1f %s\n""" % (id0,ra0s,dec0s,ra0,dec0,x,y,mag0,dmag0,fwhm0,slp,dslp,chi2,nj,len(j),highpm,fader[i],submag_ar[i],deblended_ar[i],sdis,c1mag,c1dis,c1name))

        ofile.close()

    # make light curve plots
    if (do_lc_plots and len(lc_files)>0):
        lock0 = Lock()
        pool = Pool(initializer=init, initargs=(lock0,))
        pool.map( make_lc_plot, lc_files )

    #
    # now make a summary plot
    #
    f = srcdata['fwhm'].clip(fwhm_stack/4,fwhm_stack*10)
    # read in the lightcurve data again
    data = getdata(photfile,3)[1:,srt,:]

    mag_zero = imdata['magzero']
    gm = mag_zero >= median(mag_zero)-1.
    flx, dflx = data[gm,:,0], data[gm,:,1]

    j = (dflx>0)*(flx>1.)
    jj = j.sum(axis=0)
    h = jj>=jj.max()/2
    jh = j[:,h]

    m0=m0[h]; f=f[h];
    m=25-2.5*log10(flx[:,h].clip(1))
    dm = 1.0857*dflx[:,h]/flx[:,h].clip(1)

    ma = masked_array(m,~jh)
    dma = masked_array(dm,~jh)
    ma -= m0
    r0 = masked_median(dma,axis=0).data
    ma /= dma/r0

    ng = jh.sum(axis=0)
    dm00 = ma.std(axis=0)*sqrt( ng/(ng-1) )
    s = m0.argsort()
    r0 = medfilt(r0[s],11)

    fig,ax = subplots()
    axs = ax.scatter(m0,dm00,c=f,alpha=0.2)
    ax.set_xlim((m0.min()-0.2,m0.max()+0.2)); ax.set_ylim((dm00.min()/2,dm00.max()*2))
    ax.set_xlabel(cfilter+" Mag",fontsize=14)
    ax.set_ylabel("Photometric Scatter [mag]",fontsize=14)
    ax.plot (m0[s],r0,label='Expectation',alpha=0.8)
    ax.plot (m0[s],sqrt(1.e-4+r0**2),label='Expectation (+1% sys)',alpha=0.8)
    ax.semilogy()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb=fig.colorbar(axs,cax=cax,orientation='vertical')
    cb.set_label('fwhm')
    ax.legend()
    fig.savefig("scatter_plot.jpg")


if __name__ == "__main__":
    """
    make summary photometry plots
    """
    if (len(sys.argv)<2): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==False): usage()

    outfile='source_fitting.txt'

    snr_cut=10.
    if (len(sys.argv)>2): snr_cut=float(sys.argv[2])

    lc_plots_ddoti(infile,outfile,snr_cut=snr_cut)
