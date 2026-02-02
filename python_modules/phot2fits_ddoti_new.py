#!/usr/bin/python3
"""
 phot2fits.py timefile masterfile t0 <filter>
"""
import sys,os
import astropy.io.fits as pyfits
from numpy import loadtxt,zeros,atleast_1d,log10,dot,sqrt
from ut2gps import ut2gps
from sfdmap import SFDMap

def usage():
    print (__doc__)
    sys.exit()

# dictionary of lambda, A_lambda/E[B-V]
lam_dict = {
 'B': (435., 4.1), 'g': (481., 3.69), 'gri': (617., 2.72),
 'r': (617., 2.72), 'i': (752., 2.08), 'z': (866., 1.63),
 'zy': (914., 1.49), 'y': (962., 1.37), 'J': (1235., 0.9),
 'H': (1662., 0.51)
}

sfdir = os.environ['REDUX_BASE_DIR'] + '/catalogs/sfddata'

def phot2fits(timefile,masterfile,t0,cfilter='GAIA-EDR3 G(AB)'):
    """
       allow a maximum of nmax sources through
       if present not just in the stack, must be present at least nmin times
    """
    # get the image information
    tdata = atleast_1d(loadtxt(timefile,dtype={'names': ('file','t0','t1','dt','x0','y0','magzero','maglim','mag0','dmag0','airmass','skylev','var0','fwhm'), 'formats': ('S60','S20','S20','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4')}))

    # stack file header
    hdr = pyfits.getheader(tdata['file'][0])
    tc0, tc1 = tdata['t0'][0].decode("utf-8"), tdata['t1'][0].decode("utf-8")
    if ( ('TCSTART' in hdr) and ('TCSTOP' in hdr) ):
        if (len(hdr['TCSTART'])>0 and len(hdr['TCSTOP'])>0):
            tc0 = hdr['TCSTART']
            tc1 = hdr['TCSTOP']

    n0 = len(tdata)

    stack_base = tdata['file'][0].decode("utf-8").replace('.fits','')
    sub_cat = stack_base + '.diff_dir/calibrated_catalog.fits'

    filter = cfilter.replace("_D","")
    ref_cat = 'ps1_stack_'+filter+ '_dir/calibrated_catalog.fits'

    lam, A_lam_fac = 617., 2.72  # r-band as default
    if (filter in lam_dict): lam, A_lam_fac = lam_dict[filter]

    # check that calibrated frame data exist
    good = zeros(n0,dtype='bool')
    for i in range(n0):
        base = tdata['file'][i].decode("utf-8").replace('.fits','')
        pfile = base+'_dir/calibrated_catalog.fits'
        if (os.path.exists(pfile)): good[i]=True

    tdata = tdata[good]
    n0 = len(tdata)

    t0_gps = ut2gps(t0)[0]
    gps0,gps1 = ut2gps(tdata['t0'])-t0_gps, ut2gps(tdata['t1'])-t0_gps
    cgps0,cgps1 = ut2gps(tc0)[0]-t0_gps, ut2gps(tc1)[0]-t0_gps  # common epoch times
   
    # time epochs shared by different filters 
    common = (gps1>cgps0)*(gps0<cgps1)
    common[0] = False
    exptime_common = tdata['dt'][common].sum()

    # get the stack photometry
    mdata = pyfits.getdata(masterfile)
    id = mdata['VECTOR_ASSOC']
    flx, dflx = mdata['FLUX_APER'][:,0], mdata['FLUXERR_APER'][:,0]

    m0 = len(mdata)

    # get faint source data
    if (os.path.exists('catalog_matches.fits')):
        faintdis,faintmag,faintid,hpm = pyfits.getdata('catalog_matches.fits')
        faintid = faintid.astype('int32')
        hpm = hpm.astype('int32')
        fhdr = pyfits.getheader('catalog_matches.fits')
        have_faint=True
    else:
        print ("unable to find faint source list")
        faintdis,faintmag,faintid,hpm = zeros(m0),zeros(m0),zeros(m0,dtype='int32'),zeros(m0,dtype='int32')
        have_faint=False

    if (os.path.exists('subtraction_matches.fits')):
        subdis,submag,subid,subhpm = pyfits.getdata('subtraction_matches.fits')
        subid = subid.astype('int32')
        subhpm = subhpm.astype('int32')
        shdr = pyfits.getheader('subtraction_matches.fits')
        if (os.path.exists(sub_cat)):
            sub_dat = pyfits.getdata(sub_cat)
            flx_sub, dflx_sub = sub_dat['FLUX_APER'][:,0], sub_dat['FLUXERR_APER'][:,0]
            submag[submag>0] = 25-2.5*log10(flx_sub[submag>0].clip(1))
        if (os.path.exists(ref_cat)):
            ref_dat = pyfits.getdata(ref_cat)
            flx_ref, dflx_ref = ref_dat['FLUX_APER'][:,0], ref_dat['FLUXERR_APER'][:,0]
    else:
        flx_sub, dflx_sub, flx_ref, dflx_ref = 0.*flx,0.*flx,0.*flx,0.*flx
        print ("unable to find subtraction source list")
        subdis,submag,subid,subhpm = zeros(m0),zeros(m0),zeros(m0,dtype='int32'),zeros(m0,dtype='int32')

    # gather the data from separate epochs, including the master (stack, position 0)
    data = zeros((n0,m0,4),dtype='float32')
    files = tdata['file']
    data[0,:,0] = flx
    data[0,:,1] = dflx
    data[0,:,2] = mdata['FLUX_APER'][:,1]
    data[0,:,3] = mdata['FLUXERR_APER'][:,1]
    for i in range(n0-1):
        base = files[i+1].decode("utf-8").replace('.fits','')
        pfile = base+'_dir/calibrated_catalog.fits'
        pdata = pyfits.getdata(pfile)
        data[i+1,:,0] = pdata['FLUX_APER'][:,0]
        data[i+1,:,1] = pdata['FLUXERR_APER'][:,0]
        data[i+1,:,2] = pdata['FLUX_APER'][:,1]
        data[i+1,:,3] = pdata['FLUXERR_APER'][:,1]

    # calculate common epoch time fluxes
    flxc, dflxc = 1.*flx, 1.*dflx
    if (common.sum()<n0-1):
        print ("Calculating common epoch fluxes...")
        wt = 1/tdata['var0']
        flxc = dot(wt[common],data[common,:,0])
        dflxc = dot(wt[common]**2,data[common,:,1]**2)
        norm = dot(wt[common],data[common,:,1]>0)
        h = norm>0
        flxc[h] /= norm[h]; dflxc[h] = sqrt(dflxc[h]) / norm[h]
        flxc[~h] = 0 ; dflxc[~h] = 0

    # now write this to a fits file
    prihdr = pyfits.Header()
    prihdr['CFILTER'] = cfilter
    prihdr['TRIGTM'] = t0
    prihdr['TCSTART'] = tc0
    prihdr['GCSTART'] = cgps0
    prihdr['GCSTOP'] = cgps1
    prihdr['TCSTOP'] = tc1
    prihdr['TCEXP'] = exptime_common
    prihdr['LAMBDA'] = lam
    prihdr['ALAMFAC'] = A_lam_fac
    hdr0 = pyfits.getheader(masterfile)
    try:
        prihdr['SEXX0'] = hdr0['SEXX0']
        prihdr['SEXY0'] = hdr0['SEXY0']
        prihdr['PS'] = hdr0['PS']
    except: pass

    if (have_faint):
        ncat=fhdr["NCAT"]
        prihdr["NCAT"] = ncat
        for i in range(ncat):
            cid = """CATID%d""" % (i+1)
            prihdr[cid] = fhdr[cid]

    prihdu = pyfits.PrimaryHDU(header=prihdr)

    # first extension will be image information table
    col1 = pyfits.Column(name='filename', format='60A', array=files)
    col2 = pyfits.Column(name='DATE-OBS', format='20A', array=tdata['t0'])
    col3 = pyfits.Column(name='DATE-OBE', format='20A', array=tdata['t1'])
    col4 = pyfits.Column(name='T0', format='D', array=gps0)
    col5 = pyfits.Column(name='T1', format='D', array=gps1)
    col6 = pyfits.Column(name='exptime', format='E', array=tdata['dt'])
    col7 = pyfits.Column(name='x0', format='E', array=tdata['x0'])
    col8 = pyfits.Column(name='y0', format='E', array=tdata['y0'])
    col9 = pyfits.Column(name='fwhm', format='E', array=tdata['fwhm'])
    col10 = pyfits.Column(name='airmass', format='E', array=tdata['airmass'])
    col11 = pyfits.Column(name='skylev', format='E', array=tdata['skylev'])
    col12 = pyfits.Column(name='magzero', format='E', array=tdata['magzero'])
    col13 = pyfits.Column(name='maglim', format='E', array=tdata['maglim'])
    col14 = pyfits.Column(name='mag0', format='E', array=tdata['mag0'])
    col15 = pyfits.Column(name='dmag0', format='E', array=tdata['dmag0'])
    col16 = pyfits.Column(name='var0', format='E', array=tdata['var0'])
    cols = pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16])
    tbhdu1 = pyfits.BinTableHDU.from_columns(cols)

    # get extinction data from SF2011
    sfmap = SFDMap(sfdir)
    ebv = sfmap.ebv(mdata['ALPHA_J2000'], mdata['DELTA_J2000'])
    A_lam = A_lam_fac * ebv

    # second extension will be master stack photometry
    col1 = pyfits.Column(name='SRCID', format='J', array=mdata['VECTOR_ASSOC'])
    col2 = pyfits.Column(name='RA', format='D', array=mdata['ALPHA_J2000'])
    col3 = pyfits.Column(name='DEC', format='D', array=mdata['DELTA_J2000'])
    col4 = pyfits.Column(name='flx', format='E', array=flx)
    col5 = pyfits.Column(name='dflx', format='E', array=dflx)
    col6 = pyfits.Column(name='flx_ref', format='E', array=flx_ref)
    col7 = pyfits.Column(name='dflx_ref', format='E', array=dflx_ref)
    col8 = pyfits.Column(name='flx_sub', format='E', array=flx_sub)
    col9 = pyfits.Column(name='dflx_sub', format='E', array=dflx_sub)
    col10 = pyfits.Column(name='flxc', format='E', array=flxc)
    col11 = pyfits.Column(name='dflxc', format='E', array=dflxc)
    col12 = pyfits.Column(name='x', format='E', array=mdata['X_IMAGE'])
    col13 = pyfits.Column(name='y', format='E', array=mdata['Y_IMAGE'])
    col14 = pyfits.Column(name='fwhm', format='E', array=mdata['FWHM_IMAGE'].clip(0))
    col15 = pyfits.Column(name='HIGHPM', format='J', array=hpm)
    col16 = pyfits.Column(name='CATDIS', format='E', array=faintdis)
    col17 = pyfits.Column(name='CATMAG', format='E', array=faintmag)
    col18 = pyfits.Column(name='CATID', format='J', array=faintid)
    col19 = pyfits.Column(name='SUBMAG', format='E', array=submag)
    col20 = pyfits.Column(name='A_lam', format='E', array=A_lam)
    col21 = pyfits.Column(name='FLAGS', format='J', array=mdata['FLAGS'])
    cols = pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21])
    tbhdu2 = pyfits.BinTableHDU.from_columns(cols)

    # third extension is the data image
    imghdu = pyfits.ImageHDU(data=data)

    hdulist = pyfits.HDUList([prihdu,tbhdu1,tbhdu2,imghdu])
    hdulist.writeto('photometry.fits',overwrite=True)


if __name__ == "__main__":
    """
    make a photometry fits file
    """
    if (len(sys.argv)<4): usage()

    timefile=sys.argv[1]
    if (os.path.exists(timefile)==False): usage()
    masterfile=sys.argv[2]
    if (os.path.exists(masterfile)==False): usage()
    t0=sys.argv[3]

    cfilter='USNO-R'
    if (len(sys.argv)>4): cfilter=sys.argv[4]

    phot2fits(timefile,masterfile,t0=t0,cfilter=cfilter)
