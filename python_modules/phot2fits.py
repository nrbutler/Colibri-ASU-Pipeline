#!/usr/bin/python3
"""
 phot2fits.py timefile masterfile t0 <filter>
"""
import sys,os
from astropy.io.fits import Header,PrimaryHDU,Column,ColDefs,BinTableHDU,ImageHDU,HDUList
from numpy import loadtxt,zeros,atleast_1d,arange
from ut2gps import ut2gps

def usage():
    print (__doc__)
    sys.exit()

def phot2fits(timefile,masterfile,t0,cfilter='USNO-R(AB)'):
    """
       allow a maximum of nmax sources through
       if present not just in the stack, must be present at least nmin times
    """
    # get the image information
    tdata = atleast_1d(loadtxt(timefile,dtype={'names': ('file','imid', 'imnum', 'impid', 't0','t1','dt','fwhm','magzero','maglim','mag0','dmag0'), 'formats': ('S60','i4','i4','i4','S20','S20','f4','f4','f8','f8','f8','f8')}))
    n0 = len(tdata)

    t0_gps = ut2gps(t0)
    gps0,gps1 = ut2gps(tdata['t0'])-t0_gps, ut2gps(tdata['t1'])-t0_gps

    # get the photometry master information
    mdata = atleast_1d(loadtxt(masterfile,dtype={'names': ('id','ra', 'dec', 'mag', 'dmag','fwhm','dt'), 'formats': ('i','f8','f8','f8','f8','f4','f4')}))
    m0 = len(mdata)

    i0,i1 = mdata['id'].min(), mdata['id'].max()
    ii = zeros(i1-i0+1,dtype='int32')
    ii[mdata['id']-i0] = arange(1,m0+1,dtype='int32')

    # gather the data from separate epochs, ignoring the master (stack)
    data = zeros((n0-1,m0+1,6),dtype='float32')
    files = tdata['file']
    for i in range(1,n0):
        base = files[i].decode("utf-8").replace('.fits','')
        pfile = base+'_dir/'+base+'_radec.txt.photometry.txt'
        try:
            dat = atleast_1d(loadtxt(pfile,unpack=True,usecols=(0,1,2,3,6,13,14)))
            id1 = dat[-1].astype('int32')
            good = (id1>=i0)*(id1<=i1)
            ids = ii[id1[good]-i0]
            data[i-1,ids] = dat[:-1,good].T
        except:
            sys.stderr.write("""Failed to read file %s\n""" % pfile)

    data = data[:,1:]

    # now write this to a fits file
    prihdr = Header()
    prihdr['CFILTER'] = cfilter
    prihdr['TRIGTM'] = t0
    prihdu = PrimaryHDU(header=prihdr)

    # first extension will be image information table
    col1 = Column(name='filename', format='60A', array=files)
    col2 = Column(name='IMID', format='I', array=tdata['imid'])
    col3 = Column(name='IMNUM', format='I', array=tdata['imnum'])
    col4 = Column(name='IMPID', format='I', array=tdata['impid'])
    col5 = Column(name='DATE-OBS', format='20A', array=tdata['t0'])
    col6 = Column(name='DATE-OBE', format='20A', array=tdata['t1'])
    col7 = Column(name='T0', format='D', array=gps0)
    col8 = Column(name='T1', format='D', array=gps1)
    col9 = Column(name='dt', format='E', array=tdata['dt'])
    col10 = Column(name='fwhm', format='E', array=tdata['fwhm'])
    col11 = Column(name='magzero', format='E', array=tdata['magzero'])
    col12 = Column(name='maglim', format='E', array=tdata['maglim'])
    col13 = Column(name='mag0', format='E', array=tdata['mag0'])
    col14 = Column(name='dmag0', format='E', array=tdata['dmag0'])
    cols = ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14])
    tbhdu1 = BinTableHDU.from_columns(cols)

    # second extension will be master stack photometry
    col1 = Column(name='SRCID', format='J', array=mdata['id'])
    col2 = Column(name='RA', format='D', array=mdata['ra'])
    col3 = Column(name='DEC', format='D', array=mdata['dec'])
    col4 = Column(name='mag', format='E', array=mdata['mag'])
    col5 = Column(name='dmag', format='E', array=mdata['dmag'])
    col6 = Column(name='fwhm', format='E', array=mdata['fwhm'])
    col7 = Column(name='dt', format='E', array=mdata['dt'])
    cols = ColDefs([col1,col2,col3,col4,col5,col6,col7])
    tbhdu2 = BinTableHDU.from_columns(cols)

    # third extension is the data image
    imghdu = ImageHDU(data=data)

    hdulist = HDUList([prihdu,tbhdu1,tbhdu2,imghdu])
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

    cfilter='USNO-R(AB)'
    if (len(sys.argv)>4): cfilter=sys.argv[4]

    phot2fits(timefile,masterfile,t0=t0,cfilter=cfilter)
