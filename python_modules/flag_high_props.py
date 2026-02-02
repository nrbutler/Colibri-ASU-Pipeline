#!/usr/bin/python3
"""
flag_high_props.py radecfile.txt high_prop_file.fits mp_radec.txt

   flag high proper motion stars
"""
import sys,os
from numpy import loadtxt,sqrt,abs,pi,sin,cos,where,zeros,ones,vstack,log10,median
from astropy.io.fits import getdata,writeto,getheader

def usage():
    print (__doc__)
    sys.exit()


def do_the_flagging(r,d,r0,d0,dr0,dd0,m,m0):
    """
    """
    dis0=sqrt(dd0**2+dr0**2)*3600.

    nr=len(r)
    good = ones(nr,dtype='bool')
    dis_ar = zeros(nr,dtype='float32')
    ii_ar = zeros(nr,dtype='int32')

    for i in range(nr):
        cd=cos(d[i]*pi/180.) 
        dis = sqrt( ((r[i]-r0)*cd)**2+(d[i]-d0)**2 )*3600.
        dis1= abs( dd0*(r[i]-r0)*cd-dr0*(d[i]-d0) )/sqrt( dd0**2+dr0**2 )*3600.
        h = (dis<10*dis0)*(dis1<dis0)*(abs(m[i]-m0)<0.5)
        if (h.sum()>0):
            dis = dis[h]
            j = dis.argmin()
            ii_ar[i] = j
            dis_ar[i] = dis[j]
            good[i] = False

    return good,dis_ar[~good],ii_ar[~good]


def flag_high_props(infile='radecfile.txt',propfile='high_prop_file.fits',radecfile='null.txt',ofile='catalog_matches.fits'):
    """
    """
    phot_dat=getdata(infile)
    good0 = phot_dat['VECTOR_ASSOC']<0

    if (good0.sum()>0):

        n0 = len(phot_dat)
        wg = where(good0)[0]
        phot_dat = phot_dat[good0]

        cat_id=0
        if (os.path.exists(ofile)):
            os.system("""sethead BSCALE=1 %s""" % ofile)
            res = getdata(ofile).astype('int16')
            ohdr = getheader(ofile)
            cat_id += ohdr["NCAT"]
        else: res = zeros((4,n0),dtype='int16')

        cat_id0 = cat_id
        cat_ids=[]
        cat_names=[]

        r=phot_dat['ALPHA_J2000']
        d=phot_dat['DELTA_J2000']
        m=25-2.5*log10(phot_dat['FLUX_APER'][:,0].clip(1))
        fwhm=phot_dat['FWHM_IMAGE']
        fwhm0 = median(fwhm[fwhm>1])
        fwhm = 2*fwhm.clip(fwhm0)  #   convert to arcsec

        if (os.path.exists(propfile)):
            r0,d0,m0,dm0,dr0,dd0 = getdata(propfile)

            good,dis,ii = do_the_flagging(r,d,r0,d0,dr0,dd0,m,m0)
            if (good.sum()>0):
                cat_id+=1
                cat_ids.append(cat_id)
                cat_names.append("usno")
                whg=wg[~good]
                res[0,whg] = (dis*10).round().astype('int16')
                res[1,whg] = (m0[ii]*10).round().astype('int16')
                res[2,whg] = 10*cat_id
                res[3,whg] = 10

                nr=len(r)
                nr1 = good.sum()
                sys.stderr.write("""Flagging %d detections near high-proper-motion stars...\n""" % (nr-nr1))

        if (os.path.exists(radecfile)):
            mpdata = loadtxt(radecfile,dtype={'names': ('ra','dec','vmag','name','dra','ddec','err'), 'formats': ('f4','f4','f4','S20','f4','f4','f4')},ndmin=1)
            nmp,nmp0=0,len(mpdata)
            if (nmp0>0):
                fac = 180./pi
                sa,ca = sin(r/fac),cos(r/fac)
                sd,cd = sin(d/fac),cos(d/fac)
                for i in range(nmp0):
                    ra0,dec0,vmag,name,dra0,ddec0,err0=mpdata[i]

                    sa0,ca0 = sin(ra0/fac),cos(ra0/fac)
                    sd0,cd0 = sin(dec0/fac),cos(dec0/fac)
                    cos_c = sd0*sd + cd0*cd*(ca*ca0+sa*sa0)
                    # mean, offset in arcsec
                    dx = fac*cd*(sa*ca0-ca*sa0) *3600./ cos_c
                    dy = fac*(cd0*sd-sd0*cd*(ca*ca0+sa*sa0)) *3600./ cos_c

                    # dis2 = (dx + s*dra0)**2 + (dy + s*ddec0)**2
                    #  where s goes from -0.5 to 0.5 over the observation, so min is:
                    s = - ( (dx*dra0+dy*ddec0)/(dra0**2+ddec0**2) ).clip(-0.5,0.5)
                    dis2_min = (dx + s*dra0)**2 + (dy + s*ddec0)**2

                    jj = where(dis2_min <= 100. + fwhm**2 + err0**2)[0]
                    if (len(jj)>0):
                        nmp+=len(jj)
                        whg=wg[jj]
                        res[0,whg] = (sqrt(dis2_min[jj])*10).round().astype('int16')
                        res[1,whg] = (mpdata['vmag'][i]*10).round().astype('int16')
                        res[3,whg] = -10
                        cat_id += 1
                        res[2,whg] = 10*cat_id
                        cat_ids.append(cat_id)
                        cat_names.append(mpdata['name'][i].decode('UTF-8').replace('\'',''))

                if (nmp>0): sys.stderr.write("""Flagging %d detections near minor planets...\n""" % nmp)

        if (cat_id>cat_id0):        
            if (cat_id0>0): writeto(ofile, res,ohdr,overwrite=True)
            else: writeto(ofile, res,overwrite=True)

        if (cat_id>0):
            os.system("""sethead BSCALE=0.1 NCAT=%.d %s""" % (cat_id,ofile))
            for i in range(len(cat_ids)):
                os.system("""sethead CATID%d=%s %s""" % (cat_ids[i],cat_names[i],ofile))


if __name__ == "__main__":

    if (len(sys.argv)<3): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==0): usage()
    propfile=sys.argv[2]

    radecfile='null.txt'
    if (len(sys.argv)>3): radecfile=sys.argv[3]

    flag_high_props(infile,propfile,radecfile=radecfile)
