#!/usr/bin/python3
"""
fit_wcs.py fits_file matchfile outfile
"""

def usage():
    print (__doc__)
    sys.exit()

from numpy import sqrt,sin,cos,sqrt,arctan,arctan2,arcsin,arccos,zeros,pi,loadtxt,where,sign,array,log,atleast_1d,empty,modf
import sys,os
from scipy.optimize import fmin
from astropy.io.fits import getdata,getheader,writeto

def ra2sex(ra):
    """
    """
    hfrac,hd = modf(ra/15.+1.e-8)
    min_frac,m = modf(hfrac*60)

    if (len(atleast_1d(ra))>1):
        n=len(ra)
        ras=empty(n,dtype='|S11')
        for i in range(n):
            ras[i] = """%02.0f:%02.0f:%05.2f""" % (hd[i],m[i],min_frac[i]*60)
        ras = ras.astype('str')
    else:
        ras = """%02.0f:%02.0f:%05.2f""" % (hd,m,min_frac*60)

    return ras

def dec2sex(dec):
    """
    """
    hfrac,hd = modf(dec+sign(dec)*1.e-8)
    min_frac,m = modf(abs(hfrac)*60)

    if (len(atleast_1d(dec))>1):
        n=len(dec)
        decs=empty(n,dtype='|S11')
        for i in range(n):
            decs[i] = """%+03.0f:%02.0f:%04.1f""" % (hd[i],m[i],min_frac[i]*60)
        decs = decs.astype('str')
    else:
        decs = """%+03.0f:%02.0f:%04.1f""" % (hd,m,min_frac*60)

    return decs

C0_xy0=array([513,513])
C1_xy0=array([513,513])
C2_xy0=array([1177.,1031])
C3_xy0=array([924,982])
C4_xy0=array([283,189])
par_C0 = array([30.0,0.0,0.3168,0.3171,0.68,0.,0.])
par_C1 = array([30.0,0.0,0.3171,0.3191,1.48,0.,0.])
par_C2 = array([30.0,0.0,0.2988,0.2954,-89.0,-19.60381671,-4128.15179797])
par_C3 = array([30.0,0.0,0.2983,-0.2945,89.5,-19.60381671,-4128.15179797])
par_C4 = array([30.0,0.0,0.29010212,-0.29010212,91.69,0.,0.])

hdr_C0={}
hdr_C0['CRPIX1'] = C0_xy0[0]; hdr_C0['CRPIX2'] = C0_xy0[1]
hdr_C1={}
hdr_C1['CRPIX1'] = C1_xy0[0]; hdr_C1['CRPIX2'] = C1_xy0[1]
hdr_C2={}
hdr_C2['CRPIX1'] = C2_xy0[0]; hdr_C2['CRPIX2'] = C2_xy0[1]
hdr_C3={}
hdr_C3['CRPIX1'] = C3_xy0[0]; hdr_C3['CRPIX2'] = C3_xy0[1]
hdr_C4={}
hdr_C4['CRPIX1'] = C4_xy0[0]; hdr_C4['CRPIX2'] = C4_xy0[1]

def movhdr(x,y,hdr):
    """
    move wcs reference pixel to x,y
    """
    ra,dec = xy2ad(x,y,hdr)
    hdr['CRPIX1'], hdr['CRPIX2'] = x,y
    hdr['CRVAL1'], hdr['CRVAL2'] = ra,dec


def hdr2par(hdr,distort=False):
    """
    """
    d2r=pi/180.

    if (distort):
        par=zeros(7,dtype='float32')
    else:
        par=zeros(5,dtype='float32')

    par[0:2] = hdr['CRVAL1'], hdr['CRVAL2']
    cd1_1, cd1_2 = hdr['CD1_1'], hdr['CD1_2']
    cd2_1, cd2_2 = hdr['CD2_1'], hdr['CD2_2']

    par[4] = arctan2(cd2_1,-cd1_1)
    if (abs(cos(par[4]))>0):
        s1 = -sign(cd1_1/cos(par[4]))
        s2 = sign(cd2_2/cos(par[4]))
    else:
        s1 = sign(cd2_1/sin(par[4]))
        s2 = sign(cd1_2/sin(par[4]))

    par[4] /= d2r

    par[2] = s1*sqrt(cd1_1**2+cd2_1**2)*3600.
    par[3] = s2*sqrt(cd2_2**2+cd1_2**2)*3600.

    if (distort):
        try:
            par[5] = hdr['PV1_17']
            par[6] = hdr['PV1_31']
        except:
            par[5] = hdr['PV1_7']
            par[6] = hdr['PV1_17']

    return par


def par2hdr(par,hdr,distort=False):
    """
    need to fix
    """
    d2r=pi/180.

    hdr['CRVAL1']=float(par[0])
    hdr['CRVAL2']=float(par[1])
    hdr['CD1_1'] = float(-par[2]*cos(par[4]*d2r)/3600.)
    hdr['CD2_1'] = float(par[2]*sin(par[4]*d2r)/3600.)
    hdr['CD1_2'] = float(par[3]*sin(par[4]*d2r)/3600.)
    hdr['CD2_2'] = float(par[3]*cos(par[4]*d2r)/3600.)

    if (distort):
        a=float(par[5])
        b=float(par[6])
        hdr['PV1_1']=hdr['PV2_1']=1.
        try:
            x = hdr['PV1_7']
            hdr['PV1_7']=hdr['PV1_9']=hdr['PV2_7']=hdr['PV2_9']=a
            hdr['PV1_17']=hdr['PV1_21']=hdr['PV2_17']=hdr['PV2_21']=b
            hdr['PV1_19']=hdr['PV2_19']=2*b
        except:
            hdr['PV1_17']=hdr['PV2_17']=hdr['PV1_21']=hdr['PV2_21']=a
            hdr['PV1_19']=hdr['PV2_19']=2*a
            hdr['PV1_31']=hdr['PV2_31']=hdr['PV1_37']=hdr['PV2_37']=b
            hdr['PV1_33']=hdr['PV2_33']=hdr['PV1_35']=hdr['PV2_35']=3*b


def xy2ij(x,y,hdr,distort=False):

    try:
        cd1_1, cd1_2 = hdr['CD1_1'], hdr['CD1_2']
        cd2_1, cd2_2 = hdr['CD2_1'], hdr['CD2_2']
    except: 
        cd1_1, cd1_2 = hdr['CDELT1'],0.
        cd2_2, cd2_1 = hdr['CDELT2'],0.
    x0, y0 = hdr['CRPIX1'], hdr['CRPIX2']

    i = cd1_1*(x-x0)+cd1_2*(y-y0)
    j = cd2_1*(x-x0)+cd2_2*(y-y0)

    if (distort):

        order=5
        try:
            try: a,b = hdr['PV1_17'],hdr['PV1_31']
            except:
                a,b = hdr['PV1_7'],hdr['PV1_17'] 
                order=3
        except:
            cam='X'
            if (cam=='C0'):
                a,b=par_C0[5], par_C0[6]
            elif (cam=='C1'):
                a,b=par_C1[5], par_C1[6]
            elif (cam=='C2'):
                a,b=par_C2[5], par_C2[6]
            elif (cam=='C4'):
                a,b=par_C4[5], par_C4[6]
            elif (cam=='C3'):
                a,b=par_C3[5], par_C3[6]
            else:
                a,b=0.,0.

        r2 = i*i+j*j
        if (order==5): fac=(a*r2**2 + b*r2**3)
        else: fac=(a*r2 + b*r2**2)
        i += fac*i
        j += fac*j

    return i,j


def ij2ad(i,j,hdr,sex=False):

    fac = 57.29577951

    xi,eta = i/fac,j/fac
    xi0, eta0 = hdr['CRVAL1']/fac, hdr['CRVAL2']/fac
    seta0, ceta0 = sin(eta0), cos(eta0)

    cosc = 1./sqrt(1.+xi**2.+eta**2.)
    a = fac*( xi0 + arctan2( xi, ceta0 - eta*seta0 ) )
    d = fac*( arcsin(cosc*(seta0 + eta*ceta0)) )

    try:
        if (a<0): a+=360
        elif (a>360): a-=360
    except:
     a[a<0]+=360
     a[a>360]-=360

    if (sex): a,d = ra2sex(a),dec2sex(d)

    return a,d


def xy2ad(x,y,hdr,distort=False,sex=False):
    i,j = xy2ij(x,y,hdr,distort=distort)
    return ij2ad(i,j,hdr,sex=sex)


def ad2ij(a,d,hdr):
    fac = 57.29577951

    a0, d0 = hdr['CRVAL1']/fac, hdr['CRVAL2']/fac

    cos_c = sin(d0)*sin(d/fac) + cos(d0)*cos(d/fac)*cos(a/fac-a0)
    i = fac*cos(d/fac)*sin(a/fac-a0) / cos_c
    j = fac*(cos(d0)*sin(d/fac)-sin(d0)*cos(d/fac)*cos(a/fac-a0)) / cos_c

    return i,j


def ij2xy(i,j,hdr,distort=False):

    if (distort):
        p11,p12,p14,p15,p16=hdr['PV1_1'], hdr['PV1_2'], hdr['PV1_4'], hdr['PV1_5'], hdr['PV1_6']
        p21,p22,p24,p25,p26=hdr['PV2_1'], hdr['PV2_2'], hdr['PV2_4'], hdr['PV2_5'], hdr['PV2_6']

        i1,j1 = 1.*i,1.*j
        for ii in range(3):
            di = i1*(p11-1) + j1*p12 + i1*i1*p14 + i1*j1*p15 + j1*j1*p16
            dj = j1*(p21-1) + i1*p22 + j1*j1*p24 + i1*j1*p25 + i1*i1*p26
            i1,j1 = i - di,j - dj
    else: 
        i1, j1 = i,j

    cd1_1, cd1_2 = hdr['CD1_1'], hdr['CD1_2']
    cd2_1, cd2_2 = hdr['CD2_1'], hdr['CD2_2']
    denom = cd1_1*cd2_2 - cd1_2*cd2_1

    x0, y0 = hdr['CRPIX1'], hdr['CRPIX2']

    x = x0 + (i1*cd2_2 - j1*cd1_2)/denom
    y = y0 + (j1*cd1_1 - i1*cd2_1)/denom

    return x,y


def ad2xy(a,d,hdr,distort=False):
    i,j = ad2ij(a,d,hdr)
    return ij2xy(i,j,hdr,distort=distort)


def xy2ad_par(x,y,par,x0=1025,y0=1025):

    fac = 57.29577951

    secpix1,secpix2,theta=par[2:]
    cd1_1 = secpix1*cos(theta)
    cd2_1 = secpix1*sin(theta)
    cd1_2 = -secpix2*sin(theta)
    cd2_2 = secpix2*cos(theta)

    xi = ( cd1_1*(x-x0)+cd1_2*(y-y0) )/ fac
    eta = ( cd2_1*(x-x0)+cd2_2*(y-y0) )/ fac

    xi0, eta0 = par[0:2]
    seta0, ceta0 = sin(eta0), cos(eta0)

    cosc = 1./sqrt(1.+xi**2.+eta**2.)
    a = xi0 + arctan2( xi, ceta0 - eta*seta0 )
    d = arcsin(cosc*(seta0 + eta*ceta0))

    return a*fac,d*fac


def comp_radec(x,y,ra0,dec0,par,hdr,distort=False):
    #model ~ i0-i
    par2hdr(par,hdr,distort=distort)
    i,j=xy2ij(x,y,hdr,distort=distort)
    #i,j=xy2ij(x,y,hdr,distort=True)
    i0,j0=ad2ij(ra0,dec0,hdr)
    return ( (i-i0)**2+(j-j0)**2 ).sum()


par2hdr(par_C0,hdr_C0,distort=True)
hdr_C0['CCD_NAME']='C0'
hdr_C0['NAXIS1']=1024; hdr_C0['NAXIS2']=1024
par2hdr(par_C1,hdr_C1,distort=True)
hdr_C1['CCD_NAME']='C1'
hdr_C1['NAXIS1']=1024; hdr_C1['NAXIS2']=1024
par2hdr(par_C2,hdr_C2,distort=True)
hdr_C2['CCD_NAME']='C2'
hdr_C2['NAXIS1']=2040; hdr_C2['NAXIS2']=2040
par2hdr(par_C3,hdr_C3,distort=True)
hdr_C3['CCD_NAME']='C3'
hdr_C3['NAXIS1']=2040; hdr_C3['NAXIS2']=2040
par2hdr(par_C4,hdr_C4,distort=True)
hdr_C4['CCD_NAME']='C4'
hdr_C4['NAXIS1']=546; hdr_C4['NAXIS2']=368

def main():
    """
    """
    if (len(sys.argv)<4): usage()

    fits_file=sys.argv[1];
    matchfile=sys.argv[2];
    outfile=sys.argv[3];
    cam=sys.argv[4];

    hdr=getheader(fits_file)
    if (cam=='C0'): movhdr(hdr_C0['CRPIX1'],hdr_C0['CRPIX2'],hdr)
    if (cam=='C1'): movhdr(hdr_C1['CRPIX1'],hdr_C1['CRPIX2'],hdr)
    if (cam=='C2'): movhdr(hdr_C2['CRPIX1'],hdr_C2['CRPIX2'],hdr)
    if (cam=='C3'): movhdr(hdr_C3['CRPIX1'],hdr_C3['CRPIX2'],hdr)
    if (cam=='C4'): movhdr(hdr_C4['CRPIX1'],hdr_C4['CRPIX2'],hdr)
    #par2hdr(par_C2,hdr,distort=True)
    par0=hdr2par(hdr)

    (ra0,dec0,m,ra,dec,m,dm,x,y,num)=loadtxt(matchfile,unpack=True)
    if (cam=='C0' or cam=='C1'): rad = sqrt( (x-513)**2+(y-513)**2 )
    else: rad = sqrt( (x-1025)**2+(y-1025)**2 )
    j=where(rad<800)

    def get_chi2(par):
        xx=comp_radec(x[j],y[j],ra0[j],dec0[j],par,hdr)
        #return xx/4.68021257769e-05 + abs(log(abs(par[2]/par[3]))/1.e-3)
        return xx/4.68021257769e-05

    res = fmin(get_chi2,par0)
    res = fmin(get_chi2,res)
    par2hdr(res,hdr)

    dat=getdata(fits_file)
    os.remove(fits_file)
    writeto(fits_file,dat,hdr)

    #i, j = xy2ij(x,y,hdr,distort=True)
    i, j = xy2ij(x,y,hdr)
    i0, j0 = ad2ij(ra0,dec0,hdr)

    f = open(outfile,'w')
    for k in range(len(x)):
        f.write("""%f %f %f %f\n""" % (i[k],j[k],i0[k]-i[k],j0[k]-j[k]))

    f.close()

    print ("SEC_PIX1=%.6f SECPIX2=%.6f (%.4f) Theta=%.4f" % (res[2],res[3],res[2]/res[3],res[4]*180/pi))


if __name__ == "__main__":
        main()
