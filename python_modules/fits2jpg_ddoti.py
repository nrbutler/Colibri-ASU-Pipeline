#!/usr/bin/python3
"""
 fits2jpg.py fitsfile [size] [noinvert] [make_thumbs] [thumb_size] [thumb_edge] [save_fits]
"""

import sys,os

#from PIL import Image
from PIL import Image, ImageDraw, ImageFont
from multiprocessing import Pool

from astropy.io.fits import getdata,writeto
from numpy.random import rand
from numpy import sqrt,log,zeros,arcsinh,pi,arctan2,cos,array,ceil,floor
from fits_resize import fits_resize

from fit_wcs import xy2ad

save_fits=False

def usage():
    print (__doc__)
    sys.exit()

font = ImageFont.truetype(font="arialbd.ttf", size=24)
tlabel="PS1"

def run_imsave(indata):
    #Image.fromarray(obj=indata[1][::-1], mode='F').convert('L').save(indata[0],format = 'jpeg')
    img = Image.fromarray(obj=indata[1][::-1], mode='F').convert('RGB')
    draw = ImageDraw.Draw(img)
    draw.text((10,10), tlabel, font=font, fill="green")
    img.save(indata[0],format = 'jpeg')

def fits2jpg(file0,out_size=1024,wt_min=0,invert=True,logscale=False,zscale=True,log_fac=5.,make_thumbs=False,ts=301,te=15,ij0=[],ttag="",tag_label=""):
    """
    """
    outfile=file0.replace('.fz','').replace('.fits','.jpg')
    wfile=file0.replace('.fits','.wt.fits')
    sys.stderr.write("""Creating %s from %s ...\n""" % (outfile,file0))

    global tlabel
    tlabel = tag_label

    have_wfile=False
    if (os.path.exists(wfile)):
        sys.stderr.write("""Using weight file %s\n""" % (wfile))
        have_wfile=True
        wt = getdata(wfile)

    x,hdr=getdata(file0,header=True)
    a,b=x.shape
    n=a*b

    if (make_thumbs): fits_resize(x,wt,hdr,file0,fac=15)

    try: sky_levl = hdr['SKYLEV'] + hdr['READ_NS']**2/hdr['GAIN']
    except: sky_levl=1.

    if (have_wfile):
        h = wt>wt_min
        x[~h]=0
        x[h]*=sqrt(wt[h])

    # summary jpg images
    if (zscale):
        if (n>1.e4): nsample = int(1.e-3*n)
        else: nsample=n
        xind = (rand(nsample)*(a-1)+1).astype('int32')
        yind = (rand(nsample)*(b-1)+1).astype('int32')

        x1 = x[xind,yind]
        s= x1.argsort()
        s1,s2 = s[int(0.25*nsample)],s[int(0.95*nsample)]
        x1,x2 = x1[s1],x1[s2]
    else:
        x1,x2 = x.min(),x.max()

    if (x2>x1):
        if (logscale):
            y = ( 255*(arcsinh(x*log_fac/x2)*x2/log_fac-x1)/(x2-x1) ).clip(0,255)
        elif (zscale): y = (255*(x-x1)/(x2-x1)).clip(0,255)
        else: y = (255*(x-x1)/(x2-x1)).clip(0,255)
    else: y=1.*x

    # make the black black
    blk=0.
    if (have_wfile):
        blk = y[wt<=wt_min].mean()
        y[wt<=wt_min]=0

    if (invert): y = 255.-y

    if (make_thumbs):

        #pixel is about 2"
        # thumb size (ts) 301 pixels is about 10'
        nt_x = int(ceil(a/ts))
        nt_y = int(ceil(b/ts))
        extra_x = a - nt_x*ts
        extra_y = b - nt_y*ts

        if (len(ij0)==0):
            i0,j0 = int(extra_x/2), int(extra_y/2)
        else:
            i0,j0 = ij0
        ydat=[]

        ia, ib = int(floor((te-i0)/ts)), int(floor((a-i0-te)/ts))
        ja, jb = int(floor((te-j0)/ts)), int(floor((b-j0-te)/ts))
        for i in range(ia,ib+1):
            padx1=padx2=0
            if (i==ia): padx1 = te-i0-i*ts
            if (i==ib): padx2 = te+i0+(i+1)*ts-a
            for j in range(ja,jb+1):
                pady1=pady2=0
                if (j==ja): pady1 = te-j0-j*ts
                if (j==jb): pady2 = te+j0+(j+1)*ts-b
                r0,d0 = xy2ad(1+j0+j*ts+ts//2,1+i0+i*ts+ts//2,hdr)
                r1,d1 = xy2ad(2+j0+j*ts+ts//2,1+i0+i*ts+ts//2,hdr)
                theta = arctan2(d1-d0,(r0-r1)*cos(d0*pi/180.))*180/pi
                r,d = xy2ad(1+j0+j*ts+ts//2,1+i0+i*ts+ts//2,hdr,sex=True)
                tfile = """%sthumb_%d_%d_%s_%s_%.0f_.jpg""" % (ttag,j0+j*ts,i0+i*ts,r,d,theta)
                yt = zeros((ts+2*te,ts+2*te),dtype=y.dtype)
                if (invert): yt += 255
                #yt[padx1:ts+2*te-padx2,pady1:ts+2*te-pady2] = 1*y[i0+i*ts-te+padx1:i0+(i+1)*ts+te-padx2,j0+j*ts-te+pady1:j0+(j+1)*ts+te-pady2][::-1]
                yt[padx1:ts+2*te-padx2,pady1:ts+2*te-pady2] = 1*y[i0+i*ts-te+padx1:i0+(i+1)*ts+te-padx2,j0+j*ts-te+pady1:j0+(j+1)*ts+te-pady2]
                ydat.append( (tfile,yt) )

        pool = Pool()
        pool.map( run_imsave, ydat )

    #
    # make the stack jpg
    #
    if (out_size<=0): out_size=a
    if (a!=out_size):
        nx,ny = out_size,int(round(1.*out_size*a/b))
        sys.stderr.write("""Resizing Image %s to %dx%d\n""" % (outfile,nx,ny))
        y = array(Image.fromarray(obj=y, mode='F').resize(size=(nx,ny), resample=Image.BICUBIC))

    if (invert): y = (255*y/(255-blk)).clip(0,255)
    else: y = (255*(y-blk)/(255-blk)).clip(0,255)
    Image.fromarray(obj=y[::-1], mode='F').convert('L').save(outfile,format = 'jpeg')


if __name__ == "__main__":
    """
    """
    if (len(sys.argv)<2): usage()

    file0=sys.argv[1]
    if (os.path.exists(file0)==False): usage()

    sz=1024
    if (len(sys.argv)>2): sz=int(sys.argv[2])

    invert=True
    logscale=False
    zscale=True
    if (len(sys.argv)>3):
        if(sys.argv[3]=="noinvert"): invert=False
        elif(sys.argv[3]=="linvert"): invert,logscale = True,True
        elif(sys.argv[3]=="noz"): zscale=False
        elif(sys.argv[3]=="lnoz"): logscale,zscale=True,False
        elif(sys.argv[3]=="lnoinvert"): invert,logscale = False,True

    make_thumbs=False
    if (len(sys.argv)>4):
        if (sys.argv[4]=='make_thumbs'): make_thumbs=True

    thumb_size=401
    if (len(sys.argv)>5): thumb_size=int(sys.argv[5])

    thumb_edge=30
    if (len(sys.argv)>6): thumb_edge=int(sys.argv[6])

    ij0 = []
    if (len(sys.argv)>7): ij0=eval(sys.argv[7])

    save_fits=False
    if (len(sys.argv)>8): 
        if (sys.argv[8]=="save_fits"): save_fits=True
 
    tag_label=""
    if (len(sys.argv)>9): tag_label=sys.argv[9]

    ttag=""
    if (len(sys.argv)>10): ttag=sys.argv[10]
    
    fits2jpg(file0,out_size=sz,invert=invert,logscale=logscale,zscale=zscale,make_thumbs=make_thumbs,ts=thumb_size,te=thumb_edge,ij0=ij0,ttag=ttag,tag_label=tag_label)
