#!/usr/bin/python3
"""
 backsub.py image_file wmap_file <bg_order>
"""
import os, sys
from numpy import abs,median,linspace,dot,empty
from scipy.linalg import lstsq

from astropy.io.fits import getdata,writeto

def usage():
    print (__doc__)
    sys.exit()


def backsub(image,wt=[],order=4,return_bg=False):
    """
      weighted polynomial background subtraction
    """
    if (order<0): order=0
    nf = ((1 + order)*(2+order))//2

    nx,ny=image.shape
    x,y = linspace(-1,1,nx), linspace(-1,1,ny)

    fx0 = empty((1+order,nx),dtype='float64')
    fy0 = empty((1+order,ny),dtype='float64')
    fx0[0],fy0[0] = 1,1
    if (order>0): fx0[1],fy0[1] = x,y

    # legendre polynomials
    for i in range(1,order):
        fx0[i+1] = ((2*i+1)*x*fx0[i] - i*fx0[i-1])/(i+1.)
        fy0[i+1] = ((2*i+1)*y*fy0[i] - i*fy0[i-1])/(i+1.)

    fx = empty((nf,nx),dtype='float64')
    fy = empty((nf,ny),dtype='float64')

    i0 = 0
    for i in range(order+1):
        fx[i0:i0+i+1],fy[i0:i0+i+1] = fx0[:i+1],fy0[i::-1]
        i0 += i+1

    if (len(wt)==0):
        yy = (dot(fx,image)*fy).sum(axis=1)
        vec = yy/( (fx**2).sum(axis=1)*(fy**2).sum(axis=1) )
    else:

        matr = empty((nf,nf),dtype='float64')
        if (order>0): matr[0] = (dot(fx,wt)*fy).sum(axis=1)
        else: matr[0,0] = wt.astype('float64').sum()
        for k in range(1,nf):
            matr[1:k+1,k] = (dot(fx[k]*fx[1:k+1],wt)*fy[k]*fy[1:k+1]).sum(axis=1)

        for k in range(nf): matr[k,:k] = matr[:k,k]

        yy = (dot(fx,image*wt)*fy).sum(axis=1)

        try:
            vec = lstsq(matr,yy)[0]
        except:
            print ("Matrix inversion failed")
            vec = 0.*yy

    if (return_bg): return dot(fx.T*vec,fy)
    else: image -= dot(fx.T*vec,fy)


def clean_image(img,weight):
    """
       subtract median columns and rows, using masking
    """
    nx,ny = img.shape
    imgT = img.T
    weightT = weight.T
    for i in range(ny):
        h = weightT[i]>0
        if (h.sum()>0): imgT[i] -= median(imgT[i,h])
    for i in range(nx):
        h = weight[i]>0
        if (h.sum()>0): img[i] -= median(img[i,h])
    #img[weight==0]=0


def qclean_image(img,weight):
    """
       subtract median quadrants, using masking
    """
    nx,ny = img.shape
    nx2,ny2 = nx//2,ny//2
    h = weight[:nx2,:ny2]>0
    img[:nx2,:ny2] -= median( img[:nx2,:ny2][h] )
    h = weight[nx2:,:ny2]>0
    img[nx2:,:ny2] -= median( img[nx2:,:ny2][h] )
    h = weight[:nx2,ny2:]>0
    img[:nx2,ny2:] -= median( img[:nx2,ny2:][h] )
    h = weight[nx2:,ny2:]>0
    img[nx2:,ny2:] -= median( img[nx2:,ny2:][h] )


if __name__ == '__main__':
    """
    """



if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<3): usage()

    datafile= sys.argv[1]
    if (os.path.exists(datafile)==False): usage()

    wmapfile=sys.argv[2]
    if (os.path.exists(wmapfile)==False): usage()

    bg_order=4
    if (len(sys.argv)>3): bg_order=int(sys.argv[3])

    image,hdr = getdata(datafile,header=True)
    wt = getdata(wmapfile,ignore_missing_end=True)

    if (bg_order<-1): qclean_image(image,wt)
    elif (bg_order==-1): clean_image(image,wt)
    else: backsub(image,wt,order=bg_order)

    writeto(datafile,image,hdr,overwrite=True)
