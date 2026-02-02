import os
from numpy import loadtxt,maximum,vstack,arange,zeros,where,linspace
from scipy.interpolate import interp1d

fdir = os.environ['REDUX_BASE_DIR'] + '/calfiles/filters/'

def bin_filter(lam,f,nbin=10):
    wh = where(f>0)[0]
    i0,i1 = wh[0],wh[-1]
    ll = linspace( lam[i0-1],lam[i1+1], nbin+1 )
    cf = interp1d(lam[i0-1:i1+2],f[i0-1:i1+2].cumsum())(ll)
    cf = cf[1:]-cf[:-1]
    return (ll[1:]+ll[:-1])/2., cf/cf.sum()

def get_filters(lam,nbin=10):
    """
    lam in nm, set nbin=0 to skip binning
    """
    lB,B0 = loadtxt(fdir+'B.dat',unpack=True)
    lg,g0 = loadtxt(fdir+'g.dat',unpack=True)
    lr,r0 = loadtxt(fdir+'r.dat',unpack=True)
    li,i0 = loadtxt(fdir+'i.dat',unpack=True)
    lz,z0 = loadtxt(fdir+'z.dat',unpack=True)
    ly,y0 = loadtxt(fdir+'y.dat',unpack=True)
    lJ,J0 = loadtxt(fdir+'J.dat',unpack=True)
    lH,H0 = loadtxt(fdir+'H.dat',unpack=True)

    B = interp1d(lB,B0,bounds_error=False,fill_value=0)(lam*10).clip(0)
    g = interp1d(lg,g0,bounds_error=False,fill_value=0)(lam*10).clip(0)
    r = interp1d(lr,r0,bounds_error=False,fill_value=0)(lam*10).clip(0)
    i = interp1d(li,i0,bounds_error=False,fill_value=0)(lam*10).clip(0)
    z = interp1d(lz,z0,bounds_error=False,fill_value=0)(lam*10).clip(0)
    y = interp1d(ly,y0,bounds_error=False,fill_value=0)(lam*10).clip(0)
    J = interp1d(lJ,J0,bounds_error=False,fill_value=0)(lam*10).clip(0)
    H = interp1d(lH,H0,bounds_error=False,fill_value=0)(lam*10).clip(0)

    gri = maximum(maximum(g,r),i)
    zy = maximum(z,y)

    PB = B/lam; PB /= PB.sum()
    Pg = g/lam; Pg /= Pg.sum()
    Pr = r/lam; Pr /= Pr.sum()
    Pi = i/lam; Pi /= Pi.sum()
    Pz = z/lam; Pz /= Pz.sum()
    Py = y/lam; Py /= Py.sum()
    PJ = J/lam; PJ /= PJ.sum()
    PH = H/lam; PH /= PH.sum()
    Pgri = gri/lam; Pgri /= Pgri.sum()
    Pzy = zy/lam; Pzy /= Pzy.sum()

    matr = vstack((PB, Pg, Pr, Pi, Pz, Py, PJ, PH, Pgri, Pzy))

    if (nbin>0):
        Nf=len(matr)
        lam_matr0 = zeros((Nf,nbin),dtype='float32')
        filter_matr0 = zeros((Nf,nbin),dtype='float32')
        for i in range(Nf):
            l1,f1 = bin_filter(lam,matr[i],nbin=nbin)
            lam_matr0[i] = l1
            filter_matr0[i] = f1

        return lam_matr0, filter_matr0

    return lam,matr
