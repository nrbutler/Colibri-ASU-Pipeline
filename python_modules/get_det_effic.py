#!/usr/bin/python3
"""
 get_det_effic.py
"""

from astropy.io.fits import getdata,getheader
from numpy import arange,log,exp,where,log10,zeros,array,maximum
from glob import glob
from linfit import linfit
from scipy.interpolate import interp1d

import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import subplots,ioff,close
ioff()

# dictionary of lambda, A_lambda/E[B-V]
lam_dict = {
 'B': (435., 4.1), 'g': (481., 3.69), 'gri': (617., 2.72),
 'r': (617., 2.72), 'i': (752., 2.08), 'z': (866., 1.63),
 'zy': (914., 1.49), 'y': (962., 1.37), 'J': (1235., 0.9),
 'H': (1662., 0.51)
}

def usage():
    print (__doc__)
    sys.exit()


def get_det_effic(files):
    """
    make some efficiency plots
    """

    fig,ax = subplots()
    ax.plot ([10,10],[0.,1.05],'k:',alpha=0.5)
    ax.plot ([5,5],[0.,1.05],'k:',alpha=0.5)
    ax.plot ([3,3],[0.,1.05],'k:',alpha=0.5)

    lam=zeros(len(files),dtype='float32')
    for i in range(len(files)):
        filter0 = files[i].split('_')[-1].split('/')[0]
        lam[i], Alam_fac = lam_dict[filter0]

    s= lam.argsort()
    files = list(array(files)[s])

    lims0 = " Search 90% Efficiency Minimum SNR Estimates: "
    lims1 = " Estimated Detection Efficiencies at 3, 5, and 10-sigma: "
    lims2 = " Estimated PS1 Association Efficiencies at 3, 5, and 10-sigma: "
    for file in files:
        fs = file.split('_')
        cam, filter = fs[-2],fs[-1].split('/')[0]
        sfile = 'redux_colibri_'+cam+'_'+filter+'/stack_'+cam+'.fits'
        hdr = getheader(sfile)
        maglim = hdr['MAGLIM']
        x=getdata(file,2)
        mfile = file.replace('photometry.fits','catalog_matches.fits')
        m=getdata(mfile)
        flx,dflx = x['flx'],x['dflx']
        j=dflx>0
        #j=(flx>dflx)
        snr = flx[j]/dflx[j]
        s = snr.argsort()[::-1]
        cp = (1+arange(len(s)))/len(s)

        h = (snr[s]>3)*(snr[s]<30)
        res = linfit(log(snr[s][h]),log(cp[h]))
        cp0 = maximum(cp,exp(res[0]+log(snr[s].clip(1))*res[1]))
        i0 = where((snr[s]>0.99)*(cp>0.9*cp0))[0][-1]
        snr_min0 = snr[s][i0]

        maglims = maglim - 2.5*log10(snr_min0/10)
        lims0 = lims0 + f" {snr_min0:.1f} ({filter}={maglims:.2f}) "

        res = interp1d(snr[s],cp/cp0)
        lim3 = min(1,res(3.))
        lim5 = min(1,max(lim3,res(5.)))
        lim10 = min(1,max(lim5,res(10.)))
        lims1 = lims1 + f"{lim3*100:.0f}% {lim5*100:.0f}% {lim10*100:.0f}% ({filter}) "

        p = ax.plot (snr[s],cp,label=filter)
        ax.plot (snr_min0,cp[i0],'ko')
        cat = m[1,j]>0
        cat_cp = cat[s].cumsum()/len(cat)
        res = linfit(log(snr[s][h]),log(cat_cp[h]))
        cat_cp0 = maximum(cp,exp(res[0]+log(snr[s].clip(1))*res[1]))
        res = interp1d(snr[s],cat_cp/cat_cp0)
        clim3 = min(1,res(3.))
        clim5 = min(1,max(clim3,res(5.)))
        clim10 = min(1,max(clim5,res(10.)))
        lims2 = lims2 + f"{clim3*100:.0f}% {clim5*100:.0f}% {clim10*100:.0f}% ({filter}) "
        far = 1 - cat_cp/cat_cp0
        g = snr[s]<10
        ax.plot (snr[s][g],far[g],color=p[0].get_color(),linestyle=':')
        #ax.plot (snr[s],cp0,color=p[0].get_color(),linestyle='--')
   
    print (lims0+'<BR>') 
    print (lims1+'<BR>')
    print (lims2)
    ax.set_xlim((1.,1.e3))
    ax.set_ylim((0,1.05))
    #ax.set_ylim((0.1,2.))
    #ax.loglog()
    ax.semilogx()
    ax.set_xlabel(r"Detection SNR $\sigma$",fontsize=14)
    ax.set_ylabel(r"N($>\sigma$)",fontsize=14)
    ax.set_title("Black circles are estimated 90% efficiency points")
    ax.legend()
    fig.savefig("det_effic.jpg")


if __name__ == "__main__":

    files=glob('redux_colibri_C*/photometry.fits')
    if (len(files)>0): get_det_effic(files)
