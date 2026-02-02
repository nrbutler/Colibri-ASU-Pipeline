from get_filters import get_filters
from get_90conf_prob import get_90conf_prob
from numpy import arange,zeros,dot,log,log10,sqrt,newaxis,array,exp,unravel_index,outer,abs
from numpy.random import randn

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#plt.ioff()
plt.rcParams['figure.max_open_warning'] = 100

import os
from astropy.io.fits import getdata,writeto
fdir = os.environ['REDUX_BASE_DIR'] + '/calfiles/filters/'

def flx2mag(flx,dflx):
    mag = 25 - 2.5*log10(flx.clip(dflx)); 
    dmag = 2.5/log(10)*dflx/flx.clip(dflx)
    return mag,dmag

class PhotoZ:

    def __init__(self,lam1=350,lam2=2000,filters=array([0,1,1,1,1,1,1,1,0,0],dtype='bool'),z1=0,z2=10,dz=0.01,sys_err=0.1,verbose=True):
        """
        with lam in nm:
        """
        self.sys_err=sys_err
        self.ebv = 0
        self.verbose = verbose

        self.lam1 = lam1
        self.lam2 = lam2
        self.filters = filters
        self.filter_names = array(["B","g","r","i","z","y","J","H","gri","zy"])[filters]
        lam = 1.*lam1 + arange(1+int(lam2-lam1))

        Nz = int(1 + (z2-z1)/dz)
        self.zgrid = z1 + dz*arange(Nz)
        self.Nz = Nz

        self.Av_grid = arange(51)*0.1
        self.alpha_grid = -3.5+arange(51)*0.1

        # priors
        Avz = 0.5/(1+self.zgrid)**0.7 # Hayes, et al. 2011
        self.Pr_Av = exp(-outer(1/Avz,self.Av_grid)) / ((1+self.zgrid)*Avz)[:,newaxis]
        self.Pr_al = exp(-8*(self.alpha_grid+1)**2)

        lamf = getdata(fdir+'lambda.fits')[filters]
        self.filter_matr = getdata(fdir+'filter_matr.fits')[filters]
        self.lam0 = (self.filter_matr*lamf).sum(axis=1)
        self.extin_gal = getdata(fdir+'extin_grb1.fits')[0,filters]
        self.mw=getdata(fdir+'mw_data.fits')[:,:,:,filters]
        self.lmc=getdata(fdir+'lmc_data.fits')[:,:,:,filters]
        self.smc=getdata(fdir+'smc_data.fits')[:,:,:,filters]

        self.sed_filename = "sed_plot.jpg"
        self.redshift_filename = "redshift_plot.jpg"

        plt.ioff()
        self.saveplots=True

    def plotson(self):
        plt.ion()
        mpl.use('TkAgg')
        self.saveplots=False

    def sim_data(self,i0=550,Av=0.3,alpha=-1.,Av0=0.5,Hmag=21.,gain=1.e3,gal=3):
        """
        """
        i = abs(self.Av_grid-Av).argmin()
        j = abs(self.alpha_grid-alpha).argmin()
        if (gal==1): fmdl = self.mw[i0,i,j]
        elif (gal==2): fmdl = self.lmc[i0,i,j]
        else: fmdl = self.smc[i0,i,j]

        if (Av0>0): fmdl *= (10**(-0.4*Av0*self.extin_gal)*self.filter_matr).sum(axis=1)

        mag = -2.5*log10(fmdl)
        mag += Hmag - mag[-1]
        flx = 10**(-0.4*(mag-25))
        dflx = sqrt( 1 + flx/gain )
        return flx + randn(len(flx))*dflx, dflx

    def fit_flux(self,flx,dflx,Av0=0.):
        """ fitting by flux """
        self.ebv = Av0/3.08

        fac = (10**(-0.4*Av0*self.extin_gal)*self.filter_matr).sum(axis=1)
        flx1, dflx1 = flx.clip(0)/fac,dflx/fac
        wt = 0.*flx
        h = dflx>0
        if (h.sum()<2): return False

        wt[h]=1./(dflx1[h]**2+(self.sys_err*flx1[h])**2)
        wts = wt.sum()
        norm = (self.mw**2*wt).sum(axis=3); X1 = (self.mw*flx1*wt).sum(axis=3)**2; h = norm>1.e-10*wts; X1[h]/=norm[h]; X1[~h]=0
        norm = (self.lmc**2*wt).sum(axis=3); X2 = (self.lmc*flx1*wt).sum(axis=3)**2; h = norm>1.e-10*wts; X2[h]/=norm[h]; X1[~h]=0
        norm = (self.smc**2*wt).sum(axis=3); X3 = (self.smc*flx1*wt).sum(axis=3)**2; h = norm>1.e-10*wts; X3[h]/=norm[h]; X1[~h]=0

        Xm = max(X1.max(),X2.max(),X3.max())
        P1 = exp(0.5*(X1-Xm))*self.Pr_Av[:,:,newaxis]*self.Pr_al[newaxis,newaxis,:]
        P2 = exp(0.5*(X2-Xm))*self.Pr_Av[:,:,newaxis]*self.Pr_al[newaxis,newaxis,:]
        P3 = exp(0.5*(X3-Xm))*self.Pr_Av[:,:,newaxis]*self.Pr_al[newaxis,newaxis,:]

        self.prob1 = P1.sum(axis=2).sum(axis=1)
        self.prob2 = P2.sum(axis=2).sum(axis=1)
        self.prob3 = P3.sum(axis=2).sum(axis=1)

        i1 = P1.argmax()
        i1a,i1b,i1c = unravel_index(i1,X1.shape)
        self.z0_1,self.p2_1,self.p3_1 = self.zgrid[i1a],self.Av_grid[i1b],self.alpha_grid[i1c]
        self.fmdl_1 = self.mw[i1a,i1b,i1c]
        norm = dot(self.fmdl_1**2,wt)
        if (norm>0): norm = dot(self.fmdl_1,flx1*wt) / norm
        else: norm = 0
        self.fmdl_1 *= norm
        self.chi2_1m = dot( (flx1-self.fmdl_1)**2,wt )
        self.fmdl_1 *= fac
        i2 = P2.argmax()
        i2a,i2b,i2c = unravel_index(i2,X2.shape)
        self.z0_2,self.p2_2,self.p3_2 = self.zgrid[i2a],self.Av_grid[i2b],self.alpha_grid[i2c]
        self.fmdl_2 = self.lmc[i2a,i2b,i2c]
        norm = dot(self.fmdl_2**2,wt)
        if (norm>0): norm = dot(self.fmdl_2,flx1*wt) / norm
        else: norm = 0
        self.fmdl_2 *= norm
        self.chi2_2m = dot( (flx1-self.fmdl_2)**2,wt )
        self.fmdl_2 *= fac
        i3 = P3.argmax()
        i3a,i3b,i3c = unravel_index(i3,X3.shape)
        self.z0_3,self.p2_3,self.p3_3 = self.zgrid[i3a],self.Av_grid[i3b],self.alpha_grid[i3c]
        self.fmdl_3 = self.smc[i3a,i3b,i3c]
        norm = dot(self.fmdl_3**2,wt)
        if (norm>0): norm = dot(self.fmdl_3,flx1*wt) / norm
        else: norm = 0
        self.fmdl_3 *= norm
        self.chi2_3m = dot( (flx1-self.fmdl_3)**2,wt )
        self.fmdl_3 *= fac

        self.z1_1,self.z2_1 = get_90conf_prob(self.zgrid,self.prob1)
        self.z1_2,self.z2_2 = get_90conf_prob(self.zgrid,self.prob2)
        self.z1_3,self.z2_3 = get_90conf_prob(self.zgrid,self.prob3)

        if (self.verbose):
            print (f"  MW fit: {self.z1_1:5.2f} < z < {self.z2_1:<5.2f} (90% conf.) pars=({self.p2_1:.1f},{self.p3_1:.1f}) chi2/nu={self.chi2_1m:.1f}/{len(flx)-3}")
            print (f" LMC fit: {self.z1_2:5.2f} < z < {self.z2_2:<5.2f} (90% conf.) pars=({self.p2_2:.1f},{self.p3_2:.1f}) chi2/nu={self.chi2_2m:.1f}/{len(flx)-3}")
            print (f" SMC fit: {self.z1_3:5.2f} < z < {self.z2_3:<5.2f} (90% conf.) pars=({self.p2_3:.1f},{self.p3_3:.1f}) chi2/nu={self.chi2_3m:.1f}/{len(flx)-3}")

        return True
        
    def make_plots(self,flx,dflx,jfac=10**(-0.4*1.1),normalize=False):
        fig,ax = plt.subplots()
        axs = ax.errorbar(self.lam0,jfac*flx,xerr=0*self.lam0,yerr=jfac*dflx,fmt='bo',linestyle='None',markersize=5,mew=1,capsize=0)
        ax.plot (self.lam0,jfac*flx,':',alpha=0,label=r'            $A_v$     $\beta$       $z$        $\chi^2$')
        ax.plot(self.lam0,jfac*self.fmdl_1,label="""%-5s %5.1f %5.1f %6.2f %6.2f""" % ("MW",self.p2_1,self.p3_1,self.z0_1,self.chi2_1m))
        ax.plot(self.lam0,jfac*self.fmdl_2,label="""%-5s %5.1f %5.1f %6.2f %6.2f""" % ("LMC",self.p2_2,self.p3_2,self.z0_2,self.chi2_2m))
        ax.plot(self.lam0,jfac*self.fmdl_3,label="""%-5s %5.1f %5.1f %6.2f %6.2f""" % ("SMC",self.p2_3,self.p3_3,self.z0_3,self.chi2_2m))
        ax.semilogy()
        ax.set_xlabel("Wavelength [nm]",fontsize=14)
        ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]",fontsize=14)
        ymin = max(1.e-1*jfac,jfac*(flx-dflx).min()*0.9)
        ymax = (flx+dflx).max()*jfac
        ymax1 = ymin*exp(1.1*log(ymax/ymin))
        ax.set_ylim((ymin,ymax1))
        ypos = ymin*exp(1.05*log(ymax/ymin))
        for i in range(len(self.lam0)):
            an=plt.annotate(self.filter_names[i],xy=(self.lam0[i],ypos),ha='center')
        ax.legend(loc='lower right',frameon=False,fontsize=12,labelspacing=0.1)
        ax.set_title("""E[B-V]=%.3f, sys_err=%.2f""" % (self.ebv,self.sys_err))
        fig.tight_layout()
        if (self.saveplots):
            fig.savefig(self.sed_filename)
            plt.close(fig)

        fig,ax = plt.subplots()
        p1,p2,p3 = 1.*self.prob1,1.*self.prob2,1.*self.prob3
        if (normalize):
            p1/=p1.sum()
            p2/=p2.sum()
            p3/=p3.sum()
        axs = ax.plot(self.zgrid,p1,label=f" MW {self.z1_1:5.2f} < z < {self.z2_1:<5.2f}")
        ax.plot(self.zgrid,p2,label=f"LMC {self.z1_2:5.2f} < z < {self.z2_2:<5.2f}")
        ax.plot(self.zgrid,p3,label=f"SMC {self.z1_3:5.2f} < z < {self.z2_3:<5.2f}")
        ax.set_xlabel("Redshift z",fontsize=14)
        ax.set_ylabel("dP/dz",fontsize=14)
        ax.legend()
        za = min(self.z1_1,self.z1_2,self.z1_3)
        zb = max(self.z2_1,self.z2_2,self.z2_3)
        dz = zb-za
        ax.set_xlim(( max(0,za-dz/2),min(zb+dz/2,10) ))
        ax.set_title("90% conf. limits")
        fig.tight_layout()
        if (self.saveplots):
            fig.savefig(self.redshift_filename)
            plt.close(fig)
