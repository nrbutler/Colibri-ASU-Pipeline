from pei_av import pei_av
from madau_extin import madau_extin
from get_filters import get_filters
from get_90conf_prob import get_90conf_prob
from numpy import arange,zeros,dot,log,log10,sqrt,newaxis,array,exp,abs
from numpy.random import randn
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d

#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
#plt.ioff()

import os
from astropy.io.fits import getdata,writeto
fdir = os.environ['REDUX_BASE_DIR'] + '/calfiles/filters/'

def flx2mag(flx,dflx):
    mag = 25 - 2.5*log10(flx.clip(dflx)); 
    dmag = 2.5/log(10)*dflx/flx.clip(dflx)
    return mag,dmag

class PhotoZ:

    def __init__(self,lam1=350,lam2=2000,filters=array([0,1,1,1,1,1,1,1,0,0],dtype='bool'),nbin=10,z1=0,z2=10,dz=0.01,sys_err=0.1):
        """
        with lam in nm:
        """
        self.sys_err=sys_err
        self.ebv = 0

        self.lam1 = lam1
        self.lam2 = lam2
        self.filters = filters
        lam = 1.*lam1 + arange(1+int(lam2-lam1))

        Nz = int(1 + (z2-z1)/dz)
        self.zgrid = z1 + dz*arange(Nz)
        self.Nz = Nz

        # fitting priors
        self.Pr_alpha = -1.
        self.Pr_dalpha = 0.25
        self.Pr_Av = 0.5

        self.Av_min=0
        self.Av_max=5
        self.alpha_min=-3
        self.alpha_max=1

        if (filters.sum()==10):
            self.lamf, self.filter_matr = get_filters(lam,nbin=nbin)
            Nf,Nl = self.lamf.shape
            self.extin_grb1 = zeros((Nz,Nf,Nl),dtype='float32')
            self.extin_grb2 = zeros((Nz,Nf,Nl),dtype='float32')
            self.extin_grb3 = zeros((Nz,Nf,Nl),dtype='float32')
            self.igm_atten = zeros((Nz,Nf,Nl),dtype='float32')
            for i in range(Nz):
                self.extin_grb1[i] = pei_av(self.lamf*10/(1+self.zgrid[i]),gal=1)
                self.extin_grb2[i] = pei_av(self.lamf*10/(1+self.zgrid[i]),gal=2)
                self.extin_grb3[i] = pei_av(self.lamf*10/(1+self.zgrid[i]),gal=3)
                if (nbin<100): self.igm_atten[i] = interp1d(lam,gaussian_filter(madau_extin(lam*10,z=self.zgrid[i]),100/nbin))(self.lamf)
                else: self.igm_atten[i] = madau_extin(self.lamf*10,z=self.zgrid[i])

            writeto(fdir+'lambda.fits',self.lamf,overwrite=True)
            writeto(fdir+'filter_matr.fits',self.filter_matr,overwrite=True)
            writeto(fdir+'extin_grb1.fits',self.extin_grb1,overwrite=True)
            writeto(fdir+'extin_grb2.fits',self.extin_grb2,overwrite=True)
            writeto(fdir+'extin_grb3.fits',self.extin_grb3,overwrite=True)
            writeto(fdir+'igm_atten.fits',self.igm_atten,overwrite=True)
        else:
            self.lamf = getdata(fdir+'lambda.fits')[filters]
            self.filter_matr = getdata(fdir+'filter_matr.fits')[filters]
            self.extin_grb1 = getdata(fdir+'extin_grb1.fits')[:,filters]
            self.extin_grb2 = getdata(fdir+'extin_grb2.fits')[:,filters]
            self.extin_grb3 = getdata(fdir+'extin_grb3.fits')[:,filters]
            self.igm_atten = getdata(fdir+'igm_atten.fits')[:,filters]

        self.lam0 = (self.filter_matr*self.lamf).sum(axis=1)


    def sim_data(self,i0=550,Av=0.3,alpha=-1.,Av0=0.5,Hmag=21.,gain=1.e3,gal=3):
        """
        """
        if (gal==1): extin = self.extin_grb1
        elif (gal==2): extin = self.extin_grb2
        else: extin = self.extin_grb3

        x0 = 2.5*log10(self.lamf/620.)
        fmdl = ( 10**(-0.4*(Av*extin[i0]+Av0*self.extin_grb1[0]+alpha*x0))*self.igm_atten[i0]*self.filter_matr ).sum(axis=1)
        mag = -2.5*log10(fmdl)

        mag += Hmag - mag[-1]
        flx = 10**(-0.4*(mag-25))
        dflx = sqrt( 1 + flx/gain )
        return ( flx + randn(len(flx))*dflx ).clip(0),dflx

    def fit_flux(self,flx,dflx,gal=1,Av0=0.,dpmax=1.5,p20=0.3,p30=-1.,NITER=3):
        """ fitting by flux """
        p2 = p20 + 0*self.zgrid
        p3 = p30 + 0*self.zgrid
        self.ebv = Av0/3.08

        fac = -0.4*log(10)
        if (gal==1): extin = fac*self.extin_grb1
        elif (gal==2): extin = fac*self.extin_grb2
        else: extin = fac*self.extin_grb3
        extin0 = extin.mean(axis=2).mean(axis=1)
        extin -= extin0[:,newaxis,newaxis]

        wt = 1./(dflx**2 + (self.sys_err*flx.clip(0))**2)
        x0 = -log(self.lamf/620.)

        atten = 1.*self.igm_atten
        if (Av0>0): atten *= 10**(-0.4*Av0*self.extin_grb1[0])

        ep2 = p2[:,newaxis,newaxis]*extin
        fmdl = (exp(ep2+p3[:,newaxis,newaxis]*x0[newaxis,:,:])*atten*self.filter_matr).sum(axis=2)
        norm = dot(fmdl,flx*wt) / dot(fmdl**2,wt)
        fmdl *= norm[:,newaxis]
        dmagA = norm[:,newaxis]* (ep2*exp(ep2+p3[:,newaxis,newaxis]*x0[newaxis,:,:])*atten*self.filter_matr).sum(axis=2)
        x = norm[:,newaxis]* (x0*exp(ep2+p3[:,newaxis,newaxis]*x0[newaxis,:,:])*atten*self.filter_matr).sum(axis=2)

        # X2 = Sumi (fi-A*exp(p2*extin+p3*x))^2 * wt + 2*p2/p20 + (p3-p30)^2/dp3^2
        # mi = A*exp(p2*extin+p3*x
        # -0.5*dX2/dp2 = Sumi (fi-A*exp(p2*extin+p3*x))*mi*extin * wt - 1/p20
        # -0.5*dX2/dp3 = Sumi (fi-A*exp(p2*extin+p3*x))*mi*x * wt - (p3-p30)/dp3^2

        for i in range(NITER):
            v1 = dot(fmdl,flx*wt) - dot(fmdl**2,wt)
            v2 = dot(dmagA,flx*wt) - dot(dmagA*fmdl,wt) - 1./self.Pr_Av
            v3 = dot(x,flx*wt) - dot(x*fmdl,wt) - (p3-self.Pr_alpha)/self.Pr_dalpha**2
            #m00 = dot(fmdl**2,wt)
            m11 = 2*dot(fmdl**2,wt) - dot(fmdl,flx*wt)
            m12 = 2*dot(fmdl*dmagA,wt) - dot(dmagA,flx*wt)
            m22 = 2*dot(dmagA**2,wt) - dot(dmagA**2/fmdl,flx*wt)
            m13 = 2*dot(fmdl*x,wt) - dot(x,flx*wt)
            m23 = 2*dot(dmagA*x,wt) - dot(dmagA/fmdl*x,flx*wt)
            m33 = 2*dot(x**2,wt) - dot(x**2/fmdl,flx*wt) + 1./self.Pr_dalpha**2
            #m11 = dot(fmdl**2,wt)
            #m12 = dot(fmdl*dmagA,wt)
            #m22 = dot(dmagA**2,wt)
            #m13 = dot(fmdl*x,wt)
            #m23 = dot(dmagA*x,wt)
            #m33 = dot(x**2,wt) + 1./(self.Pr_dalpha*fac)**2
            det = m11*(m22*m33-m23**2)-m12*(m12*m33-m23*m13)
            det += m13*(m12*m23-m22*m13)
            mi11,mi12 = m22*m33-m23**2, m13*m23-m12*m33
            mi13,mi22 = m12*m23-m13*m22, m11*m33-m13**2
            mi23,mi33 = m13*m12-m11*m23, m11*m22-m12**2
            h = det>0
            dp1 = (mi11*v1+mi12*v2+mi13*v3)[h]/det[h]
            dp2 = (mi12*v1+mi22*v2+mi23*v3)[h]/det[h]
            dp3 = (mi13*v1+mi23*v2+mi33*v3)[h]/det[h]
            dp0 = sqrt(dp1**2+dp2**2+dp3**2)
            g = dp0>dpmax
            dp1[g] *= dpmax/dp0[g]
            dp2[g] *= dpmax/dp0[g]
            dp3[g] *= dpmax/dp0[g]
            p2[h] *= exp(dp2); p3[h] += dp3
            ep2 = p2[:,newaxis,newaxis]*extin
            fmdl = (exp(ep2+p3[:,newaxis,newaxis]*x0[newaxis,:,:])*atten*self.filter_matr).sum(axis=2)
            norm = dot(fmdl,flx*wt) / dot(fmdl**2,wt)
            fmdl *= norm[:,newaxis]
            dmagA = norm[:,newaxis]* (ep2*exp(ep2+p3[:,newaxis,newaxis]*x0[newaxis,:,:])*atten*self.filter_matr).sum(axis=2)
            x = norm[:,newaxis]* (x0*exp(ep2+p3[:,newaxis,newaxis]*x0[newaxis,:,:])*atten*self.filter_matr).sum(axis=2)

        self.chi2 = ( (flx[newaxis,:]-fmdl)**2/dflx**2 ).sum(axis=1)
        i = self.chi2.argmin()
        self.fmdl = norm[i]*exp(ep2[i]+p3[i]*x0)*atten[i]
        self.fmdl0 = fmdl[i]
    
        return 25+log(norm)/fac+p2*extin0/fac,p2,p3

    def eval_model(self,flx,dflx,iz=350,alpha=-1,Av=0.3,Av0=0.,gal=3):
        """
        """
        x0 = 2.5*log10(self.lamf)
        fac = 0.4*log(10.) 
        p2 = -fac*Av
        p3 = -fac*alpha
        wt = 1./(dflx**2+(self.sys_err*flx.clip(0))**2)

        if (gal==1): extin = self.extin_grb1[iz]
        elif (gal==2): extin = self.extin_grb2[iz]
        else: extin = self.extin_grb3[iz]
       
        atten = 10**(-0.4*Av0*self.extin_grb1[0])*self.igm_atten[iz] 
        fmdl = (exp(p2*extin+p3*x0)*atten*self.filter_matr).sum(axis=1)
        norm = dot(fmdl,flx*wt) / dot(fmdl**2,wt)
        self.mdl_chi2 = dot((flx-norm*fmdl)**2,wt)
        return norm*fmdl

    def fit_3flux(self,flx,dflx,Av0=0.,NITER=3,p20=0.3,p30=-1.):
        p1,p2,p3 = self.fit_flux(flx,dflx,Av0=Av0,gal=1,NITER=NITER,p20=p20,p30=p30)
        i1 = self.chi2.argmin()
        self.z0_1,self.p1_1,self.p2_1,self.p3_1 = self.zgrid[i1],p1[i1],p2[i1],p3[i1]
        chi2_1, self.fmdl_1 = 1.*self.chi2, 1.*self.fmdl0
        p1,p2,p3 = self.fit_flux(flx,dflx,Av0=Av0,gal=2,NITER=NITER,p20=p20,p30=p30)
        i2 = self.chi2.argmin()
        self.z0_2,self.p1_2,self.p2_2,self.p3_2 = self.zgrid[i2],p1[i2],p2[i2],p3[i2]
        chi2_2, self.fmdl_2 = 1.*self.chi2, 1.*self.fmdl0
        p1,p2,p3 = self.fit_flux(flx,dflx,Av0=Av0,gal=3,NITER=NITER,p20=p20,p30=p30)
        i3 = self.chi2.argmin()
        self.z0_3,self.p1_3,self.p2_3,self.p3_3 = self.zgrid[i3],p1[i3],p2[i3],p3[i3]
        chi2_3, self.fmdl_3 = 1.*self.chi2, 1.*self.fmdl0
        self.chi2_1m, self.chi2_2m, self.chi2_3m = chi2_1[i1],chi2_2[i2],chi2_3[i3]
        chi0 = min(self.chi2_1m,self.chi2_2m,self.chi2_3m)
        self.prob1 = exp(-0.5*(chi2_1-chi0))
        self.prob2 = exp(-0.5*(chi2_2-chi0))
        self.prob3 = exp(-0.5*(chi2_3-chi0))
        self.z1_1, self.z2_1 = get_90conf_prob(self.zgrid,self.prob1)
        self.z1_2, self.z2_2 = get_90conf_prob(self.zgrid,self.prob2)
        self.z1_3, self.z2_3 = get_90conf_prob(self.zgrid,self.prob3)

        print (f"  MW fit: {self.z1_1:5.2f} < z < {self.z2_1:<5.2f} (90% conf.) pars=({self.p1_1:.1f},{self.p2_1:.1f},{self.p3_1:.1f}) chi2/nu={self.chi2_1m:.1f}/{len(flx)-3}")
        print (f" LMC fit: {self.z1_2:5.2f} < z < {self.z2_2:<5.2f} (90% conf.) pars=({self.p1_2:.1f},{self.p2_2:.1f},{self.p3_2:.1f}) chi2/nu={self.chi2_2m:.1f}/{len(flx)-3}")
        print (f" SMC fit: {self.z1_3:5.2f} < z < {self.z2_3:<5.2f} (90% conf.) pars=({self.p1_3:.1f},{self.p2_3:.1f},{self.p3_3:.1f}) chi2/nu={self.chi2_3m:.1f}/{len(flx)-3}")
        
    def make_plots(self,flx,dflx,jfac=10**(-0.4*1.1),normalize=False):
        fig,ax = plt.subplots()
        axs = ax.errorbar(self.lam0,jfac*flx,xerr=0*self.lam0,yerr=jfac*dflx,fmt='bo',linestyle='None',markersize=5,mew=1,capsize=0)
        ax.plot (self.lam0,jfac*flx,':',alpha=0,label=r'            $A_v$     $\beta$       $z$        $\chi^2$')
        ax.plot(self.lam0,jfac*self.fmdl_1,label="""%-5s %5.1f %5.1f %6.2f %6.2f""" % ("MW",self.p2_1,self.p3_1,self.z0_1,self.chi2_1m))
        ax.plot(self.lam0,jfac*self.fmdl_2,label="""%-5s %5.1f %5.1f %6.2f %6.2f""" % ("LMC",self.p2_2,self.p3_2,self.z0_2,self.chi2_2m))
        ax.plot(self.lam0,jfac*self.fmdl_3,label="""%-5s %5.1f %5.1f %6.2f %6.2f""" % ("SMC",self.p2_3,self.p3_3,self.z0_3,self.chi2_2m))
        ax.loglog()
        ax.set_xlabel("Wavelength [nm]",fontsize=14)
        ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]",fontsize=14)
        ax.set_ylim(( max(1.e-1*jfac,jfac*(flx-dflx).min()*0.9),(flx+dflx).max()*1.1*jfac ))
        ax.legend(loc='lower right',frameon=False,fontsize=12,labelspacing=0.1)
        ax.set_title("""E[B-V]=%.3f, sys_err=%.2f""" % (self.ebv,self.sys_err))
        fig.savefig("sed_plot.jpg")

        fig,ax = plt.subplots()
        if (normalize):
            axs = ax.plot(self.zgrid,self.prob1/self.prob1.sum(),label=f" MW {self.z1_1:5.2f} < z < {self.z2_1:<5.2f}")
            ax.plot(self.zgrid,self.prob2/self.prob2.sum(),label=f"LMC {self.z1_2:5.2f} < z < {self.z2_2:<5.2f}")
            ax.plot(self.zgrid,self.prob3/self.prob3.sum(),label=f"SMC {self.z1_3:5.2f} < z < {self.z2_3:<5.2f}")
        else:
            axs = ax.plot(self.zgrid,self.prob1,label=f" MW {self.z1_1:5.2f} < z < {self.z2_1:<5.2f}")
            ax.plot(self.zgrid,self.prob2,label=f"LMC {self.z1_2:5.2f} < z < {self.z2_2:<5.2f}")
            ax.plot(self.zgrid,self.prob3,label=f"SMC {self.z1_3:5.2f} < z < {self.z2_3:<5.2f}")
        ax.set_xlabel("Redshift z",fontsize=14)
        ax.set_ylabel("dP/dz",fontsize=14)
        ax.legend()
        za = min(self.z1_1,self.z1_2,self.z1_3)
        zb = max(self.z2_1,self.z2_2,self.z2_3)
        dz = zb-za
        ax.set_xlim(( max(0,za-dz/2),min(zb+dz/2,10) ))
        fig.savefig("redshift_plot.jpg")

    def build_gridfitter(self):
        x=PhotoZ(filters=array([1,1,1,1,1,1,1,1,1,1],dtype='bool'),nbin=300)

        Av_grid=0.1*arange(51)
        alpha_grid=-3.5+arange(51)*0.1
        Nz = len(x.zgrid)
        atn = x.igm_atten
        Nf = len(x.filter_matr)
        x0 = -log(x.lamf/620.)

        data = zeros((Nz,51,51,Nf),dtype='float32')

        extin = -0.4*log(10)*x.extin_grb1; file='mw_data.fits'
        for i in range(51):
            Av = Av_grid[i]
            for j in range(51):
                alpha = alpha_grid[j]
                data[:,i,j] = (exp(Av*extin+alpha*x0[newaxis,:,:])*atn*x.filter_matr).sum(axis=2)
        writeto(fdir+file,data,overwrite=True)

        extin = -0.4*log(10)*x.extin_grb2; file='lmc_data.fits'
        for i in range(51):
            Av = Av_grid[i]
            for j in range(51):
                alpha = alpha_grid[j]
                data[:,i,j] = (exp(Av*extin+alpha*x0[newaxis,:,:])*atn*x.filter_matr).sum(axis=2)
        writeto(fdir+file,data,overwrite=True)

        extin = -0.4*log(10)*x.extin_grb3; file='smc_data.fits'
        for i in range(51):
            Av = Av_grid[i]
            for j in range(51):
                alpha = alpha_grid[j]
                data[:,i,j] = (exp(Av*extin+alpha*x0[newaxis,:,:])*atn*x.filter_matr).sum(axis=2)
        writeto(fdir+file,data,overwrite=True)
