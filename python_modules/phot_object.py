from astropy.io.fits import getdata,getheader
from numpy import where,log10,log
from fit_wcs import ra2sex,dec2sex

class PhotObject:

    def __init__(self,file):
        self.srcdata = getdata(file,2)
        self.indx = 0

        hdr = getheader(file)
        self.t0 = hdr['TCSTART']
        self.t1 = hdr['TCSTOP']
        self.dt = hdr['TCEXP']
        self.filter = hdr['CFILTER']
        self.lam = hdr['LAMBDA']
        ncat=hdr.get('NCAT') or 0
        self.cnames = ["NONE"]
        for i in range(ncat):
            self.cnames.append( hdr["""CATID%d""" % (i+1)] )

    def setidx(self,srcid):
        self.indx = where(self.srcdata['srcid']==srcid)[0][0]

    def __str__(self):
        """
        print things nicely
        """
        x = self.srcdata[self.indx]
        ra, dec = x['RA'], x['DEC']
        ras, decs = ra2sex(ra), dec2sex(dec)

        str0 = f"Source {x['SRCID']}: {ras} {decs} ({ra:.5f}, {dec:.5f})"
        str0 += f" Catalog {self.cnames[x['catid']]},"

        str_oth = f" lam0: {self.lam:.0f} nm, A_lam: {x['A_lam']:.4f}, FWHM: {x['fwhm']:.1f}"
        strc = f"Common Epoch: {self.t0} - {self.t1} (dt={self.dt:.1f}s)"
        str_out = str0 + str_oth + '\n ' + strc

        return str_out

    def print_flx(self):
        """
        strings to print magnitude and flux information
        """
        x = self.srcdata[self.indx]
        flx, dflx = x['flx'], x['dflx']
        flxc, dflxc = x['flxc'], x['dflxc']
        flx_ref, dflx_ref = x['flx_ref'], x['dflx_ref']
        flx_sub, dflx_sub = x['flx_sub'], x['dflx_sub']

        if (flx>dflx): 
            mag, dmag = 25-2.5*log10(flx), 2.5/log(10)*dflx/flx
        else: 
            mag, dmag = 25-2.5*log10(max(0,flx)+3*dflx), 0
        if (flxc>dflxc): 
            magc, dmagc = 25-2.5*log10(flxc), 2.5/log(10)*dflxc/flxc
        else: 
            magc, dmagc = 25-2.5*log10(max(0,flxc)+3*dflxc), 0

        if (flx_sub>dflx_sub):
            mag_sub, dmag_sub = 25-2.5*log10(flx_sub), 2.5/log(10)*dflx_sub/flx_sub
        else: 
            mag_sub, dmag_sub = 25-2.5*log10(max(0,flx_sub)+3*dflx_sub), 0

        strf = f"{flx:.2f} +/- {dflx:.2f} (sci)"
        strf0 = f"{flx_ref:.2f} +/- {dflx_ref:.2f} (ref)"
        strdf = f"{flx_sub:.2f} +/- {dflx_sub:.2f} (diff)"
        strc = f"{flxc:.2f} +/- {dflxc:.2f} (common epoch)"
        flstr = f" Flux {self.filter}: "+strf+", "+strf0+", "+strdf+", "+strc

        if (dmag>0):
            magstr = f" Mag {self.filter}: {mag:.4f} +/- {dmag:.4f}"
        else:
            magstr = f" Mag {self.filter}: >{mag:.2f} (3-sigma)"
        if (dmag_sub>0):
            magstr += f", {mag_sub:.4f} +/- {dmag_sub:.4f} (diff)"
        else:
            magstr += f", >{mag_sub:.2f} (diff, 3-sigma)"
        if (dmagc>0):
            magstr += f", {magc:.4f} +/- {dmagc:.4f} (common epoch)"
        else:
            magstr += f", >{magc:.2f} (common epoch, 3-sigma)"

        str_out = flstr + '\n' + magstr
        return str_out
