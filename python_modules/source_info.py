#!/usr/bin/python3
"""
 source_info.py <candidate_list.txt>
"""

import sys,os

from astropy.io.fits import getdata,getheader
from glob import glob
from fit_wcs import ra2sex, dec2sex
from numpy import loadtxt,array,log,log10,zeros,where
from photoz import PhotoZ

from multiprocessing import Pool

# dictionary of lambda, A_lambda/E[B-V]
lam_dict = {
 'B': (435., 4.1), 'g': (481., 3.69), 'gri': (617., 2.72),
 'r': (617., 2.72), 'i': (752., 2.08), 'z': (866., 1.63),
 'zy': (914., 1.49), 'y': (962., 1.37), 'J': (1235., 0.9),
 'H': (1662., 0.51)
}
filters_available = ['B', 'g', 'r', 'i', 'z', 'y', 'J', 'H', 'gri', 'zy']

def html_style(fx,fy,sz=461):
    ostr = f"""<style>
    .container {{ width: 200px; height: 200px; overflow: hidden; position: relative; display: inline-block; border: 1px solid #000000; }}
    .container img {{ position: absolute; top: {100-(1-fy)*sz:.0f}px; left: {100-fx*sz:.0f}px; }}
    .overlay-text {{ position: absolute; top: 0%; left: 0%; color: white; font-size: 18px; text-align: center; background-color: rgba(0, 0, 0, 0.5); padding: 2px; }}
    .circle-overlay {{ position: absolute; width: 24px; height: 24px; border-radius: 30%; border: 2px solid red; top: 86px; left: 86px; }}
    </style>"""
    return ostr

file0='redux_colibri_AL/stack_AL_dir/catalog.fits'
dat0 = getdata(file0)
xpos0, ypos0 = dat0['X_IMAGE'], dat0['Y_IMAGE']
dat0=0

files=glob('redux_colibri_C?_*/photometry.fits')
Nf = len(files)

filters=[]
lam=zeros(Nf,dtype='float32')
Alam_fac=zeros(Nf,dtype='float32')
for i in range(Nf):
    i0 = files[i].rfind('_')
    i1 = files[i].find('/')
    filter0 = files[i][i0+1:i1]
    filters.append(filter0)
    lam[i], Alam_fac[i] = lam_dict[filter0]

s = lam.argsort()
filters = array(filters)
lam = lam[s]; Alam_fac = Alam_fac[s]; filters = list(filters[s]); files = list(array(files)[s])
files.insert(0,'redux_colibri_AL/photometry.fits')

filters_present = []
for filter in filters_available:
    filters_present.append( filter in filters )

filters_present = array(filters_present)
pz = PhotoZ(filters=filters_present,verbose=False)

filters.insert(0,'Chi2')

npz=0
npz_max=1000

def flx2mag(flx,dflx):
    if (dflx<=0): return 0.0, 0.0, " "
    elif (flx>2*dflx):
        return 25-2.5*log10(flx), 2.5/log(10)*dflx/flx, " "
    else:
        return 25-2.5*log10(max(0,flx)+3*dflx), 0.0, ">"

def print_mags(x0,i=0):
     """
     """
     x = x0[i]
     m, dm, m_str = flx2mag(x['flx'],x['dflx'])
     ms, dms, ms_str = flx2mag(x['flx_sub'],x['dflx_sub'])
     mc, dmc, mc_str = flx2mag(x['flxc'],x['dflxc'])
     mr, dmr, mr_str = flx2mag(x['flx_ref'],x['dflx_ref'])
     alam = x['A_lam']
     return  f"{m_str}{m:<7.4f} {dm:<7.4f} {ms_str}{ms:<7.4f} {dms:<7.4f} {mc_str}{mc:<7.4f} {dmc:<7.4f} {mr_str}{mr:<7.4f} {dmr:<7.4f} {alam:<7.4f}"


x=[]
for file in files:
    x.append( getdata(file,2) )

hdr = getheader(files[1])
ncat=hdr.get('NCAT') or 0
cnames = ["NONE"]
for i in range(ncat):
    cnames.append( hdr["""CATID%d""" % (i+1)] )

x0 = getdata(files[0],1)
t00,t11,gt00,gt11 = x0['DATE-OBS'][0],x0['DATE-OBE'][0],x0['T0'][0],x0['T1'][0]

t0,t1 = hdr['TCSTART'], hdr['TCSTOP']
gt0,gt1 = hdr['GCSTART'], hdr['GCSTOP']

i0=[]; j0=[]; ts=[]; te=[]; pdir=[]
for file in files:
    hdr = getheader(file,1)
    i0.append(hdr['TI0'])
    j0.append(hdr['TJ0'])
    ts.append(hdr['TS'])
    te.append(hdr['TE'])
    pdir.append(hdr['PDIR'])

def print_all_mags(idp,snr_cut=10,q=0):
    qstr="Bright"
    if (q==1): qstr="Faint"
    elif (q==2): qstr="Questionable"
    i = where(x[1]['srcid']==idp)[0][0]
    x0 = x[1][i]
    id,r,d,f,cid = x0['srcid'],x0['ra'],x0['dec'],x0['fwhm'],x0['catid']
    for j in range(2,Nf+1):
        cid1 = x[j][i]['catid']
        if (cid1>cid): cid=cid1
       
    ebv = x0['A_lam']/Alam_fac[0]
    rs,ds = ra2sex(r), dec2sex(d)

    style=""
    for j in range(1,Nf+1):
        onframe = x[j]['dflx'][i]>0
        bright = x[j]['flx'][i] > snr_cut*x[j]['dflx'][i]
        if (not (onframe and bright)): continue
        xi=i0[j]+int((x[j]['x'][i]-i0[j]-0.5)/ts[j])*ts[j]
        yi=j0[j]+int((x[j]['y'][i]-j0[j]-0.5)/ts[j])*ts[j]
        fx,fy = (x[j]['x'][i]-xi+te[j])/(ts[j]+2*te[j]),(x[j]['y'][i]-yi+te[j])/(ts[j]+2*te[j])
        style = html_style(fx,fy)
        break

    if (len(style)==0):
        for j in range(1,Nf+1):
            onframe = x[j]['dflx'][i]>0
            if (not onframe): continue
            xi=i0[j]+int((x[j]['x'][i]-i0[j]-0.5)/ts[j])*ts[j]
            yi=j0[j]+int((x[j]['y'][i]-j0[j]-0.5)/ts[j])*ts[j]
            fx,fy = (x[j]['x'][i]-xi+te[j])/(ts[j]+2*te[j]),(x[j]['y'][i]-yi+te[j])/(ts[j]+2*te[j])
            style = html_style(fx,fy)
            break

    html_str = f"<HTML><HEAD><TITLE>Source {idp}</TITLE> {style} </HEAD><BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\"><PRE>\n"

    # data section, str_out
    str0=f"Source {id}: {rs} {ds} ({r:.5f}, {d:.5f}) ({xpos0[i]:.1f},{ypos0[i]:.1f}) {qstr}\n Catalog {cnames[cid]} E[B-V]: {ebv:.4f}, FWHM: {f:.1f}\n"
    strt = f"Full Time Range: {t00} - {t11} ({gt00/3600:.4f} - {gt11/3600:.4f} hours)\n"
    strc = f"Common Epoch: {t0} - {t1} ({gt0/3600:.4f} - {gt1/3600:.4f} hours)\n"
    str_out = str0 + strt + strc
    str_out = str_out + ("-"*79) + '\n'
    str_out = str_out + "filter  mag            mag_subtraction  mag_common       mag_reference   A_lam" + '\n'
    for j in range(Nf):
        str_out = str_out + f" {filters[1+j]:<4s}" + print_mags(x[j+1],i) + '\n'
    str_out = str_out + ("-"*79) + '\n'

    html_str = html_str + str_out + '</PRE>\n'

    for j in range(Nf+1):
        onframe = x[j]['dflx'][i]>0
        if (not onframe): continue
        xi=i0[j]+int((x[j]['x'][i]-i0[j]-0.5)/ts[j])*ts[j]
        yi=j0[j]+int((x[j]['y'][i]-j0[j]-0.5)/ts[j])*ts[j]
        istr = f"<A HREF=\"../{pdir[j]}/thumb_{xi}_{yi}_.html\" TARGET=\"_blank\">"
        istr = istr + f"<div class=\"container\"><IMG SRC=\"../{pdir[j]}/thumb_{xi}_{yi}_.jpg\"><div class=\"overlay-text\">{filters[j]}</div><div class=\"circle-overlay\"></div></div></A>\n"
        html_str = html_str + istr
    html_str = html_str + '<BR>\n'

    for j in range(Nf+1):
        onframe = x[j]['dflx'][i]>0
        if (not onframe): continue
        xi=i0[j]+int((x[j]['x'][i]-i0[j]-0.5)/ts[j])*ts[j]
        yi=j0[j]+int((x[j]['y'][i]-j0[j]-0.5)/ts[j])*ts[j]
        istr = f"<A HREF=\"../{pdir[j]}/thumb_{xi}_{yi}_.html\" TARGET=\"_blank\">"
        istr = istr + f"<div class=\"container\"><IMG SRC=\"../{pdir[j]}/ps1_thumb_{xi}_{yi}_.jpg\"><div class=\"overlay-text\">PS1{filters[j]}</div><div class=\"circle-overlay\"></div></div></A>\n"
        html_str = html_str + istr
    html_str = html_str + '<BR>\n'

    for j in range(Nf+1):
        onframe = x[j]['dflx'][i]>0
        if (not onframe): continue
        xi=i0[j]+int((x[j]['x'][i]-i0[j]-0.5)/ts[j])*ts[j]
        yi=j0[j]+int((x[j]['y'][i]-j0[j]-0.5)/ts[j])*ts[j]
        istr = f"<A HREF=\"../{pdir[j]}/diff_thumb_{xi}_{yi}_.html\" TARGET=\"_blank\">"
        istr = istr + f"<div class=\"container\"><IMG SRC=\"../{pdir[j]}/diff_thumb_{xi}_{yi}_.jpg\"><div class=\"overlay-text\">Sub{filters[j]}</div><div class=\"circle-overlay\"></div></div></A>\n"
        html_str = html_str + istr
    html_str = html_str + '<BR>\n'

    jstr='Ascii text files: '
    for j in range(1,Nf+1):
        onframe = x[j]['dflx'][i]>0
        bright = x[j]['flx'][i] > snr_cut*x[j]['dflx'][i]
        if (not (onframe and bright)): continue
        istr = f"<A HREF=\"../{pdir[j]}/lc_{idp}.jpg\" TARGET=\"_blank\"><IMG WIDTH=300 SRC=\"../{pdir[j]}/lc_{idp}.jpg\"></A>\n"
        jstr = jstr + f"<A HREF=\"../{pdir[j]}/lc_{idp}.txt\" TARGET=\"_blank\">lc_{filters[j]:<4s}_{idp}.txt</A> &nbsp; &nbsp; "
        html_str = html_str + istr

    html_str = html_str + "<BR>" + jstr

    if (cid==0):   # new sources only

        # try to do photoz analyasis
        flx = zeros(Nf,dtype='float32')
        dflx = zeros(Nf,dtype='float32')
        for j in range(Nf):
            flx[j] = x[j+1]['flxc'][i]
            dflx[j] = x[j+1]['dflxc'][i]

        global npz
        fit=0
        if (npz<=npz_max):
            fit = pz.fit_flux(flx,dflx,Av0=ebv*3.08)
            npz += 1
        #else: print (f"Warning: {npz}>{npz_max}")

        if (fit):
            pz.sed_filename = 'new_sources/'+f"sed_{idp}.jpg"
            pz.redshift_filename = 'new_sources/'+f"redshift_{idp}.jpg"
            pz.make_plots(flx,dflx)
            istr = f"<P><A HREF=\"sed_{idp}.jpg\" TARGET=\"_blank\"><IMG WIDTH=300 SRC=\"sed_{idp}.jpg\"></A>"
            istr = istr + f"<A HREF=\"redshift_{idp}.jpg\" TARGET=\"_blank\"><IMG WIDTH=300 SRC=\"redshift_{idp}.jpg\"></A>\n"
            html_str = html_str + istr
        #else: print (f"Not making sed and photo-z plots for Source {idp}")

    html_str = html_str + '<BR></BODY></HTML>\n'

    return html_str


if __name__ == "__main__":
    """
    just run
    """
    if (len(sys.argv)<2): usage()

    clist = sys.argv[1]
    if (os.path.exists(clist)==False): usage()

    snr_cut=10.
    if (len(sys.argv)>2): snr_cut=float(sys.argv[2])

    def make_file(args):
        idp0,q0 = args
        str_out = print_all_mags(idp0,snr_cut=snr_cut,q=q0)
        ofile = open(f"new_sources/source_{idp0}.html","w")
        ofile.write(str_out)
        ofile.close()

    try: idp,q = loadtxt(clist,unpack=True,dtype='int32',ndmin=2)
    except: idp=[]

    if (len(idp)>0):
        pool = Pool()
        pool.map( make_file, zip(idp, q) )
