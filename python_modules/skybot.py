from astropy.io import ascii
import requests

def run_skybot(ra_deg, dec_deg, rad_deg, juldate, dt_hrs=0., vmax=24.):
    """
      run a skybot cone search to find minor planets
    """

    server='http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php'

    r = requests.get(server, params={'RA': ra_deg, 'DEC': dec_deg,
        'SR': rad_deg, 'EPOCH': str(juldate),
        '-observer': '679', '-output': 'basic', '-mime': 'text'},
       timeout=20)

    # motion dra=dra*cos(dec), ddec is in arcsec/hr , posunc is in arcsec

    try:
        results = ascii.read(r.text, delimiter='|', names=('number', 'name', 'ra', 'dec', 'type', 'V', 'posunc', 'd', 'dra', 'ddec', 'dg', 'dh'))
    except:
        results=[]

    nsources=len(results)
    if (nsources>0):
        ras = results['ra']
        dra = results['dra']
        decs = results['dec']
        ddec = results['ddec']
        err = results['posunc']
        vmag = results['V']
        name = results['name']
        for i in range(nsources):
            h,m,s = ras[i].split()
            ra=15.*(float(h)+float(m)/60.+float(s)/3600.)
            h,m,s = decs[i].split()
            if (float(h)>=0): dec=float(h)+float(m)/60.+float(s)/3600.
            else: dec=float(h)-float(m)/60.-float(s)/3600.
            if( vmag[i]<vmax): print (ra,dec,vmag[i],name[i].replace(' ','_'),dra[i]*dt_hrs,ddec[i]*dt_hrs,err[i])
