from astropy.table import Table
import numpy as np

# pixel scale is 0.25 arcsec
pixscale = 0.25

# table of rings info
rings = Table.read("""zone projcell nband dec dec_min dec_max xcell ycell crpix1 crpix2
13 487 72 -38.0 -39.98569922229713 -35.98748871852416 6305 6279 239.5 238.0
14 559 76 -34.0 -35.98748871852416 -31.988947843740302 6265 6274 237.5 237.5
15 635 79 -30.0 -31.988947843740302 -27.990585819085737 6274 6269 239.5 241.5
16 714 82 -26.0 -27.990585819085737 -23.991856317278298 6255 6265 241.5 240.5
17 796 84 -22.0 -23.991856317278298 -19.99333106003115 6279 6261 241.0 242.0
18 880 86 -18.0 -19.99333106003115 -15.994464100667857 6274 6258 241.5 238.0
19 966 88 -14.0 -15.994464100667857 -11.995671545879146 6242 6254 240.5 239.5
20 1054 89 -10.0 -11.995671545879146 -7.997014822849273 6248 6250 240.0 242.0
21 1143 89 -6.0 -7.997014822849273 -3.998086931795508 6291 6247 237.5 240.0
22 1232 90 -2.0 -3.998086931795508 1.3877787807814457e-17 6240 6243 240.0 242.0
23 1322 90 2.0 1.3877787807814457e-17 3.998086931795508 6240 6243 240.0 242.0
24 1412 89 6.0 3.998086931795508 7.997014822849273 6291 6247 237.5 240.0
25 1501 89 10.0 7.997014822849273 11.995671545879146 6248 6250 240.0 242.0
26 1590 88 14.0 11.995671545879146 15.994464100667857 6242 6254 240.5 239.5
27 1678 86 18.0 15.994464100667857 19.99333106003115 6274 6258 241.5 238.0
28 1764 84 22.0 19.99333106003115 23.991856317278298 6279 6261 241.0 242.0
29 1848 82 26.0 23.991856317278298 27.990585819085737 6255 6265 241.5 240.5
30 1930 79 30.0 27.990585819085737 31.988947843740302 6274 6269 239.5 241.5
31 2009 76 34.0 31.988947843740302 35.98748871852416 6265 6274 237.5 237.5
32 2085 72 38.0 35.98748871852416 39.98569922229713 6305 6279 239.5 238.0
33 2157 68 42.0 39.98569922229713 43.98384988850235 6320 6284 239.5 240.0
34 2225 64 46.0 43.98384988850235 47.98179099058205 6307 6289 238.0 242.0
35 2289 60 50.0 47.98179099058205 51.97972846465295 6261 6295 241.0 239.5
36 2349 55 54.0 51.97972846465295 55.97676677254534 6283 6302 239.0 242.0
37 2404 50 58.0 55.97676677254534 59.97411834356914 6278 6311 238.5 237.5
38 2454 45 62.0 59.97411834356914 63.97005521318337 6240 6319 240.0 241.0
39 2499 39 66.0 63.97005521318337 67.96475284198297 6307 6333 239.5 240.0
40 2538 33 70.0 67.96475284198297 71.95817764373305 6365 6350 238.5 240.0
41 2571 27 74.0 71.95817764373305 75.94974140116577 6413 6371 240.5 242.0
42 2598 21 78.0 75.94974140116577 79.93968190571101 6452 6398 240.0 242.0
43 2619 15 82.0 79.93968190571101 83.93610604798963 6481 6429 241.0 242.0
44 2634 9 86.0 83.93610604798963 87.96940975734594 6501 6419 239.0 239.5
45 2643 1 90.0 87.96940975734594 90.0001 6240 6240 239.5 240.0
""", format="ascii.csv", delimiter=" ")

dec_limit = rings['dec_min'].min()

def findskycell(ra, dec):

    """Given input arrays RA, DEC (degrees), returns an astropy table with the skycell info

    RA and Dec can be scalars or arrays.
    This returns just the best sky cell for each ra/dec (the one where the position is closest to the center).
    This uses the rings table for info on the tessellation.
    The returned table has these columns:
    Column      Value
    ra          Input position
    dec         
    index       0-based index into original array
    projcell    Projection cell (0 if outside coverage)
    subcell     Subcell (0..99)
    crval1      Sky position of reference pixel for projection cell
    crval2      
    crpix1      Reference pixel in this skycell
    crpix2      
    x           Source pixel position in this skycell
    y           
    iring       Index into rings array
    """

    if np.isscalar(ra) and np.isscalar(dec):
        return _findskycell_array(np.array([ra]),np.array([dec]))
    if len(ra) == len(dec):
        return _findskycell_array(np.asarray(ra),np.asarray(dec))
    else:
        raise ValueError("ra and dec must both be scalars or be matching length arrays")


def _findskycell_array(ra, dec):

    """Internal function: ra and dec are known to be arrays"""

    if ra.ndim != 1 or dec.ndim != 1:
        raise ValueError("ra and dec must be 1-D arrays")
    index = np.arange(len(ra),dtype=int)
    # find dec zone where rings.dec_min <= dec < rings.dec_max
    idec = np.searchsorted(rings['dec_max'], dec)

    # special handling at pole where overlap is complicated
    # do extra checks for top 2 rings
    # always start with the ring just below the pole
    nearpole = np.where(idec >= len(rings)-2)
    idec[nearpole] = len(rings)-2

    nband = rings['nband'][idec]
    # get normalized RA in range 0..360
    nra = ra % 360.0
    ira = np.rint(nra*nband/360.0).astype(int) % nband

    projcell = rings['projcell'][idec] + ira
    dec_cen = rings['dec'][idec]
    ra_cen = ira*360.0/nband

    # locate subcell(s) within the projection cell

    # use tangent project to get pixel offsets
    x, y = sky2xy_tan(nra, dec, ra_cen, dec_cen)

    pad = 480

    if len(nearpole[0]) > 0:
        # handle the points near the pole (if any)
        # we know that this "ring" has only a single field
        projcell2 = rings['projcell'][-1]
        dec_cen2 = rings['dec'][-1]
        ra_cen2 = 0.0
        x2, y2 = sky2xy_tan(nra[nearpole], dec[nearpole], ra_cen2, dec_cen2)
        # compare x,y and x2,y2 to image sizes to select best image
        # returns a Boolean array with true for values where 2nd image is better
        use2 = poleselect(x[nearpole], y[nearpole], x2, y2, rings[-2], rings[-1], pad)
        if use2.any():
            # slightly obscure syntax here makes this work even if ra, dec are multi-dimensional
            wuse2 = np.where(use2)[0]
            w2 = tuple(xx[wuse2] for xx in nearpole)
            idec[w2] = len(rings)-1
            nband[w2] = 1
            ira[w2] = 0
            projcell[w2] = projcell2
            dec_cen[w2] = dec_cen2
            ra_cen[w2] = ra_cen2
            x[w2] = x2[wuse2]
            y[w2] = y2[wuse2]

    # compute the subcell from the pixel location
    px = rings['xcell'][idec]-pad
    py = rings['ycell'][idec]-pad
    k = np.rint(4.5+x/px).astype(int).clip(0,9)
    j = np.rint(4.5+y/py).astype(int).clip(0,9)
    subcell = 10*j + k

    # get pixel coordinates within the skycell image
    crpix1 = rings['crpix1'][idec] + px*(5-k)
    crpix2 = rings['crpix2'][idec] + py*(5-j)
    ximage = x + crpix1
    yimage = y + crpix2
    iring = idec

    # insert zeros where we are below lowest dec_min
    w = np.where(dec < dec_limit)
    projcell[w] = 0
    subcell[w] = 0
    crpix1[w] = 0
    crpix2[w] = 0
    ximage[w] = 0
    yimage[w] = 0

    # return an astropy Table
    return Table([ra,dec,index,projcell,subcell,ra_cen,dec_cen,crpix1,crpix2,ximage,yimage,iring],
                 names="ra,dec,index,projcell,subcell,crval1,crval2,crpix1,crpix2,x,y,iring".split(","))


def poleselect(x1, y1, x2, y2, rings1, rings2, pad):

    """Compares x,y values from 2 images to determine which is best

    Returns boolean array with True where x2,y2 is best
    """

    nx1 = 10*(rings1['xcell']-pad)+pad
    ny1 = 10*(rings1['ycell']-pad)+pad
    nx2 = 10*(rings2['xcell']-pad)+pad
    ny2 = 10*(rings2['ycell']-pad)+pad
    # compute minimum distances to image edges
    # note negative values are off edge
    d1 = np.minimum(np.minimum(x1+nx1//2,nx1//2-1-x1), np.minimum(y1+ny1//2,ny1//2-1-y1))
    d2 = np.minimum(np.minimum(x2+nx2//2,nx2//2-1-x2), np.minimum(y2+ny2//2,ny2//2-1-y2))
    return (d1 < d2)


def sky2xy_tan(ra, dec, ra_cen, dec_cen, crpix=(0.0,0.0)):

    """Convert RA,Dec sky position (degrees) to X,Y pixel position

    ra[n], dec[n] are input arrays in degrees
    ra_cen[n], dec_cen[n] are image centers in degrees
    crpix is the reference pixel position (x,y)
    Returns tuple (x,y) where x and y are arrays with pixel position for each RA,Dec
    """

    dtor = np.pi/180
    cd00 = -pixscale*dtor/3600
    cd01 = 0.0
    cd10 = 0.0
    cd11 = -cd00
    determ = cd00*cd11-cd01*cd10
    cdinv00 =  cd11/determ
    cdinv01 = -cd01/determ
    cdinv10 = -cd10/determ
    cdinv11 =  cd00/determ

    cos_crval1 = np.cos(dtor*dec_cen)
    sin_crval1 = np.sin(dtor*dec_cen)

    radif = (ra - ra_cen)*dtor
    w = np.where(radif > np.pi)
    radif[w] -= 2*np.pi
    w = np.where(radif < -np.pi)
    radif[w] += 2*np.pi

    decrad = dec*dtor
    cos_dec = np.cos(decrad)
    sin_dec = np.sin(decrad)
    cos_radif = np.cos(radif)
    sin_radif = np.sin(radif)
    h = sin_dec*sin_crval1 + cos_dec*cos_crval1*cos_radif
    xsi = cos_dec*sin_radif/h
    eta = (sin_dec*cos_crval1 - cos_dec*sin_crval1*cos_radif)/h
    xdif = cdinv00*xsi + cdinv01*eta
    ydif = cdinv10*xsi + cdinv11*eta
    return (xdif+crpix[0], ydif+crpix[1])


def xy2sky_tan(x, y, ra_cen, dec_cen, crpix=(0.0,0.0)):

    """Convert X,Y pixel position to RA,Dec sky position (degrees)

    x[n], y[n] are input arrays in pixels
    ra_cen[n], dec_cen[n] are image centers in degrees
    crpix is the reference pixel position (x,y)
    Returns tuple (ra,dec) where ra and dec are arrays with pixel position for each x,y
    """

    dtor = np.pi/180
    cd00 = -pixscale*dtor/3600
    cd01 = 0.0
    cd10 = 0.0
    cd11 = -cd00

    cos_crval1 = np.cos(dtor*dec_cen)
    sin_crval1 = np.sin(dtor*dec_cen)

    xdif = x - crpix[0]
    ydif = y - crpix[1]
    xsi = cd00*xdif + cd01*ydif
    eta = cd10*xdif + cd11*ydif
    beta = cos_crval1 - eta*sin_crval1
    ra = np.arctan2(xsi, beta) + dtor*ra_cen
    gamma = np.sqrt(xsi**2 + beta**2)
    dec = np.arctan2(eta*cos_crval1+sin_crval1, gamma)
    return (ra/dtor, dec/dtor)
