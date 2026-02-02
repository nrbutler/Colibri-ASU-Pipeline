from numpy import array,exp
from scipy.ndimage import gaussian_filter

def madau_extin(lam,z,smooth=False):
    """
    """
    # define the lyman series lines and edge strengths
    lam0=array([1215.67,1025.72,972.537,949.743,937.803,930.748,926.226,923.150, 920.963,919.352,918.129,917.181,916.429,915.824,915.329,914.919,914.576])

    AA=array([0.0036,0.0017,0.0011846,0.0009410,0.0007960,0.0006967,0.0006236,0.0005665,0.0005200,0.0004817,0.0004487,0.0004200,0.0003947,0.000372,0.0003520,0.0003334,0.00031644])

    tau = 0.*lam
    for i in range(len(lam0)):
        h=(lam<=lam0[i]*(1+z))*(lam>=lam0[i])
        tau[h] += AA[i]*(lam[h]/lam0[i])**3.46

    # photo-electric absorption below the lyman limit
    h = lam<912.*(1+z)
    xe = 1+z
    xc = lam[h]/912.; xc[xc<1]=1
    xc[xc>xe]=xe

    tau[h] += 0.25*xc**3*(xe**0.46-xc**0.46)
    tau[h] += 9.4*xc**1.5*(xe**0.18-xc**0.18)
    tau[h] -= 0.7*xc**3*(xc**(-1.32)-xe**(-1.32))
    tau[h] -= 0.023*(xe**1.68-xc**1.68)

    mdl = exp(-tau)
    if (smooth):
        # number of lines between lyman-a and lyman-b
        nlines = 0.5*(1+z)**3.46*( 1 - (lam0[1]/lam0[0])**3.46 )
        dlam = (lam0[0]-lam0[1])*(1+z)/(lam[1]-lam[0])/nlines
        if (dlam<1): dlam=1.
        mdl = gaussian_filter(mdl,dlam)

    return mdl
