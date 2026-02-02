from numpy import abs,sqrt,std,median,atleast_1d

def robust_mean(x,dx):
    """
    ala stetson
    """
    alpha = 2.
    beta = 3.
    x0 = median(x)
    nn = len(atleast_1d(x))

    if (nn<2):
        return x0,0.
    else:

        for i in range(20):
            resid = x - x0
            weight = 1./dx**2
            resid_err = abs(resid)/dx
            weight1 = weight/(1. + (resid_err/alpha)**beta)

            x0 = (x*weight1).sum() / weight1.sum()

        dx0 = sqrt( ( (x-x0)**2*weight1 ).sum()/weight1.sum()/(nn-1.) )

        return x0,dx0
