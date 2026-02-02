from numpy import ones

def linfit(x,y,dy=[],slope_prior=0.0,slope_prior_err=0.0):
    """
    m = a+b*x
    minimize chi^2 = Sum (y-m)^2/dy^2
    """
    lx=len(x)

    if (lx==0): return (0,0,[[0,0],[0,0]])
    if (lx==1): return (y.mean(),0,[[0,0],[0,0]])

    if (len(dy)==0):
        dy = ones(lx,dtype='float32')

    wt = 1./dy**2
    ss = wt.sum()
    sx = (wt * x).sum()
    sy = (wt * y).sum()
    t =  (x - sx/ss) / dy
    b = (t * y / dy).sum()

    st2 = (t*t).sum()

    # parameter estimates
    if (st2>0):
        b /= st2
        if (slope_prior_err>0):
            b = ( b + slope_prior/(st2*slope_prior_err**2) ) / ( 1. + 1./(st2*slope_prior_err**2) )
    else: b=0.0
    a = (sy - sx * b) / ss

    # error estimates
    svara,svarb,covar=1./ss,0.,0.
    if (st2>0):
        svarb = 1. / st2
        covar = -sx/(ss*st2)
        svara = (1. + sx * sx / (ss * st2)) / ss

    covar = [[svara, covar], [covar, svarb]]

    return (a,b,covar)
