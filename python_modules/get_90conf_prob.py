from numpy import hstack,newaxis,where

def get_90conf_prob(xx,probxx,plim=0.9):
    """
    returns most compact 90% conf interval (for default plim=0.9)
    """
    s=xx.argsort(); x=xx[s]

    dx = x[1:]-x[:-1]
    dx = hstack( (dx[0],dx) )
    cprob = (probxx[s]*dx).cumsum()
    cprob = hstack( (0,cprob) )

    cprob_m = cprob[-1]

    cprob1 = cprob[:,newaxis] - cprob[newaxis,:]

    #i>j
    i,j = where( cprob1 >= plim*cprob_m )

    delta = x[i-1] - x[j]
    i0 = delta.argmin()

    return x[j[i0]],x[i[i0]-1]
