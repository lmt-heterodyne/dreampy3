import numpy
from dreampy.utils.mpfit import mpfit

def gaussfunc2gauss2d(p, x, y):
    return (p[0]*numpy.exp(-(x-p[1])**2/(2*p[3]**2)) * \
            numpy.exp(-(y-p[2])**2/(2*p[4]**2))) + \
            (p[5]*numpy.exp(-(x-p[6])**2/(2*p[8]**2)) * \
             numpy.exp(-(y-p[7])**2/(2*p[9]**2)))

def gaussfunc2gauss2d_samesigma(p, x, y):
    return (p[0]*numpy.exp(-(x-p[1])**2/(2*p[3]**2)) * \
            numpy.exp(-(y-p[2])**2/(2*p[4]**2))) + \
            (p[5]*numpy.exp(-(x-p[6])**2/(2*p[3]**2)) * \
             numpy.exp(-(y-p[7])**2/(2*p[4]**2)))

def twogauss_deviation_func(p, fjac=None, Z=None,
                            X=None, Y=None):
    """
    Z is actual Z data
    X and Y are two dimensional arcsecond offsets
    """
    
    zmodel = gaussfunc2gauss2d(p, X, Y)
    model_deviate = (Z - zmodel)
    #print model_deviate.mean()
    status = 0
    return [status, numpy.sqrt(model_deviate**2)
            ]

def twogauss_samesigma_deviation_func(p, fjac=None, Z=None,
                                      X=None, Y=None):
    """
    Z is actual Z data
    X and Y are two dimensional arcsecond offsets
    """
    
    zmodel = gaussfunc2gauss2d_samesigma(p, X, Y)
    model_deviate = (Z - zmodel)
    #print model_deviate.mean()
    status = 0
    return [status, numpy.sqrt(model_deviate**2)
            ]

def two_gaussian_fit(X, Y, Z, pinit=numpy.zeros(10, dtype='float'),
                     maxiter=1000, quiet=False):
    fa = {'Z' : Z.flatten(),
          'X' : X.flatten(),
          'Y' : Y.flatten()
          }
    m = mpfit(twogauss_deviation_func, pinit,
              functkw = fa, maxiter=maxiter, quiet=quiet)
    
    return m

def two_gaussian_samesigma_fit(X, Y, Z, pinit=numpy.zeros(8, dtype='float'),
                               maxiter=1000, quiet=False):
    fa = {'Z' : Z.flatten(),
          'X' : X.flatten(),
          'Y' : Y.flatten()
          }
    m = mpfit(twogauss_samesigma_deviation_func, pinit,
              functkw = fa, maxiter=maxiter, quiet=quiet)
    
    return m
