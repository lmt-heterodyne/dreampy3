"""
This module contains some general classes for
curve fitting that could be used by any of
dreampy's modules
Uses scipy ODR (orthogonal distance regression)
for its fitting
"""

from scipy import odr
#from scipy import optimize
import numpy
import math
from dreampy3.utils import DreampyGeneralError
from dreampy3.redshift.utils.spectrum_utils import nanmean

class OutputLeastSquares(object):
    """A pseudo output class for the
    scipy.optimize.leastsq module to look somewhat
    similar to the objects in odr.Output class
    Not all objects in ODR output are reproduced!"""
    def __init__(self, leastsq_out):
        """leastsq_out is the output of leastsq with
        full_output=1
        """
        self.leastsq_out = leastsq_out
        self.beta = self.leastsq_out[0]
        self.cov_x = self.leastsq_out[1]
        self.infodict = self.leastsq_out[2]
        self.stopreason = self.leastsq_out[3]
        self.info = self.leastsq_out[4]
        self.x = None
        self.y = None
        self.res_var = None
        self.beta_err = None
    
    def update1d(self, x=None, y=None, yfit=None):
        """Given x, y, yfit, calculates and updates
        residual variance, beta_error, etc."""
        self.x = x
        self.y = yfit
        self.res_var = ((yfit - y)**2.).sum()/(x.size-self.beta.size)
        self.sd_beta = numpy.sqrt((self.cov_x*self.res_var).diagonal())

    def update2d(self, x=None, y=None, z=None, zfit=None):
        """Given x, y, z, zfit, calculates and updates
        residual variance, beta_error, etc."""
        self.x = x
        self.y = zfit
        self.res_var = ((zfit - z)**2.).sum()/((x.size+y.size)-self.beta.size)
        self.sd_beta = numpy.sqrt((self.cov_x*self.res_var).diagonal())

    def pprint(self):
        mystr = "Beta: %s\n" % self.beta
        mystr += "Beta Std Error: %s\n" % self.sd_beta
        mystr += "Beta Covariance: %s\n" % self.cov_x
        mystr += "Residual Variance : %s\n" % self.res_var
        mystr += "Reasons for Halting :\n%s" % self.stopreason
        print(mystr)
        #return mystr

class DataFit(object):
    def __init__(self, use_odr=True):
        self.use_odr = use_odr
        self.explicit = False
        
    def _set_Data(self, x, z):
        self.Data = odr.Data(x, z)

    def _set_Model(self, func, fjacb=None, fjacd=None,
                   estimate=None, extra_args=None):
        if fjacb is not None or fjacd is not None:
            self.explicit = True
        else:
            self.explicit = False
        self.Model = odr.Model(func, fjacb=fjacb,
                               fjacd=fjacd, estimate=estimate,
                               extra_args=extra_args)

    def _run(self, beta0=None):
        self.ODR = odr.ODR(data=self.Data, model=self.Model,
                           beta0=beta0)
        if self.explicit:
            self.ODR.set_job(deriv=2)  #give explicit (user-supplied) derivatives, Jacobians
        self.out = self.ODR.run()
        #print self.out.stopreason
        return self.out

class Gauss1DFit(DataFit):
    """
    Fits 1-dimensional gaussian to two input vectors
    x and z.

    >>> import numpy
    >>> x = numpy.linspace(-10, 10, 1000)
    >>> y = 1*numpy.exp(-0.5*((x-0)/3.)**2)
    >>> gfit = Gauss1DFit(x, y)
    >>> out = gfit._run()
    >>> print out.beta
    [  1.00000000e+00  -8.63063585e-17   3.00000000e+00]
    """
    def __init__(self, x, z, peak=None,
                 mu=None, sigma=None):
        DataFit.__init__(self)
        self.peak, self.mu, self.sigma = peak, mu, sigma
        self._set_Data(x, z)
        self._set_Model(self._gaussfunc1d, fjacb=self._gaussfunc1d_fjacb,
                        fjacd=self._gaussfunc1d_fjacd,
                        estimate=self._estimate)
                        
    def _gaussfunc1d(self, p, x):
        return p[0]*numpy.exp(-0.5*((x-p[1])/p[2])**2)


    def _gaussfunc1d_fjacb(self, p, x):
        """Jacobian wrt to the p parameters"""
        #return numpy.vstack([self._gaussfunc1d(p,x)/p[0],
        #                     self._gaussfunc1d(p,x)*(x-p[1])/p[2]**2.,
        #                     self._gaussfunc1d(p,x)*((x-p[1])/p[2])**2/p[2]])
        return numpy.vstack([self._gaussfunc1d(p,x)/p[0],
                             self._gaussfunc1d(p,x)*(x-p[1])/p[2]**2,
                             self._gaussfunc1d(p,x)*(((x-p[1])/p[2])**2)/p[2]])
    
    def _gaussfunc1d_fjacd(self, p, x):
        """Jacobian wrt to x"""
        return (-(x-p[1])/p[2]**2)*self._gaussfunc1d(p,x)

    def _estimate(self, data):
        x = data.x
        z = data.y
        # guess some fit parameters
        if self.peak is None:
            self.peak = numpy.max(z)
        if self.mu is None:
            self.mu = numpy.sum(z*x)/numpy.sum(z)
        if self.sigma is None:
            self.sigma = numpy.sqrt(numpy.abs(numpy.sum(z*(x-self.mu)**2.)/numpy.sum(z)))
        if self.sigma < 0.0:
            print('No signal to estimate 2nd moment')
            raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
        #print peak, mu, math.sqrt(sigma)
        return numpy.array([self.peak, self.mu, math.sqrt(self.sigma)])

class Gauss1DFitBase(DataFit):
    """
    Fits 1-dimensional gaussian with a baseline to two input vectors
    x and z.

    >>> import numpy
    >>> x = numpy.linspace(-10, 10, 1000)
    >>> y = 1*numpy.exp(-0.5*((x-0)/3.)**2)
    >>> gfit = Gauss1DFitBase(x, y)
    >>> out = gfit._run()
    >>> print out.beta
    [  1.00000000e+00  -8.63063585e-17   3.00000000e+00]
    """
    def __init__(self, x, z, base=None, peak=None,
                 mu=None, sigma=None):
        DataFit.__init__(self)
        self.base, self.peak, self.mu, self.sigma = base, peak, mu, sigma
        self._set_Data(x, z)
        self._set_Model(self._gaussfunc1d, fjacb=self._gaussfunc1d_fjacb,
                        fjacd=self._gaussfunc1d_fjacd,
                        estimate=self._estimate)
                        
    def _gaussfunc1d(self, p, x):
        return p[0] + p[1]*numpy.exp(-0.5*((x-p[2])/p[3])**2)

    def _gaussfunc1d_fjacb(self, p, x):
        """Jacobian wrt to the p parameters"""
        #return numpy.vstack([self._gaussfunc1d(p,x)/p[0],
        #                     self._gaussfunc1d(p,x)*(x-p[1])/p[2]**2.,
        #                     self._gaussfunc1d(p,x)*((x-p[1])/p[2])**2/p[2]])
        return numpy.vstack([numpy.ones(x.size),
                             self._gaussfunc1d(p,x)/p[1],
                             self._gaussfunc1d(p,x)*(x-p[2])/p[3]**2,
                             self._gaussfunc1d(p,x)*(((x-p[2])/p[3])**2)/p[3]])
    
    def _gaussfunc1d_fjacd(self, p, x):
        """Jacobian wrt to x"""
        return (-(x-p[2])/p[3]**2)*self._gaussfunc1d(p,x)

    def _estimate(self, data):
        x = data.x
        z = data.y
        if self.base is None:
            self.base = z.mean()
        # guess some fit parameters
        if self.peak is None:
            self.peak = numpy.max(z) - self.base
        if self.mu is None:
            self.mu = numpy.sum(z*x)/numpy.sum(z)
        if self.sigma is None:
            self.sigma = numpy.sum(z*(x-self.mu)**2.)/numpy.sum(z)
        if self.sigma < 0.0:
            print('No signal to estimate 2nd moment')
            raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
        #print peak, mu, math.sqrt(sigma)
        return numpy.array([self.base, self.peak, self.mu, math.sqrt(self.sigma)])
    
    
class Gauss2DFit(DataFit):
    def __init__(self, x, y, z):
        DataFit.__init__(self)
        self.x, self.y, self.z = x, y, z
        self._set_Data(x, z)
        self._set_Model(self._gaussfunc2d, 
                        estimate=self._estimate,
                        extra_args=(y,))

    def _gaussfunc2d(self, p, x, y):
        return p[0]*numpy.exp(-0.5*((x-p[1])/p[3])**2) * \
               numpy.exp(-0.5*((y-p[2])/p[4])**2)

    def _estimate(self, data):
        # guess some fit parameters
        peak = numpy.nanmax(self.z)
        mux = (self.z*self.x).sum()/self.z.sum()
        muy = numpy.sum(self.z*self.y)/numpy.sum(self.z)
        sigmax = (self.z*(self.x-mux)**2.).sum()/numpy.abs(self.z.sum())
        sigmay = (self.z*(self.y-muy)**2.).sum()/numpy.abs(self.z.sum())
        if sigmax < 0.0:
            print('No signal to estimate 2nd x moment')
            raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
        if sigmay < 0.0:
            print('No signal to estimate 2nd y moment')
            raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
        p0 = numpy.array([peak, mux, muy, math.sqrt(sigmax), math.sqrt(sigmay)])
        print(p0)
        return p0


class Gauss2DFitAvg(DataFit):
    def __init__(self, x, y, z, explicit=False):
        DataFit.__init__(self)
        self.x, self.y, self.z = x, y, z
        self._set_Data(x, z)
        self.explicit = explicit
        if self.explicit:
            self._set_Model(self._gaussfunc2d, 
                            estimate=None,
                            extra_args=(y,))
        else:
            self._set_Model(self._gaussfunc2d, 
                            estimate=self._estimate,
                            extra_args=(y,))

    def _gaussfunc2d(self, p, x, y):
        return p[0] + p[1]*numpy.exp(-0.5*((x-p[2])/p[4])**2) * \
               numpy.exp(-0.5*((y-p[3])/p[5])**2)

    def _estimate(self, data):
        # guess some fit parameters
        avg = nanmean(self.z)
        peak = numpy.nanmax(self.z)-avg
        mux = (self.z*self.x).sum()/self.z.sum()
        muy = numpy.sum(self.z*self.y)/numpy.sum(self.z)
        # let's not mess with sigma estimation
        #sigmax = (self.z*(self.x-mux)**2.).sum()/numpy.abs(self.z.sum())
        #sigmay = (self.z*(self.y-muy)**2.).sum()/numpy.abs(self.z.sum())
        # if sigmax < 0.0:
        #     print 'No signal to estimate 2nd x moment'
        #     raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
        # if sigmay < 0.0:
        #     print 'No signal to estimate 2nd y moment'
        #     raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
        sigmax, sigmay = 100, 100
        p0 = numpy.array([avg, peak, mux, muy, math.sqrt(sigmax), math.sqrt(sigmay)])
        print(p0)
        return p0

    # def _estimate(self, data):
    #     # guess some fit parameters
    #     peak = numpy.nanmax(self.z)
    #     mux = (self.z[numpy.isfinite(self.z)]*self.x[numpy.isfinite(self.x)]).sum()/self.z[numpy.isfinite(self.z)].sum()
    #     muy = numpy.sum(self.z[numpy.isfinite(self.z)]*self.y[numpy.isfinite(self.y)])/numpy.sum(self.z[numpy.isfinite(self.z)])
    #     sigmax = (self.z[numpy.isfinite(self.z)]*(self.x[numpy.isfinite(self.x)]-mux)**2.).sum()/numpy.abs(self.z[numpy.isfinite(self.z)].sum())
    #     sigmay = (self.z[numpy.isfinite(self.z)]*(self.y[numpy.isfinite(self.y)]-muy)**2.).sum()/numpy.abs(self.z[numpy.isfinite(self.z)].sum())
    #     if sigmax < 0.0:
    #         print 'No signal to estimate 2nd x moment'
    #         raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
    #     if sigmay < 0.0:
    #         print 'No signal to estimate 2nd y moment'
    #         raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
    #     p0 = numpy.array([peak, mux, muy, math.sqrt(sigmax), math.sqrt(sigmay)])
    #     print p0
    #     return p0


class TwoGauss2DFit(DataFit):
    """Fits two gaussians in the same field"""
    def __init__(self, x, y, z, pinit=None):
        DataFit.__init__(self)
        self.x, self.y, self.z = x, y, z
        self._set_Data(x, z)
        self.pinit = pinit
        self._set_Model(self._gaussfunc2gauss2d, 
                        estimate=self._estimate,
                        extra_args=(y,))

    def _gaussfunc2gauss2d(self, p, x, y):
        return (p[0]*numpy.exp(-0.5*((x-p[1])/p[3])**2) * \
               numpy.exp(-0.5*((y-p[2])/p[4])**2)) + \
               (p[5]*numpy.exp(-0.5*((x-p[6])/p[8])**2) * \
               numpy.exp(-0.5*((y-p[7])/p[9])**2)) 

    def _estimate(self, data):
        # guess some fit parameters
        if self.pinit is not None:
            print(self.pinit)
            return self.pinit
        peak1 = self.z.max()
        mux1 = self.x[self.z == self.z.max()]
        muy1 = self.y[self.z == self.z.max()]
        #mux = (self.z*self.x).sum()/self.z.sum()
        #muy = numpy.sum(self.z*self.y)/numpy.sum(self.z)
        #sigmax = numpy.abs((self.z*(self.x-mux)**2.).sum()/numpy.abs(self.z.sum()))
        #sigmay = numpy.abs((self.z*(self.y-muy)**2.).sum()/numpy.abs(self.z.sum()))
        #if sigmax < 0.0:
        #    print 'No signal to estimate 2nd x moment'
        #    raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
        #if sigmay < 0.0:
        #    print 'No signal to estimate 2nd y moment'
        #    raise DreampyGeneralError("No Signal", "No signal to estimate 2nd moment")
        sigmax1, sigmay1 = 20., 20.0
        peak2 = self.z.min()
        mux2 = self.x[self.z == self.z.min()]
        muy2 = self.y[self.z == self.z.min()]        
        p0 = numpy.array([peak1, mux1, muy1, sigmax1, sigmay1,
                          peak2, mux2, muy2, sigmax1, sigmay1])
        print(p0)
        return p0

class ParabolicFit(DataFit):
    """
    Fits 1-dimensional parabola to two input vectors
    x and y.
    >>> vertical parabolic formula
    >>> (x-h)**2 = 4p(y-k)
    >>> import numpy
    >>> x = numpy.linspace(-10, 10, 1000)
    >>> y = x**2/4*1.5 - 2
    >>> gfit = Gauss1DFit(x, y)
    >>> out = gfit._run()
    >>> print out.beta
    [  1.00000000e+00  -8.63063585e-17   3.00000000e+00]
    """
    def __init__(self, x, z, p=None,
                 h=None, k=None, explicit=True):
        DataFit.__init__(self)
        self.p, self.h, self.k = p, h, k
        self._set_Data(x, z)
        self.explicit = explicit
        if self.explicit:
            self._set_Model(self._parabolic1d, fjacb=None,
                            fjacd=None,
                            estimate=self._estimate)
        else:
            self._set_Model(self._parabolic1d, fjacb=self._parabolic1d_fjacb,
                            fjacd=self._parabolic1d_fjacd,
                            estimate=self._estimate)
                        
    def _parabolic1d(self, p, x):
        #return p[0]*numpy.exp(-0.5*((x-p[1])/p[2])**2)
        return x**2/(4*p[0]) - 2*x*p[1]/(4*p[0]) + (p[1]**2 + 4*p[0]*p[2])/(4*p[0])

    def _parabolic1d_fjacb(self, p, x):
        """Jacobian wrt to the p parameters"""
        #return numpy.vstack([self._gaussfunc1d(p,x)/p[0],
        #                     self._gaussfunc1d(p,x)*(x-p[1])/p[2]**2.,
        #                     self._gaussfunc1d(p,x)*((x-p[1])/p[2])**2/p[2]])
        return numpy.vstack([-x**2/(4*p[0]**2) - (x*p[1])/(2*p[0]**2) - p[1]**2/(4*p[0]**2),
                              -x/(2*p[0]) + (2*p[1])/(4*p[0]),
                              numpy.ones(x.size)])

    
    def _parabolic1d_fjacd(self, p, x):
        """Jacobian wrt to x"""
        return (-x/(2*p[0])) - (p[1]/(2*p[0]))

    def _estimate(self, data):
        x = data.x
        z = data.y
        # guess some fit parameters
        if self.p is None:
            self.p = -1.0
        if self.h is None:
            self.h = x.mean()
        if self.k is None:
            self.k = z.max()
        return numpy.array([self.p, self.h, self.k])


if __name__ == '__main__':
    import doctest
    doctest.testmod()
