"""
Simple utilities to do with smoothing of
data
"""
import numpy

def rebin(specin, smooth=2):
    """Very simple boxcar average of input 1-d data"""
    #print len(specin), smooth
    if smooth == 1:
        return specin
    specin = specin.copy()
    outsize = len(specin)/smooth
    f = smooth
    return numpy.array([specin[f*i:f*i+f].mean() for i in range(outsize)],
                       dtype=specin.dtype)

