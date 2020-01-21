import numpy

def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = numpy.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx
