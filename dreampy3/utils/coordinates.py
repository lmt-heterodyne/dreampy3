#import math
#import numpy
from astropy import units as u
from astropy.coordinates import Angle
import numpy


def sixty(valrad, deg=False, ra=False):
    """returns a pretty string version of the RA or Dec value. Input
    is in radians if deg is set to False. If ra is needed, set ra=True.
    >>> sixty(numpy.radians(63.3233))
    '63:19:23'
    >>> sixty(numpy.radians(4.344*360/24.), ra=True)
    '04:20:38'
    """
    if deg is False:
        unit = 'rad'
    else:
        unit = 'deg'
    ang = Angle("%s%s" % (valrad, unit))
    valdeg = numpy.degrees(valrad)
    if ra:
        retunit = u.hour
    else:
        retunit = u.degree
    return ang.to_string(unit=retunit, sep=':')
    # if ra:
    #     valdeg *= 24./360.0
    # degs = int(valdeg)
    # absdegs = abs(degs)
    # valmins = (abs(valdeg)-absdegs)*60.0
    # mins = int(valmins)
    # secs = int((valmins-mins)*60.0)
    # return "%s%02d:%02d:%02d" % (neg, degs, mins, secs)
