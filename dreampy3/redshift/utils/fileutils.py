"""
A few simple file utils inspired by 
Pete's simplification functions
"""
from dateutil.parser import parse
from .exceptions import LMTRedshiftError
import glob

def make_filename(date, scan, chassis, subobsnum=1):
    """
    makes a filename from date, scan number, and chassis number
    date is in format 'yyyy-mm-dd' or even plain English like
    'Apr 22, 2013'
    scan is scan number (integer)
    chassis is chassis number (one of 0, 1, 2 or 3)
    """
    try:
        dt = parse(date)
    except ValueError:
        raise LMTRedshiftError("Arg Error", "Cannot parse date: %s" % date)
    if chassis not in (0, 1, 2, 3):
        raise LMTRedshiftError("Arg Error", "Chassis should be one of 0, 1, 2 or 3")
    filename = '/data_lmt/RedshiftChassis%s/RedshiftChassis%s_%s_%06d_%02d_0001.nc' % (chassis, chassis, dt.strftime('%Y-%m-%d'), scan, subobsnum)
    return filename

def make_generic_filename(scan, chassis, subobsnum=1):
    if subobsnum is None:
        glb = glob.glob('/data_lmt/RedshiftChassis%s/RedshiftChassis%s_*_%06d_*_0001.nc' % (chassis, chassis, scan))
    else:
        glb = glob.glob('/data_lmt/RedshiftChassis%s/RedshiftChassis%s_*_%06d_%02d_0001.nc' % (chassis, chassis, scan, subobsnum))
    if glb:
        return glb[0]
    
