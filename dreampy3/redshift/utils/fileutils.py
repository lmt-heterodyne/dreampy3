"""
A few simple file utils inspired by 
Pete's simplification functions
"""
from dateutil.parser import parse
from .exceptions import LMTRedshiftError
import glob
import os

def get_data_root():
    if "DATA_LMT" in os.environ.keys():
        data_root =os.environ["DATA_LMT"].strip()
        if data_root[-1] !=os.sep:
            data_root +=os.sep
    else:
        data_root="/data_lmt/"
    return data_root

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
    filename = get_data_root()+ 'RedshiftChassis%s/RedshiftChassis%s_%s_%06d_%02d_0001.nc' % (chassis, chassis, dt.strftime('%Y-%m-%d'), scan, subobsnum)
    return filename

def make_generic_filename_old(scan, chassis, subobsnum=1):
    if subobsnum is None:
        glb = glob.glob(get_data_root()+'RedshiftChassis%s/RedshiftChassis%s_*_%06d_*_0001.nc' % (chassis, chassis, scan))
    else:
        glb = glob.glob(get_data_root()+'RedshiftChassis%s/RedshiftChassis%s_*_%06d_%02d_0001.nc' % (chassis, chassis, scan, subobsnum))
    if glb:
        return glb[0]

def make_generic_filename(obsnum, chassis):
    glb = glob.glob(os.path.join(get_data_root(), 'RedshiftChassis%s/RedshiftChassis%s*%06d*.nc' % (chassis, chassis, obsnum)))
    if glb:
        return glb[0]
    
    
