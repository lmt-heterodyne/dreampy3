"""
A few utilities to manipulate corr_cal files
"""
import types
import numpy
import os
#from dreampy.dreamdb.rsr.models import CorrCal as CorrCalTable
from .exceptions import LMTRedshiftError
from dreampy3 import dreampyParams, corr_cal
from dreampy3.logging import logger
from dreampy3.redshift.netcdf.redshift_corrcal_netcdf_file import RedshiftCorrCalFile
import glob
logger.name = __name__

blank = -32767.    # like in Gildas Class, use for blanking some channels

def get_blank_indices_array(array):
    return numpy.where(array==blank, True, False)

def get_blank_indices(array):
    return numpy.where(array == blank)

def blank_to_nan(array):
    arr = array.copy()
    ind = numpy.where(arr == blank)
    arr[ind] = numpy.nan
    return arr

def blank_to_zero(array):
    arr = array.copy()
    ind = numpy.where(arr == blank)
    arr[ind] = 0.0
    return arr

def nanblankmean(array):
    """Returns the mean of an array after ignoring the 
    nan numbers and blank'ed numbers"""
    ind = numpy.logical_and(numpy.isfinite(array), array != blank)
    return array[ind].mean()

def nanmean(array):
    """Returns the mean of an array after ignoring the
    nan numbers"""
    return array[numpy.isfinite(array)].mean()

def _set_bad_lags(chassis, board, blags):
    """
    Sets the bad lags on a given board
    for a given chassis"""
    if chassis not in range(4):
        raise LMTRedshiftError("_set_bad_lags", "Chassis should be 0,1,2 or 3")
    if board not in range(6):
        raise LMTRedshiftError("_set_bad_lags", "Board should be 0,1,2,3,4 or 5")
    if type(blags) in (types.ListType, types.TupleType):
        blag_str = '/'.join(map(str, blags))
    elif type(blags) == types.IntType:
        blag_str = '%d' % blags
    elif type(blags) == types.StringType:
        blag_str = blags
    bl = []
    for bd in range(6):
        if bd != board:
            bl.append(dreampyParams['redshiftchassis']['bad_lags%d' % chassis][bd])
        else:
            bl.append(blag_str)
    dreampyParams['redshiftchassis']['bad_lags%d' % chassis] = bl
    #dreampyParams['redshiftchassis']['bad_lags%d' % chassis][board] = blag_str
    print(dreampyParams['redshiftchassis']['bad_lags%d' % chassis][board])
    
def set_bad_lags():
    """
    Interactive set_bad_lags for redshift boards
    """
    chassis = int(raw_input("Enter Chassis number (one of 0,1,2 or 3) > "))
    board = int(raw_input("Enter Board number (one of 0,1,2,3,4 or 5) > "))
    print("Current bad lags for chassis %d and board %d : %s" % (chassis, board, get_bad_lags(chassis, board)))
    blags_str = raw_input("Enter bad lags for Chassis %d, board %d in comma-separated form > " % (chassis, board))
    if blags_str:
        blags = map(int, blags_str.split(','))
    else:
        blags = ''
    _set_bad_lags(chassis, board, blags)
    logger.info("Bad Lags for chassis %d board %d set to %s" % (chassis, board, dreampyParams['redshiftchassis']['bad_lags%d' % chassis][board]))

    
def get_bad_lags(chassis, board):
    """Given a chassis and a board, returns
    badlags set in configuration rc params"""
    blag_str = dreampyParams['redshiftchassis']['bad_lags%d' % chassis][board]
    if blag_str == '':
        return []
    else:
        return list(map(int, blag_str.split('/')))

def apply_bad_lags(data, chassis, board):
    """Given a 256 point ACF, applies the bad lags for
    the given chassis and given board number"""
    dt = data.copy()
    badlags = get_bad_lags(chassis, board)
    if badlags:
        print("Found Bad lags %s in configuration for chassis %d board %d" % (badlags, chassis, board))
    for blag in badlags:
        dt[blag] = numpy.nan
    return dt

def print_bad_lags(chassis=None, board=None):
    """Lists the bad lags by chassis and board.
    If chassis and board are not specified, lists bad lags for
    all chassis and all boards"""
    if chassis is None:
        chassis = numpy.arange(4)
    elif type(chassis) in (types.ListType, types.TupleType):
        chassis = chassis
    elif isinstance(chassis, numpy.ndarray):
        chassis = chassis
    elif type(chassis) == types.IntType:
        chassis = [chassis]
    if board is None:
        board = numpy.arange(6)
    elif type(board) in (types.ListType, types.TupleType):
        board = board
    elif isinstance(board, numpy.ndarray):
        board = board
    elif type(board) == types.IntType:
        board = [board]
    for chas in chassis:
        for bd in board:
            badlags = get_bad_lags(chas, bd)
            print("Chassis: %d, board: %d, badlags=%s" % (chas, bd, badlags))


def compress_corr_cal(corr_cal):
    for i in range(256):
        if i == 0:
            cc = corr_cal[1+4*i:4*i+5].mean(axis=0)
        else:
            cc = numpy.vstack((cc, corr_cal[1+4*i:4*i+5].mean(axis=0)))
    return cc


def check_array(data):
    """checks if data is a numpy array with at least one non-zero
    value"""
    if type(data) not in (numpy.ndarray, numpy.ma.core.MaskedArray) or ((type(data) == numpy.ndarray) and not data.any()):
        return False
    else:
        return True

def rebin(specin, outsize=256):
    """Very simple boxcar average of input spectra"""
    f = len(specin)/outsize
    return numpy.array([specin[f*i:f*i+f].mean() for i in range(outsize)])

def makespectrum(data, ratio=None, g_norm=None, corr_cal=None,
                 #lowfreq=1.25, highfreq=7.9, badlags=None, numchan=256):
                 lowfreq=1.3, highfreq=7.9, badlags=None, numchan=256):
    """returns g_normed and corr_cal smoothed spectrum. Also
    does the elimination of end channels and bad channels"""

    freq = (0.0078125*2.5) + numpy.arange(256)*0.03125
    if type(data) not in (numpy.ndarray, numpy.ma.core.MaskedArray):
        print("Needs to be array data type")
        return

    if badlags is not None:
        for lag in badlags:
            data[lag] = numpy.nan
    if check_array(g_norm):
        normacf = data*g_norm
        normacf = normacf - nanmean(normacf[65:])
        normacf  = normacf/g_norm
    else:
        normacf = data
    if check_array(ratio):
        normacf = normacf*ratio

    if not check_array(corr_cal):
        print("Corr cal not present, return acf")
        return normacf
    
    spec = numpy.dot(corr_cal, numpy.nan_to_num(normacf))
    #x = 8.0*numpy.arange(len(spec))/(len(spec)-1.)
    #xx = 8.0*numpy.arange(numchan)/(numchan-1.)
    #intinst = Interpolate(x,spec)
    #return intinst(xx)
    if corr_cal.shape[0] > 256:
        spec = rebin(spec[1:])

    #spec = numpy.flipud(spec)
    ind = numpy.logical_or(freq<lowfreq, freq>highfreq)
    #spec[:lochan] = numpy.nan
    #spec[hichan:] = numpy.nan
    #spec = numpy.flipud(spec)
    spec[ind] = numpy.nan
    return spec


def get_corr_cal_matrix(header, slotno=None):
    """Given a netCDF based header value, uses the UTDate
    in the header to determine from the corr_cal database which
    corr_cal file to read, and reads it and returns the corr_cal
    matrix.
    If slotno is given it can either be a single integer or a list of
    slotnumbers. Returned value is a dictionary of corr_cal values
    with slotnumbers as keys. If slotno is None, all slots in the header
    are returned.
    """
    #from dreampy.redshift.netcdf.redshift_corrcal_netcdf_file import RedshiftCorrCalFile
    if slotno is None:
        slotno = header.SlotNumber
    elif type(slotno) == types.IntType:
        slotno = [slotno]
    elif type(slotno) in (numpy.ndarray, types.ListType, type.TupleType):
        slotno = slotno
    else:
        #something wrong
        slotno = header.SlotNumber

    chassis = int(header.ChassisNumber)
    ccaldic = {}
    corr_cal_dir = os.environ.get('CORR_CAL_DIR', '/raw/rsr/cal')
    for slot in slotno:
        try:
            ccal = CorrCalTable.objects.filter(chassis=chassis,
                                               slotno=slot,
                                               obsdate__lt=header.utdate()).latest('obsdate')
            ccal_filename = ccal.filename
        except:
            #cannot reach the internet, or the database
            glb = glob.glob(os.path.join(corr_cal_dir, 'corr_cal_%s_%s.nc' % (int(header.CorrCalID), chassis)))
            if glb:
                ccal_filename = os.path.basename(glb[0])
            else:
                raise LMTRedshiftError("get_corr_cal_matrix", "Some error in getting corr_cal table information")
        fname = os.path.join(corr_cal_dir, os.path.basename(ccal_filename))
        logger.info("Corr Cal file for Slot no %d is %s" % (slot, fname))
        if fname in corr_cal and slot in corr_cal[fname]:
            ccaldic[slot] = corr_cal[fname][slot]
            logger.info("Re-reading Corr Cal from memory for Slot no %d for %s" % (slot, fname))
        else:
            if not fname in corr_cal:
                corr_cal[fname] = {}
            if not slot in corr_cal[fname]:
                if os.path.exists(fname):
                    logger.info("Reading first corr cal for slot no %d for %s" % (slot, fname))
                    ccalnc = RedshiftCorrCalFile(fname)
                    idx = int(numpy.where(ccalnc.hdu.header.slotno == slot)[0])
                    ccaldic[slot] = ccalnc.hdu.data.corr_cal[idx, :, :]
                    corr_cal[fname][slot] = ccalnc.hdu.data.corr_cal[idx, :, :]
                else:
                    raise LMTRedshiftError("get_corr_cal_matrix", "Corr cal file %s does not exist" % fname)
    return ccaldic
