"""Contains utilities to convert redshift ascii files to
Redshift NetCDF files"""

from dreampy3.lmtnetcdf import LMTNetCDFFile
from dreampy3.redshift.utils.redshift_netcdf_variables import dimensions, header_variables_dict, data_variables_dict
from dreampy3.redshift.utils import LMTRedshiftError
from dreampy3.redshift.utils.redshift_ascii import ReadAccumFiletoDic
#from dreampy.dreamdb.rsr.models import CorrCal as CorrCalTable
#from dreampy3.utils import OrderedDict
from collections import OrderedDict
import numpy
import os
import glob
from dateutil import parser
from netCDF4 import num2date, date2num
import time
import datetime
import shutil

def create_dimensions(nc, dims, numcards=6):
    dims['numcards'] = numcards
    for k, v in dims.items():
        nc.createDimension(k, v)

def create_variables(nc, hvardict, hname='RedshiftChassis',
                     fstr="Header"):
    vdic = {}
    for name, dic in hvardict.items():
        varname = '%s.%s.%s' % (fstr, hname, name)
        #print "Creating variable %s" % varname
        vdic[varname] = nc.createVariable(varname, dic['type'],
                                          dic['dim'])
        for k,v in dic['attributes'].items():
            vdic[varname].setncattr(k, v)
    return vdic


def fill_header_variables(nc, h, chassis=0,
                          hname='RedshiftChassis',
                          numcards=6):
    for name in ('blank', 'board', 'boardid', 'corr_cal_id',
                 'gain', 'main_acf_count', 'polar_demod',
                 'polar_out', 'swta', 'swtb', 'refacc', 'slotno',
                 'Hardware_Version', 'Software_Version'):
        varname = 'Header.%s.%s' % (hname, name)
        dtyp = nc.variables[varname].dtype
        typ = dtyp.type
        val = numpy.zeros((numcards,), dtype=dtyp)
        for i, board in enumerate(h[chassis].keys()):
            #print name, h[chassis][board][name]
            val[i] = typ(h[chassis][board][name])
        nc.variables[varname][:] = val
    #now for the string variables
    for name in ('Errormes', 'Probemes', 'date'):
        varname = 'Header.%s.%s' % (hname, name)
        txt = numpy.fromstring(h[chassis][h[chassis].keys()[0]][name],
                               dtype=numpy.dtype('S1'))
        nc.variables[varname][:len(txt)] = txt
    #single scalar values
    for name in ("chassis", "error", "target_size", "storage_mode",
                 "integ_status"):
        varname = 'Header.%s.%s' % (hname, name)
        nc.variables[varname][:] = int(h[chassis][h[chassis].keys()[0]][name])

def fill_data_variables(nc, d, chassis=0,
                        hname='RedshiftChassis',
                        numcards=6):
    #varname = 'Data.%s.AccData' % hname
    for i, board in enumerate(d[chassis].keys()):
        val = d[chassis][board]
        nc.variables['Data.%s.AccData' % hname][i,:] = val[:256]
        nc.variables['Data.%s.RefData' % hname][i,:] = val[256:]
    

def convert_ascii_netcdf(asciifile, ncfile, chassis=0):
    """asciifile is the full path name for ascii file
    ncfile is new netCDF file to be created"""
    h,c,d = ReadAccumFiletoDic(asciifile)
    if h.has_key(chassis):
        numcards = len(h[chassis])
        nc = LMTNetCDFFile(ncfile, 'w')
        create_dimensions(nc, dimensions, numcards=numcards)
        hdic = create_variables(nc, header_variables_dict,
                                hname='RedshiftChassis_%d_' % chassis)
        ddic = create_variables(nc, data_variables_dict, fstr='Data',
                                hname='RedshiftChassis_%d_' % chassis)
        fill_header_variables(nc, h, chassis=0, numcards=numcards,
                              hname='RedshiftChassis_%d_' % chassis)
        fill_data_variables(nc, d, chassis=0, numcards=numcards,
                            hname='RedshiftChassis_%d_' % chassis)
        #nc._make_hdu()
        nc.close()
        #return nc

def convert_asciilist_netcdf(asciilist, ncfile, chassis=0):
    """asciilist is a list of ascii file
    ncfile is new netCDF file to be created"""
    nc = LMTNetCDFFile(ncfile, 'w')
    for asciifile in asciilist:
        name, ext = os.path.splitext(os.path.basename(asciifile))
        h,c,d = ReadAccumFiletoDic(asciifile)
        if h.has_key(chassis):
            group = nc.createGroup(name) 
            numcards = len(h[chassis])
            create_dimensions(group, dimensions, numcards=numcards)
            hdic = create_variables(group, header_variables_dict,
                                    hname='RedshiftChassis_%d_' % chassis)
            ddic = create_variables(group, data_variables_dict, fstr='Data',
                                    hname='RedshiftChassis_%d_' % chassis)
            fill_header_variables(group, h, chassis=0, numcards=numcards,
                                  hname='RedshiftChassis_%d_' % chassis)
            fill_data_variables(group, d, chassis=0, numcards=numcards,
                                hname='RedshiftChassis_%d_' % chassis)
    nc.close()

def _get_ccal_info(filename):
    fp = open(filename)
    lines = fp.readlines()
    fp.close()
    args = lines[0].strip().split(',')
    svd = int(args[0].strip('#'))
    zero_lag = int(args[1])
    return svd, zero_lag

def _get_bid(filename, chassis=0):
    h,c,d = ReadAccumFiletoDic(filename)
    lis = [(dd['slotno'], dd['boardid']) for dd in h[chassis].values()]
    slotdic = {}
    for sno, bid in lis:
        slotno = int(sno)
        boardid = int(bid)
        slotdic[slotno] = boardid
    return slotdic

def _get_obsdate(filename, chassis=0):
    h,c,d = ReadAccumFiletoDic(filename)
    return parser.parse(h[chassis].values()[-1]['date'])

def _get_reduce_date(filename):
    return parser.parse(time.asctime(time.localtime(os.path.getmtime(filename))))

def make_corr_cal_nc(ncfile=None, chassis=0, direc='.', writedb=False,
                     copy=False):
    """
    Given a directory, this function finds all the
    corr_cal files in that directory in the form of
    corr_cal_ccid_chano_slotno.dat. For each file it
    reads the corr_cal matrix and creates a new netCDF file
    with all corr_cal matrices in it.
    If writedb is True, it will also file the data into the redshiftdb
    database into the corrcal table"""
    ccal = OrderedDict()
    slotdic = _get_bid(direc+os.path.sep+'calib200.dat', chassis=chassis)    
    time_units = "hours since 0001-01-01 00:00:00.0"
    calendar = "gregorian"
    obsdate = _get_obsdate(direc+os.path.sep+'calib200.dat', chassis=chassis)
    for slotno in (5, 8, 11, 14, 17, 20):
        files = glob.glob(direc+os.path.sep+'corr_cal_*_%d_%d.dat' % (chassis, slotno))
        if files:
            fname = files[0]
            ccid = os.path.basename(fname).split('_')[2]
            if os.path.exists(fname+'.info'):
                svd, zero_lag = _get_ccal_info(fname+'.info')
            else:
                svd, zero_lag = 0, 0
            corr_cal = numpy.loadtxt(fname)
            ccal[slotno] = (slotno, slotdic[slotno], ccid, svd, zero_lag, corr_cal)
            reducedate = _get_reduce_date(fname)
    if ncfile is None:
        ncfile = direc+os.path.sep+'corr_cal_%s_%d.nc' % (ccal.values()[0][2], chassis)
    nc = LMTNetCDFFile(ncfile, 'w')
    numcards = len(ccal.keys())
    dims = {
        "numchannels" : 256,
        "scalarval" : 1,
        "ccid_slen" : 9
        }
    create_dimensions(nc, dims, numcards=numcards)
    hvars = {
        "slotno" : {"type" : numpy.dtype(int),
                    "dim" : ("numcards",),
                    "attributes" : {"long_name" : "VME Slot number for correlator card"}
                    },
        "boardid" : {"type" : numpy.dtype(int),
                     "dim" : ("numcards",),
                     "attributes" : {"long_name": "Board ID as inscribed on board, and programmed on card."}
                 },
        "chassis" : {"type" : numpy.dtype(int),
                     "dim" : ("scalarval",),
                     "attributes" : {"long_name": "Chassis Number for the board (one of 0,1,2,3)"}
                     },    
        "corr_cal_id" : {"type" : numpy.dtype(int),
                         "dim" : ("numcards",),
                         "attributes" : {"long_name": "corr cal id for each board"}
                         },    
        "svd" : {"type" : numpy.dtype(int),
                 "dim" : ("numcards",),
                 "attributes" : {"long_name": "svd cut off value used in corr_cal calculation"}
                 },    
        "zerolag" : {"type" : numpy.dtype(int),
                     "dim" : ("numcards",),
                     "attributes" : {"long_name": "zero lag of board"}
                     },
        "obsdate" : {"type" : numpy.dtype(float),
                     "dim" : ("scalarval",),
                     "attributes" : {"long_name": "observed datetime for frequency calibration",
                                     "units" : time_units,
                                     "calendar" : calendar,
                                     },
                     },
        "reducedate" : {"type" : numpy.dtype(float),
                         "dim" : ("scalarval",),
                         "attributes" : {"long_name": "reduced datetime for frequency calibration",
                                         "units" : time_units,
                                         "calendar" : calendar,
                                         },
                     },
        "createdate" : {"type" : numpy.dtype(float),
                        "dim" : ("scalarval",),
                        "attributes" : {"long_name": "creation datetime for frequency calibration netcdf file",
                                        "units" : time_units,
                                        "calendar" : calendar,
                                        },
                     },        
        
        }
    createdate = datetime.datetime.now()
    hdic = create_variables(nc, hvars,
                            hname='RedshiftChassis_%d_' % chassis)
    dvars = {
        "corr_cal" : {"type" : numpy.dtype(float),
                      "dim" : ("numcards", "numchannels", "numchannels"),
                      "attributes" : {"long_name" : "corr_cal matrix"}
                      }
        }
    ddic = create_variables(nc, dvars, fstr='Data',
                            hname='RedshiftChassis_%d_' % chassis)
    
    #print ccal
    #fill in header variables
    nc.variables['Header.RedshiftChassis_%d_.slotno' % chassis][:] = numpy.array(ccal.keys(),
                                                                                 dtype=numpy.dtype(int))
    nc.variables['Header.RedshiftChassis_%d_.boardid' % chassis][:] = numpy.array([c[1] for c in ccal.values()],
                                                                                  dtype=numpy.dtype(int))
    #print numpy.array([c[1] for c in ccal.values()], dtype=numpy.dtype(int))
    nc.variables['Header.RedshiftChassis_%d_.corr_cal_id' % chassis][:] = numpy.array([int(c[2]) for c in ccal.values()],
                                                                                  dtype=numpy.dtype(int))
    nc.variables['Header.RedshiftChassis_%d_.chassis' % chassis][:] = chassis
    nc.variables['Header.RedshiftChassis_%d_.svd' % chassis][:] = numpy.array([c[3] for c in ccal.values()],
                                                                              dtype=numpy.dtype(int))
    nc.variables['Header.RedshiftChassis_%d_.zerolag' % chassis][:] = numpy.array([c[4] for c in ccal.values()],
                                                                              dtype=numpy.dtype(int))
    nc.variables['Header.RedshiftChassis_%d_.obsdate' % chassis][:] = date2num(obsdate, units=time_units, calendar=calendar)
    nc.variables['Header.RedshiftChassis_%d_.reducedate' % chassis][:] = date2num(reducedate, units=time_units, calendar=calendar)
    nc.variables['Header.RedshiftChassis_%d_.createdate' % chassis][:] = date2num(createdate, units=time_units, calendar=calendar)    
    
    for i, val in enumerate(ccal.values()):
        nc.variables['Data.RedshiftChassis_%d_.corr_cal' % chassis][i,:] = val[5]
    nc.close()
    print(ccal[5][1])
    
    if writedb:
        #write to CorrCal Table
        for slotno in ccal.keys():
            ccaltable = CorrCalTable(obsdate=obsdate, reducedate=reducedate,
                                     createdate=createdate, chassis=chassis,
                                     boardid=int(ccal[slotno][1]),
                                     slotno=slotno,
                                     svd=ccal[slotno][3], zerolag=ccal[slotno][4],
                                     corr_cal_id="%s" % ccal[slotno][2],
                                     filename=ncfile,
                                     valid=True)
            ccaltable.save()
    if copy:
        basefile = os.path.basename(ncfile)
        try:
            shutil.copyfile(ncfile, "/raw/rsr/cal" + os.path.sep + basefile)
        except:
            raise LMTRedshiftError("make_corr_cal_nc", "Error in copy of file %s to system directory /raw/rsr/cal" % basefile)
    return ncfile

