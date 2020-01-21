"""
Class to read IFProc files for SEQUOIA, MSIP1mm receiver
and future receivers that uses IF processor
"""
from collections import OrderedDict
from dreampy3.lmtnetcdf import LMTNetCDFFile
from .ifprocdata import IFProcData
from .ifproc_header import IFProcHeader
from .ifproc_scan import IFProcScan
from dreampy3.ifproc.utils import LMTIFProcError
from netCDF4 import Dataset, Variable
import numpy
#from lmtnetcdf import LMTNetCDFFile
import netCDF4
import os
import time


class IFProcNetCDFFile(LMTNetCDFFile):

    def __init__(self, filename, mode='r'):
        if mode in ('r', 'a') and not os.path.exists(filename):
            raise LMTIFProcError('NetCDF Error', 'File %s not found' % filename)
        LMTNetCDFFile.__init__(self, filename, mode)
        self.filename = filename
        self.hdus = OrderedDict()
        if mode in ('r', 'a'):
            self._make_hdus()
        self.data_index = 0

    def _populate_headers(self, variables,
                          dimensions):
        return IFProcHeader(ncvariables=variables,
                            dimensions=dimensions)

    def _populate_data(self, variables):
        return IFProcData(variables)

    def _make_hdus(self):
        #only one single root group
        data = self._populate_data(self.nc.variables)
        header = self._populate_headers(self.nc.variables, self.nc.dimensions)
        self.hdus['rootgrp'] = IFProcScan(data=data, header=header,
                                          filename=self.filename)
        self.hdu = list(self.hdus.values())[0]

