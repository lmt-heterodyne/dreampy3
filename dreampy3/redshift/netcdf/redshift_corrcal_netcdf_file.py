from dreampy3.lmtnetcdf import LMTNetCDFFile, NetCDFFile, LMThdu
from .redshift_header import RedshiftCorrCalHeader
#from dreampy3.utils import OrderedDict
from collections import OrderedDict
from .redshift_data import RedshiftData
import os

class RedshiftCorrCalFile(LMTNetCDFFile):
    """
    A RedshiftCorrCalFile class opens up a redshift receiver
    corrcal netCDF file.
    """
    def __init__(self, filename):
        if not os.path.exists(filename):
            raise LMTRedshiftError('Redshift CorrCalFile Error', 'File %s not found' % filename)
        LMTNetCDFFile.__init__(self, filename, mode='r')
        self.filename = filename
        self.hdus = OrderedDict()
        self._make_hdus()

    def _populate_headers(self, variables,
                          dimensions):
        return RedshiftCorrCalHeader(ncvariables=variables,
                                     dimensions=dimensions)

    def _populate_data(self, variables):
        return RedshiftData(variables)

    def _make_hdus(self):
        #only one single root group
        data = self._populate_data(self.variables)
        header = self._populate_headers(self.variables, self.dimensions)
        self.hdus['rootgrp'] = LMThdu(data=data, header=header,
                                      filename=self.filename)
        self.hdu = list(self.hdus.values())[0]
