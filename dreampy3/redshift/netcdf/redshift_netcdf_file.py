"""Various classes to do with specifically
the reading and writing of LMT Redshift Receiver
netCDF files which are typically netCDF4 files.

As read from the LMT M&C data writer system, the files are
typically netCDF3. Inside dreampy's redshift module,
we like to store the intermediate files as netCDF4 files
"""

from dreampy3.lmtnetcdf import LMTNetCDFFile, NetCDFFile, LMThdu
from collections import OrderedDict
from .redshift_data import RedshiftData
from .redshift_header import RedshiftHeader, RedshiftCorrCalHeader
from .redshift_scan import RedshiftScan
from dreampy3.redshift.utils.redshift_netcdf_variables import dimensions, header_variables_dict, data_variables_dict
from dreampy3.redshift.utils import LMTRedshiftError
from dreampy3.redshift.utils.spectrum_utils import get_bad_lags
import os
import numpy
from dreampy3.logging import logger
logger.name = __name__

class RedshiftNetCDFFile(LMTNetCDFFile):
    """
    A RedshiftNetCDFFile class opens up a redshift receiver
    netCDF file. If it is version 3, it creates a list of hdus
    with just one element, that made from the root directory of
    the netCDF file.
    If it is version 4, it creates a list of hdus, populated with all
    the header and data items from each rootgroup of the netcdf4 file.
    """
    def __init__(self, filename, mode='r'):
        if mode in ('r', 'a') and not os.path.exists(filename):
            raise LMTRedshiftError('Redshift NetCDF Error', 'File %s not found' % filename)
        LMTNetCDFFile.__init__(self, filename, mode)
        self.filename = filename
        self.hdus = OrderedDict()
        if mode in ('r', 'a'):
            self._make_hdus()

    def _populate_headers(self, variables,
                          dimensions):
        return RedshiftHeader(ncvariables=variables,
                              dimensions=dimensions)

    def _populate_data(self, variables):
        return RedshiftData(variables)
    
    def _make_hdus(self):
        if len(self.groups) > 0:
            #more than one group available,
            #but don't open everything just the first one
            #so you don't run out of memory
            i = 0
            for key, group in self.groups.items():
                if i == 0:
                    data = self._populate_data(group.variables)
                    header = self._populate_headers(group.variables, group.dimensions)
                    self.hdus[key] = RedshiftScan(data=data, header=header,
                                                  filename=self.filename)
                else:
                    self.hdus[key] = None
                i += 1
        else:
            #only one single root group
            data = self._populate_data(self.variables)
            header = self._populate_headers(self.variables, self.dimensions)
            self.hdus['rootgrp'] = RedshiftScan(data=data, header=header,
                                                filename=self.filename)
        self.hdu = list(self.hdus.values())[0]

    def _get_group_hdu(self, groupname):
        if not self.groups.has_key(groupname):
            raise LMTRedshiftError('make_group_hdu', 'Group with name %s not found' % groupname)
        group = self.groups[groupname]
        data = self._populate_data(group.variables)
        header = self._populate_headers(group.variables, group.dimensions)
        return RedshiftScan(data=data, header=header,
                            filename=self.filename,
                            groupname=groupname)
        
    def copy_dimensions(self, ncin, hdu, ncout):
        """Copy dimension object from ncin to ncout"""
        for k, dim in ncin.dimensions.items():
            if k != 'Header.Spectrum.compwindows_xlen':
                ncout.createDimension(k, len(dim))
        if hasattr(hdu, 'spectrum'):
            #has processed spectra, so let us save it appropriately
            #create dimensions for badlags
            chassis = int(hdu.header.ChassisNumber)
            for board in hdu.header.BoardNumber:
                blags = get_bad_lags(chassis, board)
                if blags:
                    if 'Header.Spectrum.BadLags%d' % board not in ncout.dimensions:
                        ncout.createDimension('Header.Spectrum.BadLags%d' % board,
                                              len(blags))
            #create dimensions for Tsys and tint
            if hasattr(hdu, 'Tsys'):
                if 'Header.Spectrum.Tsys_xlen' not in ncout.dimensions:
                    ncout.createDimension('Header.Spectrum.Tsys_xlen',
                                          len(hdu.Tsys))
            if hasattr(hdu, 'tint'):
                if 'Header.Spectrum.tint' not in ncout.dimensions:
                    ncout.createDimension('Header.Spectrum.tint',
                                          1)
            #create dimensions for spectral and frequency data
            if 'Data.Spectrum.spectrum_rlen' not in ncout.dimensions:
                ncout.createDimension('Data.Spectrum.spectrum_rlen',
                                      hdu.spectrum.shape[0])
            if 'Data.Spectrum.spectrum_xlen' not in ncout.dimensions:
                ncout.createDimension('Data.Spectrum.spectrum_xlen',
                                      len(hdu.header.BoardNumber))
            if 'Data.Spectrum.spectrum_ylen' not in ncout.dimensions:
                ncout.createDimension('Data.Spectrum.spectrum_ylen',
                                      hdu.spectrum.shape[2])
        if hasattr(hdu, 'compspectrum'):
            #has processed compspectrum
            chassis = int(hdu.header.ChassisNumber)
            if 'Data.Spectrum.compspectrum_rlen' not in ncout.dimensions:
                ncout.createDimension('Data.Spectrum.compspectrum_rlen',
                                      hdu.compspectrum.shape[0])
            if 'Data.Spectrum.compspectrum_xlen' not in ncout.dimensions:
                ncout.createDimension('Data.Spectrum.compspectrum_xlen',
                                      hdu.compspectrum.shape[1])
        if hasattr(hdu, 'windows') and hdu.windows:
            #baseline possibly taken out, so we may consider
            #multiple windows in dimensions and headers
            for board in hdu.header.BoardNumber:
                if hdu.windows[board]:
                    if 'Header.Spectrum.windows%d_xlen' % board not in ncout.dimensions:
                        ncout.createDimension('Header.Spectrum.windows%d_xlen' % board,
                                              len(hdu.windows[board]))
                    if 'Header.Spectrum.windows_ylen' not in ncout.dimensions:
                        ncout.createDimension('Header.Spectrum.windows_ylen', 2)
            #create dimension for sigma
            if 'Header.Spectrum.sigma_rlen' not in ncout.dimensions:
                ncout.createDimension('Header.Spectrum.sigma_rlen',
                                      hdu.spectrum.shape[0])
            if 'Header.Spectrum.sigma_xlen' not in ncout.dimensions:
                ncout.createDimension('Header.Spectrum.sigma_xlen',
                                      len(hdu.header.BoardNumber))
        if hasattr(hdu, 'compwindows') and hdu.compwindows:
            if 'Header.Spectrum.compwindows_xlen' not in ncout.dimensions:
                ncout.createDimension('Header.Spectrum.compwindows_xlen',
                                      len(hdu.compwindows))
            if 'Header.Spectrum.compwindows_ylen' not in ncout.dimensions:
                ncout.createDimension('Header.Spectrum.compwindows_ylen',
                                      2)
            if hasattr(hdu, 'compsigma') and hdu.compsigma.any():
                if 'Header.Spectrum.compsigma_xlen' not in ncout.dimensions:
                    ncout.createDimension('Header.Spectrum.compsigma_xlen',
                                          hdu.compsigma.shape[0])
            
    def copy_header(self, ncin, hdu, ncout):
        """Copy header netCDF variables from ncin to ncout"""
        heads = [name for name in ncin.variables.keys() if name.find('Header') != -1]
        for head in heads:
            #logger.info("Creating header %s" % head)
            if head != 'Header.Spectrum.compwindows':
                var = ncout.createVariable(head,
                                           ncin.variables[head].dtype,
                                           dimensions=ncin.variables[head].dimensions)
                var[:] = ncin.variables[head][:]
                #copy the ncattrs as well
                for attr in ncin.variables[head].ncattrs():
                    var.setncattr(attr, getattr(ncin.variables[head], attr))
        if hasattr(hdu, 'spectrum'):
            #has spectral data so save some headers
            #eventually let us save baseline and window data here
            chassis = int(hdu.header.ChassisNumber)
            for board in hdu.header.BoardNumber:
                blags = get_bad_lags(chassis, board)
                if blags:
                    if not ncout.variables.has_key('Header.Spectrum.BadLags%d' % board):
                        blvar = ncout.createVariable('Header.Spectrum.BadLags%d' % board,
                                                     numpy.dtype('int32'),
                                                     dimensions='Header.Spectrum.BadLags%d' % board)
                    else:
                        blvar = ncout.variables['Header.Spectrum.BadLags%d' % board]
                    blvar[:] = numpy.array(blags)
            if hasattr(hdu, 'Tsys'):
                if not ncout.variables.has_key('Header.Spectrum.Tsys'):
                    tsysvar = ncout.createVariable('Header.Spectrum.Tsys',
                                                   numpy.dtype('float'),
                                                   dimensions='Header.Spectrum.Tsys_xlen')
                else:
                    tsysvar = ncout.variables['Header.Spectrum.Tsys']
                tsysvar[:] = hdu.Tsys
            if hasattr(hdu, 'tint'):
                if not ncout.variables.has_key('Header.Spectrum.tint'):
                    tintvar = ncout.createVariable('Header.Spectrum.tint',
                                                   numpy.dtype('float'),
                                                   dimensions='Header.Spectrum.tint')
                else:
                    tintvar = ncout.variables['Header.Spectrum.tint']
                tintvar[:] = hdu.tint
        if hasattr(hdu, 'windows') and hdu.windows:
            for board in hdu.header.BoardNumber:
                if hdu.windows[board]:
                    if 'Header.Spectrum.windows%d' % board not in ncout.variables:
                        winvar = ncout.createVariable('Header.Spectrum.windows%d' % board,
                                                      numpy.dtype('float'),
                                                      dimensions=('Header.Spectrum.windows%d_xlen' % board,
                                                                  'Header.Spectrum.windows_ylen'))
                    else:
                        winvar = ncout.variables['Header.Spectrum.windows%d' % board]
                    winvar[:] = numpy.array(hdu.windows[board])
            if hasattr(hdu, 'order'):
                if 'Header.Spectrum.order' not in ncout.variables:
                    ordervar = ncout.createVariable('Header.Spectrum.order',
                                                    numpy.dtype('int32'),
                                                    dimensions=())
                else:
                    ordervar = ncout.variables['Header.Spectrum.order']
                ordervar[:] = hdu.order
            if hasattr(hdu, 'sigma'):
                if 'Header.Spectrum.sigma' not in ncout.variables:
                    sigmavar = ncout.createVariable('Header.Spectrum.sigma',
                                                    numpy.dtype('float'),
                                                    dimensions=('Header.Spectrum.sigma_rlen',
                                                                'Header.Spectrum.sigma_xlen'))
                else:
                    sigmavar = ncout.variables['Header.Spectrum.sigma']
                #print sigmavar.shape
                #sigarr = numpy.zeros((ncin.hdu.spectrum.shape[0],
                #                      len(ncin.hdu.header.BoardNumber)),
                #                     dtype=numpy.dtype('float'))
                                 
                #for rpt in range(ncin.hdu.spectrum.shape[0]):
                #    for board in ncin.hdu.header.BoardNumber:
                #        sigarr[rpt, board] = ncin.hdu.sigma[board][0, rpt]
                sigmavar[:] = hdu.sigma
        if hasattr(hdu, 'compwindows') and hdu.compwindows:
            if 'Header.Spectrum.compwindows' not in ncout.variables:
                compwinvar = ncout.createVariable('Header.Spectrum.compwindows',
                                                  numpy.dtype('float'),
                                                  dimensions=('Header.Spectrum.compwindows_xlen',
                                                              'Header.Spectrum.compwindows_ylen'))
            else:
                compwinvar = ncout.variables['Header.Spectrum.compwindows']
            print(compwinvar.shape)
            compwinvar[:] = numpy.array(hdu.compwindows)
            if hasattr(hdu, 'compsigma') and hdu.compsigma.any():
                if 'Header.Spectrum.compsigma' not in ncout.variables:
                    compsigvar = ncout.createVariable('Header.Spectrum.compsigma',
                                                      numpy.dtype('float'),
                                                      dimensions=('Header.Spectrum.compsigma_xlen'))
                else:
                    compsigvar = ncout.variables['Header.Spectrum.compsigma']
                compsigvar[:] = hdu.compsigma
                
                                                  
    def copy_data(self, ncin, hdu, ncout, remove_data=['BufPos']):
        """Copy data netCDF variables from ncin to ncout.
        Do not include items in remove_data list in ncout"""
        datas = [name for name in ncin.variables.keys() if name.find('Data') != -1]
        dvar = {}
        for data in datas:
            dname = data.split('.')[-1]
            if dname not in remove_data:
                #print dname
                dvar[dname] = ncout.createVariable(data,
                                                   ncin.variables[data].dtype,
                                                   dimensions=ncin.variables[data].dimensions)
                dvar[dname][:] = ncin.variables[data][:]
                #copy the attrs
                for attr in ncin.variables[data].ncattrs():
                    dvar[dname].setncattr(attr,
                                          getattr(ncin.variables[data], attr))
        if hasattr(hdu, 'spectrum'):
            #has spectral data so save some data values
            if not ncout.variables.has_key('Data.Spectrum.spectrum'):
                dvar['Data.Spectrum.spectrum'] = ncout.createVariable('Data.Spectrum.spectrum',
                                                                      numpy.dtype(float),
                                                                      dimensions=('Data.Spectrum.spectrum_rlen',
                                                                                  'Data.Spectrum.spectrum_xlen',
                                                                                  'Data.Spectrum.spectrum_ylen')
                                                                      
                                                                      )
            else:
                dvar['Data.Spectrum.spectrum'] = ncout.variables['Data.Spectrum.spectrum']
            dvar['Data.Spectrum.spectrum'][:] = hdu.spectrum
            if 'Data.Spectrum.frequencies' not in ncout.variables:
                dvar['Data.Spectrum.frequencies'] = ncout.createVariable('Data.Spectrum.frequencies',
                                                                         numpy.dtype(float),
                                                                         dimensions=('Data.Spectrum.spectrum_xlen',
                                                                                     'Data.Spectrum.spectrum_ylen')
                                                                         )
            else:
                dvar['Data.Spectrum.frequencies'] = ncout.variables['Data.Spectrum.frequencies']
            dvar['Data.Spectrum.frequencies'][:] = hdu.frequencies
        if hasattr(hdu, 'compspectrum'):
            #has compspectral data so save some data values
            if 'Data.Spectrum.compspectrum' not in ncout.variables:
                dvar['Data.Spectrum.compspectrum'] = ncout.createVariable('Data.Spectrum.compspectrum',
                                                                          numpy.dtype(float),
                                                                          dimensions=('Data.Spectrum.compspectrum_rlen',
                                                                                      'Data.Spectrum.compspectrum_xlen')
                                                                          
                                                                          )
            else:
                dvar['Data.Spectrum.compspectrum'] = ncout.variables['Data.Spectrum.compspectrum']
            dvar['Data.Spectrum.compspectrum'][:] = hdu.compspectrum
            if 'Data.Spectrum.compfreq' not in ncout.variables:
                dvar['Data.Spectrum.compfreq'] = ncout.createVariable('Data.Spectrum.compfreq',
                                                                      numpy.dtype(float),
                                                                      dimensions=('Data.Spectrum.compspectrum_xlen',)
                                                                      )
            else:
                dvar['Data.Spectrum.compfreq'] = ncout.variables['Data.Spectrum.compfreq']
            dvar['Data.Spectrum.compfreq'][:] = hdu.compfreq

    def append_observation(self, ncout, groupname=None):
        """This method will append the current netCDF3 based single
        observation to the output ncout netCDF4 file.
        It is assumed that ncout is an open writeable netCDF4 file.
        A new rootgroup will be created in ncout for each netCDF3
        input observation"""
        if groupname is None:
            groupname = "%s_%d_%d_%d_C%d" % (self.hdu.header.utdate().strftime('%Y-%m-%d'),
                                             int(self.hdu.header.ObsNum),
                                             int(self.hdu.header.SubObsNum),
                                             int(self.hdu.header.ScanNum),
                                             int(self.hdu.header.ChassisNumber),
                                             )
        logger.info("Creating Groupname %s" % groupname)
        group = ncout.createGroup(groupname)
        self.copy_dimensions(self, self.hdu, group)
        self.copy_header(self, self.hdu, group)
        self.copy_data(self, self.hdu, group)
        group.sync()
        ncout.sync()

    def save_hdu(self, group, hdu, groupname=None):
        """Saves a group and associated hdu in the opened
        netcdf file"""
        if groupname is None:
            groupname = "%s_%d_%d_%d_C%d" % (hdu.header.utdate().strftime('%Y-%m-%d'),
                                             int(hdu.header.ObsNum),
                                             int(hdu.header.SubObsNum),
                                             int(hdu.header.ScanNum),
                                             int(hdu.header.ChassisNumber),
                                             )
        logger.info("Creating Groupname %s" % groupname)
        outgroup = self.createGroup(groupname)
        self.copy_dimensions(group, hdu, outgroup)
        self.copy_header(group, hdu, outgroup)
        self.copy_data(group, hdu, outgroup)
        outgroup.sync()
        self.sync()

        
    def average_all_hdus(self, weight='sigma', threshold_sigma=None):
        """Given a list of hdus in a netCDF4 file this method
        will average all the scans and returns a representative group
        and new hdu with the  average of all hdus"""
        avghdu = RedshiftScan(data=self.hdu.data, header=self.hdu.header,
                              filename=self.filename)
        avggroup = self.groups[self.groups.keys()[0]]
        firsthdu = self._get_group_hdu(self.groups.keys()[0])
        frequencies = firsthdu._get_frequencies()
        spec = firsthdu._get_spectrum()
        avgspectrum = numpy.zeros((1, spec.shape[1], spec.shape[2]),
                                  dtype=spec.dtype)
        weights = {}
        tsys = {}
        tint = 0.0
        for board in range(frequencies.shape[0]):
            weights[board] = 0.0
            tsys[board] = 0.0
        for gname in self.groups.keys():
            hdu = self._get_group_hdu(gname)
            spec = hdu._get_spectrum()
            htsys = hdu._get_tsys()
            logger.info("Processing hdu in %s" % gname)
            count = 0
            for board in range(frequencies.shape[0]):
                wt = hdu._calc_weight(0, board, weight=weight)
                if weight == 'sigma' and threshold_sigma is not None:
                    #print gname, board, wt, threshold_sigma
                    if 1/numpy.sqrt(wt) < threshold_sigma:
                        avgspectrum[0, board, :] += spec[0, board, :] * wt
                        tsys[board] += htsys[board] * wt
                        weights[board] += wt
                        count += 1
                    else:
                        logger.info("Rejecting Board %d group %s due to sigma threshold" % (board, gname))
                else:
                    avgspectrum[0, board, :] += spec[0, board, :] * wt
                    tsys[board] += htsys[board] * wt
                    weights[board] += wt
                    count += 1
            if count != 0:
                tint += hdu._get_tint()
        for board in range(frequencies.shape[0]):
            avgspectrum[0, board, :] /= weights[board]
            tsys[board] /= weights[board]
        firsthdu.spectrum = avgspectrum
        firsthdu.frequencies = frequencies
        firsthdu.tsys = tsys
        firsthdu.tint = tint
        return avggroup, firsthdu

    def average_all_hdus_compspectrum(self, weight='sigma', threshold_sigma=None):
        """Given a list of hdus in a netCDF4 file this method
        will average all the compspectrum scans and returns a representative group
        and new hdu with the  average of all hdus"""
        avghdu = RedshiftScan(data=self.hdu.data, header=self.hdu.header,
                              filename=self.filename)
        avggroup = self.groups[self.groups.keys()[0]]
        firsthdu = self._get_group_hdu(self.groups.keys()[0])
        compfreq = firsthdu._get_comp_freq()
        spec = firsthdu._get_comp_spectrum()
        avgspectrum = numpy.zeros((1, spec.shape[1]),
                                  dtype=spec.dtype)

        weights = 0.0
        for gname in self.groups.keys():
            logger.info("Processing hdu in %s" % gname)
            hdu = self._get_group_hdu(gname)
            spec = hdu._get_comp_spectrum()
            wt = hdu._calc_weight_comp(0, weight=weight)
            if weight == 'sigma' and threshold_sigma is not None:
                print(gname, wt, threshold_sigma)
                if 1/numpy.sqrt(wt) < threshold_sigma:
                    avgspectrum[0, :] += spec[0, :] * wt
                    weights += wt
                else:
                    logger.info("Rejecting group %s due to sigma threshold" % gname)
            else:
                avgspectrum[0, :] += spec[0, :] * wt
                weights += wt                    
        avgspectrum[0, :] /= weights
        firsthdu.compspectrum = avgspectrum
        firsthdu.compfreq = compfreq
        return avggroup, firsthdu

    def export_to_sdfits(self, sdfilename, obs_list=None):
        """
        Convert nc data to fits format and save as a fits file. (sdfilename.fits)
        If nc is hetergeneous (except if all 'Cal' and one other type), ask user whether to save separate fits files.
        If user says no, ask which scan type to save.
        If user does not properly specify, use first non-cal scan type.
        """
        from dreampy.redshift.io.fits.redshift_sdfits import RedshiftSDFITS
        def find_obspgms(nc):
            """
            If all scans are the same (or all 'Cal' and one other type), return the type of scan.
            If there are different types of scans, return a list with number of scans of each type.
            """
            on = []
            x = 0
            bs = []
            y = 0
            ps = []
            z = 0
            cal = []
            w = 0
            for groupname in self.hdus.keys():
                hdu = self._get_group_hdu(groupname)
                if 'On' in hdu.header.keys():
                    on.append(groupname)
                    x = 1
                if 'Bs' in hdu.header.keys():
                    bs.append(groupname)
                    y = 1
                if 'Ps' in hdu.header.keys():
                    ps.append(groupname)
                    z = 1
                if 'Cal' in hdu.header.keys():
                    cal.append(groupname)
                    w = 1
            if x + y + z == 0:
                if w == 1:
                    raise LMTRedshiftError('Redshift NetCDF Error', 'nc contains only Cal scans. Cal scan export to SDFits not supported.')
                else:
                    raise LMTRedshiftError('Redshift NetCDF Error', 'Unknown scan types in nc.')
            elif x + y + z > 1:
                print("There are %d 'On' scans, %d 'Bs' scans, %d 'Ps' scans, and %d 'Cal' scans in this nc." % (len(on),len(bs),len(ps),len(cal)))
                return [len(on),len(bs),len(ps),len(cal)]
            elif w == 1:
                if on:
                    print("All scans are 'On' scans and 'Cal' Scans. 'Cal' scan export to SDFits not supported; exporting only 'On' scans.")
                    return 'On'
                if bs:
                    print("All scans are 'Bs' scans and 'Cal' Scans. 'Cal' scan export to SDFits not supported; exporting only 'Bs' scans.")
                    return 'Bs'
                if ps:
                    print("All scans are 'Ps' scans and 'Cal' Scans. 'Cal' scan export to SDFits not supported; exporting only 'Ps' scans.")
                    return 'Ps'
            else:
                if on:
                    print("All scans are 'On' scans.")
                    return 'On'
                if bs:
                    print("All scans are 'Bs' scans.")
                    return 'Bs'
                if ps:
                    print("All scans are 'Ps' scans.")
                    return 'Ps'


        check_sdfilename = sdfilename.split('.')
        if check_sdfilename[-1] == 'fits':
            junk_variable = check_sdfilename.pop()
            sdfilename = '.'.join(check_sdfilename)
        if os.path.exists(sdfilename+'.fits'):
            raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
        if (type(obs_list) != list) and (obs_list != None):
            raise LMTRedshiftError('Redshift NetCDF Error',"obs_list must be a list.")
        if (self.file_format == 'NETCDF3_CLASSIC') or (len(self.hdus) == 1):
           sdfits = RedshiftSDFITS(self)
           sdfits.save_fits(sdfilename+'.fits')
        else:
            obspgm = find_obspgms(self)
            if obs_list != None:
                if ('Bs' in obs_list) or ('Ps' in obs_list) or ('On' in obs_list):
                    any_scans = False
                    if ('On' in obs_list) and (obspgm[0] > 0):
                        if os.path.exists(sdfilename+'_On.fits'):
                            print("%s_On.fits already exists" % sdfilename)
                            raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
                        sdfits_On = RedshiftSDFITS(self, obspgm='On')
                        sdfits_On.save_fits(sdfilename+'_On.fits')
                        print("On scans exported to %s_On.fits" % sdfilename)
                        any_scans = True
                    if ('Bs' in obs_list) and (obspgm[1] > 0):
                        if os.path.exists(sdfilename+'_Bs.fits'):
                            print("%s_Bs.fits already exists" % sdfilename)
                            raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
                        sdfits_Bs = RedshiftSDFITS(self, obspgm='Bs')
                        sdfits_Bs.save_fits(sdfilename+'_Bs.fits')
                        print("Bs scans exported to %s_Bs.fits" % sdfilename)
                        any_scans = True
                    if ('Ps' in obs_list) and (obspgm[2] > 0):
                        if os.path.exists(sdfilename+'_Ps.fits'):
                            print("%s_Ps.fits already exists" % sdfilename)
                            raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
                        sdfits_Ps = RedshiftSDFITS(self, obspgm='Ps')
                        sdfits_Ps.save_fits(sdfilename+'_Ps.fits')
                        print("Ps scans exported to %s_Ps.fits" % sdfilename)
                        any_scans = True
                    if any_scans == False:
                        print("No scans of the type(s) you wanted.")
                else:
                    raise LMTRedshiftError('Redshift NetCDF Error',"obs_list must contain at least one of the strings 'On','Bs', and 'Ps'.")
            elif type(obspgm) == str:
                if obspgm == 'Cal':
                    raise LMTRedshiftError('Redshift NetCDF Error', "'Cal' scan export to SDFits not supported.")
                else:
                    sdfits = RedshiftSDFITS(self, obspgm=obspgm)
                    sdfits.save_fits(sdfilename+'.fits')
            else:
                multiple_files = raw_input('Would you like to save each scan type as its own SDFits file? (excluding Cal) (y/n): ')
                multiple_files = multiple_files.lower()
                if multiple_files == 'y' or multiple_files == 'yes':
                    if obspgm[0]:
                        if os.path.exists(sdfilename+'_On.fits'):
                            print("%s_On.fits already exists" % sdfilename)
                            raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
                        sdfits_On = RedshiftSDFITS(self, obspgm='On')
                        sdfits_On.save_fits(sdfilename+'_On.fits')
                        print("On scans exported to %s_On.fits" % sdfilename)
                    if obspgm[1]:
                        if os.path.exists(sdfilename+'_Bs.fits'):
                            print("%s_Bs.fits already exists" % sdfilename)
                            raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
                        sdfits_Bs = RedshiftSDFITS(self, obspgm='Bs')
                        sdfits_Bs.save_fits(sdfilename+'_Bs.fits')
                        print("Bs scans exported to %s_Bs.fits" % sdfilename)
                    if obspgm[2]:
                        if os.path.exists(sdfilename+'_Ps.fits'):
                            print("%s_Ps.fits already exists" % sdfilename)
                            raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
                        sdfits_Ps = RedshiftSDFITS(self, obspgm='Ps')
                        sdfits_Ps.save_fits(sdfilename+'_Ps.fits')
                        print("Ps scans exported to %s_Ps.fits" % sdfilename)
                    if obspgm[3]:
                        print("Reminder: 'Cal' scans not saved. ('Cal' scan export to SDFits not supported.)")
                elif obspgm[0] == obspgm[1] == obspgm[2] == 1:
                    two_scans = raw_input('Would you liked to save two scan types but not the third? (y/n): ')
                    if two_scans == 'y' or two_scans == 'yes':
                        dont_save = raw_input('Which scan type would you like NOT to save? (On, Bs, Ps -- Case Sensitive): ')
                        if (dont_save != 'On') and (dont_save != 'Bs') and (dont_save != 'Ps'):
                            print("Scan type entered incorrectly. Enter only two letters. First should be uppercase and second should be lowercase.")
                            dont_save = raw_input('One more try... Which scan type would you like NOT to save? (On, Bs, Ps -- Case Sensitive): ')
                        if (dont_save != 'On') and (dont_save != 'Bs') and (dont_save != 'Ps'):
                            raise LMTRedshiftError('Redshift NetCDF Error',"Scan type entered incorrectly.")
                        if dont_save != 'On':
                            if os.path.exists(sdfilename+'_On.fits'):
                                print("%s_On.fits already exists" % sdfilename)
                                raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
                            sdfits_On = RedshiftSDFITS(self, obspgm='On')
                            sdfits_On.save_fits(sdfilename+'_On.fits')
                            print("On scans exported to %s_On.fits" % sdfilename)
                        if dont_save != 'Bs':
                            if os.path.exists(sdfilename+'_Bs.fits'):
                                print("%s_Bs.fits already exists" % sdfilename)
                                raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
                            sdfits_Bs = RedshiftSDFITS(self, obspgm='Bs')
                            sdfits_Bs.save_fits(sdfilename+'_Bs.fits')
                            print("Bs scans exported to %s_Bs.fits" % sdfilename)
                        if dont_save != 'Ps':
                            if os.path.exists(sdfilename+'_Ps.fits'):
                                print("%s_Ps.fits already exists" % sdfilename)
                                raise LMTRedshiftError('Redshift NetCDF Error', "Filename already exists.")
                            sdfits_Ps = RedshiftSDFITS(self, obspgm='Ps')
                            sdfits_Ps.save_fits(sdfilename+'_Ps.fits')
                            print("Ps scans exported to %s_Ps.fits" % sdfilename)
                    else:
                        observepgm = raw_input('Select which scan type to save. (On, Bs, Ps -- Case Sensitive): ')
                        if (observepgm == 'On') or (observepgm == 'Bs') or (observepgm == 'Ps'):
                            sdfits = RedshiftSDFITS(self, obspgm=observepgm)
                            sdfits.save_fits(sdfilename+'.fits')
                            print("%s scans exported to %s.fits" % (observepgm, sdfilename))
                        else:
                            print("Saving only scans with obspgm of first scan that isn't a cal scan.")
                            sdfits = RedshiftSDFITS(self)
                            sdfits.save_fits(sdfilename+'.fits')
                            print("%s scans exported to %s.fits" % (sdfits.obspgm, sdfilename))
                else:
                    observepgm = raw_input('Select which scan type to save. (On, Bs, Ps -- Case Sensitive): ')
                    if (observepgm == 'On') or (observepgm == 'Bs') or (observepgm == 'Ps'):
                        sdfits = RedshiftSDFITS(self, obspgm=observepgm)
                        sdfits.save_fits(sdfilename+'.fits')
                        print("%s scans exported to %s.fits" % (observepgm, sdfilename))
                    else:
                        print("Saving only scans with obspgm of first scan that isn't a cal scan.")
                        sdfits = RedshiftSDFITS(self)
                        sdfits.save_fits(sdfilename+'.fits')
                        print("%s scans exported to %s.fits" % (sdfits.obspgm, sdfilename))
