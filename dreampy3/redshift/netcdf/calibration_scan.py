"""RedshiftCalibration is a class to handle diagnostic or astronomical
CAL Scans. It contains methods and attributes that returns g_norm,
Tsys vectors, ratio vectors to decompress the data.
It takes as init arguments the ObsNum, chassis and whether to use
the db for the CAL observation"""
import glob
import os
import time
import numpy
import datetime

from .redshift_netcdf_file import RedshiftNetCDFFile
from dreampy3.redshift.utils import LMTRedshiftError
from dreampy3.redshift.utils.spectrum_utils import makespectrum, get_corr_cal_matrix, get_bad_lags
from dreampy3.redshift.utils.fileutils import make_generic_filename 
from dreampy3.logging import logger
from numpy import array

logger.name = __name__

diffsat = {'A': {1: array([  7.34152079e+00,  -1.79580840e-01,   1.24253497e-03]),
                 2: array([  7.93315138e+00,  -1.60071049e-01,   9.02031498e-04]),
                 3: array([  2.80288506e+01,  -6.53401895e-01,   3.91676190e-03]),
                 4: array([  1.22405665e+02,  -2.55689766e+00,   1.34491661e-02]),
                 5: array([  1.10510165e+01,  -1.65585867e-01,   6.59110670e-04]),
                 6: array([  3.98872506e+01,  -7.78625325e-01,   3.88134224e-03])},
           'B': {1: array([  3.24748980e+01,  -8.32696797e-01,   5.48390685e-03]),
                 2: array([  1.79916438e+01,  -3.87490619e-01,   2.18600744e-03]),
                 3: array([ -1.06681896e+01,   2.84144909e-01,  -1.73839248e-03]),
                 4: array([  3.21421487e+01,  -6.53854850e-01,   3.41765455e-03]),
                 5: array([  2.37064758e+01,  -4.17595082e-01,   1.90348556e-03]),
                 6: array([  4.63857755e+01,  -9.12545568e-01,   4.56504927e-03])},
           'C': {1: array([ -4.76731687e-01,   1.53251650e-02,   2.77181514e-05]),
                 2: array([  6.86855779e+01,  -1.50843761e+00,   8.37623012e-03]),
                 3: array([  1.90276497e+01,  -4.46213606e-01,   2.73114576e-03]),
                 4: array([  9.13517445e+01,  -1.91474237e+00,   1.01294210e-02]),
                 5: array([  6.57368094e+01,  -1.18706618e+00,   5.42288265e-03]),
                 6: array([  6.65131276e+01,  -1.31728143e+00,   6.59743826e-03])},
           'D': {1: array([  1.60296400e+01,  -4.01648845e-01,   2.65939840e-03]),
                 2: array([  2.41446468e+01,  -5.18848756e-01,   2.88828948e-03]),
                 3: array([  2.27700502e+01,  -5.33397060e-01,   3.24635137e-03]),
                 4: array([  6.13163689e+01,  -1.25404142e+00,   6.50543794e-03]),
                 5: array([ -1.84163344e+01,   3.61224637e-01,  -1.69379621e-03]),
                 6: array([  2.39658059e+01,  -4.74742570e-01,   2.43468229e-03])}
           }

def polyfunc(p, x):
    return p[0] + p[1]*x + p[2]*x**2

def nanmean(array):
    """Returns the mean of an array after ignoring the
    nan numbers"""
    return array[numpy.isfinite(array)].mean()
    
def bad_lags_acf(acf, chassis, board):
    badlags = get_bad_lags(chassis, board)
    acf = acf.copy()
    if badlags:
        logger.info('Setting bad lags for chassis %d board %d to %s' % (chassis, board, badlags))
    for blag in badlags:
        acf[blag] = numpy.nan
    return acf

def corr_acf(acf, scale_factor=1.3e6):
     acf = acf.copy()
     acf[:3] = acf[:3]*3.0
     acf = acf*(1. + (numpy.abs(acf)/scale_factor))
     acf[:3] = acf[:3]/3.0
     return acf

def ratio_blank(W, X, Y, Z, scale_factor=1.3e6):
     W_corr = corr_acf(W, scale_factor=scale_factor)
     X_corr = corr_acf(X, scale_factor=scale_factor)
     Y_corr = corr_acf(Y, scale_factor=scale_factor)
     Z_corr = corr_acf(Z, scale_factor=scale_factor)
     blank_corr = (W_corr - X_corr - Y_corr + Z_corr)/4.0
     blank = (W - X - Y + Z)/4.0
     return blank_corr/blank

def get_gnorm(W, X, Y, Z, scale_factor=1.3e6):
     W_corr = corr_acf(W, scale_factor=scale_factor)
     X_corr = corr_acf(X, scale_factor=scale_factor)
     Y_corr = corr_acf(Y, scale_factor=scale_factor)
     Z_corr = corr_acf(Z, scale_factor=scale_factor)
     blank_corr = (W_corr - X_corr - Y_corr + Z_corr)/4.0
     return nanmean(blank_corr)/blank_corr

class RedshiftCalibration(object):
    def __init__(self, ObsNum, chassis, usedb=False,
                 corr_linear=True,
                 chassis_pix_map={0 : 'A',
                                  1 : 'B',
                                  2 : 'C',
                                  3 : 'D'
                                  }):
        """
        corr_linear : correct for linearity issues in pixels
        """
        self.ObsNum = ObsNum
        self.chassis = chassis
        self.chas_str = 'RedshiftChassis_%d_' % self.chassis
        self.chassis_pix_map = chassis_pix_map
        self.corr_linear = corr_linear
        self.SingleFileCal = False
        if usedb:
            #here we do what it takes to access database and get the files
            pass
        else:
            #filenames = glob.glob('/data_lmt/RedshiftChassis%d/RedshiftChassis%d_*_%06d_00_0001.nc' % (self.chassis, self.chassis, self.ObsNum))
            filename = make_generic_filename(self.ObsNum, self.chassis)
            print("Cal filename", filename)
            #ftup = [(time.localtime(os.path.getmtime(f)), f) for f in filenames]
            #ftup.sort()
            #ftup.reverse()
            #fname = ftup[0][1]
            fname = filename
            #nc = RedshiftNetCDFFile(fname)
            nc = RedshiftNetCDFFile(fname)
            self.header = nc.hdu.header
            #eventually replace with real ambient temperature
            #this will be one of the values of Rsfend.Temperature
            #like self.Tamb = self.header.get('Rsfend.Temperature')[6]
            try:
                self.Tamb = self.header.get('Rsfend.Temperature')[6]
            except:
                self.Tamb = 279.0
            if self.Tamb < 100.0:
                #maybe GPIB card is not reading?
                self.Tamb = 279.0   #6 deg C
            if self.Tamb > 400:
                self.Tamb = 279.0
            self.datestr = self._get_datestr(fname)
            #if nc.hdu.header.get('Dcs.Receiver') != 'RedshiftReceiver':
            #    raise LMTRedshiftError('RedshiftCalibration', 'netcdf file %s not a RedshiftReceiver file' % fname)
            if nc.hdu.header.get('Dcs.ObsPgm') != 'Cal':
                raise LMTRedshiftError('RedshiftCalibration', 'netcdf file %s not a RedshiftReceiver Cal file' % fname)
            #self.CalMode = int(nc.hdu.header.get('Dcs.CalMode'))
            self.CalMode = int(nc.hdu.header.get('Cal.CalMode'))
            #self.CalMode = int(nc.hdu.header.get('%s.SwitchSettings' % self.chas_str, 0))
            self.ccal = get_corr_cal_matrix(nc.hdu.header)
            if self.corr_linear:
                nc.hdu._get_freqs()
                self.frequencies = nc.hdu.frequencies
            if nc.hdu.data.AccData.shape[0] == 6:
                self.SingleFileCal = True
                self.process_calibration(self.CalMode, nc)
            else:
                self.SingleFileCal = False
                nc.close()
                self.process_calibration(self.CalMode)

    def _get_datestr(self, fname):
        return os.path.basename(fname).split('_')[1]

    def _get_filename_for_scan(self, ScanNum):
        """returns the file name of a cal scan with a given scan number"""
        return "/raw_mnc/RedshiftChassis%d/RedshiftChassis%d_%s_%06d_00_%04d.nc" % (self.chassis, self.chassis, self.datestr, self.ObsNum, ScanNum)
    
    def _get_accdata(self, fname):
        nc = RedshiftNetCDFFile(fname)
        accdata = nc.hdu.data.AccData[0]/numpy.outer(nc.hdu.data.AccSamples[0], numpy.ones(256))  #just the first repeat
        nc.close()
        return accdata

    def _get_single_file_accdata(self, nc, index):
        accdata = nc.hdu.data.AccData[index, :, :]/numpy.outer(nc.hdu.data.AccSamples[index, :], numpy.ones(256))
        return accdata

    def _get_ratio(self, scans=(2,3,4,5)):
        W = self._get_accdata(self._get_filename_for_scan(scans[0]))
        X = self._get_accdata(self._get_filename_for_scan(scans[1]))
        Y = self._get_accdata(self._get_filename_for_scan(scans[2]))
        Z = self._get_accdata(self._get_filename_for_scan(scans[3]))
        ratio = numpy.ma.zeros(W.shape, dtype='float')
        for board in self.header.BoardNumber:
            logger.info('Processing ratio for board %d' % board)
            ratio[board, :] = ratio_blank(bad_lags_acf(W[board, :],
                                                       self.chassis, board),
                                          bad_lags_acf(X[board, :],
                                                       self.chassis, board),
                                          bad_lags_acf(Y[board, :],
                                                       self.chassis, board),
                                          bad_lags_acf(Z[board, :],
                                                       self.chassis, board),
                                          scale_factor=0.65e6)
        return ratio

    def _get_single_file_ratio(self, nc, scans=(2,3,4,5)):
        W = self._get_single_file_accdata(nc, scans[0])
        X = self._get_single_file_accdata(nc, scans[1])
        Y = self._get_single_file_accdata(nc, scans[2])
        Z = self._get_single_file_accdata(nc, scans[3])
        ratio = numpy.ma.zeros(W.shape, dtype='float')
        for board in self.header.BoardNumber:
            logger.info('Processing ratio for board %d' % board)
            ratio[board, :] = ratio_blank(bad_lags_acf(W[board, :],
                                                       self.chassis, board),
                                          bad_lags_acf(X[board, :],
                                                       self.chassis, board),
                                          bad_lags_acf(Y[board, :],
                                                       self.chassis, board),
                                          bad_lags_acf(Z[board, :],
                                                       self.chassis, board),
                                          scale_factor=0.65e6)
        return ratio

    def _get_gnorm(self, scans=(2,3,4,5)):
        W = self._get_accdata(self._get_filename_for_scan(scans[0]))
        X = self._get_accdata(self._get_filename_for_scan(scans[1]))
        Y = self._get_accdata(self._get_filename_for_scan(scans[2]))
        Z = self._get_accdata(self._get_filename_for_scan(scans[3]))
        g_norm = numpy.ma.zeros(W.shape, dtype='float')
        for board in self.header.BoardNumber:
            logger.info('Processing g_norm for board %d' % board)
            g_norm[board, :] = get_gnorm(bad_lags_acf(W[board, :],
                                                      self.chassis, board),
                                         bad_lags_acf(X[board, :],
                                                      self.chassis, board),
                                         bad_lags_acf(Y[board, :],
                                                      self.chassis, board),
                                         bad_lags_acf(Z[board, :],
                                                      self.chassis, board),
                                         scale_factor=1.3e6)
            # hack to fix faulty board 1 chassis 0 blanking 
            if int(self.header.ChassisNumber) == 0:
                if self.header.utdate().date() > datetime.date(2015, 12, 1):
                    g_norm[1, :] = numpy.ones(256, dtype='float')
        return g_norm

    def _get_single_file_gnorm(self, nc, scans=(2,3,4,5)):
        W = self._get_single_file_accdata(nc, scans[0])
        X = self._get_single_file_accdata(nc, scans[1])
        Y = self._get_single_file_accdata(nc, scans[2])
        Z = self._get_single_file_accdata(nc, scans[3])
        g_norm = numpy.ma.zeros(W.shape, dtype='float')
        for board in self.header.BoardNumber:
            logger.info('Processing g_norm for board %d' % board)
            g_norm[board, :] = get_gnorm(bad_lags_acf(W[board, :],
                                                      self.chassis, board),
                                         bad_lags_acf(X[board, :],
                                                      self.chassis, board),
                                         bad_lags_acf(Y[board, :],
                                                      self.chassis, board),
                                         bad_lags_acf(Z[board, :],
                                                      self.chassis, board),
                                         scale_factor=1.3e6)
        # hack to fix faulty board 1 chassis 0 blanking 
        if int(self.header.ChassisNumber) == 0:
            if self.header.utdate().date() > datetime.date(2015, 12, 1):
                g_norm[1, :] = g_norm[2, :]
        return g_norm

    def get_scaling(self, board):
        print("Getting scale for %d" % board)
        #print(self.frequencies.shape)
        scale = polyfunc(diffsat[self.chassis_pix_map[self.chassis]][board+1],
                         self.frequencies[board])
        return scale
    
    def process_diagnostic_calibration(self):
        #14 subscans
        logger.info('Processing diagnostic calibration')        
        hot_ratio = self._get_ratio(scans=(2,3,4,5))
        hotacf = self._get_accdata(self._get_filename_for_scan(1))
        sky_ratio = self._get_ratio(scans=(7,8,9,10))
        skyacf = self._get_accdata(self._get_filename_for_scan(6))
        self.hotspec = numpy.ma.zeros(hotacf.shape, dtype='float')
        self.skyspec = numpy.ma.zeros(hotacf.shape, dtype='float')
        self.Tsys = numpy.ma.zeros(hotacf.shape, dtype='float')
        for board in self.header.BoardNumber:
            badlags = get_bad_lags(self.chassis, board)
            self.hotspec[board, :] = makespectrum(hotacf[board, :],
                                                  ratio=hot_ratio[board, :],
                                                  badlags=badlags,
                                                  corr_cal=self.ccal[self.header.SlotNumber[board]])
            self.skyspec[board, :] = makespectrum(skyacf[board, :],
                                                  ratio=sky_ratio[board, :],
                                                  badlags=badlags,
                                                  corr_cal=self.ccal[self.header.SlotNumber[board]])
            self.Tsys[board, :] = self.Tamb*self.skyspec[board, :]/(self.hotspec[board, :] - self.skyspec[board, :])
            if self.corr_linear:
                scale = self.get_scaling(board)
                self.Tsys[board, :] = self.Tsys[board, :] * scale
        self.g_norm = self._get_gnorm(scans=(11,12,13,14))

    def process_astronomical_calibration(self, nc):
        #6 subscans
        if self.SingleFileCal:
            logger.info('Processing astronomical calibration using SingleFileCal')
        else:
            logger.info('Processing astronomical calibration using 6 files')
        if self.SingleFileCal:
            hotacf = self._get_single_file_accdata(nc, 0)
            skyacf = self._get_single_file_accdata(nc, 1)
        else:
            hotacf = self._get_accdata(self._get_filename_for_scan(1))
            skyacf = self._get_accdata(self._get_filename_for_scan(2))
        #print(type(hotacf))
        #print(hotacf)
        self.hotspec = numpy.ma.zeros(hotacf.shape, dtype='float')
        self.skyspec = numpy.ma.zeros(hotacf.shape, dtype='float')
        self.Tsys = numpy.ma.zeros(hotacf.shape, dtype='float')
        for board in self.header.BoardNumber:
            badlags = get_bad_lags(self.chassis, board)
            self.hotspec[board, :] = makespectrum(hotacf[board, :],
                                                  badlags=badlags,
                                                  corr_cal=self.ccal[self.header.SlotNumber[board]])
            #print(self.hotspec[board, :])
            self.skyspec[board, :] = makespectrum(skyacf[board, :],
                                                  badlags=badlags,
                                                  corr_cal=self.ccal[self.header.SlotNumber[board]])
            self.Tsys[board, :] = self.Tamb*self.skyspec[board, :]/(self.hotspec[board, :] - self.skyspec[board, :])
            #print(self.Tsys[board, :])
            #print(self.Tsys[board, :].shape)
            if self.corr_linear:
                scale = self.get_scaling(board)
                self.Tsys[board, :] = self.Tsys[board, :] * scale
        if self.SingleFileCal:
            self.skyratio = self._get_single_file_ratio(nc, (2, 3, 4, 5))
            self.g_norm = self._get_single_file_gnorm(nc, (2, 3, 4, 5))
            nc.nc.close()
        else:
            self.skyratio = self._get_ratio(scans=(3,4,5,6))
            self.g_norm = self._get_gnorm(scans=(3,4,5,6))
            
    def process_calibration(self, calmode, nc=None):
        if calmode == 0:
            self.process_diagnostic_calibration(nc)
        elif calmode == 1:
            self.process_astronomical_calibration(nc)
