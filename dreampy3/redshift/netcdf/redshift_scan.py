"""
The RedshiftScan class is a simple extension to the
LMThdu class"""

import numpy
from dreampy3.lmtnetcdf import LMThdu
from dreampy3.redshift.utils.spectrum_utils import makespectrum, get_corr_cal_matrix, get_bad_lags, apply_bad_lags, nanmean, blank, blank_to_zero, blank_to_nan, get_blank_indices, get_blank_indices_array
from dreampy3.redshift.utils import LMTRedshiftError
#import types
from scipy.interpolate.interpolate import interp1d
#from scipy.stats import nanmean as scipynanmean
from numpy import nanmean as scipynanmean
from dreampy3.logging import logger
from dreampy3.utils.find_closest_indices import find_closest
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
#from scipy.interpolate import griddata
#from dreampy3.utils.two_gaussian_fit import two_gaussian_fit

logger.name = __name__

#mapping between slot numbers and band names
defaultBandNames = {
    5: 1,
    8: 3,
    11: 2,
    14: 4,
    17: 6,
    20: 5
    }

Frontend_LO_LOWFREQ = 93.4   #GHz
Frontend_LO_HIGHFREQ = 112.296 #GHz
Backend_LO_LOWFREQ = 15.4 #GHz
Backend_LO_HIGHFREQ = 21.7  #GHz

TOL = 1.e-16

class RedshiftScan(LMThdu):
    """Doc string here"""
    def __init__(self, data=None, header=None,
                 filename=None, groupname=None):
        """
        Redshift netcdf files can be netcdf4 which can
        have multiple groups. groupname if provided is the
        group name of the particular redshift scan
        """
        LMThdu.__init__(self, data=data, header=header,
                        filename=filename)
        self.windows = {}
        self.groupname = groupname
        #self.sigma = {}
        self.calibrated = False
        self.blank = blank
        
    def _get_chassis_str(self):
        return "RedshiftChassis_%d_" % (int(self.header.ChassisNumber))

    def _get_xfr1_deltasign(self, bandname):
        if bandname == 1:
            return ((self.fend_LO_LowFreq-self.bend_LO_HighFreq), 1)
        if bandname == 2:
            return((self.fend_LO_LowFreq-self.bend_LO_LowFreq), 1)
        if bandname == 3:
            return ((self.fend_LO_LowFreq), -1)
        if bandname == 4:
            return ((self.fend_LO_HighFreq-self.bend_LO_HighFreq), 1)
        if bandname == 5:
            return ((self.fend_LO_HighFreq-self.bend_LO_LowFreq), 1)
        if bandname == 6:
            return ((self.fend_LO_HighFreq), -1)
        
    def _get_freqs(self):
        chas_str = self._get_chassis_str()
        self.fend_LO_LowFreq = self.header.get('Rsfend.LoLowFreq', Frontend_LO_LOWFREQ)
        if abs(self.fend_LO_LowFreq - 0.0) < TOL:
            self.fend_LO_LowFreq = Frontend_LO_LOWFREQ
        self.fend_LO_HighFreq = self.header.get('Rsfend.LoHighFreq', Frontend_LO_HIGHFREQ)
        if abs(self.fend_LO_HighFreq - 0.0) < TOL:
            self.fend_LO_HighFreq = Frontend_LO_HIGHFREQ
        self.bend_LO_LowFreq = self.header.get('%s.LoLowFreq' % chas_str,
                                                Backend_LO_LOWFREQ)
        self.bend_LO_HighFreq = self.header.get('%s.LoHighFreq' % chas_str,
                                                Backend_LO_HIGHFREQ)
        self.basefreq = (0.0078125*2.5) + numpy.arange(256)*0.03125
        self.frequencies = numpy.ma.zeros((self.header.SlotNumber.size, 256),
                                       dtype='float')
        #bandnames = self.header.get('%s.BandName' % chas_str, None)
        bandnames = self.header.BandName
        for i, slotno in enumerate(self.header.SlotNumber):
            if bandnames is not None:
                bandname = bandnames[i]
            else:
                bandname = defaultBandNames[slotno]
            #bandname = self.header.get('%s.BandName' % chas_str,
            xfr1, deltasign = self._get_xfr1_deltasign(bandname)
            self.frequencies[i, :] = xfr1 + deltasign*self.basefreq

            
    def _get_frequencies(self):
        if hasattr(self, 'frequencies'):
            return self.frequencies
        elif hasattr(self.data, 'frequencies'):
            return self.data.frequencies
        else:
            return None

    def _get_spectrum(self):
        if hasattr(self, 'spectrum'):
            return self.spectrum
        elif hasattr(self.data, 'spectrum'):
            return self.data.spectrum
        else:
            return None        

    def _get_comp_freq(self):
        if hasattr(self, 'compfreq'):
            return self.compfreq
        elif hasattr(self.data, 'compfreq'):
            return self.data.compfreq
        else:
            return None

    def _get_comp_spectrum(self):
        if hasattr(self, 'compspectrum'):
            return self.compspectrum
        elif hasattr(self.data, 'compspectrum'):
            return self.data.compspectrum
        else:
            return None        

    def _get_tsys(self):
        if hasattr(self, 'Tsys'):
            return self.Tsys
        elif hasattr(self.header, 'Tsys'):
            return self.header.Tsys
        else:
            return None

    def _get_tint(self):
        if hasattr(self, 'tint'):
            return self.tint
        elif hasattr(self.header, 'tint'):
            return self.header.tint
        else:
            return None
        
    def _get_hdu(self, nc):
        if hasattr(nc, 'hdu'):
            return nc.hdu
        else:
            return nc

    def _get_sigma(self):
        if hasattr(self, 'sigma'):
            return self.sigma
        elif hasattr(self.header, 'sigma'):
            return self.header.sigma
        else:
            return None

    def _get_compsigma(self):
        if hasattr(self, 'compsigma'):
            return self.compsigma
        elif hasattr(self.header, 'compsigma'):
            return self.header.compsigma
        else:
            return self.header.get('Spectrum.compsigma', None)

    def _set_spectrum(self, spectrum):
        if hasattr(self, 'spectrum'):
            self.spectrum = spectrum
        elif hasattr(self.data, 'spectrum'):
            self.data.spectrum = spectrum

    def get_cal(self, corr_linear=True):
        if hasattr(self, 'calibrated') and self.calibrated:
            return
        #from dreampy.redshift.netcdf.calibration_scan import RedshiftCalibration
        from .calibration_scan import RedshiftCalibration
        logger.info("Getting cal object from Cal ObsNum %d" % self.header.CalObsNum)
        self._get_freqs()
        chassis = int(self.header.ChassisNumber)
        self.cal = RedshiftCalibration(self.header.CalObsNum,
                                       chassis,
                                       corr_linear=corr_linear)
        self.calibrated = True
        
    def process_uncalibrated_scan(self, acf='AccData', ps=False):
        """Wrapper for all kinds of complicated processing of
        data"""
        ccal = get_corr_cal_matrix(self.header)
        self._get_freqs()
        self.spectrum = numpy.ma.zeros(self.data.AccData.shape,
                                    dtype=self.data.AccData.dtype)
        chassis = int(self.header.ChassisNumber)
        for rpt in numpy.arange(self.data.AccData.shape[0]):
            for board, slotno in enumerate(self.header.SlotNumber):
                data = getattr(self.data, acf)
                normacf = data[rpt, board, :]/self.data.AccSamples[0, board]
                badlags = get_bad_lags(chassis, board)
                if len(badlags) > 0:
                    logger.info("Removing bad lags in chassis %d, board %d: %s" % (chassis, board, badlags))
                self.spectrum[rpt, board, :] = makespectrum(normacf, badlags=badlags, corr_cal=ccal[slotno])
        for board in self.header.BoardNumber:
            self.windows[board] = []
        self.calibrated = False
        self.ylabel = acf
        
    def process_calibrated_on_scan(self, calibrate=True,
                                   scalar_tsys=False,
                                   corr_linear=True):
        logger.info("Processing Calibrated ON Observation with ObsNum %d" % int(self.header.ObsNum))
        ccal = get_corr_cal_matrix(self.header)
        self._get_freqs()
        self.spectrum = numpy.ma.zeros(self.data.AccData.shape,
                                    dtype=self.data.AccData.dtype)
        self.offspec = numpy.ma.zeros(self.data.AccData.shape,
                                    dtype=self.data.AccData.dtype)        
        chassis = int(self.header.ChassisNumber)
        if calibrate and not self.header.CalObsNum:
            raise LMTRedshiftError("Redshift.process_calibrated_scan",
                                   "No CalObsNum to use in scan")
        if calibrate:
            #self.cal = RedshiftCalibration(self.header.CalObsNum,
            #                               chassis)
            self.get_cal(corr_linear=corr_linear)
        for rpt in numpy.arange(self.data.AccData.shape[0]):
            logger.info("Processing repeat %d" % (rpt+1))
            for board, slotno in enumerate(self.header.SlotNumber):
                onminusoff = self.data.AccData[rpt, board, :]/self.data.AccSamples[rpt, board]
                off = self.data.RefData[rpt, board, :]/self.data.AccSamples[rpt, board]
                badlags = get_bad_lags(chassis, board)
                if len(badlags) > 0:
                    logger.info("Removing bad lags in chassis %d, board %d: %s" % (chassis, board, badlags))
                if calibrate and self.cal and hasattr(self.cal, 'Tsys'):
                    onminusoffspec = makespectrum(self.data.AccData[rpt, board, :],
                                                  badlags=badlags,
                                                  g_norm=self.cal.g_norm[board, :],
                                                  #ratio=self.cal.skyratio[board, :],
                                                  corr_cal=ccal[slotno])
                    offspec = makespectrum(self.data.RefData[rpt, board, :],
                                           badlags=badlags,
                                           #ratio=self.cal.skyratio[board, :],
                                           #g_norm=self.cal.g_norm[board, :],
                                           corr_cal=ccal[slotno])
                    self.offspec[rpt, board, :] = offspec
                    self.spectrum[rpt, board, :] = (onminusoffspec/offspec) * self.cal.Tsys[board, :]
                    #self.spectrum[rpt, board, :] = (onminusoffspec/(offspec-self.cal.skyspec[board,:])) * self.cal.Tsys[board, :]
                    #self.spectrum[rpt, board, :] = (onminusoffspec/self.cal.skyspec[board,:]) * self.cal.Tsys[board, :]
                    self.calibrated = True
                    self.ylabel = 'TA* (K)'
                else:
                    onminusoffspec = makespectrum(self.data.AccData[rpt, board, :],
                                                  badlags=badlags,
                                                  corr_cal=ccal[slotno])
                    offspec = makespectrum(self.data.RefData[rpt, board, :],
                                           badlags=badlags,
                                           corr_cal=ccal[slotno])
                    self.offspec[rpt, board, :] = offspec
                    self.spectrum[rpt, board, :] = (onminusoffspec/offspec)
                    self.calibrated = False
                    self.ylabel = '(ON-OFF)/OFF (arb)'
        for board in self.header.BoardNumber:
            self.windows[board] = []

    def process_calibrated_acf_crossscan(self, calibrate=True,
                                         scalar_tsys=False,
                                         corr_linear=True):
        """
        ACF CrossScan: What the heck is it?
        SWTA=SWTB=0, polar demod = 2, blank = 0
        refacc = 0
        for now we will ignore calibration
        """
        from dreampy.redshift.netcdf.calibration_scan import RedshiftCalibration
        #ccal = get_corr_cal_matrix(self.header)
        #self._get_freqs()
        logger.info("Processing Calibrated ACF CrossScan Observation with ObsNum %d" % int(self.header.ObsNum))
        self.continuum = numpy.ma.zeros((self.data.AccData.shape[0], #dmps
                                      self.data.AccData.shape[1] #boards
                                      ),
                                    dtype=self.data.AccData.dtype)
        chassis = int(self.header.ChassisNumber)
        if not self.header.CalObsNum and calibrate:
            raise LMTRedshiftError("Redshift.process_calibrated_acf_crossscan",
                                   "No CalObsNum to use in scan")
        if calibrate:
            #self.cal = RedshiftCalibration(self.header.CalObsNum,
            #                               chassis)
            self.get_cal(corr_linear=corr_linear)
        for dmp in numpy.arange(self.data.AccData.shape[0]):
            for board, slotno in enumerate(self.header.SlotNumber):
                normacf = self.data.AccData[dmp, board, :]/self.data.RefData[dmp, board, :]
                #badlags = get_bad_lags(chassis, board)
                #if len(badlags) > 0:
                #logger.info("Removing bad lags in chassis %d, board %d: %s" % (chassis, board, badlags))
                normacf = apply_bad_lags(normacf, chassis, board)
                if calibrate and self.cal and hasattr(self.cal, 'Tsys'):
                    self.continuum[dmp, board] = nanmean(nanmean(normacf) * self.cal.Tsys[board, :])
                    self.calibrated = True
                    self.ylabel = 'CS ACF TA* (K)'
                else:
                    self.continuum[dmp, board] = nanmean(normacf)
                    self.calibrated = False
                    self.ylabel = 'CS ACF (arb)'

    def process_calibrated_spectral_crossscan(self, calibrate=True,
                                              scalar_tsys=False,
                                              corr_linear=True):
        """
        What the heck is spectral CrossScan?
        SWTA=SWTB=2, polar demod = 2, blank = 0
        refacc = 0
        """
        from dreampy.redshift.netcdf.calibration_scan import RedshiftCalibration
        logger.info("Processing Calibrated Spectral CrossScan Observation with ObsNum %d" % int(self.header.ObsNum))        
        ccal = get_corr_cal_matrix(self.header)
        self._get_freqs()
        self.continuum = numpy.ma.zeros((self.data.AccData.shape[0], #dmps
                                      self.data.AccData.shape[1] #boards
                                      ),
                                     dtype=self.data.AccData.dtype)
        chassis = int(self.header.ChassisNumber)
        if not self.header.CalObsNum and calibrate:
            raise LMTRedshiftError("Redshift.process_calibrated_acf_crossscan",
                                   "No CalObsNum to use in scan")
        if calibrate:
            #self.cal = RedshiftCalibration(self.header.CalObsNum,
            #                               chassis)
            self.get_cal(corr_linear=corr_linear)
        for dmp in numpy.arange(self.data.AccData.shape[0]):
            for board, slotno in enumerate(self.header.SlotNumber):
                #normacf = self.data.AccData[dmp, board, :]/self.data.RefData[dmp, board, :]
                badlags = get_bad_lags(chassis, board)
                if len(badlags) > 0:
                    logger.info("Removing bad lags in chassis %d, board %d: %s" % (chassis, board, badlags))
                if calibrate and self.cal and hasattr(self.cal, 'Tsys'):
                    onminusoffspec = makespectrum(self.data.AccData[dmp, board, :],
                                                  badlags=badlags,
                                                  g_norm=self.cal.g_norm[board, :],
                                                  corr_cal=ccal[slotno])
                    offspec = makespectrum(self.data.RefData[dmp, board, :],
                                           badlags=badlags,
                                           corr_cal=ccal[slotno])
                    self.continuum[dmp, board] = ((onminusoffspec/offspec) * self.cal.Tsys[board, :])[50:245].sum()
                    self.calibrated = True
                    self.ylabel = 'CS Spectral TA* (K)'
                else:
                    onminusoffspec = makespectrum(self.data.AccData[dmp, board, :],
                                                  badlags=badlags,
                                                  corr_cal=ccal[slotno])
                    offspec = makespectrum(self.data.RefData[dmp, board, :],
                                           badlags=badlags,
                                           corr_cal=ccal[slotno])                    
                    self.continuum[dmp, board] = ((onminusoffspec/offspec))[50:245].sum()
                    self.calibrated = False
                    self.ylabel = 'CS Spectral (arb)'

    def process_calibrated_compressed_continuum_crossscan(self, calibrate=True,
                                                          scalar_tsys=False,
                                                          numoffset=10,
                                                          gain=[0.88, 0.8, 1.0, 1.0],
                                                          apply_gain=True,
                                                          corr_linear=True):
        """
        ACF Compressed Continuum CrossScan: What the heck is it?
        SWTA=SWTB=0, polar demod = 2, blank = 0
        refacc = 0
        for now we will ignore calibration
        """
        from dreampy.redshift.netcdf.calibration_scan import RedshiftCalibration
        #ccal = get_corr_cal_matrix(self.header)        #self._get_freqs()
        logger.info("Processing Compressed Continuum ACF CrossScan Observation with ObsNum %d" % int(self.header.ObsNum))
        self.continuum = numpy.ma.zeros(self.data.AccAverage.shape,
                                     dtype=self.data.AccAverage.dtype)
        chassis = int(self.header.ChassisNumber)
        if not self.header.CalObsNum and calibrate:
            raise LMTRedshiftError("Redshift.process_calibrated_acf_crossscan",
                                   "No CalObsNum to use in scan")
        if calibrate:
            #self.cal = RedshiftCalibration(self.header.CalObsNum,
            #                               chassis)
            self.get_cal(corr_linear=corr_linear)
        for dmp in numpy.arange(self.data.AccAverage.shape[0]):
            for board, slotno in enumerate(self.header.SlotNumber):
                normacf = self.data.AccAverage[dmp, board]
                if calibrate and self.cal and hasattr(self.cal, 'Tsys'):
                    self.continuum[dmp, board] = (normacf * self.cal.Tsys[board, :])
                    self.calibrated = True
                    self.ylabel = 'CS ACF TA* (K)'
                else:
                    self.continuum[dmp, board] = normacf
                    self.calibrated = False
                    self.ylabel = 'CS ACF (arb)'
        for board, slotno in enumerate(self.header.SlotNumber):
            offset = self.continuum[:, board][:numoffset].mean()
            self.continuum[:, board] -= offset
            if chassis in (1, 2):
                self.continuum[:, board] = -self.continuum[:, board]
            if apply_gain:
                self.continuum[:, board] *= gain[chassis]
            


    def process_calibrated_crossscan(self, calibrate=True, scalar_tsys=False,
                                     corr_linear=True):
        #CrossScanMode = int(self.header.get('CrossScan.CrossScanMode', -1))
        chas_str = self._get_chassis_str()

        if hasattr(self.data, 'AccAverage'):
            #this is the new compressed continuum data set
            self.process_calibrated_compressed_continuum_crossscan(calibrate=calibrate,
                                                                  scalar_tsys=scalar_tsys)
            return

        self.CrossScanMode = int(self.header.get('%s.SwitchSettings' % chas_str, -1))

        if self.CrossScanMode == -1:
            raise LMTRedshiftError('process_calibrated_crossscan', "No CrossScanMode in header")
        if self.CrossScanMode == 0:
            #Acf Continuum Mode
            self.process_calibrated_acf_crossscan(calibrate=calibrate,
                                                  scalar_tsys=scalar_tsys)
        elif self.CrossScanMode == 1:
            #Spectral Continuum Mode
            self.process_calibrated_spectral_crossscan(calibrate=calibrate,
                                                       scalar_tsys=scalar_tsys)
        else:
            raise LMTRedshiftError('process_calibrated_crossscan', "No such CrossScanMode %d" % self.CrossScanMode)

    def process_calibrated_bs_scan(self, calibrate=True, ps=False,
                                   scalar_tsys=False,
                                   corr_linear=True):
        if self.header.AccumMode == 'AccumCompressed':
            raise LMTRedshiftError('process_calibrated_bs_scan', "AccumMode was wrongly set to AccumCompressed. Put it to AccumNormal")
        if ps:
          logger.info("Processing Calibrated Ps Observation with ObsNum %d" % int(self.header.ObsNum))
        else:
            logger.info("Processing Calibrated Bs Observation with ObsNum %d" % int(self.header.ObsNum))
        chas_str = self._get_chassis_str()
        ccal = get_corr_cal_matrix(self.header)
        self._get_freqs()
        nrpt2, nbrds, nchan = self.data.AccData.shape
        self.spectrum = numpy.ma.zeros((int(nrpt2/2), nbrds, nchan),
                                    dtype=self.data.AccData.dtype)
        #self.rawspectrum = numpy.ma.zeros((nrpt2, nbrds, nchan),
        #                            dtype=self.data.AccData.dtype)
        #self.rawrefspectrum = numpy.ma.zeros((nrpt2, nbrds, nchan),
        #                                  dtype=self.data.AccData.dtype)        
        if ps:
            self.tint = 2.*int(self.header.get('Ps.NumRepeats'))*float(self.header.get('Ps.TSamp'))
        else:
            self.tint = 2.*int(self.header.get('Bs.NumRepeats'))*float(self.header.get('Bs.TSamp'))
        chassis = int(self.header.ChassisNumber)
        if calibrate and not self.header.CalObsNum:
            raise LMTRedshiftError("Redshift.process_calibrated_scan",
                                   "No CalObsNum to use in scan")
        if calibrate:
            #self.cal = RedshiftCalibration(self.header.CalObsNum,
            #                               chassis)
            self.get_cal(corr_linear=corr_linear)
            self.Tsys = scipynanmean(self.cal.Tsys, axis=1)
        for i, rpt in enumerate(range(0, nrpt2, 2)):
            logger.info("Processing repeat %d" % (rpt+1))
            for board, slotno in enumerate(self.header.SlotNumber):
                if self.data.BufPos[rpt] == 0:
                    pind = rpt
                    if self.data.BufPos[rpt+1] == 1:
                        nind = rpt+1
                    else:
                        raise LMTRedshiftError("Redshift Bs", "error in BufPos indices")
                if self.data.BufPos[rpt] == 1:
                    nind = rpt
                    if self.data.BufPos[rpt+1] == 0:
                        pind = rpt+1
                    else:
                        raise LMTRedshiftError("Redshift Bs", "error in BufPos indices")
                pnorm = numpy.outer(self.data.AccSamples[pind], numpy.ones(256))
                nnorm = numpy.outer(self.data.AccSamples[nind], numpy.ones(256))
                pos_beam = self.data.AccData[pind, board, :]/pnorm[board]
                #pos_tot = 0.5*self.data.RefData[pind, board, :]/pnorm[board]
                pos_tot = self.data.RefData[pind, board, :]/pnorm[board]
                neg_beam = self.data.AccData[nind, board, :]/nnorm[board]
                #neg_tot = 0.5*self.data.RefData[nind, board, :]/nnorm[board]
                neg_tot = self.data.RefData[nind, board, :]/nnorm[board]
                badlags = get_bad_lags(chassis, board)
                if len(badlags) > 0:
                    logger.info("Removing bad lags in chassis %d, board %d: %s" % (chassis, board, badlags))
                if calibrate and self.cal and hasattr(self.cal, 'Tsys'):
                    posbeamspec = makespectrum(pos_beam,
                                               badlags=badlags,
                                               g_norm=self.cal.g_norm[board, :],
                                               #ratio=self.cal.skyratio[board, :],
                                               corr_cal=ccal[slotno])
                    postotspec = makespectrum(pos_tot,
                                              badlags=badlags,
                                              corr_cal=ccal[slotno])
                    negbeamspec = makespectrum(neg_beam,
                                               badlags=badlags,
                                               #ratio=self.cal.skyratio[board, :],
                                               g_norm=self.cal.g_norm[board, :],
                                               corr_cal=ccal[slotno])
                    negtotspec = makespectrum(neg_tot,
                                              badlags=badlags,
                                              corr_cal=ccal[slotno])
                    #print(negtotspec)
                    if chassis in (1,2):
                        beam_sign = -1.0
                    else:
                        beam_sign = 1.0
                    if scalar_tsys:
                        self.spectrum[i, board, :] = beam_sign*((2*posbeamspec/(postotspec-posbeamspec))\
                                                                - (2*negbeamspec/(negtotspec-negbeamspec))) * self.Tsys[board]/2.  #nanmean(self.cal.Tsys[board, :])/2.
                    else:
                        self.spectrum[i, board, :] = beam_sign*((2*posbeamspec/(postotspec-posbeamspec))\
                                                                - (2*negbeamspec/(negtotspec-negbeamspec))) * self.cal.Tsys[board, :]/2.
                        #print(posbeamspec[:4], postotspec[:4], negbeamspec[:4], negtotspec[:4])
                        #print(type(posbeamspec), type(self.spectrum[i, board, :]))
                        #print(self.spectrum[i, board, :])
                    #self.rawspectrum[rpt, board, :] = posbeamspec * self.cal.Tsys[board, :]
                    #self.rawspectrum[rpt+1, board, :] = negbeamspec * self.cal.Tsys[board, :]
                    #self.rawrefspectrum[rpt, board, :] = postotspec
                    #self.rawrefspectrum[rpt+1, board, :] = negtotspec
                    self.calibrated = True
                    self.ylabel = 'TA* (K)'
                else:
                    posbeamspec = makespectrum(pos_beam,
                                               badlags=badlags,
                                               corr_cal=ccal[slotno])
                    postotspec = makespectrum(pos_tot,
                                              badlags=badlags,
                                              corr_cal=ccal[slotno])
                    negbeamspec = makespectrum(neg_beam,
                                               badlags=badlags,
                                               corr_cal=ccal[slotno])
                    negtotspec = makespectrum(neg_tot,
                                              badlags=badlags,
                                              corr_cal=ccal[slotno])                    
                    if chassis in (1,2):
                        beam_sign = -1
                    else:
                        beam_sign = 1
                    self.spectrum[rpt, board, :] = 0.5*beam_sign*((posbeamspec/postotspec) \
                                                              -(negbeamspec/negtotspec))
                    self.calibrated = False
                    self.ylabel = '(ON-OFF)/OFF (arb)'
        for board in self.header.BoardNumber:
            self.windows[board] = []

    def process_calibrated_map(self, calibrate=True, scalar_tsys=False,
                               remove_offset=True, numoff=20,
                               gain=[0.88, 0.8, 1.0, 1.0],
                               apply_gain=True, remove_spikes=True,
                               **kwargs):
        """
        processes Map Obspgm that has been made using compressed
        continuum mode. Uses a regridding algorithm
        and uses some kwargs arguments to derive output
        grid size and sampling
        """
        logger.info("Processing Compressed Continuum Map data and regridding for Observation with ObsNum %d" % int(self.header.ObsNum))
        xlength = numpy.degrees(self.header.get('Map.XLength'))*3600.0
        ylength = numpy.degrees(self.header.get('Map.YLength'))*3600.0
        ind = numpy.where(self.data.BufPos == 0)
        xpos = numpy.degrees(self.data.XPos[ind])*3600.
        ypos = numpy.degrees(self.data.YPos[ind])*3600.
        if int(self.header.ChassisNumber) in (1, 2):
            z = -self.data.AccAverage[ind[0], :]
        else:
            z = self.data.AccAverage[ind[0], :]
        if apply_gain:
            z *= gain[int(self.header.ChassisNumber)]

        sz_before = z.copy().size
        #check to see if AccSamples is available and if so only
        #use values with 3936 AccSamples
        if remove_spikes and hasattr(self.data, 'AccSamples'):
            ind = numpy.where(self.data.AccSamples.mean(axis=1) == 3936)
            xpos, ypos = xpos[ind], ypos[ind]
            z = z[ind[0], :]
            num_spikes = sz_before - z.size
            if num_spikes > 0:
                print("Removed %d number out of %d data points that failed AccSamples criteria" % (num_spikes, sz_before))
        if remove_offset:
            for board in range(6):
                z[:, board] = z[:, board] - z[:, board][:numoff].mean()
        ramp = kwargs.get('ramp', 5.)
        numpoints = kwargs.get('numpoints', 100)
        numypoints = kwargs.get('numypoints', 100)
        xlength = xlength * (1.-ramp/100.)
        ylength = ylength * (1.-ramp/100.)
        ind = numpy.logical_and(xpos > -xlength/2., xpos < xlength/2.)
        xpos, ypos = xpos[ind], ypos[ind]
        z = z[ind]
        ind = numpy.logical_and(ypos > -ylength/2., ypos < ylength/2.)
        xpos, ypos = xpos[ind], ypos[ind]
        z = z[ind]
        self.xi = numpy.linspace(-xlength/2, xlength/2, numpoints)
        self.yi = numpy.linspace(-ylength/2, ylength/2, numypoints)
        print("Making %d x %d map" % (numpoints, numypoints))
        self.BeamMap = numpy.ma.zeros((self.yi.size, self.xi.size,
                                    self.data.AccAverage.shape[1]),
                                   dtype='float')
        for i in range(self.data.AccAverage.shape[1]):
            self.BeamMap[:, :, i] = griddata(xpos, ypos, z[:, i], self.xi, self.yi)

    def process_calibrated_map_old(self, calibrate=True, scalar_tsys=False,
                                   **kwargs):
        """
        processes Map Obspgm that has been made using compressed
        continuum mode. Uses a regridding algorithm
        and uses some kwargs arguments to derive output
        grid size and sampling
        """
        logger.info("Processing Compressed Continuum Map data and regridding for Observation with ObsNum %d" % int(self.header.ObsNum))
        xlength = numpy.degrees(self.header.get('Map.XLength'))*3600.0
        ylength = numpy.degrees(self.header.get('Map.YLength'))*3600.0
        hpbw = numpy.degrees(self.header.get('Map.HPBW', default=0.0072))*3600.0
        beam_sample = kwargs.get('beam_sample', 10.)
        stepsize = hpbw/beam_sample
        self.num_grid = int(xlength/stepsize)
        # make sure it is odd
        if self.num_grid == 2 * (self.num_grid/2):
            self.num_grid += 1
        logger.info("Regridding into output grid of %d x %d points" % (self.num_grid, self.num_grid))
        logger.info("Regridding to output size of %.2f to %.2f arcseconds" % (-xlength/2., xlength/2.))
        self.grid = numpy.linspace(-xlength/2., xlength/2.,
                                   self.num_grid)
        self.numpoints = numpy.ma.zeros((self.num_grid, self.num_grid),
                                     dtype='int')
        self.BeamMap = numpy.ma.zeros((self.num_grid, self.num_grid,
                                    self.data.AccAverage.shape[1]),
                                   dtype='float')
        ind = numpy.where(self.data.BufPos == 0)
        xpos = numpy.degrees(self.data.XPos[ind])*3600.
        ypos = numpy.degrees(self.data.YPos[ind])*3600.
        z = self.data.AccAverage[ind[0], :]

        #ind is indices that will fall in output grid
        ind = numpy.logical_and(xpos >= self.grid[0], xpos <= self.grid[-1])
        xpos, ypos = xpos[ind], ypos[ind]
        z = z[ind, :]
        xind = find_closest(self.grid, xpos)
        yind = find_closest(self.grid, ypos)
        for i in range(xind.size):
            self.numpoints[yind[i], xind[i]] += 1
            self.BeamMap[yind[i], xind[i], :] += z[i, :]
        ind = numpy.where(self.numpoints != 0)
        #normalize beam_map by average number of points
        for i in range(self.BeamMap.shape[2]):
            bm = self.BeamMap[:, :, i]
            bm[ind] = bm[ind]/self.numpoints[ind]
            self.BeamMap[:, :, i] = bm
        logger.info("Output beam data in BeamMap with shape (%d, %d, %d)" % (self.BeamMap.shape[0], self.BeamMap.shape[1], self.BeamMap.shape[2]))
        
    def process_scan(self, calibrate=True, raw=False, rawdata='AccData',
                     scalar_tsys=False, corr_linear=True, **kwargs):
        """wrapper for processing any redshift scan"""
        if self.header.ObsPgm == 'On':
            if not raw:
                self.process_calibrated_on_scan(calibrate=calibrate, scalar_tsys=scalar_tsys, corr_linear=corr_linear)
            else:
                self.process_uncalibrated_scan(acf=rawdata)
        elif self.header.ObsPgm == 'Bs':
            #Beam-switched observation
            if not raw:
                self.process_calibrated_bs_scan(calibrate=calibrate, scalar_tsys=scalar_tsys, corr_linear=corr_linear)
            else:
                self.process_uncalibrated_scan(acf=rawdata)
        elif self.header.ObsPgm == 'Ps' and int(self.header.get('Ps.RefSwitch')) == 2:
            #exactly analogous to Beam-switched observations
            if not raw:
                self.process_calibrated_bs_scan(calibrate=calibrate, ps=True, scalar_tsys=scalar_tsys, corr_linear=corr_linear)
            else:
                self.process_uncalibrated_scan(acf=rawdata)            
        elif self.header.ObsPgm == 'CrossScan':
            if not raw:
                self.process_calibrated_crossscan(calibrate=calibrate, scalar_tsys=scalar_tsys, corr_linear=corr_linear)
            else:
                self.process_uncalibrated_scan(acf=rawdata)
        elif self.header.ObsPgm == 'Map':
            if not raw:
                self.process_calibrated_map(calibrate=calibrate, scalar_tsys=scalar_tsys,
                                            **kwargs)
        else:
            logger.warning("Processing unknown data type as simply uncalibrated scan. ObsPgm = %s" % self.header.ObsPgm)
            self.process_uncalibrated_scan(acf=rawdata)
        
    def baseline(self, windows={}, order=0, subtract=False):
        """applies a baseline over given windows. windows is a
        dictionary of list of two-tuples when given. If subtract=True,
        remove the baseline from the data"""
        frequencies = self._get_frequencies().filled(fill_value=numpy.nan)
        blankind = get_blank_indices_array(self._get_spectrum())
        spectrum = blank_to_nan(self._get_spectrum().filled(fill_value=numpy.nan))
        if not windows:
            if self.windows:
                for board, wins in self.windows.items():
                    if not wins:
                        c1, c2 = frequencies[board][0], frequencies[board][-1]
                        self.windows[board] = [(c1, c2)]
                windows = self.windows
            else:
                windows = {}
                for i in self.header.BoardNumber:
                    windows[i] = [(frequencies[i][0], frequencies[i][-1])]
                self.windows = windows
                logger.info("Computing automatic windows: %s" % self.windows)
        else:
            self.windows = windows
        if not isinstance(self.windows, dict):
            raise LMTRedshiftError("RedshiftScan.baseline", "windows should be a dictionary type")
        if order<0 or order>3:
            raise LMTRedshiftError("RedshiftScan.baseline", "Order should be 0,1,2 or 3")
        self.order = order
        self.sigma = numpy.zeros((spectrum.shape[0], spectrum.shape[1]),
                                 dtype=spectrum.dtype)
        #self.sigma = {}
        for rpt in numpy.arange(spectrum.shape[0]):
            for board in self.header.BoardNumber:
                ct = 0
                for win in windows[board]:  
                    if isinstance(win, list) or isinstance(win, tuple) and len(win) == 2:
                        c1, c2 = win
                        c1,c2 = sorted([c1,c2])
                    else:
                        c1, c2 = frequencies[board][0], frequencies[board][-1]
                    ind = numpy.logical_and(frequencies[board] >= c1,
                                            frequencies[board] <= c2)
                    if ct == 0:
                        finalind = numpy.logical_or(ind,ind)
                    else:
                        finalind = numpy.logical_or(finalind,ind)
                    ct += 1 
                finalind = numpy.logical_and(finalind,
                                             numpy.isfinite(spectrum[rpt, board, :]))
                if finalind.any():
                    # only do a fit if there is any part of spectrum that is 
                    # not blanked out
                    ind = numpy.where(finalind)
                    p = numpy.polyfit(frequencies[board][ind],
                                      spectrum[rpt, board, :][ind],
                                      order)
                    spectrum[rpt, board, :] = spectrum[rpt, board, :] - \
                        numpy.polyval(p, frequencies[board])
                    sig = spectrum[rpt, board, :][ind].std()
                else:
                    # since there are no channels available to fit out
                    # baseline, sigma is effectively put to a large value
                    sig = 1e38
                # if rpt == 0:
                #     self.sigma[board] = numpy.array(sig)
                #     self.sigma[board].shape = (1,1)
                # else:
                #     self.sigma[board] = numpy.column_stack((self.sigma[board],
                #                                             numpy.array(sig)))
                self.sigma[rpt, board] = sig
                if not subtract and finalind.any():
                    spectrum[rpt, board, :] = spectrum[rpt, board, :] + \
                                              numpy.polyval(p, frequencies[board])
                logger.info("Chassis %d, Board %d, rpt %d, sigma = %.3f" % \
                            (int(self.header.ChassisNumber), board, rpt, sig))
        spectrum[blankind] = blank
        #self.spectrum = spectrum
        self.spectrum = numpy.ma.masked_array(spectrum, mask=numpy.isnan(spectrum))
        
    def baseline_compspectrum(self, windows=[], order=0, subtract=False):
        """applies a baseline over given windows. windows is a
        list of two-tuples when given. If subtract=True,
        remove the baseline from the data"""
        compfreq = self._get_comp_freq()
        blankind = get_blank_indices(self._get_comp_spectrum())
        compspectrum = blank_to_nan(self._get_comp_spectrum())
        if not windows:
            if hasattr(self, 'compwindows') and self.compwindows:
                windows = self.compwindows
            else:
                c1, c2 = compfreq[0], compfreq[-1]
                windows = [(c1, c2)]
                logger.info("Using automatic windows")
        if isinstance(windows, list) or isinstance(windows, tuple):
            raise LMTRedshiftError("RedshiftScan.baseline", "windows should be a list or tuple type")
        if order<0 or order>3:
            raise LMTRedshiftError("RedshiftScan.baseline_compspectrum", "Order should be 0,1,2 or 3")
        self.order = order
        self.compsigma = numpy.ma.zeros((compspectrum.shape[0], ),
                                     dtype=compspectrum.dtype)
        self.compwindows = windows
        #self.sigma = {}
        for rpt in numpy.arange(compspectrum.shape[0]):
            ct = 0
            for win in windows:  
                c1, c2 = win
                c1, c2 = sorted([c1,c2])
                ind = numpy.logical_and(compfreq >= c1,
                                        compfreq <= c2)
                if ct == 0:
                    finalind = numpy.logical_or(ind,ind)
                else:
                    finalind = numpy.logical_or(finalind,ind)
                ct += 1
            #print finalind.shape, finalind
            #print rpt
            #print numpy.isfinite(compspectrum[rpt, :]).shape
            finalind = numpy.logical_and(finalind,
                                         numpy.isfinite(compspectrum[rpt, :]))
            ind = numpy.where(finalind)
            p = numpy.polyfit(compfreq[ind],
                              compspectrum[rpt, :][ind],
                              order)
            compspectrum[rpt, :] = compspectrum[rpt, :] - \
                                   numpy.polyval(p, compfreq)
            sig = compspectrum[rpt, :][ind].std()
            # if rpt == 0:
            #     self.sigma[board] = numpy.array(sig)
            #     self.sigma[board].shape = (1,1)
            # else:
            #     self.sigma[board] = numpy.column_stack((self.sigma[board],
            #                                             numpy.array(sig)))
            self.compsigma[rpt] = sig
            if not subtract:
                compspectrum[rpt, :] = compspectrum[rpt, :] + \
                                       numpy.polyval(p, compfreq)
            logger.info("Chassis %d, rpt %d, sigma = %.3f" % \
                        (int(self.header.ChassisNumber), rpt, sig))
        compspectrum[blankind] = blank
        self.compspectrum = compspectrum

    def make_composite_scan(self):
        frequencies = self._get_frequencies()
        spectrum = blank_to_nan(self._get_spectrum())
        #tsys = self._get_tsys()
        if frequencies is None or spectrum is None:
            raise LMTRedshiftError('make_composite_scan', "ACF data here. First run process_scan on the scan first")
        fmin, fmax = frequencies.min(), frequencies.max()
        self.nchan = int(abs((fmax-fmin)/0.03125))
        self.compfreq = fmin + numpy.arange(self.nchan)*0.03125
        self.compspectrum = numpy.zeros((spectrum.shape[0], self.nchan), dtype='float')
        #if tsys is not None:
        #    self.Tsys = tsys.mean()
        wt = numpy.ma.zeros((spectrum.shape[0], self.nchan), dtype='float')
        for rpt in numpy.arange(spectrum.shape[0]):       
            for board in range(frequencies.shape[0]):
                xd = frequencies[board, :].copy()
                diff = numpy.abs(self.compfreq - xd.min()) 
                x1 = diff.argmin()
                x2 = x1 + 256
                if x2 > self.nchan-1:
                    x2 = self.nchan
                sc_fmin = self.compfreq[x1]
                rebin_freq = sc_fmin + numpy.arange(x2-x1)*0.03125
                yd = spectrum[rpt, board, :].copy()
                if xd[0] > xd[-1]:
                    xd = numpy.flipud(xd)
                    yd = numpy.flipud(yd)
                interp = interp1d(xd, yd.filled(fill_value=numpy.nan), bounds_error=False)
                dt = interp(rebin_freq)
                weights = numpy.ma.ones(x2-x1)
                idx = numpy.isnan(dt)
                dt[idx] = 0.0
                weights[idx] = 0.0
                self.compspectrum[rpt, x1:x2] += dt
                wt[rpt, x1:x2] += weights
        for rpt in numpy.arange(spectrum.shape[0]):
            self.compspectrum[rpt, :] = self.compspectrum[rpt, :]/wt[rpt, :]

    def write_composite_scan_to_ascii(self, filename):
        if not hasattr(self, 'compspectrum'):
            raise LMTRedshiftError('write_composite_scan_to_ascii', "ACF data here. First run make_composite_scan on the scan first")
        x = numpy.column_stack((self.compfreq, self.compspectrum.mean(axis=0)))
        logger.info("Saving Composite Spectrum to file %s" % filename)
        numpy.savetxt(filename, x)

    def write_spectra_to_ascii(self, filename):
        if not hasattr(self, 'spectrum'):
            raise LMTRedshiftError('write_spectra_to_ascii', "ACF data here. First run process_scan on the scan first")
        for board in range(self.frequencies.shape[0]):        
            if board == 0:
                x = numpy.column_stack((self.frequencies[board,:], self.spectrum.mean(axis=0)[board, :]))
            else:
                x = numpy.column_stack((x, self.frequencies[board,:]))
                x = numpy.column_stack((x, self.spectrum.mean(axis=0)[board, :]))
        logger.info("Saving spectra to file %s" % filename)
        numpy.savetxt(filename, x)
    
    def make_ps_scan(self, refnc):
        """Given a reference netcdf file, this method
        uses the current netcdf file data and the reference data to
        construct a traditional position-switched spectra"""
        frequencies = self._get_frequencies()
        spectrum = self._get_spectrum()
        refnc.hdu.process_scan()
        for rpt in numpy.arange(spectrum.shape[0]):       
            for board in range(frequencies.shape[0]):
                spectrum[rpt, board, :] = self.cal.Tsys[board,:] * (self.offspec[rpt, board, :] - \
                                           refnc.hdu.offspec[rpt, board, :])/ \
                                           refnc.hdu.offspec[rpt, board, :]
    

    def make_dbs_scan(self, negnc):
        """Given a negative beam netcdf instance, this method will use the
        current netcdf file data as the positive beam and then
        construct a Dual-beam-switched data set"""
        if not hasattr(self, 'frequencies'):
            self.process_scan()
        frequencies = self._get_frequencies()
        spectrum = self._get_spectrum()
        negnc.hdu.process_scan()
        for rpt in numpy.arange(spectrum.shape[0]):
            for board in range(frequencies.shape[0]):
                #spectrum[rpt, board, :] = self.cal.Tsys[board, :] * (self.spectrum[rpt, board, :] - negnc.hdu.spectrum[rpt, board, :])
                spectrum[rpt, board, :] = (self.spectrum[rpt, board, :] - negnc.hdu.spectrum[rpt, board, :])


    def smooth(self, nchan=2):
        """
        boxcar average the data given the number of channel
        """
        spectrum = self._get_spectrum()
        #spec = spectrum.copy()
        frequencies = self._get_frequencies()
        if nchan < 2:
            raise LMTRedshiftError("smooth", "at least 2 channel smooth")
        if nchan > 256:
            raise LMTRedshiftError("smooth", "has to be smaller than 256")
        for rpt in range(spectrum.shape[0]):
            for board in range(frequencies.shape[0]):
                #xmin, xmax = frequencies[board].min(), frequencies[board].max()
                spectrum[rpt][board] = numpy.convolve(spectrum[rpt][board],
                                                      numpy.ones(nchan)/nchan, 'same')

        
    def smooth_windows(self, window_len=10, window='hanning'):
        """
        smooth the data using a window with requested size.
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        input:
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming',
        'bartlett', 'blackman'
        flat window will produce a moving average smoothing.
        
        output:
        the smoothed signal

        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
    
        """
        if window_len<3:
            return
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise LMTRedshiftError("smooth", "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
        spectrum = self._get_spectrum()
        #spec = spectrum.copy()
        frequencies = self._get_frequencies()
        for rpt in xrange(spectrum.shape[0]):
            for board in xrange(frequencies.shape[0]):
                #xmin, xmax = frequencies[board].min(), frequencies[board].max()
                if spectrum[rpt][board].size < window_len:
                    raise LMTRedshiftError("smooth", "data vector needs to be bigger than window size.")
                x = spectrum[rpt][board]
                s = numpy.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
                if window == 'flat': #moving average
                    w = numpy.ones(window_len,'d')
                else:
                    w = eval('numpy.'+window+'(window_len)')
                y = numpy.convolve(w/w.sum(), s, mode='same')
                spectrum[rpt][board] = y[window_len-1:-window_len+1]
                    
    def _calc_weight(self, rpt, board, weight='sigma'):
        if weight not in ('sigma', 'equal', 'tsys'):
            raise LMTRedshiftError('_calc_weight', "weight should be one of 'sigma', 'equal', 'tsys'")
        if weight == 'sigma':
            if hasattr(self.header, 'sigma') or hasattr(self, 'sigma'):
                weight = 'sigma'
            else:
                logger.warning("No Baselines have been calculated. Defaulting to equal weighting")
                weight = 'equal'
        if weight == 'sigma':
            sigma = self._get_sigma()
            logger.info("Weighting by sigma %s for board %d" % (sigma[rpt, board], board))            
            #return 1./self.sigma[rpt, board]**2.
            return 1./sigma[rpt, board]**2.
        elif weight == 'tsys':
            tsys = self._get_tsys()
            tint = self._get_tint()
            if tsys is None or tint is None:
                logger.warning("Weighting by tsys will not work as tint or tsys is not available. Resetting to equal weighting")
                return 1.0
            return tint/tsys[board]**2.
        else:
            return 1.0

    def _calc_weight_comp(self, rpt, weight='sigma'):
        if weight not in ('sigma', 'equal', 'tsys'):
            raise LMTRedshiftError('_calc_weight_comp', "weight should be one of 'sigma', 'equal', 'tsys'")
        if weight == 'sigma':
            if hasattr(self.header, 'compsigma') or hasattr(self, 'compsigma'):
                weight = 'sigma'
            else:
                logger.warning("No Baselines have been calculated. Defaulting to equal weighting")
                weight = 'equal'
        if weight == 'sigma':
            sigma = self._get_compsigma()
            if isinstance(sigma, numpy.ndarray):
                logger.info("Weighting by sigma %s for rpt %d" % (sigma[rpt], rpt))
                #return 1./self.sigma[rpt, board]**2.
                return 1./sigma[rpt]**2.
            else:
                logger.info("Weighting by sigma %s for rpt %d" % (sigma, rpt))
                #return 1./self.sigma[rpt, board]**2.
                return 1./sigma**2.
        elif weight == 'tsys':
            tsys = self._get_tsys()
            tint = self._get_tint()
            if tsys is None or tint is None:
                logger.warning("Weighting by tsys will not work as tint or tsys is not available. Resetting to equal weighting")
                return 1.0
            return tint/tsys.mean()**2.
        else:
            return 1.0
        
    def average_all_repeats(self, weight='sigma', threshold_sigma=None):
        if not hasattr(self, 'spectrum'):
            raise LMTRedshiftError('average_all_repeats', "ACF data here. First run process_scan on the scan first")
        if not hasattr(self, 'sigma'):
            weight='equal'
        frequencies = self._get_frequencies()            
        blankind = get_blank_indices_array(self._get_spectrum())
        spectrum = blank_to_zero(self._get_spectrum())
        spec = spectrum.copy()
        self.spectrum = numpy.ma.zeros((1, spec.shape[1], spec.shape[2]),
                                    dtype=spec.dtype)

        for board in range(frequencies.shape[0]):
            weights = 0.0
            for rpt in range(spec.shape[0]):
                wt = self._calc_weight(rpt, board, weight=weight)
                if weight == 'sigma' and threshold_sigma is not None:
                    if 1/numpy.sqrt(wt) < threshold_sigma:
                        self.spectrum[0, board, :] += spec[rpt, board, :] * wt
                        weights += wt
                    else:
                        logger.info("Rejecting Board %s repeat %s due to sigma threshold" % (board, rpt))
                else:
                    if not spec.mask[rpt,board,:].all():
                        self.spectrum[0, board, :] += spec[rpt, board, :] * wt
                        weights += wt     

            if weights != 0.0:
                self.spectrum[0, board, :] /= weights
            else:
                self.spectrum[0, board, :] = numpy.nan
            self.spectrum[0,board,:]=numpy.ma.masked_invalid(self.spectrum[0,board,:])

        #blankarr = blankind.mean(axis=0).astype(numpy.bool)
        blankarr = blankind.mean(axis=0).astype(bool)
        blankarr.shape = (1, blankarr.shape[0], blankarr.shape[1])
        self.spectrum[blankarr] = self.blank
        self.baseline()

    def average_all_repeats_compspectrum(self, weight='sigma',
                                         threshold_sigma=None):
        if not hasattr(self, 'compspectrum'):
            raise LMTRedshiftError('average_all_repeats_compspectrum', "No composite spectrum available. First run make_composite_scan on the scan first")
        if not hasattr(self, 'compsigma'):
            weight='equal'
        compfreq = self._get_comp_freq()
        compspectrum = blank_to_nan(self._get_comp_spectrum())
        spec = compspectrum.copy()
        self.compspectrum = numpy.ma.zeros((1, spec.shape[1]),
                                        dtype=spec.dtype)

        weights = 0.0
        for rpt in range(spec.shape[0]):
            wt = self._calc_weight_comp(rpt, weight=weight)
            if weight == 'sigma' and threshold_sigma is not None:
                if 1/numpy.sqrt(wt) < threshold_sigma:
                    self.compspectrum[0, :] += spec[rpt, :] * wt
                    weights += wt
                else:
                    logger.info("Rejecting repeat %s due to sigma threshoalse, False], fill_value=999999)ld" % rpt)
            else:
                self.compspectrum[0, :] += spec[rpt, :] * wt
                weights += wt                
        self.compspectrum[0, :] /= weights
        self.baseline_compspectrum()
            
    def average_scans(self, nclist, weight='sigma', threshold_sigma=None):
        """Given a list of additional nc instances as a list
        averages current spectrum with the list of other ncs
        """
        frequencies = self._get_frequencies()
        #spectrum = numpy.nan_to_num(self._get_spectrum())
        blankind = get_blank_indices_array(self._get_spectrum())
        spectrum = blank_to_zero(self._get_spectrum())
        spec = spectrum.copy()
        avg_spectrum = numpy.ma.zeros((1, spec.shape[1], spec.shape[2]),
                                    dtype=spec.dtype)
        #self.spectrum = numpy.ma.zeros((1, spec.shape[1], spec.shape[2]),
        #                            dtype=spec.dtype)

        weights = {}
        for board in range(frequencies.shape[0]):
            weights[board] = numpy.zeros(frequencies.shape[1])
        nclist.insert(0, self._get_hdu(self))
        for nc in nclist:
            hdu = self._get_hdu(nc)
            #spec = numpy.nan_to_num(hdu._get_spectrum())
            spec = blank_to_zero(hdu._get_spectrum())
            blank_nind = numpy.logical_not(
                         get_blank_indices_array(hdu._get_spectrum()))
            for board in range(frequencies.shape[0]):
                wt = hdu._calc_weight(0, board, weight=weight)
                if weight == 'sigma' and threshold_sigma is not None:
                    print(nc, board, wt, threshold_sigma)
                    if 1/numpy.sqrt(wt) < threshold_sigma:
                        avg_spectrum[0, board, :] += spec[0, board, :] * wt * \
                                                       blank_nind[0,board,:]
                        weights[board] += wt*blank_nind[0,board,:]*numpy.isfinite(spec[0,board,:].data)
                    else:
                        logger.info("Rejecting Board %d nc %s due to sigma threshold" % (board, nc))
                else:
                    if not spec.mask[0,board,:].all() and wt > 1e-40:
                        avg_spectrum[0, board, :] += spec[0, board, :] * wt * \
                                                   blank_nind[0,board,:]
                        weights[board] += wt * blank_nind[0,board,:]*numpy.isfinite(spec[0,board,:].data)                   
                #spectrum[0, board, :] += nc.hdu.spectrum[:, board, :].mean(axis=0)
        for board in range(frequencies.shape[0]):
            self.spectrum[0, board, :] =avg_spectrum[0,board,:]/weights[board]
            self.spectrum.mask[0, board,:] = weights[board] == 0
        self.baseline()
        #blankarr = blankind.mean(axis=0).astype(numpy.bool)
        #blankarr.shape = (1, blankarr.shape[0], blankarr.shape[1])
        #self.spectrum[blankarr] = self.blank
        #spectrum /= float(len(nclist)+1)
        
        # for board in range(frequencies.shape[0]):
        #     weights = 0.0
        #     for nc in nclist:
        #         hdu = self._get_hdu(nc)
        #         spec = hdu._get_spectrum()
        #         wt = hdu._calc_weight(0, board, weight=weight)
        #         if weight == 'sigma' and threshold_sigma is not None:
        #             print nc, board, wt, threshold_sigma
        #             if 1/numpy.sqrt(wt) < threshold_sigma:
        #                 self.spectrum[0, board, :] += spec[0, board, :] * wt
        #                 weights += wt
        #             else:
        #                 logger.info("Rejecting Board %d nc %s due to sigma threshold" % (board, nc))
        #         else:
        #             self.spectrum[0, board, :] += spec[0, board, :] * wt
        #             weights += wt                    
        #         #spectrum[0, board, :] += nc.hdu.spectrum[:, board, :].mean(axis=0)
        #     self.spectrum[0, board, :] /= weights
        # #spectrum /= float(len(nclist)+1)

    # def fit_double_gaussian(self, board=5, pinit=None,
    #                         maxiter=1000, quiet=False):
    #     if self.header.ObsPgm != 'Map':
    #         raise LMTRedshiftError('fit_double_gaussian', "can only do this for Map data")
    #     if self.header.AccumMode != 'AccumCompressed':
    #         raise LMTRedshiftError('fit_double_gaussian', "can only do this for Compressed Continuum Data")
    #     if self.header.SwitchSettings != 'Acf':
    #         raise LMTRedshiftError('fit_double_gaussian', "can only do this for CC data with Switch Settings : ACF")
    #     ind = numpy.where(self.data.BufPos == 0)
    #     xi = numpy.degrees(self.data.XPos[ind])*3600
    #     yi = numpy.degrees(self.data.YPos[ind])*3600
    #     if int(self.header.ChassisNumber) in (1, 2):
    #         z = -self.data.AccAverage[ind[0], board]
    #     else:
    #         z = self.data.AccAverage[ind[0], board]
    #     zi = z - z[:20].mean()  # subtract out average
    #     if pinit is None:
    #         # let's make some guesses
    #         indmax = numpy.where(zi == zi.max())
    #         indmin = numpy.where(zi == zi.min())
    #         pinit = [zi.max(), xi[indmax], yi[indmax], 20., 20.0,
    #                  zi.min(), xi[indmin], yi[indmin], 20.0, 20.0]
    #     m = two_gaussian_fit(xi, yi, zi, pinit=pinit,
    #                          maxiter=maxiter, quiet=quiet)
    #     if m.status > 0 and m.status != 4:
    #         logger.info("Looks like minimization may have converged. Reporting details below")
    #         return m
        
    def blank_frequencies(self, freqrange):
        """
        Given a dictionary (with board numbers as keys) of frequency tuples
        this function sets those frequency channels to the blank value
        For eg. set freqrange to:
        freqrange = {3: [(95.6, 97.3),],
                     4: [(104.3, 104.5), (109.1, 110)]
                     }
        will blank channels between 95.6 to 97.3 GHz in board 3 and
        between 104.3 to 104.5 GHz and 109.1 to 110 GHz in board 4
        """
        if not isinstance(freqrange, dict):
            raise LMTRedshiftError("RedshiftScan.blank_frequencies", "freqrane should be a dictionary with board numbers as keys and a list (or tuples) of frequency tuples as values")
        for board, freqlist in freqrange.items():
            if board not in range(6):
                raise LMTRedshiftError("RedshiftScan.blank_frequencies", "Not a valid board number: %d" % board)
            for i, (fmin, fmax) in enumerate(freqlist):
                fmin, fmax = numpy.sort((fmin, fmax))
                if i == 0:
                    ind = numpy.logical_and(self.frequencies[board] >= fmin,
                                            self.frequencies[board] <= fmax)
                else:
                    ind = numpy.logical_or(ind,
                                           self.numpy.logical_and(self.frequencies[board] >= fmin,
                                                                  self.frequencies[board] <= fmax)
                                           )
                logger.info("Blanking board %d freq range: %s GHz to %s GHz" % (board, fmin, fmax))
            self.spectrum[:, board, ind] = self.blank

        
    def area(self, fmin, fmax):
        """
        Computes the area under the curve for a given window.  The
        units of fmin and fmax are in GHz.

        A new area object and area_sigma object is filled into the hdu
        object. Area is in K.km/s and uncertainty in same units.
        The area object and area_sigma object are also returned by this call
        """
        C = 2.99792458e+05
        if not hasattr(self, 'spectrum'):
            raise LMTRedshiftError('area', "No spectrum available. First run process_scan")
        if not hasattr(self, 'sigma'):
            raise LMTRedshiftError('area', "No sigma available. Make sure to remove baseline before computing area")
        if hasattr(self, 'compspectrum'):
            #we'll use this instead
            comp = True
            freq = self.compfreq
        else:
            comp = False
            freq = self.frequencies
        fmin, fmax = numpy.sort([fmin, fmax])
        # boards = []
        # for board in self.hdu.header.BoardNumber:
        #     if fmin > self.hdu.frequencies[board].min() and fmin <  self.hdu.frequencies[board].max():
        #         boards.append(board)
        #     elif 

        if comp:
            find = numpy.logical_and(freq>=fmin, freq<=fmax)
            #if window is None:
            bandid = None
            for i, f in enumerate(self.frequencies):
                ind = numpy.logical_and(f>=fmin, f<=fmax)
                if ind.any():
                    bandid = i
        else:
            bandid = None
            for i, f in enumerate(freq):
                find = numpy.logical_and(f>=fmin, f<=fmax)
                if find.any():
                    bandid = i
                    break
        if comp:
            dv = C*0.03125/freq[find]
            tdv = self.compspectrum[:, find].mean(axis=0) * dv
            integ_inten = tdv.sum()
            dv_square = dv**2
            integ_inten_rms = self.sigma[:, bandid].mean()*numpy.sqrt(dv_square.sum())
        else:
            dv = C*0.03125/freq[bandid, find]
            tdv = self.spectrum[:, bandid, find].mean(axis=0) * dv
            integ_inten = tdv.sum()
            dv_square = dv**2
            integ_inten_rms = self.sigma[:, bandid].mean()*numpy.sqrt(dv_square.sum())
        logger.info("Integ Intensity:  %f +/- %f K.km/s" % (integ_inten, integ_inten_rms))
        return integ_inten, integ_inten_rms

            
