import numpy
from dreampy3.lmtnetcdf import LMThdu
from dreampy3.ifproc.utils import LMTIFProcError
from dreampy3.logging import logger
from scipy.interpolate.interpolate import interp1d
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
#from scipy.interpolate import griddata as sci_griddata
import scipy.interpolate
logger.name = __name__

def process_signal_2(bb_level, chop, window=3, harm=2):
    npts = len(chop)
    nchannels = numpy.shape(bb_level)[1]

    ww = 2*window+1

    # create array of second harmonic values
    tcos = numpy.cos(harm*chop/8000*2.*numpy.pi)
    tsin = numpy.sin(harm*chop/8000*2.*numpy.pi)
    
    cresult = numpy.zeros((npts,nchannels))
    sresult = numpy.zeros((npts,nchannels))
    aresult = numpy.zeros((npts,nchannels))
    presult = numpy.zeros((npts,nchannels))
                    
    for i in range(nchannels):
        c_level = bb_level[:,i]*tcos
        s_level = bb_level[:,i]*tsin

        for j in range(window, npts-window):
            cresult[j,i] = numpy.sum(c_level[j-window:j+window+1])
            sresult[j,i] = numpy.sum(s_level[j-window:j+window+1])

        cresult[:window,i] = cresult[window,i]*numpy.ones(window)
        cresult[npts-window:,i] = cresult[npts-window-1,i]*numpy.ones(window)
        sresult[:window,i] = sresult[window,i]*numpy.ones(window)
        sresult[npts-window:,i] = sresult[npts-window-1,i]*numpy.ones(window)
        
        aresult[:,i] = 2/ww*numpy.sqrt(cresult[:,i]**2+sresult[:,i]**2)
        presult[:,i] = numpy.arctan2(sresult[:,i],cresult[:,i])
                                

    return(aresult,presult)

def process_chopped_encoder(chop, thresholds=[[15,45,181],[110,135]]):
    # create array of indices for main and ref based on chop array
    ang = (chop/8000*360)%180

    midx = numpy.where(numpy.logical_or(numpy.logical_and(ang > thresholds[0][0], ang < thresholds[0][1]),numpy.logical_and(ang>thresholds[0][2],ang<=180)))[0]
    ridx = numpy.where(numpy.logical_and(ang > thresholds[1][0], ang < thresholds[1][1]))[0]
    return midx, ridx

def process_chopped_signal(bb_level, chop, window=6,
                           thresholds=[[15,45,181],[110,135]]):
    '''
    gated chopper signal processor
    inputs:
         bb_level is npts by nchannels 2D array with baseband if data samples
         chop is chopper wheel position 0 to 8000 corresponds to 0 to 360 degrees
         window defines smoothing window of 2*window + 1 points.  The total smoothing
          window must span at least one chop cycle
         thresholds are positions for including data points in the main and ref
           thresholds[0] elements give main limits in degrees from 0 to 180
                         data are included if between thresholds[0][0] and thresholds[0][1]
                         OR if greater than thresholds[0][2]
           thresholds[1] elements give reference limits
                         data are included if between thresholds[1][0] and thresholds[1][1]
    output:
         result is a 2D array with npts samples to match input arrays and nchannels.
    '''

    # look at the shape of the arrays to determine if super sampled and reshape
    s1 = numpy.shape(bb_level)
    s2 = numpy.shape(chop)
    if len(s1) == 3 and len(s2) == 2:
        if s1[0] != s2[0] or s1[1] != s2[1]:
            return None
        bb_level = bb_level.reshape(s1[0]*s1[1], s1[2])
        chop = chop.reshape(s2[0]*s2[1])
        super_sample = s1[1]
        window = window*super_sample
    else:
        super_sample = 1


    npts = len(chop)
    nchannels = numpy.shape(bb_level)[-1]

    # define the smoothing window
    ww = 2*window+1

    # find indices where encoder value are within a range
    midx, ridx = process_chopped_encoder(chop, thresholds=thresholds)

    msig = numpy.zeros(npts)
    msig[midx] = 1
    rsig = numpy.zeros(npts)
    rsig[ridx] = 1

    result = numpy.zeros((npts,nchannels))

    for i in range(nchannels):
        channel_level = bb_level[:,i] # gets rid of "masked array

        # create a rolling sum of the main points
        msum = numpy.cumsum(numpy.insert(msig*channel_level,0,0))
        mrollsum = msum[ww:]-msum[:-ww]

        # to do this accurately we also need a rolling sum for normalization
        mnorm = numpy.cumsum(numpy.insert(msig,0,0))
        mrollnorm = mnorm[ww:]-mnorm[:-ww]

        # same procedure for reference points
        rsum = numpy.cumsum(numpy.insert(rsig*channel_level,0,0))
        rrollsum = rsum[ww:]-rsum[:-ww]

        # same normalization procedure for reference points
        rnorm = numpy.cumsum(numpy.insert(rsig,0,0))
        rrollnorm = rnorm[ww:]-rnorm[:-ww]

        # now compute difference between main and ref for all points 
        result[window:npts-window,i] = mrollsum/mrollnorm - rrollsum/rrollnorm
        result[:window,i] = result[window,i]*numpy.ones(window)
        result[npts-window:,i] = result[npts-window-1,i]*numpy.ones(window)

    # average the arrays back down if super sampled
    if super_sample > 1:
        result = numpy.mean(result.reshape(-1, super_sample, nchannels), axis=1)
        msig = numpy.mean(msig.reshape(-1, super_sample), axis=1)
        rsig = numpy.mean(rsig.reshape(-1, super_sample), axis=1)

    return(result)

# def process_chopped_signal(bb_level, chop, window=3, threshold=0.5,
#                            shift=0):
#     npts = len(chop)
#     nchannels = numpy.shape(bb_level)[1]
    
#     # define the smoothing window
#     ww = 2*window + 1

#     # create array of indices for main, ref, blank based on chop array
#     print('chop shift',  shift)
#     switch = numpy.cos(2.*(chop-shift)/8000*2.*numpy.pi)

#     # find indices where cos of encoder value exceeds a threshold
#     midx = numpy.where(switch > threshold)[0]
#     ridx = numpy.where(switch < -threshold)[0]

#     # create arrays for main and ref's identified by indices
#     msig = numpy.zeros(npts)
#     rsig = numpy.zeros(npts)

#     # load arrays with 1 where indices for main and ref identified
#     msig[midx] = 1.
#     rsig[ridx] = 1.

#     result = numpy.zeros((npts, nchannels))
#     main = numpy.zeros((npts, nchannels))
#     ref = numpy.zeros((npts, nchannels))
    
#     for i in range(nchannels):
#         channel_level = bb_level[:,i] # gets rid of "masked array

#         # create a rolling sum of the main points
#         msum = numpy.cumsum(numpy.insert(msig*channel_level,0,0))
#         mrollsum = msum[ww:]-msum[:-ww]

#         # to do this accurately we also need a rolling sum for normalization
#         mnorm = numpy.cumsum(numpy.insert(msig,0,0))
#         mrollnorm = mnorm[ww:]-mnorm[:-ww]
#         if not 0 in mrollnorm:
#             mrollsum /= mrollnorm

#         # same procedure for reference points
#         rsum = numpy.cumsum(numpy.insert(rsig*channel_level,0,0))
#         rrollsum = rsum[ww:]-rsum[:-ww]

#         # same normalization procedure for reference points
#         rnorm = numpy.cumsum(numpy.insert(rsig,0,0))
#         rrollnorm = rnorm[ww:]-rnorm[:-ww]
#         if not 0 in rrollnorm:
#             rrollsum /= rrollnorm

#         # now compute difference between main and ref for all points 
#         result[window:npts-window,i] = mrollsum - rrollsum
#         result[:window,i] = result[window,i]*numpy.ones(window)
#         result[npts-window:,i] = result[npts-window-1,i]*numpy.ones(window)

#         main[window:npts-window,i] = mrollsum
#         main[:window,i] = main[window,i]*numpy.ones(window)
#         main[npts-window:,i] = main[npts-window-1,i]*numpy.ones(window)

#         ref[window:npts-window,i] = rrollsum
#         ref[:window,i] = ref[window,i]*numpy.ones(window)
#         ref[npts-window:,i] = ref[npts-window-1,i]*numpy.ones(window)        

        
#     return(result, main, ref)


class IFProcScan(LMThdu):
    def __init__(self, data=None, header=None,
                 filename=None, groupname=None):
        LMThdu.__init__(self, data=data, header=header,
                        filename=filename)
        self.groupname = groupname
        self.FS = 100.0
        self.nyquist = self.FS/2
        
    def process_cross_scan(self,
                           stepsize=1.0, # in arcseconds
                           apply_Tsys=False,
                           Tsys=[],
                           baseline=True,
                          ):
        """
        Cross scan
        If apply_Tsys is True,
        pass a list or dictionary for Tsys with
        Tsys values.
        If baseline is True, removes the off-source value 
        from the overall measurements
        """
        if self.header.get('Dcs.ObsPgm') != 'CrossScan':
            logger.error("Not a cross-scan datafile")
            return
        self.numpixels = self.data.BasebandLevel.shape[1]
        logger.info("Processing IFProc Continuum CrossScan Observation with ObsNum %d" % int(self.header.get('Dcs.ObsNum')))
        ind = numpy.where(self.data.BufPos == 0)
        self.off_source = {}
        for pixel in range(self.numpixels):
            self.off_source[pixel] = numpy.histogram(self.data.BasebandLevel[ind, pixel].flatten())[1][:4].mean()
        if self.header.get('CrossScan.MapCoord') == 'Az':
            xpos = numpy.degrees(self.data.TelAzMap[ind])*3600.
        else:
            xpos = numpy.degrees(self.data.TelElMap[ind])*3600.
        mapsize = 2* int(((xpos[-1] - xpos[0])*0.95)/2.) # slightly smaller
        numsteps = int(mapsize/stepsize + 1)
        self.grid = numpy.linspace(-mapsize/2, mapsize/2, numsteps)
        self.continuum = numpy.zeros((self.grid.size, self.numpixels),
                                     dtype=self.data.BasebandLevel.dtype)
        cont = {}
        if baseline:
            for pixel in range(self.numpixels):
                cont[pixel] = (self.data.BasebandLevel[ind, pixel] - self.off_source[pixel]).flatten()
        else:
            for pixel in range(self.numpixels):
                cont[pixel] = self.data.BasebandLevel[ind, pixel].flatten()
        for pixel in range(self.numpixels):
            f = interp1d(xpos, cont[pixel])
            self.continuum[:, pixel] = f(self.grid)
        return self.continuum

    def process_cal_scan(self, msip1mm=True):
        """
        Process cal scan
        """
        idxh = numpy.where(self.data.BufPos == 3)
        idxc = numpy.where(self.data.BufPos == 2)
        self.numpixels = self.data.BasebandLevel.shape[1]
        self.Tsys = numpy.zeros(self.numpixels)
        self.Ph = numpy.zeros(self.numpixels)
        self.Pc = numpy.zeros(self.numpixels)
        for pixel in range(self.numpixels):
            self.Ph[pixel] = self.data.BasebandLevel[idxh, pixel].mean()
            self.Pc[pixel] = self.data.BasebandLevel[idxc, pixel].mean()
            self.Tsys[pixel] = 280 * self.Pc[pixel]/(self.Ph[pixel] - self.Pc[pixel])
        if msip1mm:
            print("Tsys: P0 USB: %.3f K" % self.Tsys[0])
            print("Tsys: P0 LSB: %.3f K" % self.Tsys[1])
            print("Tsys: P1 USB: %.3f K" % self.Tsys[3])
            print("Tsys: P1 LSB: %.3f K" % self.Tsys[2])
        return self.Tsys

    def process_chopped_cal_scan(self, msip1mm=True):
        """
        Process cal scan
        """
        if self.header.get('Msip1mm.BeamChopperActState', 0) != 3:
            raise LMTIFProcError('IFProc Error', 'Scan is not a chopped scan')
        bb_level = self.data.BasebandLevel
        chop = self.data.BeamChopperActPos
        shift = 1000
        window = 12
        threshold = 0.5
        signal, main, ref = process_chopped_signal(bb_level, chop,
                                                   window=window, threshold=threshold,
                                                   shift=shift)
        idxh = numpy.where(self.data.BufPos == 3)
        idxc = numpy.where(self.data.BufPos == 2)
        self.numpixels = self.data.BasebandLevel.shape[1]
        self.Tsys = numpy.zeros(self.numpixels)
        self.Ph = numpy.zeros(self.numpixels)
        self.Pc = numpy.zeros(self.numpixels)
        for pixel in range(self.numpixels):
            self.Ph[pixel] = main[idxh, pixel].mean()
            self.Pc[pixel] = ref[idxc, pixel].mean()
            self.Tsys[pixel] = 280 * self.Pc[pixel]/(self.Ph[pixel] - self.Pc[pixel])
        if msip1mm:
            print("Tsys: P0 USB: %.3f K" % self.Tsys[0])
            print("Tsys: P0 LSB: %.3f K" % self.Tsys[1])
            print("Tsys: P1 USB: %.3f K" % self.Tsys[3])
            print("Tsys: P1 LSB: %.3f K" % self.Tsys[2])
        return self.Tsys

    def process_chopped_on_scan(self, msip1mm=True, window=12, 
                                harm=2):
        """
        Process ON scan
        """
        if self.header.get('Msip1mm.BeamChopperActState', 0) != 3:
            raise LMTIFProcError('IFProc Error', 'Scan is not a chopped scan')
        bb_level = self.data.BasebandLevel
        chop = self.data.BeamChopperActPos
        #shift = 780
        #window = 12
        #threshold = 0.5
        #signal, main, ref = process_chopped_signal(bb_level, chop,
        #                                           window=window, threshold=threshold,
        #                                          shift=shift)
        aresult, presult = process_signal_2(bb_level, chop, window=12, harm=2)
        #return (signal, main, ref)
        return (aresult, presult)
        
    def process_map(self, remove_offset=False, numoff=20,
                    scigrid=False,
                    **kwargs):
        """
        processes Map Obspgm that has been made using compressed
        continuum mode. Uses a regridding algorithm
        and uses some kwargs arguments to derive output
        grid size and sampling
        """
        if self.header.ObsPgm not in ('Map', 'Lissajous'):
            logger.error("Not a Map datafile")
            return
        else:
            maptype = self.header.ObsPgm
        logger.info("Processing MSIP 1mm Continuum Map data and regridding for Observation with ObsNum %d" % int(self.header.ObsNum))
        self.numpixels = self.data.BasebandLevel.shape[1]
        xlength = numpy.degrees(self.header.get('%s.XLength' % maptype))*3600.0
        ylength = numpy.degrees(self.header.get('%s.YLength' % maptype))*3600.0
        if maptype == 'Lissajous':
            xlength = xlength/numpy.cos(numpy.radians(45))
            xlength = ylength/numpy.cos(numpy.radians(45))
        ind = numpy.where(self.data.BufPos == 0)
        xpos = numpy.degrees(self.data.TelAzMap[ind])*3600.
        ypos = numpy.degrees(self.data.TelElMap[ind])*3600.
        rows = self.header.get('Map.RowsPerScan')
        z = {}
        self.off_source = {}
        for chan in range(self.numpixels):
            z[chan] = self.data.BasebandLevel[ind, chan].flatten()
            self.off_source[chan] = numpy.histogram(self.data.BasebandLevel[ind, chan].flatten())[1][:4].mean()
            if remove_offset:
                z[chan] = z[chan] - self.off_source[chan]
            print(z[chan].shape)
        ramp = kwargs.get('ramp', 5.)
        numpoints = kwargs.get('numpoints', 100)
        numypoints = kwargs.get('numypoints', 100)
        xlength = xlength * (1.-ramp/100.)
        ylength = ylength * (1.-ramp/100.)
        ind = numpy.logical_and(xpos > -xlength/2., xpos < xlength/2.)
        xpos, ypos = xpos[ind], ypos[ind]
        # add a tiny random number to stop griddata from crashing when two pixels are same
        xpos = xpos + numpy.random.random(xpos.size)*1e-6
        ypos = ypos + numpy.random.random(ypos.size)*1e-6
        for chan in range(self.numpixels):
            z[chan] = z[chan][ind]
        ind = numpy.logical_and(ypos > -ylength/2., ypos < ylength/2.)
        xpos, ypos = xpos[ind], ypos[ind]
        for chan in range(self.numpixels):
            z[chan] = z[chan][ind]
        self.xi = numpy.linspace(-xlength/2, xlength/2, numpoints)
        self.yi = numpy.linspace(-ylength/2, ylength/2, numypoints)
        print("Making %d x %d map" % (numpoints, numypoints))
        #self.z = z
        #self.xpos = xpos
        #self.ypos = ypos
        self.BeamMap = numpy.zeros((self.yi.size, self.xi.size, self.numpixels))
        for chan in range(self.numpixels):
            #self.BeamMap[chan] = numpy.zeros((self.yi.size, self.xi.size),
            #                                 dtype='float')
            if scigrid:
                self.BeamMap[:, :, chan] = scipy.interpolate.griddata(xpos, ypos, z[chan], (self.xi, self.yi),
                                                                      method='cubic')
            else:
                self.BeamMap[:, :, chan] = griddata(xpos, ypos, z[chan], self.xi, self.yi)


    def interpolate_scan(self, **kwargs):
        """
        interpolates scan and writes out new 
        Baseband Level data for interpolated data
        """
        self.numpixels = self.data.BasebandLevel.shape[1]
        time = self.data.BasebandTime
        t0 = time[0]
        time = time - t0
        tnew = numpy.arange(time[0], time[-1]+1e-6, 1/self.FS)
        self.IFProcTime = tnew
        self.BasebandLevel = numpy.zeros((tnew.size, self.numpixels))
        for chan in range(self.numpixels):
            self.BasebandLevel[:, chan] = interp1d(time, self.data.BasebandLevel[:, chan])(tnew)
        self.XPos = interp1d(time, self.data.TelAzMap)(tnew)
        self.YPos = interp1d(time, self.data.TelElMap)(tnew)
        self.BufPos = interp1d(time, self.data.BufPos.astype(numpy.bool))(tnew).astype(numpy.bool)


    # def process_detrended_filtered_map(self, nfwhm=2, lpf_freq=5.0,
    #                                    **kwargs):
    #     """
    #     Interpolates scan, filters based on known frequencies,
    #     detrends, cleans map, and then puts it into regular 
    #     grid
    #     Based on Lindy Blackburn's implementation
    #     """
    #     self.interpolate_scan()
    #     maptype = self.header.ObsPgm
    #     t = self.IFProcTime
    #     power = numpy.zeros(self.BasebandLevel.shape)
    #     for chan in range(self.numchannels):
    #         power[:, chan] = self.BasebandLevel[:, chan] - self.BasebandLevel[:, chan].mean()
    #     # if self.header.utdate().year < 2016:
    #     #     f1 = 1.25
    #     #     f2 = 1.55
    #     #     f3 = 11.47
    #     #     cfreqs = numpy.array([f1, f2, 2*f1, f3, f3-f2, f3-2*f2, f3+f2])
    #     # elif self.header.utdate().year == 2016 and self.header.utdate().month <= 3:
    #     #     f1 = 1.17
    #     #     #cfreqs = numpy.array([f1, f2, 2*f1])
    #     #     cfreqs = numpy.array([f1, 2*f1])
    #     # else:
    #     #     f1 = 1.20
    #     #     f3 = 10.00
    #     #     cfreqs = numpy.array([f1, 2*f1, f3, f3-f1])
    #     # power_p = {}
    #     # for key in ('APower', 'BPower'):
    #     #     power_p[key] = power[key].copy()
    #     power_p = power.copy()
    #     (b, a) = butter(6, lpf_freq/self.nyquist, "lowpass")
    #     for key in ('APower', 'BPower'):
    #         power_p[key] = lfilter(b, a, lfilter(b, a, power_p[key])[::-1])[::-1]
    #     # w = 0.3
    #     # for f in cfreqs:
    #     #     (b, a) = butter(3, [(f-w)/self.nyquist, (f+w)/self.nyquist], 'bandstop')
    #     #     for key in ('APower', 'BPower'):
    #     #         power_p[key] = lfilter(b, a, lfilter(b, a, power_p[key])[::-1])[::-1]
    #     # moving minimum filter
    #     tfwhm = nfwhm * 10. / numpy.degrees(self.header.get('%s.ScanRate' % maptype)*3600) # beam size over tracking time
    #     (b, a) = butter(3, [0.02/self.nyquist, 2./self.nyquist], 'bandpass')
    #     for chan in range(self.numchannels):
    #         y2 = lfilter(b, a, lfilter(b, a, power_p[:, chan])[::-1])[::-1]
    #         bottom = minimum_filter1d(power_p[:, chan], int(tfwhm * self.FS), mode='nearest')
    #         bottom2 = minimum_filter1d(y2, int(tfwhm * self.FS), mode='nearest')
    #         imin = numpy.nonzero((power_p[:, chan] == bottom) | (y2 == bottom2))[0]
    #         tmin = numpy.hstack((t[0]-.1, t[imin], t[-1]+.1))
    #         ymin = numpy.hstack((power_p[key][imin[0]], power_p[key][imin], power_p[key][imin[-1]]))
    #         ymints = interp1d(tmin, ymin)(t)
    #         (b, a) = butter(3, 0.2/self.nyquist, 'lowpass')
    #         yrefl = numpy.hstack(((2*ymints[0]-ymints[1:])[::-1], ymints, (2*ymints[-1]-ymints[:-1])[::-1]))
    #         ybaseline = lfilter(b, a, lfilter(b, a, yrefl[::-1])[::-1])[len(ymints):2*len(ymints)]
    #         ybaseline += numpy.median(power_p[key] - ybaseline)
    #         #setattr(self.data, key, power_p[key] - ybaseline)
    #     for val in ('XPos', 'YPos', 'Vlbi1mmTpmTime', 'BufPos'):
    #         if val == 'BufPos':
    #             self.data.BufPos = (~self.BufPos).astype(numpy.int)
    #         else:
    #             setattr(self.data, val, getattr(self, val))
        
    
    
    def process_scan(self, detrend=True, **kwargs):
        """wrapper for processing any onemm scan"""
        obspgm = self.header.get('Dcs.ObsPgm')
        if obspgm == 'Cal':
            return self.process_cal_scan(**kwargs)
        elif obspgm == 'CrossScan':
            return self.process_cross_scan(**kwargs)
        elif obspgm in ('Map', 'Lissajous'):
            if detrend:
                self.process_detrended_filtered_map(**kwargs)
            return self.process_map(**kwargs)
    
