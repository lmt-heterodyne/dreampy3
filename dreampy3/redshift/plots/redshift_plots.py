"""
Extension to DreamPlot
for redshift visualizations
"""
from dreampy3.plots import DreamPlot
from dreampy3.plots.dream_plots import DreamPlotBase, DreamPlotChart
from dreampy3.redshift.netcdf import RedshiftNetCDFFile, RedshiftScan
from dreampy3.redshift.utils import LMTRedshiftError
from dreampy3.redshift.utils.linetable import freqdict as linefreqdict
from dreampy3.redshift.utils.spectrum_utils import apply_bad_lags, nanmean, blank_to_nan, blank
from dreampy3.utils.smoothing import rebin
from dreampy3.utils.curve_fitting import Gauss1DFit, Gauss2DFit, Gauss1DFitBase
import numpy
from matplotlib import cm
from matplotlib import patches
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
from matplotlib.font_manager import FontProperties
import math
import os
import types
from dreampy3.logging import logger
logger.name = __name__

class RedshiftPlotBase(DreamPlotBase):
    """
    Main class to visualize Redshift data and scans
    Do not use directly. Instead use extended class RedshiftPlot
    or RedshiftPlotChart
    """
    def __init__(self, **kwargs):
        #super(RedshiftPlotBase, self).__init__(**kwargs)
        #self.set_title('Redshift Receiver Plot Window')
        self.image = None
        self.panel_plots = []
        self.subpanel_plots = []
        self.basewinlines = []
                
    def _get_hdu(self, nc):
        if isinstance(nc, RedshiftNetCDFFile):
            hdr = nc.hdu.header
            data = nc.hdu.data
            fname = nc.filename
        elif isinstance(nc, RedshiftScan):
            hdr = nc.header
            data = nc.data
            fname = nc.filename
        else:
            raise LMTRedshiftError("Argument Error", "Argument nc should be of type RedshiftNetCDFFile or RedshiftScan")
        return hdr, data, fname

    def _get_repeats(self, nc):
        hdr, data, fname = self._get_hdu(nc)
        repeats = data.AccData.shape[0]
        if hdr.ObsPgm == 'Bs':
            return repeats/2
        else:
            return repeats

    def plot_data(self, nc, board=None, datatype='AccData',
                  rpt=None, xtype='Chan', norm=True, hold=False,
                  dosubplots=True, smooth=1,
                  *args, **kwargs):
        """
        Plots one of datatype = 'AccData', or 'RefData'
        against xtype.
        xtype can be one of 'Chan' or 'Freq'
        None (just uses index numbers).
        The smooth parameter allows a simple boxcar average of input
        data. If smooth=1, raw data is plotted, smooth=10 does a boxcar
        average of 10 adjacent samples.
        """
        self.check_alive()        
        if not hold:
            self.clear()
        if datatype not in ('AccData', 'RefData'):
            raise LMTRedshiftError("Argument Error", """Datatype should be one of 'AccData', 'RefData'""")
        if xtype not in ('chan', 'Chan', 'c', 'Freq', 'freq', 'f', None):
            raise LMTRedshiftError("Argument Error", """xtype should be one of 'Chan', 'Freq', None""")
        hdr, data, fname = self._get_hdu(nc)
        if xtype is None:
            xlabel = 'Index'
            x = numpy.arange(data.AccData.shape[2])
        elif xtype.lower() in ('chan', 'c'):
            xlabel = 'Chan'
            x = numpy.arange(data.AccData.shape[2])
        xlabel = xtype
        if board is None:
            #do all boards
            boards = hdr.BoardNumber
        else:
            if board < 0 or board > 5:
                raise LMTRedshiftError("Argument Error",
                                       """board should be between 0 and 5""")
            boards = [board]
        if rpt is None:
            #do all repeats
            #repeats = xrange(data.AccData.shape[0])
            repeats = self._get_repeats(nc)
            rpt = 1
        else:
            repeats = [rpt-1]
        for repeat in xrange(repeats):
            for board in boards:
                if dosubplots:
                    self.add_subplot(2, 3, board+1)
                dt = getattr(data, datatype).copy()[repeat, board, :]
                dt = apply_bad_lags(dt, int(hdr.ChassisNumber), board)
                if norm:
                    dt = dt/float(data.AccSamples[repeat, board])
                self.plot(rebin(x, smooth=smooth), rebin(dt, smooth=smooth),
                          drawstyle='steps',
                          label="%s C%dB%d r%d" % (datatype, int(hdr.ChassisNumber), board+1, repeat),
                          *args, **kwargs)
                self.set_legend(loc='best', prop=FontProperties(size="xx-small"))
                self.set_xlabel(xlabel)
                self.set_ylabel(datatype)
        self.set_figtitle('%s: Source %s' % (hdr.obs_string(), hdr.SourceName))

    def plot_accdata(self, nc, board=None, rpt=None,
                     xtype='Chan', norm=True, hold=False,
                     dosubplots=True, smooth=1, *args, **kwargs):
        """
        Plots redshift AccData in a given scan (pass the nc
        instance) against xtype.
        xtype can be one of
        """
        return self.plot_data(nc, board=board, rpt=rpt, xtype=xtype,
                              norm=norm, hold=hold, dosubplots=dosubplots,
                              smooth=smooth, datatype='AccData', 
                              *args, **kwargs)

    def plot_refdata(self, nc, board=None, rpt=None,
                     xtype='Chan', norm=True, hold=False,
                     dosubplots=True, smooth=1, *args, **kwargs):
        """
        Plots redshift RefData in a given scan (pass the nc
        instance) against xtype.
        xtype can be one of
        """
        return self.plot_data(nc, board=board, rpt=rpt, xtype=xtype,
                              norm=norm, hold=hold, dosubplots=dosubplots,
                              smooth=smooth, datatype='RefData', 
                              *args, **kwargs)    

    def plot_spectra(self, nc, board=None, avgrepeats=True, rpt=None,
                     hold=False, *args, **kwargs):
        """Plots the redshift spectral data. This is to be obtained
        from the ACF in a prescribed way. If the spectral data is not
        available in the scan, it will produce an error"""
        if isinstance(nc, RedshiftNetCDFFile):
            hdu = nc.hdu
        elif isinstance(nc, RedshiftScan):
            hdu = nc
        else:
            raise LMTRedshiftError("plot spectra", "Argument nc should be of type RedshiftNetCDFFile or RedshiftScan")
        if hasattr(hdu, 'spectrum'):
            spectrum = blank_to_nan(hdu.spectrum)
            frequencies = hdu.frequencies
        elif hasattr(hdu.data, 'spectrum'):
            spectrum = blank_to_nan(hdu.data.spectrum)
            frequencies = hdu.data.frequencies
        else:
            raise LMTRedshiftError("plot spectra", """There are no spectra in this hdu. Make sure to run process_scan on the HDU first""")
        if hasattr(hdu, 'ylabel'):
            ylabel = hdu.ylabel
        elif hasattr(hdu.header, 'ylabel'):
            ylabel = hdu.header.ylabel
        else:
            ylabel = ''
        self.check_alive()        
        if not hold:
            self.clear()
        if board is None:
            boards = hdu.header.BoardNumber
        else:
            boards = [board]
        if rpt is None:
            #do all repeats
            #repeats = xrange(hdu.data.AccData.shape[0])
            repeats = self._get_repeats(nc)
            rpt = 1
        else:
            repeats = [rpt-1]
        if avgrepeats:
            for board in boards:
                self.plot(frequencies[board, :], spectrum[:, board, :].mean(axis=0),
                          drawstyle='steps-mid', label='%d.%d' % (int(hdu.header.ChassisNumber), board))
        else:
            for repeat in xrange(repeats):
                for board in boards:
                    self.plot(frequencies[board, :], spectrum[repeat, board, :],
                              drawstyle='steps-mid', label='%d.%d r%d' % (int(hdu.header.ChassisNumber), board, repeat))
        self.set_xlabel("Frequency (GHz)")
        self.set_legend(loc='best', prop=FontProperties(size="xx-small"))
        self.set_ylabel(ylabel)
        self.set_subplot_title('%s: Source %s' % (hdu.header.obs_string(), hdu.header.SourceName))        

    def reset_basewinlines(self):
        for line in self.basewinlines:
            line.set_visible(False)
        self.redraw_plot()
        self.basewinlines = []
        self.basewin = []
        self.winpair = []
        print("Press 'n' to set window limits and 'e' to exit, and 'r' to reset")
        
    def keypress_event(self, event):
        if not event.inaxes: return
        if event.key.lower() == 'n':
            if len(self.winpair) % 2 == 0:
                #start new window pair
                self.winpair = []
                self.winpair.append(event.xdata)
                print("xwin1 = %.3f. Press 'n' to mark other end of window. 'e' to exit; 'r' to reset" % event.xdata)
                ymin, ymax = self.get_ylims()
                self.basewinlines.append(self.vline(event.xdata, ymin=ymin,
                                                    ymax=ymax, linestyle='--',
                                                    color='r', linewidth=2,
                                                    label='_nolegend_'))
            else:
                #finish existing pair
                self.winpair.append(event.xdata)
                self.basewin.append(self.winpair)
                print("xwin1 = %.3f, xlim2 = %.3f. Press 'n' to start marking next window., Press 'e' to  exit, 'r' to reset" % (self.winpair[0], self.winpair[1]))
                ymin, ymax = self.get_ylims()
                self.basewinlines.append(self.vline(event.xdata, ymin=ymin,
                                                    ymax=ymax, linestyle='--',
                                                    color='r', linewidth=2,
                                                    label='_nolegend_'))
        elif event.key.lower() == 'e':
            if len(self.winpair) != 2:
                print("Last window was not completed. It is being dropped. Full baseline window: %s" % self.basewin)
                return
            if hasattr(self.nc.hdu, 'windows') and type(self.nc.hdu.windows) == types.DictType:
                self.nc.hdu.windows[self.board] = self.basewin
            elif hasattr(self.nc.hdu.data, 'windows') and type(self.nc.hdu.data.windows) == types.DictType:
                self.nc.hdu.data.windows[self.board] = self.basewin
            else:
                raise LMTRedshiftError("make_windows", "netcdf instance has no windows attribute")
            self.disconnect_event(self.basewinsig)
            print("Exiting make_windows. Final baseline windows: %s" % self.basewin)
        elif event.key.lower() == 'r':
            self.reset_basewinlines()


    def make_windows(self, nc, board=0):
        self.nc = nc
        self.board = board
        self.reset_basewinlines()
        self.basewinsig = self.connect_event('key_press_event',
                                             self.keypress_event)

    def plot_tsys(self, nc, board=None, hold=False,
                  set_subplot_title=True,
                  *args, **kwargs):
        """Plots the redshift Tsys data. This is to be obtained
        from the CalObsNum scan in a prescribed way. If the Cal data is not
        available in the scan, it will produce an error"""
        self.check_alive()        
        if not hold:
            self.clear()
        if isinstance(nc, RedshiftNetCDFFile):
            hdu = nc.hdu
        elif isinstance(nc, RedshiftScan):
            hdu = nc
        else:
            raise LMTRedshiftError("plot tsys", "Argument nc should be of type RedshiftNetCDFFile or RedshiftScan")
        if not hasattr(hdu, 'cal'):
            raise LMTRedshiftError("plot tsys", "hdu does not have cal object. Process first with get_cal()")
        if board is None:
            boards = hdu.header.BoardNumber
        else:
            boards = [board]
        for board in boards:
            self.plot(hdu.frequencies[board, :], hdu.cal.Tsys[board, :],
                      drawstyle='steps-mid', label='%d.%d' % (int(hdu.header.ChassisNumber), board))
        self.set_xlabel("Frequency (GHz)")
        self.set_legend(loc='best', prop=FontProperties(size="xx-small"))
        self.set_ylabel("Tsys (K)")
        if set_subplot_title:
            self.set_subplot_title('%s: Source %s' % (hdu.header.obs_string(), hdu.header.SourceName))        

    def plot_crosscan(self, nc, board=None, hold=False,
                      fitgauss=False, 
                      *args, **kwargs):
        """Plots the redshift crossscan data. """
        self.check_alive()        
        if not hold:
            self.clear()
        if isinstance(nc, RedshiftNetCDFFile):
            hdu = nc.hdu
        elif isinstance(nc, RedshiftScan):
            hdu = nc
        else:
            raise LMTRedshiftError("plot tsys", "Argument nc should be of type RedshiftNetCDFFile or RedshiftScan")
        if not hasattr(hdu, 'continuum'):
            raise LMTRedshiftError("plot crosscan", "hdu does not have continuum object. Process first with process_scan()")
        if board is None:
            boards = hdu.header.BoardNumber
        else:
            boards = [board]
        if hdu.header.get('CrossScan.MapCoord') == 'Az':
            x = numpy.degrees(hdu.data.XPos)*3600
            xlabel = 'Az Offset (")'
        elif hdu.header.get('CrossScan.MapCoord') == 'El':
            x = numpy.degrees(hdu.data.YPos)*3600
            xlabel = 'El Offset (")'
        for board in boards:
            # if align_sign:
            #     if int(hdu.header.ChassisNumber) in (1, 2):
            #         sgn = -1.0
            #     else:
            #         sgn = 1.0
            # else:
            #     sgn = 1.0
            self.plot(x, hdu.continuum[:, board],
                      drawstyle='steps-mid', label='%d.%d' % (int(hdu.header.ChassisNumber), board))
            if fitgauss:
                fit = Gauss1DFit(x, hdu.continuum[:, board])
                out = fit._run()
                print(out.pprint())
                self.plot(x, out.y, 'r-')
                                 
        self.set_xlabel(xlabel)
        self.set_legend(loc='best', prop=FontProperties(size="xx-small"))
        self.set_ylabel("TA*(K)")
        self.set_subplot_title('%s: Source %s' % (hdu.header.obs_string(), hdu.header.SourceName))

    def plot_composite_spectrum(self, nc, hold=False,
                                *args, **kwargs):
        """Plots the redshift composite data. """
        self.check_alive()        
        if not hold:
            self.clear()
        if isinstance(nc, RedshiftNetCDFFile):
            hdu = nc.hdu
        elif isinstance(nc, RedshiftScan):
            hdu = nc
        else:
            raise LMTRedshiftError("plot composite spectrum", "Argument nc should be of type RedshiftNetCDFFile or RedshiftScan")
        if hasattr(hdu, 'compspectrum'):
            spectrum = blank_to_nan(hdu.compspectrum)
            frequencies = hdu.compfreq
        elif hasattr(hdu.data, 'compspectrum'):
            spectrum = blank_to_nan(hdu.data.compspectrum)
            frequencies = hdu.data.compfreq
        else:
            raise LMTRedshiftError("plot composite spectra", """There are no composite spectra in this hdu. Make sure to run make_composite_scan on the HDU first""")
        self.plot(frequencies, spectrum.mean(axis=0), drawstyle='steps-mid',
                  label="CompSpectrum")
        self.set_xlabel("Frequency (GHz)")
        self.set_legend(loc='best', prop=FontProperties(size="xx-small"))
        self.set_ylabel("TA*(K)")
        self.set_subplot_title('%s: Source %s' % (hdu.header.obs_string(), hdu.header.SourceName))

    def plot_line_frequencies(self, nc, yloc=None, elim=[],
                              minint=0.0, z=0.0, txtsize=6.0, hold=True,
                              *args, **kwargs):
        self.check_alive()        
        if not hold:
            self.clear()
        if isinstance(nc, RedshiftNetCDFFile):
            hdu = nc.hdu
        elif isinstance(nc, RedshiftScan):
            hdu = nc
        else:
            raise LMTRedshiftError("plot spectrum", "Argument nc should be of type RedshiftNetCDFFile or RedshiftScan")
        fmin, fmax = hdu.frequencies.min(), hdu.frequencies.max()
        freqs = numpy.array(linefreqdict.keys())
        redfreqs = freqs/(z+1.)
        ind = numpy.logical_and(redfreqs>=fmin, redfreqs<=fmax)
        if yloc is None:
            ymin, ymax = self.get_ylims()
            yloc = ymin + (ymax-ymin)*0.8
        for line in freqs[ind]:
            if linefreqdict[line][2] >= minint:
                if linefreqdict[line][0] not in elim:
                    freqtext = "--"+linefreqdict[line][0]
                    freq = line/(z+1)
                    self.set_text(freq, yloc, freqtext,
                                  rotation='vertical',
                                  horizontalalignment='center',
                                  verticalalignment='center',
                                  size=txtsize)

    def plot_xypos(self,  nc, bufpos=0, hold=True, *args, **kwargs):
        """
        Plots the xy positions in the nc file. If bufpos is
        set to any negative number all XPos and YPos are plotted
        """
        self.check_alive()        
        if not hold:
            self.clear()
        hdr, data, fname = self._get_hdu(nc)
        if bufpos < 0:
            xpos = data.XPos
            ypos = data.YPos
        else:
            ind = numpy.where(data.BufPos == bufpos)
            xpos = numpy.degrees(data.XPos[ind])*3600.
            ypos = numpy.degrees(data.YPos[ind])*3600.
        units = "(arcsec)"
        self.plot(xpos, ypos, *args, **kwargs)
        self.set_xlabel('XPos %s' % units)
        self.set_ylabel('YPos %s' % units)
        self.set_subplot_title('%s: Source %s' % (os.path.basename(fname), hdr.get('Source.SourceName', '')))

    def plot_map_contour(self, nc,
                         nlevels=10,
                         cmap=cm.get_cmap('Spectral'),
                         board=1,
                         hold=False):
        self.check_alive()
        if not hold:
            self.clear()
        if isinstance(nc, RedshiftNetCDFFile):
            hdu = nc.hdu
        elif isinstance(nc, RedshiftScan):
            hdu = nc
        else:
            raise LMTRedshiftError("plot map", "Argument nc should be of type RedshiftNetCDFFile or RedshiftScan")
        if hasattr(hdu, 'BeamMap'):       
            beammap = hdu.BeamMap
        else:
            raise LMTRedshiftError("plot map", "No BeamMap attribute found. Please run process_scan() on the hdu first")
        hdr, data, fname = self._get_hdu(nc)        
        xi = hdu.xi
        yi = hdu.yi
        stepsize = xi[1] - xi[0]
        zi = beammap[:, :, board]
        X, Y = numpy.meshgrid(xi,yi)
        
        cs = self.contour(xi, yi, zi, nlevels, linewidths=0.5, colors='k')
        cs = self.contourf(xi, yi, zi, nlevels, cmap=cmap)
        self.colorbar()

        def format_coord(x, y):
            ind = numpy.logical_and(numpy.abs(X-x)<stepsize,
                                    numpy.abs(Y-y)<stepsize)
            zzz = nanmean(zi[ind])
            return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x, y, zzz)

        self.set_xlim(xi[0], xi[-1])
        self.set_ylim(yi[0], yi[-1])
        ax, kw = self.plotobj._get_current_axes()
        ax.format_coord = format_coord
        self.redraw_plot()
        self.set_subplot_title('%s: Source %s; Board:%d' % (os.path.basename(fname), hdr.get('Source.SourceName', ''), board), size=10)
        
        
    def implot_map(self, nc, board=1,vmin=None,
                   vmax=None, extents=None,
                   colorbar=True, hold=False, 
                   **kwargs):
        self.check_alive()
        if not hold:
            self.clear()
        if isinstance(nc, RedshiftNetCDFFile):
            hdu = nc.hdu
        elif isinstance(nc, RedshiftScan):
            hdu = nc
        else:
            raise LMTRedshiftError("implot map", "Argument nc should be of type RedshiftNetCDFFile or RedshiftScan")
        if hasattr(hdu, 'BeamMap'):       
            beammap = hdu.BeamMap
        else:
            raise LMTRedshiftError("implot map", "No BeamMap attribute found. Please run process_scan() on the hdu first")
        hdr, data, fname = self._get_hdu(nc)
        xi = hdu.xi
        yi = hdu.yi
        stepsize = xi[1] - xi[0]
        zi = beammap[:, :, board]
        if extents is None:
            extents = [xi[0], xi[-1], yi[0], yi[-1]]
        
        print(extents)
        self.image = self.imshow(zi, cmap=cm.get_cmap('Spectral'),
                                 interpolation='bilinear', vmin=vmin, vmax=vmax,
                                 origin='lower', extent=extents,
                                 **kwargs)
        X, Y = numpy.meshgrid(xi,yi)
        def format_coord(x, y):
            ind = numpy.logical_and(numpy.abs(X-x)<stepsize,
                                    numpy.abs(Y-y)<stepsize)
            zzz = nanmean(zi[ind])
            return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x, y, zzz)


        ax, kw = self.plotobj._get_current_axes()
        ax.format_coord = format_coord
        self.set_subplot_title('%s: Source %s; Board:%d' % (os.path.basename(fname), hdr.get('Source.SourceName', ''), board), size=10)
        if colorbar:
            self.colorbar()
        return self.image


class RedshiftPlot(RedshiftPlotBase, DreamPlot):
    def __init__(self, **kwargs):
        #super(RedshiftPlot, self).__init__(**kwargs)
        DreamPlot.__init__(self, **kwargs)
        RedshiftPlotBase.__init__(self, **kwargs)
        #self.set_title('Redshift Receiver Plot Window')

class RedshiftPlotChart(RedshiftPlotBase, DreamPlotChart):
    def __init__(self, chart_prop=None, **kwargs):
        DreamPlotChart.__init__(self, chart_prop, **kwargs)
        RedshiftPlotBase.__init__(self, **kwargs)
