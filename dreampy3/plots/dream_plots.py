#import gtk
import os
import numpy
import datetime 
#from dreampy import dreampyParams
from dreampy3.utils import DreampyGeneralError
from dreampy3.utils.lmtcat_to_ephem import read_lmt_cat, source_uptimes, planetary_uptimes
from dreampy3.logging import logger
from dateutil import parser
from matplotlib.dates import  DateFormatter 
logger.name = __name__
from dreampy3.plots.plots import MainView, BasePlotView, ChartView, ChartProperties
from dreampy3.observing.locations import location, LMT, GBT, Amherst

class DreamPlotBase:
    """This is the base class for dreampy plotting
    Do not use this class directly. Instead use the extended
    class DreamPlot or DreamPlotChart"""
    def __init__(self, **kwargs):
        self.plotobj = BasePlotView(**kwargs)
        self.alive = True

    def check_alive(self):
        """A simple routine that checks if the plot window is still alive
        """
        if not self.alive:
            raise DreampyGeneralError("PlotWindow", "Plot window has been killed, restart plot window")
        
    def add_subplot(self, nrows, ncols, plotnum, name=None, 
                    **kwargs):
        """Adds a subplot where nrows is number of rows in grid,
        and ncols is number of columns in grid, and plotnum is the
        plotnumber that always starts at 1.
        name is an optional string that you want to associate with this
        particular subplot. If name is not given, it will default to
        plotnum"""
        self.check_alive()
        return self.plotobj.add_subplot(nrows, ncols, plotnum, name,
                                        **kwargs)

    def plot(self, *args, **kwargs):
        """Plots given parameters such as data to to the current axes
        unless specified by current_axes_name='axes_name'.  If no subplot has been created a default (1,1,1) axes will be created
        NOTE: ALL AXES FUNCTIONS have latter property!!! """
        self.check_alive()
        return self.plotobj.plot(*args, **kwargs)

    def plot_uptimes(self, catalog_fname=None, delimiter=',', day=None,
                     add_planets=True, sources=[],
                     loc=LMT, **kwargs):
        #utcoffset=5.0, **kwargs):
        """plot uptimes for a catalog filename that is based on
        the LMT catalog format. If catalog filename is None,
        then only planets are plotted. If day is None, today's date is assumed.
        day can be a datetime.date instance or a string like 'May 21, 2011',
        or '2011-05-30', etc.
        If sources is a non-empty list, only objects in that list are plotted
        """
        from matplotlib.pyplot import setp

        self.clear()
        if catalog_fname is None:
            only_planets = True
        else:
            only_planets = False
            if not os.path.exists(catalog_fname):
                raise DreampyGeneralError("Plot Uptimes", "File %s does not exist" % catalog_fname)
            sdic = read_lmt_cat(catalog_fname, delimiter=delimiter)
            if sources:
                rsdic = {}
                for source in sources:
                    if sdic.has_key(source):
                        rsdic[source] = sdic[source]
                    else:
                        print("Source %s is not in catalog %s. Skipping.." % (source, catalog_fname))
                if not rsdic:
                    #empty
                    raise DreampyGeneralError("Plot Uptimes", "No sources in %s found in catalog %s" % (sources, catalog_fname))
                else:
                    print("Plotting only sources: %s" % rsdic.keys())
                sdic = rsdic
                        
        if day is not None:
            try:
                day = parser.parse(day)
            except:
                day = None
        if day is None:
            daydate = datetime.datetime.today()
        else:
            daydate = day
        if only_planets:
            localtime, uttime, selev = planetary_uptimes(day=day,
                                                         add_planets=add_planets,
                                                         loc=loc)
            #utcoffset=utcoffset)
        else:
            localtime, uttime, selev = source_uptimes(sdic, day=day,
                                                      add_planets=add_planets,
                                                      loc=loc)

        formatter = DateFormatter('%H:%M:%S') 
        for i, src in enumerate(selev.keys()):
            yval = (i+1)*2
            source_yval = yval*numpy.ones(len(localtime))
            ind = numpy.where(selev[src] < 20.)
            source_yval[ind] = numpy.nan
            idx = numpy.where(selev[src] == selev[src].max())
            self.plot_date(localtime, source_yval, '-', label=src)
            self.set_text(localtime[idx[0][0]], yval+0.5, src, fontdict={'size': "xx-small"})

        ax,kw = self.plotobj._get_current_axes(**kwargs)
        ax.xaxis.set_major_formatter(formatter) 
        labels = ax.get_xticklabels() 
        setp(labels, rotation=30, fontsize=8) 
        ymax = (((len(selev.keys())+1)*2)/10)*10 + 10
        self.set_ylim(0, ymax)
        #self.plotobj.f.autofmt_xdate()
        if only_planets:
            self.set_subplot_title('Planet Uptimes at %s for Date: %s' % (loc.name, daydate.strftime('%m-%d-%y')))
        else:
            self.set_subplot_title('Source Uptimes at %s for %s, Date: %s' % (loc.name, catalog_fname, daydate.strftime('%m-%d-%y')))
        self.set_xlabel('Localtime')
        self.redraw_plot()

    def plot_planet_uptimes(self, day=None, add_planets=True,
                            loc=LMT, **kwargs):
        self.plot_uptimes(catalog_fname=None, day=day, add_planets=add_planets,
                          loc=loc)


    def plot_horizon(self, catalog_fname=None, delimiter=',', day=None,
                     add_planets=True, sources=[],
                     loc=LMT, **kwargs):
        #utcoffset=5.0, **kwargs):
        """plot horizon plots for a catalog filename that is based on
        the LMT catalog format. If catalog filename is None,
        then only planets are plotted. If day is None, today's date is assumed.
        day can be a datetime.date instance or a string like 'May 21, 2011',
        or '2011-05-30', etc.
        If sources is a non-empty list, only objects in that list are plotted
        """
        from matplotlib.pyplot import setp
        self.clear()
        if catalog_fname is None:
            only_planets = True
        else:
            only_planets = False
            if not os.path.exists(catalog_fname):
                raise DreampyGeneralError("Plot Uptimes", "File %s does not exist" % catalog_fname)
            sdic = read_lmt_cat(catalog_fname, delimiter=delimiter)
            if sources:
                rsdic = {}
                for source in sources:
                    if sdic.has_key(source):
                        rsdic[source] = sdic[source]
                    else:
                        print("Source %s is not in catalog %s. Skipping.." % (source, catalog_fname))
                if not rsdic:
                    #empty
                    raise DreampyGeneralError("Plot Uptimes", "No sources in %s found in catalog %s" % (sources, catalog_fname))
                else:
                    print("Plotting only sources: %s" % rsdic.keys())
                sdic = rsdic

        if day is not None:
            try:
                day = parser.parse(day)
            except:
                day = None
        if day is None:
            daydate = datetime.datetime.today()
        else:
            daydate = day
        if only_planets:
            localtime, uttime, selev = planetary_uptimes(day=day,
                                                         add_planets=add_planets,
                                                         loc=loc)
        else:
            localtime, uttime, selev = source_uptimes(sdic, day=day,
                                                      add_planets=add_planets,
                                                      loc=loc)
        formatter = DateFormatter('%H:%M:%S') 
        for i, src in enumerate(selev.keys()):
            idx = numpy.where(selev[src] == selev[src].max())
            self.plot_date(localtime, selev[src], '-', label=src)
            self.set_text(localtime[idx[0][0]], selev[src].max()+0.5, src, fontdict={'size': "xx-small", 'horizontalalignment': 'center'})

        ax,kw = self.plotobj._get_current_axes(**kwargs)
        ax.xaxis.set_major_formatter(formatter) 
        labels = ax.get_xticklabels() 
        setp(labels, rotation=30, fontsize=8) 
        self.set_ylim(0, 90)
        #self.plotobj.f.autofmt_xdate()
        if only_planets:
            self.set_subplot_title('Horizon Plot at %s for Planets for Date: %s' % (loc.name, daydate.strftime('%m-%d-%y')))
        else:
            self.set_subplot_title('Horizon Plot at %s for Sources in %s, Date: %s' % (loc.name, catalog_fname, daydate.strftime('%m-%d-%y')))
        self.set_xlabel('Localtime')
        self.set_ylabel('Elev (deg)')
        self.redraw_plot()

    def plot_planet_horizon(self, day=None, add_planets=True,
                            loc=LMT, **kwargs):
        self.plot_horizon(catalog_fname=None, day=day, add_planets=add_planets,
                          loc=loc)
        
                     
    def plot_date(self, *args, **kwargs):
        self.check_alive()
        return self.plotobj.plot_date(*args, **kwargs)
    
    def hide_lines(self, lines):
        self.check_alive()
        self.plotobj.hide_lines(lines)

    def show_lines(self, lines):
        self.check_alive()
        self.plotobj.show_lines(lines)
        
    def refresh_canvas(self):
        """Clears figure, can also use L{clf} or L{clear}"""
        self.check_alive()        
        self.plotobj.refresh_canvas()

    def clf(self):
        """Clears figure, can also use L{clear} or L{refresh_canvas}"""
        self.check_alive()        
        self.refresh_canvas()

    def clear(self):
        """Clears figure, can also use L{clf} or L{refresh_canvas}"""
        self.check_alive()        
        self.refresh_canvas()

    def redraw_plot(self):
        self.check_alive()        
        self.plotobj.redraw_plot()
        
    def set_figtitle(self, title, *args, **kwargs):
        """Adds a title 'title' to the Figure. Text features
        can be added to kwargs, see matplotlib documentation on U{suptitle<http://matplotlib.sourceforge.net/api/figure_api.html?highlight=suptitle#matplotlib.figure.Figure.suptitle>}"""
        self.check_alive()        
        self.plotobj.set_figtitle(title, *args, **kwargs)

    def set_subplot_title(self, title, **kwargs):
        """Adds a title 'title' to the subplot.
        Text features can be added to kwargs, see matplotlib documentation
        on U{suptitle<http://matplotlib.sourceforge.net/api/pyplot_api.html?highlight=suptitle#matplotlib.pyplot.suptitle>}"""
        self.check_alive()        
        self.plotobj.set_subplot_title(title, **kwargs)

    def set_xlabel(self, xlabel, fontdict=None, labelpad=None, **kwargs):
        """Creates a xlabel 'xlabel'. Fondict is an optional dictionary
        which can override default text properties,
        fontdict=None sets to default properties.
        labelpad is the spacing from the xlabel to the x-axis in points. """
        self.check_alive()        
        self.plotobj.set_xlabel(xlabel, fontdict=fontdict,
                                labelpad=labelpad, **kwargs)

    def set_ylabel(self, ylabel, fontdict=None, labelpad=None, **kwargs):
        """Creates a ylabel 'ylabel'. Fondict is an optional dictionary
        which can override default text properties,
        fontdict=None sets to default properties.
        labelpad is the spacing from the ylabel to the y-axis in points. """
        self.check_alive()
        self.plotobj.set_ylabel(ylabel, fontdict=fontdict,
                                labelpad=labelpad, **kwargs)

    def annotate(self, label, xy, xytext=None, xycoords='data',
                 textcoords='data', arrowprops=None, **kwargs):
        """Annotates the point xy with text s at xytext.  If xytext is set
        to None the text appears at the point.
        arrowprops=None is the default properties set for the arrow
        which will connect the text to the data point"""
        self.check_alive()        
        self.plotobj.annotate(label, xy, xytext=xytext,
                              xycoords=xycoords, textcoords=textcoords,
                              arrowprops=arrowprops, **kwargs)

    def set_text(self, x, y, s, fontdict=None, **kwargs):
        """Adds text s at a position (x,y)
        fontdict is an optional dictionary in which you can override default
        text properties,
        fontdict=None set to default properties"""
        self.check_alive()        
        self.plotobj.set_text(x, y, s, fontdict=fontdict, **kwargs)

    def grid(self, b=None, **kwargs):
        """Adds a grid to subplot. b=None is a boolean expression,
        setting b=False for example would create the grid but not show it.
        It is assumed if b=None that the grid is to be shown"""
        self.check_alive()        
        self.plotobj.grid(b=b, **kwargs)

    def set_legend(self, *args, **kwargs):
        """Adds a legend to a subplot.  When using plot() it is useful to say
        label='label1' and so on. Adding loc=number(1-10) provides a convenient
        spot for the legend, find documentation for U{legend<http://matplotlib.sourceforge.net/api/axes_api.html?highlight=legend#matplotlib.axes.Axes.legend>} to see number locatins"""
        self.check_alive()        
        self.plotobj.set_legend(*args, **kwargs)

    def imshow(self, X, **kwargs):
        """Displays an image X to the current axes. X may be an image,
        a two- dimensional array etc.  Everything is optional from there
        but can be useful see documentaion on matplotlib U{imshow<http://matplotlib.sourceforge.net/api/axes_api.html?highlight=imshow#matplotlib.axes.Axes.imshow>}"""
        self.check_alive()
        return self.plotobj.imshow(X, **kwargs)

    def colorbar(self, **kwargs):
        """Displays colorbar for an image in current axis.
        current axis can also be passed using current_axes_name keyword"""
        self.check_alive()        
        self.plotobj.colorbar(**kwargs)
        
    def get_current_axis(self):
        """Returns the name of the current axis"""
        self.check_alive()
        return self.plotobj.get_current_axis()

    def hist(self,x, bins=10, range=None, normed=False,
             weights=None, cumulative=False, bottom=None,
             histtype='bar', align='mid', orientation='vertical',
             rwidth=None, log=False, **kwargs):
        """Takes a dataset x and creates a histogram of the data.
        Bins is the number of stacks the data will fit in.
        All other kwargs are optional but useful, see matplotlib
        documentation on U{hist<http://matplotlib.sourceforge.net/api/axes_api.html?highlight=hist#matplotlib.axes.Axes.hist>}"""
        self.check_alive()
        return self.plotobj.hist(x, bins=bins, range=range,
                                 normed=normed, weights=weights,
                                 cumulative=cumulative, bottom=bottom,
                                 histtype=histtype, align=align,
                                 orientation=orientation, rwidth=rwidth,
                                 log=log, **kwargs)

    def set_axis_off(self,**kwargs):
        """Turns axis off but keeps any plots,hists etc."""
        self.check_alive()
        self.plotobj.set_axis_off()

    def set_axis_on(self,**kwargs):
        """Turns axis on if off"""
        self.check_alive()
        self.plotobj.set_axis_on()

    def set_xlim(self,xmin=None, xmax=None, emit=True, **kwargs):
        """Sets limits on x-axis. emit=True notifys observer of limit changes"""
        self.check_alive()
        self.plotobj.set_xlim(xmin=xmin, xmax=xmax, emit=emit, **kwargs)
        self.redraw_plot()
        
    def set_ylim(self,ymin=None, ymax=None, emit=True, **kwargs):
        """Sets limits on x-axis. emit=True notifys observer of limit changes"""
        self.check_alive()
        self.plotobj.set_ylim(ymin=ymin, ymax=ymax, emit=emit, **kwargs)
        self.redraw_plot()
        
    def set_xscale(self,value, **kwargs):
        """Sets xscale to something other than linear such as logarithmic.  value can be 'linear' or  'log' or  'symlog'. see matplotlib documentation on U{set_xscale<http://matplotlib.sourceforge.net/api/axes_api.html?highlight=set_xscale#matplotlib.axes.Axes.set_xscale>} for usefull kwargs"""
        self.check_alive()
        self.plotobj.set_xscale(value, **kwargs)

    def set_yscale(self,value, **kwargs):
        """Sets yscale to something other than linear such as logarithmic.  value can be 'linear' or  'log' or  'symlog'. see matplotlib documentation on U{set_xscale<http://matplotlib.sourceforge.net/api/axes_api.html?highlight=set_xscale#matplotlib.axes.Axes.set_xscale>} for usefull kwargs"""
        self.check_alive()        
        self.plotobj.set_yscale(value, **kwargs)

    def set_xticks(self,ticks, minor=False, **kwargs):
        """Accepts an array or list of ticks and overrides default ticks"""
        self.check_alive()
        self.plotobj.set_xticks(ticks, minor=minor)

    def set_yticks(self,ticks, minor=False, **kwargs):
        """Accepts an array or list of ticks and overrides default ticks"""
        self.check_alive()
        self.plotobj.set_yticks(ticks, minor=minor, **kwargs)

    def stem(self,x, y, linefmt='b-', markerfmt='bo', basefmt='r-', **kwargs):
        """Creates a stem plot for x and y """
        self.check_alive()
        self.plotobj.stem(x, y, linefmt=linefmt,
                          markerfmt=markerfmt, basefmt=basefmt,
                          **kwargs)

    def set_position(self, pos, which='both', **kwargs):
        """Sets the axis position with a pos list [left, bottom, width, height] in relative coordinates you would use to create a subplot"""
        self.check_alive()
        self.plotobj.set_position(pos, which=which,**kwargs)

    def set_axis_bgcolor(self, color, **kwargs):
        """Sets the current axis' background color"""
        self.check_alive()
        self.plotobj.set_axis_bgcolor(color, **kwargs)

    def semilogx(self, *args, **kwargs):
        """Sets semilog scale of x"""
        self.check_alive()
        self.plotobj.semilogx(*args, **kwargs)

    def semilogy(self, *args, **kwargs):
        """Sets semilog scale of y"""
        self.check_alive()
        self.plotobj.semilogy(*args, **kwargs)

    def scatter(self,x, y, s=20, c='b', marker='o', cmap=None,
                norm=None, vmin=None, vmax=None, alpha=1.0,
                linewidths=None, faceted=True, verts=None, **kwargs):
        """Creates a scatter plot of xdata versus ydata. s is size in points^2, c is a color, and different markers can be found in the matplotlib documentation of U{scatter<http://matplotlib.sourceforge.net/api/axes_api.html?highlight=scatter#matplotlib.axes.Axes.scatter>} """
        self.check_alive()        
        self.plotobj.scatter(x, y, s=s, c=c, marker=marker,
                             cmap=cmap, norm=norm, vmin=vmin,
                             vmax=vmax, alpha=alpha,
                             linewidths=linewidths, faceted=faceted,
                             verts=verts, **kwargs)

    def pie(self,x, explode=None, labels=None, colors=None,
            autopct=None, pctdistance=0.6, shadow=False,
            labeldistance=1.1):
        """Creates a pie chart of an array x, the area of each wedge is x/sum(x).  See matplotlib documentation on U{pie<http://matplotlib.sourceforge.net/api/axes_api.html?highlight=pie#matplotlib.axes.Axes.pie>} for optional parameters"""
        self.check_alive()
        self.plotobj.pie(x, explode=explode, labels=labels,
                         colors=colors, autopct=autopct,
                         pctdistance=pctdistance, shadow=sahdow,
                         labeldistance=labeldistance)
        
    def minorticks_on(self, **kwargs):
        """Creates autoscaled minor ticks to the axes"""
        self.check_alive()
        self.plotobj.minorticks_on(**kwargs)

    def minorticks_off(self, **kwargs):
        """Removes autoscaled minor ticks to the axes"""
        self.check_alive()
        self.plotobj.minorticks_off(**kwargs)

    def invert_xaxis(self,**kwargs):
        """Inverts the x-axis"""
        self.check_alive()
        self.plotobj.invert_xaxis(**kwargs)

    def invert_yaxis(self, **kwargs):
        """Inverts the y-axis"""
        self.check_alive()
        self.plotobj.invert_yaxis(**kwargs)

    def hlines(self,y, xmin, xmax, colors='k',
               linestyles='solid', label='', **kwargs):
        """Plots horizontal lines at each y from xmin to xmax """
        self.check_alive()
        self.plotobj.hlines(y, xmin, xmax, colors=colors,
                            linestyles=linestyles, label=label, **kwargs)

    def fill(self, *args, **kwargs):
        """Fills regions of x arrays and y arrays specified"""
        self.check_alive()
        self.plotobj.fill(*args, **kwargs)

    def errorbar(self, x, y, yerr=None, xerr=None, fmt='-',
                 ecolor=None, elinewidth=None, capsize=3,
                 barsabove=False, lolims=False, uplims=False,
                 xlolims=False, xuplims=False, **kwargs):
        """Plots x versus y with error deltas in yerr and xerr. Vertical errorbars are plotted if yerr is not None. Horizontal errorbars are plotted if xerr is not None.  See matplotlib documentation on U{errorbar's<http://matplotlib.sourceforge.net/api/axes_api.html?highlight=error#matplotlib.axes.Axes.errorbar>} optional paramaters"""
        self.check_alive()
        self.plotobj.errorbar(x, y, yerr=yerr, xerr=xerr, fmt=fmt,
                              ecolor=ecolor, elinewidth=elinewidth,
                              capsize=capsize, barsabove=barsabove,
                              lolims=lolims, uplims=uplims,
                              xlolims=xlolims, xuplims=xuplims, **kwargs)

    def clear_axis(self, **kwargs):
        """Clears current or specified axis"""
        self.check_alive()
        self.plotobj.clear_axis(**kwargs)
        
    def boxplot(self, x, notch=0, sym='b+', vert=1,
                whis=1.5, positions=None, widths=None):
        """Creates a boxplot of the data in x. Notch = 0 (default) produces a rectangular box plot, notch = 1 will produce a notched box plot.  Vert=1 (default) makes a vertical boxplot vert=0 makes a horizontal boxplot.  See matplotlib documentation on U{boxplot<http://matplotlib.sourceforge.net/api/axes_api.html?highlight=boxplot#matplotlib.axes.Axes.boxplot>} for optional parameters"""
        self.check_alive()
        self.plotobj.boxplot(x, notch=notch, sym=sym, vert=vert,
                             whis=whis, positions=positions, widths=widths)

    def bar(self,left, height, width=0.8, bottom=None, color=None,
            edgecolor=None, linewidth=None, yerr=None, xerr=None,
            ecolor=None, capsize=3, align='edge',
            orientation='vertical', log=False, **kwargs):
        """Creates a bar chart of the data in x, see matplotlib documentation on bar for details on optional parameters"""
        self.check_alive()
        self.plotobj.bar(left, height, width=width, bottom=bottom,
                         color=color, edgecolor=edgecolor,
                         linewidth=linewidth, yerr=yerr, xerr=xerr,
                         ecolor=ecolor, capsize=capsize, align=align,
                         orientation=orientation, log=log, **kwargs)

    def vline(self,x=0, ymin=0, ymax=1, **kwargs):
        """Draws a vertical line at x from ymin to ymax """
        self.check_alive()
        return self.plotobj.vline(x=x, ymin=ymin, ymax=ymax, **kwargs)


    def contour(self, *args, **kwargs):
        """Contour plot of data"""
        self.check_alive()
        self.plotobj.contour(*args, **kwargs)

    def contourf(self, *args, **kwargs):
        """Filled contour plot of data"""
        self.check_alive()
        self.plotobj.contourf(*args, **kwargs)

    def get_xlims(self, **kwargs):
        return self.plotobj._get_xlims(**kwargs)

    def get_ylims(self, **kwargs):
        return self.plotobj._get_ylims(**kwargs)

    def print_figure(self, format='PNG'):
        return self.plotobj.print_figure(format=format)
    
class DreamPlot(DreamPlotBase):  #, gtk.Window):
    """
    The base class for dreampy interactive plotting.
    This class provides a handy plotting tool based on
    python's matplotlib library.

    >>> from dreampy.plots import DreamPlot
    >>> pl = DreamPlot()
    >>> x = numpy.arange(10)
    >>> y = x**2
    >>> pl.plot(x, y, 'bo-')
    
    """
    def __init__(self, figsize=(10, 8), **kwargs):
        super().__init__()
        #gobject.threads_init()
        #gtk.Window.__init__(self)
        #self.set_title("Dreampy Plot Window")
        #x, y = dreampyParams['plot']['figsize']
        #self.set_default_size(x, y)
        #self.connect("destroy", lambda x: gtk.main_quit())
        #self.connect("delete_event", self.delete_event)
        #self.vbox = gtk.VBox(False, 0)
        #self.add(self.vbox)
        self.plotobj = MainView(figsize=figsize, **kwargs)
        #self.add(self.plotobj)
        #self.alive = True
        #self.show_all()
        #self.gtk_catchup()

    # def gtk_catchup(self):
    #     if ENABLE_GTK:
    #         from IPython.lib.inputhook import enable_gui
    #         enable_gui(gui='gtk')
        
    def delete_event(self, widget, event, data=None):
        #return True means will not generate destroy event next
        #print "Type Ctrl-D or exit() to exit from main terminal window"
        logger.info("Killing plot window")
        self.alive = False
        return False
    
    def disconnect_event(self, sig):
        self.plotobj.disconnect_event(sig)

    def connect_event(self, eventname, eventfunc):
        return self.plotobj.connect_event(eventname, eventfunc)

class DreamPlotChart(DreamPlotBase):
    """
    A chart widget for dreampy for non-interactive plots
    """
    def __init__(self, chart_prop=None, **kwargs):
        if chart_prop is None:
            chart_prop = ChartProperties()
        chart_prop.get_properties(kwargs)
        self.plotobj = ChartView(chart_prop, **kwargs)
        self.alive = True


    def savefig(self, filename, dpi=None, facecolor='w',
                edgecolor='w', orientation='portrait',
                format='png', transparent=False, bbox_inches=None,
                pad_inches=0.1):
        """
        saves the figure into a file
        """
        logger.info("Saving figure to file %s" % filename)
        self.plotobj.savefig(filename, dpi=dpi, facecolor=facecolor,
                             edgecolor=edgecolor, orientation=orientation,
                             format=format, transparent=transparent,
                             bbox_inches=bbox_inches, pad_inches=pad_inches)
