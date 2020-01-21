import numpy
import ephem
from dreampy3.observing.timezone import Eastern, Central, Central2
from dreampy3.utils import DreampyGeneralError
import datetime
import math
import matplotlib.dates as mdates
from dateutil.parser import parse
import pytz
import os

class location(ephem.Observer):
    def __init__(self, name, loc_lat, loc_long, loc_elev, loc_tzone):
        self.name = name
        self.lat = loc_lat
        self.long = loc_long
        self.elev = loc_elev
        self.tzone = loc_tzone
        self.bodies = {}
        self._init_planets()
        
    def set_time_now(self):
        self.date = ephem.now()

    def set_time(self, dt):
        """
        dt should be a datetime object
        if not timezone aware we will assume that it is in
        the local timezone of the location and then convert to utc
        """
        if type(dt) != datetime.datetime:
            # use a parser
            try:
                dt = parse(dt)
            except ValueError:
                raise DreampyGeneralError("Malformed Datetime", "Cannot parse datetime %s" % dt)
        if dt.tzinfo is None:
            #make it localtime first
            dt = dt.replace(tzinfo=self.tzone)
        self.date = dt.astimezone(pytz.utc)

    def _init_planets(self):
        for planet in ('Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter',
                       'Saturn', 'Uranus', 'Neptune'):
            self.bodies[planet] = getattr(ephem, planet)()

    def read_catalog(self, filename, delimiter=',',
                     append=True):
        """
        Given a filename where the first line is a descriptor
        of the columns for eg:
         NAME,        RA,         Dec,     RA_REF,     Dec_REF,    VLSR,  EPOCH
        this function will read all the lines and convert to XEphem
        format catalog format, and append to bodies object if append=True
        If append=False, will create a new dictionary of bodies
        """
        if not append:
            self.bodies = {}
        if not os.path.exists(filename):
            raise DreampyGeneralError("No Such File", "Cannot find catalog filename %s" % filename)
        #now open and read first line
        fp = open(filename, 'r')
        first = fp.readline()
        colnames = [col.strip().lower() for col in first.strip().split(delimiter)]
        colname_index = {}
        for col in ('name', 'ra', 'dec', 'epoch'):
            try:
                colname_index[col] = colnames.index(col)
            except ValueError:
                raise DreampyGeneralError("No col %s" % col, "Catalog Filename %s does not contain column name %s; it is needed!" % (filename, col))
        numbodies = 0
        for line in fp.readlines():
            #print line
            if line[0] == '#' or line[0] == ' ':
                #comment line or empty line
                continue
            try:
                columns = [col.strip() for col in line.strip().split(delimiter)]
                name = columns[colname_index['name']]
                ra = columns[colname_index['ra']]
                dec = columns[colname_index['dec']]
                epoch = columns[colname_index['epoch']]
                # see http://www.mmto.org/obscats/edb.html for XEphem format
                body = ephem.readdb('%s,f|J,%s,%s,,%s' % (name, ra, dec, epoch))
                numbodies += 1
                self.bodies[name] = body
                print("Read in %s" % name)
            except:
                continue
        print("Read %d bodies from catalog file %s" % (numbodies, filename))
        fp.close()
        

    def report_body(self, name):
        """
        Given a source name which has already been read in using read_catalog
        report_body will produce a simple summary of current az and elevation
        and rise_time, set_time and transit_time all in localtime
        """
        if not self.bodies.has_key(name):
            raise DreampyGeneralError("No such Source", "The current catalog does not have source with name %s" % name)
        body = self.bodies[name]
        self.set_time_now()
        body.compute(self)
        if type(body) == ephem.FixedBody:
            rstr = 'Source: %s; Ra: %s, Dec: %s\n' % (name, body._ra, body._dec)
        else:
            rstr = 'Source: %s; Ra: %s, Dec: %s\n' % (name, body.ra, body.dec)
        rstr += 'Current time: %s\n' % (datetime.datetime.now())
        rstr += 'Current Az: %s, El: %s\n' % (body.az, body.alt)
        rstr += 'Prev Rise Time: %s\n' % (ephem.localtime(self.previous_rising(body)))
        rstr += 'Next Rise Time: %s\n' % (ephem.localtime(self.next_rising(body)))        
        rstr += 'Next Transit Time: %s\n' % (ephem.localtime(self.next_transit(body)))
        rstr += 'Next Set Time: %s\n' % (ephem.localtime(self.next_setting(body)))
        print(rstr)
            

Amherst = location(name='Amherst', loc_lat='42.38028', loc_long='-72.52361',
                   loc_elev=100, loc_tzone=Eastern)

LMT = location(name='LMT', loc_lat='18.9859', loc_long='-97.31458',
               loc_elev=4640, loc_tzone=Central)

LMT2 = location(name='LMT', loc_lat='18.9859', loc_long='-97.31458',
               loc_elev=4640, loc_tzone=Central2)

GBT = location(name='GBT', loc_lat='38.4331211', loc_long='-79.839835',
               loc_elev=806, loc_tzone=Eastern)
