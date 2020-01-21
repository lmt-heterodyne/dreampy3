import os
import ephem
from matplotlib.font_manager import FontProperties
import datetime
import numpy
from dreampy3.utils import DreampyGeneralError
from dreampy3.observing.locations import location, LMT, GBT, Amherst

class Source(object):
    def __init__(self):
        pass

    def ephem_db(self):
        return self.__repr__()
    
    def __repr__(self):
        return "%s,f|M|F7,%s,%s,%s,%s" % (self.name, self.RA, self.Dec, self.VLSR, self.epoch)
        
# class location(ephem.Observer):
#     def __init__(self, name, loc_lat, loc_long, loc_elev): #, loc_tzone):
#         self.name = name
#         self.lat = loc_lat
#         self.long = loc_long
#         self.elev = loc_elev
#         self.epoch = '2000'
#         #self.tzone = loc_tzone        

# LMT = location(name='LMT', loc_lat='18.9859', loc_long='-97.31458',
#                loc_elev=4640)#, loc_tzone=Central)


class Planets(object):
    def __init__(self): #, utcoffset=6.0):
        self.planet_list = []
        #self.utcoffset = utcoffset
        for planet in ('Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter',
                       'Saturn', 'Uranus', 'Neptune'):
            setattr(self, planet, getattr(ephem, planet)())
            self.planet_list.append(planet)
            
    def get_position(self, planet, time=None, loc=LMT):
        self.utcoffset = abs(loc.tzone.stdoffset.total_seconds()/3600.)
        if planet.title() not in self.planet_list:
            print("Planet %s not found in list" % planet.title())
        src = getattr(self, planet.title())
        if time is None:
            time = datetime.datetime.utcnow()
        localtime = time - datetime.timedelta(hours=self.utcoffset)
        LMT.date = time
        src.compute(LMT)
        return "%s at UTC %s (localtime %s) is Az: %s, Elev: %s" % (planet.title(), time.strftime("%Y-%m-%d %H:%M:%S"), localtime.strftime("%Y-%m-%d %H:%M:%S"), src.az, src.alt)


def read_lmt_cat_old(filename):
    """Reads an LMT like catalog file
    and returns a source dictionary. Each element of the
    source dictionary is a Source instance"""
    if os.path.exists(filename):
       fp = open(filename, 'r')
       sdic = {}
       for i, line in enumerate(fp.readlines()):
           if i == 0:
               #skip first line
               continue
           if line[0] == '#' or line[0] == ' ':
               #skip empty lines or comments
               continue
           elements = line.strip().split(',')
           source = Source()
           #print elements
           source.name = elements[0].strip(' ')
           source.RA = elements[1].strip(' ')
           source.Dec = elements[2].strip(' ')
           source.RARef = elements[3].strip(' ')
           source.DecRef = elements[4].strip(' ')
           source.VLSR = elements[5].strip(' ')
           source.epoch = elements[6].strip(' ')
           sdic[source.name] = source.ephem_db()
       return sdic

def read_lmt_cat(filename, delimiter=','):
    """
    Given a filename where the first line is a descriptor
    of the columns for eg:
    NAME,        RA,         Dec,     RA_REF,     Dec_REF,    VLSR,  EPOCH
    this function will read all the lines and convert to XEphem
    format catalog format and return source dictionary
    """
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
    sdic = {}
    for line in fp.readlines():
        #print line
        if line[0] == '#' or line[0] == ' ':
            #comment line or empty line
            continue
        try:
            columns = [col.strip() for col in line.strip().split(delimiter)]
            source = Source()
            source.name = columns[colname_index['name']]
            source.RA = columns[colname_index['ra']]
            source.Dec = columns[colname_index['dec']]
            source.VLSR = 0.0
            source.epoch = columns[colname_index['epoch']]
            sdic[source.name] = source.ephem_db()
            numbodies += 1
        except:
            continue
    fp.close()
    print("Read %d bodies from catalog file %s" % (numbodies, filename))
    return sdic

def source_uptimes(sdic, day=None, add_planets=True,
                   loc=LMT):
    #utcoffset=6):
    """Given a source dictionary, and a day, returns a list of local
    times and UT Times, and elevation angles of source.
    day is given as datetime.date
    """
    utcoffset = abs(loc.tzone.stdoffset.total_seconds()/3600.)
    localtime = []
    uttime = []
    if day is None:
        day = datetime.datetime.today()
    for i in range(288):
        ltime = datetime.datetime(day.year, day.month, day.day, 0,0,0) + \
                datetime.timedelta(seconds=5*60*i)
        utc = ltime + datetime.timedelta(hours=utcoffset)
        localtime.append(ltime)
        uttime.append(utc)
    source_elev = {}
    for source in sdic.keys():
        src = ephem.readdb(sdic[source])
        source_elev[source] = []
        for utime in uttime:
            loc.date = utime
            src.compute(loc)
            source_elev[source].append(src.alt)
        source_elev[source] = numpy.degrees(source_elev[source])
    if add_planets:
        for planet in ('Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter',
                       'Saturn', 'Uranus', 'Neptune'):
            src = getattr(ephem, planet)()
            source_elev[planet] = []
            for utime in uttime:
                loc.date = utime
                src.compute(loc)
                source_elev[planet].append(src.alt)
            source_elev[planet] = numpy.degrees(source_elev[planet])
    return localtime, uttime, source_elev

def planetary_uptimes(day=None, add_planets=True,
                      loc=LMT):
    #utcoffset=6):
    utcoffset = abs(loc.tzone.stdoffset.total_seconds()/3600.)
    localtime = []
    uttime = []
    if day is None:
        day = datetime.datetime.today()
    for i in range(288):
        ltime = datetime.datetime(day.year, day.month, day.day, 0,0,0) + \
                datetime.timedelta(seconds=5*60*i)
        utc = ltime + datetime.timedelta(hours=utcoffset)
        localtime.append(ltime)
        uttime.append(utc)
    source_elev = {}
    for planet in ('Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter',
                   'Saturn', 'Uranus', 'Neptune'):
        src = getattr(ephem, planet)()
        source_elev[planet] = []
        for utime in uttime:
            loc.date = utime
            src.compute(loc)
            source_elev[planet].append(src.alt)
        source_elev[planet] = numpy.degrees(source_elev[planet])
    return localtime, uttime, source_elev
    
if __name__ == '__main__':
    local, utc, sdic = source_uptimes('XGAL_Line_Survey.cat')
    for i, src in enumerate(sdic.keys()):
        yval = (i+1)*2
        source_yval = yval*numpy.ones(len(local))
        ind = numpy.where(sdic[src] < 20.)
        source_yval[ind] = numpy.nan
        idx = numpy.where(sdic[src] == sdic[src].max())
        plot_date(local, source_yval, '-', label=src)
        text(local[idx[0]], yval+0.5, src, fontdict={'size': "xx-small"})
    ylim(0, 130)
    gcf().autofmt_xdate()
    draw()
