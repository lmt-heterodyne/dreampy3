"""
This is the LMT DREAMPY data analysis package
"""
import os, sys, shutil, platform
import atexit
#if os.path.exists(os.path.join(os.environ["HOME"], '.dreampy')):
#first-time setup - don't run logger yet
from dreampy3.logging import logger
logger.name = __name__


#@atexit.register
#def _ext1():
#    print 'ext1'
    
#@atexit.register
#def _ext2():
#    print 'ext2'

revString = "$Rev: 278 $: Last Commit"

def version():
    print("Dreampy: %s" % revString)
    return revString

@atexit.register
def _save_defaults():
    logger.info('Saving dreampy defaults')
    if cfgobj:
        cfgobj.save_config(dreampyParams)


def first_time_setup():
    """
    Check to see if the user has setup their
    environment correctly
    """
    dreamdata = '/usr/share/dreampy'
    #if os.environ.has_key("DREAMDATA"):
    #    if os.path.exists(os.environ["DREAMDATA"]):
    #        dreamdata = os.environ["DREAMDATA"]
    # set up user space
    userdir = os.environ["HOME"]+os.path.sep+".dreampy"
    if not os.path.exists(userdir):
        print('First time DREAMPY use. Setting up ~/.dreampy')
        os.mkdir(userdir)
        #shutil.copyfile(dreamdata+"/data/ipythonrc-dreampy", userdir+"/ipythonrc-dreampy")
        for fname in ('ipythonrc-dreampy', 'ipy_user_conf.py'):
            shutil.copyfile(os.path.join(dreamdata, fname),
                            os.path.join(userdir, fname))
        f = file(userdir+"/ipythonrc", "w")
        f.close()
    else:
        if not os.path.exists(os.path.join(userdir, 'upgraded')):
            logger.info("Upgrading ipy_user_conf.py")
            shutil.copyfile(os.path.join(dreamdata, 'ipy_user_conf.py'),
                            os.path.join(userdir, 'ipy_user_conf.py'))
            f = file(os.path.join(userdir, 'upgraded'), "w")
            f.close()
        #from django.core.management import setup_environ, syncdb
        #import spadb.settings
        #loc = setup_environ(spadb.settings)
        #syncdb(interactive=False)

def corr_cal_setup():
    # Allow user defined data location
    dreamdata = '/usr/share/dreampy'
    default_corr_cal_dir = '/raw/rsr/cal'
    #if os.environ.has_key('CORR_CAL_DIR'):
    if 'CORR_CAL_DIR' in os.environ:
        corr_cal_dir = os.environ['CORR_CAL_DIR']
        if not os.path.exists(corr_cal_dir):
            logger.info("corr_cal_dir %s does not exist. Setting to current directory %s. Please check environment variable CORR_CAL_DIR" % (corr_cal_dir, os.curdir))
            if os.path.exists(default_corr_cal_dir):
                corr_cal_dir = default_corr_cal_dir
            else:
                corr_cal_dir = os.curdir
        else:
            logger.info("Setting corr_cal_dir to directory %s" % corr_cal_dir)
    else:
        if os.path.exists(default_corr_cal_dir):
            corr_cal_dir = default_corr_cal_dir
        else:    
            corr_cal_dir = os.curdir
        logger.info("Setting corr_cal_dir to current directory %s" % corr_cal_dir)

    # remove from namespace
    #if dreamdata:
    #    del dreamdata
    #del userdir, shutil, platform
    return corr_cal_dir

# workaround for ipython, which redirects this if banner=0 in ipythonrc
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__

#import dreampy
#from dreampy.plots import DreamPlot
from dreampy3.utils.configuration import Configuration, default_config, validate_dictionary

first_time_setup()
corr_cal_dir = corr_cal_setup()
version()

#from dreampy.holography.plots import HologPlot
#from dreampy.holography.holoIO import HoloFile
#from dreampy.utils import OrderedDict

# from spapy.spa import spa_tfile, SpaScan, SpaAverage, spa_tempfile
# from spadb.spaheaders.models import spafile, spaheader
# from spapy.spaplots import SpaPlot
# from spapy.selection import Selection
# from spadb.get_correlation_lines import default_species as species

#def plotgui(**kwargs):
#    return DreamPlot(**kwargs)

def _dreampy_fname():
    """
    Returns the path to the rc file
    Search order:

    * current working dir
    * environ var DREAMPYRC
    * HOME/.dreampy/dreampyrc
    """

    fname = os.path.join( os.getcwd(), '.dreampyrc')
    if os.path.exists(fname): return fname

    #if os.environ.has_key('DREAMPYRC'):
    if 'DREAMPYRC' in os.environ:
        path =  os.environ['DREAMPYRC']
        if os.path.exists(path):
            fname = os.path.join(path, '.dreampyrc')
            if os.path.exists(fname):
                return fname
    #if os.environ.has_key('HOME'):
    if 'HOME' in os.environ:
        home =  os.environ['HOME']
        fname = os.path.join(home, '.dreampy', 'dreampyrc')
        if os.path.exists(fname):
            return fname
    return None


dreampyDefaultParams = default_config()

def get_configobj():
    """Return the default configuration object in the rc file"""
    fname = _dreampy_fname()
    if fname is None or not os.path.exists(fname):
        logger.info("could not find rc file; creating file with defaults")
        home =  os.environ['HOME']
        fname = os.path.join(home, '.dreampy', 'dreampyrc')
    return Configuration(fname)


cfgobj = get_configobj()

def dreampy_get_configuration():
    """Returns default configuration values updated with values
    from the rc file"""
    params = dreampyDefaultParams.copy()
    params.update(cfgobj.cfg)
    return params

dreampyParams = dreampy_get_configuration()

def setParam(group, **kwargs):
    """
    Set the current dreampyParams.  Group is the grouping for the rc, eg
    for plot.figsize the group is 'plot', For redshift.chassis, the
    group is 'redshift', and so on.  kwargs is a list of attribute
    name/value pairs, eg

      setParam('plot', figsize=[1024, 768])

    sets the current rc params and is equivalent to

      dreampyParams['plot']['figsize'] = [1024, 768]

    Use dreamdefaults to restore the default rc params after changes.
    """
    #if not dreampyDefaultParams.has_key(group):
    if group not in dreampyDefaultParams:
        raise KeyError('Unrecognized group "%s"' %  group)
    newParams = {}
    newParams[group] = {}
    for k,v in kwargs.items():
        #if not dreampyDefaultParams[group].has_key(k):
        if k not in dreampyDefaultParams[group]:
            raise KeyError('Unrecognized key "%s" in group "%s"' % (k, group))
        newParams[group][k] = v
        if not validate_dictionary(newParams):
            #something wrong with this value remove
            newParams[group].pop(k)
    dreampyParams[group].update(newParams[group])
    #dreampyParams.update(newParams)

def resetdefaultParams():
    """
    Restore the default dreampy params - the ones that were created at
    dreampy load time
    """
    dreampyParams.update(dreampyDefaultParams)


fileloglevel = dreampyParams['dreampy'].get('fileloglevel', 10)
consoleloglevel = dreampyParams['dreampy'].get('consoleloglevel', 10)

logger.handlers[1].setLevel(fileloglevel)
logger.handlers[2].setLevel(consoleloglevel)

corr_cal = {}  # a persistent corr_cal dictionary. Key is filename, and value
               # of dictionary is corr_cal dictionary
# def search_unique_files(name):
#     """Search for all filenames in the database that contain the given
#     source name"""
#     sel = Selection(name=name)
#     sc = sel._filter_scans().distinct('fileid').values('fileid')
#     sc = [s['fileid'] for s in sc]
#     sf = spafile.objects.in_bulk(sc)
#     print "Filename Number_of_Scans"
#     for fileobj in sf.values():
#         print "%s %d" % (fileobj.filename, fileobj.spaheader_set.filter(name__icontains=name).count())

# def search_scans(name=None, scan1=1, scan2=31000, utmin=None,
#                  utmax=None, oper=None, otype=None, sfile=None):
#     """Search for scans in the database for sourcename = name"""
#     sel = Selection(name=name, scan1=scan1, scan2=scan2,
#                     utmin=utmin, utmax=utmax, oper=oper, otype=otype,
#                     sfile=sfile)
#     sc = sel._filter_scans().order_by('fileid', 'utd')
#     ct = sc.count()
#     print "Filename ScanNo SourceName Operator UTDate OType"
#     for s in sc:
#         print "%s %d %s %s %s %s" % (os.path.basename(s.fileid.filename), \
#                                   s.iscan, s.name, s.oper, s.utdate, s.otype_str())
#     print "%d scans in database for Source Name containing %s" % (ct, name)

__date__ = '$Date: 2015-03-20 03:15:17 -0400 (Fri, 20 Mar 2015) $'.split()[1]
__version__  = '$Rev: 278 $'.split()[1]

def is_ipython():
    return '__IP' in dir(sys.modules["__main__"])

if is_ipython():
    def version(): logger.info("dreampy %s(%s)"% (__version__, __date__))

def welcome():
    return """Welcome to dreampy v%s (%s) - the LMT Data Reduction and Analysis Package

Please report any bugs to:
Gopal Narayanan <gopal@astro.umass.edu>

""" % (__version__, __date__)

def myexit():
    print("Exiting dreampy")
    raise SystemExit



#_ip.IP.ask_exit = myexit

if __name__ == '__main__':
    first_time_setup()
