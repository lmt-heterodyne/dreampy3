"""
Read a default configuration file for default values
and save it upon each exit. Next time we fire up, use those
default values.

This uses python's configobj library.
"""

from configobj import ConfigObj,flatten_errors
from validate import Validator

import os
#import types

from dreampy3.logging import logger
logger.name = __name__

config_spec_text = """
#dreampy general items
[dreampy]
#loglevel can be one of these. If loglevel is set to 20 (say)
#then only messages with a loglevel higher than 20 will be emitted
#default loglevel is 10
# List of loglevels:
# Name    Level
# DEBUG0 : 5
# DEBUG1 : 6
# DEBUG2 : 7
# DEBUG3 : 8
# DEBUG4 : 9
# DEBUG  : 10
# INFO   : 20
# WARNING: 30
# ERROR  : 40
# CRITICAL: 50
# fileloglevel is for dreampy.log in .dreampyrc
# consoleloglevel is for console printout
fileloglevel  = integer(5, 50, default=10)
consoleloglevel  = integer(5, 50, default=10)

#general plot specific items
[plot]
#figsize is figure size in pixels (x,y)
figsize = int_list(min=2, max=2, default=list(800, 600))
#more plot configs to follow later

[holography]
plot_panels = boolean(default=True)
plot_sub_panels = boolean(default=False)

[redshiftchassis]
#This object contains configuration items specific to redshift chassis
#each bad_lags item is chassis specific, each bad_lags item has six entries
#one for each board in the chassis
#Bad Lags: slash-separated list like '220/234'
bad_lags0 = string_list(min=6, max=6, default=list('', '', '', '', '', ''))
bad_lags1 = string_list(min=6, max=6, default=list('', '', '', '', '', ''))
bad_lags2 = string_list(min=6, max=6, default=list('', '', '', '', '', ''))
bad_lags3 = string_list(min=6, max=6, default=list('', '', '', '', '', ''))
"""

def default_config():
    cfg_spec = ConfigObj(config_spec_text.splitlines(), list_values=False)
    valid = Validator()
    cfg = ConfigObj(configspec=cfg_spec, stringify=True, list_values=True)
    test = cfg.validate(valid, copy=True)
    if test:
        return cfg

def validate_dictionary(cdic):
    """This function validates a dictionary against the config spec here"""
    cfg_spec = ConfigObj(config_spec_text.splitlines(), list_values=False)
    valid = Validator()
    cfg = ConfigObj(cdic, configspec=cfg_spec)
    rtn = cfg.validate(valid, preserve_errors=True)
    if isinstance(rtn, bool) and rtn:
        return True
    else:
        res = flatten_errors(cfg, rtn)
        errortxt = ''
        for row in res:
            errortxt += 'In Section %s, key %s has error: %s' % (row[0], row[1], row[2])
            logger.error(errortxt)
        return False
    
class Configuration:
    def __init__(self, filename):
        """Initializes a config file if does not exist. If exists, uses
        it to validate the file, and setup default initial parameters"""
        self.cfg_spec = ConfigObj(config_spec_text.splitlines(), list_values=False)
        self.cfg_filename = filename
        valid = Validator()
        if not os.path.exists(self.cfg_filename):
            #no .dreampyrc file found
            cfg = ConfigObj(configspec=self.cfg_spec, stringify=True, list_values=True)
            cfg.filename = self.cfg_filename
            test = cfg.validate(valid, copy=True)
            cfg.write()
        self.cfg = ConfigObj(self.cfg_filename, configspec=self.cfg_spec)
        rtn = self.cfg.validate(valid, preserve_errors=True)
        if isinstance(rtn, bool) and rtn:
            logger.info("Config file validated")
            self.tested = True
        else:
            self.tested = False
            res = flatten_errors(self.cfg, rtn)
            self.errortxt = ''
            for row in res:
                self.errortxt += 'In Section %s, key %s has error: %s' % (row[0], row[1], row[2])
            logger.error(self.errortxt)

    def save_config(self, new_config):
        """Writes the config file upon exiting the program"""
        self.cfg.update(new_config)
        self.cfg.filename = self.cfg_filename
        self.cfg.write()
