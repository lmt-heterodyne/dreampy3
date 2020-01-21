"""Utilities to read a Redshift ASCII file
Primarily used to read Ascii files generated by
redgui in the lab"""

import re
import numpy
from dreampy3.redshift.utils import LMTRedshiftError
import os

def countLines(file):
    """Count number of lines in file"""
    cur_offset = file.tell()
    numlines = len(file.readlines())
    file.seek(cur_offset)
    return numlines

def countChars(file):
    """Count Number of chars in file"""
    cur_offset = file.tell()        # remmember current position
    file.seek(0, 2)                 # seek to end of file
    file_len = file.tell()          # the end-position (number of bytes)
    file.seek(cur_offset)           # seek back to where we were
    return file_len

def parseheader(header):
    """Parse 17 lines in header of accum file and return a header
    dictionary"""
    import re
    headerdict = {}
    dashes = re.compile('^-+')
    dictpat = re.compile('\s*(\S+\s*\S+)\s*[=:]\s*(\S+)')
    datepat = re.compile('^Mon|^Tue|^Wed|^Thu|^Fri|^Sat|^Sun')
    for line in header.splitlines():
        if not dashes.match(line):
            if datepat.match(line):
                headerdict['date'] = line
            mat = re.split(',|;', line)
            if mat:
                for m in mat:
                    res = dictpat.match(m)
                    if res:
                        (key, val) = res.groups()
                        key = re.sub('Hardware Version', 'Hardware_Version', key)
                        key = re.sub('Software Version', 'Software_Version', key)
                        headerdict[key] = val
    
    (jnk,buf,sbuf, jnk) = headerdict['Errormes'].split('_')
    if not headerdict.has_key('Probemes'):
        headerdict['Probemes'] = ''
    headerdict['buffer'] = buf
    headerdict['subbuffer'] = sbuf
    return headerdict

def ReadAccumFiletoDic(filename, numhead=17, numchan=512):
    """Read accum.dat and return dictionary of header, chan and data arrays
    The dictionary keys are chassis first. Each chassis can have six board keys
    """
    if not os.path.exists(filename):
        raise LMTRedshiftError("ReadAccumFiletoDic",
                               "File %s does not exist" % filename)        
    try:
        acc = open(filename, "r")
    except:
        raise LMTRedshiftError("ReadAccumFiletoDic",
                               "File IO error: %s" % filename)
    TotalLines = countLines(acc)
    numlines_per_board = numchan+numhead
    numboards = TotalLines/numlines_per_board
    headers = {}
    chans = {}
    datas = {}
    for boards in range(numboards):
        head = ''
        for line in range(numhead):
            head += acc.readline()
        headerdic = parseheader(head)
        chassis, board = int(headerdic.get('chassis',0)), int(headerdic.get('board',0))
        if headers.has_key(chassis):
            headers[chassis][board] = headerdic
        else:
            headers[chassis] = { board: headerdic}
        ch = []
        da = []
        for line in range(numchan):
            (c,d) = acc.readline().split()
            ch.append(float(c))
            da.append(float(d))
        ch = numpy.array(ch)
        da = numpy.array(da)
        if chans.has_key(chassis):
            chans[chassis][board] = ch
        else:
            chans[chassis] = {board: ch}
        if datas.has_key(chassis):
            datas[chassis][board] = da
        else:
            datas[chassis] = {board: da}
    acc.close()
    return (headers, chans, datas)
