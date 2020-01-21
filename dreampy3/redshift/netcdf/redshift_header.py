"""The RedshiftHeader class is an extension
of the LMTHeader class"""

#from lmtdcs.lmtheader import LMTHeader

from dreampy3.lmtheader import LMTHeader
from dreampy3.utils.coordinates import sixty
import numpy
from netCDF4 import num2date, date2num
import datetime

switch_settings = { 0: 'Acf',
                    1: 'Spectral',
                    2: 'Blank'
                    }

class RedshiftHeader(LMTHeader):
    def __init__(self, ncvariables=None, dimensions=None):
        """Initializes a redshift header instance.
        @param ncvariables: A L{netCDF4} variables ordered dictionary
        @type ncvariables: netcdf4 variables type
        @param dimensions: Dimensions object from L{netCDF4}
        @type dimensions: netCDF4 dimensions type
        """
        LMTHeader.__init__(self, ncvariables=ncvariables,
                           dimensions=dimensions,
                           fromold=False)
        self.calculate_attributes()

    def _identify_chassis(self):
        self.chassis = None
        for key in self.keys():
            if key.find('RedshiftChassis') != -1:
                chas_string = key
        if chas_string:
            self.chassis = int(chas_string.split('_')[1])
            
    def obs_string(self):
        return "%d.%d.%d" % (int(self.ObsNum), int(self.SubObsNum),
                             int(self.ScanNum))
    
    def _get_Spectrum_attributes(self):
        if 'Spectrum' in self:
            for key in self['Spectrum'].keys():
                setattr(self, key, self.get('Spectrum.%s' % key))

    def calculate_attributes(self):
       """Given a conventional LMTHeader object,
        creates redshift specific header attributes"""
       self._identify_chassis()
       chas_str = 'RedshiftChassis_%d_' % self.chassis
       self.SwtA = self.get(chas_str+'.SwtA', None)
       if self.SwtA is None:
           self.SwtA = self.get(chas_str+'.swta')
       self.SwtB = self.get(chas_str+'.SwtB', None)
       if self.SwtB is None:
           self.SwtB = self.get(chas_str+'.swtb')
       self.Blank = self.get(chas_str+'.Blank', None)
       if self.Blank is None:
           self.Blank = self.get(chas_str+'.blank')           
       self.Polar = self.get(chas_str+'.Polar', None)
       if self.Polar is not None:
           self.PolarOut = numpy.zeros(self.Polar.size, dtype=self.Polar.dtype)
           self.PolarDemod = numpy.zeros(self.Polar.size, dtype=self.Polar.dtype)
           for i, polar in enumerate(self.Polar):
               self.PolarOut[i] = polar & 0x03
               self.PolarDemod[i] = (polar >> 2) & 0x03
       #if self.PolarOut is None:
       #    self.PolarOut = self.get(chas_str+'.polar_out')
       #self.PolarDemod = self.get(chas_str+'.PolarDemod', None)
       #if self.PolarDemod is None:
       #    self.PolarDemod = self.get(chas_str+'.polar_demod')
       self.BoardNumber = self.get(chas_str+'.BoardNumber', None)
       if self.BoardNumber is None:
           self.BoardNumber = self.get(chas_str+'.board')
       self.SlotNumber = self.get(chas_str+'.SlotNumber', None)
       if self.SlotNumber is None:
           self.SlotNumber = self.get(chas_str+'.slotno')
       self.BoardID = self.get(chas_str+'.DipBid', None)
       if self.BoardID is None:
           self.BoardID = self.get(chas_str+'.boardid')
       self.BandName = self.get(chas_str+'.BandName', None)
       self.Tamb = self.get(chas_str+'.Tamb', 290.)
       self.Bandwidth = self.get(chas_str+'.Bandwidth', None)
       self.NumPixels = self.get(chas_str+'.NumPixels', None)
       self.NumBands = self.get(chas_str+'.NumBands', None)
       self.NumChannels = self.get(chas_str+'.NumChannels', None)
       self.IfFreq = self.get(chas_str+'.IfFreq', None)
       self.IfConfig = self.get(chas_str+'.IfConfig', None)
       self.Gain = self.get(chas_str+'.Gain', None)
       for i, gain in enumerate(self.Gain):
           self.Gain[i] = self.Gain[i] & 0x0f
       self.RefAcc = self.get(chas_str+'.RefAcc', None)
       self.DiagCtl = self.get(chas_str+'.DiagCtl', None)
       self.DiagData = self.get(chas_str+'.DiagData', None)
       self.AccDesiredValue = self.get(chas_str+'.AccDesiredValue', None)
       self.RawCntValue = self.get(chas_str+'.RawCntValue', None)
       self.DipBid = self.get(chas_str+'.DipBid', None)
       self.ChassisNumber = self.get(chas_str+'.ChassisNumber', None)
       self.HwVersion = self.get(chas_str+'.HwVersion', None)
       self.SwVersion = self.get(chas_str+'.SwVersion', None)
       self.DiagStatus = self.get(chas_str+'.DiagStatus', None)
       self.Valid = self.get(chas_str+'.Valid', None)
       self.Dumpmode = self.get(chas_str+'.Dumpmode', None)
       self.AcquireMode = self.get(chas_str+'.AcquireMode', None)
       self.AccumMode = self.get(chas_str+'.AccumMode', None)
       try:
           sw_settings = self.get(chas_str+'.SwitchSettings', [0])[0]
       except IndexError:
           # diff python versions return this differently
           sw_settings = self.get(chas_str+'.SwitchSettings', 0)
       self.SwitchSettings = switch_settings[sw_settings]
       self.TargetTime = self.get(chas_str+'.TargetTime', None)
       self.TargetCount = self.get(chas_str+'.TargetCount', None)
       self.ActualCount = self.get(chas_str+'.ActualCount', None)
       self.CorrCalID = self.get(chas_str+'.CorrCalId', None)
       CalObsNum = self.get(chas_str+'.CalObsNum', None)
       if CalObsNum is not None:
           self.CalObsNum = int(CalObsNum)
       else:
           self.CalObsNum = None
       self.SourceName = self.get('Source.SourceName', None)
       self.Epoch = self.get('Source.Epoch', None)
       self.Az = self.get('Source.Az', None)
       self.El = self.get('Source.El', None)
       AzReq = self.get('Sky.AzReq', None)
       ElReq = self.get('Sky.ElReq', None)
       if AzReq is not None:
           if AzReq.any():
               self.AzReq = numpy.degrees(AzReq[0])
       else:
           self.AzReq = None
       if ElReq is not None:
           if ElReq.any():
               self.ElReq = numpy.degrees(ElReq[0])
       else:
           self.ElReq = None
       self.Ra = self.get('Source.Ra', None)
       self.Dec = self.get('Source.Dec', None)
       if self.Ra is not None:
           if self.Ra.any():
               self.RA = sixty(self.Ra[0], ra=True)
       if self.Dec is not None:
           if self.Dec.any():
               self.DEC = sixty(self.Dec[0])
       self.L = self.get('Source.L', None)
       self.B = self.get('Source.B', None)
       self.CoordSys = self.get('Source.CoordSys', None)
       self.Velocity = self.get('Source.Velocity', None)
       self.VelSys = self.get('Source.VelSys', None)
       self.AzOff = self.get('Sky.AzOff', None)
       self.ElOff = self.get('Sky.ElOff', None)
       self.RaOff = self.get('Sky.RaOff', None)
       self.DecOff = self.get('Sky.DecOff', None)
       self.LOff = self.get('Sky.LOff', None)
       self.BOff = self.get('Sky.BOff', None)
       self.ObsVel = self.get('Sky.ObsVel', None)
       self.BaryVel = self.get('Sky.BaryVel', None)
       self.ParAng = self.get('Sky.ParAng', None)
       self.Operator = self.get('Telescope.Operator', None)
       self.AzTotalPoff = self.get('Telescope.AzTotalPoff', None)
       self.ElTotalPoff = self.get('Telescope.ElTotalPoff', None)
       self.AzPcor = self.get('Telescope.AzPcor', None)
       self.ElPcor = self.get('Telescope.ElPcor', None)
       self.XAct = self.get('M2.XAct', None)
       self.YAct = self.get('M2.YAct', None)
       self.ZAct = self.get('M2.ZAct', None)
       self.TipAct = self.get('M2.TipAct', None)
       self.TiltAct = self.get('M2.TiltAct', None)
       self.Wvp = self.get('Environment.Wvp', None)
       self.Temp = self.get('Environment.Temp', None)
       self.Pressure = self.get('Environment.Pressure', None)
       self.LST = self.get('TimePlace.LST', None)
       self.UTDate = self.get('TimePlace.UTDate', None)
       self.UT1 = self.get('TimePlace.UT1', None)
       self.ModRev = self.get('PointModel.ModRev', None)
       self.ObsPgm = self.get('Dcs.ObsPgm', None)
       self.Receiver = self.get('Dcs.Receiver', None)
       self.ObsNum = self.get('Dcs.ObsNum', None)
       self.SubObsNum = self.get('Dcs.SubObsNum', None)
       self.ScanNum = self.get('Dcs.ScanNum', None)
       self.ObsType = self.get('Dcs.ObsType', None)
       self.ObsMode = self.get('Dcs.ObsMode', None)
       self.CalMode = self.get('Dcs.CalMode', None)
       self.Timer = self.get('Dcs.Timer', None)
       self.Valid = self.get('ScanFile.Valid',None)
       self._get_Spectrum_attributes()
        

class RedshiftCorrCalHeader(RedshiftHeader):
    def __init__(self, ncvariables=None, dimensions=None):
        """Initializes a redshift corrcal header instance.
        @param ncvariables: A L{netCDF4} variables ordered dictionary
        @type ncvariables: netcdf4 variables type
        @param dimensions: Dimensions object from L{netCDF4}
        @type dimensions: netCDF4 dimensions type
        """
        RedshiftHeader.__init__(self, ncvariables=ncvariables,
                                dimensions=dimensions)

    def calculate_attributes(self):
       """Given a conventional LMTHeader object,
       creates redshift specific header attributes"""
       self._identify_chassis()
       chas_str = 'RedshiftChassis_%d_' % self.chassis
       self.zerolag = self.get(chas_str+'.zerolag', None)
       self.boardid = self.get(chas_str+'.boardid', None)
       self.svd = self.get(chas_str+'.svd', None)
       self.slotno = self.get(chas_str+'.slotno', None)
       self.corr_cal_id = self.get(chas_str+'.corr_cal_id', None)
       obsdate = self[chas_str]['obsdate']
       self.createdate = num2date(self.get(chas_str+'.createdate', 0), obsdate.units,
                                  obsdate.calendar)
       self.reducedate = num2date(self.get(chas_str+'.reducedate', 0), obsdate.units,
                                  obsdate.calendar)
       self.obsdate = num2date(self.get(chas_str+'.obsdate', 0), obsdate.units,
                               obsdate.calendar)
       
