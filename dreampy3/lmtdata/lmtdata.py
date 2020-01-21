
class LMTData(dict):
    """The most generic LMT Data class. This class is not usually
    instantiated by the user. It is part of the initialization of
    opening an LMT NetCDF file."""
    def __init__(self, ncvariables):
        #self.ncvariables = ncvariables
        self.make_data(ncvariables)

    def make_data(self, ncvariables):
        """Make the actual data"""
        datas = [name for name in ncvariables.keys() if name.find('Data') == 0]
        for d in datas:
            data_type = d.split('.')[1]
            self[data_type] = []
        for data_type in self.keys():
            keys = [name.split('.')[-1] for name in ncvariables.keys() if name.find('Data.%s.' % data_type) == 0]
            for key in keys:
                try:
                    #self.__setattr__(key, ncvariables['Data.%s.%s' % (data_type, key)].get())
                    if hasattr(self, key):
                        # there  is already a data attribute of same name so let's
                        # change how we write this
                        #print "Key %s already present" % key
                        self.__setattr__("%s%s" % (data_type, key), ncvariables['Data.%s.%s' % (data_type, key)][:])
                    else:
                        #print "Adding %s to data" % key
                        self.__setattr__(key, ncvariables['Data.%s.%s' % (data_type, key)][:])
                    self[data_type].append(key)
                except:
                    self.__setattr__(key, None)
        
    def _make_data_keys(self, ncvariables):
        dnames = [name for name in ncvariables.keys() if name.find('Data') != -1]
        for dname in dnames:
            dtype =dname.split('.')[1]
            self[dtype] = {}
        for key in self.keys():
            for subdata in [n for n in ncvariables.keys() if n.startswith('Data.%s.' % key)]:
                sdata = subdata.split('.')[-1]
                self[key][sdata] = ncvariables[subdata]

    # def make_data(self, ncvariables):
    #     self._make_data_keys(ncvariables)
