from dreampy3.utils import DreampyGeneralError

class LMTIFProcError(DreampyGeneralError):
    def __init__(self, errorname, reason):
        super(LMTIFProcError, self).__init__(errorname, reason)
