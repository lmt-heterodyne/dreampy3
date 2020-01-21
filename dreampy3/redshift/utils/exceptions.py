from dreampy3.utils import DreampyGeneralError

class LMTRedshiftError(DreampyGeneralError):
    def __init__(self, errorname, reason):
        super(LMTRedshiftError, self).__init__(errorname, reason)
