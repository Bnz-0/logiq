#### Qerrors.py
#
# This file contains all type of custom exception that can be raised in logiq
#
####


class_name = lambda obj: '<'+obj.__class__.__name__+'>'

class LogiqError(Exception):

    def __init__(self, mex, prev=None):
        self.mex = mex
        self.prev = prev
    

    def __str__(self):
        return self.mex + ('' if self.prev is None else ', raised by ' + class_name(self.prev) + ': ' + str(self.prev))


class GenericLogiqError(LogiqError):
    "Generically raised when something goes wrong"
    pass

class InitializationError(LogiqError):
    "Raised when an object failed to initialize itself"
    pass

class DimensionError(LogiqError):
    "Raised when a linear object's dimension isn't compatible"
    pass

class NotAllowError(LogiqError):
    "Raised when the operation is not allowed"
    pass

class IllegalOperationError(LogiqError):
    "Raised when the simulator system can't do something"
    pass

class IncomprehensibleStatusError(LogiqError):
    "Raised when a status of a quantum state, described in some way, contain a \"sintax\" error"
    pass
