'''
Base module for all of the exceptions classes used internally.

Created on 10 Dec 2013

@author: alex
'''


class PheException(Exception):
    '''
    This is the top level class that EVERYTHING must be derived from. In particular,
    this class contains an abstract property called 'phe_return_code'. This property
    must be implemented and the individual implementation will have it's own
    exit code. which will be propogated to the calling functions, if needs be.
    
    PheException must not be passed as is.
    '''

    def __init__(self, msg, cause, phe_return_code=255):
        '''
        Constructor
        '''
        super(Exception, self).__init__(msg)
        self._phe_return_code = phe_return_code
        self._cause = cause
        
    @property
    def phe_return_code(self):
        '''
        Read-only attribute that holds the return status that should be exited with.
        '''
        return self._phe_return_code
    
    @property
    def cause(self):
        '''
        Read-only attribute that indicates the root cause of the exception raised.
        '''
        return self._cause
        
class PheExternalError(PheException):
    '''
    Exception class designed to be raised when an external command/process
        fails. Instead of falling over quietly, this exception can be raised. The
        exception includes the message to be put into the logs and the cause of
        the exception. In this case, the cause should generally be subprocess.CallerProcessError.
        The particulars of the failed command can be found inside the cause.
        If the catcher of this exception choses to exit the code, 'phe_return_code'
        should be used to indicate the cause of it all.
    '''        
    def __init__(self, msg, cause):
        '''
        Constructor for the PheExternalError
        
        @param msg: Message to be displayed with the exception.
        @type msg: str.
        @param cause: Cause of this exception, usually subprocess.CalledProcessError.
        @type cause: class.
        '''
        super(PheExternalError, self).__init__(msg, cause, 55)

