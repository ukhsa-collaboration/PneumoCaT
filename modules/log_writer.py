import logging
import os
import time


def setup_logging(logs_dirname, time_stamp=None, project_basename=None, logger_name=None):
    '''
    setup_logging(logs_dirname, time_stamp=None, project_basename=None, logger_name=None)

    Wrapper to setup_logger which abstracts the details & provides for time-stamping.

    Parameters
    ----------
    logs_dirname : str
        Path to which logs/ dir will be appended
    time_stamp : None or str
        If None, a stamp-stamp in the format 'yyyy-mm-dd-hm' is generated, otherwise
        the timestamp passed is used to name the output log files.
    project_basename : None or str
        Basename of the calling project; used to name the logger unless a logger_name is passed
    logger_name : None or str
        Default is None, in which case the logger_name is project_basename_logger
        Otherwise the passed logger_name is used.

    Return
    ------
    logger : An instance of the logger class

    Note
    ----
    When project_basename is passed the logger_name is project_basename_logger &
    can be recovered elsewhere in the project using the statement
        logger = logging.getlogger( logger_name )
    otherwise the passed logger_name is used.

    When project_basename is passed the log files are written to log_dirname/logs/ as
    project_basename_timestamp.stdout & project_basename_timestamp.stderr
    otherwise the passed logger_name is used to write to files
    logger_name_timestamp.stdout & logger_name_timestamp.stdout 

    MGGoulden
    2014OCT01
    '''

    logs_path = os.path.join(logs_dirname, 'logs/')
    print 'logs_path:', logs_path
    if not os.path.exists(logs_path):
        os.makedirs(logs_path)

    if not time_stamp:
        time_stamp = time.strftime('%Y-%m-%d_%H%M%S')

    if not logger_name:
        logger_name = '_'.join([ project_basename, 'logger' ])
        stdout, stderr = [ ''.join([ logs_path, project_basename, '_', time_stamp, log ]) for log in ['.stdout', '.stderr'] ]
    else:
        stdout, stderr = [ ''.join([ logs_path, logger_name, '_', time_stamp, log ]) for log in ['.stdout', '.stderr'] ]

    print 'stdout:', stdout  # tmp
    print 'stderr:', stderr  # tmp
    print 'logger_name:', logger_name  # tmp

    logger = setup_logger(info_file=stdout, error_file=stderr, logger_name=logger_name)
    return logger




def setup_logger(info_file="stdout.log", error_file="stderr.log", logger_name='stdout_stderr_logger'):
    """
    A function to set up a logger for writing out log files

    Parameters
    ----------
    info_file : String
        Path to info level log file (default: stdout.log)
    error_file : String
        Path to error level log file (default: stderr.log)
    logger_name : String
        Name for logger (default: 'stdout_stderr_logger')

    Returns
    -------
    logger : A instance of the Logger class
    
    Note
    ----
    In a code context into which 'logging' has been imported, the logger_name
    parameter allows a logger instantiated elsewhere to be recovered using
    logger = logging.getLogger( 'stdout_stderr_logger' ), obviating the need
    to pass logger as an arg to any function.
     
    """

    class LogLevelFilter(object):
        def __init__(self, level):
            self.__level = level

        def filter(self, logRecord):
            return logRecord.levelno <= self.__level

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s\n%(message)s')

    handler_stderr = logging.FileHandler(error_file)
    handler_stderr.setLevel(logging.ERROR)
    handler_stderr.setFormatter(formatter)
    handler_stderr.addFilter(LogLevelFilter(logging.ERROR))
    logger.addHandler(handler_stderr)

    handler_stdout = logging.FileHandler(info_file)
    handler_stdout.setLevel(logging.INFO)
    handler_stdout.setFormatter(formatter)
    handler_stdout.addFilter(LogLevelFilter(logging.INFO))
    logger.addHandler(handler_stdout)

    return logger

def write_log(logger, log_text, log_level):
    """
    Writes text to a logger at a particular log level

    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    log_text : String
        The text to be written to the log file
    log_level : String
        The level of logging to which the text should be written (either 'info' or 'error')

    """

    if log_level == "error":
        logger.error(log_text)
    elif log_level == "info":
        logger.info(log_text)

def log_process(logger, process, log_info_to="info", log_error_to="error", limit_logging=0):
    """
    A function to log the output of a subprocess.Popen call 

    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    process : process pipe
        A process created by subprocess.Popen
    log_info_to: String
        The level at which to log info level logs into (default 'info')
    log_error_to: String
        The level at which to log error level logs into (default 'error')
    limit_logging: Integer
        Limit logging to either stdout (1) or stderr (2), log both if 0 (default 0)

    Returns
    -------
    stdout : String
        The stdout from the process
    stderr : String
        The stderr from the process
    """
    # if the process gives a exit status greater than 1, i.e. an genuine error, then log_error_to 'error'.
    if process.returncode > 0:
        log_error_to = "error"

    stdout = ""
    stderr = ""

    if limit_logging != 2:
        if not process.stdout == None:
            stdout = process.stdout.read()
        if len(stdout) > 0:
            write_log(logger, stdout, log_info_to)
    if limit_logging != 1:
        if not process.stderr == None:
            stderr = process.stderr.read()
        if len(stderr) > 0:
            write_log(logger, stderr, log_error_to)

    return stdout, stderr

def write_header_to_log(logger, header_text, log_level):
    """
    A utility function to write some header text bound by asterisks to the log
    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    header_text : String
        The header text to embed within asterisks
    log_level : String
        The level at which to log
    """
    write_log(logger, "******** " + header_text + " ********", log_level)

def error_header(logger, header_text):
    """
    A utility function to write some header text bound by asterisks to the error log_text
    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    header_text : String
        The header text to embed within asterisks
    """
    write_header_to_log(logger, header_text, "error")

def info_header(logger, header_text):
    """
    A utility function to write some header text bound by asterisks to the info log_text
    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    header_text : String
        The header text to embed within asterisks
    """
    write_header_to_log(logger, header_text, "info")


def get_logger_path(ARGS, dir_name='logs'):
    ''' 
    Returns a path for logger output according to params in the return from
    parser.parse_args(). By default returns output_dir/logs, else returns
    input_dir/logs, else fastqfile_dir/logs

    Args:
        ARGS, Namespace object : the return from parser.parse_args() 
        dir_name, string : terminal dir appended to derived path

    Returns:
        output_dir, string : fully-specified output path terminating
        with /dir_name to which logs can be written

    Side effect:
        creates path if necessary   

    MGGoulden 20130709
    amended 20130927 & moved into log_writer module
    amended 20140109 to make more generic

    '''

    log_path = None
    # if an output_dir has been passed check it & use it
    if not ARGS.output_dir is None:
        if not os.path.isdir(ARGS.output_dir):
            print 'the output_dir passed (' + ARGS.output_dir + ') does not exist'
            print 'making dir: ' + ARGS.output_dir
            os.makedirs(ARGS.output_dir)
        log_path = ARGS.output_dir

    # otherwise use any input_dir
    elif not ARGS.input_dir is None:
        if not os.path.isdir(ARGS.input_dir):
            print 'ERROR: the input_dir passed (' + ARGS.input_dir + ') is not valid'
        else:
            log_path = ARGS.input_dir

    # with neither input_dir nor output_dir passed, extract a path from fastq files passed
    elif ARGS.fastq_1 and ARGS.fastq_2:
        fqs = [ ARGS.fastq_1, ARGS.fastq_2 ]
        fq_path = set([ os.path.split(fq)[0] for fq in fqs ])
        if not len(fq_path) == 1:
            print 'ERROR: fastq files passed (' + str(fqs) + ') have different paths'
        else:
            (fq_path,) = fq_path
            if not os.path.isdir(fq_path):
                print ''.join(['fastq files passed have invalid path: ', fq_path ])
            else:
                log_path = fq_path

    assert log_path is not None, 'a valid path for logger could not be derived from the passed ARGS'

    full_log_path = ''.join([ log_path, '/', dir_name ])
    if not os.path.isdir(full_log_path):
        os.mkdir(full_log_path)
    return full_log_path
