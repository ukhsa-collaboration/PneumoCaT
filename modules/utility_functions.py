'''
TODO: Needs module definition, kinds of functions to be found here.
More description.
'''

import logging
import os
import sys, inspect

module_folder = os.path.dirname(os.path.abspath(__file__))
if module_folder not in sys.path:
    sys.path.insert(1, module_folder)

import phe_exceptions

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------


def try_and_except(error_filepath, function, *parameters, **named_parameters):
    """
    This wraps a function in try and except clause. If an error is caught this will trigger
    1) reporting of the error to stdout
    2) writing of the error into an error file
    3) a sys.exit with code 1

    Parameters
    ----------
    error_file : String
        path to log file to capture error in
        (a useful default = logger.handlers[0].baseFilename, which returns 
            the stderr.log FileHandler added first by log_writer.setup_logger)
    function: Function
        the function
    parameters : all non-named parameters for the function
    named_parameters : all named paramaters for the function

    Notes
    -----
    Returns the returns from the function called

    Examples
    --------

    assuming a function

    def my_func(a,b,c = None)
        ......
        return d

    This function would be called as follows:

    return_val = try_and_except("stderr.log", my_func, 1, 2, c = 3)
    """
    import traceback

    try:
        return function(*parameters, **named_parameters)
    except phe_exceptions.PheException as phe_e:
        # This exception is created by the 'call_external', when exit code != 0

        logger = logging.getLogger("stdout_stderr_logger")
        logger.exception(function.__name__ + " has raised an exception: \n" + str(phe_e))

        # Exit with the returncode specified in the called process.
        # TODO: Should it be a different return code? E.g. ranged for traceback.
        sys.exit(phe_e.phe_return_code)
    except Exception:
        error_string = "There was an error in the function '" + function.__name__ + "'"
        error_divider = "_" * 60
        print error_string
        print error_divider
        traceback.print_exc()
        print error_divider

        error_file = open(error_filepath, "a")
        error_file.write(error_string + "\n")
        error_file.write(error_divider + "\n")
        traceback.print_exc(file=error_file)
        error_file.write(error_divider + "\n")
        error_file.close()
        sys.exit(1)


# -----------------------------------------------------------------------------


def check_file_exists(filepath, file_description):
    """
    A function to check if a file exists.
    It will print out an error message and exit if the file is not found

    Parameters
    ----------
    filepath : String
        the path to the file to be checked
    file_description : String
        a description of the file to be checked e.g "config file"
    """
    if not os.path.exists(filepath):
        print("The " + file_description + " (" + filepath + ") does not exist")
        sys.exit(1)


# -----------------------------------------------------------------------------


def write_component_complete(output_dir):
    '''
    Creates marker of complete module run by writing empty file
    'module_complete.txt' to result/module/. Intended to be
    called as final call from module.

    Parameters:
    -----------
        output_dir, str : full path to the component output dir

    Returns:
    --------
        none

    Side-effects:
    -------------
        writes empty file 'ComponentComplete.txt' into output_dir
    '''
    cc = '/'.join([ output_dir, 'ComponentComplete.txt' ])
    with open(cc, 'w') as cc_out:
        cc_out.write('')


# -----------------------------------------------------------------------------


