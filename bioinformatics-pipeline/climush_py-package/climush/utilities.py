import subprocess, re, sys, pathlib, shutil, json, tomlkit, gzip, psutil, math, zipfile
from pathlib import Path
from datetime import datetime
from functools import wraps
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from climush.constants import *

# tried to order by category but some required specific order due to dependencies on other functions

#######################
# MISC UTILITIES ######
#######################

# check if the provided file path is a pathlib Path object
def is_pathclass(file_path, exit_if_false=True):
    '''
    Checks if the input is a Path class.

    :param file_path: path in question
    :param exit_if_false: whether to exit the script if the
    file path is found to not be a valid Path object; defaults
    to true, as it tends to mess up downstream processes.
    :return: True/error message
    '''
    if isinstance(file_path, pathlib.PurePath):
        return True
    else:
        if exit_if_false:
            msg = f'The provided path, {file_path}, is not a valid Path object.\n'
            exit_process(message=msg)
            return None
        else:
            return False

# alert if multiple files matching the search pattern in the file path are returned
# prompt for user input if multiple or no files matching the pattern are returned
# moved here because used in continue_to_next()
def flag_multiple_files(file_path, search_for, auto_respond=False):
    assert is_pathclass(file_path, exit_if_false=False)

    result = list(file_path.glob(search_for))

    if len(result) == 1:
        return result[0]
    elif len(result) == 0:
        all_files = file_path.glob('*')
        print(f'No files matching the pattern \'{search_for}\' were detected in the file path: {file_path}. Do '
              f'you mean any of these files in this directory?')
        which_file = prompt_print_options([all_files, 'none of these (exit)'], auto_respond=auto_respond)
        if which_file in all_files:
            return which_file
        else:
            exit_process(message=f'The response \'none of these\' was chosen when searching for the correct'
                                 f'file matching the pattern: {search_for}')
            return None
    else:
        print(f'{len(result)} files matching the pattern \'{search_for}\' were detected in the file path: '
              f'{file_path}. Please type the number corresponding to the correct file to use:')
        which_file = prompt_print_options([result, 'none of these (exit)'], auto_respond=auto_respond)
        if which_file == 'none of these (exit)':
            exit_process(message=f'The response \'none of these\' was chosen when searching for the correct'
                                 f'file matching the pattern: {search_for}.')
            return None
        else:
            return which_file

def check_dir_exists(dir_path, auto_respond = False):
    '''
    Creates a Path object and checks if the directory exists.

    Takes the input path and checks whether the path is a Path object. If it is not,
    it converts it to a Path object. Then checks whether the input path is an existing
    directory. If it is not, it will prompt the user for a different file path. This
    function is recursive, in that it will rerun the function if a new file path is
    provided via the command line prompt.
    :param dir_path: path to the directory to check; can be a Path object or a string
    :param auto_respond: if set to True, and the provided path does not exist, then the
    function will cause the script it is used within to exit if the dir_path is not an
    existing file path
    :return: a Path object of an existing file path
    '''

    # create a path object from the input, if not already a Path object
    if isinstance(dir_path, pathlib.PurePath):
        if dir_path.is_absolute():  # check if absolute; needs to be to get full name of parent, etc.
            pass
        else:
            dir_path = dir_path.resolve()
    else:
        dir_path = Path(dir_path).resolve()  # also check if absolute here

    # check if the path exists, which is required for this script
    if dir_path.is_dir():  # if the input path is an existing directory...
        return dir_path  # return the Path as is

    else:  # if the input path is NOT an existing directory...

        # trigger yes/no/quit prompt for adding a different path
        msg = f'The provided input path is not an existing directory:\n '\
              f'\t{dir_path}\n'\
              f'Would you like to provide a different directory?'  # do not include last line break, prompt_ will add

        if auto_respond:
            print(msg)
            print(f'\tauto response: quit\n')
            sys.exit(1)
        else:
            prompt_yes_no_quit(message=msg)  # will exit here if no/quit is selected, continue to next line if 'yes'

            # if yes, then prompt for updated path
            print(f'Please provide the new file path to use below:\n')
            new_path = input('>  ')

            # recursively run function to ensure that this new path exists
            return check_dir_exists(new_path)

# exit current script due to error, save script and timestamp of where error occurred
def exit_process(message, config_section='error.message'):
    script_name = sys.argv[0]  # unsure if will get name of script it is executed in or the one it is compiled in
    exit_time = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    # write_to_config(config_section=config_section, script=script_name, timestamp=exit_time, details=message)
    print(message)
    print(f'Exiting {script_name}...\n')
    return sys.exit(1)

# recursive function that will seek out yes/no/quit response continuously until achieved
def prompt_yes_no_quit(message, auto_respond=False):
    '''
    User prompt accepting specific responses.

    Takes the input of a print message that shows just prior to a user prompt. The message
    should specify the potential responses of yes/y, no/n, or quit/exit (case insensitive).
    If the user inputs a response that is not one of these potential responses, it will provide
    a notification that the input response was not recognize, and relist the acceptable responses.
    It then will immediately reprint the initial prompt to allow the user another chance to enter
    and acceptable input string.
    :param message: any message describing what is wanted from the user; see the return for how
    to set up the message so that the desired output is achieved
    :return: if response is yes/y, it will return nothing and then next line of code in the
    script will be run; if response is no/n/quit/exit, then it will exit the current script; if the
    response is unfamiliar, it will print an alert for invalid input, and restart the prompt
    again until an acceptable answer is received.
    '''

    print(f'{message} [yes/no/quit]\n')

    if auto_respond:
        return print(f'\tauto response: yes\n')

    response = input()

    if re.search(YES_RE, response, re.I):
        return None
    elif re.search(NO_RE, response, re.I) or re.search(QUIT_RE, response, re.I):
        print(f'Exiting...\n')
        return sys.exit()
    else:
        print(f'Your response was not recognized out of the list of potential responses. Please respond with '
              f'\'yes\', \'no\', or \'quit\'.\n')
        return prompt_yes_no_quit(message)

# print out list of options for option-based prompts
def prompt_print_options(option_list, auto_respond=False):

    option_list.sort()  # sort option list for easier user experience
    option_list.append('quit')  # add option to quit, listed last

    if auto_respond:
        print(f'options: {option_list}\n')
        print(f'\tauto response: quit\n')
        return sys.exit()

    for o,opt in enumerate(option_list):
        if (o+1) < len(option_list):
            print(f'\t[{o+1}]: {opt}')  # print out all options from option list, excluding final option
        else:
            selected_opt = input(f'\t[{o+1}]: {opt}\n')  # print out last option from option list as prompt

    result = option_list[int(selected_opt)-1]  # return the item from option list matching the user input

    if result == 'quit':
        return sys.exit()
    else:
        return result

# print a CLI prompt when multiple files are detected when only one is expected
def prompt_multiple_files(file_path, auto_respond=False):
    if isinstance(file_path, list):
        print(f'\nWARNING: Multiple files were detected when only '
              f'one was expected. Please type the number '
              f'corresponding to the correct file to use:')
        return prompt_print_options(file_path, auto_respond=auto_respond)
    else:
        if is_pathclass(file_path, exit_if_false=False):
            if file_path.is_dir():
                print(f'\nWARNING: Multiple files were detected in the '
                      f'{file_path.stem} folder. Please type the number '
                      f'corresponding to the correct file to use:')
                file_list = [file.stem for file in file_path.glob('*') if not re.search(HIDDEN_FILE_REGEX, file.stem)]
                return prompt_print_options(file_list, auto_respond=auto_respond)

# provide options for the sequencing platform if it cannot be otherwise detected
def prompt_sequencing_platform(sample_id, auto_respond=False):
    print(f'\nWARNING: The sequencing platform could not be inferred from the sample: {sample_id}. Please type '
          f'the number corresponding to the correct sequencing platform from the options below. ')
    platform = prompt_print_options([SEQ_PLATFORM_OPTS, 'multiple'],
                                    auto_respond=auto_respond)
    if platform == 'multiple':
        print(f'If you have a combination of sequences from multiple platforms, you will need to either:\n'
              f'(1) manually sort the files into their sequencing platform directories\n'
              f'(2) use the file renaming script, rename_sequence_files.py, to rename the files to match '
              f'the file naming convention.\n')
        return

# prompt user for generic input, with ability to assign an automated response depending on config settings
def prompt_generic(message, auto_respond=False, **kwargs):
    '''
    For use with generic user prompts, with option to customize the auto response using
    a keyword argument.
    :param message: prompt message to print through command line
    :param auto_respond: will check the configuration file to determine whether to use the
    automatic response (i.e., if settings['automate']['auto_respond'] == True)
    :param kwargs: currently only accepts a keyword pair for auto_response, where key needs to
    be auto_response and its value is the automatic response to use for this prompt
    :return: prompt response value
    '''

    # check for any keyword arguments
    try:
        auto_response = kwargs['auto_response']
    except KeyError:
        auto_response = None

    # print the response message regardless of possible error below, so that error is easier to trace if
    #  autoresponse is missing but required
    print(message)

    # if auto_respond is set to True, require an auto_response and exit if not provided
    if auto_respond and (auto_response is None):
        print(f'ERROR. If the pipeline is configured to run automatically, an auto-response is required '
              f'in order to continue. Please include an automated response for this prompt, and then rerun.\n'
              f'Exiting...\n')
        return sys.exit()

    # if auto_respond is set to True and an auto_response is provided, return the auto_response
    elif auto_respond and (auto_response is not None):
        return auto_response

    # if auto_respond is not set to True, then prompt user for input via command line, and return value
    else:
        user_response = input()
        return user_response

# print indented list
def print_indented_list(print_list):

    # CONFIRM INPUT IS A LIST

    # if input is already a list, continue to formatting
    if isinstance(print_list, list):
        pass

    # if not a list, but can be converted to a list, then converty type to list
    elif isinstance(print_list, set) or isinstance(print_list, tuple) or isinstance(print_list, np.ndarray):
        print_list = list(print_list)
    elif isinstance(print_list, pd.Series):
        print_list = print_list.to_list()

    # if not a list, and cannot obviously be converted to a list, then print unformatted
    else:
        return print(print_list)

    # add leading tab for first item to be printed in the list; does not lead with a line break
    print_list[0] = '\t' + print_list[0]

    # FORMAT INPUT LIST FOR PRINTING

    # add a line break and tab to the rest of the items in the
    # list (join will add after each existing list element)
    formatted_list = '\n\t'.join(print_list)

    # PRINT LIST FORMATTED WITH INDENTS AND LINE BREAKS

    # do not return an object, just print the formatted list (\n\t will be last thing to print)
    return print(formatted_list)

# run shell command and save stdout and stderr to file
def run_subprocess(cli_command_list, dest_dir, run_name, program=None, separate_sample_output=True, auto_respond=False):
    '''
    Run a subprocess, saving the output and error to a log file.

    Takes the command line (CL) argument assembled as a list, and runs
    the subprocess, saving the standard in (stdin), standard out (stdout),
    and standard error (stderr) from the process. The output of stderr and
    stdout are written to a log file in the destination directory
    :param cli_command_list: the list of the components of the command line
    argument to run through subprocess; each argument should be separated
    out into a list, rather than a single string.
    :param dest_dir: directory to write the file to; should be a Path
    object
    :param run_name: the name of the bioinformatics run describing the group of sequences
    being processed by this function; typically read in from the pipeline settings .toml file
    :param separate_sample_output: True/False; whether to include a separator in the output file
    between the standard output for each sample; in most cases, a separator is helpful
    when reviewing the output and seeing a clear distinction between samples, but for some
    programs like cutadapt, it prevents summary calculations from being performed
    :param auto_respond: True/False; whether to automatically respond to any prompts that
    might arise from this function
    :return: None
    '''
    assert is_pathclass(dest_dir)

    assert isinstance(cli_command_list, list) and len(cli_command_list) > 1, print(f'Command line process to run '
                                                                             f'should be in a list of its components. '
                                                                             f'Try shlex.split() if unsure how to '
                                                                             f'compose the list. ')

    # run_cmd = subprocess.run(cli_command_list, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    # out, err = run_cmd.communicate()
    if program is None:
        program = cli_command_list[0]

        if re.search(r'.\..', program):
            program = program.split('.')[0]

    # add the run name to the output file name, with the program name
    output_filename = '_'.join([program, run_name])

    run_cmd = subprocess.run(cli_command_list, capture_output=True)

    out_path = dest_dir / f'{output_filename}.out'
    with open(out_path, 'at') as fout:
        out_as_str = run_cmd.stdout.decode('utf-8')
        if separate_sample_output:
            if not int(out_path.stat().st_size) == 0:
                fout.write(f'-----------------------------\n')
            if len(out_as_str) > 0:
                fout.write(f'{out_as_str}\n')
        else:
            fout.write(f'{out_as_str}\n')

    err_path = dest_dir / f'{output_filename}.err'
    with open(err_path, 'at') as fout:
        err_as_str = run_cmd.stderr.decode('utf-8')
        if separate_sample_output:
            if not int(err_path.stat().st_size) == 0:
                fout.write(f'-----------------------------\n')
            if len(err_as_str) > 0:
                fout.write(f'{err_as_str}\n')
        else:
            fout.write(f'{err_as_str}\n')

    temp_file = dest_dir / f'{output_filename}.temp'

    if len(run_cmd.stderr) == 0:
        pass
    else:
        if not temp_file.is_file():  # if the temp_file does not exist (i.e., prompt already responded to)
            print(f'Running {program} produced an error. Please review the output in {err_path.name}. '
                  f'Would like to continue despite this error? [yes/no]')
            if auto_respond:
                print(f'\tauto response: yes\n')
                with open(temp_file, 'wt') as fout:  # still want to write to this file so it stops asking
                    fout.write('stop_prompt')
            else:
                continue_ok = input()
                if re.search(continue_ok, AFFIRM_REGEX, re.I):
                    with open(temp_file, 'wt') as fout:
                        fout.write('stop_prompt')
                else:
                    exit_process(message=f'Running {program} produced an error. See {err_path.name} for details.')
                    return None

    return None

# add a series of commands to a subprocess command list
def append_subprocess(cli_command_list, options_to_add, position, return_copy=False):

    ## CONFIRM IF LIST

    # check type of input received for options_to_add; if valid, create list
    if isinstance(options_to_add, str):
        options_to_add = [options_to_add]
    elif isinstance(options_to_add, list):
        pass
    else:
        print(f'ERROR. The input for command line options to add using {append_subprocess.__name__} must be either '
              f'a list or a string. The input provided is {type(options_to_add)}.\n')
        return sys.exit(1)

    ## MAKE POSITION POSITIVE INDEX

    # if a negative index is provided, convert to positive; makes everything below easier
    if position < 0:
        position = len(cli_command_list) + (position+1)

    ## CONFIRM POSITION IN BOUNDS

    # confirm that the list index provided in position is valid
    if isinstance(position, int):  # must be integer

        # get length of input list to know bounds; save as var, used several times below
        len_input_list = len(cli_command_list)

        # confirm that the index position is in bounds
        if position <= len_input_list + 1:
            pass

        # otherwise not valid
        else:
            print(f'ERROR. The provided index to add new command line options to, {position}, is outside of the '
                  f'bounds of the original command line list, which has a length of {len_input_list}.\n')
            return sys.exit(1)
    else:
        print(f'ERROR. The input value for the position to add new command line options in '
              f'{append_subprocess.__name__} must be an integer, but a {type(position)} was provided.\n')
        return sys.exit(1)

    # if you don't want to overwrite the input list, then create a copy before continuing
    if return_copy:
        cli_command_list = cli_command_list.copy()
    else:
        pass

    ## INSERT NEW COMMANDS AT RESPECTIVE CONSECUTIVE POSITIONS

    # go through each value in the options_to_add list
    for i in range(len(options_to_add)):
        position += i
        cli_command_list.insert(position, options_to_add[i])

    # RETURN DEPENDING ON OPTION

    if return_copy:
        return cli_command_list
    else:
        return None

# t = [1,2,3,4,5,6]
# t_in = ['A', 'B']
# append_subprocess(t, t_in, position=-8)

def func_timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = datetime.now()
        func_output = func(*args, **kwargs)
        end_time = datetime.now()
        runtime = str(end_time - start_time).split('.')[0]  # round seconds down
        print(f"{func.__name__} was executed in {runtime}.\n")
        return func_output
    return wrapper()

def get_seq_platform(fastx_file, delim='_'):
    if isinstance(fastx_file, str):
        sample_id = fastx_file
    elif is_pathclass(fastx_file, exit_if_false=False):
        sample_id = fastx_file.name
    else:
        print(f'ERROR. Could not read the file name when trying to detect the sequencing platform '
              f'from the file name: {fastx_file}\n')
        return None

    # check the file name components for any of the expected platforms: pacbio, illumina, sanger
    platform = []
    for v in sample_id.split(delim):
        if re.search(PLATFORM_ANYWHERE_RE, v, re.I):
            platform.append(v)
        else:
            continue

    # confirm that only one platform was encountered in the file name
    if len(platform) == 1:
        return platform[0]
    elif len(platform) > 1:
        print(f'There were multiple sequencing platforms detected:\n'
              f'\tfile name: {sample_id}\n'
              f'\tplatforms detected: {platform}\n')
    else:
        print(f'No sequence platforms were detected:\n'
              f'\tfile name: {sample_id}\n')

    print(f'May not be able to continue without a clear sequencing platform.\n'
          f'Exiting...\n')
    return None

def import_mapping_df(df_path, ignore_tabs=True, filter=True, auto_respond=False):
    '''
    Import .csv, .txt, or .xlsx table as a dictionary.

    Read in a mapping file for demultiplexing. Can accomodate the file formats
    .xlsx, .csv, and .txt. Will output a dictionary, where the key is the name
    of the tab in the dataframe and the value is the dataframe in that tab. Will
    always return a dictionary, although only .xlsx files will have tabs. Formats
    the name of the keys in the output dictionary to be 'pool01' or 'pool02'; the
    number of the pool is inferred from the name of the tab (.xlsx) or the name
    of the file (.csv, .txt). Returned as dictionary so that the output can be
    handled in the same way, regardless of file type (e.g., loop through tabs
    even if a .csv, which will have a single tab).
    :param df_path: path to the dataframe
    :param ignore_tabs: True/False; ignore any tabs in the mapping file that do not
    match 1 or 2 when looking for pool numbers; useful when mapping files have more
    than two pools (PacBio Revio) but you only want to use pool 1 and pool 2
    :param filter: True/False; whether to remove all tabs that are not detected as
    a pool tab (pool01/pool02); set to the default True because I am not sure that
    the demultiplex() function will work if more than just pool01 and pool02
    :param auto_respond: True/False; whether to automatically respond to the prompt
    within or not; this value should be provided from the settings read in from the
    configuration file, within ['automate']['auto_respond']
    :return: dictionary, where the key is 'pool1' or 'pool2', and the value is the
    dataframe belonging to that tab, or in the case of a .csv or .txt file, the
    entirety of that dataframe file.
    '''

    # if the mapping file is an excel file...
    if re.search(r'^\.x', df_path.suffix):

        # read in all tabs from the table
        mapping_tabs = pd.read_excel(df_path, sheet_name=None)  # need to set sheet_name to None to get all tabs
        old_tab_names = list(mapping_tabs.keys())  # make list of old names, otherwise cannot update key names in loop

        # create a set to add pool numbers to when they are found from tabs in .xlsx file
        # if there are multiple tabs, some of which are unrelated to pools, will need to reference this
        found_pools = {'pool1':'',
                       'pool2':''}
        miscell_tabs = []

        for tab in old_tab_names:

            try:  # format the name of the tab to be uniform across all tabs/dataframes
                pool_num = re.search(POOL_NUM_RE, tab, re.I).group(0).strip()

                # add the found pool number to the found_pools set
                if re.search(r'1', pool_num):
                    found_pools['pool1'] += tab
                elif re.search(r'2', pool_num):
                    found_pools['pool2'] += tab
                else:
                    if ignore_tabs:
                        continue
                    else:
                        print(f'ERROR. A pool number was detected in the tab {tab} in the mapping file {df_path.name}, '
                              f'but it could not be determined whether this was the number for pool 01 or pool 02.\n')
                        return sys.exit()

            except AttributeError:  # if this tab does not contain pool information, or at least not detected as such...

                # add this tab name to the miscell_tabs list to check at end whether might be a pool tab or not
                miscell_tabs.append(tab)

        for pool_num, pool_tab in found_pools.items():

            # if no tab was found matching this pool...
            if pool_tab == '':

                # print a warning and prompt user to chose the tab that corresponds to this pool
                print(f'WARNING. The pool number for {pool_num} could not be inferred from tabs in the '
                      f'mapping file {df_path.name}. Please enter the number of the tab that corresponds to '
                      f'{pool_num}')

                # will return the name of the tab to use, so update var pool_tab to match
                pool_tab = prompt_print_options(miscell_tabs, auto_respond=auto_respond)

            mapping_tabs[pool_num] = mapping_tabs.pop(pool_tab)  # replace old tab (key) with reformatted one


    # if the mapping file is a .csv or .txt file...
    elif re.search(r'^\.c|^\.txt$', df_path.suffix):  # I think you can read in .txt and .csv files the same way?

        try:  # try to get the pool number from the file name

            pool_num = re.search(POOL_NUM_RE, df_path.name, re.I).group(0)

        except AttributeError:  # if there's no detected pool number in the name, prompt user to specify one

            print(f'The pool number could not be inferred from the file name of the mapping file {df_path.name}. '
                  f'Please type the correct pool number for this file: ')
            pool_num = prompt_print_options(['1', '2'], auto_respond=auto_respond)  # choose from 1 or 2 (or quit, built into function)

        # make sure that the pool number tab is formatted correctly (single number, no leading zero, etc.)
        if re.search(r'1', pool_num):
            formatted_pool_num = '1'
        elif re.search(r'2', pool_num):
            formatted_pool_num = '2'
        else:
            if ignore_tabs:
                pass
            else:
                print(f'ERROR. A pool number was detected in the tab {tab} in the mapping file {df_path.name}, '
                      f'but it could not be determined whether this was the number for pool 01 or pool 02.\n')
                return sys.exit()

        mapping_tabs = {f'pool{formatted_pool_num}': pd.read_csv(df_path)}  # make dict to match format from .xlsx

    # if the mapping file doesn't appear to be .xlsx, .csv, or .txt...
    else:

        # can't use this file with the current version of the pipeline functions
        print(f'ERROR. The file format {df_path.suffix} of the mapping file {df_path.name} is not a recognized '
              f'file type. Accepted file types are: \'.xlsx\', \'.csv\', and \'.txt\'.\n')
        sys.exit()

    # filter out any other tabs in the input file if they are not the pool01 or pool02 tab
    if filter:
        filtered_tabs = {}

        # go through the tabs in the input dataframe (which now has updated pool1 and pool2 tabs)...
        for tab_name,content in mapping_tabs.items():

            # if you want to ignore any tabs that are not pool1 or pool2...
            if ignore_tabs:

                # only add the pool1 and pool2 tabs to the output dataframe
                 if tab_name in found_pools.keys():
                     filtered_tabs.update({tab_name: content})

            # if you want to include any tabs that weren't found to be 'pool' tabs...
            ## THIS PART DOESN'T MAKE SENSE ##
            else:

                if tab_name in miscell_tabs:
                    filtered_tabs.update({tab_name: content})
                else:
                    filtered_tabs.update({tab_name: content})

        return filtered_tabs

    # otherwise, return all tabs from input Excel table
    else:
        return mapping_tabs

# determine and run the next step in the pipeline
def continue_to_next(current_script, config_dict):
    """
    Continue to the next step of the pipeline.

    Determines which script should run next in the pipeline based on
    the settings in the configuration file. If the settings indicate
    that nothing should be run automatically, the function causes the
    input script to silently ext.
    :param current_script: the file path of the currently running script;
    input is typically __file__, which will pull the absolute path of
    the file (script) in which it is run
    :return: None, assembles and runs a shell command.
    """

    # if the input is not a Path class, make it a Path class
    if is_pathclass(current_script, exit_if_false=False):
        pass
    else:
        current_script = Path(current_script)

    # read in automation details from the configuration file
    automate = config_dict['automate']['run_all']  # whether to automatically run all steps
    steps_to_run = config_dict['automate']['run_some']  # which steps should be run automatically

    # set a variable so that a single print message can be used at the end for multiple situations
    no_action_needed = False  # set to False; if no next step to run, will switch to True

    # if automation is turned on, figure out which script should run next...
    if automate:

        ## GET NUMBER OF CURRENT SCRIPT

        # get the number of the current script that just completed
        try:
            current_script_number = int(re.search(PIPE_SCRIPT_NUM_RE, current_script.name).group(0))

        # if it cannot be located in the current script's file name, print ERROR and exit
        except AttributeError:
            msg = f'ERROR. The pipeline script number cannot be located in the following file name:\n'\
                  f'  {current_script.name}\n'\
                  f'The pipeline script number should contain at least two digits at the start of the '\
                  f'file name, indicating the order in which the script is run in the pipeline. If the '\
                  f'pipeline step is <10, it should be preceded by a 0, e.g., 03_remove-primers.py.\n'
            exit_process(msg)
            return None  # not really necessary here, but it is giving me warnings below that are annoying

        ## DETERMINE THE NUMBER OF THE NEXT SCRIPT TO RUN, BASED ON USER CONFIGURATION SETTINGS

        # if a list of scripts to run was provided, prioritize these scripts
        if len(steps_to_run) > 0:

            ## UPDATE THE PRIORITY STEPS TO RUN LIST, IF NECESSARY

            # determine whether the current script was included in this list, which it typically will be
            if current_script_number in steps_to_run:

                # if this current script was among these priority pipeline steps, do nothing to this list
                pass

            # if the script was not included in this priority list, print a warning and add it to this list
            else:

                # print a warning that the function will default to the next highest numbered script in the run list
                print(f'WARNING. The pipeline script that is currently running,\n'
                      f'  {current_script.name}\n'
                      f'is not listed in the specific list of pipeline scripts to run from the configuration '
                      f'file settings. Using the logical next script in the specified list of pipeline scripts '
                      f'to run, which is the next highest numbered script.\n'
                      f'\n'
                      f'This warning usually occurs if you run a specific pipeline script from the command line '
                      f'and do not update your configuration file list of pipeline scripts to run, run_some, to '
                      f'include this script.\n')

                # append the current script's number to the list of steps to run, and sort the list
                steps_to_run.append(current_script_number)  # add current script's number to list
                steps_to_run.sort()  # numerically sort the list, ascending


            ## FIND THE NEXT SCRIPT TO RUN BASED ON THE CURRENT SCRIPT'S POSTION IN THE PRIORITY RUN LIST

            # get the index (location) of th current script's number in the steps_to_run list
            current_script_runlist_i = steps_to_run.index(current_script_number)

            # if it is the last integer in this list (and therefore last script to run)...
            if current_script_runlist_i == (len(steps_to_run) - 1):

                # switch the no_action_needed Boolean to True, which will result in same output as no automation
                no_action_needed = True
                no_action_reason = ('No further pipeline steps are set to run, based on the run_some list from '
                                    'the configuration file')
                next_int = None

            # if the current script is not the last script that needs to run...
            else:

                # get the integer value of the next script to run
                next_int = steps_to_run[current_script_runlist_i + 1]

        # if a list of scripts to run was not provided, then run the pipeline scripts in numeric sequence
        else:
            next_int = current_script_number + 1


        ## FIND THE FILE NAME OF THE NEXT SCRIPT TO RUN

        # next_int set to None if there is not another script to run next
        if next_int is None:
            pass
        else:
            # format the next script integer as a string, w/ (0-9) or w/o (10+) leading zero
            if next_int >= 10:  # if in double-digits, don't add a leading zero
                next_num_str = str(next_int)
            else:  # if single digit, add leading zero
                next_num_str = '0' + str(next_int)

    # if set not to automatically run any other scripts...
    else:

        # set the no_action_needed variable to true, as no further scripts will be run
        no_action_needed = True
        no_action_reason = 'Pipeline is configured to run manually'


    ## DETERMINE FUNCTION RETURNED BASED ON ACQUIRED INFORMATION

    # if pipeline has completed or if it has been set to run only manually...
    if no_action_needed:

        # print message, in case confusion about a script terminating, and exit with success code
        print(f'{no_action_reason}. {current_script.name} has completed, exiting pipeline...\n')
        return sys.exit(0)

    # if another script needs to run...
    else:

        # find the file name/path of the next script, based on its number prefix determined above
        next_script = flag_multiple_files(file_path=current_script.parent, search_for=f'{next_num_str}*',
                                          auto_respond=config_dict['automate']['auto_respond'])

        print(f'Running next step, {next_script.name}...\n')
        return subprocess.run(['python3', next_script])

def file_finder(reference_dir, search_glob, multiple_matches=False, max_depth=5000):
    '''
    Locate a file using a glob string via a directory walk.

    :param reference_dir: a Path object to a directory from which the directory
    walk will begin
    :param search_glob: a glob-formatted search string that is used to find the
    target file
    :param multiple_matches: True/False; whether to return multiple matches as a list (True) or
    print all matching file paths and require a choice of a single file path (False)
    :param max_depth: the maximum file depth to search for a matching file
    :return: if successful, returns a Path object path to the matching file; if multiple
    matches are found, the user is prompted to select the correct target path; if no
    matches are found, the script will exit with a descriptive error message
    '''

    # if the input directory is not a Path object, create a Path object from it
    if is_pathclass(reference_dir, exit_if_false=False):
        pass
    else:
        reference_dir = Path(reference_dir)

    # if the reference_dir is instead a file, the path walk won't work
    # if a file path is provided, use its parent instead as the reference_dir
    if reference_dir.is_dir():
        pass
    else:
        reference_dir = reference_dir.parent

    # set a boolean variable to keep track of whether at least a single file was located
    file_not_found = True

    # create a list of file Paths that were found in the walk that match the search_glob provided
    files_matching_search = []

    # start with a top-down search of the file tree from the reference directory provided
    for root, dirs, files in reference_dir.walk(top_down=True):

        # check for match to the search_glob among the files
        file_matches = [file for file in root.glob(search_glob)]

        # if there are matches, append these to the file match list, so all matches are gathered
        if len(file_matches) > 0:
            file_not_found = False
            files_matching_search = [*files_matching_search, *file_matches]

        # if there are no matches, check child directories
        else:
            continue

    # if no file was located through a top-down search, try a bottom-up search
    if file_not_found:

        # keep track of the depth searched, do not exceed the max_depth param
        depth_counter = 0

        # walk the parent of the input reference directory from bottom-up
        for root, dirs, files in reference_dir.parent.walk(top_down=False):

            # only search up to a certain depth, as defined in params
            if depth_counter > max_depth:
                break
            else:
                depth_counter +=1

            # check for match to the search_glob among the files
            file_matches = [file for file in root.glob(search_glob)]

            # if there are matches, stop the walk and move on to file check
            if len(file_matches) > 0:
                file_not_found = False
                files_matching_search = [*files_matching_search, *file_matches]

            # if there are no matches, check parent directories
            else:
                continue

    # if still no matching file was located, print an error and exit with code 1
    if file_not_found:
        err_msg = (f'A file matching the glob {search_glob} was not located through a path walk starting at:\n'
                   f'   {reference_dir}\n'
                   f'when searched both top-down and bottom-up (with a maximum bottom-up depth of {max_depth}). '
                   f'Check that the file exists in relative proximity to the reference directory and retry.\n')
        return exit_process(err_msg)

    # check how many file matches there are, since there should only be 1

    # if one file match is found, return this single path
    if len(files_matching_search) == 1:
        return files_matching_search[0]

    # if multiple file matches are found, prompt user to select the correct match
    else:

        # check if multiple files is okay and should be returned as list
        if multiple_matches:
            return files_matching_search

        # if multiple matches not okay, then prompt user to chose one from list of all matches
        else:
            return prompt_multiple_files(files_matching_search)

def import_filepath(arg_value, must_exist, make_absolute=True):
    '''
    Import a file path from parsed command line options as a Path object.

    Takes the value returned from a parsed argparse argument, typically as
    a dictionary element. Should also accept an argparse Namespace class as
    long as it has already been indexed into a string (makes no difference).
    :param arg_value: value from the parsed command line argument that expects
    a path-like string or Path object (e.g., for default values from config files)
    :param
    :param make_absolute: whether to convert the Path object to an absolute
    file path if a relative file path is provided
    :return: the input path as a Path object
    '''

    # convert input to a Path class object
    if is_pathclass(arg_value, exit_if_false=False):
        input_path = arg_value
    else:
        input_path = Path(arg_value)

    # check through the rest of the settings for the path (must_exist, make_absolute)

    # if the path must exist and want to make it absolute...
    if make_absolute and must_exist:  # [T/T]
        try:  # using strict=True will both check that the path exists and make it absolute
            return input_path.resolve(strict=True)
        except FileNotFoundError:  # if it doesn't exist, then FileNotFoundError will be triggered
            # in this case, print an error message with the error-causing path, and exit function
            # because there is one other section with same error message, just pass and will print error + exit at end
            pass

    # if only want to make the path absolute...
    elif make_absolute:  # [T/F] if it were both, it would have already execute first if statement
        return input_path.resolve(strict=False)

    # if I only want to ensure that the path exists...
    elif must_exist:  # [F/T] if make_absolute were true, it would have executed by now
        if input_path.is_dir() or input_path.is_file():  # accept either directories or files
            return input_path
        else:
            # because there is one other section with same error message, just pass and will print error + exit at end
            pass

    # if I don't want either setting; I don't want to require that path exists nor that the path be absolute
    else:
        return input_path

    # print collected ERROR message, if any, and exit script if this point is reached
    # only errors will not have returned anything by now
    msg = (f'ERROR. The following path provided does not exist: \n'
           f' {input_path} \n'
           f'A path to an existing file or directory is required for this function. '
           f'Please check this path and then rerun.\n')
    print(msg)
    return sys.exit()

def filter_empty_files(dir_path, keep_empty=True):
    empty_files = 0
    empty_folders = 0
    for file in dir_path.glob('*'):
        if file.stat().st_size == 0:
            empty_files += 1
            if keep_empty:
                # REPLACE W/ MKDIR EXIST OKAY
                dest_dir = dir_path / 'empty-files'
                dest_dir.mkdir(exist_ok=True)
                dest_path = dest_dir / file.name
                file.rename(dest_path)
                fate = f'moved to the directory {dest_dir}'
            else:
                file.unlink()
                fate = 'deleted'
        elif file.is_dir():
            try:  # empty dirs not always size 0 (hidden files); this will only rm dir if it is empty (but not 0)
                file.rmdir()
                fate = 'deleted'
                empty_folders += 1
            except:
                fate = ''
                continue
        else:
            fate = ''
            continue

    return print(f'{empty_files} files and {empty_folders} folders were {fate}.\n')

# check whether a path exists; return path if it does or raise error if it does not
def flag_if_not_path_exists(file_path, absolute=True, exit_if_false=True):
    '''
    Create a Path class using the provided filepath string.

    If the filepath does not exist, send an error message,
    and will exit the script.
    :param filepath: string of absolute or relative path.
    :param absolute: default True; if True, will return the
    Path as an absolute path; otherwise will remain relative
    (which may cause issues when joining paths later on)
    :param exit_if_false: whether to exit the Python script if
    the provided path does not exist; defaults to True
    :return: filepath as PosixPath if it exists, otherwise
    an error will exit the script.
    '''
    try_path = Path(file_path)

    if try_path.exists():
        if absolute:
            return try_path.resolve()
        else:
            return try_path
    else:
        msg = f'Cannot locate the filepath {file_path}. Please '\
              f'make sure that either the original file structure '\
              f'of the container is unchanged, or use the --seq-path/-sp '\
              f'flags while running this script.\n'
        print(msg)
        if exit_if_false:
            exit_process(message=msg)
            return None
        else:
            return None

# create a log file to write details of pipeline step to
def make_log_file(file_name, dest_path, log_suffix=LOG_SUFFIX):
    '''
    Create file path for a log file.

    Files with multiple file extensions (e.g., fastq.gz) are
    not easily manipulated with most standard libraries like
    pathlib. Only '.gz' will be recognized as the file suffix/
    extension. This function deals with the different approaches
    required when the file extension is complex (.fastq.gz) and
    when it is simple (.fasta)
    :param file_name: name of the file to use to base
    the log file name off of
    :param log_suffix: the suffix/file extension to use for the log
    file; defaults to '.log'
    :param dest_path: directory to create the file
    :return: path to log file
    '''
    assert is_pathclass(file_name)
    assert is_pathclass(dest_path)

    if not log_suffix.startswith('.'):  # for user-provided suffix, check there is a '.'
        log_suffix = '.' + log_suffix

    if re.search(GZIP_REGEX, file_name.name):
        log_file_name = re.sub(GZIP_REGEX, log_suffix, file_name.name)
    else:
        log_file_name = file_name.with_suffix(log_suffix)

    return Path(dest_path / log_file_name)

# create directory if it doesn't exist and return path, or return the directory path if it already exists
def mkdir_exist_ok(new_dir, parent_dir=None):
    '''
    Make a new directory if one does not exist. Bypass making
    the directory if it does exist.
    :param new_dir: name of new directory as string
    :param parent_dir: path object to new directory; if none
    is specified, the new directory must be an absolute path.
    :return: Path object of directory path
    '''
    if not (isinstance(new_dir, str)) and (new_dir.is_absolute()):
        new_dir.mkdir(exist_ok=True, parents=True)
        return new_dir
    elif (parent_dir is None):
        return print(f'\nThe provided file path of the new directory is not absolute. Please '
                           f'provide a parent_dir to specify the path of the new directory, and '
                           f'retry.\n')
    else:
        new_path = Path(parent_dir / new_dir)
        new_path.mkdir(exist_ok=True, parents=True)
        return new_path

# move files from current location (source_path) to new location (dest_path)
def move_file(source_path, dest_path):
    assert isinstance(source_path, pathlib.PurePath), f'The provided source path is not a Path object: {source_path}'
    assert isinstance(dest_path, pathlib.PurePath), f'The provided destination path is not a Path object: {dest_path}'

    full_dest_path = dest_path / source_path.name

    source_path.rename(full_dest_path)

    if full_dest_path.is_file():
        return None
    else:
        return print(f'The file {source_path.name} was not properly copied to its new destination: {dest_path}.')

# count the number of files in a directory
def count_files(file_path, search_for='*'):
    '''
    Counts the number of files in a file path, and returns as int.

    :param file_path: path to directory to count files in.
    :param search_for: the .glob string that describes the file types
    to search for; default is '*', which is all files
    :return: int of the number of files in directory
    '''
    # check if it is a file path that exists
    flag_if_not_path_exists(file_path, absolute=False, exit_if_false=False)

    return len(list(file_path.glob(search_for)))

def check_for_input(file_dir, config_dict, path_must_exist=True, file_identifier=None, file_ext=SEQ_FILE_GLOB, file_prefix=None):
    '''
    Checks if there are input files for the process.

    :param file_dir: the directory path to check for files and to check
    whether the path exists
    :param config_dict: dictionary of imported settings from the pipeline
    configuration .toml file; required for inputing the auto respond param in
    check_dir_exists
    :param path_must_exist: True/False; whether or not the input directory must
    exist; default True means that the input path must exist, and if it does not
    exist, an exit code of 1 will result
    :param file_identifier: regex string or list of regex strings; a file identifier to look for
    within a file or directory name, that is as specific as the regex composed (i.e., if just a
    regular string and not a formatted regex, it will look for this string anywhere in the file
    name); typically the prefix of the sequence file that describes the sequencing platform
    (e.g., illumina) or sequenced gene region (e.g., its1, 18s); defaults to None, meaning this
    function will return True if any non-empty directory or file is located
    :param file_ext: the expected file extension of the files to search for,
    using asterisk '*' for wildcard (used in glob); if the desired return list
    is a list of subdirectories and not files, set this argument to None
    :param file_prefix: string; a file prefix to search for at the very start of a sequence file
    name; will only work if it is the start of the file name, which is different from file_identifier
    which can occur anywhere within the file name
    :return: two items; [1] Boolean T/F, whether there is a directory in this
    path that contains items matching the file_ext; [2] a list of the files matching the
    file extension in the directory, if present
    '''

    # if the path must exist to continue, then...
    if path_must_exist:

        # confirm that the path exists, which will return an absolute Path object,
        #   even if input is str and relative
        input_dir = check_dir_exists(file_dir, auto_respond=config_dict['automate']['auto_respond'])

        # once the path is confirmed to exist, check that it contains the desired content

        # if a file extension isn't provided, then the function looks for any subdirectories
        if file_ext is None:

            # get a list of subdirectories in the input directory
            input_subdirs = [d for d in input_dir.glob('*') if d.is_dir()]

            # if a specific prefix is required, filter out anything that doesn't have this prefix
            if file_prefix is None:
                pass
            else:
                input_subdirs = [d for d in input_subdirs if d.name.startswith(file_prefix)]

            # if there are subdirectories, return True and a list of the subdirectory paths
            if len(input_subdirs) > 0:
                return True, input_subdirs

        # if a file extension is provided, check that there are files within this directory
        #   that contain this type of file
        elif count_files(file_path=input_dir, search_for=file_ext) > 0:

            # then create a regex to match the sequencing platform files that you want returned in the file list
            if file_identifier is None:
                seq_re = r'.+?'  # will return all files in directory
            elif isinstance(file_identifier, list):
                seq_re = '|'.join(file_identifier)  # returns only those in provided in seq platform list
            else:
                seq_re = file_identifier  # returns only those of a single platform

            # use the platform-specific regex to create a list of relevant file paths
            file_list = [f for f in input_dir.glob(file_ext) if re.search(seq_re, f.name, re.I)]

            # if a file prefix is provided, further filter the file list for a files with specific prefix
            if file_prefix is None:
                pass
            else:
                file_list = [f for f in file_list if f.name.startswith(file_prefix)]

            # if at least 1 file was recovered matching the sequence platform(s)...
            if len(file_list) > 0:

                # return True for input files found, and the list of files matching the sequence platform(s)
                return True, file_list

            # if there are no files matching the wanted sequencing platform(s)...
            else:

                # return False for input files found, and the (empty) list of files
                return False, file_list

        # if the input directory is empty, pass on to the False and empty list return at the end
        else:
            pass

    # if the path does not need to exist to continue...
    else:

        # print a warning and continue to the False / empty list return
        print(f'The input directory\n'
              f'   {file_dir}\n'
              f'does not exist. Continuing...\n')

    # if directory is non-existent or empty, return False (no matching files/directories) and an empty file list
    file_list = []
    return False, file_list

# get name of previous script
def compress_data(input_path, output_path=None, compress_fmt='gzip', keep_input=False):
    '''
    Compress or decompress .zip or .gz files.

    :param input_path: path to the file or directory or directory of files that contain either
    compressed or uncompressed files that should be decompressed / compressed
    :param output_path: path to write the decompressed / compressed output files to
    this directory; if None, the files are written to the input_path
    :param compress_fmt: ['gzip','zip']; the compression algorithm of the input compressed
    files or the desired output compressed files / directories
    :param keep_input: True/False; if False (default), the input files are removed and only
    the output files are retained
    :return: a list of the output file paths
    '''

    ## CHECK FILE PATHS ##

    # confirm that the input path is a path object
    is_pathclass(input_path, exit_if_false=True)

    # if an output path was not provided...
    if output_path is None:

        # if the input path is to a file, use the parent directory as the output directory
        if input_path.is_file():
            output_path = input_path.parent

        # if the input path is a directory, use the input path as the output path
        else:
            output_path = input_path

    # if any input path was provided...
    else:

        # confirm that the output path is a path object
        is_pathclass(output_path, exit_if_false=True)

    ## CHECK COMPRESSION FORMAT ##

    # create a list of the compression formats currently working in this function
    accepted_compress_fmts = {
        'gzip': {
            'open_function': gzip.open,   # how to
            'magic_number': b'\x1f\x8b',
            'file_extension': '.gz',
        },
        'zip': {
                'open_function': zipfile.ZipFile,     # this doesn't work, need to create a ZipFile class object, then use .
                'magic_number': b'\x50\x4b\x03\x04',
                'file_extension': '.zip',
            },
    }

    # confirm that the argument for compress_fmt is a valid compression format
    if compress_fmt in accepted_compress_fmts:
        pass
    else:
        err_msg = (f'The compression format provided, {compress_fmt}, is not one of the accepted compression '
                   f'file formats for this function at this time. The accepted formats are:\n'
                   f'   {accepted_compress_fmts}\n')
        raise KeyError(err_msg)


    # create a list of the files in the input file path; a list will be created even if the input is a path to a file
    input_file_list = create_file_list(input_path)

    # create an empty list to add the input names to and make updates in cases where a name conflict w/ output file occurs
    input_files_corrected = []

    # keep track of the output file names, to return a list of the output file paths
    output_file_list = []

    # iterate through each file in the input directory (or if input is a file, just that file)
    for input_file in input_file_list:

        ## IS THE INPUT FILE COMPRESSED? ##

        # check if the input file is zip compressed
        if zipfile.is_zipfile(input_file):
            input_compressed = True
            input_fmt = accepted_compress_fmts['zip']['file_extension']
            output_fmt = ''

        # if not zip compressed...
        else:

            # open file and inspect the first few lines...
            with open(input_file, 'rb') as file_in:

                # if the first two lines match the gzip magic number, it is compressed
                if file_in.read(2) == accepted_compress_fmts['gzip']['magic_number']:  # compare the first two bytes
                    input_compressed = True
                    input_fmt = accepted_compress_fmts['gzip']['file_extension']  # .gz
                    output_fmt = input_file.suffixes[0]                           # whichever file ext precedes .gz

                # if the first two lines don't match the gzip magic number, it is not compressed
                else:

                    input_compressed = False

                    # check if the input file has a erroneous compressed file format (issue w/vsearch output)
                    if len(input_file.suffixes) > 1:

                        # create a new file name for the input file so it isn't named the same as the output file
                        input_file_new = input_path / input_file.name.replace(input_file.suffixes[1], '')  # remove compress fmt

                        # rename the input file using this re-formatted file name
                        input_file = input_file.replace(input_file_new)

                    # if the input file name only has an uncompressed file suffix...
                    else:
                        # do nothing here
                        pass

                    # add the input file, as it originally was or corrected, to the input_files_corrected list
                    input_files_corrected.append(input_file)

                    # create the output file format based on the single input file suffix and the output compressed file ext
                    output_fmt = input_file.suffix + accepted_compress_fmts[compress_fmt]['file_extension']


        ## CREATE OUTPUT FILE PATH ##

        # get the input file name without any file extensions included
        output_basename = input_file.name.replace(''.join(input_file.suffixes), '')

        # create a file path by including the output file extension and path
        output_file = (output_path / output_basename).with_suffix(output_fmt)

        # add the output file to the output file list
        output_file_list.append(output_file)


        ## DECOMPRESS INPUT FILE IF IT IS COMPRESSED ##

        # if input file is compressed...
        if input_compressed:

            # use the compression type's open() function to read in the information within it
            with accepted_compress_fmts[compress_fmt]['open_function'](input_file, 'rb') as compressed_in:

                # open the uncompressed output file to write the compressed content into...
                with open(output_file, 'wb') as decompressed_out:

                    # copy the content from the compressed input file to the decompressed output file
                    shutil.copyfileobj(compressed_in, decompressed_out)


        ## COMPRESS INPUT FILE IF IT IS DECOMPRESSED ##

        # if input file is uncompressed...
        else:

            # open the uncompressed input file to copy content into a compressed file...
            with open(input_file, 'rb') as decompressed_in:

                # use the compression type's open() function to open the file to write the uncompressed data into
                with accepted_compress_fmts[compress_fmt]['open_function'](output_file, 'wb') as compressed_out:

                    # copy the content from the uncompressed input file to the compressed output file
                    shutil.copyfileobj(decompressed_in, compressed_out)


        ## CONFIRM OUTPUT FILE WAS CREATED ##

        # confirm that the output file was created
        if output_file.exists():
            continue
        else:
            raise OSError(f'The output file was not properly written to the file system:\n'
                          f'   {output_file}\n')


    ## OPTIONAL; REMOVE INPUT FILES ##

    # if the input files should be kept still...
    if keep_input:

        # do nothing here
        pass

    # if the default option is used, where keep_input=False...
    else:

        # check if input_files_corrected has files paths included...
        if len(input_files_corrected) > 0:

            # go through each of the input files from the corrected input file list...
            for input_file in input_files_corrected:

                # delete each input file
                input_file.unlink()

        # if there weren't any instances of renamed / corrected input file paths...
        else:

            # iterate through the original input file paths...
            for input_file in input_file_list:

                # delete each input file
                input_file.unlink()

    return output_file_list

# add prefix to file name
def add_prefix(file_path, prefix, dest_dir, action='rename', f_delim='_', output_compressed=False, replace_prefix=True):
    '''
    Add prefix to file name.

    Adds a prefix to the file name that indicates which step of the
    pipeline the given file is produced from. Typically, will use one
    of the constants for the prefix, but any string can be provided.
    Prefix constants only exist for the files that persist to the
    next step of the pipeline. To create prefixes for the non-persistent
    files (e.g., reads/samples that are filtered out), use the flip_prefix
    function from the textmanip module.
    :param file_path: either a directory containing files to add a prefix to
    or a single file to add a prefix to
    :param prefix: prefix to add to the file
    :param f_delim: file name separator to use between the prefix and the rest
    of the file name; defaults to an underscore
    :param output_compressed: True/False; whether to retain the compressed file
    extension if the file extension of the input file (`file_path`) is also a
    compressed file extension (e.g., .gz)
    :param replace_prefix: True/False; if True, the input file name in the input file
    path (`file_path`) has a file name prefix that should be replaced with the input
    string `prefix`; if False, the input file name in the input file path (`file_path`)
    does not have a file name prefix, or if it does, the prefix should be retained and not
    replaced by the `prefix` string
    :return: Path object with new file name
    '''

    ## CHECK FUNCTION INPUT ##

    # input file path must be a Path object
    assert is_pathclass(file_path, exit_if_false=False)

    # the action must be an accepted value
    if action in ['mkdir', 'rename', 'copy', None]:
        pass
    else:
        err_msg = (f'The input provided to the add_prefix() parameter `action`, {action}, is not a valid input. '
                   f'Use one of the following valid input options for '
                   f'`action` = `mkdir`, `rename`, `copy`, or None.\n')
        return exit_process(err_msg)

    # subfunction; at this time, this is the entire function of add_prefix() really
    def prefix_single(file_name):

        ## GET OLD FILE NAME ##

        # if output is a directory...
        if action == 'mkdir':

            # do not include the file extension in the old name
            old_name = file_name.stem

        # if the output is NOT a directory (i.e., is a file)...
        else:

            # if the input file has more than one file extension, it is likely compressed
            if len(file_name.suffixes) > 1:

                # if it is okay that the output file also has a compressed file extension...
                if output_compressed:

                    # use both the original name and file extension(s) of the input file path
                    old_name = file_name.name

                # if the output file should NOT have a compressed file extension...
                else:

                    # drop the compressed part of the file extension, retaining only the
                    #    first file extension of the input file
                    old_name = Path(file_name.stem).with_suffix(file_name.suffixes[0]).name

            # if the input file only has one file extension, no need to evaluate whether output should be compressed
            else:

                # use both the original name and file extension(s) of the input file path
                old_name = file_name.name


        ## CREATE OUTPUT FILE NAME WITH NEW PREFIX ##

        # if the input prefix should be replaced...
        if replace_prefix:

            # determine where to replace the old prefix with the new prefix by using the sequence
            #   platform as a reference point

            # use regex to try to locate the platform in the input file name
            platform_present = re.search(PLATFORM_ANYWHERE_RE, old_name, re.I)

            # if a match to a sequence platform is located in the input file name...
            if platform_present:

                # get location of sequence platform in the input file name
                location = platform_present.span()[0]

                # add the file prefix just before the sequence platform to create the output file name
                # only keep the part of the input file name from the sequence platform through the end, removing old prefix
                new_name = prefix + f_delim + old_name[location:]

            # if a match to a sequence platform is not located, then just replace whatever the first label of the
            #  input file path is
            else:
                # this was originally written to work with denoise()
                old_name_no_prefix = f_delim.join(old_name.split(f_delim)[1:])
                new_name = f_delim.join([prefix, old_name_no_prefix])

        # if there is not a prefix to be replaced in the input file name...
        else:

            # add to the existing input file name without replacing anything
            new_name = prefix + f_delim + old_name


        ## CREATE FUNCTION OUTPUT ##

        # rename
        if action == 'rename':
            # replace the input file name with the output file name
            new_path = file_name.rename(dest_dir / new_name)

        # copy
        elif action == 'copy':
            # make a copy of the original input file, renaming the copy with the new file prefix
            new_path = dest_dir / new_name
            shutil.copy(file_path, new_path)

        # make directory (mkdir)
        elif action == 'mkdir':
            # create a directory using the updated directory name with the prefix added
            new_path = dest_dir / new_name
            mkdir_exist_ok(new_dir=new_path, parent_dir=dest_dir)

        # other (None)
        else:
            # create a file path; does not make any changes to or within the file system
            new_path = dest_dir / new_name

        return new_path

    return prefix_single(file_path)

    # if file_path.is_dir():
    #     new_path_list = []
    #     for file in file_path.glob('*'):
    #         new_path_list.append(prefix_single(file))
    #     return new_path_list
    # elif file_path.is_file():
    #     return prefix_single(file_path)
    # else:
    #     msg = (f'FAILURE. Unrecognized '
    #            f' file; cannot recognize as either a directory or file.\n')
    #     exit_process(msg)
    #     return None

    # filename_start = file_name.name.split('_')[0]
    # if re.search(ANY_PLATFORM_REGEX, filename_start, re.I):
    #     return dest_path / f'{prefix}_{file_name.name}'
    # elif re.search(ANY_PLATFORM_REGEX, file_name.name.split('_')[1], re.I):
    #     return dest_path / file_name.name.replace(filename_start, prefix)
    # else:
    #     msg = f'\nERROR: Cannot locate start of the original filename to add prefix to for '\
    #           f'file: {file_name}. Please make sure that the file name either starts with one '\
    #           f'of the prefixes expected given the naming convention (i.e., sanger, illumina, pacbio), '\
    #           f'or one of these expected prefixes is the second component of the filename.\n'
    #     print(msg)
    #     return exit_process(message=msg)
    #     return exit_process(message=msg)

def script_name_as_dir(script_name, parent, remove_num=False, suffix=None):
    '''
    Create a directory matching the name of the script from which it is executed.

    :param script_name: typically can just use __file__, but if I were to write that
    into the function, it would be the name of the module script for the climush
    Python package, not the file the function was used in
    :param suffix: suffix to add to end of directory name; may be useful
    to reflect the user-settings for that file
    :param parent: the parent directory under which to nest the directory; will
    default to the output directory of the pipeline, as that's what this function
    tends to be used for.
    :return: Path object to resulting directory
    '''
    assert is_pathclass(parent)

    if not parent.is_dir():
        parent = mkdir_exist_ok(new_dir=parent)

    script_num = Path(script_name).stem.split('_')[0]
    script_name_nonum = re.sub(r'^\d+_(?=\w)',r'', Path(script_name).stem)
    completed_name_nonum = OUTPUT_DIRS[script_name_nonum]
    completed_name = f'{script_num}_{completed_name_nonum}'

    if remove_num:
        completed_name = re.sub(r'^\d+_(?=\w)',r'', completed_name)

    if suffix:
        out_path = parent / f'{completed_name}_{suffix}'
    else:
        out_path = parent / completed_name

    out_path.mkdir(exist_ok=True)

    return out_path

def flip_prefix(prefix_const):
    # not sure why I would check if the prefix was a path? I don't think I ever use a path, just a constant string?
    # was throwing error when checking for PhiX in prefilter step, though was creating directory, error creating file
    # if is_pathclass(prefix_const, exit_if_false=False):
    #     prefix_const = prefix_const.stem
    # else:
    #     pass

    if re.search(r'^no-', prefix_const, re.I):  # if provided prefix starts with no, remove no
        return re.sub(r'^no-', r'', prefix_const, re.I)
    else:  # if it does not start with no, add no
        return 'no-' + prefix_const

def escape_path_special(file_path):
    if is_pathclass(file_path):
        file_path = str(file_path)

    if re.search(r'\(', file_path):
        file_path = re.sub(r'\(', r'\(', file_path)

    if re.search(r'\)', file_path):
        file_path = re.sub(r'\)', r'\)', file_path)

    return file_path

def get_sample_id(file_path, log_error=False, platform=None):
    '''
    Return the sample ID as a string from a file path.

    :param file_path: Path object, whose .name attribute contains
    the sample ID
    :param platform: illumina, pacbio, sanger, or a custom platform
    prefix used in the file name
    :return: the sample ID of the input file, as a string
    '''

    # if at start, the platform is not provided...
    if platform is None:
        # try to get it from the file name
        platform = get_seq_platform(file_path)
        # if it is still None...
        if platform is None:
            platform_re = SAMPLE_ID_RE  # search for sample ID for any of the seq platforms
        else:  # if it did find a platform, then use in regex
            platform_re = f'{platform.lower()}' + r'_.+?(?=\.)'
    else:  # if a platform was provided, then use this in the regex
        platform_re = f'{platform.lower()}' + r'.+?(?=\.)'

    # try to locate the sample ID in the file path
    try:

        # if the sample ID is successfully retrieved, return the sample ID
        sample_id = re.search(platform_re, file_path.name, re.I).group(0)
        return sample_id

    # if there is an error finding the sample ID in the file path...
    except AttributeError:

        # if logging the error file and not exiting script...
        if log_error:

            # use the parent directory's name to create the log file name
            log_file = (file_path.parent / f'{file_path.parent.name}_sample-id-errors').with_suffix('.log')

            # if the log file has already been created, then append the new error file name to the log file
            if log_file.is_file():
                with open(log_file, 'at') as log_out:
                    log_out.write(f'{file_path.name}\n')

            # if the file doesn't exist yet, write this file name out, along with a description header
            else:
                with open(log_file, 'wt') as log_out:
                    log_out.write(file_path.name)

            return None

        # if NOT logging the error but instead exiting the script...
        else:

            # assemble error message
            msg = (f'ERROR. The provided file name, {file_path.name}, does not follow the expected '
                   f'naming convention, so the sample ID cannot be inferred from the file name. '
                   f'This will affect the logging of summary data to output files.\n')

            # print error and exit
            exit_process(message=msg)
            return None  # won't print if exit_process used as return?

def get_read_orient(file_path):
    '''
    Return the read orientation in file name as string (R1/R2).

    :param file_path: Path object, whose .name attribute contains
    the read orientation (R1 or R2)
    :return: the sample ID of the input file, as a string
    '''

    try:
        read_orient = re.search(ORIENT_RE, file_path.name, re.I).group(0)
        return read_orient
    except AttributeError:
        msg = (f'ERROR. The provided file name, {file_path.name}, does not follow the expected '
               f'naming convention, so the read direction cannot be inferred from the file name. '
               f'This will affect the logging of summary data to output files.\n')
        exit_process(message=msg)
        return None

def copy_original_files(directory, copy_directory='original_files', compress=True):
    '''
    Creates a copy of original sequence files to zipped subdirectory.

    Copies all sequence files in the provided directory
    :param directory: the path to the directory containing the
    sequence files to be renamed.
    :param copy_directory: name or path of the subdirectory to
    copy the original sequence files to; default is a directory
    called 'original_files' within the provided directory
    :param compress: whether to compress the folder containing the
    copies of the original sequence files; default is True and will
    zip the copy directory.
    :return: print confirmation that all files were successfully
    copied to the copy directory
    '''

    # convert input to a Path class object
    if is_pathclass(directory, exit_if_false=False):
        dir_path = directory
    else:
        dir_path = Path(directory)

    # create the output directory (copy) in same directory as the input directory
    if is_pathclass(copy_directory, exit_if_false=False):
        copy_path = copy_directory
    else:
        copy_path = dir_path / copy_directory

    # create the copy directory

    # if the directory already exists, prompt user for input
    if copy_path.is_dir():

        # prompt for user input
        msg = f'A directory with the name {copy_path.name} already exists in: \n'\
              f'\t{dir_path}\n'\
              f'Do you want to delete this directory and continue? If you do not want to delete this directory, '\
              f'then the program will exit.'

        # if response is yes, will continue; if no/quit, will exit here
        prompt_yes_no_quit(message=msg)

        # remove the pre-existing directory that matches the copy directory name
        shutil.rmtree(copy_path)

    # make the new copy path directory
    mkdir_exist_ok(new_dir = copy_path)

    # create a list of files to copy, excluding any non-sequencing files
    files_to_copy = [ file for file in dir_path.glob(GZIP_GLOB) ]

    # copy each sequence file to the new copy directory
    for file in files_to_copy:
        shutil.copy2(file, copy_path)  # to preserve extra metadata, use copy2

    # before archiving copy, get count of sequence files from source and destination
    num_copied_files = len(list(copy_path.glob('*')))
    num_original_files = len(files_to_copy)

    # confirm all files were copied to the copy directory
    if num_copied_files == num_original_files:

        # if all files copied, then compress output
        if compress:
            # zip the copy directory
            shutil.make_archive(copy_path, 'zip', copy_path)

            # remove the unzipped copy directory
            shutil.rmtree(copy_path)

        return print(f'\nA copy of {num_copied_files} original sequence files was successfully made in:\n'
                     f'\t{copy_path}')

    else:
        num_missing = abs(num_copied_files - num_original_files)
        return print(f'ERROR: {num_missing}/{num_original_files} sequence files from the source directory:\n'
                     f'\t{dir_path}\n'
                     f'were not copied to the destination directory:\n'
                     f'\t{copy_path}')

def rename_read_header(input_dir, run_name, file_format='.fasta', unique_headers=False, no_copy=True, append_sample_str=False):

    ## CHECK INPUT FILES ARE FASTA FORMAT FILES

    # this function has only be tested with .fasta formatted files, so confirm the file extension of the input

    # get a list of all fasta files in the input directory
    fasta_files = list(input_dir.glob(f'*{file_format}*'))

    # if there are no fasta files, immediately print error and exit
    if len(fasta_files) == 0:
        err_msg = (f'ERROR. No files in the input directory:\n'
                   f'   {input_dir}\n'
                   f'match the {file_format} file format. This function, {rename_read_header.__name__}, only works '
                   f'on seuqences files with the {file_format} format.')
        exit_process(message=err_msg)

    # if there are at least some fasta files in the input file path...
    else:

        # check if *all* input files are fasta files
        if len(fasta_files) == len(list(input_dir.glob(SEQ_FILE_GLOB))):  # compare strictly fasta files to all seq files
            pass

        # if not all of the input files are fasta files, print a warning that only the reads in the fasta files in
        #   the input path will be renamed
        else:

            # determine which files are sequencing files (e.g., fastq, fastx, fasta), but NOT fasta files
            non_fasta_files = list(set(fasta_files).difference(set(input_dir.glob(SEQ_FILE_GLOB))))

            # print warning, and the files that are excluded, but continue, renaming only fasta files
            print(f'WARNING. {len(non_fasta_files)} files in the input path are sequencing files, but not in '
                  f'the {file_format} that this function, {rename_read_header.__name__}, requires as input. '
                  f'These files, listed below, will not be renamed:')
            print_indented_list(non_fasta_files)

    ## CREATE COPY OF ORIGINAL INPUT FILES

    # if default setting is used, then no copy of the original files is created
    if no_copy:
        pass

    # if no_copy is set to True, then a copy of all original input files is created, and compressed to .zip format
    else:
        # use the copy_original_files custom function to copy all files and compress them
        copy_original_files(directory=input_dir,
                            copy_directory=(input_dir.parent / f'{input_dir.name}_original-headers'))


    ## COLLECT SAMPLE IDS AND ORIGINAL READ HEADERS

    # create a conversion dictionary, in case there's an error in renaming the headers it can be reversed
    #   key = sample ID, value = key/value pair of old read header, new read header
    conversion_dict = {}

    # go through each fasta file in the input directory
    for file in fasta_files:

        # get the sample ID from the input file name; log file name if function doesn't work
        sample_id = get_sample_id(file, log_error=True, platform='illumina')

        # if the sample ID cannot be found in the file name, continue to next file (file name will be logged)
        if sample_id is None:
            continue

        # if the sample ID was found in the file name, add the sample ID to convert dict then continue to records
        else:

            # add to dictionary with empty dict as value, to add key/value pairs of old and new read headers
            conversion_dict.update({sample_id: {}})

        # create an empty list to add updated sequence records to for this sample
        new_sample_records = []

        # go through each read in this fasta file...
        for record in SeqIO.parse(file, file_format.replace('.','')):

            # get the original read header to add to the conversion dict
            old_read_header = record.description

            if unique_headers:

                # get the last part of the Illumina sequencer read ID from the old read header, to use in new read ID
                unique_read_num = re.search(ILLUMINA_SEQ_OG_RE, old_read_header, re.I).group(0)

                # add this unique read number to the end of the sample ID to use as new read header
                unique_read_id = sample_id + '_' + unique_read_num

                # replace the old read ID for the new read ID in the read header
                new_read_header = re.sub(ILLUMINA_READ_ID_OG, unique_read_id, record.description)

            else:

                # just use the sample ID as the read header, with the read's size included
                new_read_header = re.sub(ILLUMINA_READ_ID_OG, sample_id, record.description)


            # before adding the new read header to the conversion dict, check if sample= should be appended to front
            if append_sample_str:
                new_read_header = 'sample=' + new_read_header
            else:
                pass

            # add the old and new read headers to the conversion dictionary
            conversion_dict[sample_id].update({old_read_header: new_read_header})

            # I don't think you can overwrite info in a record, I think you have to create a new one
            updated_record = SeqRecord(record.seq,
                                       name='',
                                       id=new_read_header,
                                       description='')

            # add this updated record to a list of updated records for this sample
            new_sample_records.append(updated_record)



        ## CHECK NEW READ HEADERS ARE UNIQUE

        if unique_headers:

            # compare whether the length of all read IDs and the unique read IDs are same (i.e., all are unique)
            # before writing the new sequencing file with the updated read headers, check that the new read IDs are
            #  unique within the sample

            # get a list of all new read IDs
            read_headers_all = list(conversion_dict[sample_id].values())

            # create a set from this list of new read IDs, which will contain only unique read IDs
            read_headers_unique = set(read_headers_all)

            # if all new read headers are unique, continue without doing anything
            if len(read_headers_all) == len(read_headers_unique):
                pass  # continue to writing out file with new read headers if they are all unique

            # if the new read headers are not unique...
            else:

                # compose an error message
                err_msg = f'ERROR. The last digits of the original read ID provided by the Illumina sequencer are not '\
                          f'unique. If it is important that they are unique, consider using a larger portion of the '\
                          f'original Illumina sequencer read ID.\n'

                # write out the sample name, and the original read headers that caused duplication of new read headers
                duplicates_outpath = input_dir.parent / 'new-read-ids_duplicates.txt'

                # for each unique read header...
                for read_id in read_headers_unique:

                    # if it occurs more than once in the total list of read headers (i.e., is a duplicate, triplicate, etc.)
                    if read_headers_all.count(read_id) > 1:

                        # name this read ID as the non-unique new read ID
                        nonunique_new_header = read_id  # new read ID that is duplicated across reads
                        nonunique_count = read_headers_all.count(read_id)  # number of times new read ID is duplicated

                        # create a list of the original read IDs that match this duplicated read ID
                        nonunique_original_headers = [old_header for old_header, new_header in conversion_dict[sample_id].items() if new_header == nonunique_new_header]

                        # if the duplicates file already exists, define mode as append text, to add to this file
                        if duplicates_outpath.is_file():
                            mode = 'at'

                        # if the duplicates file does not yet exist, define mode as write text to create this file
                        else:
                            mode = 'wt'

                        with open(duplicates_outpath, mode) as fout:
                            fout.write(f'{nonunique_new_header}------------------\n')
                            fout.write(f'  non-unique header count = {nonunique_count}\n')
                            for id in nonunique_original_headers:
                                fout.write(f'  {id}\n')

                # write out the conversion_dict that has been created so far, to a .json file w/ tag _incomplete
                with open(input_dir.parent / f'{run_name}_read-header-conversion_incomplete.json', 'wt') as json_out:
                    json.dump(conversion_dict, json_out)

                # exit script, printing composed error message
                exit_process(message=err_msg)

        # if it isn't required that all read headers within a sample are unique, continue to writing new records w/o action
        else:
            pass

        # after going through each read in this sample, write out the new records, with updated read IDs, to a new
        #    sequence file for this sample, using the original sample file name so that it overwrites the old one
        SeqIO.write(new_sample_records, file, format=file_format.replace('.',''))

    # write out the conversion dictionary to a .json file so that there is a record of the conversion
    with open(input_dir.parent / f'{run_name}_read-header-conversion.json', 'wt') as json_out:
        json.dump(conversion_dict, json_out)

    return None

# define function to simplify import of configuration files
def get_settings(reference_dir, default_only=True):
    '''
    ***TEMPORARILY ONLY ACCEPTS THE DEFAULT CONFIG FILE.***

    Imports user and default settings from .toml to a dictionary.

    Looks for a default and a user configuration file in the config directory of
    the bioinformatics pipeline. Reads in the information from the user configuration
    file with priority. If a given field is empty in the user configuration, indicating
    no user preference was provided for that setting, then the value for that field will
    be searched for in the default configuration file.
    :param reference_dir: the directory from which to start searching for the configuration files;
    typically, when used within the climush pipeline scripts, this will be the path of the script
    in which the function is called (e.g., Path(__file__))
    :param default_only: True; inserted this *TEMPORARY* parameter so that it will
    only read in the default configuration file, rather than compiling the user
    and default configuration
    :return: dictionary of settings, with keys as configuration file settings and
    values as parameters chosen by the user (priority) or default settings.
    '''

    # construct glob for default configuration file
    default_settings_glob = '*' + DEFAULT_SETTINGS_HANDLE + CONFIG_FILETYPE


    if default_only:

        # locate the default configuration file by walking file paths
        config_path = file_finder(reference_dir, search_glob=default_settings_glob)

        # read in the default configuration file as a dictionary
        with open(config_path, 'rt') as config_in:
            compiled_config = tomlkit.parse(config_in.read())

        return compiled_config

    else:
        err_msg = f'ERROR. Compilation of pipeline settings from both the user and default configuration files is '\
                  f'still under construction. If you see message, ensure that default_only=True for the function '\
                  f'get_settings().\n'
        return exit_process(err_msg)

# sort the sequences provided in the 'sequences' input directory
def sort_input_files(filepath_dict, to_sort='main'):
    '''
    Sort the input sequencing files in their respective sequence
    platform subdirectories, or into subdirectories for required
    pre-processing.
    :param input_file_path: Path object to the directory containing
    the input sequences.
    :return: will print a count of how many files were sorted and
    into which folders, and then prompt to next steps if in
    interactive, or automatically carry out next steps if in
    automated.
    '''

    # UPDATE TO CHECK FOR ALL POSSIBLE DELIMITERS?
    def follows_naming_convention(filepath):
        if re.search(ANY_PLATFORM_REGEX, filepath.stem, re.I):
            return True
        else:
            return False

    def is_demultiplexed(filepath):
        if follows_naming_convention(filepath) or re.search(r'^[A-Z]+', filepath.stem):
            return True
        else:
            return False

    def is_pacbio(filepath):
        if re.search(r'^pacbio|^\d{4}', filepath.stem):
            return True
        else:
            return False

    def is_illumina(filepath):
        if re.search(r'^illumina|^[A-Z]+', filepath.stem):
            return True
        else:
            return False

    def is_sanger(filepath):
        if re.search(r'^sanger', filepath.stem):
            return True
        else:
            return False

    # I guess converting a generator to a list (like in print statement below) exhausts it? so then loop won't work
    files_to_sort = list(filepath_dict['sequences'][to_sort].glob('*.fast*'))

    print(f'Sorting {len(files_to_sort)} sequences files...\n')

    for file in files_to_sort:
        if is_pacbio(file):
            if is_demultiplexed(file):
                if follows_naming_convention(file):
                    move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['pacbio']))
                else:
                    move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['rename']))
            else:
                move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['demux']))
        elif is_illumina(file):
            if is_demultiplexed(file):
                if follows_naming_convention(file):
                    move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['illumina']))
                else:
                    move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['rename']))
            else:
                move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['demux']))
        elif is_sanger(file):
            if is_demultiplexed(file):
                if follows_naming_convention(file):
                    move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['sanger']))
                else:
                    move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['rename']))
            else:
                move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['demux']))
        else:
            move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=filepath_dict['sequences']['unclear']))

    return None

# locate relevant mapping file for demultiplexing
# moved down because it requires the import_config_as_dict function
# def find_mapping_file(path_to_mapping, path_to_plex, path_to_config, **kwargs):
#     '''
#     Locate the barcode mapping file(s) for demultiplexing.
#
#     Attempts to automatically detect the location of the mapping file
#     corresponding to the multiplexed files located in the sorted
#     directory of sequences to demultiplex. It inspects the name of the files
#     that were sorted into the needs_demux folder of the sequences folder, and
#     identifies the name of the sequencing run. It then will go into the mapping
#     file directory in the config folder, mapping-files, and looks for the mapping
#     file that shares that same sequencing run name. If multiple mapping files
#     are detected for a single sequencing run, the user will be prompted to select
#     the file that should be used. If no matching mapping file is detected, it
#     will throw an error and exit the pipeline.
#     :param path_to_mapping: defaults to expected location within structured
#     file directories; otherwise can specify the path as a Path object
#     :param path_to_plex: defaults to expected location within structured
#     file directories; otherwise can specify the path to the mapping file(s) as
#     a Path object
#     :param kwargs: specify the bioinfo-settings configuration file dictionary
#     if one is already loaded into the environment; if none is provided, the
#     function will open the configuration file to create a new instance of the
#     dictionary; adding a dictionary, if available, will likely be quicker so it
#     is recommended to do so if possible
#     :return: if found, returns mapping file(s) as a list
#     '''
#     assert is_pathclass(path_to_mapping)
#     assert is_pathclass(path_to_plex)
#
#     # get the multiplexed file's queue ID, which the UO sequencing core uses as the seq run identifier in the filename
#     multiplex_ids = [file.stem.split('.')[0] for file in path_to_plex.glob('*.fast*')]
#
#     # look up the multiplex IDs in the bioinfo-settings configuration file for its corresponding CliMush run name
#     kwargs_keys = list(kwargs.keys())
#     if kwargs_keys > 0:
#         if kwargs_keys > 1:
#             msg = f'ERROR: Too many keyword arguments (**kwargs) used. Only one kwarg is accepted and should have the '\
#                   f'name of the configuration dictionary as the value to a key with any name.\n'
#             print(msg)
#             return exit_process(message=msg)
#         else:
#             config_multiplex_dict = kwargs[kwargs_keys[0]]['pacbio-multiplex-ids']
#     else:
#         config_multiplex_dict = get_settings(file_handle=path_to_config)['pacbio-multiplex-ids']
#
#     # get the CliMush run name(s) and create regex to search for any of these at start of mapping file name
#     climush_ids = [config_multiplex_dict[q] for q in multiplex_ids]
#
#     # find the corresponding barcode mapping files for the run name(s) in the barcode-mapping dir in config
#     mapping_file_names = []
#     for id in climush_ids:
#         match = list(path_to_mapping.glob(f'*{id}*'))
#         if len(match) == 1:
#             mapping_file_names.append(match[0])
#         elif len(match) == 0:
#             msg = f'ERROR: No mapping files matching the CliMush run name {id} were located in {path_to_mapping}. ' \
#                   f'Check the mapping-files directory in the config folder, then retry.\n'
#             print(msg)
#             return exit_process(message=msg)
#         else:
#             msg = f'ERROR: {len(match)} mapping files were located for CliMush run name {id}. Please select the number ' \
#                   f'of the correctly corresponding barcode mapping file for this sample:\n'
#             print(msg)
#             mapping_file_names.append(prompt_print_options(match))
#
#     return mapping_file_names


def rescale_percent_identity(input_df, column='percent_match'):
    '''
    Standardize how percent identity is scaled.

    Percent identity scores from taxonomic assignments can be scaled from 0-100 or 0-1, depending on
    the algorithm used. This function will look for a column named 'percent_match' (default) in the
    input table, and scale any values between (and including) 0 and 1 and convert them to a percentage
    scale between 0 and 100. If the values are already between 1 and 100, they will not be altered.
    :param input_df: dataframe containing percent identity values to be scaled
    :param column: the name of the column in the input dataframe that contains the percent identity values.
    :return: a copy of the dataframe, with the scaled percent match values replacing the input
    percent mach values
    '''

    # make sure that the column dtype is float, tends to be imported as a string
    output_df = input_df.astype({column:'float'}, copy=True)

    # create a list of values where the percent_match values are all percents and not decimals
    percent_vals = []
    for val in output_df[column]:

        # if the percent match value is less than or equal to 1, it is likely a decimal form, so convert to percent
        if val <= 1.0:
            percent_vals.append(val*100)

        # if the percent match value is greater than 1 but less than 100, it should already be a percent
        #  I don't think you would get 1%, so it seems more likely to say a value of 1.0 is 100% decimal, not 1% percent
        elif 1.0 < val <= 100.0:
            percent_vals.append(val)

        else:
            raise KeyboardInterrupt(f'The value {val} is outside of the range of values that are evaluated in this '
                                    f'loop. Add another contingency before continuing.')

    # with this list of percent values, replace the mixed decimal/percent column with this list of percent only
    output_df['percent_match'] = percent_vals

    return output_df

def sort_taxonomy_info(input_df, drop_col=True):
    '''
    Split and format taxonomic information from one column to several.

    This function has only been testing on taxonomy information formatted by amptk taxonomy. When
    provided an input dataframe, the column name 'Taxonomy' is searched for. The string in this
    Taxonomy column is then split and manipulated so that 'clean' versions of each piece of information
    in the taxonomy string is sorted into its own column. This includes not only the taxonomic
    information (e.g., species, genus) but also the reference match information.
    :param input_df: a dataframe, typically an OTU table, that has a Taxonomy column to be formatted
    :param drop_col: True/False; whether to drop the Taxonomy column from the dataframe after extracting
    its information, as the column should now be redundant
    :return: a copy of the input dataframe with the split and formatted information from taxonomic
    identification
    '''

    # create a copy of the input dataframe, as to not change the original version
    dataframe = input_df.copy()

    # create a list of new columns to add to the input dataframe, separated by taxonomy info and reference match info

    # create a list of column names for splitting the taxonomy columns
    tax_cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    tax_codes = ['k', 'p', 'c', 'o', 'f', 'g', 's']

    # create a dictionary of the one-letter codes in the amptk taxonomy column with the full string version
    tax_key = {k: v for k, v in zip(tax_codes, tax_cols)}

    # new column names for the reference match info
    info_cols = ['amptk_method', 'percent_match', 'genbank_ref', 'unite_ref']

    # create a dictionary to store values for each of the new columns
    new_tax_cols = {col: [] for col in (info_cols + tax_cols)}

    # go through each ASV and sort out the taxonomy
    for tax_str in dataframe['Taxonomy']:

        # make a asv-level dict, so it is clearer which of the new columns needs an empty string before adding to megadict
        # fill all values with empty string; will either be replaced with a value, or left an empty str if info missing
        asv_tax_info = {col: '' for col in (info_cols + tax_cols)}

        ## SORT TAXONOMIC LEVELS ##

        # split string into the taxonomic levels
        tax_info = tax_str.split(';')[1].split(',')

        # keep track of which taxonomic levels had values that were handled; otherwise, will add empty strings
        tax_lvl_sorted = set()

        # some of the ASVs have 100% match and genbank/UNITE refs, but no taxonomy assigned?
        if tax_info == [''] or tax_info == ['No hit']:

            # do nothing here, tax_lvl_sorted will be empty, and new columns handled below
            pass

        # if at least some of the taxonomic info is there, go through each string
        else:

            # sort the taxonomy into the different taxonomic levels
            for t in tax_info:
                # get the taxonomic level string based on the one-letter code in the amptk output info
                tax_lvl_abbr = t.split(':')[0]
                tax_lvl = tax_key[tax_lvl_abbr]

                # add the taxonomy value to the new columns
                asv_tax_info.update({tax_lvl: t.split(':')[1]})

                # add the taxonomic lvl to the set that keeps track of which taxonomic levels have been dealt with
                tax_lvl_sorted.add(tax_lvl_abbr)

        ## SORT REFERENCE DB INFO ##

        # split the taxonomy string for a list of references related to the taxonomy match
        ref_match_info = tax_str.split(';')[0].split('|')

        # sort the reference db match info into the new columns

        # create a regular expression for three of four pieces of info in the ref db match info
        ref_regex = {r'^[A-Z]{2,3}$': 'amptk_method',
                     r'^([0-9]{1,3}\.[0-9]{1})|(0\.\d{4})$': 'percent_match',
                     r'^SH': 'unite_ref'}

        # genbank ref number is more variable, so anything not matching a regex in the ref_regex dict will be the genbank ref

        # go through each regular expression...
        for regex, new_col in ref_regex.items():

            # and try to match to an item in the ref match info
            match_found = [info for info in ref_match_info if re.search(regex, info)]

            # check if and how many matches were found for this regex/new column
            if len(match_found) == 1:

                # if a match is found, add it to the new column dictionary
                asv_tax_info.update({new_col: match_found[0]})

                # and remove it from the list, to keep track of which values have made it into the new columns
                ref_match_info.remove(match_found[0])

            # if multiple matches are found, then the regex is likely not specific enough and should be reworked
            elif len(match_found) > 1:
                raise KeyboardInterrupt('Refine the regex used to locate the column values, as the current version '
                                        'matches multiple values rather than just one.')

            # if no match found, then move to next regex
            else:
                continue

        # after using the regex to sort the values, confirm there is one more unsorted value, which should
        #  be the genbank ref number
        if len(ref_match_info) == 1:
            asv_tax_info.update({'genbank_ref': ref_match_info[0]})
        elif len(ref_match_info) == 0:
            raise KeyboardInterrupt('There are no remaining reference db info strings for the GenBank reference '
                                    'number; confirm that no reference number exists before continuing.')
        else:
            raise KeyboardInterrupt('Multiple values remain after sorting the amptk method, percent match, and unite '
                                    'reference number; check the values to see what is left unsorted.')

        # update the big dictionary with info for this asv
        for col, val in asv_tax_info.items():
            new_tax_cols[col].append(val)

    # add the taxonomy info from the dictionary into the input table, just after the ASV ID
    for c, col in enumerate(new_tax_cols):
        dataframe.insert(loc=(c + 1), column=col, value=new_tax_cols[col])

    # rescale the percent identity column, so all values are on the scale of 0-100, not a mix of 0-1/0-100
    dataframe = rescale_percent_identity(input_df=dataframe)

    # if drop_col=True, return the dataframe without the Taxonomy column
    if drop_col:
        return dataframe.drop('Taxonomy', axis=1)

    # if drop_col=False, return the dataframe with the Taxonomy column included
    else:
        return dataframe

def strip_file_ext(file_path, common_suffixes=KNOWN_FILE_SUFFIXES, verbose=False):
    '''
    Remove suffix or suffixes from a file name or file path.

    This function will remove one or more suffixes (i.e., file extensions)
    from a file name and return the file name without any suffix. This is an
    adaptation using the pathlib library so that file paths with periods in
    the file name and file paths with two suffixes can still be properly trimmed
    of their suffix(es).
    :param file_path: Path object to a file name that needs it suffix(es) removed
    :param common_suffixes: a list of common  file extentions that are
    likely to be encountered if the input file has more than one suffix; will help
    to determine if multiple recovered suffixes are true file extension components
    or part of a file name that includes periods.
    :param verbose: will print a warning if an input file path did not contain
    a file suffix
    :return: a string of the file name with the suffix(es) removed
    '''

    # first, confirm that the input file_path is a Path object
    is_pathclass(file_path, exit_if_false=False)

    # create a list of any suffixes detected by pathlib
    file_suffixes = file_path.suffixes

    # check for how many suffixes were found in the file name

    # if no suffixes were found, then the input was likely a directory, not a file
    if len(file_suffixes) == 0:

        # create a warning message
        warning_msg = [f'The input file path:\n   {file_path}\n']

        # update warning message based on whether input file_path is an existing directory or not
        if file_path.is_dir():
            warning_msg.append('is a directory and therefore does not contain a suffix.\n')
        else:
            warning_msg.append('did not contain any detectable file suffixes and does not appear '
                               'to be an existing directory.\n')

        # regardless of whether an existing directory or not, just return the string
        warning_msg.append('The name of the input file/directory was returned as a string.')

        # only print this warning if the verbose param is set to True (False by default)
        if verbose:
            print(''.join(warning_msg))
        else:
            pass

        return file_path.name

    # if a single suffix was found, then use the .stem attribute to return the file name
    elif len(file_suffixes) == 1:
        return file_path.stem

    # if there was more than one suffix found in the file name...
    else:

        # reduce the list of suffixes to valid file suffixes only
        valid_file_suffixes = [suff for suff in file_suffixes if suff in common_suffixes]

        # check what was filtered out (i.e., not true file suffixes)
        removed_file_suffixes = set(file_suffixes).difference(set(valid_file_suffixes))
        if len(removed_file_suffixes) > 0:
            if verbose:
                print(f'The following strings were erroneously detected as file suffixes in the file '
                      f'name {file_path.name}:\n'
                      f'   {removed_file_suffixes}\n'
                      f'These were removed from the list of valid file suffixes for this file.')
            else:
                pass
        else:
            pass

        # remove the valid file suffixes from the input file name
        valid_file_suffix_re = '|'.join(valid_file_suffixes)  # bit worried this won't fly in 3.12 (not raw string)
        file_without_suffixes = re.sub(valid_file_suffix_re, '', file_path.name)

        # return the file name with the valid file suffixes removed
        return file_without_suffixes

def replace_file_ext(file_path, output_ext, create_file=False, replace_file=False, output_dir=None):
    '''
    Replace the file's file extension with a different file extension. Does not
    reformat the file, it will only change its suffix.

    :param file_path: path to the file that should have its file extension
    replaced with the value provided to output_ext
    :param output_ext: the file extension to use to replace the file
    extension of the input file provided to file_path
    :param create_file: True / False; if False (default), will return a file path
    object with the output path and file extension, but will not create or alter
    any existing files
    :param replace_file: True / False; if False (default) the input file provided to
    the file_path parameter will not be replaced by the output file path with the
    output_ext file extension; option is only valid if create_file=True
    :param output_dir: the output directory to which to write the output file
    with the output_ext file extension; if None, the output file will be written
    to the same location as the input file
    :return: output file path with the output_ext file extension
    '''

    ## DEFINE INPUT / OUTPUT FILE PATHS ##

    assert is_pathclass(file_path, exit_if_false=True)

    if output_dir is None:
        output_dir = file_path.parent
    else:
        assert is_pathclass(output_dir, exit_if_false=True)

    ## CHECK FOR LEADING . IN FILE EXTENSION ##

    # if the output_ext arg already has a leading ., pass
    if output_ext.startswith('.'):
        pass

    # if the output_ext arg doesn't have a leading ., add one
    else:
        # must be formatted this way to use with pathlib's .with_suffix()
        output_ext = '.' + output_ext

    ## STRIP INPUT FILE NAME'S FILE EXTENSION ##

    # use utilities.py function strip_file_ext() to remove file extension, even
    #   if the input file has multiple suffixes (e.g., fastq.gz)
    input_file_basename = strip_file_ext(file_path=file_path)

    ## CREATE OUTPUT FILE PATH WITH OUTPUT_EXT AS SUFFIX ##

    output_path = (output_dir / input_file_basename).with_suffix(output_ext)

    ## IF A FILE SHOULD BE CREATED USING THE OUTPUT PATH W/ OUTPUT FILE EXTENSION ##

    if create_file:

        # REPLACE INPUT FILE #

        # if the output file should replace the input file...
        if replace_file:

            # rename the input file with the output file path
            file_path.replace(output_path)

        # CREATE COPY OF INPUT FILE WITH OUTPUT FILE NAME #

        # if the output file shouldn't replace the input file...
        else:

            # copy the content from the input file to this new file path
            with open(file_path, 'rt') as f_in:
                with open(output_path, 'wt') as f_out:
                    shutil.copyfileobj(f_in, f_out)

    ## OTHERWISE, DO NOTHING HERE ##
    else:
        pass

    ## RETURN NEW FILE PATH ##
    return output_path


def create_file_list(file_input, file_regex=[SEQ_FILE_RE, GZIP_REGEX]):
    '''
    Converts a directory into a list of its file contents. If the input provided is
    already a list of files, then it simply returns the input unchanged.

    :param file_input: item received by the parent function as the input of files;
    either as a directory from which this function will extract a list of files or
    a list of files, in which case nothing will be done and the list will be returned
    without changes.
    :param file_regex: the regular expression to use to filter out any files; this
    is primarly to exclude files that are not sequence files, so by default it is a
    sequence file regular expression; if a list is provided, the regex in the list
    will be joined so that it will match either (|) of the regex
    :return: a list of files
    '''

    # create a single regex if list is provided
    if isinstance(file_regex, list):
        file_filter_regex = '|'.join(file_regex)
    else:
        file_filter_regex = file_regex

    # is the input of files a list or a PosixPath object?

    # if the input is a Path object, create a (filtered) list of its contents
    if is_pathclass(file_input, exit_if_false=False):

        # if the input file is a single file, just create a list with one item
        if file_input.is_file():
            return [file_input]

        # if the input is a directory, create a list of its contents that match the input regex
        else:
            output_file_list = [file for file in file_input.glob('*') if re.search(file_filter_regex, file.name, re.I)]

    # if the input is not a Path object but a list, ensure it is filtered and return its contents
    else:
        output_file_list = [file for file in file_input if re.search(file_filter_regex, file.name, re.I)]

    # return the filtered file list; which is a list of PosixPath objects matching the provided file regex
    return output_file_list

# create a function that will detect the available system memory
def get_available_memory(memory_units='GB'):
    '''
    Get the available system memory of the current system.

    :param memory_units: the units of measurement that the returned memory will be in
    :return: the available system memory, in memory_units, rounded down to the
    nearest whole number
    '''

    # check input value for memory_units to ensure its an accepted value
    accepted_memory_units = {
        'bytes': 'B',
        'megabytes': 'MB',
        'gigabytes': 'GB',
    }
    if memory_units in accepted_memory_units:
        pass
    else:
        err_msg = (f'The input value for the parameter memory_units in the '
                   f'function {get_available_memory.__name__} is not among '
                   f'the accepted values for this parameter:\n'
                   f'   accepted values     = {list(accepted_memory_units.values())}\n'
                   f'   invalid input value = {memory_units}\n')
        exit_process(err_msg)

    # define a function to convert the input units of memory to the desired output units, memory_units
    def convert_system_memory(sys_mem, units_in, units_out):

        # convert from bytes...
        if units_in == 'B':
            if units_out == 'B':     # to bytes
                return sys_mem
            elif units_out == 'MB':  # to megabytes
                return sys_mem*1e-6
            elif units_out == 'GB':  # to gigabytes
                return sys_mem*1e-9
            else:                    # return error if not bytes, megabytes, or gigabytes
                pass

        # convert from megabytes...
        elif units_in == 'MB':
            if units_out == 'B':     # to bytes
                return sys_mem/1e-6
            elif units_out == 'MB':  # to megabytes
                return sys_mem
            elif units_out == 'GB':  # to gigabytes
                return sys_mem/1000
            else:                    # return error if not bytes, megabytes, or gigabytes
                pass

        # convert from gigabytes...
        elif units_in == 'GB':
            if units_out == 'B':     # to bytes
                return sys_mem*1e-9
            elif units_out == 'MB':  # to megabytes
                return sys_mem*1000
            elif units_out == 'GB':  # to gigabytes
                return sys_mem
            else:                    # return error if not bytes, megabytes, or gigabytes
                pass

        # invalid unit of measurement, print error below
        else:
            pass

        err_msg = (f'Invalid input units provided to the function {convert_system_memory.__name__}, '
                   f'a subfunction within {get_available_memory.__name__}\n'
                   f'   accepted values     = {list(accepted_memory_units.values())}\n'
                   f'   invalid input value = {units_in}\n')
        exit_process(err_msg)

    # convert the available system memory in bytes to the unit of memory_units
    avail_sys_mem = convert_system_memory(
        sys_mem=psutil.virtual_memory().available,  # get available system memory
        units_in='B',                               # psutil will always return in bytes
        units_out=memory_units,                     # return value as memory_units
    )

    bbduk_fmt_unit = memory_units.lower()[0]

    return bbduk_fmt_unit, math.floor(avail_sys_mem)
    return bbduk_fmt_unit, math.floor(avail_sys_mem)