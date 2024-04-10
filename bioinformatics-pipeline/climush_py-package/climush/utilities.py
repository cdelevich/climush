import subprocess, re, sys, pathlib, shutil, json, tomlkit
from pathlib import Path
from datetime import datetime
from functools import wraps
import pandas as pd
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
        msg = f'The provided path, {file_path}, is not a valid Path object.\n'
        print(msg)
        if exit_if_false:
            return exit_process(message=msg)
        else:
            return False

# alert if multiple files matching the search pattern in the file path are returned
# prompt for user input if multiple or no files matching the pattern are returned
# moved here because used in continue_to_next()
def flag_multiple_files(file_path, search_for):
    assert is_pathclass(file_path, exit_if_false=False)

    result = list(file_path.glob(search_for))

    if len(result) == 1:
        return result[0]
    elif len(result) == 0:
        all_files = file_path.glob('*')
        print(f'No files matching the pattern \'{search_for}\' were detected in the file path: {file_path}. Do '
              f'you mean any of these files in this directory?')
        which_file = prompt_print_options([all_files, 'none of these (exit)'])
        if which_file in all_files:
            return which_file
        else:
            return exit_process(message=f'The response \'none of these\' was chosen when searching for the correct'
                                        f'file matching the pattern: {search_for}')
    else:
        print(f'{len(result)} files matching the pattern \'{search_for}\' were detected in the file path: '
              f'{file_path}. Please type the number corresponding to the correct file to use:')
        which_file = prompt_print_options([result, 'none of these (exit)'])
        if which_file == 'none of these (exit)':
            return exit_process(message=f'The response \'none of these\' was chosen when searching for the correct'
                                        f'file matching the pattern: {search_for}.')
        else:
            return which_file

# activate next script in python
def continue_to_next(this_script, config_dict):
    '''
    Continue to the next step of the pipeline.

    Will activate the script from the CLI with default
    parameters. If you do not want to use the default
    parameters, opt out of automated continuation by
    responding to the CLI prompt that precedes the use
    of this function.
    :param this_script: use __file__ always
    :return: None, assembles and runs a shell command.
    '''
    automate_dict = config_dict['automate']
    auto = automate_dict['run_all']
    to_run = automate_dict['run_some']
    to_exclude = automate_dict['exclude']

    current_script = Path(this_script).stem
    current_num = int(current_script.split('_')[0])

    def format_num_for_search(current_script_number):
        if current_script_number > 9:
            next_num_search = str(current_script_number + 1)
        else:
            next_num_search = '0' + str(current_script_number + 1)
        return next_num_search

    next_num_search = format_num_for_search(current_num)
    next_script = flag_multiple_files(file_path=this_script.parent, search_for=f'{next_num_search}*')

    if (auto) or ((current_num + 1) in to_run):
        print(f'\n\nRunning next step, {next_script.name}...\n')
        subprocess.run(['python3', next_script])
        return None
    elif (current_num + 1) in to_exclude:
        i = 2
        while (current_num + i) in to_exclude:
            i += 1
        next_num_search = format_num_for_search(current_num + 1)
        next_wanted = flag_multiple_files(file_path=this_script.parent, search_for=f'{next_num_search}')
        print(f'\n\nRunning next step, {next_wanted.name}...\n')
        subprocess.run(['python3'], next_wanted)
        return None
    else:
        to_next = input(f'The script {current_script} has completed. Would you like to continue '
                        f'to the next step in the pipeline, {next_script.stem}?[Y/N]\n'
                        f'\n')
        if re.search(AFFIRM_REGEX, to_next, re.I):
            print(f'\n\nRunning next step, {next_script.name}...\n')
            subprocess.run(['python3', next_script])
            return None
        else:
            return print(f'\n\nExiting the completed step, {current_script}.\n')

# log progress of bioinformatics pipeline
def log_progress(file_map, run_name):
    log_file = file_map['pipeline-output']['summary'] / f'log_{run_name}.json'

    log_dict = {'run_name': run_name,
                'error': {'script':'04_quality-filtering',
                          'function':''}
                }

    with open(log_file, 'at') as log_out:
        log_out.write(json.dumps(log_dict))



# exit current script due to error, save script and timestamp of where error occurred
def exit_process(message, config_section='error.message'):
    script_name = sys.argv[0]  # unsure if will get name of script it is executed in or the one it is compiled in
    exit_time = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    # write_to_config(config_section=config_section, script=script_name, timestamp=exit_time, details=message)
    print(f'Exiting {script_name}...\n')
    return sys.exit()

# recursive function that will seek out yes/no/quit response continuously until achieved
def prompt_yes_no_quit(message):
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
def prompt_print_options(option_list):
    option_list.sort()  # sort option list for easier user experience
    option_list.append('quit')  # add option to quit, listed last

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
def prompt_multiple_files(file_path):
    print(f'\nWARNING: Multiple files were detected in the '
          f'{file_path.stem} folder. Please type the number '
          f'corresponding to the correct file to use:')
    file_list = [file.stem for file in file_path.glob('*') if not re.search(HIDDEN_FILE_REGEX, file.stem)]
    return prompt_print_options(file_list)

# provide options for the sequencing platform if it cannot be otherwise detected
def prompt_sequencing_platform(sample_id):
    print(f'\nWARNING: The sequencing platform could not be inferred from the sample: {sample_id}. Please type '
          f'the number corresponding to the correct sequencing platform from the options below. ')
    platform = prompt_print_options([SEQ_PLATFORM_OPTS, 'multiple'])
    if platform == 'multiple':
        print(f'If you have a combination of sequences from multiple platforms, you will need to either:\n'
              f'(1) manually sort the files into their sequencing platform directories\n'
              f'(2) use the file renaming script, rename_sequence_files.py, to rename the files to match '
              f'the file naming convention.\n')
        return

# print indented list
def print_indented_list(print_list):
    print_list[0] = '\t' + print_list[0] # add leading tab for first item printed from list, otherwise adds after line break
    formatted_list = '\n\t'.join(print_list)
    return print(formatted_list)

# run shell command and save stdout and stderr to file
def run_subprocess(cli_command_list, dest_dir):
    '''
    Run a subprocess, saving the output and error to a log file.

    Takes the command line (CL) argument assembled as a list, and runs
    the subprocess, saving the standard in (stdin), standard out (stdout),
    and standard error (stderr) from the process. The output of stderr and
    stdout are written to a log file in the destination directory
    :param cli_command_list: the list of the components of the command line
    argument to run through subprocess; each argument should be separated
    out into a list, rather than a single string.
    :param dest_path: directory to write the file to; should be a Path
    object
    :return: None
    '''
    assert is_pathclass(dest_dir)

    assert isinstance(cli_command_list, list) and len(cli_command_list) > 1, print(f'Command line process to run '
                                                                             f'should be in a list of its components. '
                                                                             f'Try shlex.split() if unsure how to '
                                                                             f'compose the list. ')

    # run_cmd = subprocess.run(cli_command_list, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    # out, err = run_cmd.communicate()
    program = cli_command_list[0]

    if re.search('.\..', program):
        program = program.split('.')[0]

    run_cmd = subprocess.run(cli_command_list, capture_output=True)

    out_path = dest_dir / f'{program}.out'
    with open(out_path, 'at') as fout:
        as_str = run_cmd.stdout.decode('utf-8')
        if int(out_path.stat().st_size) == 0:
            fout.write(as_str)
        else:
            fout.write(as_str.split('\n')[1] + '\n')

    err_path = dest_dir / f'{program}.err'
    with open(err_path, 'ab') as fout:
        fout.write(run_cmd.stderr)

    temp_file = dest_dir / f'{program}.temp'

    if len(run_cmd.stderr) == 0:
        pass
    else:
        if not temp_file.is_file():
            continue_ok = input(f'Running {program} produced an error. Please review the output in {err_path.name}. '
                                f'Would like to continue despite this error? [Y/N]')
            if re.search(continue_ok, AFFIRM_REGEX, re.I):
                with open(temp_file, 'wt') as fout:
                    fout.write('stop_prompt')
                return None
            else:
                return exit_process(message=f'Running {program} produced an error. See {err_path.name} for details.')

    return None

def func_timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = datetime.datetime.now()
        func_output = func(*args, **kwargs)
        end_time = datetime.datetime.now()
        runtime = str(end_time - start_time).split('.')[0]  # round seconds down
        print(f"{func.__name__} was executed in {runtime}.\n")
        return func_output
    return wrapper()

def get_seq_platform(fastx_file, delim):
    if isinstance(fastx_file, str):
        sample_id = fastx_file
    elif is_pathclass(fastx_file):
        sample_id = fastx_file.stem
    else:
        print(f'ERROR. An error occurred when trying to detect the sequencing platform from the '
              f'sample ID: {fastx_file}\n')


    platform = sample_id.split(delim)[0]
    if platform in SEQ_PLATFORM_OPTS:
        return platform
    else:
        pass

def import_mapping_df(df_path):
    '''
    Import .csv, .txt, or .xlsx table as a dictionary.

    Read in a mapping file for demultiplexing. Can accomodate the file formats
    .xlsx, .csv, and .txt. Will output a dictionary, where the key is the name
    of the tab in the dataframe and the value is the dataframe in that tab. Will
    always return a dictionary, although only .xlsx files will have tabs. Formats
    the name of the keys in the output dictionary to be 'pool1' or 'pool2'; the
    number of the pool is inferred from the name of the tab (.xlsx) or the name
    of the file (.csv, .txt). Returned as dictionary so that the output can be
    handled in the same way, regardless of file type (e.g., loop through tabs
    even if a .csv, which will have a single tab).
    :param df_path: path to the dataframe
    :return: dictionary, where the key is 'pool1' or 'pool2', and the value is the
    dataframe belonging to that tab, or in the case of a .csv or .txt file, the
    entirety of that dataframe file.
    '''

    POOL_NUM_RE = '(?<=pool).?(\d)'

    if re.search('^\.x', df_path.suffix):  # if an excel file
        mapping_tabs = pd.read_excel(df_path, sheet_name=None)  # need to set sheet_name to None to get all tabs
        old_tab_names = list(mapping_tabs.keys())  # make list of old names, otherwise cannot update key names in loop
        for tab in old_tab_names:
            try:  # format the name of the tab to be uniform across all tabs/dataframes
                pool_num = re.search(POOL_NUM_RE, tab, re.I).group(0)
            except AttributeError:
                print(f'The pool number could not be inferred from tab {tab} in the mapping file {df_path.name}. '
                      f'Please type the correct pool number for this file: ')
                pool_num = prompt_print_options(['1', '2'])  # choose from 1 or 2 (or quit, built into function)
            mapping_tabs[f'pool{pool_num}'] = mapping_tabs.pop(tab)  # replace old tab (key) with reformatted one
    elif re.search('^\.c|^\.txt$', df_path.suffix):  # I think you can read in .txt and .csv files the same way?
        try:  # try to get the pool number from the file name
            pool_num = re.search(POOL_NUM_RE, df_path.name, re.I).group(0)
        except AttributeError:  # if there's no detected pool number in the name, prompt user to specify one
            print(f'The pool number could not be inferred from the file name of the mapping file {df_path.name}. '
                  f'Please type the correct pool number for this file: ')
            pool_num = prompt_print_options(['1', '2'])  # choose from 1 or 2 (or quit, built into function)
        mapping_tabs = {f'pool{pool_num}': pd.read_csv(df_path)}  # make dict to match format from .xlsx
    else:
        print(f'ERROR. The file format {df_path.suffix} of the mapping file {df_path.name} is not a recognized '
              f'file type. Accepted file types are: \'.xlsx\', \'.csv\', and \'.txt\'.\n')
        sys.exit()

    return mapping_tabs


#######################
# FILE PATHS ##########
#######################

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
                continue
                fate = ''
        else:
            continue
            fate = ''

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
            return exit_process(message=msg)
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
        new_dir.mkdir(exist_ok=True)
        return new_dir
    elif (parent_dir is None):
        return print(f'\nThe provided file path of the new directory is not absolute. Please '
                           f'provide a parent_dir to specify the path of the new directory, and '
                           f'retry.\n')
    else:
        new_path = Path(parent_dir / new_dir)
        new_path.mkdir(exist_ok=True)
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

# sort the sequences provided in the 'sequences' input directory

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


def check_for_input(file_dir, seq_platform=None, file_ext=SEQ_FILE_GLOB):
    '''
    Checks if there are input files for the process.

    :param file_dir: the directory path to check for files and to check
    whether the path exists
    :param seq_platform: the type of sequencing files to look for ['illumina',
    'sanger', 'pacbio']; defaults to None, meaning it will return True if any
    non-empty directory or file is located
    :param file_ext: the expected file extension of the files to search for,
    using asterisk '*' for wildcard (used in glob)
    :return: two items; [1] Boolean T/F, whether there is a directory in this
    path that contains sequencing files; [2] a list of the files matching the
    file extension in the directory, if present
    '''
    assert is_pathclass(file_dir)

    # I decided to leave this out so you could put any regex in (i.e., multiplexed pacbio files will start with \d{4}
    # accepted_seq_input = [None, 'illumina', 'sanger', 'pacbio']
    # if not seq_platform in accepted_seq_input:
    #     print(f'The provided sequencing platform, {seq_platform}, is not one of the accepted '
    #           f'values: {accepted_seq_input}. Please retry with one input value from this list '
    #           f'as a string, or multiple values as a list of strings.\n')
    #     return sys.exit()

    if seq_platform is None:
        seq_re = '.+?'  # will return all files in directory
    elif isinstance(seq_platform, list):
        seq_re = '|'.join(seq_platform)  # returns only those in provided list
    else:
        seq_re = seq_platform  # returns only those of the provided platform

    if file_dir.is_dir():  # check that input is a directory
        if count_files(file_path=file_dir, search_for=file_ext) > 0:  # confirm it is not empty
            file_list = [f for f in file_dir.glob('*') if re.search(seq_re, f.name, re.I)]  # check for specific files
            if len(file_list) > 0:
                return True, file_list  # return true if files matching criteria are found
            else:
                print(f'The directory {file_dir} contains files, but none that match the platform-specific search '
                      f'criteria: {seq_platform}.')
                return False, file_list  # return false if not
        else:
            file_list = []  # if input is a directory, but it is empty...
            print(f'The directory {file_dir} exists, but is empty.\n')
            return False, file_list  # return false and empty list
    else:
        file_list = []  # if input is not a directory
        print(f'The directory {file_dir} does not exist.\n')
        return False, file_list  # return false and empty list

# get name of previous script

#######################
# TEXT MANIPULATION ###
#######################

# add prefix to file name
def add_prefix(file_path, prefix, dest_dir, action='rename', f_delim='_'):
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
    :return: Path object with new file name
    '''

    assert is_pathclass(file_path, exit_if_false=False)

    def prefix_single(file_name):
        old_name = file_name.name
        platform_present = re.search('illumina|pacbio|sanger', old_name, re.I)
        try:
            # if there's a match, put prefix before the platform, remove all else
            location = platform_present.span()[0]
            new_name = prefix + f_delim + old_name[location:]
        except:
            # if there's no match, just add prefix to start of file name
            new_name = prefix + f_delim + old_name

        # replace old name with new name
        if action == 'rename':
            new_path = file_name.rename(dest_dir / new_name)
        elif action == 'copy':
            new_path = dest_dir / new_name
            shutil.copy(file_path, new_path)
        else:
            new_path = dest_dir / new_name

        return new_path


    if file_path.is_dir():
        new_path_list = []
        for file in file_path.glob('*'):
            new_path_list.append(prefix_single(file))
        return new_path_list
    elif file_path.is_file():
        return prefix_single(file_path)
    else:
        msg = f'FAILURE. Unrecognized input file; cannot recognize as either a directory or file.\n'
        return exit_process(msg)

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
    script_name_nonum = re.sub('^\d+_(?=\w)','',Path(script_name).stem)
    completed_name_nonum = OUTPUT_DIRS[script_name_nonum]
    completed_name = f'{script_num}_{completed_name_nonum}'

    if remove_num:
        completed_name = re.sub('^\d+_(?=\w)','', completed_name)

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

    if re.search('^no-', prefix_const, re.I):  # if provided prefix starts with no, remove no
        return re.sub('^no-', '', prefix_const, re.I)
    else:  # if it does not start with no, add no
        return 'no-' + prefix_const

def rename_read_header(fasta_file, header_delim=';'):
    # fasta_file = next(post_itsx.glob('*.fast*'))
    parent_dir = fasta_file.parent.name

    if re.search('itsx', parent_dir, re.I):  # if they are post-itsx reads...

        READ_COUNT_RE = '(?<=size=)[0-9]{1,}'
        READ_LEN_RE = '[0-9]{1,}(?=\sbp)'
        READ_ID_RE = '[0-9]{1,}(?=\/ccs)'
        READ_REGION_RE = '(?<=\|\w\|).+(?=\sExtracted)'
        SAMPLE_ID_RE = '(?<=\w_)(pacbio|sanger|illumina).+?(?=\.)'

        # REPLACE WITH REMOVE PREFIX??
        sample_id = re.search(SAMPLE_ID_RE, fasta_file.name).group(0)

        updated_records = []
        no_update_count = 0
        for record in SeqIO.parse(fasta_file, 'fasta'):
            header = record.description
            try:
                read_id = re.search(READ_ID_RE, header).group(0)
                region = re.search(READ_REGION_RE, header).group(0)
                read_count = re.search(READ_COUNT_RE, header).group(0)
                read_length = re.search(READ_LEN_RE, header).group(0)

                updated_header = header_delim.join([f'{sample_id}_{read_id}', region, f'region_len={read_length}bp',
                                                    f'full-len_copies={read_count}'])

                record.id = updated_header
                record.name = f'{sample_id}_{read_id}'
                record.description = ''
                updated_records.append(record)
            except:
                no_update_count += 1
                updated_records.append(record)

        SeqIO.write(updated_records, fasta_file, 'fasta')

        if no_update_count == 0:
            return print(f'SUCCESS. All reads ({len(updated_records)}) were renamed '
                         f'in {fasta_file.name}.\n')
        elif no_update_count == len(updated_records):
            return print(f'FAILURE. None of the reads were renamed in {fasta_file.name}. This can '
                         f'happen if the fasta file has already been renamed by this function, '
                         f'so open the fasta file to check.\n')
        else:  # if some reads were renamed, but not all
            return print(f'FAILURE. {no_update_count} sequence headers could not be renamed due to '
                         f'a missing data field for reads in {fasta_file.name}.\n')
    else:
        return print(f'UNDER CONSTRUCTION. I haven\'t yet updated this function to work with any other '
                     f'fasta files than those produced after ITSx.\n')

def escape_path_special(file_path):
    if is_pathclass(file_path):
        file_path = str(file_path)

    if re.search('\(', file_path):
        file_path = re.sub('\(', '\(', file_path)

    if re.search('\)', file_path):
        file_path = re.sub('\)', '\)', file_path)

    return file_path

#######################
# CONFIGURATION #######
#######################

# define function to simplify import of configuration files
def import_config_as_dict(file_path, file_handle, config_section='all'):
    '''
    Import a configuration file and return as a dictionary.

    Checks the input path to ensure that the filepath exists; if
    not, it will stop the process and exit the entire script.
    Uses a custom parser to preserve the case of the text in the
    configuration file (i.e., uppercase letters remain uppercase,
    lowercase remain lowercase). Creates a dictionary where the
    sections of the configuration file are the primary keys, if it
    is a nested dictionary. If not a nested dictionary, will be
    returned as an un-nested dictionary where the keys are the
    section's keys.
    :param file_handle: standard file handle used for the config file; will
    default to the constants that work if the file structure of the container
    is maintained
    :param config_section: option to output only one section of the configuration
    file.
    :return: dictionary of values in the configuration file
    '''
    path_to_file_abs = flag_if_not_path_exists(file_path, absolute=True)

    config_file_list = list(path_to_file_abs.glob(f'*{file_handle}*'))

    if len(config_file_list) > 1:
        msg = f'\nERROR. Multiple configuration files matching the provided file handle, {file_handle}, were located '\
              f'in this directory:\n'\
              f'\t{[file.name for file in config_file_list]}\n'\
              f'Please include more information in the file handle input, then retry.\n'
        print(msg)
        # return exit_process(message=msg)
    else:
        config_file = config_file_list[0]

    with open(config_file, 'rt') as fin:
        config_dict = tomlkit.parse(fin.read())

    if config_section == 'all':
        return config_dict
    else:
        if config_section in config_dict.keys():
            return config_dict[config_section]
        else:
            msg = f'The provided configuration file section header, {config_section}, is not a valid header for ' \
                  f'the {file_handle} configuration file. Please enter the number of the correct header that you ' \
                  f'want:'
            print(msg)
            correct_header = prompt_print_options(list(config_dict.keys()))
            return config_dict[correct_header]

# write or update settings to the configuration file
# NOT WORKING RIGHT NOW
# def write_to_config(config_section, output_path, update_obj=False, overwrite=True, **kwargs):
#     # import the pipeline configuration file directly
#     config_section = 'error_message'
#     config_dict = import_config_as_dict(output_path, file_handle)
#
#     # confirm that the provided section is in the configuration file
#     if not config_section in config_dict.keys():
#         msg = f'ERROR: Could not locate the provided section \'{config_section}\' in the pipeline configuration '\
#               f'file.'
#         print(msg)
#         return exit_process(message=msg)
#
#     # update the configuration file with the provided key-value pairs
#     for k, v in kwargs.items():
#         if k in config_dict[config_section].keys():  # if there's already a key for this kwarg...
#             if (config_dict[config_section][k] == '') or overwrite:  # if it has no value or if we want to overwrite...
#                 config_dict[config_section][k] == v  # update/overwrite current key with new key value pair
#             else:  # if there is something already there and we do not want to overwrite it...
#                 continue  # continue to next key-value pair
#         else:  # if there's no key matching in this section, add key-value pair without further inspection
#             config_dict[config_section][k] == v
#
#     # write out this updated configuration dict to replace the previous configuration file
#     config_filename = flag_multiple_files(output_path, search_for=f'*{PIPELINE_CONFIG_HANDLE}*{CONFIG_FILETYPE}')
#
#     with open(config_filename, 'wt') as fout:
#         fout.write(tomlkit.dumps(config_dict))
#
#     return None

# locate relevant mapping file for demultiplexing
# moved down because it requires the import_config_as_dict function
def find_mapping_file(path_to_mapping, path_to_plex, path_to_config, **kwargs):
    '''
    Locate the barcode mapping file(s) for demultiplexing.

    Attempts to automatically detect the location of the mapping file
    corresponding to the multiplexed files located in the sorted
    directory of sequences to demultiplex. It inspects the name of the files
    that were sorted into the needs_demux folder of the sequences folder, and
    identifies the name of the sequencing run. It then will go into the mapping
    file directory in the config folder, mapping-files, and looks for the mapping
    file that shares that same sequencing run name. If multiple mapping files
    are detected for a single sequencing run, the user will be prompted to select
    the file that should be used. If no matching mapping file is detected, it
    will throw an error and exit the pipeline.
    :param path_to_mapping: defaults to expected location within structured
    file directories; otherwise can specify the path as a Path object
    :param path_to_plex: defaults to expected location within structured
    file directories; otherwise can specify the path to the mapping file(s) as
    a Path object
    :param kwargs: specify the bioinfo-settings configuration file dictionary
    if one is already loaded into the environment; if none is provided, the
    function will open the configuration file to create a new instance of the
    dictionary; adding a dictionary, if available, will likely be quicker so it
    is recommended to do so if possible
    :return: if found, returns mapping file(s) as a list
    '''
    assert is_pathclass(path_to_mapping)
    assert is_pathclass(path_to_plex)

    # get the multiplexed file's queue ID, which the UO sequencing core uses as the seq run identifier in the filename
    multiplex_ids = [file.stem.split('.')[0] for file in path_to_plex.glob('*.fast*')]

    # look up the multiplex IDs in the bioinfo-settings configuration file for its corresponding CliMush run name
    kwargs_keys = list(kwargs.keys())
    if kwargs_keys > 0:
        if kwargs_keys > 1:
            msg = f'ERROR: Too many keyword arguments (**kwargs) used. Only one kwarg is accepted and should have the '\
                  f'name of the configuration dictionary as the value to a key with any name.\n'
            print(msg)
            return exit_process(message=msg)
        else:
            config_multiplex_dict = kwargs[kwargs_keys[0]]['pacbio-multiplex-ids']
    else:
        config_multiplex_dict = import_config_as_dict(file_handle=path_to_config)['pacbio-multiplex-ids']

    # get the CliMush run name(s) and create regex to search for any of these at start of mapping file name
    climush_ids = [config_multiplex_dict[q] for q in multiplex_ids]

    # find the corresponding barcode mapping files for the run name(s) in the barcode-mapping dir in config
    mapping_file_names = []
    for id in climush_ids:
        match = list(path_to_mapping.glob(f'*{id}*'))
        if len(match) == 1:
            mapping_file_names.append(match[0])
        elif len(match) == 0:
            msg = f'ERROR: No mapping files matching the CliMush run name {id} were located in {path_to_mapping}. ' \
                  f'Check the mapping-files directory in the config folder, then retry.\n'
            print(msg)
            return exit_process(message=msg)
        else:
            msg = f'ERROR: {len(match)} mapping files were located for CliMush run name {id}. Please select the number ' \
                  f'of the correctly corresponding barcode mapping file for this sample:\n'
            print(msg)
            mapping_file_names.append(prompt_print_options(match))

    return mapping_file_names

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
        if follows_naming_convention(filepath) or re.search('^[A-Z]+', filepath.stem):
            return True
        else:
            return False

    def is_pacbio(filepath):
        if re.search('^pacbio|^\d{4}', filepath.stem):
            return True
        else:
            return False

    def is_illumina(filepath):
        if re.search('^illumina|^[A-Z]+', filepath.stem):
            return True
        else:
            return False

    def is_sanger(filepath):
        if re.search('^sanger', filepath.stem):
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

