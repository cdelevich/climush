import argparse, configparser, itertools, os, shutil, re, sys
from pathlib import Path
import pandas as pd
from climush.config import import_config_as_dict

# command line arguments
parser = argparse.ArgumentParser(description='Rename CliMush raw sequence files.')

path_help = 'Path to the directory containing the sequences to rename, relative to the location of ' \
            'rename_sequence_files.py script.'
parser.add_argument('path', help=path_help)  # positional

inipath_help = 'Relative or absolute path to the .ini file for this set of sequence files, which should be updated ' \
               'by the user prior to running this script; if a path is not provided, it will look for a file with the ' \
               'extension \'.ini\' in the same directory as the input path provided.'
parser.add_argument('-inipath', help=inipath_help, required=False)

tableout_help = 'Whether to save a conversion table showing the old and new file names. Output table will not ' \
                'be saved unless argument used in command line. '
parser.add_argument('--table_out', help=tableout_help, required=False, action=argparse.BooleanOptionalAction)

nocopy_help = 'Whether to save a copy of the files with their original file names in a zip file ' \
              'before renaming. Tends to be time consuming for Illumina sequencing files. If you do not ' \
              'want copies made, include the argument in the command line. Program will confirm either way ' \
              'before proceeding.'
parser.add_argument('--no_copy', help=nocopy_help, required=False, action=argparse.BooleanOptionalAction)

ext_help = 'File extension of the sequence files; if not included, will default to fasta.gz'
parser.add_argument('--ext', default='.fastq.gz', help=ext_help)

colldate_help = 'The general collection date used for this group of sequences, formatted as YYYY-MM; for spring collections, ' \
                'use 202#-05 and for fall use 202#-10. if a collection date is not provided, the collection date will be ' \
                'inferred from the name of the configuration (.ini) file as long as it follows the standard naming ' \
                'convention for configuration files (i.e., YYYY-MM date is included in the file name). Error will be ' \
                'raised if the program cannot locate a date with this format in the configuration file name, and none is ' \
                'provided with the \'--coll_date\' CLI argument.'
parser.add_argument('--coll_date', required=False, help=colldate_help)

args = parser.parse_args()

def positive_response(input_response):
    if re.search('^Y', input_response, re.I):
        return True
    else:
        return False
def check_input(cli_input):
    '''
    Confirms that the format of the input from the command line is
    formatted correctly. If it is a directory, it must have a trailing
    '/' and if it is a file extension, it must have a leading '.'
    :param input: input variable
    :return: if not formatted correctly, returns reformatted variable;
    otherwise will return original input variable
    '''

    if cli_input is None:
        return cli_input
    else:
        if (os.path.isdir(cli_input)) and not (cli_input.endswith('/')):
            cli_input = cli_input + '/'
        elif re.search('^\d', cli_input):
            if not (re.search('\d{4}-\d{2}', cli_input)) and not (os.path.isdir(cli_input)):
                print(f'\nCollection date provided in the command line is not '
                      f'formatted correctly. Must be YYYY-MM.')
                sys.exit()
        else:
            if not cli_input.startswith('.'):
                cli_input = '.' + cli_input

    return cli_input

table_out = args.table_out
colldate_group = check_input(args.coll_date)
path_input = check_input(args.path)
file_ext = check_input(args.ext)

if args.inipath is None:
    ini_types = list(Path(path_input).glob('*.ini'))
    if len(ini_types) > 1:
        print(f'\nERROR: Multiple .ini files ({len(ini_types)}) were located in the provided input directory '
              f'{path_input}. Please specify the correct configuration file in the '
              f'command line with argument \'-inipath\' or remove all other .ini files in this directory.\n')
        sys.exit()
    else:
        path_init = ini_types[0]
        print(f'\nUsing configuration file: {path_init}.')
else:
    path_init = args.inipath

if args.no_copy is True:
    cont_input = input('A copy of the files with their original file names will *not* be created, so it '
                       'is advised that you have your own backup of these files with their original file names. '
                       'Do you wish to continue? [Y/N]')
else:
    cont_input = input('A copy of the files with their original file names will be created, and may '
                       'be time intensive. Do you wish to continue? [Y/N]')

if positive_response(cont_input):
    no_copy = args.no_copy
else:
    print('To save a copy of the files with their original file names, remove \'--no_copy\' '
          'from the command line input.\n'
          'To prevent a copy of the files with their original file names from being saved, '
          'include \'--no_copy\' from the command line input.')
    sys.exit()

## FOR TESTING/TROUBLESHOOTING
# table_out = True
# colldate_group = None
# path_input = './illumina_renamed-files/illumina_2023-05_test/'
# path_init = f'{path_input}!file-rename_config_2023-05_soil-litter.ini'
# file_ext = '.fastq.gz'
# no_copy = True

def is_valid_entry(string):
    '''
    Tests if a string contains data.

    Uses the VALID_ENTRY regex to search an input string for usable
    data. Checks whether the word 'NA' or 'None' is detected as a
    whole word, and not a part of another string, case insensitive.
    If NA/None or empty string detected, returns False as it is not
    a valid entry. If NA/None or empty string is *not* detected, it
    is assumed that string is a valid piece of data (e.g., from the
    configuration file), and returns True.
    :param string: string to test for valid data entry
    :return: True, if string is a valid entry; otherwise False.
    '''
    r = re.compile(VALID_ENTRY_REGEX, re.I)
    results = list(filter(r.search, string))

    if len(results) > 0:
        return True
    else:
        return False

# read in data from the initialization file into a dict (config_dict)
# THIS WILL THROW ERROR BECAUSE PATH_INIT IS FILE NOT PATH
config_dict = import_config_as_dict(path_to_file=path_init, file_handle=)

config_parser = configparser.RawConfigParser()
config_parser.optionxform = lambda option: option  # use custom config to preserve case of config file
config_parser.read(path_init)
config_dict = {sect: dict(config_parser.items(sect)) for sect in config_parser.sections()}

# define constants
SEQ_FILE_REGEX = r'fast.'
VALID_ENTRY_REGEX = r'[^(^NA$)|^(^None$)]'
MOCK_REGEX = r'^mock'
NTC_REGEX = r'^ntc'
UNDET_REGEX = r'^undet'

DELIMITER = config_dict['delimiter']['delimiter']

# helpful dictionaries
filename_components = {'seq_platform': ['sanger', 'illumina', 'pacbio'],
                       'compartment': ['litter', 'soil', 'spore', 'sporocarp-f', 'sporocarp-a', 'leaf-sp01',
                                       'leaf-sp02', 'seed-sp01', 'seed-sp02', 'root-sp01', 'root-sp02'],
                       'coll_date': 'YYYY-MM',
                       'ecoregion': ['D01', 'D03', 'D05', 'D06', 'D13', 'D14', 'D16', 'D19'],
                       'treatment': ['UC', 'UO', 'UG', 'BC', 'BO', 'BG'],
                       'subplot': ['01', '02', '03', '04', '05', '06', '07', '08', '09']}

# create dict of NEON domains, their associated states, and room to add unique site names from original labels
neon_domains = {'AK': {'D19': []},
                'AZ': {'D14': []},
                'CO': {'D13': []},
                'FL': {'D03': []},
                'KS': {'D06': []},
                'MA': {'D01': []},
                'MN': {'D05': []},
                'OR': {'D16': []}}

# get values from a dictionary; works for nested and not-nested
def get_dict_values(dictionary):
  for val in dictionary.values():
    if isinstance(val, dict):
      yield from get_dict_values(val)
    else:
      yield val

# extracted nested portion of a nested dictionary
def extract_nested(dictionary):
    output_dict = {}

    for val in dictionary.values():
        if isinstance(val, dict):
            output_dict.update(val)
        else:
            return print('Provided dictionary does not appear to be nested.')

    return output_dict

# function that will update the neon_domains dictionary if there is a new value for the domain
def update_domain_dict(domain_dict, site_section='site_name', configuration_dict=config_dict):
    '''
    Updates the NEON domains dictionary with new site values.

    Looks through the site names given by the configuration
    file submitted by the user and checks these values against
    the values in the neon_domains dictionary. If a site name
    provided in the initialization file is not already a value
    in the neon_domains dictionary, it will add it to the
    dictionary to the corresponding domain.
    :param domain_dict: a nested dictionary where the primary
    key is the state abbreviation, the secondary key is the
    NEON domain, and the value of the NEON domain is an empty
    list.
    :param site_section: the name of the section that contains
    information on the site name used in the old file names;
    defaults to 'site_name', which is the name given in the
    original configuration file template, but can be changed
    if altered in the configuration file.
    :param configuration_dict: dictionary read in from the user-
    provided configuration file; by default, it should be called
    config_dict, but able to change if not.
    :return: None; will update the input dictionary.
    '''

    site_dict = configuration_dict[site_section]
    state_abbr = list(site_dict.keys())

    for state in state_abbr:
        # get list of config sites to add to neon domain dict
        config_site_name = site_dict[state]

        domain = list(domain_dict[state].keys())[0]

        domain_dict[state][domain] += config_site_name

    return None

update_domain_dict(neon_domains)

# function that differentiates Illumina versus PacBio raw reads based on start of file name
def determine_sequence_platform(filename):
    '''
    Get sequence platform based on file name format.

    Determines the sequence platform that the sequence
    file originated from based on the format of the raw
    sequence file before renaming. If it is a PacBio
    dataset, it will start with four digits; otherwise,
    it is an Illumina sequence file.
    :param filename: original raw read file name as string
    :return: sequencing platform as string; either 'pacbio'
    or 'illumina'
    '''
    if re.search('^\d{4}', filename):
        seq_platform = 'pacbio'
    else:
        seq_platform = 'illumina'

    return seq_platform

# function that will drop the file extension of a file name
def drop_file_extension(filename, add_path='', add_extension='', keep_path=False):
    '''
    Get base file name of a sequencing file.

    Removes the sequence file extension of the input filename and
    returns a string of the file name with the extension removed.
    Function looks specifically for a sequence file extension such
    as fasta, fastq, fastq.gz, fastx to remove. Options to add a
    new file path or a new file extension to the base name of the
    input file. Helper function to create dummy/test files and when
    renaming sequencing files to updated naming convention.
    :param filename: name of the sequencing file with its original
    file extension; should be some kind of fastx file (O.K. if it is
    a fastq.gz file); an error will be raised if a part of the filename
    that is not the file extension contains 'fast*'.
    :param add_path: file path to add to the beginning of the base file
    name; optional, will default to nothing added.
    :param add_extension: file extension to add to the beginning of the
    base file name; optional, will default to nothing added.
    :param keep_path: whether to keep the path on the file name or remove
    it; default is False, which means the file path is not kept
    :return: base file name of the input file name with original file
    extension removed; may also contain a file path prefix (if provided)
    and/or a new file extension (if provided).
    '''

    # remove the input file path prefix if it is not wanted
    if keep_path:
        pass
    else:  # default
        filename = os.path.basename(filename)

    # ensure the output path has a trailing slash so filepath is created correctly
    if (not add_path.endswith('/')) and (add_path != ''):
        add_path = add_path + '/'

    # ensure the file type has a leading '.' so empty files created correctly
    if (not add_extension.startswith('.')) and (add_extension != ''):
        add_extension = '.' + add_extension

    # split at '.', which should at minimum split away the file extension
    filename_components = filename.split('.')

    # find part (index) of the file name from split component list that matches the sequence file type
    index_to_drop = [c for c, comp in enumerate(filename_components) if re.search(SEQ_FILE_REGEX, comp)]

    if len(index_to_drop) == 1:  # should only be one part of the file name that matches fast*
        # keep only part of original filename that comes before the sequence file type
        components_keep = filename_components[:index_to_drop[0]]
        if len(components_keep) > 1:  # if the filename contains a '.' outside its file extension...
            final_out = add_path + '.'.join(components_keep) + add_extension  # ...join the components back with a '.'
        else:  # if the original filename does not otherwise contain a '.'...
            final_out = add_path + components_keep[0] + add_extension  # ...no joining required
    else:  # raise an error if there's confusion about what is the file extension and what is part of the filename
        print(f"\nERROR: more than one component of the file name {filename} contains text matching 'fast*'. Rename "
              f"this section of the file name and retry.\n")
        sys.exit()

    return final_out

# convert input from config file to regular expression
def convert_config_regex(config_input):
    '''
    Convert simple text with asterisk wildcard to regex.

    Config file accepts wildcards in the common form of an
    asterisk (*) to make it more user-friendly. This function
    will take that input and convert to a regular expression.
    If there are multiple strings in the input for this function,
    it will return a list of regex strings, which can be later
    joined to kept separate.
    :param config_input: string or list of strings where asterisk
    wildcard is used; wildcard not necessary, but if not present, will
    build regex for exact match to input (i.e., not part of other
    word).
    :return: list of regex if input list, single regex string if input
    string.
    '''
    result = []

    for search_str in config_input:  # assemble regex based on config string/list provided

        search_str = search_str.replace('.', '\.')  # if there's a '.' as part of original label, need to escape for regex

        if search_str.startswith('*'):  # if wildcard at start, search for string at end of label
            result.append(search_str.replace('*', '.*') + '$')
        elif search_str.endswith('*'):  # if wildcard at end, search for string at beginning of label
            result.append(
                '^' + search_str.replace('*', '.*'))  # last asterisk, can be anything or nothing after
        else:  # otherwise, search for exact match
            result.append('^' + search_str + '$')

    return result

# get tuple of original filenames; tuple so that no changes can be made to these original file names
original_filenames = tuple([drop_file_extension(filename=file) for file in os.listdir(path_input) if re.search(SEQ_FILE_REGEX, file)])
count_undet = len([f for f in original_filenames if re.search('undetermined', f, re.I)])

# create a dictionary of the original file names and their sizes, as a check after renaming
def file_size_dict(directory, filesize_unit='MB'):
    '''
    Create a dictionary of original file names and their
    corresponding file sizes.

    For each of the original file names, the size of the file
    is calculated and saved into a dictionary. This will serve
    as a test after renaming files that the renamed files are
    corresponding to the correct original files.
    :param directory: the path to the directory containing
    the sequence files to rename; can be absolute or relative to
    the location of this script.
    :param filesize_unit: the memory unit to use for the output,
    default is megabytes (MB); options: bytes (B), kilobytes (KB),
    megabytes (MB), and gigabytes (GB)
    :return: dictionary where keys are names of the original files
    and values are the size of the original files in (bites?)
    '''
    # ensure the output path has a trailing slash so filepath is created correctly
    if not directory.endswith('/'):
        directory = directory + '/'

    # get path of all original sequence files in the input directory
    path_seq_files = [(directory + file) for file in os.listdir(directory) if re.search(SEQ_FILE_REGEX, file)]

    # get the conversion to use for output file size
    conversion_dict = {'B': 0, 'KB': 1024, 'MB': 1024**2, 'GB': 1024**3}  # conversion relative to bytes
    accepted_unit_regex = '|'.join(list(conversion_dict.keys()))
    if re.search(accepted_unit_regex, filesize_unit, re.I):
        conversion = conversion_dict[filesize_unit]
        seq_file_dict = {drop_file_extension(k): {'file_size': ((os.stat(k).st_size)/conversion),
                                                  'unit': filesize_unit.upper()}
                         for k in path_seq_files}
    else:
        return print(f'ERROR: provided file size unit ({filesize_unit}) is not a valid unit. Valid '
                     f'file size units are: {list(conversion_dict.keys())}. Choose one of the accepted units and '
                     f'rerun.')

    return seq_file_dict

original_file_sizes = file_size_dict(directory=path_input)

# confirm that the original filenames in the file size dict are exactly the same as in original_filenames
assert tuple(original_file_sizes.keys()) == original_filenames, 'Original file names do not match the file names ' \
                                                                'in the file size dictionary. The final renamed file ' \
                                                                'check will be inaccurate if not fixed. '

# create dictionary of original label components; helper function in get_old_new_labels()
def get_original_labels(filenames, configuration_dict=config_dict):
    '''
    Get file name label components of original Illumina file names.

    Substitution of old labels with new labels from the updated
    naming convention requires knowledge of all possible old
    labels used. In some cases, there are differences in the
    use of punctuation in label components (i.e., UB-C and UBC)
    so gathering all possible label components is crucial to
    the relabeling process. This function only works with raw
    reads from Illumina sequencing, as PacBio raw reads do not
    contain any relevant information. This function also requires
    the configuration file to contain entries for all component
    positions in the original label, and the delimiter used in the
    original label.
    :param filenames: path of directory containing sequences
    to rename.
    :param configuration_dict: the name of the dictionary that was
    built from the configuration file; default is config_dict, and
    is unlikely to be different unless changed by user.
    :return: a dictionary with keys being the components of the
    new sequence file name and the values being the old components
    used for those categories.
    '''

    # create an empty set for each relevant file name component
    filename_comps = list(filename_components.keys())
    filename_comps.remove('coll_date')  # collection date is handled separately, do not include in this function
    set_list = []  # used at end when created output dictionary
    for comp in filename_comps:
        set_name = f'og_{comp}'  # create name of set based on component
        set_list.append(set_name)  # add name of set to list for use when creating output dict
        exec(f'{set_name} = set()', globals())  # execute creation of set using this set name

    missing_components = []
    # helper function to find label for each component in a file name
    def find_original_components(component_name, filename):
        '''
        Find the original label used for a given component.

        Helper function for get_original_labels. Looks for the label used
        for each component present in the original file name. Skips over
        components that are missing in the original file name.
        :param component_name: the name of the component of the label (str), such as
        seq_platform, compartment, etc.; see keys of filename_components
        dictionary for all options; will raise an error if it is not one
        of the component names in this dictionary.
        :return: none; will add the label for the component to the set of
        file name components created in get_original_labels function.
        '''
        comp_dict = configuration_dict['component_positions']
        component_list = filename.split(DELIMITER)

        comp_options = list(comp_dict.keys())
        if not component_name in comp_options:
            print(f'\nERROR: component provided ({component_name}) not a valid option. '
                  f'Valid names for components are: {comp_options}')
            return False

        # get component from file name, if component present
        comp_value = comp_dict[component_name][0]
        try:  # if an integer is provided for the component's index...
            comp_label = component_list[int(comp_value)]  # index the component list by this integer to get label
        except:  # if component index cannot be converted to an integer (i.e., is None, NA, '', etc.)
            # print(f'\nERROR: the component {component_name} was not found in the original file name.\n')
            missing_components.append(component_name)
            return True  # end function, component not in original filename

        # get the name of the set to add component label to
        og_set = f'og_{component_name}'

        # add component label to the original component set if it follows expected naming convention
        if component_name == 'ecoregion':  # ecoregion is handled differently to check it is added to NEON dictionary
            if comp_label in list(
                    itertools.chain(*list(get_dict_values(neon_domains)))):  # look for it in neon_domain nested values
                exec(f'{og_set}.add(comp_label)')  # add ecoregion to set once confirmed its in NEON domain dict
            else:
                print(
                    f'\nERROR: Ecoregion/site used in original file name {comp_label} not properly added to the NEON domains '
                    f'dictionary. Check that you updated your configuration file before starting, and rerun.')
                return False
        else:
            if component_name == 'compartment':
                search_list = [val for val in itertools.chain(*list(config_dict['compartment_labels'].values())) if
                               is_valid_entry(val)]  # get all possible label values
            elif (component_name == 'treatment') or (component_name == 'subplot'):
                search_list = [val for val in config_dict['accepted_labels'][component_name] if is_valid_entry(val)]
            else:
                print(f'Unrecognized file component name: {component_name}. Accepted component names are: '
                      f'{list(filename_components.keys())}')
                return False

            search_re_list = convert_config_regex(search_list)

            search_regex = '|'.join(search_re_list)  # search for all possible labels for this component

            if re.search(search_regex, comp_label, re.I):
                exec(f'{og_set}.add(comp_label)')
            else:
                print(f"\nERROR: File {filename} cannot be renamed using the automated program because it does not "
                      f"follow the naming convention provided in the configuration file.")
                return False

        return True

    # generator to produce results only if component items are valid
    def search_until_error(component_list, filenames):
        invalid_filenames = []
        for file in filenames:
            for comp in component_list:
                if not find_original_components(comp, file):
                    invalid_filenames.append(file)
                    break
                else:
                    continue

        yield invalid_filenames

    filenames_use = [file for file in filenames if (determine_sequence_platform(file) == 'illumina') and
                     not (re.search('undetermined', file, re.I))]

    invalid_filenames = tuple(search_until_error(filename_comps, filenames_use))[0]


    print(f'\n{len(missing_components)} out of {len(filenames)-count_undet} original file names were missing the '
          f'following file name component(s): {set(missing_components)}. This information will be gathered '
          f'automatically from either the configuration file or the program\'s auto-detect functions.\n')

    # create dictionary where keys are component names and values are all possible values used in original filenames
    output_dict = {k: [] for k in filename_comps}
    for k in output_dict.keys():
        output_set = f'og_{k}'
        exec(f'output_dict[k] += {output_set}')

    return output_dict, invalid_filenames

original_labels, invalid_filenames = get_original_labels(original_filenames)

# create a list of only valid filenames that follow the naming convention from config file
valid_filenames = tuple(set(original_filenames).difference(set(invalid_filenames)))

# remove any punctuation that might be in the old label components; helper function in get_old_new_labels
def remove_punctuation(old_labels, replace='', remove_duplicates=True):
    '''
    Remove punctuation from old labels.

    Before creating the oldnew dictionaries within the function
    get_old_new_labels, the old labels will be screened for
    punctuation, and if found, the punctuation will be removed.
    If there are duplicate labels after punctuation is removed,
    it will also remove duplicate labels. Removing punctuation
    is necessary in order for future renaming steps to work
    correctly.
    :param old_labels: a list of the old labels used in the
    original file names.
    :param replace: option to replace the punctuation with a custom
    string; default is set to remove punctuation and replace with
    nothing.
    :return: list of the old labels with all punctuation removed.
    '''
    s = '\W'  # search regex

    if isinstance(old_labels, list):
        no_punctuation = []
        non_strings = []
        for label in old_labels:
            if isinstance(label, str):
                no_punctuation.append(re.sub(s, replace, label))
            else:
                no_punctuation.append(label)
                non_strings.append(label)
        if len(non_strings) > 0:
            print(f'Some elements of the input list were not strings ({non_strings}), so they were not '
                  f'assessed for punctuation, but remain in the output list.')
        if remove_duplicates:
            return list(set(no_punctuation))
        else:
            return no_punctuation
    elif isinstance(old_labels, str):
        return re.sub(s, replace, old_labels)
    else:
        return print(f'The provided input data type ({type(old_labels)}) is not a valid input.'
                     f'This function only accepts strings or lists of strings as input.')

# create dictionary for each component of old labels and new labels (2 keys each dictionary, 'old_labels', 'new_labels')
def get_old_new_labels(filenames, original_labels):
    '''
    Create nested dictionary of old and new labels by file name component.

    Goes through the original file names to access the original labels
    of the file names to be updated, using the get_original_labels
    function. Then will create a nested dictionary of all file
    name components that are present in the old file names. The primary
    keys of this output dictionary is the name of the component, and then
    each of the components will have two secondary keys: 'old_labels'
    and 'new_labels'. These keys have values of all old labels used for
    that component and all new labels that will be used to replace the
    old ones.
    :param filenames: a list or tuple of the original file names
    that will be renamed with the updated file naming convention.
    :return: a nested dictionary where the primary keys are
    the different file name components and each primary key has the secondary
    key 'old_labels' and 'new_labels'.
    '''

    output_dict = {k: {'old_labels': [], 'new_labels': []} for k in original_labels.keys()}

    for k in original_labels.keys():
        old_labels = remove_punctuation(original_labels[k])
        new_labels = list(filename_components[k])

        if len(old_labels) == len(new_labels):
            output_dict[k]['old_labels'] += remove_punctuation(original_labels[k])
            output_dict[k]['new_labels'] += list(filename_components[k])
        else:
            if k == 'seq_platform':
                output_dict[k]['old_labels'] += old_labels
                output_dict[k]['new_labels'] += set([determine_sequence_platform(file) for file in filenames])
            elif k == 'compartment':
                valid_new_labels = set()

                for lab in config_dict['compartment_labels'].keys():
                    values = config_dict['compartment_labels'][lab]
                    for v in values:
                        if is_valid_entry(v):
                            valid_new_labels.add(lab)
                        else:
                            continue

                output_dict[k]['old_labels'] += old_labels
                output_dict[k]['new_labels'] += valid_new_labels
            elif k == 'ecoregion':
                valid_new_labels = set()
                label_dict = extract_nested(neon_domains)

                for lab in label_dict.keys():
                    values = label_dict[lab]
                    for v in values:
                        if is_valid_entry(v):
                            valid_new_labels.add(lab)
                        else:
                            pass

                output_dict[k]['old_labels'] += old_labels
                output_dict[k]['new_labels'] += valid_new_labels
            elif k == 'subplot':
                old_numbers = '|'.join([re.search(r'([1-9])', ol)[0] for ol in old_labels])
                valid_new_labels = [nl for nl in new_labels if re.search(old_numbers, nl)]

                output_dict[k]['old_labels'] += old_labels
                output_dict[k]['new_labels'] += valid_new_labels
            else:
                print(f'I did not yet write a solution for this scenario, sorry!')

    return output_dict

oldnew_dict = get_old_new_labels(valid_filenames, original_labels)

# before renaming files, create a copy of all files to be renamed and send to zipped subdirectory (for safety)
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
    # ensure the output path has a trailing slash so filepath is created correctly
    if not directory.endswith('/'):
        directory = directory + '/'

    # assemble the path of the copy directory
    copy_path = directory + copy_directory + '/'

    # create the copy directory
    if os.path.isdir(copy_path):  # if a directory already exists with this name...
        if len(os.listdir(copy_path)) > 0:  # if the directory isn't empty...
            ask_to_delete = f'\nA directory with the name {copy_directory} already exists in {directory}. Do you ' \
                            f'want to delete this directory and continue? If you choose no (\'N\') and do not want ' \
                            f'to delete the directory, then the program will exit. (Y/N)\n'
            response = input(ask_to_delete)  # ask if you want to delete it
            if positive_response(response):  # if CLI response yes, delete old directory and create new
                shutil.rmtree(copy_path)
                os.mkdir(copy_path)
            else:  # otherwise, exit program
                sys.exit()
        else:
            shutil.rmtree(copy_path)  # remove existing (empty) directory
            os.mkdir(copy_path)  # create new directory
    else:  # if the directory doesn't exist...
        os.mkdir(copy_path)  # create it

    # get a list of the sequence files to copy to the copy directory
    files_to_copy = [(directory + file) for file in os.listdir(directory) if re.search(SEQ_FILE_REGEX, file)]

    # copy each sequence file to the copy directory
    for file in files_to_copy:
        shutil.copy(file, copy_path)  # to preserve extra metadata, use copy2

    # before archiving copy, get count of sequence files from source and destination
    num_copied_files = len(os.listdir(copy_path))
    num_original_files = len(files_to_copy)

    if compress:
        # zip the copy directory
        shutil.make_archive(copy_path, 'zip', copy_path)
        # remove the unzipped copy directory
        shutil.rmtree(copy_path)

    # confirm all files were copied to the copy directory
    if num_copied_files == num_original_files:
        return print(f'\n{num_copied_files} original sequence files were successfully copied to {copy_directory}.')
    else:
        num_missing = abs(num_copied_files - num_original_files)
        return print(f'ERROR: {num_missing}/{num_original_files} sequence files from the source directory {directory} '
                     f'were not copied to the destination directory {copy_directory}.')

if no_copy:
    pass
else:
    copy_original_files(directory=path_input)

# create a conversion dictionary with pairs of old and new labels for each component
def create_conversion_dict(filenames, old_new_dictionary=oldnew_dict):
    '''
    Create a dictionary with key value pairs of old and new file name components.
    :param filenames: list of the original file names for which to create
    the conversion dictionary
    :param old_new_dictionary: dictionary created using the function get_old_new_labels,
    which creates a dictionary of old labels and new labels for each component of
    the new file name; by default, will look for dictionary called 'oldnew_dict'
    :return: nested dictionary, where primary keys are the file name components,
    nested keys are the original (old) file name labels and values are the
    corresponding labels for the new file naming convention.
    '''

    # helper functions
    def in_original_filename():
        '''
        Checks whether provided component is in the original filename.

        Uses the old_new_dictionary to assess, on a collection-wide basis,
        whether the provided component was a part of the original (old)
        file name. Does not assess on a per-file name basis, but by using
        the old_new_dictionary, it looks whether the original file name
        structure should contain the given component.
        :return: True/False
        '''
        old_labels = component_dict['old_labels']

        if len(old_labels) == 0:  # if no old labels are listed for given component
            return False  # return False; component not in old file names
        else:
            return True  # return True; component is in old file names

    def pair_oldnew_labels(component):

        if component == 'compartment':
            config_section = 'compartment_labels'
        else:
            config_section = 'accepted_labels'

        for l in component_dict['new_labels']:
            if component == 'treatment':  # treatment has two pieces of info collapsed into one string
                burn, habitat = list(l)
                new_label_regex = f'^{burn}.*{habitat}$'
                old_labels = [l for l in list(component_dict['old_labels']) if re.search(new_label_regex, l, re.I)]
                conversion_dict[component][l] = old_labels
            elif component == 'subplot':
                subplot_re = re.compile(f"[{re.search('[1-9]', l)[0]}]")
                old_labels = list(filter(subplot_re.search, component_dict['old_labels']))
                conversion_dict[component][l] = old_labels
            elif component == 'compartment':
                old_label_regex = '|'.join(convert_config_regex(config_dict[config_section][l]))
                old_labels = [l for l in list(component_dict['old_labels']) if re.search(old_label_regex, l, re.I)]
                conversion_dict[component][l] = old_labels
            else:
                print(f'This function was not built to work with {component}.')

        return None

    # create empty dictionary with primary keys for new file name components
    conversion_dict = {k:{} for k in filename_components.keys() if not k == 'coll_date'}

    # go through each component of the new file name and add key value pairs of old new labels
    for comp in conversion_dict.keys():
        component_dict = old_new_dictionary[comp]  # used in both helper functions

        if in_original_filename():  # is there a conversion to be made?
            if comp == 'ecoregion':
                conversion_dict[comp] = extract_nested(neon_domains)
            else:
                pair_oldnew_labels(comp)
        else:  # if no conversion, get info from config file or other means
            if comp == 'seq_platform':
                conversion_dict[comp][determine_sequence_platform(filenames[0])] = ''
            else:
                print('IDK what to do here yet.')

    return conversion_dict

conversion_dict = create_conversion_dict(filenames=valid_filenames)

def find_colldate(config_path):
    coll_dates = re.findall('\d{4}-\d{2}', str(path_init).split('/')[-1])
    if len(coll_dates) == 1:
        return coll_dates[0]
    else:
        print(f'\nMultiple collection dates ({coll_dates}) inferred from the name of the provided '
              f'configuration file (.ini). Please rename the file so that it only contains one '
              f'collection date with the format YYYY-MM (typically 202#-05 for spring collections, 202#-10 '
              f'for fall collections), or specify the collection date for this group of sequences '
              f'with the command line argument \'--coll_date\' following this same date format.')
        return sys.exit()

def rename_invalid_files(filename, conversion_dictionary=conversion_dict):
    if not 'id_index' in globals():
        id_index = 0

    unique_id = filename.split(DELIMITER)[id_index]

    # format known/excepted invalid file name identifiers
    if re.search(MOCK_REGEX, unique_id, re.I):
        unique_id_lower = unique_id.lower()
        mock_break = re.search(MOCK_REGEX, unique_id_lower, re.I).span()[1]
        unique_id = unique_id_lower[:mock_break] + '-' + unique_id_lower[mock_break:]
    elif re.search(UNDET_REGEX, unique_id, re.I):
        unique_id = unique_id.lower()
    elif re.search(NTC_REGEX, unique_id, re.I):
        ntc_break = re.search(NTC_REGEX, unique_id, re.I).span()[1]
        ntc_num = unique_id[ntc_break:]
        if int(ntc_num) < 10:
            unique_id = unique_id[:ntc_break] + '-0' + unique_id[ntc_break:]
        else:
            unique_id = unique_id[:ntc_break] + '-' + unique_id[ntc_break:]
    else:
        pass

    compartment = '-'.join(list(conversion_dictionary['compartment'].keys()))
    new_name = {k: v for k, v in zip(['platform', 'compartment', 'coll_date', 'description'],
                                     ['illumina', compartment, '', unique_id])}
    if colldate_group is None:
        new_name['coll_date'] = find_colldate(config_path=path_init)
    else:
        new_name['coll_date'] = colldate_group

    return new_name

# main function that creates new file name
def rename_sequence_files(filenames, conversion_dictionary=conversion_dict):

    filename_change = {}

    positions = config_dict['component_positions']

    for old_filename in filenames:
        # haven't figured out pacbio yet so just filter everything that isn't illumina
        if (determine_sequence_platform(old_filename) == 'illumina'):
            components = old_filename.split(DELIMITER)
            if (re.search(UNDET_REGEX, old_filename, re.I)) or (old_filename in invalid_filenames):
                new_name = rename_invalid_files(old_filename)
            else:
                new_name = {k:'' for k in filename_components.keys()}

                for c in new_name.keys():  # search for each new component in old file name
                    try:
                        if is_valid_entry(positions[c]):
                            old_value = remove_punctuation(components[int(positions[c][0])])
                            new_value = [k for k in conversion_dict[c].keys() if old_value in conversion_dict[c][k]]
                            new_name[c] = new_value[0]
                    except:
                        if c == 'coll_date':
                            site = components[int(positions['ecoregion'][0])]
                            domain = \
                                [k for k in extract_nested(neon_domains).keys() if site in extract_nested(neon_domains)[k]][0]
                            state = [k for k in neon_domains.keys() if domain in list(neon_domains[k].keys())][0]
                            coll_dates = config_dict['collection_date'][state]
                            if len(coll_dates) > 1:
                                coll_date_index = config_dict['site_name'][state].index(site)
                                coll_date = coll_dates[coll_date_index]
                            else:
                                coll_date = config_dict['collection_date'][state][0]
                            new_name[c] = coll_date
                    else:
                        if c == 'seq_platform':
                            new_name[c] = list(conversion_dict[c].keys())[0]

            read_id_regex = re.compile('^R1$|^R2$', re.I)
            read_id = list(filter(read_id_regex.search, components))
            new_name['read_id'] = read_id[0]

            # join together new file name, add to file name dictionary
            new_filename = '_'.join(list(get_dict_values(new_name)))
            filename_change[old_filename] = new_filename

    # check for duplicates in the new file names
    def check_duplicate_filenames(filename_change):
        def find_duplicates(dictionary):
            checked = set()
            duplicate_names = {}

            for pair in dictionary.items():
                if pair[1] in checked:
                    duplicate_names.update({pair})
                else:
                    checked.add(pair[1])

            # ensures that all duplicate file names, even first encountered instance, included in duplicates
            for pair in dictionary.items():
                if (pair[1] in get_dict_values(duplicate_names)) and not (pair[0] in duplicate_names.keys()):
                    duplicate_names.update({pair})
                else:
                    continue

            return duplicate_names

        wc_compartment = [l[1] for l in config_dict['compartment_labels'].items() if re.search('\*', l[1][0])][0]
        compart_regex = convert_config_regex(wc_compartment)[0].replace('.*',
                                                                        '\w+')  # only find labels with S + another letter
        compart_pos = int(positions['compartment'][0])

        def get_original_basename(filename, include_added_value=False, return_max=False):
            '''
            Extracts the base sample name from original file name.

            Separates the original filename input by the delimiter designated in the configuration
            file. Uses the positions sections of the configuration file to determine the maximum
            index of the file name components. Returns the sample's base name as a string. This is
            helpful for removing any additional information that the sequencer might add onto a
            sample. Also required for fixing issues of duplicate new names.
            :param filename: original filename as a string
            :param include_added_value: if to include an additional value in the original filename
            that is not accounted for in the configuration file; this is helpful for renaming
            duplicates, as new names that are duplicates are often because there is additional information
            that is not accounted for in the standard file name composition of the old file names;
            default False, will only return base file name as outlined by configuration file.
            :param return_max: whether to return the maximum index used to find the base file name;
            default is False, but may be helpful to return at times. If set to True, it will be the
            second item returned by the function.
            :return: the file's base name
            '''
            max_index = max([int(val[0]) for val in get_dict_values(positions) if is_valid_entry(val)])
            if include_added_value:
                max_index += 2
            else:
                max_index += 1  # add one to include up to max_index in slice

            basename = DELIMITER.join(filename.split(DELIMITER)[:max_index])

            if return_max is True:
                return basename, max_index
            else:
                return basename

        def find_similar_names(filename, search_dict, include_added_value=False, return_max=False):
            basename = get_original_basename(filename, include_added_value=include_added_value, return_max=return_max)
            base = re.compile(f'^({basename})', re.I)
            similar = list(filter(base.search, search_dict.keys()))
            return [s for s in similar if s != filename]

        def add_subplot_letter(duplicate_new_name):
            subplot_pos = list(filename_components.keys()).index('subplot')
            subplot = duplicate_new_name.split(DELIMITER)[subplot_pos]
            name = duplicate_new_name

            ascii_start = ord('A')  # start with A
            while name in get_dict_values(duplicate_names):  # while it remains a duplicate
                name = re.sub(f'_{subplot}_', f'_{subplot}{chr(ascii_start)}_', duplicate_new_name)
                ascii_start += 1  # move through alphabet
                if ascii_start == ord('Z'):  # stop if end of alphabet reached
                    return print('Reached end of alphabet and ran out of letters; probably issues elsewhere in code.')
            else:
                return name

        def are_read_pairs(file01, file02):

            def check_input_type(file):
                if isinstance(file, str):
                    return file
                elif isinstance(file, list) and len(file) == 1:
                    return file[0]
                else:
                    return print('Function only handles strings or lists of single string value.')

            file01 = check_input_type(file01)
            file02 = check_input_type(file02)

            parts1 = file01.split(DELIMITER)
            parts2 = file02.split(DELIMITER)

            match = 0
            for p1, p2 in zip(parts1, parts2):  # confirm they're exactly same at each location except that R1 == R2
                if (p1 == p2) or ((p1 in ['R1', 'R2']) and (p2 in ['R1', 'R2']) and (p1 != p2)):
                    match += 1
                    continue
                else:
                    match += 0

            if match == len(parts2):
                return True
            else:
                return False

        def read_pair_name_match(file01, file02):
            subplot_index = list(filename_components.keys()).index('subplot')
            f1 = DELIMITER.join(file01.split(DELIMITER)[:subplot_index + 1])
            f2 = DELIMITER.join(file02.split(DELIMITER)[:subplot_index + 1])
            if f1 == f2:
                return True
            else:
                return print(f'File renaming issue. Paired read files do not match:\n'
                             f'\t{file01}\n'
                             f'\t{file02}')

        duplicate_names = find_duplicates(filename_change)
        encountered = set()
        for pair in duplicate_names.items():
            old = pair[0]
            new = pair[1]
            compartment_old = old.split(DELIMITER)[compart_pos]

            if old in encountered:  # because some names changed pairwise, don't try renaming same file twice
                continue

            if re.search(compart_regex, compartment_old, re.I):
                similar = find_similar_names(old, duplicate_names)
            else:
                similar = find_similar_names(old, duplicate_names, include_added_value=True)

            if are_read_pairs(old, similar):
                fixed_name = add_subplot_letter(new)
                fixed_paired_name = add_subplot_letter(duplicate_names[similar[0]])
                if read_pair_name_match(fixed_name, fixed_paired_name):
                    duplicate_names[old] = fixed_name
                    duplicate_names[similar[0]] = fixed_paired_name

            encountered.add(similar[0])
            encountered.add(old)

        if len(list(find_duplicates(duplicate_names).keys())) == 0:
            for old in list(duplicate_names.keys()):
                filename_change[old] = duplicate_names[old]
            return print(
                f'\n{len(list(duplicate_names.keys()))} files were detected to have duplicate new '
                f'names, and have been renamed so that they are unique. A letter has been added to the subplot.')
        else:
            print(f'ERROR: Did not solve all instances of duplicate new name values. Need to check manually.\n'
                  f'duplicate sets:\n\t {find_duplicates(duplicate_names)}')
            return sys.exit()

    check_duplicate_filenames(filename_change)

    renamed_file_count = 0
    for i in range(len(list(filename_change.items()))):  # go through each old-new name pair and change file name
        renamed_file_count += 1
        old, new = list(filename_change.items())[i]

        original_filename = Path((path_input + old + file_ext))  # need file extension?
        updated_filename = Path((path_input + new + file_ext))  # don't need file extension?

        # rename file
        original_filename.replace(updated_filename)

    if table_out:
        combo_compart = '-'.join(list(conversion_dictionary['compartment'].keys()))
        out_table_path = f'{path_input}illumina_{find_colldate(path_init)}_{combo_compart}_file-rename-conversion.csv'

        if os.path.isfile(out_table_path):
            current_df = pd.read_csv(out_table_path)
            added_df = pd.DataFrame.from_dict(filename_change, orient='index').reset_index()
            added_df.columns = ['old_name', 'new_name']
            output_df = pd.concat([current_df, added_df])
        else:
            output_df = pd.DataFrame.from_dict(filename_change, orient='index').reset_index()
            output_df.columns = ['old_name', 'new_name']

        output_df.to_csv(out_table_path, index=False)

    return print(f'\n{renamed_file_count} files were successfully renamed.\n')

rename_sequence_files(valid_filenames)

ntc_example = [f for f in invalid_filenames if re.search(NTC_REGEX, f, re.I)]
mock_example = [f for f in invalid_filenames if re.search(MOCK_REGEX, f, re.I)]
fix_error_files = input(f'{len(invalid_filenames)} of {len(original_filenames)} files in this directory '
                        f'were not renamed because they did not follow the naming convention of the old '
                        f'file names provided in the configuration file.\n'
                        f'If these sequences are related to the CliMush project, they can still be renamed '
                        f'using the first component of the old file name for these files (e.g., '
                        f'{ntc_example[0]}, {mock_example[0]}).\n'
                        f'Do you wish to continue with renaming these files this way? [Y/N]\n')

def move_unnammed(destination, files):
    print(f'Moving all remaining files that have not been renamed to {dest_path} in this current directory.\n')

    if os.path.isdir(dest_path):  # if a directory already exists with this name...
        if len(os.listdir(dest_path)) > 0:  # if the directory isn't empty...
            ask_to_delete = f'\nA directory with the name {dest_path} already exists in {path_input}. Do you ' \
                            f'want to delete this directory and continue? If you choose no (\'N\') and do not want ' \
                            f'to delete the directory, then the program will exit. (Y/N)\n'
            response = input(ask_to_delete)  # ask if you want to delete it
            if positive_response(response):  # if CLI response yes, delete old directory and create new
                shutil.rmtree(dest_path)
                os.mkdir(dest_path)
            else:  # otherwise, exit program
                sys.exit()
        else:
            shutil.rmtree(dest_path)  # remove existing (empty) directory
            os.mkdir(dest_path)  # create new directory
    else:  # if the directory doesn't exist...
        os.mkdir(dest_path)  # create it

    # copy each sequence file to the copy directory
    for file in files:
        shutil.move((path_input+file+file_ext), dest_path)
    sys.exit()

dest_path = f'{path_input}non-climush_sequences'
if positive_response(fix_error_files):
    able_to_rename = [f for f in invalid_filenames if re.search(f'{MOCK_REGEX}|{NTC_REGEX}|{UNDET_REGEX}', f, re.I)]
    cannot_rename = list(set(invalid_filenames).difference(able_to_rename))
    rename_sequence_files(able_to_rename)
    if len(cannot_rename) == 0:
        print('All invalid file names were successfully renamed.\n')
    else:
        move_unnammed(dest_path, cannot_rename)
else:
    move_unnammed(dest_path, invalid_filenames)

# try_other = input(f'If you do not want to continue with renaming these files, you can specify a different '
#                   f'index to get the identifying component from the invalid filename. If you\'d like to '
#                   f'provide an index, do so now (zero-indexing). If not, enter \'N\' to exit. This will '
#                   f'move all invalid file names to a sub directory called \'non-climush_sequences\'.\n')
# if isinstance(try_other, int):
#     id_index = try_other
#     rename_sequence_files(invalid_filenames)
# else:
#     move_unnammed(dest_path, invalid_filenames)
