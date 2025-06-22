import argparse, pathlib
from pathlib import Path
from climush.utilities import get_settings, check_for_input, continue_to_next, mkdir_exist_ok
from climush.bioinfo import demultiplex

## IMPORT PIPELINE CONFIGURATION #######################################################################################

# create a reference directory path for file-finding functions
ref_dir=Path(__file__).parent

# import the settings for the bioinformatics configuration
settings = get_settings(ref_dir)

########################################################################################################################


## COMMAND LINE ARGUMENTS ##############################################################################################

## INSTANTIATE PARSER ##

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Demultiplex PacBio reads using PacBio\'s Lima demultiplexing '
                                             'tool.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

## FILE PATHS ##

# input directory containing the files to demultiplex
parser.add_argument('-i', '--input',
                    type=pathlib.PosixPath,
                    help='The path to a directory containing the sequencing files to demultiplex; '
                         'must be an absolute path.')

# output directory to send the demultiplexed reads to
parser.add_argument('-o', '--output',
                    type=pathlib.PosixPath,
                    help='The path to write out the demultiplexed sequence files to; must be an absolute path.')


## BARCODE MATCH SETTINGS ##

# if option provided, will use this minimum score instead of configuration settings
parser.add_argument('--min-score', nargs='?',
                    default=settings['demultiplex']['min_precision'],
                    help='The minimum precision of detecting barcodes in the multiplexed reads. A minimum score '
                         'of 80 yields a precision of >99.99%. The higher the score, the lower the level of'
                         'contaminant, but this will also lead to a lower post-demultiplexing read yield. The '
                         'default minimum precision score is 93, as recommended by Tedersoo et al. 2021.')

## STANDARD OUTPUT ##

# if option provided, will print out all warnings
parser.add_argument('--verbose',
                    action='store_true',
                    default=settings['automate']['verbose'],
                    help='Whether to print out all errors and warnings in the terminal; if True, will print all '
                         'errors and warnings; if False, will only print out fatal errors.')

## PARSE OPTIONS INTO DICTIONARY ##

args = vars(parser.parse_args())

########################################################################################################################

# I don't have a way to accept a non-default minimum precision score yet (--min-score) so print warning if provided one
if not isinstance(args['min_score'], dict):
    default_precision = settings['demultiplex']['min_precision']
    print(f'WARNING. This version of {Path(__file__).name} is not yet configured to accept values for '
          f'--min-score that differ from the default values from the configuration file. Using the default '
          f'values of:\n'
          f'   {default_precision.items()}')

# set variable that will only switch to True if files requiring demultiplexing are detected
files_demuxed = False

#####################
# ILLUMINA ##########
#####################
# platform = 'illumina'
#
# # check if there are Illumina reads that need to be demultiplexed
# is_input, illumina_files = = check_for_input(
#     args['input'],
#     config_dict=settings,
#     file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform]
# )
#
# if is_input:
#     msg = f'WARNING. Currently, there is no script that can demultiplex Illumina reads. You can continue with the ' \
#           f'pipeline, but these {len(illumina_files)} multiplexed Illumina sequencing files will be ignored. Do you ' \
#           f'wish to continue without these sequences?'
#     prompt_yes_no_quit(message = msg)
# else:
#     pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# check if there are PacBio reads that need to be demultiplexed


# check for multiplexed sequence files in the input path using info from the config file

# create a regex that will search for any of the 4-digit GC3F sequence run codes
gc3f_codes = '|'.join([ ('^' + num_code) for num_code in settings['demultiplex']['multiplex'].keys()])

# use this combined regex to search of the 4-digit GC3F sequence run code at the start of a file name when
#  searching for input in the input directory
is_input, pacbio_files = check_for_input(
    args['input'],
    config_dict=settings,
    file_identifier=gc3f_codes,
)

if is_input:
    files_demuxed = True
    print(f'{len(pacbio_files)} PacBio sequencing files were detected that require demultiplexing...')

    # create output directory for demultiplexed samples
    mkdir_exist_ok(new_dir = args['output'])

    # demultiplex input files
    demultiplex(
        output_dir=args['output'],
        reference_dir=Path(__file__),
        verbose=args['verbose'],
        multiplexed_files=pacbio_files,
    )

else:
    pass

#####################
# SANGER ############
#####################
# platform = 'sanger'
#
# # check if there are Sanger reads that need to be demultiplexed
# is_input, sanger_files = = check_for_input(
#     args['input'],
#     config_dict=settings,
#     file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform]
# )
#
# if is_input:
#     msg = f'WARNING. Currently, there is no script that can demultiplex Sanger reads. You can continue with the ' \
#           f'pipeline, but these {len(sanger_files)} multiplexed Sanger sequencing files will be ignored. Do you ' \
#           f'wish to continue without these sequences?'
#     prompt_yes_no_quit(message = msg)
# else:
#     pass

#########################

# when all are demultiplexed (if possible), continue to the next script
if not files_demuxed:  # print if no demultiplexing was carried out
    print(f'No sequencing files were detected in the {args["input"]} that require demultiplexing.\n')
continue_to_next(__file__, settings)