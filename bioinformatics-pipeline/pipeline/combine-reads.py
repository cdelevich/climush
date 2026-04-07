import argparse, pathlib
from pathlib import Path

from climush.bioinfo import combine_reads
from climush.utilities import get_settings, check_for_input, continue_to_next

# set a location to start looking for pipeline configuration file
ref_dir = Path(__file__).parent

settings = get_settings(ref_dir)
run_name = settings['run_details']['run_name']

# set up command line options
parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Rename read headers and combine reads.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# input directory containing the files to check for chimeras; no default, since depends on platform
parser.add_argument('-i', '--input',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The path to a directory containing the sequence files with reads that will be combined.')

parser.add_argument('-o', '--output',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The path to the directory to which the output combined sequence file will be written.')

parser.add_argument(
    '--rename',
    action='store_true',
    help='Boolean flag that, when used, with update the read headers before combining so that they have the '
         'sample ID and a unique read ID in the header, along with the size information from dereplication if '
         'present.'
)

# parse command line options and defaults into a dictionary
args = vars(parser.parse_args())

#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# check if there are ITS1 sequences in the input directory
is_input, illumina_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
    )

if is_input:
    combine_reads(
        input_dir=illumina_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        platform=platform,
        rename_headers=args['rename'],
    )
else:
    pass

#####################
# SANGER ############
#####################
platform = 'sanger'

# check if there are 18S sequences in the input directory
is_input, sanger_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
    )

if is_input:
    combine_reads(
        input_dir=sanger_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        platform=platform,
        rename_headers=args['rename'],
    )
else:
    pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# check if there are pacbio sequences in the input directory
is_input, pacbio_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
    file_ext=None,
    )

if is_input:
    combine_reads(
        input_dir=pacbio_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        platform=platform,
        rename_headers=args['rename'],
    )
else:
    pass


# continue to next
continue_to_next(__file__, settings)