import argparse, pathlib
from pathlib import Path
from climush.bioinfo import dereplicate
from climush.utilities import get_settings, check_for_input, continue_to_next

## IMPORT PIPELINE CONFIGURATION #######################################################################################

# create a reference directory path for file-finding functions
ref_dir=Path(__file__).parent

# import the settings for the bioinformatics configuration
settings = get_settings(ref_dir)

########################################################################################################################


## COMMAND LINE ARGUMENTS ##############################################################################################

## INSTANTIATE PARSER ##

parser = argparse.ArgumentParser(
    prog=Path(__file__).stem,
    description='Dereplicate full-length reads.',
    epilog='This script is part of the CliMush bioinformatics pipeline.',
)


## FILE PATHS ##

# input directory containing the files to dereplicate
parser.add_argument(
    '-i', '--input',
    required=True,
    type=pathlib.PosixPath,
    help='The path to a directory containing sequencing files to dereplicate. If nothing provided, '
         'will default to the location that is expected in the Docker container\'s native file '
         'structure, detailed in pipeline/mapping.py.',
)

# output directory for the dereplicated files
parser.add_argument(
    '-o', '--output',
    required=True,
    type=pathlib.PosixPath,
    help='The main output directory in which to create the bioinformatics run subfolder for '
         'dereplicated sequences.')

## DEREPLICATION SETTINGS ##

# minimum number of identical reads needed to pass derepliations
parser.add_argument(
    '--min-count',
    type=int,
    default=settings['dereplicate']['min_count']['derep01'],
    help='The minimum number of identical reads that is required to be written to the dereplicated output file.',
)

## OUTPUT FILE OPTIONS ##

parser.add_argument(
    '--log',
    action='store_true',
    help='Flag that, when used, will write any optional output from vsearch dereplication to a .log file.',
)

## PARSE OPTIONS INTO DICTIONARY ##

args = vars(parser.parse_args())

########################################################################################################################


#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# check that there are Illumina reads to dereplicate
is_input, illumina_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

if is_input:

    print(f'Dereplicating {len(illumina_files)} {platform} reads...\n')

    dereplicate(
        input_files=illumina_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        min_count=args['min_count'][platform],
        derep_step=1,
        keep_log=args['log'],
    )

else:
    pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# check that there are Illumina reads to dereplicate
is_input, pacbio_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

if is_input:

    print(f'Dereplicating {len(pacbio_files)} {platform} reads...\n')

    dereplicate(
        input_files=pacbio_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        min_count=args['min_count'][platform],
        derep_step=1,
        keep_log=args['log'],
    )

else:
    pass

#####################
# SANGER ############
#####################
platform = 'sanger'

# check that there are Illumina reads to dereplicate
is_input, sanger_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

if is_input:

    print(f'Dereplicating {len(sanger_files)} {platform} reads...\n')

    dereplicate(
        input_files=sanger_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        min_count=args['min_count'][platform],
        derep_step=1,
        keep_log=args['log'],
    )

else:
    pass


# when all are seqs are dereplicated, continue to next
continue_to_next(__file__, settings)