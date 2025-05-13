import argparse, pathlib
from pathlib import Path
from climush.bioinfo import quality_filter
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
    description='Read merging and quality filtering using VSEARCH.',
    epilog='This script is part of the CliMush bioinformatics pipeline.',
)

## FILE PATHS ##

# REQUIRED; input file path to the parent directory containing the sequence files to quality filter
parser.add_argument(
    '-i', '--input',
    default=None,
    required=True,
    type=pathlib.PosixPath,
    action='store',
    help='The path to the directory containing the paired-end or single-end sequence '
         'files to merge and/or quality filter.',
)

# REQUIRED; output file path to the parent directory that the final quality filtered sequences will be written to
parser.add_argument(
    '-o', '--output',
    default=None,
    required=True,
    type=pathlib.PosixPath,
    action='store',
    help='The path to the directory in which the merged reads and/or quality filtered output '
         'directories will be written to.',
)


## READ MERGING FOR PAIRED-END SEQUENCES ##

# if flag is used, the Illumina paired-end sequences in the input directory will be merge prior to quality filtering
parser.add_argument(
    '--merge-reads',
    # default=settings['quality_filtering']['merge_reads'],
    required=False,
    action='store_true',
    help='If flag is included, then the sequence files in the input directory will be merged prior to quality '
         'filtering, if there are Illumina sequences in the input directory. Otherwise, if the sequence files '
         'are not Illumina sequences, no merging will occur and quality filtering will continue without merging.',
)

## INPUT SEQUENCE QUALITY SCORES ##

# accepted minimum quality score of input sequence files
parser.add_argument(
    '--min-qscore',
    default=settings['quality_filtering']['min_qscore'],
    required=False,
    type=int,
    action='store',
    help='The minimum quality score accepted in the input sequence files. This setting does '
         'not affect the minimum quality score of the output quality filtered sequences, '
         'which are filtered based on maximum expected error (--max_error). '
)

# accepted maximum quality score of input sequence files
parser.add_argument(
    '--max-qscore',
    default=settings['quality_filtering']['max_qscore'],
    required=False,
    type=int,
    action='store',
    help='The maximum quality score accepted in the input sequence files. This setting does '
         'not affect the maximum quality score of the output quality filtered sequences, '
         'which are filtered based on maximum expected error (--max_error). '
)


## QUALITY FILTERING - SEQUENCE LENGTH ##

# if option provided, uses this minimum read length filter instead of configuration settings
parser.add_argument(
    '--min-length',
    default=settings['quality_filtering']['min_len'],
    required=False,
    type=int,
    action='store',
    help='The minimum sequence length required for sequences to be written to the output quality '
         'filtered sequences.',
)

# if option provided, uses this maximum read length filter instead of configuration settings
parser.add_argument(
    '--max-length',
    default=settings['quality_filtering']['max_len'],
    required=False,
    type=int,
    action='store',
    help='The maximum sequence length required for sequences to be written to the output quality '
         'filtered sequences.',
)

## QUALITY FILTERING - EXPECTED ERROR ##

# if option provided, uses this maximum expected error filter instead of configuration settings
parser.add_argument(
    '--max-error',
    default=settings['quality_filtering']['max_error'],
    required=False,
    type=float,
    action='store',
    help='The maximum expected errors, calculated as the sum of the error probability (see Edgar & Flybjerg 2007). '
         'A small value will result in higher quality reads and a low value will result in lower quality reads.',
)


## OUTPUT FILES ##

parser.add_argument(
    '--keep',
    action='store_true',
    help='Flag that, when used, will result in all discarded reads from merging and/or quality '
         'filtering being written out to the file system, rather than being completely '
         'discarded.',
)

parser.add_argument(
    '--log',
    action='store_true',
    help='Flag that, when used, will write any optional output from vsearch merging or quality '
         'filtering out to a .log file.',
)

parser.add_argument(
    '--no-compress',
    action='store_false',
    help='Flag that, when used, will prevent the output files from merging and/or quality '
         'filtering to be compressed to the gzip compression format.',
)


## PARSE OPTIONS INTO DICTIONARY ##

args = vars(parser.parse_args())

#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# check if there are illumina sequences to process in the input directory
is_input, illumina_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

# if there are illumina files detected in the input directory...
if is_input:

    # use the quality_filter() function from bioinfo.py to merge (if merge=True) and quality filter input sequences
    illumina_output_files = quality_filter(
        input_files=illumina_files,
        output_dir=args['output'],
        platform=platform,
        reference_dir=ref_dir,
        min_qscore=args['min_qscore'][platform],
        max_qscore=args['max_qscore'][platform],
        min_len=args['min_length'][platform],
        max_len=args['max_length'][platform],
        max_error=args['max_error'][platform],
        merge=args['merge_reads'],
        keep_removed_seqs=args['keep'],
        keep_log=args['log'],
        compress_output=args['no_compress'],
    )

# if there are no illumina files detected in the input directory, continue without doing anything
else:
    pass


#####################
# SANGER ############
#####################
platform = 'sanger'

# check if there are sanger sequences to process in the input directory
is_input, sanger_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

# if there are sanger files detected in the input directory...
if is_input:

    # use the quality_filter() function from bioinfo.py to quality filter input sequences (no merge for sanger)
    sanger_output_files = quality_filter(
        input_files=sanger_files,
        output_dir=args['output'],
        platform=platform,
        reference_dir=ref_dir,
        min_qscore=args['min_qscore'][platform],
        max_qscore=args['max_qscore'][platform],
        min_len=args['min_length'][platform],
        max_len=args['max_length'][platform],
        max_error=args['max_error'][platform],
        merge=False,
        keep_removed_seqs=args['keep'],
        keep_log=args['log'],
        compress_output=args['no_compress'],
    )

# if there are no sanger files detected in the input directory, continue without doing anything
else:
    pass


#####################
# PACBIO ############
#####################
platform = 'pacbio'

# when all are prefiltered, continue to next
continue_to_next(__file__, settings)