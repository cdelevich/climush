import argparse, pathlib
from pathlib import Path
from climush.bioinfo import filter_out_phix, prefilter_fastx
from climush.utilities import get_settings, check_for_input, continue_to_next

## IMPORT PIPELINE CONFIGURATION #######################################################################################

# create a reference directory path for file-finding functions
ref_dir=Path(__file__).parent.parent

# import the settings for the bioinformatics configuration
settings = get_settings(ref_dir)

########################################################################################################################


## COMMAND LINE ARGUMENTS ##############################################################################################

## INSTANTIATE PARSER ##

parser = argparse.ArgumentParser(
    prog=Path(__file__).stem,
    description='Remove PhiX and reads with ambiguous bases',
    epilog='This script is part of the CliMush bioinformatics pipeline.',
)

## FILE PATHS ##

# input directory
parser.add_argument(
    '-i', '--input',
    type=pathlib.PosixPath,
    help='The path to the input sequencing files that require prefiltering; path must '
         'be absolute.'
)

# output directory
parser.add_argument(
    '-o', '--output',
    type=pathlib.PosixPath,
    help='A path to write the prefiltered sequence files to; path must be absolute.',
)

## PHIX FILTER OPTIONS ##

# kmers for filtering out phix with bbduk
parser.add_argument(
    '-k', '--kmer',
    required=False,
    default=settings['quality_filtering']['prefilter']['phix_kmer'],
    type=int,
    help='The number of kmers to use in the search for PhiX spike-in reads for filtering, using bbduk. [max=31]',
)

# hamming distance for phix filtering with bbduk
parser.add_argument(
    '--hdist',
    required=False,
    default=settings['quality_filtering']['prefilter']['phix_hammingdist'],
    type=int,
    help='The hamming distance to use in the search for PhiX spike-in reads for filtering, using bbduk.',
)

# the available memory that bbduk can use for its phix filter
parser.add_argument(
    '--max_mem',
    required=False,
    default=None,
    type=int,
    help='The maximum amount of memory that the PhiX filtering function, bbduk, can use for filtering out '
         'PhiX spike-in reads. If a value is not provided, the maximum available system memory is detected and used. '
         'Accepted units are bytes, megabytes, or gigabytes and the unit should not be included in the value. '
         'All values, regardless of unit, must be rounded to the nearest whole number.',
)

## PREFILTER OPTIONS ##

# maximum number of ambiguous base calls per read
parser.add_argument(
    '--maxn',
    required=False,
    default=settings['quality_filtering']['prefilter']['maxn'],
    type=int,
    help='The maximum number of allowed ambiguous base calls in a read in order for it to pass prefiltering.',
)

# maximum read quality score of illumina sequences
parser.add_argument(
    '--qscore_max',
    required=False,
    default=settings['quality_filtering']['max_qscore']['illumina'],
    type=int,
    help='The maximum allowed read quality score that cannot be surpassed in order for Illumina reads to '
         'pass prefiltering.'
)

## FUNCTION OUTPUT OPTIONS ##

# whether to keep the sequences that were filtered out due to them being PhiX or having ambiguous base calls
parser.add_argument(
    '--keep',
    action='store_true',
    help='Flag that, when used, will retain all of the sequences that are filtered out of the input sequences '
         'in a separate discard directory. This will produce two additional directories, one for PhiX spike-in reads '
         'that were filtered out and one for reads with ambiguous bases that were filtered out.',
)

# whether to produce a log file from bbduk
parser.add_argument(
    '--log',
    action='store_true',
    help='Flag that when used will write a log summary as produced by bbduk that summarizes the search '
         'and filtering out of PhiX spike-in reads.'
)

## PARSE OPTIONS INTO DICTIONARY ##

args = vars(parser.parse_args())

########################################################################################################################


## FIND FILES + PREFILTER READS ########################################################################################


#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# check if there are ITS1 illumina reads to pre-filter
is_input, illumina_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

# if ITS1 illumina sequences were located, prefilter and remove PhiX from these sequences
if is_input:

    # detect and remove PhiX spike-in sequences from input ITS1 illumina sequence files
    illumina_nophix_files = filter_out_phix(
        input_files=illumina_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        kmer=args['kmer'],
        hdist=args['hdist'],
        keep_removed_seqs=args['keep'],
        keep_log=args['log'],
        bbduk_mem=args['max_mem'],
    )

    # detect and remove sequences with ambiguous base calls exceeding value set by maxn (default=0) from ITS1
    #   illumina sequence files that have previously been filtered to remove PhiX spike-in
    illumina_noambig_path = prefilter_fastx(
        input_files=illumina_nophix_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        maxn=args['maxn'],
        qmax=args['qscore_max'],
        keep_removed_seqs=args['keep'],
        keep_log=args['log'],
    )

# if no ITS1 illumina sequence files were located, do nothing
else:
    pass

#####################
# SANGER ############
#####################
platform = 'sanger'

# check if there are 18S illumina sequences to pre-filter
is_input, sanger_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

# if 18S illumina sequences were located, prefilter and remove PhiX from these sequences
if is_input:

    # detect and remove PhiX spike-in sequences from input 18S illumina sequence files
    sanger_nophix_files = filter_out_phix(
        input_files=sanger_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        kmer=args['kmer'],
        hdist=args['hdist'],
        keep_removed_seqs=args['keep'],
        keep_log=args['log'],
        bbduk_mem=args['max_mem'],
    )

    # detect and remove sequences with ambiguous base calls exceeding value set by maxn (default=0) from 18S
    #   illumina sequence files that have previously been filtered to remove PhiX spike-in
    sanger_noambig_path = prefilter_fastx(
        input_files=sanger_nophix_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        maxn=args['maxn'],
        qmax=args['qscore_max'],
        keep_removed_seqs=args['keep'],
        keep_log=args['log'],
    )

# if no 18S illumina sequence files were located, do nothing
else:
    pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# no pre-filtering for PacBio reads

########################################################################################################################

# when all are prefiltered, continue to next
continue_to_next(__file__, settings)


