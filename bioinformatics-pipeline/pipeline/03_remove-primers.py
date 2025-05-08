import argparse, pathlib
from pathlib import Path
from climush.bioinfo import remove_primers, confirm_no_primers
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
    description='Remove primers using cutadapt.',
    epilog='This script is part of the CliMush bioinformatics pipeline.',
)

## FILE PATHS ##

# REQUIRED; input file path to the parent directory containing the sequence files to trim primers from
parser.add_argument(
    '-i', '--input',
    required=True,
    type=pathlib.PosixPath,
    help='Path to the sequence files that require primer removal; path must be absolute.',
)

# REQUIRED; output file path to the directory that the trimmed sequence output should be written to
parser.add_argument(
    '-o', '--output',
    required=True,
    type=pathlib.PosixPath,
    help='File path to the directory to which the trimmed sequence files and any other summary files '
         'will be written.',
)

## ONLY CHECK FOR PRIMERS, NO TRIM ##

# the action is what will occur if the flag is used
# if check_only flag used, then it will be True; if not used, then False
parser.add_argument(
    '--check-only',
    action='store_true',
    help='Optional flag; if used, the sequences in the input file path will only be checked for the '
         'presence of primers in the forward and reverse complement orientation in the sequences. No '
         'primer trimming will occur.',
)

## STANDARD OUTPUT ##

# if option provided, will print out all warnings
parser.add_argument(
    '--verbose',
    action='store_true',
    default=settings['automate']['verbose'],
    help='Optional flag; if used, all errors and warnings from cutadapt will be printed to the '
         'standard output file or Terminal. If not used, then only fatal errors are shown.',
)

## PARSE OPTIONS INTO DICTIONARY ##

args = vars(parser.parse_args())

########################################################################################################################


#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# check if there are Illumina reads to trim primers from in the input directory
is_input, illumina_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

# if illumina files are detected...
if is_input:

    # determine whether primers should be located and trimmed, or just located

    # primers will just be located and reported on
    if args['check_only']:

        # warn that these settings mean that primers are only looked for in the sequences, not removed
        print(f'Only checking whether primers are still present in the input sequences, without trimming primers...\n')

        # only check if there are primers in the illumina sequences in the input path
        # confirm_no_primers(
        #     input_files=illumina_files,
        #     reference_dir=ref_dir,
        #     platform=platform,
        # )

    # primers will be located and trimmed from the sequences
    else:

        # notify that cutadapt is running to find and remove primers from the input illumina sequences
        print(f'Running cutadapt to find and remove primers...\n')

        # first run cutadapt to locate and remove primers from the input illumina sequences
        illumina_trimmed_path = remove_primers(
            input_files=illumina_files,
            output_dir=args['output'],
            reference_dir=ref_dir,
            platform=platform,
            paired_end=True,
            keep_removed_seqs=True,
        )

        # then use the primer check function to confirm that all primers in forward and reverse complement orientation
        #  are not detected in the trimmed sequence files produced in the previous cutadapt wrapper function
        # confirm_no_primers(
        #     input_files=illumina_trimmed_path,
        #     reference_dir=ref_dir,
        #     platform=platform,
        # )

# if no illumina files are detected, do nothing
else:
    pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# check if there are pacbio reads to trim primers from in the input directory
is_input, pacbio_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

# if pacbio files are detected...
if is_input:

    # determine whether primers should be located and trimmed, or just located

    # primers will just be located and reported on
    if args['check_only']:

        # warn that these settings mean that primers are only looked for in the sequences, not removed
        print(f'Only checking whether primers are still present in the input sequences, without trimming primers...\n')

        # only check if there are primers in the pacbio sequences in the input path
        # confirm_no_primers(
        #     input_files=pacbio_files,
        #     reference_dir=ref_dir,
        #     platform=platform,
        # )

    # primers will be located and trimmed from the sequences
    else:

        # notify that cutadapt is running to find and remove primers from the input pacbio sequences
        print(f'Running cutadapt to find and remove primers...\n')

        # first run cutadapt to locate and remove primers from the input pacbio sequences
        pacbio_trimmed_path = remove_primers(
            input_files=pacbio_files,
            output_dir=args['output'],
            reference_dir=ref_dir,
            platform=platform,
            paired_end=False,
            keep_removed_seqs=True,
        )

        # then use the primer check function to confirm that all primers in forward and reverse complement orientation
        #  are not detected in the trimmed sequence files produced in the previous cutadapt wrapper function
        # confirm_no_primers(
        #     input_files=pacbio_trimmed_path,
        #     reference_dir=ref_dir,
        #     platform=platform,
        # )

# if no pacbio files are detected, do nothing
else:
    pass

#####################
# SANGER ############
#####################
platform = 'sanger'

# check if there are pacbio reads to trim primers from in the input directory
is_input, sanger_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

# if sanger files are detected...
if is_input:

    # determine whether primers should be located and trimmed, or just located

    # primers will just be located and reported on
    if args['check_only']:

        # warn that these settings mean that primers are only looked for in the sequences, not removed
        print(f'Only checking whether primers are still present in the input sequences, without trimming primers...\n')

        # only check if there are primers in the sanger sequences in the input path
        # confirm_no_primers(
        #     input_files=sanger_files,
        #     reference_dir=ref_dir,
        #     platform=platform,
        # )

    # primers will be located and trimmed from the sequences
    else:

        # notify that cutadapt is running to find and remove primers from the input sanger sequences
        print(f'Running cutadapt to find and remove primers...\n')

        # first run cutadapt to locate and remove primers from the input sanger sequences
        sanger_trimmed_path = remove_primers(
            input_files=sanger_files,
            output_dir=args['output'],
            reference_dir=ref_dir,
            platform=platform,
            paired_end=False,
            keep_removed_seqs=True,
        )

        # then use the primer check function to confirm that all primers in forward and reverse complement orientation
        #  are not detected in the trimmed sequence files produced in the previous cutadapt wrapper function
        # confirm_no_primers(
        #     input_files=sanger_trimmed_path,
        #     reference_dir=ref_dir,
        #     platform=platform,
        # )

# if no sanger files are detected, do nothing
else:
    pass


# when all are primers trimmed, continue to next
continue_to_next(__file__, settings)