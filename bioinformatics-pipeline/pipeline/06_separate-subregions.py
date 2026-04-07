import argparse, pathlib
from pathlib import Path
from climush.bioinfo import separate_subregions, concat_regions
from climush.utilities import check_for_input, get_settings, continue_to_next

## IMPORT PIPELINE CONFIGURATION #######################################################################################

# create a reference directory path for file-finding functions
ref_dir=Path(__file__).parent

# import the settings for the bioinformatics configuration
settings = get_settings(ref_dir)

########################################################################################################################


## COMMAND LINE ARGUMENTS ##############################################################################################

## INSTANTIATE PARSER ##

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Identify and separate the ITS/LSU subregions.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

## FILE PATHS ##

# path to the directory containing the sequencing files to run through itsx
parser.add_argument('-i', '--input',
                    required=True,
                    type=pathlib.PosixPath,
                    help='The path to the directory containing the sequence files that need to be separated into '
                         'subregions by ITSx, if the --concat-only flag is not used. If the --concat-only flag '
                         'is used here, then this is the path to the sample directories created from ITSx that '
                         'contains the per-subregion sequence files to concatenate into longer reads.')

parser.add_argument('-o', '--output',
                    required=True,
                    type=pathlib.PosixPath,
                    help='The path to the directory in which the output files and directories will be written.')

## DO NOT RUN ITSX ##

parser.add_argument('-c', '--concat-only',
                    action='store_true',
                    help='If this flag is used, ITSx will not run but the ITSx output files will be concatenated '
                         'into: (i) full-length reads by combining all re-oriented subregions and (ii) full ITS '
                         'regions by concatenating the ITS1, 5.8S, and ITS2 subregions. When this flag is used, '
                         'the required -i / --input parameter should be the path to the directory that contains '
                         'the output of the post-ITSx sample directories.')

## PARSE OPTIONS INTO DICTIONARY ##

args = vars(parser.parse_args())

########################################################################################################################


#####################
# ILLUMINA ##########
#####################

# may want to use this to confirm that reads are ITS1? could check orientation as well?
# I don't think there's a use for concatenating regions for Illumina? at least for CliMush use

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# if only concatenating sequences...
if args['concat_only']:

    # use the path provided to -i / --input as location of seqs to concatenate
    for itsx_sample in args['input'].glob('*'):

        # print which sample is currently being processed
        print(f'\n{itsx_sample.stem}\n')

        # concatenate full ITS sequence
        concat_regions(
            dir_path=itsx_sample,
            reference_dir=ref_dir,
            platform=platform,
            regions_to_concat=['ITS1', '5_8S', 'ITS2'],
        )

        # concatenate full-length ITS-LSU sequence
        concat_regions(
            dir_path=itsx_sample,
            reference_dir=ref_dir,
            platform=platform,
            regions_to_concat=['ITS1', '5_8S', 'ITS2', 'LSU'],
        )

        # check_concat_output(itsx_dir=concat_path, full_len_dir=args['input'], num_bp_compare=50)

# if running itsx prior to concatenating sequences...
else:

    # check for ITS-LSU sequences in the input directory
    is_input, pacbio_files = check_for_input(
        file_dir=args['input'],
        config_dict=settings,
        file_identifier=platform,
    )

    # if ITS-LSU sequences are located...
    if is_input:

        # run ITSx on the sequences
        itsx_out_path = separate_subregions(
            input_files=pacbio_files,
            output_dir=args['output'],
            reference_dir=ref_dir,
            verbose=True,
        )

        # using the output path where ITSx wrote files, concatenate the ITSx sequence output

        # concatenate full ITS sequence
        concat_regions(
            dir_path=itsx_out_path,
            reference_dir=ref_dir,
            platform=platform,
            regions_to_concat=['ITS1', '5_8S', 'ITS2'],
        )

        # concatenate full-length ITS-LSU sequence
        concat_regions(
            dir_path=itsx_out_path,
            reference_dir=ref_dir,
            platform=platform,
            regions_to_concat=['ITS1', '5_8S', 'ITS2', 'LSU'],
        )

        # check_concat_output(itsx_dir=itsx_out_path, full_len_dir=args['input'], num_bp_compare=50)

    else:
        pass

#####################
# SANGER ############
#####################


# when all are seqs are dereplicated, continue to next
continue_to_next(__file__, settings)