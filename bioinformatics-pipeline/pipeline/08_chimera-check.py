import argparse, pathlib
from pathlib import Path

from climush.bioinfo import check_chimeras
from climush.utilities import get_settings, check_for_input, continue_to_next

# set a location to start looking for pipeline configuration file
ref_dir = Path(__file__).parent.parent

settings = get_settings(ref_dir)
run_name = settings['run_details']['run_name']

# set up command line options
parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Detect and remove chimeras.',
                                 epilog='This script is part of the SPUN bioinformatics pipeline.')

# input directory containing the files to check for chimeras; no default, since depends on platform
parser.add_argument('-i', '--input',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to chimera check.')

parser.add_argument('-o', '--output',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The file path that chimera-checked sequences will be written out to.')

# method to use for chimera detection
parser.add_argument('--method',
                    default=settings['chimera_check']['method'],  # depends on previous pipeline steps, settings
                    choices=['denovo', 'reference'],
                    help='The method of chimera detection; either denovo or reference. If denovo, Illumina reads will '
                         'first be denoised using the UNOISE algorithm. If reference-based, sequences will be compared '
                         'against the UNITE reference dataset for UCHIME.')

# the minimum abundance skew
parser.add_argument('--alpha',
                    default=None,  # if not provided, depends on method for default value
                    help='The minimum abundance skew; the factor by which read abundance of centroid exceeds chimera. '
                         'For de novo chimera detection using UCHIME3, the default is 16.0. For all other chimera '
                         'detection methods, the default is 2.0.')

# optional; whether to keep the chimera reads and write out to chimera output directory
parser.add_argument('--keep-chimeras',
                    default=settings['chimera_check']['keep_chimeras'],
                    action='store_true',
                    help='Whether to keep the chimera sequences and write them to their own output directory. If '
                         'flag is used, will write output files for chimera sequences. Files are created regardless '
                         'of whether chimeras are detected in a given sample. If flag is not used, will default to '
                         'whatever the setting is in the bioinformatics configuration file.')

# parse command line options and defaults into a dictionary
args = vars(parser.parse_args())


# determine the default alpha value to use if one is not provided to argparse
if args['alpha'] is None:
    default_alpha = settings['chimera_check']['alpha'][args['method']]
    args.update({'alpha': default_alpha})
    print(f'Using a default abundance skew of {default_alpha} for {args["method"]} chimera detection.\n')

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
    check_chimeras(
        input_files=illumina_files,
        reference_dir=ref_dir,
        output_dir=args['output'],
        method=args['method'],
        alpha=args['alpha'],
        keep_chimeras=args['keep_chimeras'],
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
    check_chimeras(
        input_files=sanger_files,
        reference_dir=ref_dir,
        output_dir=args['output'],
        method=args['method'],
        alpha=args['alpha'],
        keep_chimeras=args['keep_chimeras'],
    )
else:
    pass


#####################
# PACBIO ############
#####################
platform = 'pacbio'

# check if there are ITS-LSU sequences in the input directory
is_input, pacbio_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
    file_ext=None,
    )

if is_input:
    check_chimeras(
        input_files=pacbio_files,
        reference_dir=ref_dir,
        output_dir=args['output'],
        method=args['method'],
        alpha=args['alpha'],
        keep_chimeras=args['keep_chimeras'],
    )
else:
    pass


# continue to next
continue_to_next(__file__, settings)