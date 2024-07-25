from mapping import filepath_map as fpm

import argparse, sys, subprocess, pathlib
from pathlib import Path
from climush.constants import *
from climush.bioinfo import filter_out_phix, prefilter_fastx
from climush.utilities import *

settings = get_settings(fpm)

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Remove PhiX and reads with ambiguous bases',
                                 epilog='')

parser.add_argument('-i', '--input', default=fpm['sequences'],
                    type=pathlib.PosixPath,
                    help='The path to the sequencing files. Will default to the location that is '
                         'expected with the Docker container\'s native file structure.')

args = vars(parser.parse_args())

# divided into sections instead of looping through each platform to maintain a similar structure among all scripts, as
# sometimes (like in prefiltering) there are different processes based on the platform

#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# check that there are Illumina reads to pre-filter
# is_input, illumina_files = check_for_input(args['input'] / platform)

# print(f'WARNING. Hard-coded this to gather all soil/litter Illumina samples, not reproducible at the moment.\n')

# all_dirs = [Path('/home/cdelevic/roylab/globus_climush-sequences/illumina/illumina_soil-litter/illumina_soil-litter_2022-05/illumina_soil-litter_2022-05_raw-reads'),
#             Path('/home/cdelevic/roylab/globus_climush-sequences/illumina/illumina_soil-litter/illumina_soil-litter_2022-10/illumina_soil-litter_2022-10_raw-reads'),
#             Path('/home/cdelevic/roylab/globus_climush-sequences/illumina/illumina_soil-litter/illumina_soil-litter_2022-10/non-climush_sequences'),
#             Path('/home/cdelevic/roylab/globus_climush-sequences/illumina/illumina_soil-litter/illumina_soil-litter_2023-05/illumina_soil-litter_2023-05_raw-reads'),
#             Path('/home/cdelevic/roylab/globus_climush-sequences/illumina/illumina_soil-litter/illumina_soil-litter_2023-10/illumina_soil-litter_2023-10_raw-reads'),
#             Path('/home/cdelevic/roylab/globus_climush-sequences/illumina/illumina_soil-litter/illumina_soil-litter_2023-10/non-climush_sequences/horseshoe')]

print(f'WARNING. Hard-coded this to gather all Illumina samples from illumina_soil-litter_2023-10, '
      f'not reproducible at the moment.\n')

base_dir = Path('/home/cdelevic/roylab/globus_climush-sequences/illumina/illumina_soil-litter/illumina_soil-litter_2023-10/')
all_dirs = [base_dir / 'illumina_soil-litter_2023-10_raw-reads',
            base_dir / 'non-climush_sequences/horseshoe']

is_input = True
illumina_files = []
for dir in all_dirs:
    for file in dir.glob(GZIP_GLOB):
        if not file.name.startswith(platform):
            new_name = 'illumina_' + file.name
            new_path = dir / new_name
            file.rename(new_path)

        illumina_files.append(file)

print(f'Processing sequences from {len(illumina_files)} samples...\n')

if is_input:
    nophix_path = filter_out_phix(input_files=illumina_files, file_map=fpm)
    nophix_files = list(nophix_path.glob(SEQ_FILE_GLOB))
    noambig_path = prefilter_fastx(input_files=nophix_files, file_map=fpm, maxn=0)
else:
    pass

#####################
# PACBIO ############
#####################

# no pre-filtering for PacBio reads

#####################
# SANGER ############
#####################
platform = 'sanger'
# unfamiliar with whether any prefiltering necessary for Sanger reads, any PhiX spike-in?

is_input, sanger_files = check_for_input(args['input'] / platform)

if is_input:
    nophix_path = filter_out_phix(input_files=sanger_files, file_map=fpm)
    nophix_files = list(nophix_path.glob(SEQ_FILE_GLOB))
    noambig_path = prefilter_fastx(input_files=nophix_files, file_map=fpm, maxn=0)
else:
    pass

# when all are prefiltered, continue to next
continue_to_next(__file__, settings)


