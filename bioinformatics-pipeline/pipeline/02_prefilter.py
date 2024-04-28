from mapping import filepath_map as fpm

import argparse, sys, subprocess
from pathlib import Path
from climush.constants import *
from climush.bioinfo import filter_out_phix, prefilter_fastx
from climush.utilities import *

settings = get_settings(fpm)

# divided into sections instead of looping through each platform to maintain a similar structure among all scripts, as
# sometimes (like in prefiltering) there are different processes based on the platform

#####################
# ILLUMINA ##########
#####################

# check that there are Illumina reads to pre-filter
is_input, illumina_files = check_for_input(fpm['sequences']['illumina'])

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

# unfamiliar with whether any prefiltering necessary for Sanger reads, any PhiX spike-in?

is_input, sanger_files = check_for_input(fpm['sequences']['sanger'])

if is_input:
    nophix_path = filter_out_phix(input_files=sanger_files, file_map=fpm)
    nophix_files = list(nophix_path.glob(SEQ_FILE_GLOB))
    noambig_path = prefilter_fastx(input_files=nophix_files, file_map=fpm, maxn=0)
else:
    pass

# when all are prefiltered, continue to next
continue_to_next(Path(__file__), settings)


