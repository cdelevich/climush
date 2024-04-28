from mapping import filepath_map as fpm

import argparse, sys, subprocess
from pathlib import Path
from climush.constants import *
from climush.bioinfo import create_blast_db, assign_taxonomy
from climush.utilities import *

settings = get_settings(fpm)

#####################
# ILLUMINA ##########
#####################

# check that there are Illumina reads to dereplicate
last_output = [dir for dir in fpm['pipeline-output']['otus-clustered'].glob('*') if
               re.search(f'^{CLUSTER_PREFIX}', dir.stem, re.I)][0]
is_input, illumina_files = check_for_input(last_output)

if is_input:

    create_blast_db(config_dict=settings, file_map=fpm, taxa_list=None)
else:
    pass
#####################
# PACBIO ############
#####################

# check that there are PacBio reads to dereplicate
# last_output = [dir for dir in fpm['pipeline-output']['otus-clustered'].glob('*') if
#                re.search(f'^{DEREP_PREFIX}02', dir.stem, re.I)][0]
# is_input, illumina_files = check_for_input(last_output)
#
# if is_input:
#     create_blast_db(config_dict, file_map, taxa_list=None)
# else:
#     pass

#####################
# SANGER ############
#####################


# when all are seqs are dereplicated, continue to next
continue_to_next(Path(__file__), settings)
