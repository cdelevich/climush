import argparse, pathlib
from pathlib import Path
from climush.bioinfo import cluster_reads
from climush.utilities import get_settings, continue_to_next, check_for_input


# create a reference directory path for file-finding functions
ref_dir=Path(__file__).parent

# import the settings for the bioinformatics configuration
settings = get_settings(ref_dir)
run_name = settings['run_details']['run_name']


## DEFINE + PARSE COMMAND LINE OPTIONS

# set up command line options
parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Cluster reads either within a sample or among a group of samples.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# input directory containing the files with sequences to cluster; no default, since depends on platform
parser.add_argument('-i', '--input',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to cluster. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py. This depends on the sequencing platform.')

# output directory to write clustered sequence files to
parser.add_argument('-o', '--output',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The path to a directory to write clustered sequences to. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# pairwise alignment clustering threshold
parser.add_argument('--clust_threshold',
                    default=settings['otu_clustering']['min_threshold'],
                    type=float,
                    help='')

# sequence sorting method to use prior to clustering
parser.add_argument('--method',
                    default=settings['otu_clustering']['method'],
                    choices=[1, 2, 3, 4],
                    help='The method of clustering to use. The methods vary in the way that they sort the sequences '
                         'before clustering (i.e., length sorts by sequence length, size by sequence abundance); '
                         'length_fast can only be used if reads are already sorted by sequence length.')

# how to gather sequences before clustering
# parser.add_argument('--group_by',
#                     default=settings['clustering']['group_by'],
#                     choices=['community', 'sample'],
#                     help='The method of gathering together reads prior to clustering. Reads can be clustered within '
#                          'a sample, which is typical for PacBio sporocarp sequencing, or clustered across all samples '
#                          'within a provided directory of sequencing files, which is typical for Illumina sequencing.')

# # how to chose the representative read for an OTU cluster
# parser.add_argument('--choose_rep_by',
#                     default=settings['clustering']['choose_rep_by'],
#                     choices=['read count', 'read length'],
#                     help='How to determine which sequence to use as the representative read for an OTU. It can be '
#                          'by abundance (read count) or longest read (read length); if nothing is provided, then the '
#                          'standard protocol of choosing the representative read for an OTU is used.')
#
# # how to chose the representative OTU for a pacbio sporocarp sample
# parser.add_argument('--choose_top_by',
#                     default=settings['clustering']['choose_top_by'],
#                     choices=['read count', 'read length'],
#                     help='How to determine which OTU to use as the representative OTU for a sample. It can be '
#                          'by abundance (read count) or longest read (read length); if nothing is provided, then the '
#                          'standard protocol of choosing the representative OTU for a sample is used. This is only '
#                          'relevant for PacBio sporocarp sequences.')
#
#
# # option to write each cluster out to a separate .fasta file
# parser.add_argument('--separate_clusters',
#                     default=settings['clustering']['separate_clusters'],
#                     store_value='true',
#                     help='If flag used, the OTUs formed from clustering will be written out to separate sequence '
#                          'files, using the sample ID as the prefix.')

# parse command line options and defaults into a dictionary
args = vars(parser.parse_args())


#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# check that there are Illumina reads to cluster
is_input, illumina_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

if is_input:
    cluster_reads(
        input_files=illumina_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        clust_threshold=args['clust_threshold'],
        clust_method=args['method'],
        group=False,
    )
else:
    pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# check that there are pacbio reads to cluster
is_input, pacbio_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)


if is_input:
    cluster_reads(
        input_files=pacbio_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        clust_threshold=args['clust_threshold'],
        clust_method=args['method'],
        group=True,
    )
else:
    pass

#####################
# SANGER ############
#####################
platform = 'sanger'

# check that there are Illumina reads to cluster
is_input, sanger_files = check_for_input(
    file_dir=args['input'],
    config_dict=settings,
    file_identifier=platform,
)

if is_input:
    cluster_reads(
        input_files=sanger_files,
        output_dir=args['output'],
        reference_dir=ref_dir,
        clust_threshold=args['clust_threshold'],
        clust_method=args['method'],
        group=True,
    )
else:
    pass

# when all are seqs are dereplicated, continue to next
continue_to_next(__file__, settings)