from mapping import filepath_map as fpm

import argparse
import sys
import re
from pathlib import Path
##REMOVE AFTER PACKAGE TESTING#######
sys.path.insert(0, str(fpm['package']['main']))
#####################################
from climush.constants import *
from climush.utilities import import_config_as_dict, count_files, continue_to_next, check_for_input

config_dict = import_config_as_dict(fpm['config']['main'], file_handle='pipeline-settings', config_section='all')

# check for sequences that need to be demultiplexed
demux_dir = fpm['sequences']['demux']
is_input, file_list = check_for_input(demux_dir, file_ext='.fast*')
if is_input:
    print(f'{len(file_list)} files require demultiplexing...\n')
else:
    print(f'No sequences requiring demultiplexing were detected.\n')
    continue_to_next(Path(__file__), config_dict)

#####################
# ILLUMINA ##########
#####################

# check that there are Illumina reads to pre-filter
is_input, illumina_files = check_for_input(fpm['sequences']['demux'])

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


#####################################################
# VARIABLES AND PATHS REQUIRED BY THE DEMUX FUNCTIONS
#####################################################

# LIKELY GOOD TO GO
run_name = config_dict['run_details']['run_name']
quality_score = str(config_dict['quality_filtering']['pacbio']['qscore'])
mapping_files = fpm['config']['bc_mapping'].glob('*')
output_path = fpm['pipeline-output']['demultiplexed']
raw_read_path = fpm['sequences']['demux']

# REQUIRE ATTENTION
barcode_path = ''
final_demux_path = ''
primer_path = config_dict['primers']
# primer_path = "../miscell/primers/"

###################################
### FROM 01_CREATE-BARCODE-FASTA.PY
###################################
import pandas as pd
from pathlib import Path

mapping_tabs = pd.read_excel(mapping_file, sheet_name=None)
for tab in mapping_tabs.keys():
    tab_name_fixed = re.search('(pool\d)', tab, re.IGNORECASE)[0].lower() + '_df'
    exec(f'{tab_name_fixed} = pd.read_excel(mapping_file, sheet_name=tab)')

# define single function to carry out barcode file creation
def create_barcode_fasta(barcode_csv, pool_num):
    """
    Create a fasta file of unique forward and reverse barcodes from an input mapping file that
    should contain a column for the forward barcode, reverse barcode, and sample ID. Column naming
    is flexible, as function will look for any column name with 'fwd' or 'forward', 'rev' or
    'reverse', and 'sample' and 'id', 'no', 'num', with or without space or underscore, and case
    insensitive. Barcodes may have primers attached in these tables, which will be removed, if present,
    in this function.
    :param barcode_csv: input barcode file already read in
    :param pool_num: if sequencing run was split into multiple pools, include the pool number
    :return: no direct output, but will save a [1] fasta file of the barcodes and [2] a csv of
    sample IDs and their barcode combinations for in further steps of demultiplexing.
    """
    fwd_bc_regex = re.compile('(fwd)|(forward)', re.IGNORECASE)
    rev_bc_regex = re.compile('(rev)|(reverse)', re.IGNORECASE)
    sample_regex = re.compile('(sample)(_|\s)(id|no|num)', re.IGNORECASE)
    forward_bc_col = list(filter(fwd_bc_regex.search, barcode_csv.columns))[0]
    reverse_bc_col = list(filter(rev_bc_regex.search, barcode_csv.columns))[0]
    sample_col = list(filter(sample_regex.search, barcode_csv.columns))[0]

    # get the number of unique barcodes and the number of unique samples
    num_unique_fwd = len(barcode_csv[forward_bc_col].unique())
    num_unique_rev = len(barcode_csv[reverse_bc_col].unique())
    num_unique_samples = len(barcode_csv[sample_col].unique())

    # rename column names
    barcode_csv.rename(columns={sample_col: 'sample_id',
                                forward_bc_col: 'forward',
                                reverse_bc_col: 'reverse'},
                       inplace=True)

    # remove primers from barcodes
    with open(primer_path + "primer_fwd.fasta", 'rt') as fin:
        fwd_primer = fin.readlines()[1]
    with open(primer_path + "primer_rev.fasta", 'rt') as fin:
        rev_primer = fin.readlines()[1]
    barcode_csv['forward'] = [i.split(fwd_primer)[0] for i in barcode_csv['forward']]
    barcode_csv['reverse'] = [i.split(rev_primer)[0] for i in barcode_csv['reverse']]

    # create df that contains the sample IDs with the sequences of the barcodes
    wanted_cols = ['sample_id','forward','reverse']
    sampleid_df = pd.melt(barcode_csv[wanted_cols], id_vars='sample_id', var_name='barcode_direction', value_name='sequence')
    assert (len(barcode_csv['sample_id']) == len(sampleid_df['sample_id']) / 2), 'Dataframe likely contains extra columns that' \
                                                                                 'were not be removed from the dataframe (though ' \
                                                                                 'the code should do this).'

    # create a set of barcodes (set is list of unique values)
    fwd_bc = list(set(barcode_csv['forward']))
    rev_bc = list(set(barcode_csv['reverse']))

    barcode_df = sampleid_df.drop('sample_id', axis=1).drop_duplicates(subset='sequence').sort_values(
        by='sequence').reset_index().drop('index', axis=1)

    # make sure each sequence is unique in df
    assert len(fwd_bc) + len(rev_bc) == barcode_df.shape[0]
    assert len(barcode_df['sequence'].unique()) == len(barcode_df['sequence'])

    # create headers for barcode fasta file using the indices of the *sorted* dataframe, barcode_df
    # must be sorted to produce the same file each time, which is crucial for the way Lima works
    barcode_df['barcode_header'] = None

    for i in range(barcode_df.shape[0]):
        if barcode_df['barcode_direction'][i] == 'forward':
            barcode_df['barcode_header'][i] = "bc_fwd_" + str(barcode_df.index[i])
        else:
            barcode_df['barcode_header'][i] = "bc_rev_" + str(barcode_df.index[i])

    # write out barcodes and headers as fasta file
    fasta_output_name = output_path + run_name + "pool" + str(pool_num) + "_barcodes.fasta"
    with open(fasta_output_name, 'wt') as fout:
        for i in range(barcode_df.shape[0]):
            fout.write(">" + barcode_df['barcode_header'][i] + '\n')
            fout.write(barcode_df['sequence'][i] + '\n')

    if os.path.isfile(fasta_output_name) == True:
        print(f"\nBarcode fasta file for pool {pool_num}, {os.path.basename(fasta_output_name)}, has been created.\n")

    # merge sample id dataframe and barcode id dataframe by sequence
    merged_df = pd.merge(barcode_df, sampleid_df)
    wider_df = pd.pivot(data=merged_df, index='sample_id', columns='barcode_direction', values='barcode_header').reset_index()
    wider_df['fwd_index'] = [wider_df['forward'][i].split('_')[2] for i in range(wider_df.shape[0])]
    wider_df['rev_index'] = [wider_df['reverse'][i].split('_')[2] for i in range(wider_df.shape[0])]

    csv_output_name = output_path + run_name + "pool" + str(pool_num) + "_merged.csv"
    wider_df.to_csv(csv_output_name, index=False)

    assert len(wider_df['forward'].unique()) == num_unique_fwd
    assert len(wider_df['reverse'].unique()) == num_unique_rev
    assert len(wider_df.index.unique()) == num_unique_samples

    if os.path.isfile(csv_output_name) == True:
        print(f"Merged sample ID and barcode name file for pool {pool_num}, {csv_output_name}, has been created.\n")

# run command for however many pools you might have (can loop if several)
create_barcode_fasta(barcode_csv=pool1_df, pool_num=1)
create_barcode_fasta(barcode_csv=pool2_df, pool_num=2)

###################################
### FROM 02_RUN-LIMA.PY
###################################

# write timer functions
def start_runtime():
    global start_time
    start_time = datetime.now()
    return None

def end_runtime():
    global end_time
    end_time = datetime.now()
    return None

def print_runtime(custom_text=None):
    runtime_total_seconds = (end_time - start_time).total_seconds()
    runtime_mins = int(runtime_total_seconds/60)
    runtime_secs = round(runtime_total_seconds%60, 2)
    if custom_text is None:
        print(f"The program ran in {runtime_mins} min and {runtime_secs} sec.\n")
    else:
        print(f"{custom_text} in {runtime_mins} min and {runtime_secs} sec.\n")
    return None

# END SOURCE FILE INFO ###########


print(f'\nUsing the barcode fasta files located in {barcode_path} and the raw '
      f'sequence files located in {raw_read_path}. If you\'d '
      f'like to use another folder containing barcodes or raw reads, please see options by using the '
      f'--help function after the name of this script (i.e., python3 02_run-lima.py --help)\n')

def get_rawread_filename(pool_num):
    filename_regex = re.compile(('^(\d{4}\.Pool' + str(pool_num) + ')\S+(\.fastq\.gz)$'))
    rawread_file = list(filter(filename_regex.search, os.listdir(raw_read_path)))
    assert len(rawread_file) == 1, f'{len(rawread_file)} .fastq.gz files were located in the ' \
                                   f'provided filepath: {rawread_file} for pool {pool_num}. Requires exactly 1 ' \
                                   f'file per pool.'
    return raw_read_path + rawread_file[0]

def execute_lima_demux(pool_num):
    '''
    Composes the command line script for Lima demultiplexing, and executes it.
    :param run_prefix: four-digit prefix on sequencing run files, provided by sequencer
    :param pool_num: sequencing pool number, if split into more than one pool
    :return:
    '''
    start_runtime()
    rawread_file = get_rawread_filename(pool_num)
    barcode_file = barcode_path + run_name + "pool" + str(pool_num) + "_barcodes.fasta"
    output = barcode_path + run_name + "pool" + str(pool_num) + "_demux.fasta"

    lima_cmd = "lima " + rawread_file + " " + barcode_file + " " + output + \
               " --hifi-preset ASYMMETRIC" # asymmetric means forward and reverse barcodes are unique

    os.system(lima_cmd)
    end_runtime()
    print_runtime(custom_text=f'Lima demultiplexed samples from pool {pool_num}')

execute_lima_demux(pool_num = 1)
execute_lima_demux(pool_num = 2)

###################################
### FROM 03_DEMULTIPLEX.PY
###################################

final_demux_path = f"./final-demux_{run_name}/"  # new file path (created next lines)

def get_barcode_names(pool_num):
    '''
    Get the name of the barcodes from the header of the unique
    barcode fasta files used in Lima demultiplexing.
    :param pool_num: the number of the PacBio sequencing pool (1/2)
    :return: list of the name of the barcodes
    '''
    barcode_names = []

    barcode_file = f"{output_path}{run_name}pool{str(pool_num)}_barcodes.fasta"

    with open(barcode_file) as fin:
        for line in fin:
            if line.startswith(">"):
                barcode_names.append(line.strip(">").strip("\n"))
    fin.close()

    return barcode_names

def create_post_lima_df(pool_num):
    '''
    Identifies barcode pairs based on index values provided by Lima, reorients the barcodes,
    and returns a dataframe of the read IDs, their sequences, and their barcodes.
    :param pool_num: sequencing pool number
    :return: dataframe of read IDs and their barcode combinations, ready for demultiplexing
    '''
    pool_num = str(pool_num)
    print(f"\nCreating a dataframe of the Lima output for pool {pool_num}...")
    start_runtime()
    barcode_list = get_barcode_names(pool_num) # get list of barcodes

    seq_sample_id = []  # read IDs
    bc_index_01 = []  # first barcode detected by Lima (returns index value)
    bc_index_02 = []  # second barcode detected by Lima (returns index value)
    seq = []  # sequence belonging to read ID, with barcodes removed by Lima

    # output of script 02_run-lima.py
    demux_fasta = f"{output_path}{run_name}pool{pool_num}_demux.fasta"

    with open(demux_fasta) as fin:
        for line in fin.readlines():
            if line.startswith(">"):
                split_line = line.split(" ")
                seq_sample_id.append(split_line[0][1:])
                bc_index_01.append(split_line[1].split(',')[0][3:])
                bc_index_02.append(split_line[1].split(',')[1])
            else:
                seq.append(line)
    fin.close()

    post_lima = pd.DataFrame({'seq_sample_id': seq_sample_id, # see notes for lists above
                              'bc_index_01': bc_index_01,
                              'bc_index_02': bc_index_02,
                              'dna_seq': seq})

    post_lima['bc_01_id'] = None  # create empty column to add barcode ID to
    post_lima['bc_02_id'] = None

    for i in range(post_lima.shape[0]):
        index_01 = post_lima['bc_index_01'][i]  # get the index of the barcode in the fasta file given to Lima
        index_02 = post_lima['bc_index_02'][i]
        for j in barcode_list:  # for each name in the barcode list (headers in fasta file given to Lima)
            bc_end = j.split("_")[2]  # only take the index (last part of the header string)
            if bc_end == index_01:  # if this index matches Lima's index assignment
                post_lima['bc_01_id'][i] = j  # add the full barcode name to the 'fwd' index id column
            if bc_end == index_02:  # if it makes the second index...
                post_lima['bc_02_id'][i] = j  # add to second index column (may be both sometimes, so two if statements, no elif)

    # some barcodes are oriented incorrectly (e.g., first barcode is reverse, second is forward)
    # to complete demultiplexing, need to reorient them, into new columns in df
    post_lima['bc_fwd_id'] = None
    post_lima['bc_rev_id'] = None
    post_lima['bc_index_fwd'] = None
    post_lima['bc_index_rev'] = None

    # reorient barcodes
    mismatch_count = 0
    for i in range(post_lima.shape[0]):
        bc_id_01 = post_lima['bc_01_id'][i]
        bc_id_01_dir = bc_id_01.split("_")[1]
        bc_id_02 = post_lima['bc_02_id'][i]
        bc_id_02_dir = bc_id_02.split("_")[1]
        if bc_id_01_dir == bc_id_02_dir:  # remove the ones that have two fwd or two rev barcodes
            mismatch_count += 1
            post_lima['bc_fwd_id'][i] = "NA"
            post_lima['bc_rev_id'][i] = "NA"
        else:
            if bc_id_01.split("_")[1] == "fwd":
                post_lima['bc_fwd_id'][i] = bc_id_01
                post_lima['bc_index_fwd'][i] = bc_id_01.split("_")[2]
                post_lima['bc_rev_id'][i] = bc_id_02
                post_lima['bc_index_rev'][i] = bc_id_02.split("_")[2]
            else:
                post_lima['bc_fwd_id'][i] = bc_id_02
                post_lima['bc_index_fwd'][i] = bc_id_02.split("_")[2]
                post_lima['bc_rev_id'][i] = bc_id_01
                post_lima['bc_index_rev'][i] = bc_id_01.split("_")[2]

    end_runtime()
    print_runtime(custom_text=f"The post-Lima dataframe for pool {pool_num} was created")

    return post_lima

post_lima_pool1 = create_post_lima_df(pool_num=1)
post_lima_pool2 = create_post_lima_df(pool_num=2)

# import the counts post-lima info to confirm that all reads were loaded into post_lima df
def check_counts(post_lima_df, pool_num):
    pool_num = str(pool_num)
    demux_lima_counts = f"{output_path}{run_name}pool{pool_num}_demux.lima.counts"
    post_lima_counts = pd.read_csv(demux_lima_counts, sep='\t')['Counts'].sum()
    if (post_lima_counts == post_lima_df.shape[0]):
        print(f"\nAll reads are loaded into the post lima dataframe for pool {pool_num}.")
    else:
        print(f"\n**ERROR**: Not all reads are loaded into the post lima dataframe for "
              f"pooll {pool_num}. Exiting...\n")
        sys.exit()

check_counts(post_lima_pool1, pool_num=1)
check_counts(post_lima_pool2, pool_num=2)

def create_pre_lima_df(pool_num):
    pool_num = str(pool_num)

    print(f"\nCreating a dataframe of the Lima input barcodes for {pool_num}, which allows "
          f"use to reorient barcode combinations that may be flipped...")
    start_runtime()

    merged_csv = f"{output_path}{run_name}pool{pool_num}_merged.csv"
    merged_df = pd.read_csv(merged_csv)
    merged_df.rename(columns = {'fwd_index': 'bc_index_fwd',
                                'rev_index': 'bc_index_rev'}, inplace=True)

    end_runtime()
    print_runtime(custom_text=f"The pre-Lima dataframe for pool {pool_num} was created")

    return merged_df

pre_lima_pool1 = create_pre_lima_df(pool_num=1)
pre_lima_pool2 = create_pre_lima_df(pool_num=2)

def get_bc_combo(df): # fwd and rev index columns must be named 'bc_index_fwd'
    '''
    Make a combination of the forward and reverse barcodes using the index number of the barcodes
    given to Lima in a fasta file, required in order to merge sample IDs with Lima output.
    :param df: dataframe containing barcodes
    :return: none, will add column to the input dataframe
    '''
    df['bc_combo'] = [ str(df['bc_index_fwd'][i]) + "_" + str(df['bc_index_rev'][i]) for i in range(df.shape[0])]
    return None

get_bc_combo(pre_lima_pool1)
get_bc_combo(post_lima_pool1)

get_bc_combo(pre_lima_pool2)
get_bc_combo(post_lima_pool2)

def final_demux(pre_lima_df, post_lima_df, pool_num):
    pool_num = str(pool_num)
    print(f"\nDemultiplexing samples in pool {pool_num}...")
    start_runtime()

    lima_merge_df = pre_lima_df.merge(post_lima_df, how='inner', on='bc_combo')

    sample_ids = set(lima_merge_df['sample_id'])
    read_for_count = []
    sample_for_count = []

    for i,s in enumerate(sample_ids):
        print(f"\r    sample {s} ({i+1}/{len(sample_ids)})--------", end='')
        sample_fasta_name = f"{final_demux_path}{run_name}{s}.fasta"
        read_count = 0
        for i in range(lima_merge_df.shape[0]):
            if lima_merge_df['sample_id'][i] == s:
                read_count += 1
                if os.path.isfile(sample_fasta_name) == True:
                    file_mode = "a"
                else:
                    file_mode = "w"
                with open(sample_fasta_name, file_mode) as fasta_out:
                    fasta_out.write(">" + lima_merge_df['seq_sample_id'][i] + "\n")
                    fasta_out.write(lima_merge_df['dna_seq'][i])
        read_for_count.append(read_count)
        sample_for_count.append(s)

    count_summary = pd.DataFrame({'seq_pool': pool_num,
                                  'sample_id': sample_for_count,
                                  'read_count': read_for_count})

    count_summary.to_csv(f"{final_demux_path}{run_name}pool{pool_num}_read-counts.csv", index=False)
    end_runtime()
    print_runtime(custom_text=f'\nSamples in pool {pool_num} have been demultiplexed')

final_demux(pre_lima_pool1, post_lima_pool1, pool_num=1)
final_demux(pre_lima_pool2, post_lima_pool2, pool_num=2)



