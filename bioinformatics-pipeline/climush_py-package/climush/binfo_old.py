import re
import os
import datetime
from functools import wraps

# get VSEARCH version
# subprocess run
vsearch_version = os.popen('vsearch --version').readlines()[0].split(' ')[1].split('_')[0]

# get location of the mothur executable file
mothur_path = os.popen('which mothur').read().strip()
mothur_version = os.popen('mothur --version').readlines()[1].split('=')[-1]

def func_timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = datetime.datetime.now()
        func_output = func(*args, **kwargs)
        end_time = datetime.datetime.now()
        runtime = str(end_time - start_time).split('.')[0]  # round seconds down
        print(f"{func.__name__} was executed in {runtime}.\n")
        return func_output
    return wrapper()

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

def get_bc_combo(df): # fwd and rev index columns must be named 'bc_index_fwd'
    '''
    Make a combination of the forward and reverse barcodes using the index number of the barcodes
    given to Lima in a fasta file, required in order to merge sample IDs with Lima output.
    :param df: dataframe containing barcodes
    :return: none, will add column to the input dataframe
    '''
    df['bc_combo'] = [ str(df['bc_index_fwd'][i]) + "_" + str(df['bc_index_rev'][i]) for i in range(df.shape[0])]
    return None

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

# check which samples didn't make it to this final pool post-processing
def get_missing_samples(pool_num):
    '''
    Check which samples were missing from the final demultiplexing for a
    given quality score folder. Will provide a print out of the missing
    samples, and also saves a table of the missing samples.
    :param pool_num: number of the sequencing pool
    :return: saves a text file of missing samples in the pool in the
    final-demux folder
    '''
    pool_num = str(pool_num)

    print(f"\nChecking which (if any) samples are excluded from the Q{quality_score} "
          f"reads from pool {pool_num} (due to quality filtering)...\n")

    input_samples = pd.read_csv(f"{output_path}{run_name}pool{pool_num}_merged.csv")['sample_id'].tolist()
    output_samples = pd.read_csv(f"{final_demux_path}{run_name}pool{pool_num}_read-counts.csv")['sample_id'].tolist()

    missing_samples = set(input_samples).difference(set(output_samples))

    if len(missing_samples) > 0:
        print(f"    The following samples are missing from the Q{quality_score} samples in pool {pool_num}:")
        for sample in missing_samples:
            print(f"        - {sample}")

    missing_samples_out = f"{output_path}{run_name}pool{pool_num}_missing-samples.txt"
    with open(missing_samples_out, 'wt') as fout:
        for sample in missing_samples:
            fout.write(f"sample\n")
    fout.close()

    print(f"\n    The missing sample IDs have been saved to: {missing_samples_out}")

def execute_cutadapt(demux_fasta):
    demux_filepath = final_demux_path + demux_fasta
    trim_output = f"{output_path}trim_{demux_fasta}"
    untrim_output = f"{notrim_output_path}untrimmed_{demux_fasta}"

    cutadapt_cmd = f"cutadapt -g ^{fwd_primer}...{rev_revcomp_primer} -g ^{rev_primer}...{fwd_revcomp_primer} " \
                   f"-n 2 --quiet --untrimmed-output={untrim_output} " \
                   f"-o {trim_output} {demux_filepath}"

    os.system(cutadapt_cmd)

def check_read():
    '''
    Checks to see if a read (string of sequences from a single read)
    contains the forward and/or reverse primer in their forward, reverse,
    and reverse complement orientation.
    :return: prints error if a primer is detected
    '''
    with open((output_path + file), 'rt') as fin:
        lines = fin.readlines()
        for l,line in enumerate(lines):
            if line.startswith('>'):
                continue
            else:
                read = line.strip()
                for p,name in zip(primers, primer_names):
                    result = re.search(p, read)
                    if result is not None:
                        global primer_detect_count
                        primer_detect_count += 1
                        with open(primer_check_summary, 'a') as fout:
                            fout.write(f"\nSAMPLE {sample_id}____________________________\n")
                            fout.write(f"The {name} ({p}) was found starting at position {result.span()[0]} in read {lines[l-1]}")

# INCOMPLETE FUNCTION
# def get_read_lens(fasta_path):
#     '''
#     Get a dictionary of the sample Id as a key and
#     the length of each read belonging to that sample
#     as the values, as a list.
#     :param fasta_path: path to fasta files to get the
#     read lengths from
#     :return: dictionary of sample IDs and read lengths
#     '''
#     sample_ids = get_sample_id(fasta_path)
#
#     output_dict = {x:[] for x in sample_ids}
#
#     #for sample in sample_ids:
#     return None

def vsearch_dereplicate(trimmed_file):
    sample_id = trimmed_file.split('_')[-1].split('.')[0]

    input = trimmed_path + trimmed_file
    output_derep = f"{derep_path}{trimmed_file.replace('trim','derep')}"
    output_summary = f"{derep_summary_path}{trimmed_file.replace('trim','derep-summary').replace('fasta','txt')}"
    vsearch_derep_cmd = "vsearch --derep_fulllength " + input +\
        " --output " + output_derep + " --sizeout " +\
        " --uc " + output_summary + " --quiet"
    print(f"\r    sample {sample_id} ({f + 1}/{len(trimmed_files)})--------", end='')
    os.system(vsearch_derep_cmd)

def get_full_derep_summary(file_path_list):

    sample_ids = []
    num_unique_reads = []
    num_total_reads = []
    mean_cluster_size = []
    std_cluster_size = []
    min_cluster_size = []
    max_cluster_size = []

    for file in file_path_list:
        summary1 = pd.read_csv(file, delimiter="\t",
                              names=['record_type','cluster','cluster_size',
                                     'percent_centroid_similarity','match_orientation',
                                     'NA1','NA2','NA3','query_seq','centroid_seq']).drop(['NA1','NA2','NA3'], axis = 1)
        summary = summary1[(summary1['record_type'] != 'C')]
        sample_ids.append(file.split("/")[2].strip(".txt").replace('derep-summary_',""))
        num_unique_reads.append(len(summary['cluster'].unique()))
        num_total_reads.append(summary.shape[0])
        mean_cluster_size.append(round(np.mean(summary['cluster_size']),2))
        std_cluster_size.append(round(np.std(summary['cluster_size']),2))
        min_cluster_size.append(np.min(summary['cluster_size']))
        max_cluster_size.append(np.max(summary['cluster_size']))

    summary_df = pd.DataFrame(list(zip(sample_ids, num_unique_reads, num_total_reads, mean_cluster_size, std_cluster_size,
                                       min_cluster_size, max_cluster_size)),
                              columns=['sample_id','num_unique_reads','num_total_reads','mean_cluster_size',
                                       'std_cluster_size','min_cluster_size','max_cluster_size'])

    summary_df['percent_duplicated'] = round(((1 - (summary_df['num_unique_reads'] / summary_df['num_total_reads']))*100),2)

    summary_df.to_csv(df_out, index = False)

    return None

def vsearch_chimeras(derep_file):
    sample_id = derep_file.split('_')[-1].split('.')[0]
    input = derep_path + derep_file
    output_nochim = f"{nonchim_path}{derep_file.replace('derep','nochim')}"
    output_chim = f"{chim_path}{derep_file.replace('derep','chim')}"
    output_summary = f"{chim_summary_path}{derep_file.replace('derep','uchime-summary').replace('fasta','txt')}"
    vsearch_chimera_cmd = f'vsearch --uchime_denovo {input} --chimeras {output_chim} ' \
                          f'--nonchimeras {output_nochim} --uchimeout {output_summary} ' \
                          f'--quiet'
    print(f"\r    sample {sample_id} ({f + 1}/{len(derep_files)})--------", end='')
    os.system(vsearch_chimera_cmd)

def get_size(header):
    size = header.split(";")[-1].replace("size=","")
    return int(size)

def vsearch_filter(nonchim_file):
    sample_id = nonchim_file.split('_')[-1].split('.')[0]
    input = nonchim_path + nonchim_file
    output_filter = f"{filtered_path}{nonchim_file.replace('nochim','filt')}"
    discarded_filter = f"{discard_path}{nonchim_file.replace('nochim','discard')}"
    vsearch_filter_cmd = f"vsearch --fastx_filter {input} --fastaout {output_filter} " \
                         f"--minsize {minsize_threshold} --quiet --fastaout_discarded {discarded_filter}"
    print(f"\r    sample {sample_id} ({f + 1}/{len(nonchim_files)})--------", end='')
    os.system(vsearch_filter_cmd)

def mothur_align_all(filtered_file):
    sample_id = filtered_file.split('_')[-1].split('.')[0]
    input = f"{filtered_path}{filtered_file}"  # if run from CLI, no need to escape '-' characters

    namefile_cmd = f"unique.seqs(fasta={input}, outputdir={count_output_path})"
    align_cmd = f"pairwise.seqs(fasta={input}, cutoff={cutoff}, countends=F, outputdir={distmat_output_path}, processors={num_processors})"
    align_log_cmd = f"set.logfile(name={distmat_log}, append=T)"
    names_log_cmd = f"set.logfile(name={count_log}, append=T)"

    print(f"\r    sample {sample_id} ({f + 1}/{len(filtered_files)})--------", end='')
    mothur_run_cmd = f"{mothur_path} '#{align_log_cmd};{align_cmd};{names_log_cmd};{namefile_cmd};' -q" # no matter what i do it won't be quiet
    os.system(mothur_run_cmd)

def mothur_cluster(distmat_file):
    sample_id = distmat_file.split('_')[-1].split('.')[0]
    distmat_input = f"{distmat_path}{distmat_file}"
    count_input = f"{count_path}{distmat_file.replace('dist','count_table')}"

    cluster_log_cmd = f"set.logfile(name={cluster_log}, append=T)"
    mothur_cluster_cmd = f"{mothur_path} '#{cluster_log_cmd};cluster(column={distmat_input}, " \
                         f"outputdir={cluster_path}, count={count_input}, cutoff={cutoff}, " \
                         f"method=average)'"
    print(f"\r    sample {sample_id} ({f + 1}/{len(distmat_files)})--------", end='')
    os.system(mothur_cluster_cmd)

def cast_lab_as_float(input_table):
    input_table02 = input_table[input_table['label'] != 'unique']
    input_table02['label'] = input_table02['label'].astype(dtype=float)
    return input_table02
