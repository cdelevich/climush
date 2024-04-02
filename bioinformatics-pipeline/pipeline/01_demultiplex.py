from mapping import filepath_map as fpm

import argparse
import sys
import re
from pathlib import Path
import json
import pandas as pd
##REMOVE AFTER PACKAGE TESTING#######
sys.path.insert(0, str(fpm['package']['main']))
#####################################
from climush.constants import *
from climush.utilities import *

settings = import_config_as_dict(fpm['config']['main'], file_handle='pipeline-settings', config_section='all')

# set variable that will only switch to True if files requiring demultiplexing are detected
files_demuxed = False

#####################
# ILLUMINA ##########
#####################

# check if there are Illumina reads that need to be demultiplexed
is_input, illumina_files = check_for_input(fpm['sequences']['demux'], seq_platform='illumina')

if is_input:
    msg = f'WARNING. Currently, there is no script that can demultiplex Illumina reads. You can continue with the ' \
          f'pipeline, but these {len(illumina_files)} multiplexed Illumina sequencing files will be ignored. Do you ' \
          f'wish to continue without these sequences?'
    prompt_yes_no_quit(message = msg)
else:
    pass

#####################
# PACBIO ############
#####################

# check if there are PacBio reads that need to be demultiplexed
is_input, pacbio_files = check_for_input(fpm['sequences']['demux'], seq_platform='\d{4}')

if is_input:
    files_demuxed = True
    print(f'{len(pacbio_files)} PacBio sequencing files were detected that require demultiplexing...')

    # create output directory for demultiplexed samples
    mkdir_exist_ok(new_dir = fpm['pipeline-output']['demultiplexed'])


else:
    pass

#####################
# SANGER ############
#####################

# check if there are Sanger reads that need to be demultiplexed
is_input, sanger_files = check_for_input(fpm['sequences']['demux'], seq_platform='sanger')

if is_input:
    msg = f'WARNING. Currently, there is no script that can demultiplex Sanger reads. You can continue with the ' \
          f'pipeline, but these {len(sanger_files)} multiplexed Sanger sequencing files will be ignored. Do you ' \
          f'wish to continue without these sequences?'
    prompt_yes_no_quit(message = msg)
else:
    pass



#########################

# when all are demultiplexed (if possible), continue to the next script
# if not files_demuxed:  # print if no demultiplexing was carried out
#     print(f'No sequencing files were detected in the /sequences directory that require demultiplexing.\n')
# continue_to_next(Path(__file__), settings)


#####################################################
# VARIABLES AND PATHS REQUIRED BY THE DEMUX FUNCTIONS
#####################################################

# LIKELY GOOD TO GO
# run_name = settings['run_details']['run_name']
# quality_score = str(settings['quality_filtering']['pacbio']['qscore'])
# output_path = fpm['pipeline-output']['demultiplexed']
# raw_read_path = fpm['sequences']['demux']
#
# # REQUIRE ATTENTION
# barcode_path = ''
# final_demux_path = ''
# primer_path = settings['primers']  # moved up to the PacBio section
# primer_path = "../miscell/primers/"

###################################
### UPDATE TO BIOINFO MODULE
###################################

def demultiplex(file_map, multiplexed_files, settings_dict, seq_platform='pacbio'):

    # REMOVE AFTER TESTING ############
    # DON'T FORGET TO UNCOMMENT THE RUN_SUBPROCESS LINE WHERE LIMA ACTUALLY RUNS!!! #################
    file_map = fpm.copy()
    multiplexed_files = pacbio_files
    settings_dict = settings.copy()
    seq_platform = 'pacbio'

    excel_df = next(file_map['config']['bc_mapping'].glob('pacbio*.xlsx'))
    csv_df = next(file_map['config']['bc_mapping'].glob('*.csv'))
    ###################################

    # DEFINE RELEVANT FILE PATHS ####################################################################
    # get the paths needed for demux
    mapping_dir = file_map['config']['bc_mapping']  # directory containing all mapping files
    mapping_files = list(mapping_dir.glob('*'))  # get mapping files, as list (will reuse + exhaust)

    # ACCESS USER SETTINGS FROM CONFIG ##############################################################
    run_name = settings['run_details']['run_name']  # name to use for this sequence run

    # CONSTRUCT REGEX CONSTANTS #####################################################################
    # column name regex; column names vary among mapping files
    FWD_COL_RE = '(fwd)|(forward)'
    REV_COL_RE = '(rev)|(reverse)'
    SAMPLE_ID_RE = '^sample'

    # CONFIRM THAT PLATFORM IS PACBIO ###############################################################
    # check that the sequence platform is set to pacbio; not suitable for others platforms right now
    if not seq_platform == 'pacbio':
        print(f'Currently, this function only accepts PacBio sequences. Not suited for demultiplexing Sanger or '
              f'Illumina sequences.\n')
        sys.exit()

    # LOCATE MAPPING FILE FOR EACH MULTIPLEXED SEQUENCE FILE #########################################
    # get all unique run queue IDs from files needing demux; will use later on
    queue_ids = set()

    # link multiplexed sequence files in needs_demux to their corresponding mapping file in config directory
    demux_mapping = {}  # key is mapping file path, value(s) is sequence file path in needs_demux dir

    for mp_file in multiplexed_files:

        # get the sequencing core's queue ID from multiplexed seq file and corresponding climush sample ID for this queue ID
        qid = re.search('^(\d{4})', mp_file.name).group(0)  # get queue ID from file name
        queue_ids.add(qid)  # add the qid to the set of all queue ids requiring demux
        cid = settings_dict['pacbio_demultiplexing']['multiplex'][qid]  # get climush - queue ID pairing from config
        qid_files = [f for f in mapping_files if
                     re.search(f'^{cid}', f.name)]  # all files in mapping files dir that match climush id

        if len(qid_files) == 1:  # should only be one mapping file; because only 1, will also be key in dict
            if qid_files[0] in demux_mapping.keys():  # if already in dict (e.g., Pool1/Pool2)
                demux_mapping[qid_files[0]].append(mp_file)  # update list of corresponding multiplexed files w/ str
            else:
                demux_mapping.update({qid_files[0]: [mp_file]})  # add queue ID and first multiplexed file as list
        elif len(qid_files) == 0:  # if there's no corresponding mapping file..
            print(f'No mapping file matching the multiplexed file {mp_file.name} was located in '
                  f'{mapping_dir.parent.name + "/" + mapping_dir.name}. Please check that a mapping file for collection '
                  f'{cid} is present in this directory, and that the CliMush collection ID ({cid}) is at the start of the '
                  f'file name (e.g., {cid + "_mapping-file.xlsx"}).\n')
        else:  # if multiple mapping files for a given queue ID
            print(f'{len(qid_files)} mapping files were located in {mapping_dir.parent.name + "/" + mapping_dir.name}. '
                  f'Please select the number corresponding to the correct file to use for the multiplexed sequencing '
                  f'file {cid}, or type \'quit\'.\n')
            right_map = prompt_print_options(qid_files)  # get right mapping file from user input
            if right_map in demux_mapping.keys():  # add to dict
                demux_mapping[right_map].append(mp_file)  # update list of corresponding multiplexed files w/ str
            else:
                demux_mapping.update({right_map: [mp_file]})  # add queue ID and first multiplexed file as list

    # CREATE A FASTA FILE OF ALL UNIQUE BARCODES USED TO MULTIPLEX ###################################
    # define subfunctions to help locate and prepare barcodes for demultiplexing
    def get_col_name(pattern, df):

        # search for columns that match the provided pattern
        matches = [c for c in df.columns if re.search(pattern, c, re.I)]

        # output depends on the number of matches returned (ideally 1)
        if len(matches) == 1:  # if one match, return that column name
            return matches[0]

        else:  # if multiple or no matches, rerun with new user-input regex

            # check number of matches, will print out different prompt, but otherwise executes same thing
            if len(matches) == 0:  # if no matches
                print(f'No column matching the pattern \'{pattern}\' was located in the input dataframe. Please try '
                      f'another regular expression that better matches a single column out of these columns in the '
                      f'dataframe: ')
            else:
                print(f'{len(matches)} column names match the pattern \'{pattern}\' in the input dataframe. Please try '
                      f'another regular expression that better matches a single column out of these columns in the '
                      f'dataframe: ')

            print_indented_list(matches)  # print out the columns of the dataframe (formatted to be indented)
            retry_pattern = input()
            return get_col_name(pattern=retry_pattern, df=df)
    def get_barcodes(map_path, primers):

        # create an empty list to add fwd/rev barcodes to; output from function then added to set from all mapping df
        fwd_barcodes = []
        rev_barcodes = []

        # read in mapping file based on file type (.xlsx, .csv, .txt)
        mapping_tabs = import_mapping_df(df_path=map_path)

        # go through each tab in the df and get the barcodes
        for tab in mapping_tabs:  # go through each tab...
            if re.search('pool', tab, re.I):  # ...only if relevant to mapping/demux

                # get the name used for the fwd/rev barcode columns in this mapping file
                fwd_col = get_col_name(FWD_COL_RE, mapping_tabs[tab])
                rev_col = get_col_name(REV_COL_RE, mapping_tabs[tab])

                # get list of values for barcodes and sample IDs
                fwd_bc = mapping_tabs[tab][fwd_col].to_list()
                rev_bc = mapping_tabs[tab][rev_col].to_list()

                # add all barcodes in the column to the total barcode list
                fwd_barcodes.extend(fwd_bc)  # add bc to the total list
                rev_barcodes.extend(rev_bc)  # use extend to prevent sublists

        # make sure that the barcodes don't contain the primer sequences; if they do, remove primer sequence from bc
        fwd_bc_noprimer = [re.sub(f'{primers["fwd_primer"]}$', '', fwd) for fwd in fwd_barcodes]
        rev_bc_noprimer = [re.sub(f'{primers["rev_primer"]}$', '', rev) for rev in rev_barcodes]

        # confirm that primers were properly removed from the barcodes
        def primers_removed_from_bc(pre_remove_list, post_remove_list, primer):

            # OBSERVED
            # get the length of each barcode after removing the primer; keep only unique length values
            len_bc_noprimer = list(set(len(f) for f in post_remove_list))

            # EXPECTED
            # get expected length of barcode w/o primer by subtracting primer length from original barcode length

            # get the length of each barcode prior to removing the primer; keep only unique length values
            len_bc_primer = list(set(len(f) for f in pre_remove_list))

            # confirm that all original barcodes (w/ primers) are the same length (should be), then calc expected
            if len(len_bc_primer) == 1:  # if all barcodes the same length prior to removing primer...
                len_bc_expected = len_bc_primer[0] - len(primer)  # ...calc single expected length
            else:  # throw error and return False if original barcodes not all same length
                print(f'ERROR. Barcodes in the mapping file should be the same length, but are not.')
                print(f'number of unique lengths of barcodes w/ primers: {len(len_bc_primer)}')
                print(f'unique lengths of barcodes w/ primers:           {len_bc_primer}\n')
                return False

            # COMPARISON
            # if all barcodes are the same length after removing primer, and they match the expected length, return True
            if (len(len_bc_noprimer) == 1) and (int(len_bc_noprimer[0]) == len_bc_expected):
                return True
            else:
                print(f'ERROR. The barcodes after removing the primer are either different lengths or do not match '
                      f'the expected length after removing the primer.')
                print(f'number of unique barcode lengths: {len(len_bc_noprimer)}')
                print(f'length of barcode w/ primer:      {len_bc_expected}')
                print(f'length(s) of barcode w/o primer:  {len_bc}\n')
                return False

        fwd_bc_ready = primers_removed_from_bc(pre_remove_list = fwd_barcodes,
                                               post_remove_list = fwd_bc_noprimer,
                                               primer = primers["fwd_primer"])
        rev_bc_ready = primers_removed_from_bc(pre_remove_list = rev_barcodes,
                                               post_remove_list = rev_bc_noprimer,
                                               primer = primers["rev_primer"])

        # if both fwd and rev barcodes are confirmed free from primers...
        if fwd_bc_ready and rev_bc_ready:
            bc_dict = {'fwd_bc': fwd_bc_noprimer,
                       'rev_bc': rev_bc_noprimer}
            return bc_dict  # export a dictionary of primer-free barcodes
        else:  # if either had an issue..
            print(f'ERROR. There was an issue removing primers from the barcode sequencing for the mapping file: '
                  f'{map_path.name}. Demultiplexing was therefore not completed. See below to determine whether the '
                  f'forward or reverse barcodes (or both) had this issue:')
            print(f'\tprimer properly removed from...')
            print(f'\t  forward barcodes: {fwd_bc_ready}')
            print(f'\t  reverse barcodes: {rev_bc_ready}')
            return sys.exit()  # exit

        # create dict of these fwd/rev barcode lists (easier to keep track of which is fwd/rev this way)

    # get the fwd/rev primers, since they are sometimes included in the barcode sequence and need to be removed
    primer_dict = {}
    primer_dict['fwd_primer'] = settings_dict['primers']['fwd']['sequence']['pacbio']
    primer_dict['rev_primer'] = settings_dict['primers']['rev']['sequence']['pacbio']

    # create a dict to add all barcodes to; will be used to create the unique barcode fasta
    unique_barcodes = {'fwd_bc': set(),
                       'rev_bc': set()}

    # aggregate all barcodes from the mapping files that belong to the multiplexed sequences that need demux
    for dxmap in demux_mapping:

        # get all barcodes for this seq run's mapping file, with primers removed
        seq_run_bc = get_barcodes(map_path = dxmap, primers = primer_dict)

        # add barcodes from this seq run to the overall dict of unique barcodes
        for b in ['fwd_bc', 'rev_bc']:

            # before adding, count the number of barcodes already included in the unique barcodes dict
            num_bc_overall = len(unique_barcodes[b])

            # count the number of unique barcodes that are in this current mapping file
            num_bc_mapping = len(set(seq_run_bc[b]))

            # the number of unique barcodes from this mapping file should match the number in the overall bc dict
            if num_bc_overall == 0:  # if previously empty (first addition to dict), nothing to compare
                pass  # continue to adding this mapping file's barcodes to the overall dict
            else:
                if num_bc_overall < num_bc_mapping:  # if there are additional barcodes in this mapping file...
                    extra_bc = set(seq_run_bc[b]).difference(unique_barcodes[b])  # find extra bc
                    print(f'WARNING. {len(extra_bc)} new barcode(s) added to the overall list of unique barcodes:\n'
                          f'\textra {b.split("_")[0]} barcode(s): {extra_bc}\n'
                          f'All mapping files should contain the same set of unique forward and reverse barcodes. '
                          f'Continuing with demultiplexing, but if there\'s a future error, see the mapping file'
                          f'that triggered this warning: {map.name}\n')
                elif num_bc_overall > num_bc_mapping:  # if not all possible barcodes are in this mapping file...
                    missing_bc = unique_barcodes[b].difference(set(seq_run_bc[b]))  # find missing bc
                    print(f'WARNING. {len(missing_bc)} barcode(s) were missing from this mapping file\'s barcodes:\n'
                          f'\tmissing {b.split("_")[0]} barcode(s): {missing_bc}\n'
                          f'All mapping files should contain the same set of unique forward and reverse barcodes. '
                          f'Continuing with demultiplexing, but if there\'s a future error, see the mapping file '
                          f'that triggered this warning: {map.name}\n')
                else:
                    pass  # they are equal, as expected

            # regardless of above outcome, add barcodes from this mapping file to the overall dict of unique barcodes
            unique_barcodes[b].update(seq_run_bc[b])  # use update to add list to set

    # write out fasta file of unique barcodes; must be ordered so same numbering is used every time this is run
    # create new path in pipeline output, and use this path to create fasta file name
    fasta_out_dir = mkdir_exist_ok(new_dir = 'barcodes', parent_dir = file_map['pipeline-output']['demultiplexed'])
    barcode_fasta = (fasta_out_dir / f'barcodes_{run_name}').with_suffix('.fasta')

    # create a dictionary with the barcode seq as the key and the value is the ID name (header) in the output fasta file
    unique_barcode_labels = {}
    bc_name_index = {}

    fasta_index = 0
    with open(barcode_fasta, 'wt') as fout:
        for bc_dir in unique_barcodes:  # for fwd/rev unique barcodes...

            bc_list = list(unique_barcodes[bc_dir])  # convert set to list, so it can be sorted
            bc_list.sort()  # sort barcode list, so barcode headers in fasta are numbered the same way every time

            for b,bc in enumerate(bc_list):  # for each unique barcode in the sorted list...

                bc_name = f'{bc_dir}_{b+1}'  # create barcode id w/ barcode number; +1 to account for zero-index

                # write barcode name and its sequence to the fasta file
                fout.write(f'>{bc_name}\n')  # add header
                fout.write(f'{bc}\n')  # add barcode sequence

                # add barcode name and its sequence to the unique barcode labels dictionary
                unique_barcode_labels.update({bc: bc_name})

                # add barcode name and its index to the bc name index dictionary
                bc_name_index.update({fasta_index: bc_name})
                fasta_index += 1  # have to use counter because pulling from two sources, will reset to 0 for rev bc

    # DEMULTIPLEX USING LIMA #############################################################################
    # create a directory for all lima output
    lima_output = mkdir_exist_ok(new_dir='lima_output', parent_dir=file_map['pipeline-output']['demultiplexed'])
    lima_subdirs = []  # create list of output directories to read lima output from later on

    for qid in queue_ids:

        # get the climush sequencing run ID that corresponds to this queue ID
        seq_run = settings_dict['pacbio_demultiplexing']['multiplex'][qid]

        # create an output directory for this sequencing run using the climush sequencing run ID
        seq_run_output = mkdir_exist_ok(new_dir=seq_run, parent_dir=lima_output) # create an output dir for each seq run
        lima_subdirs.append(seq_run_output)

        # located all multiplexed sequencing files that match this queue ID
        pools_to_demux = [f for f in multiplexed_files if re.search(f'^{qid}', f.name)]

        # alert if only one pool is detected
        if len(pools_to_demux) == 1:
            msg = f'WARNING. There is only one multiplexed sequencing file ({pools_to_demux[0].name}) for sequencing ' \
                  f'run {seq_run} when two are typically expected: one for each pool of sequences that were ' \
                  f'sequenced simultaneously on separate SMRT cells. Do you wish to continue with only half of the ' \
                  f'sequences from this sequencing run?'  # yes/no/quit automatically added by function
            prompt_yes_no_quit(msg)  # if yes, auto-continues; if no or quit, exits

        for pool in pools_to_demux:

            # detect pool number from the file name; remove queue ID first to avoid any confusion
            try:  # should be able to do automatically
                pool_num = re.search('\d', re.sub(f'^{qid}\.', '', pool.name)).group(0)
            except AttributeError:  # but if re.search does not return a match
                print(f'The sequencing pool number for the multiplexed file {pool.name} could not be detected in ' \
                      f'the file name. Please enter the pool number corresponding to this file as a single digit'
                      f'(i.e., no leading zeros or decimals): ')
                pool_num = input()

            # define prefix to use for lima output files; include path to the output directory for this seq run
            out_prefix = seq_run_output / f'lima-demux_{seq_run}_pool{pool_num}'

            # assemble list of commands required to run lima demultiplexing
            lima_cmd = ['lima', pool, barcode_fasta, out_prefix, '--min-score', '93', '--hifi-preset', 'ASYMMETRIC']

            # run_subprocess(cli_command_list = lima_cmd, dest_dir = lima_output)

    # CREATE DICTIONARY OF BARCODE PAIRS FOR ALL SAMPLES
    sample_barcodes = {q:{'pool1': {},
                          'pool2': {}} for q in queue_ids}

    for dxmap in demux_mapping:

        # get the queue ID from the name of the multiplexed sequencing file; need it to find key in sample_barcodes
        qid = re.search('^(\d{4})', list(demux_mapping[dxmap])[0].name).group(0)

        # import the mapping file dataframe
        mapping_df = import_mapping_df(df_path=dxmap)

        # loop through each tab (= pool1/pool2) in the df to get the sample ID and its barcode combo
        for tab in mapping_df:

            # get the name of the column used for the sample ID, fwd barcode, rev barcode
            smp_col = get_col_name(pattern=SAMPLE_ID_RE, df=mapping_df[tab])
            fwd_col = get_col_name(pattern=FWD_COL_RE, df=mapping_df[tab])
            rev_col = get_col_name(pattern=REV_COL_RE, df=mapping_df[tab])

            for i in range(mapping_df[tab].shape[0]):

                # from the dataframe, get the sample ID, fwd barcode sequence, and rev barcode sequence
                smp_id = mapping_df[tab][smp_col][i]
                fwd_full_bc = mapping_df[tab][fwd_col][i]  # pull full barcode from df (may contain primer still)
                rev_full_bc = mapping_df[tab][rev_col][i]  # may also container primer

                # get the barcode matching the full barcode from the df; keep the short (no-primer) version
                fwd_bc = [f for f in unique_barcode_labels if re.search(f'^{f}', fwd_full_bc, re.I)][0]
                fwd_id = unique_barcode_labels[fwd_bc]  # use non-primer barcode to get the barcode ID
                rev_bc = [r for r in unique_barcode_labels if re.search(f'^{r}', rev_full_bc, re.I)][0]
                rev_id = unique_barcode_labels[rev_bc]

                sample_barcodes[qid][tab].update({smp_id:{}})  # add sample key
                sample_barcodes[qid][tab][smp_id].update({fwd_id:fwd_bc, rev_id:rev_bc})  # add bc IDs and seqs

    # READ IN INFO FROM LIMA OUTPUT
    # set this dictionary up the same way it was done for sample_barcodes, to distinguish different seq runs, pools, etc.
    read_barcodes = {q:{'pool1': {},
                        'pool2': {}} for q in queue_ids}

    for subdir in lima_subdirs:

        # get the queue ID from the name of the multiplexed sequencing file; need it to find key in sample_barcodes
        cid = re.search('^pacbio.+', subdir.name).group(0)  # get climush ID, which will be in name of file
        mp_dict = settings_dict['pacbio_demultiplexing']['multiplex']  # make shorter for list comp below
        qid = [k for k in mp_dict if cid in mp_dict[k]][0]  # convert to qid to match dict

        for p in ['pool1', 'pool2']:

            # find the fasta file for this pool, which has the read ID and the barcode combination
            lima_fasta = [f for f in subdir.glob('*.*') if re.search(f'{subdir.name}_{p}\.fasta', f.name, re.I)][0]

            # open the fasta file and pull the indices of the forward and reverse barcodes
            with open(lima_fasta, 'r') as fin:
                for r in fin.readlines():
                    if r.startswith('>'):
                        read_id = re.findall('(?<=>).+?(?=\sbc)', r, re.I)[0]
                        read_barcode_i = list(map(int, re.search('(?<=bc=)\d{1,2},\d{1,2}', r, re.I).group(0).split(',')))
                        read_bc_names = [bc_name_index[i] for i in read_barcode_i]  # convert index to name
                        read_barcodes[qid][p].update({read_id: read_bc_names})
                    else:
                        continue

    # COMPARE LIMA BARCODE COMBINATIONS TO MAPPING FILE BARCODE COMBINATIONS
    final_demux = {q:{'pool1': {},
                      'pool2': {}} for q in queue_ids}

    for qid in read_barcodes:
        for p in ['pool1', 'pool2']:
            sample_subdict = sample_barcodes[qid][p]
            for r in read_barcodes[qid][p]:
                # REMOVE AFTER TESTING ##############
                r = 'm64047_230523_213434/197456/ccs'
                #####################################

                read_id = r
                bc_list = read_barcodes[qid][p][r]

                # test = [list(sample_subdict[s]) for s in sample_subdict]
                sample_id = [s for s in sample_subdict if len(set(sample_subdict[s].keys()).difference(set(bc_list))) == 0]

                if len(sample_id) == 1:
                    read_seq =
                    final_demux[qid][p].update({sample_id[0]:})

    def get_sample_id():





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



