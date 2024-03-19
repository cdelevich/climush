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
run_name = settings['run_details']['run_name']
quality_score = str(settings['quality_filtering']['pacbio']['qscore'])
output_path = fpm['pipeline-output']['demultiplexed']
raw_read_path = fpm['sequences']['demux']

# REQUIRE ATTENTION
barcode_path = ''
final_demux_path = ''
# primer_path = settings['primers']  # moved up to the PacBio section
# primer_path = "../miscell/primers/"

###################################
### UPDATE TO BIOINFO MODULE
###################################

def demultiplex(file_map, multiplexed_files, settings_dict, seq_platform='pacbio'):

    # REMOVE AFTER TESTING ############
    file_map = fpm.copy()
    multiplexed_files = pacbio_files
    settings_dict = settings.copy()
    seq_platform = 'pacbio'

    excel_df = next(file_map['config']['bc_mapping'].glob('pacbio*.xlsx'))
    csv_df = next(file_map['config']['bc_mapping'].glob('*.csv'))
    ###################################

    # DEFINE RELEVANT FILE PATHS
    # get the paths needed for demux
    mapping_dir = file_map['config']['bc_mapping']  # directory containing all mapping files
    mapping_files = list(mapping_dir.glob('*'))  # get mapping files, as list (will reuse + exhaust)

    # CONFIRM THAT PLATFORM IS PACBIO
    # check that the sequence platform is set to pacbio; not suitable for others platforms right now
    if not seq_platform == 'pacbio':
        print(f'Currently, this function only accepts PacBio sequences. Not suited for demultiplexing Sanger or '
              f'Illumina sequences.\n')
        sys.exit()

    # LOCATE MAPPING FILE FOR EACH MULTIPLEXED SEQUENCE FILE
    # link multiplexed sequence files in needs_demux to their corresponding mapping file in config directory
    demux_mapping = {}  # key is mapping file path, value(s) is sequence file path in needs_demux dir
    for mp_file in multiplexed_files:

        # get the sequencing core's queue ID from multiplexed seq file and corresponding climush sample ID for this queue ID
        qid = re.search('^(\d{4})', mp_file.name).group(0)  # get queue ID from file name
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

    # CREATE A FASTA FILE OF ALL UNIQUE BARCODES USED TO MULTIPLEX
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

        # column name regex; column names vary among mapping files
        FWD_COL_RE = '(fwd)|(forward)'
        REV_COL_RE = '(rev)|(reverse)'

        # create an empty list to add fwd/rev barcodes to; output from function then added to set from all mapping df
        fwd_barcodes = []
        rev_barcodes = []

        # read in mapping file based on file type (.xlsx, .csv, .txt)
        if re.search('^\.x', map_path.suffix):  # if an excel file
            mapping_tabs = pd.read_excel(map_path, sheet_name=None)  # need to set sheet_name to None to get all tabs
        elif re.search('^\.c|^\.txt$', map_path.suffix):  # I think you can read in .txt and .csv files the same way?
            mapping_tabs = {'pool#': pd.read_csv(map_path)}  # make same format as xlsx dict so next part works for all
        else:
            print(f'ERROR. The file format {map_path.suffix} of the mapping file {map_path.name} is not a recognized '
                  f'file type. Accepted file types are: \'.xlsx\', \'.csv\', and \'.txt\'.\n')
            sys.exit()

        # go through each tab in the df and get the barcodes
        for tab in mapping_tabs:  # go through each tab...
            if re.search('pool', tab, re.I):  # ...only if relevant to mapping/demux

                fwd_col = get_col_name(FWD_COL_RE, mapping_tabs[tab])  # get name of bc column used in df
                rev_col = get_col_name(REV_COL_RE, mapping_tabs[tab])

                fwd_barcodes.extend(mapping_tabs[tab][fwd_col].to_list())  # add bc to the total list
                rev_barcodes.extend(mapping_tabs[tab][rev_col].to_list())  # use extend to prevent sublists

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
    for map in demux_mapping:

        # get all barcodes for this seq run's mapping file, with primers removed
        seq_run_bc = get_barcodes(map_path = map, primers = primer_dict)

        # add barcodes from this seq run to the overall dict of unique barcodes
        for b in ['fwd_bc', 'rev_bc']:

            # before adding, count the number of barcodes already included in the unique barcodes dict
            num_bc_overall = len(unique_barcodes[b])

            # count the number of unique barcodes that are in this current mapping file
            num_bc_mapping = len(set(seq_run_bc[b]))

            # the number of unique barcodes from this mapping file should match the number in the overall bc dict
            if num_bc_before == 0:  # if previously empty (first addition to dict), nothing to compare
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




###################################
### FROM 01_CREATE-BARCODE-FASTA.PY
###################################


# define single function to carry out barcode file creation
def create_barcode_fasta(barcode_df, pool_num, settings_dict, file_map):
    """
    Create a fasta file of unique forward and reverse barcodes from an input mapping file that
    should contain a column for the forward barcode, reverse barcode, and sample ID. Column naming
    is flexible, as function will look for any column name with 'fwd' or 'forward', 'rev' or
    'reverse', and 'sample' and 'id', 'no', 'num', with or without space or underscore, and case
    insensitive. Barcodes may have primers attached in these tables, which will be removed, if present,
    in this function.
    :param barcode_df: input barcode file already read in
    :param pool_num: if sequencing run was split into multiple pools, include the pool number
    :return: no direct output, but will save a [1] fasta file of the barcodes and [2] a csv of
    sample IDs and their barcode combinations for in further steps of demultiplexing.
    """
    # REMOVE AFTER TESTING ###############
    barcode_df = pool1_df.copy()
    pool_num = 1
    settings_dict = settings.copy()
    ######################################

    # create a dict to collect info for fwd and rev barcodes
    bc_info = {'forward':{'col_name_regex': '(fwd)|(forward)',
                          'primer_seq': '',  # get later when accessing primer seq from settings/config
                          'w_primer_col': 'fwd_w_bc',
                          'col_name': 'forward'},
               'reverse':{'col_name_regex': '(rev)|(reverse)',
                          'primer_seq': '',
                          'w_primer_col':'rev_w_bc',
                          'col_name': 'reverse'}}

    # get column names for barcodes and sample IDs
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

    fwd_col = get_col_name(pattern = bc_info['forward']['col_name_regex'], df = barcode_df)
    rev_col = get_col_name(pattern = bc_info['reverse']['col_name_regex'], df = barcode_df)
    smp_col = get_col_name(pattern = '(sample)(_|\s)(id|no|num)', df = barcode_df)

    # rename columns
    barcode_df.rename(columns={smp_col: 'sample_id',
                               fwd_col: bc_info['forward']['w_primer_col'],
                               rev_col: bc_info['reverse']['w_primer_col']},
                      inplace=True)

    # remove primers from barcodes, if present
    # get primer str from configuration file
    primer_fwd = settings_dict['primers']['fwd']['sequence']['pacbio']
    primer_rev = settings_dict['primers']['rev']['sequence']['pacbio']

    # add primer seqs to the bc_info dictionary
    bc_info['forward']['primer_seq'] = primer_fwd
    bc_info['reverse']['primer_seq'] = primer_rev

    # create new column with the primers removed from the barcodes
    barcode_df['forward'] = [re.sub(f'{primer_fwd}$', '', fwd) for fwd in barcode_df[bc_info['forward']['w_primer_col']]]
    barcode_df['reverse'] = [re.sub(f'{primer_rev}$', '', rev) for rev in barcode_df[bc_info['reverse']['w_primer_col']]]

    # confirm that primers removed from the barcodes, then remove old columns w/ primers still attached
    def primers_removed_from_bc(df, bc_direction, barcode_info):

        # get the length of the barcode without the primer
        len_bc = list(set(len(f) for f in df[barcode_info[bc_direction]['col_name']]))

        # get expected length of barcode w/o primer by subtracting primer length from original barcode length
        len_og_bc = list(set(len(f) for f in df[barcode_info[bc_direction]['w_primer_col']]))
        if len(len_og_bc) == 1:
            expected_len = len_og_bc[0] - len(barcode_info[bc_direction]['primer_seq'])
        else:
            print(f'ERROR.')
            print(f'number of unique lengths of barcodes w/ primers: {len(len_og_bc)}')
            print(f'unique lengths of barcodes w/ primers:            {len_og_bc}')
            return False

        if (len(len_bc) == 1) and (int(len_bc[0]) == expected_len):  # all same length and length is size of bc w/ primer minus len primer
            return True
        else:
            print(f'ERROR.')
            print(f'number of unique barcode lengths: {len(len_bc)}')
            print(f'length of barcode w/ primer:      {expected_len}')
            print(f'length(s) of barcode w/o primer:  {len_bc}')
            return False

    for d in bc_info.keys():
        if not primers_removed_from_bc(df = barcode_df, bc_direction = d, barcode_info = bc_info):
            print(f'Something went wrong with removing the primers from the barcodes in the '
                  f'mapping file. Exiting...\n')  # NEED MORE INFORMATIVE PRINT STATEMENT HERE
            # sys.exit()

    # create df that contains the sample IDs with the sequences of the barcodes
    # create dictionary with sample IDs as keys, and values for forward and reverse barcodes for each
    fwd_series = barcode_df[bc_info['forward']['col_name']]
    rev_series = barcode_df[bc_info['reverse']['col_name']]

    barcode_dict = {k:{'fwd_bc':'', 'rev_bc':''} for k in barcode_df['sample_id']}

    for smp in barcode_dict:
        barcode_dict[smp]['fwd_bc'] = barcode_df[barcode_df['sample_id'] == smp][bc_info['forward']['col_name']].reset_index(drop=True)[0]
        barcode_dict[smp]['rev_bc'] = barcode_df[barcode_df['sample_id'] == smp][bc_info['reverse']['col_name']].reset_index(drop=True)[0]

    # wanted_cols = ['sample_id','forward','reverse']
    # sampleid_df = pd.melt(barcode_csv[wanted_cols], id_vars='sample_id', var_name='barcode_direction', value_name='sequence')
    # assert (len(barcode_csv['sample_id']) == len(sampleid_df['sample_id']) / 2), 'Dataframe likely contains extra columns that' \
    #                                                                              'were not be removed from the dataframe (though ' \
    #                                                                              'the code should do this).'

    # create a set of forward and reverse barcodes
    fwd_unique_bc = list(set(barcode_df[bc_info['forward']['col_name']]))
    rev_unique_bc = list(set(barcode_df[bc_info['reverse']['col_name']]))

    # get the number of unique barcodes, add to bc_info, and the number of unique samples
    bc_info['forward']['num_unique_bc'] = len(fwd_unique_bc)
    bc_info['reverse']['num_unique_bc'] = len(rev_unique_bc)
    num_unique_samples = len(barcode_df[smp_col].unique())

    # barcode_df = sampleid_df.drop('sample_id', axis=1).drop_duplicates(subset='sequence').sort_values(
    #     by='sequence').reset_index().drop('index', axis=1)
    #
    # # make sure each sequence is unique in df
    # assert len(fwd_bc) + len(rev_bc) == barcode_df.shape[0]
    # assert len(barcode_df['sequence'].unique()) == len(barcode_df['sequence'])

    # add sorted list of unique barcodes to the bc_info dict
    fwd_unique_bc.sort()  # must sort prior to adding to dict; lima must receive barcodes in same order each time
    bc_info['forward']['unique_bc'] = fwd_unique_bc
    rev_unique_bc.sort()  # must sort prior to adding to dict; lima must receive barcodes in same order each time
    bc_info['reverse']['unique_bc'] = rev_unique_bc

    # # create headers for barcode fasta file using the indices of the *sorted* dataframe, barcode_df
    # # must be sorted to produce the same file each time, which is crucial for the way Lima works
    # barcode_df['barcode_header'] = None
    #
    # for i in range(barcode_df.shape[0]):
    #     if barcode_df['barcode_direction'][i] == 'forward':
    #         barcode_df['barcode_header'][i] = "bc_fwd_" + str(barcode_df.index[i])
    #     else:
    #         barcode_df['barcode_header'][i] = "bc_rev_" + str(barcode_df.index[i])

    fasta_out_dir = mkdir_exist_ok(new_dir = 'barcodes', parent_dir = fpm['pipeline-output']['demultiplexed'])
    fasta_out = (fasta_out_dir / f'pool0{str(pool_num)}_barcodes').with_suffix('.fasta')
    with open(fasta_out, 'wt') as fout:

        for bc_dir in bc_info:  # for fwd/rev

            for b,bc in enumerate(bc_info[bc_dir]['unique_bc']):  # for each unique bc

                if bc_dir == 'forward':
                    fout.write(f'>bc_fwd_{b}\n')  # add header
                else:
                    fout.write(f'>bc_rev_{b}\n')

                fout.write(f'{bc}\n')  # add sequence

    # write out the barcode info dictionary to JSON
    bc_info_out = (fasta_out_dir / f'pool0{str(pool_num)}_barcode-summary').with_suffix('.json')

    # # write out barcodes and headers as fasta file
    # fasta_output_name = output_path + run_name + "pool" + str(pool_num) + "_barcodes.fasta"
    # with open(fasta_output_name, 'wt') as fout:
    #     for i in range(barcode_df.shape[0]):
    #         fout.write(">" + barcode_df['barcode_header'][i] + '\n')
    #         fout.write(barcode_df['sequence'][i] + '\n')

    # if os.path.isfile(fasta_output_name) == True:
    #     print(f"\nBarcode fasta file for pool {pool_num}, {os.path.basename(fasta_output_name)}, has been created.\n")

    # merge sample id dataframe and barcode id dataframe by sequence
    # merged_df = pd.merge(barcode_df, sampleid_df)
    # wider_df = pd.pivot(data=merged_df, index='sample_id', columns='barcode_direction', values='barcode_header').reset_index()
    # wider_df['fwd_index'] = [wider_df['forward'][i].split('_')[2] for i in range(wider_df.shape[0])]
    # wider_df['rev_index'] = [wider_df['reverse'][i].split('_')[2] for i in range(wider_df.shape[0])]
    #
    # csv_output_name = output_path + run_name + "pool" + str(pool_num) + "_merged.csv"
    # wider_df.to_csv(csv_output_name, index=False)

    # assert len(wider_df['forward'].unique()) == num_unique_fwd
    # assert len(wider_df['reverse'].unique()) == num_unique_rev
    # assert len(wider_df.index.unique()) == num_unique_samples

    # if os.path.isfile(csv_output_name) == True:
    #     print(f"Merged sample ID and barcode name file for pool {pool_num}, {csv_output_name}, has been created.\n")

    return None

# run command for however many pools you might have (can loop if several)
create_barcode_fasta(barcode_csv=pool1_df, pool_num=1, settings_dict=settings, file_map=fpm)
create_barcode_fasta(barcode_csv=pool2_df, pool_num=2, settings_dict=settings, file_map=fpm)
sys.exit()
###################################
### FROM 02_RUN-LIMA.PY
###################################

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



