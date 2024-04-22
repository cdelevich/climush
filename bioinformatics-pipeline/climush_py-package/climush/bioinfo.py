from pathlib import Path
from Bio.Seq import Seq
from datetime import datetime
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from climush.constants import NOPHIX_PREFIX, SEQ_FILE_GLOB, TRIMMED_PREFIX, DEREP_PREFIX, QUALFILT_PREFIX
from climush.utilities import *

def demultiplex(file_map, multiplexed_files, settings_dict, seq_platform='pacbio'):

    # REMOVE AFTER TESTING ############
    # DON'T FORGET TO UNCOMMENT THE RUN_SUBPROCESS LINE WHERE LIMA ACTUALLY RUNS!!! #################
    # file_map = fpm.copy()
    # multiplexed_files = pacbio_files
    # settings_dict = settings.copy()
    # seq_platform = 'pacbio'
    #
    # excel_df = next(file_map['config']['bc_mapping'].glob('pacbio*.xlsx'))
    # csv_df = next(file_map['config']['bc_mapping'].glob('*.csv'))
    ###################################

    # DEFINE RELEVANT FILE PATHS ####################################################################
    # get the paths needed for demux
    mapping_dir = file_map['config']['bc_mapping']  # directory containing all mapping files
    mapping_files = list(mapping_dir.glob('*'))  # get mapping files, as list (will reuse + exhaust)

    # ACCESS USER SETTINGS FROM CONFIG ##############################################################
    run_name = settings_dict['run_details']['run_name']  # name to use for this bioinformatics run

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
    # create a directory for the output files that lima produces
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

            run_subprocess(cli_command_list = lima_cmd, dest_dir = lima_output)

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

    # SORT INFO FROM LIMA OUTPUT
    # set these dictionaries up the same way it was done for sample_barcodes, to distinguish different seq runs, pools, etc.
    # dictionary for read barcode combinations
    read_barcodes = {q:{'pool1': {},
                        'pool2': {}} for q in queue_ids}
    # dictionary for read sequences w/ barcodes removed
    read_seqs = {q:{'pool1': {},
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
                fasta_txt = fin.readlines()
                for l, line in enumerate(fasta_txt):
                    if line.startswith('>'):
                        read_id = re.findall('(?<=>).+?(?=\sbc)', line, re.I)[0]
                        read_barcode_i = list(map(int, re.search('(?<=bc=)\d{1,2},\d{1,2}', line, re.I).group(0).split(',')))
                        read_bc_names = [bc_name_index[i] for i in read_barcode_i]  # convert index to name
                        read_barcodes[qid][p].update({read_id: read_bc_names})  # add read's bc names to dict
                        read_seqs[qid][p].update({read_id: fasta_txt[l+1]})  # add read sequence to dict
                    else:
                        continue

    # COMPARE LIMA BARCODE COMBINATIONS TO MAPPING FILE BARCODE COMBINATIONS

    # create an output subdirectory in demultiplexed dir for all demultiplexed reads in this bioinformatics run
    demux_dir = mkdir_exist_ok(new_dir = f'demux_{run_name}', parent_dir=file_map['pipeline-output']['demultiplexed'])

    # for each sequencing run...
    for qid in read_barcodes:

        # get the name of the sequencing run for labeling the output and read headers
        seq_run = settings['pacbio_demultiplexing']['multiplex'][qid]

        # for each pool in this sequencing run...
        for p in ['pool1', 'pool2']:

            # subset dictionary for this seq run/pool, makes more succinct code below
            sample_subdict = sample_barcodes[qid][p]

            # subset dictionary for this seq run/pool of read ID and barcode-free sequences from lima fasta output
            seq_subdict = read_seqs[qid][p]

            # for each read in this sequencing run's pool...
            for r in read_barcodes[qid][p]:

                # get a list of the barcodes found by lima in this read
                bc_list = read_barcodes[qid][p][r]

                # look for the sample ID from the mapping file that has the same pair of barcodes
                sample_ids = [s for s in sample_subdict if len(set(sample_subdict[s].keys()).difference(set(bc_list))) == 0]

                # confirm that the read's barcode combination matches only one sample ID, as it should
                if len(sample_ids) == 1:

                    # get sample ID as string, not list (now that it is confirmed there's only one)
                    sample_id = sample_ids[0]

                    # create read ID header from sequencing run, sample ID, and last part of original read header
                    full_sample_id = '_'.join([seq_run, sample_id])
                    read_id = '_'.join([full_sample_id, r.split('/')[-2]])

                    # get the barcode-free sequence for this read
                    read_seq = seq_subdict[r]

                    # write the read ID and read sequence to the sample ID demultiplexed fasta file
                    sample_fasta = (demux_dir / f'demux_{full_sample_id}').with_suffix('.fasta')
                    with open(sample_fasta, 'a') as fout:
                        fout.write(f'>{read_id}\n')
                        fout.write(f'{read_seq}\n')

                else:
                    print(f'ERROR. The barcode combination detected in read {r} from {p} of {seq_run} matched to '
                          f'{len(sample_ids)} sample IDs, when it should only match to one. This read was not '
                          f'sorted, and this error was recorded to the error table.\n')
                    barcode_error = (file_map['pipeline-output']['demultiplexed'] / 'barcode_errors').with_suffix('.tsv')

                    # if the error file doesn't yet exist, start with writing the header
                    if not barcode_error.is_file():
                        with open(barcode_error, 'a') as fout:
                            fout.write(f'sequencing_run\tsequencing_pool\tsample_ids\tread_id\tbarcode_combination\n')

                    # write out details of the read that triggered this error
                    with open(barcode_error, 'a') as fout:
                        fout.write(f'{seq_run}\t{p}\t{", ".join(sample_ids)}\t{r}\t{", ".join(bc_list)}\n')

def filter_out_phix(input_files, file_map, kmer=31, hdist=1, keep_log=True, keep_removed_seqs=True):

    prefilt_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output']['prefiltered'])

    run_name = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings',
                                     config_section='run_details')['run_name']

    nophix_path = mkdir_exist_ok(new_dir=f'./{NOPHIX_PREFIX}_{run_name}', parent_dir=prefilt_parent)
    phix_path = mkdir_exist_ok(new_dir=f'./{flip_prefix(NOPHIX_PREFIX)}_{run_name}', parent_dir=prefilt_parent)

    bbduk_log = make_log_file(file_name=Path('bbduk'), dest_path=prefilt_parent)

    for file in input_files:

        nophix_out = add_prefix(file_path=file, prefix=NOPHIX_PREFIX, action=None, dest_dir=nophix_path)
        phix_out = add_prefix(file_path=file, prefix=flip_prefix(NOPHIX_PREFIX), action=None, dest_dir=phix_path)

        # bbduk does not like special characters in file path
        # nophix_out_str = escape_path_special(nophix_out)
        # phix_out_str = escape_path_special(phix_out)
        # and this didn't help it work either...

        bbduk_cmd = ['bbduk.sh', f'in={file}', f'out={nophix_out}', f'out={phix_out}',
                     'ref=phix.fa', f'k={kmer}', f'hdist={hdist}', f'stats={bbduk_log}']

        try:
            run_subprocess(bbduk_cmd, dest_dir=prefilt_parent)
        except:
            continue


    if not check_for_input(nophix_path)[0]:  # this is a temporary work around to make copies that look like output while there's a curernt issue with Java
        for file in input_files:
            add_prefix(file_path=file, prefix=NOPHIX_PREFIX, action='copy', dest_dir=nophix_path)

    return nophix_path

def pair_reads(input_files):
    pairs_dict = {}
    rev_reads = []
    for file in input_files:
        if re.search('R1', file.stem, re.I):
            pairs_dict.update({file:''})
        else:
            rev_reads.append(file)

    for rev in rev_reads:
        sample_id = re.search('.+?(?=_R2)', rev.stem).group(0)
        for k in pairs_dict.keys():
            if re.search(sample_id, k.stem, re.I):
                pairs_dict[k] = rev

    return pairs_dict

def prefilter_fastx(input_files, file_map, maxn=0):
    prefilt_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output']['prefiltered'])

    run_name = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings',
                                     config_section='run_details')['run_name']

    noambig_path = mkdir_exist_ok(new_dir=f'./{NOAMBIG_PREFIX}_{run_name}', parent_dir=prefilt_parent)
    ambig_path = mkdir_exist_ok(new_dir=f'./{flip_prefix(NOAMBIG_PREFIX)}_{run_name}', parent_dir=prefilt_parent)

    fastq_maxns_log = make_log_file(file_name=Path('fastq_maxns'), dest_path=prefilt_parent)

    # must provide pair of forward and reverse:
    pairs_dict = pair_reads(input_files)

    for file in pairs_dict.keys():

        noambig_out_fwd = add_prefix(file_path=file, prefix=NOAMBIG_PREFIX, action=None, dest_dir=noambig_path)
        noambig_out_rev = add_prefix(file_path=pairs_dict[file], prefix=NOAMBIG_PREFIX, action=None,
                                     dest_dir=noambig_path)
        ambig_out_fwd = add_prefix(file_path=file, prefix=flip_prefix(NOAMBIG_PREFIX), action=None,
                                   dest_dir=ambig_path)
        ambig_out_rev = add_prefix(file_path=pairs_dict[file], prefix=flip_prefix(NOAMBIG_PREFIX), action=None,
                                   dest_dir=ambig_path)

        vsearch_cmd = ['vsearch', '--fastq_filter', file, '--reverse', pairs_dict[file], '--fastq_maxns', str(maxn),
                       '--fastqout', noambig_out_fwd, '--fastqout_rev', noambig_out_rev, '--fastqout_discarded',
                       ambig_out_fwd, '--fastqout_discarded_rev', ambig_out_rev]

        run_subprocess(vsearch_cmd, dest_dir=prefilt_parent)

    return noambig_path

def identify_primers(platform, config_dict):
    primer_dict = config_dict['primers']

    target_primers = {}
    primer_names = []
    for d in primer_dict.keys():
        target_primers[d] = primer_dict[d]['sequence'][platform]
        primer_names.append(primer_dict[d]['name'][platform])

    print(f'Searching for {primer_names[0]} and {primer_names[1]} in {platform.title()} reads...\n')

    return target_primers

def remove_primers(input_files, file_map, platform, paired_end=True):
    trim_primers_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output']['primers-trimmed'])

    settings = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings')
    run_name = settings['run_details']['run_name']

    # get the forward and reverse primers
    primer_dict = identify_primers(platform, config_dict=settings)
    fwd_primer = primer_dict['fwd']
    rev_primer = primer_dict['rev']

    fwd_revcomp_primer = str(Seq(fwd_primer).reverse_complement())
    rev_revcomp_primer = str(Seq(rev_primer).reverse_complement())

    trim_path = mkdir_exist_ok(new_dir=f'./{TRIMMED_PREFIX}_{run_name}', parent_dir=trim_primers_parent)
    notrim_path = mkdir_exist_ok(new_dir=f'./{flip_prefix(TRIMMED_PREFIX)}_{run_name}', parent_dir=trim_primers_parent)

    cutadapt_settings = settings['cutadapt']
    fwd_max_err = int(cutadapt_settings['max_error_rate']['fwd'])
    rev_max_err = int(cutadapt_settings['max_error_rate']['rev'])
    if fwd_max_err == rev_max_err:
        max_err = fwd_max_err
    max_untrimmed = int(cutadapt_settings['max_untrimmed'])

    if paired_end:
        paired_dict = pair_reads(input_files)

        for file in paired_dict.keys():

            fwd_read_out = add_prefix(file_path=file, prefix=TRIMMED_PREFIX, dest_dir=trim_path, action=None)
            rev_read_out = add_prefix(file_path=paired_dict[file], prefix=TRIMMED_PREFIX, dest_dir=trim_path, action=None)

            cutadapt_cmd = ['cutadapt', '--report=minimal',
                            '-a', fwd_primer, '-A', rev_primer,
                            '-n', '2', '-e', str(max_err),
                            '-o', fwd_read_out, '-p', rev_read_out,
                            file, paired_dict[file]]

            run_subprocess(cutadapt_cmd, dest_dir=trim_primers_parent)

    else:

        for file in input_files:

            trim_out = add_prefix(file_path=file, prefix=TRIMMED_PREFIX, dest_dir=trim_path, action=None)

            if cutadapt_settings['keep_untrimmed']:
                notrim_out = add_prefix(file_path=file, prefix=flip_prefix(TRIMMED_PREFIX),
                                        dest_dir=notrim_path, action=None)
                untrim_cmd = f'--untrimmed-output={notrim_out}'
            else:
                untrim_cmd = ''


            cutadapt_cmd = ['cutadapt', '--report=minimal',
                            '-a', f'^{fwd_primer}...{rev_revcomp_primer}$',
                            '-a', f'^{rev_primer}...{fwd_revcomp_primer}$',
                            '-n', '2', '-e', str(max_err), untrim_cmd, '-o', trim_out, file]

            # create or append to cutadapt stderr/stdout, all output goes to one file (not one per sample)
            run_subprocess(cutadapt_cmd, dest_dir=trim_primers_parent)

    # quantify proportion untrimmed
    # with open((trim_primers_parent / 'cutadapt.out'), 'r') as fin:
    #     cutadapt_df = pd.read_table(fin)
    #     sum_in = cutadapt_df['in_reads'].sum(0)
    #     sum_out = cutadapt_df['out_reads'].sum(0)
    #     percent_lost = ((sum_in - sum_out) / (sum_in)) * 100
    #     if percent_lost > max_untrimmed:
    #         msg = f'After primer trimming, {percent_lost:.2f}% of the input reads were lost, which is ' \
    #               f'above the user-defined maximum threshold of {max_untrimmed}%.\n'
    #         exit_process(message=msg)
    #     else:
    #         print(f'{percent_lost:.2f}% of input reads were lost to primer trimming. This is above the user-provided '
    #               f'input of {max_untrimmed}%, so proceeding to next step...\n')

    return None

def merge_reads(input_files, file_map):
    settings = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings')
    run_name = settings['run_details']['run_name']

    if settings['quality_filtering']['illumina']['merge_reads']:

        qfilt_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output']['quality-filtered'])
        merge_path = mkdir_exist_ok(new_dir=f'./{MERGED_PREFIX}_{run_name}', parent_dir=qfilt_parent)
        nomerge_path = mkdir_exist_ok(new_dir=f'./{flip_prefix(MERGED_PREFIX)}_{run_name}', parent_dir=qfilt_parent)

        merge_summary = qfilt_parent / 'vsearch_merged.log'

        paired_dict = pair_reads(input_files)

        output_list = []
        for file in paired_dict.keys():

            nomerge_fwd_out = add_prefix(file_path=file, prefix=MERGED_PREFIX,
                                         dest_dir=nomerge_path, action=None)
            nomerge_rev_out = add_prefix(file_path=paired_dict[file], prefix=flip_prefix(MERGED_PREFIX),
                                         dest_dir=nomerge_path, action=None)

            merge_out_base = add_prefix(file_path=file, prefix=MERGED_PREFIX, dest_dir=merge_path, action=None)
            sample_id = re.search('.+?(?=_R1)', merge_out_base.stem).group(0)
            merge_output = (merge_out_base.parent / sample_id).with_suffix('.fastq')
            output_list.append(merge_output)

            vsearch_merge_cmd = ['vsearch', '--fastq_mergepairs', file, '--reverse', paired_dict[file],
                                 '--fastqout', merge_output,
                                 '--fastqout_notmerged_fwd', nomerge_fwd_out,
                                 '--fastqout_notmerged_rev', nomerge_rev_out,
                                 '--eetabbedout', merge_summary]

            run_subprocess(vsearch_merge_cmd, dest_dir=qfilt_parent)

        return output_list
    else:
        return input_files

def quality_filter(input_files, platform, file_map):

    settings = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings')
    run_name = settings['run_details']['run_name']
    qfilt_dict = settings['quality_filtering'][platform]
    min_len = qfilt_dict['min_len']
    max_len = qfilt_dict['max_len']
    premerged = settings['quality_filtering']['illumina']['merge_reads']

    if platform == 'illumina':
        max_error = qfilt_dict['max_error']
    else:
        max_error = ''

    qfilt_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output']['quality-filtered'])

    qfilt_path = mkdir_exist_ok(new_dir=f'./{QUALFILT_PREFIX}_{run_name}', parent_dir=qfilt_parent)
    nofilt_path = mkdir_exist_ok(new_dir=f'./{flip_prefix(QUALFILT_PREFIX)}_{run_name}', parent_dir=qfilt_parent)

    if premerged:
        for file in input_files:
            qfilt_out = add_prefix(file_path=file, prefix=QUALFILT_PREFIX,
                                           dest_dir=qfilt_path, action=None)
            nofilt_out = add_prefix(file_path=file, prefix=flip_prefix(QUALFILT_PREFIX),
                                        dest_dir=nofilt_path, action=None)

            vsearch_filt_cmd = ['vsearch', '--fastq_filter', file,
                                '--fastqout', qfilt_out,
                                '--fastqout_discarded', nofilt_out,
                                '-fastq_maxee', str(max_error), '--fastq_maxlen', str(max_len),
                                '--fastq_minlen', str(min_len), '--sizeout']

            run_subprocess(vsearch_filt_cmd, dest_dir=qfilt_parent)
    else:
        paired_dict = pair_reads(input_files)
        for file in paired_dict.keys():

            qfilt_fwd_out = add_prefix(file_path=file, prefix=QUALFILT_PREFIX,
                                       dest_dir=qfilt_path, action=None)
            nofilt_fwd_out = add_prefix(file_path=file, prefix=flip_prefix(QUALFILT_PREFIX),
                                        dest_dir=nofilt_path, action=None)

            qfilt_rev_out = add_prefix(file_path=paired_dict[file], prefix=QUALFILT_PREFIX,
                                       dest_dir=qfilt_path, action=None)
            nofilt_rev_out = add_prefix(file_path=paired_dict[file], prefix=flip_prefix(QUALFILT_PREFIX),
                                        dest_dir=nofilt_path, action=None)

            vsearch_filt_cmd = ['vsearch', '--fastq_filter', file, '--reverse', paired_dict[file],
                                '--fastqout', qfilt_fwd_out, '--fastqout_rev', qfilt_rev_out,
                                '--fastqout_discarded', nofilt_fwd_out, '--fastqout_discarded_rev', nofilt_rev_out,
                                '-fastq_maxee', str(max_error), '--fastq_maxlen', str(max_len),
                                '--fastq_minlen', str(min_len), '--sizeout']

            run_subprocess(vsearch_filt_cmd, dest_dir=qfilt_parent)

    return None

def dereplicate(input_files, derep_step, platform, file_map):

    derep_prefix = DEREP_PREFIX + '0' + str(derep_step)

    if derep_step == 1:
        out_tag = 'derep-full-length'
    else:
        out_tag = 'derep-subregions'

    settings = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings')
    run_name = settings['run_details']['run_name']
    min_unique = settings['quality_filtering'][platform]['min_count'][f'derep0{derep_step}']
    premerged = settings['quality_filtering']['illumina']['merge_reads']

    derep_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output'][out_tag])

    derep_path = mkdir_exist_ok(new_dir=f'./{derep_prefix}_{run_name}', parent_dir=derep_parent)

    # derep_summary = mkdir_exist_ok(new_dir=f'./{derep_prefix}_{run_name}_summaries', parent_dir=derep_parent)
    if premerged:
        pass
    else:
        input_files = [f for f in input_files if re.search('R1', f, re.I)]  # only take fwd if not merged
        # rename files to be just sample ID or keep R1 to be clear they're forward reads only?

    for file in input_files:

        # derep always outputs a fasta format, but need to provide file name as fasta format or will appear as fastq
        # even though it isn't really fastq when you open it up
        derep_output = add_prefix(file_path=file, prefix=derep_prefix, dest_dir=derep_path, action=None).with_suffix('.fasta')

        vsearch_derep_cmd = ['vsearch', '--derep_fulllength', file,
                             '--output', derep_output,'--minuniquesize', str(min_unique),
                             '--sizeout']

        run_subprocess(vsearch_derep_cmd, dest_dir=derep_parent)

    return None

def separate_subregions(input_files, file_map):

    settings = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings')
    run_name = settings['run_details']['run_name']

    itsx_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output']['separate-subregions'])

    itsx_path = mkdir_exist_ok(new_dir=f'./{ITSX_PREFIX}_{run_name}',  parent_dir=itsx_parent)

    for file in input_files:

        itsx_output_basename = add_prefix(file_path=file, prefix=ITSX_PREFIX, dest_dir=itsx_path, action=None).with_suffix('')

        itsx_command = ['ITSx', '-i', file, '-o', itsx_output_basename, '-t', 'F', '--multi_thread', 'T',
                        '--save_regions', 'all']

        run_subprocess(itsx_command, dest_dir=itsx_parent)

    return None

def concat_regions(dir_path, regions_to_concat=['ITS1', '5_8S', 'ITS2'], header_delim=';', **kwargs):
    '''
    Concatenates provided regions for each read.

    Looks in the provided directory for files containing the substring of the
    regions provided in the regions_to_concat list. For each region in the list,
    it will store read ID, region, read length, and full-length read count in a
    dictionary, where each read ID is the primary key with details stored for each
    region (read length, full-length read count). Then will combine this information
    from the regions for each read ID. The sequences from each region are concatenated
    in the order they appear in the provided regions_to_concat list. The read length of
    the final concatenated reads is the sum of the read lengths of the regions. It will
    check that the full-length read count is the same for each of the regions, which
    further confirms that the regions are from the same full-length read. Then will write
    these concatenated reads and updated headers to a new output fasta file.
    :param dir_path: directory containing the fasta files to concatenate.
    :param regions_to_concat: an ORDERED list of regions to concatenated, with the 5' region
    first and the 3' region last; defaults to ITS regions
    :param header_delim: delimiter to use for separating the components of the header in
    the concatenated fasta output; defaults to semicolon
    :return: None, but will write out a fasta file
    '''

    # REMOVE AFTER TESTING
    # dir_path = Path('./04_separated-subregions/itsx_pacbio_sporocarp-f_2023-12_Q40')
    # regions_to_concat = ['ITS1', '5_8S', 'ITS2']
    # header_delim = ';'
    ######################
    for kw in kwargs.keys():
        if re.search('suffix', kw, re.I):
            file_suffix = f'_{kwargs[kw]}'
        else:
            file_suffix = ''

    # create dictionary for acceptable input
    accepted_regions = {'(^S.+?(?=\bs)|SSU)': 'SSU',  # could be better but not going to obsess
                        'ITS.?1': 'ITS1',  # checks for possible punctuation within
                        '5.?8S': '5.8S',  # checks for possible punctuation within
                        'ITS.?2': 'ITS2',  # checks for possible punctuation within
                        '(^L.+?(?=\bs)|LSU)': 'LSU',  # similar to SSU
                        '(full)?.?ITS(?!.?\d)': 'ITS'}  # may include full, but can't have number after (w or w/o punct)

    # check input list for valid entries, and then format string to consistent format
    # returns dict now with formatted string as key, regex as value
    def create_region_dict(region_list):
        output_dict = {}
        for region_str in region_list:
            for r in accepted_regions.keys():  # compare input region to each regex in dict
                if re.search(r, region_str, re.I):  # if it does match, stop loop and return the formatted version
                    output_dict.update({accepted_regions[r]:r})
                    break
                else:
                    continue
            else:
                return print(f'FAILURE. Did not recognize the input region {region_str}. Please choose from the '
                             f'recognized regions: {list(accepted_regions.values())}\n')
        return output_dict

    region_search_dict = create_region_dict(regions_to_concat)

    # make nested dict: read ID > region > [sequence, read len, derep01_reads(?)]
    # ADD CHECK FOR MULTIPLE FILES MATCHING THE REGION IN THIS DIRECTORY
    def search_path_with_regex(dir_path, regex, return_all=False):
        found_file_paths = []
        for file in dir_path.glob('*'):
            if re.search(regex, file.name, re.I):
                found_file_paths.append(file)
            else:
                continue

        if len(found_file_paths) > 1:
            if return_all:
                return found_file_paths
            else:
                print(f'WARNING. Multiple files matching the provided regex {regex} were located in '
                      f'{dir_path.name}. Returning only the first match.\n')  # REPLACE WITH OPTIONS?
                return found_file_paths[0]
        elif len(found_file_paths) == 0:
            return print(f'FAILURE. No files matching the provided regex {regex} were located in '
                         f'{dir_path.name}. Exiting...\n')
        else:
            return found_file_paths[0]

    fastas_to_concat = [search_path_with_regex(dir_path, regex=r) for r in region_search_dict.values()]

    SAMPLE_ID_RE = '(?<=\w_)(pacbio|sanger|illumina).+?(?=\.)'
    PREFIX_RE = '^\w.+?(?=_pacbio|_sanger|_illumina)'

    sample_name = re.search(SAMPLE_ID_RE, fastas_to_concat[0].name, re.I).group(0)
    file_prefix = re.search(PREFIX_RE, fastas_to_concat[0].name, re.I).group(0)

    if list(region_search_dict.keys()) == ['ITS1', '5.8S', 'ITS2']:
        concat_name = 'full_ITS'
    else:
        concat_name = input(f'Please provide a name to assign to this combination of subregions (use '
                            f'underscores instead of spaces):\n')

    # renamed header search regex
    # WAY TO UPDATE REGEX GLOBALLY WHEN A RENAMING OCCURS?
    READ_ID_RE = '^.+?(?=;)'
    READ_REGION_RE = '(?<=\d;)\w.+?(?=;)'
    READ_LEN_RE = '(?<=len=)[0-9]{1,}(?=bp)'
    READ_COUNT_RE = '(?<=full\-len_copies=)[0-9]{1,}'

    # ADD THESE OPTIONS IN THE CONFIG FILE?
    details_to_include = ['sequence', 'region_len', 'full-len_copies']

    concat_dict = {}
    for fasta in fastas_to_concat:
        for record in SeqIO.parse(fasta, 'fasta'):
            try:
                read_id = re.search(READ_ID_RE, record.id).group(0)
                region = re.search(READ_REGION_RE, record.id, re.I).group(0)

                length = re.search(READ_LEN_RE, record.id, re.I).group(0)
                derep01_reads = re.search(READ_COUNT_RE, record.id, re.I).group(0)
                details = {k:v for k,v in zip(details_to_include, [record.seq, length, derep01_reads])}

                if read_id in concat_dict:  # if this sample is already in the dictionary...
                    concat_dict[read_id].update({region: details})  # append the new region info to it
                else:  # if it isn't yet in the dictionary...
                    concat_dict[read_id] = {region: details}  # add a key value pair for this sample
            except:
                print(f'There was an issue extracting information from the fasta header of the input file '
                      f'{fasta}. This will occur if the sequence files produced by ITSx have not yet '
                      f'been reformatted. Please run the rename_read_header function and retry.\n')
                return sys.exit()

    # concat sequences by read
    concatenated_records = []
    for read in concat_dict.keys():
        # create lists to append to each region's records to
        record_details = {k:[] for k in details_to_include}

        # pull info from each region for this read
        for region in region_search_dict.keys():
            region_dict = concat_dict[read][region]
            for det in details_to_include:
                if det in record_details:
                    record_details[det].append(region_dict[det])  # append before sum to ensure all regions there
                else:
                    record_details[det] = region_dict[det]  # append before sum to ensure all regions there

        # confirm all regions data was acquired, then combine
        def join_if_all_present(pulled_data_dict):

            joined_data = {'seq':'',
                           'id':[read, concat_name]}

            for k in pulled_data_dict.keys():
                if len(pulled_data_dict[k]) == len(regions_to_concat):
                    if any(isinstance(i, Seq) for i in pulled_data_dict[k]):
                        joined_data['seq'] = Seq('').join(pulled_data_dict[k])  # must join with empty seq
                    else:
                        if re.search('len$', k, re.I):
                            val = sum([int(x) for x in pulled_data_dict[k]])
                            joined_data['id'].append(f'{k}={val}bp')
                        elif re.search('^full\-len|count', k, re.I):
                            if len(set(pulled_data_dict[k])) == 1:
                                derep = int(pulled_data_dict[k][0])  # if all same, doesn't matter which you chose
                                joined_data['id'].append(f'{k}={derep}')  # make int so no '' surrounding num in header
                        else:
                            print(f'Could not recognized the input data type, {k}:\n'
                                  f'\t{pulled_data_dict[k]}\n')
                            return None
                else:
                    msg = f'Did not successfully pull all {k}s from the regions to concatenate ' \
                          f'({len(pulled_data_dict[k])} out of {len(regions_to_concat)} were pulled). ' \
                          f'Debug before continuing...\n'
                    print(msg)
                    return None  # REPLACE W/ sys.exit AFTER DEBUG

            return joined_data

        combined_data = join_if_all_present(record_details)
        concat_record = SeqRecord(combined_data['seq'],
                                  name = read,
                                  id = header_delim.join(combined_data['id']),
                                  description = '')
        concatenated_records.append(concat_record)

    fasta_out = dir_path / f'{file_prefix}_{sample_name}.{concat_name}{file_suffix}.fasta'
    SeqIO.write(concatenated_records, fasta_out, 'fasta')

    return None

def check_itsx_output(itsx_dir, full_len_dir, num_bp_compare, write_to_log=True, same_threshold=99):
    SAMPLE_ID_RE = '(?<=\w_)(pacbio|sanger|illumina).+?(?=\.)'
    READ_ID_RE = '[0-9]{1,}(?=\/ccs)'

    bp = int(num_bp_compare)

    # import records as dictionary, where the read ID is the key
    def name_as_key(record):
        header = record.id
        sample_read_id = header.split(';')[0]
        return sample_read_id

    # import full ITS records and LSU (for revcomp comparison)
    full_its_records = SeqIO.to_dict(SeqIO.parse(next(itsx_dir.glob('*full_ITS*')), 'fasta'), key_function=name_as_key)
    lsu_records = SeqIO.to_dict(SeqIO.parse(next(itsx_dir.glob('*LSU*')), 'fasta'), key_function=name_as_key)

    # collect the full length _ bp sequences at start, full ITS _ bp sequences at start, and LSU _ bp seqs at end
    post_itsx_diffs = {'read_id':[],
                       'full_start':[],
                       'its_start':[],
                       'lsu_end':[]}
    for record in SeqIO.parse(full_len_dir, 'fasta'):

        sample_id = re.search(SAMPLE_ID_RE, full_len_dir.name).group(0)
        read_num = re.search(READ_ID_RE, record.id).group(0)
        read_id = '_'.join([sample_id, read_num])

        renamed_record = SeqRecord(record.seq, id=read_id, description='')

        full_len_seq = renamed_record.seq
        full_its_seq = full_its_records[renamed_record.id].seq
        try:
            lsu_seq = lsu_records[renamed_record.id].seq
        except:
            lsu_seq = 'N'*bp


        post_itsx_diffs['read_id'].append(renamed_record.id)
        post_itsx_diffs['full_start'].append(str(full_len_seq[:bp]))
        post_itsx_diffs['its_start'].append(str(full_its_seq[:bp]))
        post_itsx_diffs['lsu_end'].append(str(lsu_seq[-bp:]))

    # create df from this dictionary
    diff_df = pd.DataFrame.from_dict(post_itsx_diffs)

    # compare the read slices first as full-ITS fwd, then full-LSU revcomp
    comparison_result = []  # True if same
    equal_as = []  # str describing nature of match: fwd, revcomp, no LSU, unresolved (error)
    for i in range(diff_df.shape[0]):
        full_fwd = Seq(diff_df['full_start'][i])
        its_fwd = Seq(diff_df['its_start'][i])
        lsu_end = Seq(diff_df['lsu_end'][i])
        if full_fwd == its_fwd:
            comparison_result.append(True)
            equal_as.append('forward_strand')
            continue
        else:
            if 'N' in lsu_end:
                comparison_result.append(False)
                equal_as.append('no_lsu')
            else:
                its_revcomp = lsu_end.reverse_complement()
                if full_fwd == its_revcomp:
                    comparison_result.append(True)
                    equal_as.append('reverse_comp')
                    continue
                else:
                    comparison_result.append(False)
                    equal_as.append('unresolved')

    diff_df['comparison'] = comparison_result
    diff_df['equal_as'] = equal_as

    percent_same = (sum(diff_df['comparison']) / diff_df.shape[0])*100
    revcomp_detect = ((diff_df[diff_df['equal_as'] == 'reverse_comp'].shape[0]) /sum(diff_df['comparison']))*100
    forward_detect = ((diff_df[diff_df['equal_as'] == 'forward_strand'].shape[0]) /sum(diff_df['comparison']))*100
    missing_lsu = ((diff_df[diff_df['equal_as'] == 'no_lsu'].shape[0]) /sum(diff_df['comparison']))*100

    if write_to_log:
        log_dir = itsx_dir / 'itsx_inspection_logs'
        log_dir.mkdir(exist_ok=True)
        diff_df.to_csv(log_dir / f'itsx_inspection_comparison.csv', index=False)
        with open((log_dir / 'itsx_inspection_summary.log'), 'wt') as fout:
            fout.write(f'sample: {re.search(SAMPLE_ID_RE, full_len_dir.name).group(0)}\n')
            fout.write(f'date: {datetime.today().strftime("%Y-%M-%d")}\n\n')
            fout.write(f'All full-length reads were compared to their ITS region after running ITSx. The first'
                       f' {bp} bp of the full-length reads were aligned to the first {bp} bp of their corresponding '
                       f'full ITS read that was concatenated after running ITSx. If the first {bp} bp did not '
                       f'match, then the reverse complement of the last {bp} bp of their LSU region was compared '
                       f'to the first {bp} bp of the full-length read. This is because reads are not reoriented '
                       f'until done so by ITSx. \n\n')
            fout.write(f'For {percent_same:.2f}% of these comparisons (out of {diff_df.shape[0]} total comparisons), '
                       f'the first {bp} bp of the full length read matched the first {bp} bp of either the forward '
                       f'ITS or reverse complement of the end of the LSU for that read.\n\n')
            fout.write(f'Out of the {percent_same:.2f}% that matched, {forward_detect:.2f}% were matched by comparing '
                       f'the full-length read to the start of the ITS, while {revcomp_detect:.2f}% were matched by '
                       f'comparing the full-length read to the reverse complement of the end of the LSU region.\n\n')
            if missing_lsu > 0:
                fout.write(f'{int(missing_lsu*diff_df.shape[0])} reads ({missing_lsu:.2f}% of total) did not have '
                           f'an LSU region to compare to, since ITSx could not locate the LSU for the read(s). This '
                           f'accounts for the proportion of reads that could not be matched.\n')
    else:
        if percent_same < same_threshold:
            print(f'WARNING. The percent of reads that matched to their pre-ITSx full-length read was above '
                  f'{same_threshold}%, at {percent_same}.\n')

    return None

def check_chimeras(input_files, file_map, ref=None):
    chim_ref_files = file_map['config']['reference-db']['chimera'].glob(SEQ_FILE_GLOB)

    settings = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings')
    run_name = settings['run_details']['run_name']

    chim_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output']['chimera-checked'])

    nochim_path = mkdir_exist_ok(new_dir=f'./{NOCHIM_PREFIX}_{run_name}', parent_dir=chim_parent)
    chim_path = mkdir_exist_ok(new_dir=f'./{flip_prefix(NOCHIM_PREFIX)}_{run_name}', parent_dir=chim_parent)

    uchime_log = chim_parent / 'uchime.log'

    if ref is None:
        for file in input_files:
            nochim_out = add_prefix(file_path=file, prefix=NOCHIM_PREFIX,
                                    dest_dir=nochim_path, action=None)
            chim_out = add_prefix(file_path=file, prefix=flip_prefix(NOCHIM_PREFIX),
                                  dest_dir=chim_path, action=None)

            vsearch_denovo_cmd = ['vsearch', '--uchime3_denovo', file,
                                  '--chimeras', chim_out,
                                  '--nonchimeras', nochim_out,
                                  '--uchimeout', uchime_log]

            run_subprocess(vsearch_denovo_cmd, dest_dir=chim_parent)
    else:
        refs_to_use = []
        for r in chim_ref_files:
            result = re.search(ref, r, re.I)
            if result:
                refs_to_use.append(result.group(0))
            else:
                continue

        for chim_ref in refs_to_use:
            for file in input_files:
                nochim_out = add_prefix(file_path=file, prefix=NOCHIM_PREFIX,
                                        dest_dir=chim_path, action=None)
                chim_out = add_prefix(file_path=file, prefix=flip_prefix(NOCHIM_PREFIX),
                                      dest_dir=chim_path, action=None)

                vsearch_ref_cmd = ['vsearch', '--uchime_ref', file,
                                   '--chimeras', chim_out,
                                   '--nonchimeras', nochim_out,
                                   '--uchimeout', uchime_log,
                                   '--db', chim_ref]

                run_subprocess(vsearch_ref_cmd, dest_dir=chim_parent)

    return None

def combine_all_reads(input_files, file_map):

    settings = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings')
    run_name = settings['run_details']['run_name']

    total_record_list = []
    for file in input_files:
        rename_read_header(file, header_delim=';')
        for record in SeqIO.parse(file, 'fasta'):
            total_record_list.append(record)

    combined_out = file_map['pipeline-output']['otus-clustered'] / f'combined_otus_{run_name}.fasta'
    SeqIO.write(total_record_list, combined_out, 'fasta')

    return combined_out

def cluster_reads(input_files, file_map):
    settings = import_config_as_dict(file_path=file_map['config']['main'], file_handle='pipeline-settings')
    run_name = settings['run_details']['run_name']

    min_threshold = settings['clustering']['min_threshold']

    clust_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output']['otus-clustered'])

    # rename headers for each file to retain source of sample
    # rename_read_headers()

    # combine all reads into one fasta file
    combined_reads = combine_all_reads(input_files, file_map)

    clust_out = add_prefix(file_path=combined_reads, prefix=CLUSTER_PREFIX,
                           dest_dir=clust_parent, action=None)
    # cluster
    vsearch_clust_cmd = ['vsearch', '--cluster_fast', combined_reads, '--clusters',
                         clust_out, '--id', str(min_threshold)]

    run_subprocess(vsearch_clust_cmd, dest_dir=clust_parent)

    return None


def create_blast_db(config_dict, file_map, taxa_list=None):
    '''
    Create a reference dataset for BLAST+ search.

    :param config_dict: configuration file dictionary; required to get the list
    of reference databases to use, as defined by user in the configuration file.
    :param file_map: file mapping, from mapping.py
    :param taxa_list: option; list of taxa to include if wanting to limit search
    to a specific taxonomic group
    :return: returns path to the newly created db for blastn search; creates a custom
    database using BLAST+ makeblastdb from the command line
    '''
    tax_settings = config_dict['taxonomy']
    db_dir = file_map['config']['reference-db']

    include_genbank = tax_settings['refdb']['genbank']['include']
    include_unite = tax_settings['refdb']['unite']['include']
    include_custom = tax_settings['refdb']['custom']['include']
    include_maarjam = tax_settings['refdb']['maarjam']['include']

    included_tag = ''
    db_include_list = []
    if include_genbank:
        genbank_fasta = db_dir.glob('*genbank*')
        db_include_list.append(genbank_fasta)
        included_tag += 'G'

    if include_unite:
        unite_fasta = db_dir.glob('*unite*')
        db_include_list.append(unite_fasta)
        included_tag += 'U'

    if include_custom:
        custom_fasta = db_dir.glob('*custom*')
        db_include_list.append(custom_fasta)
        included_tag += 'C'

    if include_maarjam:
        maarjam_fasta = db_dir.glob('*maarjam*')
        db_include_list.append(maarjam_fasta)
        included_tag += 'M'

    custom_search_records = []
    for db in db_include_list:
        for record in SeqIO.parse(db, 'fasta'):
            if taxa_list is None:
                custom_search_records.append(record)
            else:
                pass  # NEED TO SEE HOW TO GET ONLY SPECIFIC TAXONOMY, NOT SURE HOW FORMATTED

    tax_out_dir = file_map['pipeline-output']['taxonomy']

    search_date = datetime.today.strftime('%Y-%M-%d')
    output_basename = f'custom-ref-{included_tag}_{search_date}'
    output_path = (tax_out_dir / output_basename).with_suffix('.fasta')

    SeqIO.write(custom_search_records, output_path, 'fasta')

    output_db = tax_out_dir / f'{output_basename}_blast'
    blast_cmd = ['makeblastdb', '-in', output_path, '-title', output_db, '-dbtype', 'nucl',
                 output_db]

    run_subprocess(blast_cmd, dest_dir = output_db)

    return output_db

def assign_taxonomy(config_dict, file_map, taxa_list=None):
    '''
    Assign taxonomy to sequences.

    :param config_dict:
    :param file_map:
    :return:
    '''
    tax_output = file_map['pipeline-output']['taxonomy']

    run_name = config_dict['run_details']['run_name']

    if config_dict['taxonomy']['algorithm'] == 'blastn':
        ref_db = create_blast_db(config_dict, file_map, taxa_list=None)

        blast_out = (tax_output / f'{run_name}').with_suffix('.txt')
        blast_cmd = ['blastn', '-query', ref_db, '-out', blast_out]

        run_subprocess(blast_cmd, dest_dir = tax_output)

    elif config_dict['taxonomy']['algorithm'] == 'rdp':
        print(f'Not sure how to set up RDP yet...')
        rdp_cmd = []

        run_subprocess(rbp_cmd, dest_dir = tax_output)

    else:
        print(f'I don\'t yet have an option for that algorithm.')

    return None