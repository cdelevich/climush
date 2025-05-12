import json, re, os, warnings
from multiprocessing import Pool, cpu_count
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
from datetime import datetime
from climush.constants import *
from climush.utilities import *

def demultiplex(output_dir, reference_dir, multiplexed_files, verbose):

    # DEFINE RELEVANT FILE PATHS ####################################################################
    # get the paths needed for demux
    mapping_files = file_finder(
        reference_dir=reference_dir,
        search_glob='pacbio*mapping*',
        multiple_matches=True,
    )

    # ACCESS USER SETTINGS FROM CONFIG ##############################################################
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']  # name to use for this bioinformatics run

    # LOCATE MAPPING FILE FOR EACH MULTIPLEXED SEQUENCE FILE #########################################
    # get all unique run queue IDs from files needing demux; will use later on
    queue_ids = set()

    # link multiplexed sequence files to their corresponding mapping file in the config directory
    demux_mapping = {}  # key is mapping file path, value(s) is sequence file path in needs_demux dir

    for mp_file in multiplexed_files:

        # get the sequencing core's queue ID from multiplexed seq file and corresponding climush sample ID for this queue ID
        qid = re.search(r'^(\d{4})', mp_file.name).group(0)  # get queue ID from file name
        queue_ids.add(qid)  # add the qid to the set of all queue ids requiring demux
        cid = settings['demultiplex']['multiplex'][qid]  # get climush - queue ID pairing from config
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
            right_map = prompt_print_options(qid_files,
                                             auto_respond=settings['automate']['auto_respond'])  # get right mapping file from user input
            if right_map in demux_mapping.keys():  # add to dict
                demux_mapping[right_map].append(mp_file)  # update list of corresponding multiplexed files w/ str
            else:
                demux_mapping.update({right_map: [mp_file]})  # add queue ID and first multiplexed file as list

    # CREATE A FASTA FILE OF ALL UNIQUE BARCODES USED TO MULTIPLEX ###################################
    # define subfunctions to help locate and prepare barcodes for demultiplexing
    def get_col_name(pattern, df, auto_respond=False):

        # search for columns that match the provided pattern
        matches = [c for c in df.columns if re.search(pattern, c, re.I)]

        # output depends on the number of matches returned (ideally 1)
        if len(matches) == 1:  # if one match, return that column name
            return matches[0]

        else:  # if multiple or no matches, rerun with new user-input regex

            # check number of matches, will print out different prompt, but otherwise executes same thing
            if len(matches) == 0:  # if no matches
                error_msg = f'No column matching the pattern \'{pattern}\' was located in the input dataframe. Please try '\
                            f'another regular expression that better matches a single column out of these columns in the '\
                            f'dataframe: '
            else:
                error_msg = f'{len(matches)} column names match the pattern \'{pattern}\' in the input dataframe. Please try '\
                            f'another regular expression that better matches a single column out of these columns in the '\
                            f'dataframe: '

            print(error_msg)

            print_indented_list(matches)  # print out the columns of the dataframe (formatted to be indented)

            if auto_respond:
                pass  # INCOMPLETE
            else:
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
            if re.search(r'pool', tab, re.I):  # ...only if relevant to mapping/demux

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
        fwd_bc_noprimer = [re.sub(f'{primers["fwd_primer"]}$', r'', fwd) for fwd in fwd_barcodes]
        rev_bc_noprimer = [re.sub(f'{primers["rev_primer"]}$', r'', rev) for rev in rev_barcodes]

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
                if verbose:
                    print(f'ERROR. Barcodes in the mapping file should be the same length, but are not.')
                    print(f'number of unique lengths of barcodes w/ primers: {len(len_bc_primer)}')
                    print(f'unique lengths of barcodes w/ primers:           {len_bc_primer}\n')
                return False

            # COMPARISON
            # if all barcodes are the same length after removing primer, and they match the expected length, return True
            if (len(len_bc_noprimer) == 1) and (int(len_bc_noprimer[0]) == len_bc_expected):
                return True
            else:
                if verbose:
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
    # primer_dict['fwd_primer'] = settings['primers']['fwd']['sequence']['pacbio']
    primer_dict['fwd_primer'] = list(settings['primers']['fwd']['pacbio'].values())[0]
    primer_dict['rev_primer'] = list(settings['primers']['rev']['pacbio'].values())[0]

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
                    err_msg = f'WARNING. {len(extra_bc)} new barcode(s) added to the overall list of unique barcodes:\n'\
                              f'\textra {b.split("_")[0]} barcode(s): {extra_bc}\n'\
                              f'All mapping files should contain the same set of unique forward and reverse barcodes. '\
                              f'Continuing with demultiplexing, but if there\'s a future error, see the mapping file '\
                              f'that triggered this warning: {map.name}\n'
                elif num_bc_overall > num_bc_mapping:  # if not all possible barcodes are in this mapping file...
                    missing_bc = unique_barcodes[b].difference(set(seq_run_bc[b]))  # find missing bc
                    err_msg = f'WARNING. {len(missing_bc)} barcode(s) were missing from this mapping file\'s barcodes:\n'\
                              f'\tmissing {b.split("_")[0]} barcode(s): {missing_bc}\n'\
                              f'All mapping files should contain the same set of unique forward and reverse barcodes. '\
                              f'Continuing with demultiplexing, but if there\'s a future error, see the mapping file '\
                              f'that triggered this warning: {map.name}\n'
                else:
                    pass  # they are equal, as expected

                if verbose:
                    print(err_msg)

            # regardless of above outcome, add barcodes from this mapping file to the overall dict of unique barcodes
            unique_barcodes[b].update(seq_run_bc[b])  # use update to add list to set

    # write out fasta file of unique barcodes; must be ordered so same numbering is used every time this is run
    # create new path in pipeline output, and use this path to create fasta file name
    fasta_out_dir = mkdir_exist_ok(new_dir = 'barcodes', parent_dir = output_dir)
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
    lima_output = mkdir_exist_ok(new_dir='lima_output', parent_dir=output_dir)
    lima_subdirs = []  # create list of output directories to read lima output from later on

    for qid in queue_ids:

        # get the climush sequencing run ID that corresponds to this queue ID
        seq_run = settings['demultiplex']['multiplex'][qid]

        # create an output directory for this sequencing run using the climush sequencing run ID
        seq_run_output = mkdir_exist_ok(new_dir=seq_run, parent_dir=lima_output) # create an output dir for each seq run
        lima_subdirs.append(seq_run_output)

        # locate all multiplexed sequencing files that match this queue ID
        pools_to_demux = [f for f in multiplexed_files if re.search(f'^{qid}', f.name)]

        # alert if only one pool is detected
        if len(pools_to_demux) == 1:
            msg = f'WARNING. There is only one multiplexed sequencing file ({pools_to_demux[0].name}) for sequencing ' \
                  f'run {seq_run} when two are typically expected: one for each pool of sequences that were ' \
                  f'sequenced simultaneously on separate SMRT cells. Do you wish to continue with only half of the ' \
                  f'sequences from this sequencing run?'  # yes/no/quit automatically added by function
            prompt_yes_no_quit(msg, auto_respond=settings['automate']['auto_respond'])  # if yes, auto-continues; if no or quit, exits

        for pool in pools_to_demux:

            # detect pool number from the file name; remove queue ID first to avoid any confusion
            try:  # should be able to do automatically
                pool_num = re.search(r'\d', re.sub(f'^{qid}' + r'\.', r'', pool.name)).group(0)
            except AttributeError:  # but if re.search does not return a match
                print(f'The sequencing pool number for the multiplexed file {pool.name} could not be detected in ' \
                      f'the file name. Please enter the pool number corresponding to this file as a single digit'
                      f'(i.e., no leading zeros or decimals): ')
                pool_num = input()

            # define prefix to use for lima output files; include path to the output directory for this seq run
            out_prefix = (seq_run_output / f'lima-demux_{seq_run}_pool{pool_num}').with_suffix('.fasta')

            # assemble list of commands required to run lima demultiplexing
            lima_cmd = ['lima', pool, barcode_fasta, out_prefix, '--min-score', '93', '--hifi-preset', 'ASYMMETRIC']

            run_subprocess(cli_command_list=lima_cmd, dest_dir=lima_output, run_name=run_name,
                           auto_respond=settings['automate']['auto_respond'])

    # CREATE DICTIONARY OF BARCODE PAIRS FOR ALL SAMPLES
    sample_barcodes = {q:{'pool1': {},
                          'pool2': {}} for q in queue_ids}

    for dxmap in demux_mapping:

        # get the queue ID from the name of the multiplexed sequencing file; need it to find key in sample_barcodes
        qid = re.search(r'^(\d{4})', list(demux_mapping[dxmap])[0].name).group(0)

        # import the mapping file dataframe
        mapping_df = import_mapping_df(df_path=dxmap, auto_respond=settings['automate']['auto_respond'])

        # loop through each tab (= pool1/pool2) in the df to get the sample ID and its barcode combo
        for tab in mapping_df:

            # get the name of the column used for the sample ID, fwd barcode, rev barcode
            smp_col = get_col_name(pattern=SAMPLE_COL_RE, df=mapping_df[tab])
            fwd_col = get_col_name(pattern=FWD_COL_RE, df=mapping_df[tab])
            rev_col = get_col_name(pattern=REV_COL_RE, df=mapping_df[tab])

            for i in range(mapping_df[tab].shape[0]):

                # from the dataframe, get the sample ID, fwd barcode sequence, and rev barcode sequence
                smp_id = mapping_df[tab][smp_col][i].strip()  # make sure whitespace is stripped from sample ID
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
        cid = re.search(r'^pacbio.+', subdir.name).group(0)  # get climush ID, which will be in name of file
        mp_dict = settings['demultiplex']['multiplex']  # make shorter for list comp below
        qid = [k for k in mp_dict if cid in mp_dict[k]][0]  # convert to qid to match dict

        for p in ['pool1', 'pool2']:

            # find the fasta file for this pool, which has the read ID and the barcode combination
            lima_fasta = [f for f in subdir.glob('*.*') if re.search(f'{subdir.name}_{p}' + r'\.fasta', f.name, re.I)][0]

            # open the fasta file and pull the indices of the forward and reverse barcodes
            with open(lima_fasta, 'r') as fin:
                fasta_txt = fin.readlines()
                for l, line in enumerate(fasta_txt):
                    if line.startswith('>'):
                        read_id = re.findall(r'(?<=>).+?(?=\sbc)', line, re.I)[0]
                        read_barcode_i = list(map(int, re.search(r'(?<=bc=)\d{1,2},\d{1,2}', line, re.I).group(0).split(',')))
                        read_bc_names = [bc_name_index[i] for i in read_barcode_i]  # convert index to name
                        read_barcodes[qid][p].update({read_id: read_bc_names})  # add read's bc names to dict
                        read_seqs[qid][p].update({read_id: fasta_txt[l+1]})  # add read sequence to dict
                    else:
                        continue

    # COMPARE LIMA BARCODE COMBINATIONS TO MAPPING FILE BARCODE COMBINATIONS

    # create an output subdirectory in demultiplexed dir for all demultiplexed reads in this bioinformatics run
    demux_dir = mkdir_exist_ok(new_dir = f'{DEMUX_PREFIX}_{run_name}', parent_dir=output_dir)

    # for each sequencing run...
    for qid in read_barcodes:

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

                # get the number of samples that match to this combination of barcodes
                num_matches = len(sample_ids)

                # confirm that the read's barcode combination matches only one sample ID, as it should
                if num_matches == 1:

                    # get sample ID as string, not list (now that it is confirmed there's only one)
                    sample_id = sample_ids[0]

                    # create read ID header from sequencing run, sample ID, and last part of original read header
                    full_sample_id = '_'.join([seq_run, str(sample_id)])  # sample ID is int (?) debugging...
                    read_id = '_'.join([full_sample_id, r.split('/')[-2]])

                    # get the barcode-free sequence for this read
                    read_seq = seq_subdict[r]

                    # write the read ID and read sequence to the sample ID demultiplexed fasta file
                    sample_fasta = (demux_dir / f'{DEMUX_PREFIX}_{full_sample_id}').with_suffix('.fasta')
                    with open(sample_fasta, 'a') as fout:
                        fout.write(f'>{read_id}\n')
                        fout.write(f'{read_seq}\n')

                else:

                    # change the error message based on number of matches detected
                    if num_matches > 1:
                        err_msg = f'WARNING. The barcode combination detected in read {r} from {p} of {seq_run} matched to '\
                                  f'{len(sample_ids)} sample IDs, when it should only match to one. This read was not '\
                                  f'sorted, and this error was recorded to the error table.\n'
                        sample_ids_err = ", ".join(sample_ids)
                    else:
                        err_msg = f'WARNING. The barcode combination detected in read {r} from {p} of {seq_run} matched to '\
                                  f'{len(sample_ids)} sample IDs, when it should only match to one. This read was not '\
                                  f'sorted, and this error was recorded to the error table.\n'
                        sample_ids_err = "no matches"

                    if verbose:
                        print(err_msg)

                    barcode_error = (output_dir / 'barcode_errors').with_suffix('.tsv')

                    # if the error file doesn't yet exist, start with writing the header
                    if not barcode_error.is_file():
                        with open(barcode_error, 'a') as fout:
                            fout.write(f'sequencing_run\t\tsequencing_pool\tsample_ids\tread_id\t\t\t\tbarcode_combination\n')

                    # write out details of the read that triggered this error
                    with open(barcode_error, 'a') as fout:
                        fout.write(f'{seq_run}\t{p}\t{sample_ids_err}\t{r}\t{", ".join(bc_list)}\n')

def pair_reads(input_files):

    # accept either list of files or directory as input
    input_files = create_file_list(input_files)

    pairs_dict = {}
    rev_reads = []

    for file in input_files:
        # if the file is an R1 read file...
        if re.search(r'R1', file.stem, re.I):
            pairs_dict.update({file:''})  # add to dict, with empty str as value
        # if the file is an R2 read file...
        elif re.search(r'R2', file.stem, re.I):
            rev_reads.append(file)  # add to a list, will sort to match its R1 file next
        # if an R1/R2 tag is not detected
        else:
            pairs_dict.update({file:''})  # add to dict with an empty str as value, will remain w/o paired file

    # if any reverse reads are detected...
    if len(rev_reads) > 0:
        for rev in rev_reads:
            # find their associated R1 file, and match these in the dict
            sample_id = re.search(r'.+?(?=_R2)', rev.stem).group(0)
            for k in pairs_dict.keys():
                if re.search(sample_id, k.stem, re.I):
                    pairs_dict[k] = rev

    return pairs_dict

def filter_out_phix(input_files, output_dir, reference_dir, kmer, hdist, keep_log, keep_removed_seqs, bbduk_mem=None):

    # create a directory to send all output files and directories to
    phix_parent = mkdir_exist_ok(new_dir=output_dir)

    # load in the configuration settings
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    # create an output directory for reads that pass the PhiX filter
    nophix_path = mkdir_exist_ok(
        new_dir=f'./{NOPHIX_PREFIX}_{run_name}',
        parent_dir=phix_parent,
    )

    # create a log file
    bbduk_log = phix_parent / f'bbduk_{run_name}.log'

    #
    # # if no memory size is provided as argument, auto-detect the maximum available memory
    # if bbduk_mem is None:
    #     bbduk_mem_unit, bbduk_mem_size = get_available_memory(memory_units='GB')
    #
    # # if memory size is provided as an argument...
    # else:
    #
    #     # ensure that bbduk_mem is an integer, and if not, then round down to nearly whole number
    #     bbduk_mem = math.floor(bbduk_mem)
    #
    #     # auto-detect the units of the input value for bbduk_mem based on likely cut-off values
    #     if bbduk_mem < 20:
    #         bbduk_mem_unit = 'g'
    #     elif 1000 < bbduk_mem <= 20:
    #         bbduk_mem_unit = 'm'
    #     else:
    #         bbduk_mem_unit = 'b'
    #
    #     # check that the size specified by bbduk_mem is available on the system
    #     if bbduk_mem <= len(os.sched_getaffinity(0)):
    #         bbduk_mem_size = bbduk_mem
    #
    #     # if it isn't, then warn user that it exceeds what is available
    #     else:
    #         bbduk_mem_size = len(os.sched_getaffinity(0))
    #         warning_msg = (f'The requested available memory to use for bbduk ({bbduk_mem}) exceeds the available '
    #                        f'memory. Instead using the maximum available memory for this '
    #                        f'job ({bbduk_mem_size}).\n')
    #         print(warning_msg)
    #
    # # combine the memory size and unit together with the bbduk memory flag (-Xmx) into a single command string
    # bbduk_mem_str = '-Xmx' + str(bbduk_mem_size) + bbduk_mem_unit.lower()[0]
    #
    # print(f'bbduk is using {bbduk_mem_size} {bbduk_mem_unit} of memory.')


    # create pairs of R1/R2 reads in order to detect PhiX pair-wise (required)
    pairs_dict = pair_reads(input_files)

    for input_fwd, input_rev in pairs_dict.items():

        # create output file paths based on the input file basename
        # fwd (R1)
        nophix_fwd_out = add_prefix(
            file_path=input_fwd,
            prefix=NOPHIX_PREFIX,
            action=None,
            dest_dir=nophix_path,
            f_delim=settings['formatting']['filename_delim'],
            output_compressed=True,
            replace_prefix=True,
        )
        # rev (R2)
        nophix_rev_out = add_prefix(
            file_path=input_rev,
            prefix=NOPHIX_PREFIX,
            action=None,
            dest_dir=nophix_path,
            f_delim=settings['formatting']['filename_delim'],
            output_compressed=True,
            replace_prefix=True,
        )

        # create a standard bbduk command list, without any extra options like log files or keeping discarded reads
        bbduk_cmd = ['bbduk.sh',                # shell script that runs bbduk
                     f'in1={input_fwd}',        # input file - forward reads
                     f'out1={nophix_fwd_out}',  # output file - forward reads
                     f'in2={input_rev}',        # input file - reverse reads
                     f'out2={nophix_rev_out}',  # output file - reverse reads
                     'ref=phix',                # compare against the bbduk phix reference dataset
                     f'k={str(kmer)}',          # kmers to compare to; must be added as string
                     f'hdist={str(hdist)}',     # hamming distance for phix sequence comparison; must be added as string
                     # bbduk_mem_str,             # maximum available memory that bbduk can use to filter out phix
                     ]

        # check whether reads that don't pass the filter should be retained
        if keep_removed_seqs:

            # create an output directory for sequences that do *not* pass the PhiX filter
            phix_path = mkdir_exist_ok(
                new_dir=f'./{flip_prefix(NOPHIX_PREFIX)}_{run_name}',
                parent_dir=phix_parent,
            )

            # create output file paths for the PhiX reads based on the input file basename
            # fwd (R1)
            phix_fwd_out = add_prefix(
                file_path=input_fwd,
                prefix=flip_prefix(NOPHIX_PREFIX),
                action=None,
                dest_dir=phix_path,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=True,
                replace_prefix=True,
            )
            # rev (R2)
            phix_rev_out = add_prefix(
                file_path=input_rev,
                prefix=flip_prefix(NOPHIX_PREFIX),
                action=None,
                dest_dir=phix_path,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=True,
                replace_prefix=True,
            )

            # create a list of values that need to be added to the original bbduk command in order for the
            #  PhiX reads to be retained
            bbduk_keepseqs = [
                f'outm1={phix_fwd_out}',  # keep removed PhiX sequences - forward read output
                f'outm2={phix_rev_out}',  # keep removed PhiX sequences - reverse read output
            ]

            # add the outm1 / outm2 commands and output file paths for PhiX reads to the original bbduk command
            append_subprocess(
                cli_command_list=bbduk_cmd,
                options_to_add=bbduk_keepseqs,
                position=5,
                return_copy=False,
            )

        # don't alter the original bbduk command here if PhiX reads should be thrown out completely
        else:
            pass


        # if keep_log=True, then append the arguments for all bbduk summary files available
        if keep_log:
            append_subprocess(
                cli_command_list=bbduk_cmd,
                options_to_add = f'stats={bbduk_log}',
                position=-1,
                return_copy=False,
            )

        # otherwise, don't make further changes to the bbduk command
        else:
            pass

        # execute the compiled bbduk PhiX filtering command
        run_subprocess(
            bbduk_cmd,
            dest_dir=phix_parent,
            run_name=run_name,
            auto_respond=settings['automate']['auto_respond'],
        )

    return create_file_list(nophix_path)

def prefilter_fastx(input_files, output_dir, reference_dir, maxn, qmax, keep_log, keep_removed_seqs):

    ## LOAD SETTINGS ##

    # load in the configuration settings
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    ## CREATE OUTPUT FILE PATHS ##

    # create a directory to send all output files and directories to
    noambig_parent = mkdir_exist_ok(new_dir=output_dir)

    # create a parent directory to sort reads w/o ambiguous base calls into
    noambig_path = mkdir_exist_ok(
        new_dir=f'./{NOAMBIG_PREFIX}_{run_name}',
        parent_dir=noambig_parent,
    )


    ## FILTER INPUT SEQUENCES IN PAIRS ##

    # create pairs of R1/R2 reads in order to process reads in paired files (required)
    pairs_dict = pair_reads(input_files)

    for input_fwd, input_rev in pairs_dict.items():

        # determine the file format of the input sequence files; pull just the difference
        #    between .fasta and .fastq (i.e., the a or q from the string)
        input_filefmt = input_fwd.suffixes[0][-1]

        # create an output file name for sequences without ambiguous base call counts exceeding maxn
        # fwd (R1)
        noambig_fwd_out = add_prefix(
            file_path=input_fwd,
            prefix=NOAMBIG_PREFIX,
            action=None,
            dest_dir=noambig_path,
            f_delim=settings['formatting']['filename_delim'],
            output_compressed=False,
            replace_prefix=True,
        )
        # rev (R2)
        noambig_rev_out = add_prefix(
            file_path=input_rev,
            prefix=NOAMBIG_PREFIX,
            action=None,
            dest_dir=noambig_path,
            f_delim=settings['formatting']['filename_delim'],
            output_compressed=False,
            replace_prefix=True,
        )

        # compile the standard VSEARCH filtering function and options
        vsearch_filter_cmd = [
            'vsearch',                                         # call vsearch
            '--fastx_filter', input_fwd,                       # filter command - fwd reads
            '--reverse', input_rev,                            # filter command - rev reads
            f'--fast{input_filefmt}out', noambig_fwd_out,      # output file - fwd reads
            f'--fast{input_filefmt}out_rev', noambig_rev_out,  # output file - rev reads
            f'--fast{input_filefmt}_maxns', str(maxn),         # maximum number of ambiguous base calls allowed
        ]

        # if the input files are .fastq files, will need to increase the maximum quality score allowed for input
        #  files in order to avoid an error with the default maximum quality score of 41
        if input_filefmt == 'q':
            append_subprocess(
                cli_command_list=vsearch_filter_cmd,
                options_to_add=['--fastq_qmax', str(qmax)],
                position=-1,  # add to the end of the standard vsearch filter command, just after --fastx_maxns
                return_copy=False,
            )
        # if the input files are .fasta files, then no need to worry about the quality scores
        else:
            pass

        # check whether reads that don't pass the filter should be retained in a separate directory
        if keep_removed_seqs:

            # create an output directory for sequences with ambiguous base calls
            ambig_path = mkdir_exist_ok(
                new_dir=f'./{flip_prefix(NOAMBIG_PREFIX)}_{run_name}',
                parent_dir=noambig_parent,
            )

            # create an output file name for sequences *with* ambiguous base call counts exceeding maxn
            # fwd (R1)
            ambig_fwd_out = add_prefix(
                file_path=input_fwd,
                prefix=flip_prefix(NOAMBIG_PREFIX),
                action=None,
                dest_dir=ambig_path,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=False,
                replace_prefix=True,
            )
            # rev (R2)
            ambig_rev_out = add_prefix(
                file_path=input_rev,
                prefix=flip_prefix(NOAMBIG_PREFIX),
                action=None,
                dest_dir=ambig_path,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=False,
                replace_prefix=True,
            )

            # compile the list of commands to append to the filtering command so that discarded reads are saved
            vsearch_filter_keepseqs = [
                f'--fast{input_filefmt}out_discarded', ambig_fwd_out,     # output file - discarded fwd reads
                f'--fast{input_filefmt}out_discarded_rev', ambig_rev_out, # output file - discarded rev reads
            ]

            # add these additional parameters and arguments to the basic vsearch filtering command
            append_subprocess(
                cli_command_list=vsearch_filter_cmd,
                options_to_add=vsearch_filter_keepseqs,
                position=9,  # just after the output paths for the fwd/rev reads that pass the filter
                return_copy=False,
            )

        # make no additional changes to the vsearch filtering function if keep_removed_seqs=False
        else:
            pass

        # check whether a log file should be written out, which will include quality scores and read counts
        if keep_log:

            # create an output log file path
            vsearch_filter_log = noambig_parent / f'vsearch_prefilter_{run_name}.log'

            # regardless of input file format, add the --log option to the vsearch prefilter command list
            append_subprocess(
                cli_command_list=vsearch_filter_cmd,
                options_to_add=['--log', vsearch_filter_log],
                position=-1,
                return_copy=False,
            )

            # if the input files are .fastq, write out additional stats only available for .fastq files
            if input_filefmt == 'q':

                append_subprocess(
                    cli_command_list=vsearch_filter_cmd,
                    options_to_add=['--fastq_stats', vsearch_filter_log],
                    position=-1,
                    return_copy=False,
                )

            # do not add any additional flags if input files are .fasta
            else:
                pass

        # if no log is to be written out, do nothing
        else:
            pass

        # execute the vsearch filtering command that has been compiled
        run_subprocess(
            vsearch_filter_cmd,
            dest_dir=noambig_parent,
            run_name=run_name,
            program='vsearch-prefilter',
            separate_sample_output=True,
            auto_respond=settings['automate']['auto_respond'],
        )

    ## COMPRESS OUTPUT SEQUENCE FILES ##

    # compress the output written by vsearch, removing the uncompressed copies
    compress_data(
        input_path=noambig_path,    # reads without ambiguous base calls
        output_path=None,           # write to the same directory as noambig_path
        compress_fmt='gzip',        # compress using the gzip algorithm
        keep_input=False,           # remove the uncompressed input files, retaining only the .gz compressed versions
    )

    # also compress the ambiguous base call sequence files, if written to the file system
    if keep_removed_seqs:
        compress_data(
            input_path=ambig_path,  # reads with ambiguous base calls, if keep_removed_seqs=True
            output_path=None,       # write to the same directory as ambig_path
            compress_fmt='gzip',    # comrpess using the gzip algorithm
            keep_input=False,       # remove the uncompressed input files, retaining only the .gz compressed versions
        )
    else:
        pass

    return noambig_path

def identify_primers(platform, config_dict, verbose=True):
    primer_dict = config_dict['primers']

    target_primers = {}
    primer_names = []
    # for d in primer_dict.keys():
    #     target_primers[d] = primer_dict[d]['sequence'][platform]
    #     primer_names.append(primer_dict[d]['name'][platform])
    for orient in ['fwd', 'rev']:
        for name, seq in primer_dict[orient][platform].items():
            primer_names.append(name)
            target_primers.update({orient:seq})

    # add the reverse complement of each primer to this dict
    target_primers['fwd_rc'] = str(Seq(target_primers['fwd']).reverse_complement())
    target_primers['rev_rc'] = str(Seq(target_primers['rev']).reverse_complement())

    if verbose:
        print(f'Searching for {primer_names[0]} and {primer_names[1]} in {platform.title()} reads...\n')

    return target_primers

def confirm_no_primers(input_files, reference_dir, platform):

    # if directory path provided, create generator of file paths
    if input_files.is_dir():
        input_files = input_files.glob(SEQ_FILE_GLOB)

    # read in settings from the configuration file
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    # create a regex from fwd, rev, fwd_rc, and rev_rc primer sequences, to confirm no primer remains
    primer_dict = identify_primers(platform, config_dict=settings, verbose=False)  # read in primer dict for this platform
    # join the primer seq strings together with the 'or' pipe for regex search
    primer_re = '|'.join(list(primer_dict.values()))  # will search for any of these primers/orients

    # record read lengths to an output .json file
    read_info_json = (file_map['pipeline-output']['primers-trimmed'] / f'{run_name}_read-info').with_suffix('.json')

    # record primer detections to a .json log file
    primer_detected_json = (file_map['pipeline-output']['primers-trimmed'] / f'{run_name}_detected-primers').with_suffix('.json')

    # create an empty dictionary, whose contents will later be writen out to .json log
    read_info_dict = {}  # for recording read count and read lengths per sample
    p_detect_dict = {}  # for recording detection of primers in seqs
    primers_detected = 0  # keep track of how many reads had primers detected
    err_samples = set()  # keep track of how many samples had primers detected
    total_samples = 0  # counter for total number of samples, since input_files is a generator, can't get len later
    # go through each input file
    for file in input_files:
        total_samples += 1

        # get the sample ID and add it to the .json dict as key, with empty values for read count + len
        sample_id = get_sample_id(file, platform=platform)

        # get the read orientation, i.e., R1 is fwd and R2 is rev
        orient = get_read_orient(file)

        # check if there are entries yet for these values; if not, add
        if sample_id in read_info_dict.keys():  # if the sample_id is already in dict...
            if orient in read_info_dict[sample_id]:  # and the orientation is too...
                update = True  # set update var to True to inform how to add info for each read
            else:  # if the sample_id is in dict, but the orient is not...
                # add a nested dict for that orient
                read_info_dict[sample_id].update({orient: {}})
        else:  # if the sample_id is NOT in dict...
            # add it to the dict, with the subdict for the orientation
            read_info_dict[sample_id] = {orient: {}}

        # create an empty list to append read lengths to for this sample
        # easier to append vals to an empty list than an empty array
        read_lens = []

        # go through each read in this sample...
        for seq in SeqIO.parse(file, 'fastq'):

            # append the sequence length to the .json output
            read_lens.append(len(seq))

            # confirm that the primers in all orientations are not detected in the read
            found = re.search(primer_re, str(seq), re.I)

            if found:  # if a primer is detected..

                primers_detected += 1  # add to error reads counter
                err_samples.add(sample_id)  # add error sample id to err sample set

                # first, create dict for just this read
                err_read_dict = {seq.id.split(':')[-1]:{'seq_id': seq.id,  # full read ID
                                                        'primer_id': [v[0] for v in primer_dict.items() if v[1] == found.group()][0],  # whether fwd/rev/fwd_rc/rev_rc primer
                                                        'primer_seq': found.group(),  # sequence of detected primer
                                                        'pos_start': found.span()[0],  # pos of primer start in seq
                                                        'pos_end': found.span()[1],  # pos of primer end in seq
                                                        'read_len': len(seq)}  # total len of read
                                 }

                # then check whether to create new entry, or update existing one
                # if this sample already has an entry in the primer detect dict...
                if sample_id in p_detect_dict.keys():
                    # check if this read orientation (R1/R2) is already in primer detect dict...
                    if orient in p_detect_dict[sample_id]:
                        # if it is, update with the error info for this read
                        p_detect_dict[sample_id][orient].update(err_read_dict)
                    else:
                        # if not, then add as new entry
                        p_detect_dict[sample_id][orient] = err_read_dict
                # if sample not yet in dict, add sample_id and orient together
                else:
                    p_detect_dict[sample_id] = {orient: err_read_dict}

        # calculate summary information for each sample based on the lengths of each read
        read_len_arr = np.array(read_lens)

        # calculate all read metrics for this sample based on this array
        # cannot have np values written to JSON, so convert to regular non-numpy int/float
        sample_read_summary = {'read_count': int(read_len_arr.shape[0]),
                               'read_len_mean': float(np.mean(read_len_arr).round(2)),
                               'read_len_std': float(np.std(read_len_arr).round(2)),
                               'read_len_min': int(np.min(read_len_arr)),
                               'read_len_max': int(np.max(read_len_arr)),
                               'read_len_q25': float(np.quantile(read_len_arr, 0.25).round(2)),
                               'read_len_q50': float(np.quantile(read_len_arr, 0.50).round(2)),
                               'read_len_q75': float(np.quantile(read_len_arr, 0.75).round(2)),
                               'read_len_q100': float(np.quantile(read_len_arr, 1).round(2))
                               }

        # add this to the overall read dictionary
        read_info_dict[sample_id][orient] = sample_read_summary

    # dump the contents of both dictionaries to their output files
    # always write out the sequence length
    with open(read_info_json, 'w+') as info_out:
        json.dump(read_info_dict, info_out)

    # check whether there are errors to write out, then write out and print summary
    if len(p_detect_dict) > 0:
        perc_samples = (len(err_samples) / (total_samples/2))*100  # divide by 2 to account for R1/R2
        print(f'WARNING. {primers_detected} reads from {len(err_samples)} samples ({perc_samples:.1}% of all samples) '
              f'contained a primer sequence after trimming primers with cutadapt. Please consult the '
              f'details in the summary file {primer_detected_json} for details.\n')
        with open(primer_detected_json, 'w+') as err_out:
            json.dump(p_detect_dict, err_out)

    return None

def remove_primers(input_files, output_dir, reference_dir, platform, paired_end, keep_removed_seqs, max_error=None, max_untrimmed=None, linked_adapters=None, require_adapters=None, keep_log=True):

    ## IMPORT CONFIGURATION SETTINGS ##########################################

    # read in settings from the configuration file
    settings = get_settings(reference_dir)

    # get the bioinformatics run name for tagging output directories and files
    run_name = settings['run_details']['run_name']

    ## INPUT / OUTPUT DIRECTORIES #############################################

    # make the main directory for output for primer trimming
    trim_parent = mkdir_exist_ok(new_dir=output_dir)

    # create a subdirectory for the primer-trimmed sequence output
    trim_path = mkdir_exist_ok(
        new_dir=f'./{TRIMMED_PREFIX}_{run_name}',
        parent_dir=trim_parent,
    )

    ## PRIMER SEQUENCES #######################################################

    # get the forward and reverse primers
    primer_dict = identify_primers(platform, config_dict=settings)
    fwd_primer = primer_dict['fwd']
    rev_primer = primer_dict['rev']


    ## PAIRED-END PRIMER TRIMMING #############################################

    # compose and execute the command line command to run cutadapt for paired-end reads
    if paired_end:

        # create a dictionary that will link the forward (R1) and reverse (R2) sequence files for each given sample
        paired_dict = pair_reads(input_files)

        for fwd_seqs_in, rev_seqs_in in paired_dict.items():

            # forward (R1) read: create an output path for the primer-trimmed sequences from this sample
            fwd_seqs_out = add_prefix(
                file_path=fwd_seqs_in,
                prefix=TRIMMED_PREFIX,
                dest_dir=trim_path,
                action=None,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=True,
                replace_prefix=True,
            )

            # reverse (R2) read: create an output path for the primer-trimmed sequences from this sample
            rev_seqs_out = add_prefix(
                file_path=rev_seqs_in,
                prefix=TRIMMED_PREFIX,
                dest_dir=trim_path,
                action=None,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=True,
                replace_prefix=True,
            )
            # compose the standard cutadapt command to use, before incorporating additional options
            cutadapt_cmd = [
                'cutadapt',
                '-a', fwd_primer,                   # forward primer: search for fwd primer on 3' end of forward (R1) reads
                '-A', rev_primer,                   # reverse primer: search for rev primer on 3' end of reverse (R2) reads
                '--revcomp',                        # search each read's reverse complement as well
                '--cores', str(0),                  # auto-detect the number of available CPUs to use
                '--output', fwd_seqs_out,           # forward (R1) output sequence file
                '--paired-output', rev_seqs_out,    # reverse (R2) output sequence file
                fwd_seqs_in,                        # forward (R1) input sequence file
                rev_seqs_in,                        # reverse (R2) input sequence file
            ]

            # if the untrimmed / discarded reads are to be retained...
            if keep_removed_seqs:

                # make an output directory for all untrimmed / discarded reads
                notrim_path = mkdir_exist_ok(
                    new_dir=f'./{flip_prefix(TRIMMED_PREFIX)}_{run_name}',
                    parent_dir=trim_parent,
                )

                # forward (R1) read: create an output path for the untrimmed / discarded reads from this sample
                untrim_fwd_seqs_out = add_prefix(
                    file_path=fwd_seqs_in,
                    prefix=flip_prefix(TRIMMED_PREFIX),
                    dest_dir=notrim_path,
                    action=None,
                    f_delim=settings['formatting']['filename_delim'],
                    output_compressed=True,
                    replace_prefix=True,
                )

                # reverse (R2) read: create an output path for the untrimmed / discarded reads from this sample
                untrim_rev_seqs_out = add_prefix(
                    file_path=rev_seqs_in,
                    prefix=flip_prefix(TRIMMED_PREFIX),
                    dest_dir=notrim_path,
                    action=None,
                    f_delim=settings['formatting']['filename_delim'],
                    output_compressed=True,
                    replace_prefix=True,
                )

                # add the cutadapt option and output file paths to the standard cutadapt command
                append_subprocess(
                    cli_command_list=cutadapt_cmd,
                    options_to_add=[
                        '--untrimmed-output', untrim_fwd_seqs_out,          # untrimmed forward (R1) reads
                        '--untrimmed-paired-output', untrim_rev_seqs_out,   # untrimmed reverse (R2) reads
                    ],
                    position=-1,
                    return_copy=False,
                )

            # if the discarded / untrimmed reads are to be removed completely...
            else:

                # you need to specify that they are not wanted, or else will be written to the trimmed output too
                append_subprocess(
                    cli_command_list=cutadapt_cmd,
                    options_to_add=['--discard-untrimmed'],
                    position=-1,
                    return_copy=False,
                )


            run_subprocess(
                cutadapt_cmd,
                dest_dir=trim_parent,
                separate_sample_output=True,
                run_name=run_name,
                auto_respond=settings['automate']['auto_respond'],
            )

    ## SINGLE READ PRIMER TRIMMING ############################################

    # compose and execute the command line command to run cutadapt for single-end reads
    else:

        # for unpaired sequences, process one sequence file at a time
        for seqs_in in input_files:

            # create an output path for the primer-trimmed sequences
            seqs_out = add_prefix(
                file_path=seqs_in,
                prefix=TRIMMED_PREFIX,
                dest_dir=trim_path,
                action=None,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=True,
                replace_prefix=True,
            )

            # compose the standard cutadapt command to use, before incorporating additional options
            cutadapt_cmd = [
                'cutadapt',                                         # call cutadapt
                '-a', f'^{fwd_primer}...{rev_revcomp_primer}$',     # search 5' end for anchored, linked primers
                '-a', f'^{rev_primer}...{fwd_revcomp_primer}$',     # search 3' end for anchored, linked primers
                '--revcomp',                                        # search the reverse complement of each read as well
                '-n', '2',                                          # do two passes over each read
                '-e', str(max_err),                                 # maximum expected error
                '--cores', str(0),                                  # auto-detect the number of available CPUs to use
                '-o', seqs_out,                                     # output file path
                seqs_in,                                            # input file path
            ]

            # if the untrimmed / discarded reads are to be retained...
            if keep_removed_seqs:

                # make an output directory for all untrimmed / discarded reads
                notrim_path = mkdir_exist_ok(
                    new_dir=f'./{flip_prefix(TRIMMED_PREFIX)}_{run_name}',
                    parent_dir=trim_parent,
                )

                # create an output path for the untrimmed / discarded reads from this sample
                untrim_seqs_out = add_prefix(
                    file_path=seqs_in,
                    prefix=flip_prefix(TRIMMED_PREFIX),
                    dest_dir=notrim_path,
                    action=None,
                    f_delim=settings['formatting']['filename_delim'],
                    output_compressed=True,
                    replace_prefix=True,
                )

                # add the cutadapt option and output file paths to the standard cutadapt command
                append_subprocess(
                    cli_command_list=cutadapt_cmd,
                    options_to_add=['--untrimmed-output', untrim_seqs_out],
                    position=-1,
                    return_copy=False,
                )

            # if the discarded / untrimmed reads are to be removed completely...
            else:

                # you need to specify that they are not wanted, or else will be written to the trimmed output too
                append_subprocess(
                    cli_command_list=cutadapt_cmd,
                    options_to_add=['--discard-untrimmed'],
                    position=-1,
                    return_copy=False,
                )


            # create or append to cutadapt stderr/stdout, all output goes to one file (not one per sample)
            run_subprocess(
                cutadapt_cmd,
                dest_dir=trim_parent,
                separate_sample_output=True,
                run_name=run_name,
                auto_respond=settings['automate']['auto_respond'],
            )

    return trim_path

def merge_reads(input_files, output_dir, reference_dir, compress_output, keep_removed_seqs, keep_log):

    ## LOAD SETTINGS ##

    # load in the configuration settings
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    ## CREATE OUTPUT FILE PATHS ##

    # create a directory to send all quality filtered output files and directories to
    merged_parent = mkdir_exist_ok(
        new_dir=output_dir,
    )

    # create a parent directory to sort reads into that pass quality filtering
    merged_path = mkdir_exist_ok(
        new_dir=f'./{MERGED_PREFIX}_{run_name}',
        parent_dir=merged_parent,
    )

    ## PAIR INPUT FILES ##

    # if the input files is already a paired dictionary
    if isinstance(input_files, dict):

        # rename the variable to paired_dict
        paired_dict = input_files

    # if the input files are not paired by forward (R1) and reverse (R2) reads...
    else:

        # use the paired_reads() function to pair the forward and reverse reads sequence file
        paired_dict = pair_reads(input_files)

    ## COMPILE MERGED SEQUENCE FILE PATHS ##

    # create an empty list to add merged output sequence file paths to
    merged_output_files = []

    ## MERGE PAIRED-END SEQUENCES ##

    # for each pair of forward (R1) and reverse (R2) reads
    for input_fwd, input_rev in paired_dict.items():

        ## CREATE OUTPUT MERGED FILE ##

        # get the basename of the input file (i.e., without R1 / R2 tag)
        input_basename = input_fwd.parent / re.sub(r'(?<=\d{2})_(R1|R2)(?=\.)', r'', input_fwd.name)

        # create an output file for this sample's merged reads
        output_merged = add_prefix(
            file_path=input_basename,
            prefix=MERGED_PREFIX,
            dest_dir=merged_path,
            action=None,
            f_delim=settings['formatting']['filename_delim'],
            output_compressed=False,
            replace_prefix=True,
        )

        # add the output file path to the list of output file paths for all samples
        merged_output_files.append(output_merged)

        ## BASIC VSEARCH MERGE COMMAND ##

        # assemble to basic vsearch merging command list
        vsearch_merge_cmd = [
            'vsearch',                          # call the vsearch program
            '--fastq_mergepairs', input_fwd,    # the input forward (R1) read sequence file
            '--reverse', input_rev,             # the input reverse (R2) read sequence file corresponding to input fwd (R1)
            '--fastqout', output_merged,        # the output file path of the merged forward (R1) and reverse (R2) reads
        ]

        ## KEEP UNMERGED READS? ##

        # assess whether the reads that could not be merged should be retained
        if keep_removed_seqs:

            # create an output directory for unmerged sequence file output across all samples
            unmerged_path = mkdir_exist_ok(
                new_dir=f'./{flip_prefix(MERGED_PREFIX)}_{run_name}',
                parent_dir=merged_parent,
            )

            # create an output file path for unmerge forward (R1) reads
            output_unmerged_fwd = add_prefix(
                file_path=input_fwd,
                prefix=flip_prefix(MERGED_PREFIX),
                dest_dir=unmerged_path,
                action=None,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=False,
                replace_prefix=True,
            )

            # create an output file path for unmerge forward (R1) reads
            output_unmerged_rev = add_prefix(
                file_path=input_rev,
                prefix=flip_prefix(MERGED_PREFIX),
                dest_dir=unmerged_path,
                action=None,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=False,
                replace_prefix=True,
            )

            append_subprocess(
                cli_command_list=vsearch_merge_cmd,
                options_to_add=[
                    '--fastqout_notmerged_fwd', output_unmerged_fwd,  # output file for forward (R1) unmerged reads
                    '--fastqout_notmerged_rev', output_unmerged_rev,  # output file for reverse (R2) unmerged reads
                ],
                position=-1,
                return_copy=False,
            )

        # do not add any additional vsearch merge options at this time
        else:
            pass

        ## WRITE VSEARCH LOG OUT TO FILE? ##

        # if a log file should be written...
        if keep_log:

            # create an output file path for the vsearch merge log file
            vsearch_log_file = merged_parent / f'vsearch_merge_{run_name}.log'

            # append the vsearch log option to the output file
            append_subprocess(
                cli_command_list=vsearch_merge_cmd,
                options_to_add=['--eetabbedout', vsearch_log_file],
                position=-1,
                return_copy=False,
            )

        # if a log file should not be written, do nothing
        else:
            pass

        ## EXECUTE FINAL VSEARCH MERGE COMMAND ##

        # run the final vsearch merge command
        run_subprocess(
            cli_command_list=vsearch_merge_cmd,
            dest_dir=merged_parent,
            run_name=run_name,
            program='vsearch-merged',
            separate_sample_output=True,
            auto_respond=settings['automate']['auto_respond'],
        )

    ## COMPRESS OUTPUT FILES? ##

    # if the output files should be compressed...
    if compress_output:

        # compress the uncompressed output files and return the list of compressed file paths
        output_compressed = compress_data(
            input_path=merged_path,
            output_path=None,
            compress_fmt='gzip',
            keep_input=False,
        )

        # if discarded sequences were also written out to the file system...
        if keep_removed_seqs:

            # compress these sequences as well
            output_unmerged_compressed = compress_data(
                input_path=unmerged_path,
                output_path=None,
                compress_fmt='gzip',
                keep_input=False,
            )

        # if there are no discarded sequences, no additional compression needed
        else:
            pass

        # return only the compressed output file paths
        return output_compressed

    # if the output files are not to be compressed...
    else:

        # return just the original list of files
        return merged_output_files

def quality_filter(input_files, output_dir, platform, reference_dir, max_qscore, min_qscore, min_len, max_len, max_error, merge, keep_removed_seqs, keep_log, compress_output):


    ## LOAD SETTINGS ##

    # load in the configuration settings
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']


    ## CREATE OUTPUT FILE PATHS ##

    # create a directory to send all quality filtered output files and directories to
    qfilt_parent = mkdir_exist_ok(
        new_dir=output_dir,
    )

    # create a parent directory to sort reads into that pass quality filtering
    qfilt_path = mkdir_exist_ok(
        new_dir=f'./{QUALFILT_PREFIX}_{run_name}',
        parent_dir=qfilt_parent,
    )


    ## CREATE LIST OF INPUT FILES FROM PAIRED-END READS ##

    if platform == 'illumina':

        # pair togetherthe R1/R2 input files for each sample from the input files
        paired_dict = pair_reads(input_files)

        ## MERGE PAIRED-END READS ##

        if merge:

            # use the merge_reads() function to first merge forward (R1) and reverse (R2) input sequences
            qfilt_seqs_in = merge_reads(
                input_files=paired_dict,
                output_dir=qfilt_parent,
                reference_dir=reference_dir,
                compress_output=compress_output,
                keep_removed_seqs=keep_removed_seqs,
                keep_log=keep_log,
            )

        ## QUALITY FILTER ONLY FORWARD READS ##

        # if these paired-end illumina sequence files are not going to be merged...
        else:

            # create a list of input sequence files to process that are only the forward reads
            qfilt_seqs_in = list(paired_dict.keys())


    ## CREATE LIST OF INPUT FILES FROM SINGLE-END READS ##

    # if the platform isn't illumina, then it is sanger or pacbio, both of which are single-end reads
    else:

        # just rename the input file list so that it matches the name used for illumina sequences after pre-processing
        qfilt_seqs_in = input_files


    ## QUALITY FILTER INPUT FILES ##

    # keep track of how many input files were .fasta formatted and couldn't be filtered beyond lengths
    input_is_fasta = 0

    # add all output files paths to a dictionary of output file paths, sorted by quality filtered or unfiltered (if any)
    output_files = {dest:[] for dest in ['quality_filtered', 'unfiltered']}

    # iterate through each file in the input file list
    for seq_in in qfilt_seqs_in:

        ## CREATE OUTPUT FILE PATH ##

        # determine the file format of the input sequence files; pull just the difference
        #    between .fasta and .fastq (i.e., the a or q from the string)
        input_filefmt = seq_in.suffixes[0][-1]

        # create an output file path for the quality-filtered version of this input file
        seq_out = add_prefix(
            file_path=seq_in,
            prefix=QUALFILT_PREFIX,
            dest_dir=qfilt_path,
            action=None,
            f_delim=settings['formatting']['filename_delim'],
            output_compressed=False,
            replace_prefix=True,
        )

        # add this output file path to the output file dictionary in the quality filtered list of files
        output_files['quality_filtered'].append(seq_out)

        ## ASSEMBLE STANDARD VSEARCH COMMAND ##

        # assemble a standard vsearch --fastx_filter command without any additional options from params
        # only add options that are available for both .fastq and .fasta input files
        vsearch_qfilt_cmd = [
            'vsearch',                                      # call on vsearch
            '--fastx_filter', seq_in,                       # where to read the unfiltered input sequeces from
            f'--fast{input_filefmt}out', seq_out,           # where to write the filtered output sequences to
            f'--fast{input_filefmt}_maxlen', str(max_len),  # filter out any input reads above this length (bp)
            f'--fast{input_filefmt}_minlen', str(min_len),  # filter out any input reads below this length (bp)
        ]

        ## ADD ADDITIONAL VSEARCH OPTIONS IF QUALITY SCORES PRESENT ##

        # add any vsearch options available only if input files are .fastq files
        if input_filefmt == 'q':

            # .fastq options to append to the standard vsearch filtering command
            vsearch_qfilt_fastq_opts = [
                '--fastq_qmax', str(max_qscore),  # filter out any input sequences above this quality score
                '--fastq_qmin', str(min_qscore),  # filter out any input sequences below this quality score
                '--fastq_maxee', str(max_error),    # filter out any input sequences with maximum expected errors above this value
            ]

            # append the additional .fastq options to the standard vsearch filtering command
            append_subprocess(
                cli_command_list=vsearch_qfilt_cmd,
                options_to_add=vsearch_qfilt_fastq_opts,
                position=-1,
                return_copy=False,
            )


        # if the input file is a .fasta file, add to a counter keeping track of how many input files lacked quality score info
        #   and therefore could not be filtered beyond read length
        else:
            input_is_fasta += 1

        ## KEEP UNFILTERED READS ##

        if keep_removed_seqs:

            # create an output directory for sample sequence files for reads that don't pass the quality filter
            unfilt_path = mkdir_exist_ok(
                new_dir=f'./{flip_prefix(QUALFILT_PREFIX)}_{run_name}',
                parent_dir=qfilt_parent,
            )

            # create an output sequence file path for this sample's unfiltered reads
            seq_unfilt_out = add_prefix(
                file_path=seq_in,
                prefix=flip_prefix(QUALFILT_PREFIX),
                dest_dir=unfilt_path,
                action=None,
                f_delim=settings['formatting']['filename_delim'],
                output_compressed=False,
                replace_prefix=True,
            )

            # add this output file path to the output file dictionary in the unfiltered list of files
            output_files['unfiltered'].append(seq_unfilt_out)

            # add the option to write the discarded reads that don't pass the quality filtering steps to a file
            append_subprocess(
                cli_command_list=vsearch_qfilt_cmd,
                options_to_add=[f'--fast{input_filefmt}out_discarded', seq_unfilt_out],
                position=3,
                return_copy=False
            )

        ## DISCARD UNFILTERED READS ##

        else:
            pass

        ## EXECUTE FINAL VSEARCH QUALITY FILTER COMMAND FOR THIS SAMPLE ##

        run_subprocess(
            cli_command_list=vsearch_qfilt_cmd,
            dest_dir=qfilt_parent,
            run_name=run_name,
            program='vsearch-qualfilt',
            separate_sample_output=True,
            auto_respond=settings['automate']['auto_respond'],
        )

    ## REPORT NUMBER OF INPUT .FASTA FILES ##

    # if any of the input files were .fasta files...
    if input_is_fasta > 0:

        # issue a warning regarding the filtering of input files
        warnings.warn(
            f'{input_is_fasta} input files in the file path:\n'
            f'   {input_files[0].parent}\n'
            f'are .fasta files. Due to the lack of quality score information, these input files'
            f'could only be filtered based on the maximum and minimum sequence length parameters.\n'
        )

    # otherwise, issue no print
    else:
        pass

    ## COMPRESS VSEARCH OUTPUT TO GZIP FORMAT ##

    if compress_output:

        # compress the filtered output files
        seq_filt_compressed = compress_data(
            input_path=qfilt_path,
            output_path=None,
            compress_fmt='gzip',
            keep_input=False,
        )

        # if unfiltered sequence files were also created...
        if keep_removed_seqs:

            # compress the unfiltered output files
            seq_unfilt_compressed = compress_data(
                input_path=unfilt_path,
                output_path=None,
                compress_fmt='gzip',
                keep_input=False,
            )

        # otherwise, create an empty list for an unfiltered sequence file placeholder
        else:
            seq_unfilt_compressed = []

        # create a new output file dictionary with these compressed file paths
        output_files_compressed = {
            f_type:f_list for f_type, f_list in zip(
                ['quality_filtered', 'unfiltered'],
                [seq_filt_compressed, seq_unfilt_compressed]
            )
        }

        return output_files_compressed

    ## DO NOT COMPRESS VSEARCH OUTPUT TO GZIP FORMAT ##

    # do not compress any output files, leave a .fastq / .fasta
    else:
        pass

    # return a dictionary of the output file paths (filtered and unfiltered separated)
    return output_files

def dereplicate(input_files, output_dir, min_count, derep_step, platform, reference_dir):

    # create the file prefix name based on whether it is full-length (1) or region-specific (2) dereplication
    derep_prefix = DEREP_PREFIX + '0' + str(derep_step)

    if derep_step == 1:
        out_tag = 'derep-full-length'
    else:
        out_tag = 'derep-subregions'

    # read in settings from the configuration file
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    # create main file paths for dereplicated read output

    # main output directory, provided as the output directory in the pipeline/ script
    derep_parent = mkdir_exist_ok(new_dir=output_dir)

    # subfolder within the main directory for this specific bioinformatics run
    derep_path = mkdir_exist_ok(
        new_dir=f'./{out_tag}_{run_name}',
        parent_dir=derep_parent,
    )

    # for illumina sequences, if reads were not merged, then dereplicate only the forward sequence files (R1 suffix)
    if (platform == 'illumina') or (platform == 'its1') or (platform == '18s'):

        # check if the sequences were merged previously
        premerged = settings['quality_filtering']['merge_reads'][platform]

        # if the sequences weren't merged, then the input files need to be filtered to only use the R1 (fwd) seqs
        if not premerged:
            input_files = [f for f in input_files if re.search(r'R1', f, re.I)]
            # rename files to be just sample ID or keep R1 to be clear they're forward reads only?

        # if already merged, then can just continue with the pool of merged seq files
        else:
            pass

    # if not illumina sequences, then no reason to check if merged
    else:
        pass

    # for each of the files (typically derep01) or directories (typically derep02)...
    for file in input_files:

        # the derep02 for pacbio sequences will be directories, one per sample, if/when run immediately after itsx
        if file.is_dir():

            # create a regex that will match any of the possible output region / subregion file suffixes, post-itsx
            subregion_suffix_re = '|'.join(list(POST_ITSX_SUFFIXES.values()))

            # create a list of region / subregion files to dereplicate based on this regex
            subregion_files = [f for f in file.glob('*') if re.search(subregion_suffix_re, f.name, re.I)]

            # create a main output derep directory for this sample
            derep_output = add_prefix(
                file_path=file,
                prefix=derep_prefix,
                dest_dir=derep_path,
                action='mkdir',
            )

            # go through the list of region / subregion files for this sample and dereplicate
            for region_seqs in subregion_files:

                # create a path for the region / subregion derep sequences within the main sample output dir
                derep_region_output = add_prefix(
                    file_path=region_seqs,
                    prefix=derep_prefix,
                    dest_dir=derep_output,
                    action=None,
                )

                # compile the vsearch dereplication command for this region's post-itsx sequence file
                vsearch_derep_cmd = [
                    'vsearch',
                    '--derep_fulllength', region_seqs,
                    '--output', derep_region_output,
                    '--minuniquesize', str(min_count),
                    '--sizeout',
                ]

                run_subprocess(
                    vsearch_derep_cmd,
                    dest_dir=derep_parent,
                    run_name=run_name,
                    program=f'vsearch-derep{derep_prefix}',
                    auto_respond=settings['automate']['auto_respond'],
                )

        # derep01 in all cases will (should) be a list of sequence files, not directories
        else:

            # derep always outputs a fasta format, but need to provide file name as fasta format or will appear as fastq
            #   even though it isn't really fastq when you open it up
            derep_output = add_prefix(
                file_path=file,
                prefix=derep_prefix,
                dest_dir=derep_path,
                action=None,
            ).with_suffix('.fasta')

            vsearch_derep_cmd = [
                'vsearch',
                '--derep_fulllength', file,
                '--output', derep_output,
                '--minuniquesize', str(min_count),
                '--sizeout',
            ]

            run_subprocess(
                vsearch_derep_cmd,
                dest_dir=derep_parent,
                run_name=run_name,
                program=f'vsearch-derep{derep_prefix}',
                auto_respond=settings['automate']['auto_respond'],
            )

    return None

def separate_subregions(input_files, output_dir, reference_dir, cpus=None, verbose=False):

    # import configuration settings
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    # create a directory for all ITSx output, if one does not exist
    itsx_parent = mkdir_exist_ok(new_dir=output_dir)

    # create a directory within the main ITSx output for this particular pipeline run
    itsx_path = mkdir_exist_ok(new_dir=f'./{ITSX_PREFIX}_{run_name}',  parent_dir=itsx_parent)

    # if no cpu count provided as argument, auto-detect cpus available
    if cpus is None:
        # auto-detect the number of available CPUs that itsx can use (otherwise will default to 1)
        cpu_use_count = len(os.sched_getaffinity(0))  # i thought os.process_cpu_count() was the correct one here, but not attribute
    else:
        # check that the specified number of cpus is available
        if cpus <= len(os.sched_getaffinity(0)):
            cpu_use_count = cpus
        else:
            cpu_use_count = len(os.sched_getaffinity(0))
            warning_msg = (f'The requested number of CPUs to use for ITSx ({cpus}) exceeds the available '
                           f'number of CPUs. Instead using the maximum available CPUs for this '
                           f'job ({cpu_use_count}).\n')
            print(warning_msg)

    print(f'ITSx is using {cpu_use_count} CPUs.')

    # run ITSx for each input file
    for file in input_files:

        # create a new directory for each sample, since multiple files per sample are produced
        itsx_sample_dir = add_prefix(file_path=file, prefix=ITSX_PREFIX, dest_dir=itsx_path, action='mkdir')

        # construct a base name that ITSx will use for the output files
        itsx_output_basename = add_prefix(file_path=file, prefix=ITSX_PREFIX,
                                          dest_dir=itsx_sample_dir, action=None).with_suffix('')

        # construct the ITSx command; will produce more output if verbose=True (defaults to False)
        if verbose:
            itsx_command = ['ITSx', '-i', file, '-o', itsx_output_basename, '-t', 'fungi', '--multi_thread', 'T',
                            '--save_regions', 'all', '--cpu', str(cpu_use_count)]
        else:
            itsx_command = ['ITSx', '-i', file, '-o', itsx_output_basename, '-t', 'fungi', '--multi_thread', 'T',
                            '--graphical', 'F', '--positions', 'F', '--silent', 'T',
                            '--save_regions', '{ITS1,5.8S,ITS2,LSU}', '--cpu', str(cpu_use_count)]

        # run the ITSx command
        run_subprocess(itsx_command, dest_dir=itsx_parent, run_name=run_name,
                       auto_respond=settings['automate']['auto_respond'])

    # WITH MULTI-PROCESSING
    # itsx_command_list = []
    # for file in input_files:
    #
    #     # create a new directory for each sample, since multiple files per sample are produced
    #     itsx_sample_dir = add_prefix(file_path=file, prefix=ITSX_PREFIX, dest_dir=itsx_path, action='mkdir')
    #
    #     # construct a base name that ITSx will use for the output files
    #     itsx_output_basename = add_prefix(file_path=file, prefix=ITSX_PREFIX,
    #                                       dest_dir=itsx_sample_dir, action=None).with_suffix('')
    #
    #     # construct the ITSx command; will produce more output if verbose=True (defaults to False)
    #     if verbose:
    #         itsx_command = ['ITSx', '-i', file, '-o', itsx_output_basename, '-t', 'fungi', '--multi_thread', 'T',
    #                         '--save_regions', 'all']
    #     else:
    #         itsx_command = ['ITSx', '-i', file, '-o', itsx_output_basename, '-t', 'fungi', '--multi_thread', 'T',
    #                         '--graphical', 'F', '--positions', 'F', '--silent', 'T',
    #                         '--save_regions', '{ITS1,5.8S,ITS2,LSU}']
    #
    #     itsx_command_list.append(itsx_command)
    #
    # # set variables for the run_subprocess function here, don't kknow how to add multiple to .map
    # num_samples = len(itsx_command_list)
    # itsx_parent_repeats = [itsx_parent] * num_samples
    # auto_responses = [settings['automate']['auto_respond']] * num_samples
    #
    # if __name__ == '__main__':
    #     print(f'Number of cores available: {cpu_count()}')
    #
    #     p = Process(target=run_subprocess,
    #                 args=(itsx_parent_repeats, auto_responses, itsx_command_list))
    #     p.start()
    #     p.join()

    # return the ITSx output path for this sequencing run
    return itsx_path

def concat_regions(dir_path, reference_dir, platform=None, regions_to_concat=['ITS1', '5_8S', 'ITS2'], verbose=True, header_delim=None):
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
    :param reference_dir: directory from which to start searching for the bioinformatics settings
    .toml file, typically the absolute file path of the script in which this function is called
    :param platform: the sequencing platform used in the sample file name or otherwise the region
    identifier used in the sample's file name
    :param regions_to_concat: an ORDERED list of regions to concatenated, with the 5' region
    first and the 3' region last; defaults to ITS regions
    :param header_delim: delimiter to use for separating the components of the header in
    the concatenated fasta output; defaults to semicolon
    :return: None, but will write out a fasta file
    '''

    ## IMPORT CONFIGURATION SETTINGS

    # import configuration settings as a dictionary
    settings = get_settings(reference_dir)

    # get bioinformatics run name
    run_name = settings['run_details']['run_name']

    # check whether to use header_delim from settings or directly from function parameter

    # if header_delim set to default...
    if header_delim is None:

        # then import from the configuration file
        header_delim = settings['formatting']['header_delim']

    # if header_delim is not the default, use whatever value is assigned to it
    else:
        pass

    # created a variable to switched to True if a region cannot be located without an ITSx issue
    #  moved this further down, I'm not sure if its alright down there but we'll see
    # error_occurred = False

    # create an output file using the run name for this bioinformatics run; do so here so you can remove
    err_log_path = (dir_path.parent / f'{run_name}').with_suffix('.no_concat.err')


    # determine whether the input directory contains sample subdirectories or if its a sample directory

    # create list of subdirectories within input directory
    child_dirs = [child.is_dir() for child in dir_path.glob(f'{ITSX_PREFIX}*')]

    # if the input directory contains only subdirectories, then create list of subdirectories to loop through
    if all(child_dirs):
        input_dirs = list(dir_path.glob('*'))

    # if none of the contents of input directory are directories, then its a sample dir and needs to be cast as
    #  list so that you can still use the loop method to work through
    elif not all(child_dirs):
        input_dirs = [dir_path]

    # if there's a mix of directories and not directories within the input path, then assume that the input provided
    #  is a sample directory, and look for ITSx output within the input directory but warn that errors might follow
    else:
        print(f'WARNING. Encountered a mix of directories and non-directories in the input for the '
              f'concat_regions() function from climush.bioinfo. Proceeding by searching the input '
              f'directory for the fasta files to concatenate but may encounter an error because the '
              f'function does not know whether subdirectories contain per-sample fasta files or if '
              f'the subdirectories are unrelated.\n')
        input_dirs = [dir_path]

    # create dictionary for acceptable input regions
    accepted_regions = {r'(^S.+?(?=\bs))|(?<!-)SSU': 'SSU',  # could be better but not going to obsess
                        r'ITS.?1': 'ITS1',  # checks for possible punctuation within
                        r'5.?8S': '5.8S',  # checks for possible punctuation within
                        r'ITS.?2': 'ITS2',  # checks for possible punctuation within
                        r'(^L.+?(?=\bs))|(?<!-)LSU|(?<=\.)LSU(?=\.)': 'LSU',  # similar to SSU
                        r'(full)?.?ITS(?!.?\d)': 'ITS'}  # may include full, but can't have number after (w or w/o punct)

    # check input list for valid entries, and then format string to consistent format
    # returns dict now with formatted string as key, regex as value
    def create_region_dict(region_list):
        output_dict = {}
        for region_str in region_list:
            for r in accepted_regions.keys():  # compare input region to each regex in dict
                if re.search(r, region_str, re.I):  # if it does match, stop loop and return the formatted version
                    output_dict.update({accepted_regions[r]: r})
                    break
                else:
                    continue
            else:
                return print(f'FAILURE. Did not recognize the input region {region_str}. Please choose from the '
                             f'recognized regions: {list(accepted_regions.values())}\n')
        return output_dict

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
                for r in regex.split('|'):
                    indiv_file_match = [file for file in dir_path.glob('*') if re.search(r, file.name, re.I)]
                    if len(indiv_file_match) == 1:
                        return indiv_file_match[0]
                    else:
                        continue
                print(f'WARNING. Multiple files matching the provided regex {regex} were located in '
                      f'{dir_path.name}. Returning only the first match.\n')  # REPLACE WITH OPTIONS?
                return found_file_paths[0]
        elif len(found_file_paths) == 0:
            # print(f'FAILURE. No files matching the provided regex {regex} were located in '
            #       f'{dir_path.name}.\n')
            return None
        else:
            return found_file_paths[0]

    # confirm all regions data was acquired, then combine (used below in loop)
    def join_if_all_present(pulled_data_dict):

        # create an empty dictionary to add read information to once it is checked
        joined_data = {'seq':'',
                       'id':[]}

        # for each value in details_to_include from the configuration file settings...
        for k in pulled_data_dict.keys():

            # check if the number of values for this info group is equal to the number of regions in regions_to_concat
            if len(pulled_data_dict[k]) == len(regions_to_concat):

                # if this is the key for the Seq object (i.e, a list of sequences from each region to concat)...
                if any(isinstance(i, Seq) for i in pulled_data_dict[k]):

                    # then join all the sequences in this list (i.e., actually concatenate the regions to one read)
                    joined_data['seq'] = Seq('').join(pulled_data_dict[k])  # must join with empty seq

                # if this is not the key for a Seq object (i.e., not a sequence)...
                else:

                    # format the combined information for each category based on which category it is;
                    #  this format will directly affect how the information will appear in the fasta file headers

                    # if it is the read length, then sum all of the values and confirm the concat length is correct
                    if re.search(r'len$', k, re.I):

                        # compare expected and observed read length of the concatenated read before adding info to header
                        exp_rl = sum([int(x) for x in pulled_data_dict[k]])  # calc expected length based on sum of each region length
                        obs_rl = len(joined_data['seq'])  # calc observed length based on the length of the concatenated read

                        # if the expected and observed read lengths for the concatenated read match...
                        if exp_rl == obs_rl:

                            # then create a header component for read length
                            joined_data['id'].append(f'{concat_name} {obs_rl}bp')

                        else:
                            # if these values don't match, then print an error and return None (i.e., don't concat)
                            msg = f'ERROR. The sum of the values in the headers for the post-ITSx region files does not '\
                                  f'match the length of the concatenated read. Please check the ITSx output for the '\
                                  f'following regions:\n'\
                                  f'   read ID = {read}\n' \
                                  f'   concat regions = {regions_to_concat}\n' \
                                  f'   expected length (nt) = {exp_rl}\n' \
                                  f'   concat length (nt) =   {obs_rl}\n'
                            return print(msg)

                    # if this key is the read count for full-length reads (number of identical full-length reads/sample)
                    elif re.search(r'^full\-len|count', k, re.I):

                        # if all read counts are identical among the regions (they should be)...
                        if len(set(pulled_data_dict[k])) == 1:

                            # then create the header component for read count
                            joined_data['id'].append(f'{k}={int(pulled_data_dict[k][0])}')  # if all same, doesn't matter which you chose

                        # there's likely some larger error at play if a given read doesn't match in read count among regions
                        else:
                            msg = f'ERROR. The read counts for the regions of this read are different lengths:\n'\
                                  f'  {pulled_data_dict[k]}\n'
                            return print(msg)

                    # if this key is not recognized as a default header key, then make generic version to add
                    else:
                        joined_data['id'].append(f'{k} = {pulled_data_dict[k]}')

            # if there is missing information for at least one of the regions to concatenate...
            else:
                # print an error message, and return Non
                # technically, this function should never end up running if there are missing regions for a read
                #   so if this error message prints, there is likely an issue elsewhere
                msg = f'Did not successfully pull all {k} from the regions to concatenate ' \
                      f'({len(pulled_data_dict[k])} out of {len(regions_to_concat)} were pulled). ' \
                      f'See the .no_concat.txt file for more details. \n'
                return print(msg)

        return joined_data

    # if a region is missing for a read, determine whether ITSx couldn't find or if another unknown issue
    def explain_missing_region(read_id, err_region, sample_dir):
        '''
        If a read is missing from its sample's subregion fasta file, determine whether the missing subregion
        was not detected by ITSx (an ITSx error) or if an error independent of ITSx occurred.
        :param read_id: the read ID of the sample that was missing a subregion, which is the first element
        in the fasta file header
        :param err_region: the region of the sequence that could not be located (e.g., ITS1, ITS2)
        :param sample_dir: directory containing the output from ITSx for this sample
        :return: if an ITSx error, returns 'itsx' and logs error to common .no_concat.txt file for
        the bioinformatics run that this sample is processed with and does not concatenate any
        regions for this read; if not an ITSx error, returns 'other' and prints error message; sys.exit() has
        to be carried out independent of this function in order to interrupt concat_regions()
        '''

        # determine whether an ITSx error by checking the <sample-id>.problematic.txt output file from ITSx

        # get the path to the <sample-id>.problematic.txt output from ITSx
        itsx_detect_error_path = next(sample_dir.glob('*problematic*'))

        # create generator to find line that contains this read in the problematic.txt file, yield its error message
        def find_problematic_read(err_file, read_id):
            for err in err_file.readlines():
                if err.startswith(read_id):
                    yield err.split('\t')[-1].strip()
                else:
                    continue

        # get the error message created by ITSx for this read
        with open(itsx_detect_error_path, 'rt') as err_records:
            errors = next(find_problematic_read(err_file=err_records, read_id=read_id))

        # if the missing region is located in the ITSx error, then ITSx couldn't detect the region

        # get the regex used to search files, etc. produced by ITSx for this region
        region_regex = region_search_dict[err_region]

        # search for the region causing the error in the error message for this read produced by ITSx
        if re.search(region_regex, errors, re.I):
            return 'itsx'
        # if the region is not detected in the error message from ITSx for this read, then its some other unknown issue
        else:
            return 'other'

    # create a dictionary of regions to concat, where keys are formatted regions and
    #   values are regex to use for file search
    region_search_dict = create_region_dict(regions_to_concat)

    ## PER SAMPLE ########################################################################################

    # go through each sample's directory...
    for sample_dir in input_dirs:

        print(f'sample dir = {sample_dir}\n')

        # log progress by printing the name of the sample that is being processed; will append SUCCESS to string
        #   if process completes successfully
        print(f'\n{sample_dir.name} running...', end='\r')

        # create a list of the fasta files to concat based on input regions
        fastas_to_concat = [search_path_with_regex(sample_dir, regex=r) for r in region_search_dict.values()]

        # get the sample ID from any of these fasta files
        sample_name = get_sample_id(file_path = fastas_to_concat[0], platform=platform)

        # get the file prefix used for the ITSx output files, will be used as prefix for the concatenated file output
        if platform is None:
            file_prefix = re.search(PREFIX_RE, fastas_to_concat[0].name, re.I).group(0)
        else:
            file_prefix_re = f'.+?(?=_{platform})'
            file_prefix = re.search(file_prefix_re, fastas_to_concat[0].name, re.I).group(0)

        # create a list of errors when encountered in reads; formatted to write this list to an error log file
        #  as soon as this sample is completed
        # PATCH - I keep getting the same read printing over and over, so I think once the value that checks whether
        #  an error has occurred sets to True, it isn't reseting such that every read after the error read triggers
        #  the write of the same read information to the file
        reads_with_missing_regions = set()

        # collect read IDs of reads in this sample that itsx flagged as chimeric; will skip over completely later on
        chimera_fasta = search_path_with_regex(sample_dir, regex=r'\.chimeric\.')

        # if no itsx output for chimeric reads are located for this sample, then create an empty list of chimeras
        if chimera_fasta is None:
            chimer_read_ids = []
        # if the itsx chimeric output fasta is located, extract the read IDs from this file
        else:
            with open(chimera_fasta, 'rt') as chimer_in:
                chimer_headers = [line for line in chimer_in.readlines() if line.startswith('>')]
                chimer_read_ids = [read.split(';')[0].replace('>', '') for read in chimer_headers]

        # summarize information on errors and missing regions for each sample
        reads_with_nonitsx_errors = []
        reads_without_issue = []
        reads_with_itsx_errors = []
        error_occurred = False
        regions_missing_in_sample = False
        full_read_count = 0

        # auto-assign name of joined regions for common concats, otherwise prompt user for input
        formatted_input_regions = list(region_search_dict.keys())
        if formatted_input_regions == ['ITS1', '5.8S', 'ITS2']:
            concat_name = 'full-ITS'
        elif formatted_input_regions == ['ITS1', '5.8S', 'ITS2', 'LSU']:
            concat_name = 'ITS-LSU'
        else:
            msg = f'Please provide a name to assign to this combination of subregions (use '\
                  f'hyphens instead of spaces):\n'
            default_name = '-'.join(formatted_input_regions)
            concat_name = prompt_generic(message=msg, auto_respond=settings['automate']['auto_respond'],
                                         auto_response=default_name)

        # create a list of components to include in the new fasta files within read header
        #  added these to the default settings, but not yet in the user settings
        details_to_include = settings['formatting']['fasta_headers']

        ## PER REGION ########################################################################################

        # for this sample, go through each of the itsx output files for the regions to concat, collecting info
        concat_dict = {}
        for fasta in fastas_to_concat:

            ## PER READ #######################################################################################

            for record in SeqIO.parse(fasta, 'fasta'):

                read_id = re.search(READ_ID_RENAMED_RE, record.description).group(0)

                # if this read ID matches one in ITSx's chimeric read output file, then skip to the next read
                if read_id in chimer_read_ids:
                    continue

                # get information from the post-itsx read header, which depends to an extent on a certain header format
                try:
                    region = re.search(READ_REGION_RENAMED_RE, record.description, re.I).group(0)
                    length = re.search(READ_LEN_RENAMED_RE, record.description, re.I).group(0)
                    derep01_reads = re.search(READ_COUNT_RENAMED_RE, record.description, re.I).group(0)
                    details = {k: v for k, v in zip(details_to_include, [record.seq, length, derep01_reads])}

                    if read_id in concat_dict:  # if this sample is already in the dictionary...
                        concat_dict[read_id].update({region: details})  # append the new region info to it
                    else:  # if it isn't yet in the dictionary...
                        concat_dict[read_id] = {region: details}  # add a key value pair for this sample

                # if the itsx output isn't formatted as expected, can't get important info, so print error and exit
                except AttributeError:
                    print(f'{sample_dir.name} running... ERROR!')
                    print(f'  There was an issue extracting information from the fasta header of the input file:\n '
                          f'     {fasta} \n'
                          f'  This issue occurred for read:\n'
                          f'     {read_id} \n'
                          f'  This will occur if the sequence files produced by ITSx have not yet \n'
                          f'  been reformatted. Please run the rename_read_header function and retry.\n')
                    return sys.exit(1)

        ## PER READ ############################################################################################

        # create BioSeq records for each read for each region to concatenate
        concatenated_records = []

        for read in concat_dict.keys():

            # keep track of the number of reads in this sample
            full_read_count += 1

            regions_missing_in_read = False

            # create key-value pairs where key is the info to include in BioSeq record, and empty value list to add to
            record_details = {k:[] for k in details_to_include}

            # missing region
            missing_regions = []

            ## PER REGION #######################################################################################

            # pull info from each region for this read
            for region in region_search_dict.keys():

                # first try to index the given region for this read...
                try:
                    # get the sequence and read description details for this region
                    region_dict = concat_dict[read][region]

                    # add the read description to the record_details dictionary for this read
                    for det in details_to_include:
                        if det in record_details:
                            record_details[det].append(region_dict[det])  # append before sum to ensure all regions there
                        else:
                            record_details[det] = region_dict[det]  # append before sum to ensure all regions there

                # if the region doesn't exist for this read, determine whether an error from ITSx or unknown source
                except KeyError:

                    # determine the cause of the missing region for this read; will log accordingly within function
                    cause = explain_missing_region(read_id=read, err_region=region, sample_dir=sample_dir)

                    # if the region is missing because ITSx could not detect it within the full-length read...
                    if cause == 'itsx':

                        # change this variable to True, used for formatting the error log file
                        regions_missing_in_read = True

                        # add this region to a list of missing regions for this read; these regions will be written to
                        #  a log file once all regions to concatenate for this read have been searched
                        missing_regions.append(region)

                    # if ITSx is not the reason that this region is missing, and the reason is therefore unknown..
                    else:

                        # if ITSx can't find the end of the 5.8S region, it will not produce the ITS2 either;
                        #   however, it will not log this error in the <sample-id>.problematic.txt log file
                        #   so, if an itsx issue with 5.8S has been found for this read, chances are ITSx also could
                        #   not find the ITS2 and it isn't some extraneous error in locating the ITS2 seq for this read
                        if regions_missing_in_read:
                            missing_regions.append(region)
                        else:
                            error_occurred = True
                            print(f'{sample_dir.name} running... ERROR!')
                            error_msg = (f'ERROR. The {region} region could not be located and ITSx did not log an error '
                                         f'with this read, so some unknown issue is preventing the detection of the '
                                         f'{region} for this read. The regions were not concatenated. Please check the '
                                         f'ITSx output file for this sample to troubleshoot.\n'
                                         f'   sample: {sample_name}\n'
                                         f'   read:   {read}\n')
                            print(f'   {error_msg}')
                            # break out of this loop, i.e., stop going through regions for this read and move to next read
                            #  breaking here will prevent any record from being made for this read, concat won't occur
                            # break

            ## PER READ #######################################################################################

            # if there were any regions not detected by ITSx, then write this information out to a log file
            if len(missing_regions) > 0:

                reads_with_itsx_errors.append(read)

                # format the sample ID, read ID, and non-detected regions to add to list that will append to log file
                reads_with_missing_regions.add(f'   {read} = {missing_regions}\n')

            # if a non-ITSx error occurred...
            elif error_occurred:
                reads_with_nonitsx_errors.append(read)

            # if there weren't any issues detecting the regions required to concatenate...
            else:

                combined_data = join_if_all_present(record_details)

                if (len(combined_data['id']) == 0) or (combined_data is None):
                    reads_with_nonitsx_errors.append(read)
                else:
                    reads_without_issue.append(read)
                    concat_record = SeqRecord(combined_data['seq'],
                                              name=read,
                                              id=header_delim.join([read] + combined_data['id']),
                                              description='')
                    concatenated_records.append(concat_record)

        error_reads = reads_with_nonitsx_errors + reads_with_itsx_errors
        total_errs = len(error_reads)

        if verbose:
            print(f'Total reads with errors: {total_errs} ({(total_errs / full_read_count)*100:.2f}%)\n'
                  f'   ITSx errors:     {len(reads_with_itsx_errors)}\n'
                  f'   non-ITSx errors: {len(reads_with_nonitsx_errors)}\n'
                  f'Total reads without errors: {len(reads_without_issue)}\n')

        ## PER SAMPLE #######################################################################################

        # write out the concatenated read records to a new concatenated reads fasta file
        fasta_out = sample_dir / f'{file_prefix}_{sample_name}.{concat_name}.fasta'

        SeqIO.write(concatenated_records, fasta_out, 'fasta')

        # write out the information of the reads that had regions missing and weren't concatenated
        if len(reads_with_missing_regions) > 0:
            with open(err_log_path, 'at') as fout:

                # create a sample header to better distinguish between samples in the error log
                fout.write(f'\n{sample_name}  -----------------------------\n')

                # write the name of the read and the list of missing sequence regions for this read
                for read_error in reads_with_missing_regions:
                    fout.write(read_error)

            print(f'{sample_dir.name} running... WARNING!')
            print(f'  discarded reads in which ITSx could not detect subregions.')
        else:
            if fasta_out.is_file():
                print(f'{sample_dir.name} running... success!')
            else:
                print(f'{sample_dir.name} running... ERROR!\n')
                print(f'  concatenated output fasta file not detected.')

    return None

def check_chimeras(input_files, reference_dir, output_dir, method, alpha, keep_chimeras):

    ## GET SETTINGS ###################################

    # read in the configuration settings, get bioinformatics run name from settings
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']


    ## CREATE OUTPUT DIRECTORIES ######################

    # make the main output directory for chimera check output
    chim_parent = mkdir_exist_ok(new_dir=output_dir)

    # within main chimera check output directory, create a subdirectory for the chimeric and non-chimeric reads
    nochim_path = mkdir_exist_ok(new_dir=f'./{NOCHIM_PREFIX}_{run_name}', parent_dir=chim_parent)
    chim_path = mkdir_exist_ok(new_dir=f'./{flip_prefix(NOCHIM_PREFIX)}_{run_name}', parent_dir=chim_parent)

    # create a log file for any UCHIME summary content produce by vsearch flag --uchimeout
    uchime_log = chim_parent / f'uchime_{run_name}.log'


    ## DE NOVO CHIMERA DETECTION #######################

    # if the method is set to denovo, then do de novo chimera detection
    if method == 'denovo':

        # function for de novo chimera detection depends on whether sequences have been denoised or not

        # sort input files into whether they have been denoised (have denoise_ file prefix) or not
        input_sorted = {p: [] for p in [DENOISE_PREFIX, flip_prefix(DENOISE_PREFIX)]}

        # go through each input file and sort by whether sequences in file are denoised or not
        for file in create_file_list(input_files):

            # get the file prefix of the file
            file_prefix = file.name.split('_')[0]

            # add the file path to the sorted dictionary, based on prefix
            if file_prefix == DENOISE_PREFIX:
                input_sorted[DENOISE_PREFIX].append(file)
            else:
                input_sorted[flip_prefix(DENOISE_PREFIX)].append(file)

        # run UCHIME3 de novo chimera detection on denoised sequences
        for denoised_file in input_sorted[DENOISE_PREFIX]:

            # create an output file path for the non-chimeric read output
            nochim_out = add_prefix(file_path=denoised_file, prefix=NOCHIM_PREFIX,
                                    dest_dir=nochim_path, action=None)

            # assemble the vsearch UCHIME3 command for the command line
            vsearch_denovo_cmd = ['vsearch', '--uchime3_denovo', denoised_file,
                                  '--abskew', str(alpha),
                                  '--nonchimeras', nochim_out,
                                  '--uchimeout', uchime_log]

            if keep_chimeras:

                # if configured to keep chimeras, create an output file path for chimeric reads
                chim_out = add_prefix(file_path=denoised_file, prefix=flip_prefix(NOCHIM_PREFIX),
                                      dest_dir=chim_path, action=None)

                # insert the command to keep chimeras into the uchime command
                append_subprocess(
                    cli_command_list=vsearch_denovo_cmd,
                    options_to_add=['--chimeras', chim_out],
                    position=-2,
                    return_copy=False,
                )

            else:
                pass

            # execute the vsearch UCHIME command
            run_subprocess(vsearch_denovo_cmd, dest_dir=chim_parent, run_name=run_name, program='uchime3-denovo',
                           auto_respond=settings['automate']['auto_respond'])

        # run UCHIME de novo chimera detection on sequences that have not been denoised
        for undenoised_file in input_sorted[flip_prefix(DENOISE_PREFIX)]:

            # create an output file path for the non-chimeric read output
            nochim_out = add_prefix(file_path=undenoised_file, prefix=NOCHIM_PREFIX,
                                    dest_dir=nochim_path, action=None)

            # assemble the vsearch UCHIME3 command for the command line
            vsearch_denovo_cmd = ['vsearch', '--uchime_denovo', undenoised_file,
                                  '--abskew', str(alpha),
                                  '--nonchimeras', nochim_out,
                                  '--uchimeout', uchime_log]

            if keep_chimeras:

                # if configured to keep chimeras, create an output file path for chimeric reads
                chim_out = add_prefix(file_path=undenoised_file, prefix=flip_prefix(NOCHIM_PREFIX),
                                      dest_dir=chim_path, action=None)

                # insert the command to keep chimeras into the uchime command
                append_subprocess(
                    cli_command_list=vsearch_denovo_cmd,
                    options_to_add=['--chimeras', chim_out],
                    position=-2,
                    return_copy=False,
                )

            else:
                pass

            # execute the vsearch de novo UCHIME command
            run_subprocess(vsearch_denovo_cmd, dest_dir=chim_parent, run_name=run_name, program='uchime-denovo',
                           auto_respond=settings['automate']['auto_respond'])

    ## REFERENCE-BASED CHIMERA DETECTION #################

    # if the method is set to reference-based, then do reference-based chimera detection
    elif method == 'reference':

        # get the path to the directory with the chimera reference datasets
        chim_ref_dir = file_finder(
            reference_dir=reference_dir,
            search_glob='reference-sequences/chimera-check',
        )

        # if the directory with the chimera reference datasets contains a single directory, then replace
        #  the file path with the path to this child directory

        # create a list of contents of chim_ref_dir, ignoring hidden files like .DS_Store
        chim_ref_dir_contents = [child for child in chim_ref_dir.iterdir() if not child.name.startswith('.')]

        # if only a single item and this item is a directory, replace path variable with this child directory path
        if (len(chim_ref_dir_contents) == 1) and (chim_ref_dir_contents[0].is_dir()):
            chim_ref_dir = chim_ref_dir_contents[0]
        else:
            pass


        # sort the input files based on the DNA region; e.g., ITS1 and ITS2 have distinct reference datasets

        # create a dictionary where the key is the DNA region that will match the file tags and values are empty list
        ref_chim_regions = [region_tag.lower() for region_tag in POST_ITSX_SUFFIXES.values() if
                            not (region_tag in ['5_8S', 'LSU'])]
        input_files_sorted = {region: [] for region in ref_chim_regions}

        # create a regex that will search for any of the regions
        ref_chim_regions_re = '|'.join(ref_chim_regions)

        # if input files are directories (pacbio), create a list of all sequence files within all input directories
        updated_file_list = []
        for file_in in input_files:
            if file_in.is_dir():
                for fasta_file in file_in.glob(SEQ_FILE_GLOB):
                    updated_file_list.append(fasta_file)
            else:
                updated_file_list.append(file_in)

        # go through the input files and sort files by region
        for file_in in updated_file_list:
            wanted_region_found = re.search(ref_chim_regions_re, file_in.name, re.I)
            if wanted_region_found:
                input_files_sorted[wanted_region_found.group(0).lower()].append(file_in)
            else:
                continue

        # keep track of whether multiple regions are included in the input files; later use this to decide whether to
        #   create multiple output subdirectories, one for each region (only do this if multiple regions processed)
        input_region_count = 0
        for dna_region, region_file_list in input_files_sorted.items():
            if len(region_file_list) > 0:
                input_region_count += 1
            else:
                pass

        # create a dictionary with the same keys as the input_files_sorted list, to add paths to the correct
        #   reference files to use for each region
        chim_ref_by_region = {r:'' for r in input_files_sorted.keys()}

        # go through each DNA region
        for dna_region in chim_ref_by_region:

            # create a list of all file paths that match this DNA region in the chimera reference dir path
            region_ref_files = []

            # locate all matching directories or files for this DNA region
            for ref_file in chim_ref_dir.iterdir():

                # if a file or directory matches this region...
                match_found = re.search(dna_region, ref_file.name, re.I)
                if match_found:

                    # if the matching path is a directory, look inside directory for a sequence file
                    if ref_file.is_dir():
                        match_found_inside = [f for f in ref_file.glob(SEQ_FILE_GLOB) if re.search(dna_region, f.name, re.I)]

                        # if a single sequence file is located that matches the region, add this to the list of ref files
                        if len(match_found_inside) == 1:
                            region_ref_files.append(match_found_inside[0])

                        # if multiple sequence files match the region inside this directory, add all of them
                        elif len(match_found_inside) > 1:
                            region_ref_files.append(*match_found_inside)

                        # if no matching files are found in this directory, pass over it
                        else:
                            pass

                    # if the matching path is a sequence file, add it to the matching file list
                    elif re.search(SEQ_FILE_RE, ref_file.suffix, re.I):
                        region_ref_files.append(ref_file)

                    # if the matching path isn't a directory nor a sequence file, skip over it (don't add to list)
                    else:
                        pass

                else:
                    continue

            # if a single file / directory is located for this region...
            if len(region_ref_files) == 1:
                # add this as the path location of the chimera ref file for this region
                chim_ref_by_region.update({dna_region: region_ref_files[0]})

            # if multiple are located for this region...
            elif len(region_ref_files) > 1:
                # print an error; can't proceed with multiple matches
                err_msg = (f'Multiple chimera reference files were detected for the DNA region {dna_region} based '
                           f'on matching the region string to a file name in the directory: \n'
                           f'   {chim_ref_dir}')
                return exit_process(err_msg)

            # if no reference datasets exactly match this region...
            else:

                # likely indicates that the general-use reference file should be used for this region
                if dna_region.lower() in ['full-its', 'its-lsu']:
                    general_ref = [fasta_file for fasta_file in chim_ref_dir.glob(SEQ_FILE_GLOB)][0]
                    chim_ref_by_region.update({dna_region: general_ref})
                else:
                    err_msg = (f'A reference chimera dataset for the {dna_region} DNA region should be available '
                               f'for vsearch to use, but one was not detected.')
                    return exit_process(err_msg)


        # go through the list of input files by region...
        for dna_region, region_file_list in input_files_sorted.items():

            # if there aren't any files for this region, skip over it
            if len(region_file_list) == 0:
                continue

            else:

                # get the reference dataset to use based on the DNA region of the input files
                chim_ref_file = chim_ref_by_region[dna_region]

                # create a subdirectory for this DNA region, only if multiple regions are represented in input files
                if input_region_count > 1:
                    region_nonchim_out = mkdir_exist_ok(
                        new_dir=dna_region,
                        parent_dir=nochim_path,
                    )
                # otherwise, put directory into the main output directory
                else:
                    region_nonchim_out = nochim_path

                # process one input file at a time from this region file list
                for input_file in region_file_list:

                    # file name of the non-chimeric sequences for this sample
                    nochim_out = add_prefix(file_path=input_file, prefix=NOCHIM_PREFIX,
                                            dest_dir=region_nonchim_out, action=None)

                    vsearch_ref_cmd = ['vsearch', '--uchime_ref', input_file,
                                       '--nonchimeras', nochim_out,
                                       '--uchimeout', uchime_log,
                                       '--db', chim_ref_file]

                    if keep_chimeras:

                        # create a subdirectory for this DNA region, only if multiple regions are represented in input files
                        if input_region_count > 1:
                            region_chim_out = mkdir_exist_ok(
                                new_dir=dna_region,
                                parent_dir=chim_path,
                            )
                        # otherwise, put directory into the main output directory
                        else:
                            region_chim_out = chim_path

                        # file name of the chimeric sequences for this sample (if keep_chimeras=True)
                        chim_out = add_prefix(file_path=input_file, prefix=flip_prefix(NOCHIM_PREFIX),
                                              dest_dir=region_chim_out, action=None)

                        # insert the command to keep chimeras into the uchime command
                        append_subprocess(
                            cli_command_list=vsearch_ref_cmd,
                            options_to_add=['--chimeras', chim_out],
                            position=5,
                            return_copy=False,
                        )

                    else:
                        pass

                    # execute the chimera detection vsearch command for this sample sequence file
                    run_subprocess(vsearch_ref_cmd, dest_dir=chim_parent, run_name=run_name, program='uchime-ref',
                                   auto_respond=settings['automate']['auto_respond'])

    else:
        pass

    ## OUTPUT SUMMARY FILE WITH TABLE OF SAMPLES WITHOUT NON-CHIMERA READS

    # create an empty dictionary to add sample IDs and read counts of empty files only (no sequences)
    empty_nonchim = {}

    # go through each non-chimeric file that was just created
    for nonchim_file in nochim_path.glob(f'{NOCHIM_PREFIX}*fasta'):

        # for each non-chimeric file, count the number of sequences (read count)
        read_count = 0
        with open(nonchim_file) as fasta_in:
            for record in SeqIO.parse(fasta_in, 'fasta'):
                read_count += 1

        # if there are no sequences in the non-chim file for this sample...
        if read_count == 0:

            # get the sample ID
            sample_id = get_sample_id(file_path=nonchim_file)

            # append the sample ID and read count to the empty non-chim dictionary
            empty_nonchim.update({sample_id: read_count})

    # after going through each non-chimeric file, write out the empty non-chimeric samples to a summary file
    empty_nonchim_out = chim_parent / f'{run_name}_no-nonchim-reads.txt'
    with open(empty_nonchim_out, 'wt') as fout:
        fout.write(f'read count\tsample id\n')
        for sample_id, read_count in empty_nonchim.items():
            fout.write(f'{read_count}\t{sample_id}\n')

    return None


# def concat_regions(dir_path, file_map, regions_to_concat=['ITS1', '5_8S', 'ITS2'], verbose=True, header_delim=None):
#     '''
#     Concatenates provided regions for each read.
#
#     Looks in the provided directory for files containing the substring of the
#     regions provided in the regions_to_concat list. For each region in the list,
#     it will store read ID, region, read length, and full-length read count in a
#     dictionary, where each read ID is the primary key with details stored for each
#     region (read length, full-length read count). Then will combine this information
#     from the regions for each read ID. The sequences from each region are concatenated
#     in the order they appear in the provided regions_to_concat list. The read length of
#     the final concatenated reads is the sum of the read lengths of the regions. It will
#     check that the full-length read count is the same for each of the regions, which
#     further confirms that the regions are from the same full-length read. Then will write
#     these concatenated reads and updated headers to a new output fasta file.
#     :param dir_path: directory containing the fasta files to concatenate.
#     :param regions_to_concat: an ORDERED list of regions to concatenated, with the 5' region
#     first and the 3' region last; defaults to ITS regions
#     :param header_delim: delimiter to use for separating the components of the header in
#     the concatenated fasta output; defaults to semicolon
#     :return: None, but will write out a fasta file
#     '''
#
#     ## IMPORT CONFIGURATION SETTINGS
#
#     # import configuration settings as dictionary
#     settings = get_settings(file_map)
#
#     # get bioinformatics run name
#     run_name = settings['run_details']['run_name']
#
#     # check whether to use header_delim from settings or directly from function parameter
#
#     # if header_delim set to default...
#     if header_delim is None:
#
#         # then import from the configuration file
#         header_delim = settings['formatting']['header_delim']
#
#     # if header_delim is not the default, use whatever value is assigned to it
#     else:
#         pass
#
#     # created a variable to switched to True if a region cannot be located without an ITSx issue
#     #  moved this further down, I'm not sure if its alright down there but we'll see
#     # error_occurred = False
#
#     # create an output file using the run name for this bioinformatics run; do so here so you can remove
#     err_log_path = (dir_path.parent / f'{run_name}').with_suffix('.no_concat.err')
#
#
#     # determine whether the input directory contains sample subdirectories or if its a sample directory
#
#     # create list of subdirectories within input directory
#     child_dirs = [child.is_dir() for child in dir_path.glob(f'{ITSX_PREFIX}*')]
#
#     # if the input directory contains only subdirectories, then create list of subdirectories to loop through
#     if all(child_dirs):
#         input_dirs = list(dir_path.glob('*'))
#
#     # if none of the contents of input directory are directories, then its a sample dir and needs to be cast as
#     #  list so that you can still use the loop method to work through
#     elif ~all(child_dirs):
#         input_dirs = [dir_path]
#
#     # if there's a mix of directories and not directories within the input path, then assume that the input provided
#     #  is a sample directory, and look for ITSx output within the input directory but warn that errors might follow
#     else:
#         print(f'WARNING. Encountered a mix of directories and non-directories in the input for the '
#               f'concat_regions() function from climush.bioinfo. Proceeding by searching the input '
#               f'directory for the fasta files to concatenate but may encounter an error because the '
#               f'function does not know whether subdirectories contain per-sample fasta files or if '
#               f'the subdirectories are unrelated.\n')
#         input_dirs = [dir_path]
#
#     # create dictionary for acceptable input regions
#     accepted_regions = {r'(^S.+?(?=\bs))|(?<!-)SSU': 'SSU',  # could be better but not going to obsess
#                         r'ITS.?1': 'ITS1',  # checks for possible punctuation within
#                         r'5.?8S': '5.8S',  # checks for possible punctuation within
#                         r'ITS.?2': 'ITS2',  # checks for possible punctuation within
#                         r'(^L.+?(?=\bs))|(?<!-)LSU': 'LSU',  # similar to SSU
#                         r'(full)?.?ITS(?!.?\d)': 'ITS'}  # may include full, but can't have number after (w or w/o punct)
#
#     # check input list for valid entries, and then format string to consistent format
#     # returns dict now with formatted string as key, regex as value
#     def create_region_dict(region_list):
#         output_dict = {}
#         for region_str in region_list:
#             for r in accepted_regions.keys():  # compare input region to each regex in dict
#                 if re.search(r, region_str, re.I):  # if it does match, stop loop and return the formatted version
#                     output_dict.update({accepted_regions[r]: r})
#                     break
#                 else:
#                     continue
#             else:
#                 return print(f'FAILURE. Did not recognize the input region {region_str}. Please choose from the '
#                              f'recognized regions: {list(accepted_regions.values())}\n')
#         return output_dict
#
#     # make nested dict: read ID > region > [sequence, read len, derep01_reads(?)]
#     # ADD CHECK FOR MULTIPLE FILES MATCHING THE REGION IN THIS DIRECTORY
#     def search_path_with_regex(dir_path, regex, return_all=False):
#         found_file_paths = []
#         for file in dir_path.glob('*'):
#             if re.search(regex, file.name, re.I):
#                 found_file_paths.append(file)
#             else:
#                 continue
#
#         if len(found_file_paths) > 1:
#             if return_all:
#                 return found_file_paths
#             else:
#                 print(f'WARNING. Multiple files matching the provided regex {regex} were located in '
#                       f'{dir_path.name}. Returning only the first match.\n')  # REPLACE WITH OPTIONS?
#                 return found_file_paths[0]
#         elif len(found_file_paths) == 0:
#             # print(f'FAILURE. No files matching the provided regex {regex} were located in '
#             #       f'{dir_path.name}.\n')
#             return None
#         else:
#             return found_file_paths[0]
#
#     # confirm all regions data was acquired, then combine (used below in loop)
#     def join_if_all_present(pulled_data_dict):
#
#         # create an empty dictionary to add read information to once it is checked
#         joined_data = {'seq':'',
#                        'id':[]}
#
#         # for each value in details_to_include from the configuration file settings...
#         for k in pulled_data_dict.keys():
#
#             # check if the number of values for this info group is equal to the number of regions in regions_to_concat
#             if len(pulled_data_dict[k]) == len(regions_to_concat):
#
#                 # if this is the key for the Seq object (i.e, a list of sequences from each region to concat)...
#                 if any(isinstance(i, Seq) for i in pulled_data_dict[k]):
#
#                     # then join all the sequences in this list (i.e., actually concatenate the regions to one read)
#                     joined_data['seq'] = Seq('').join(pulled_data_dict[k])  # must join with empty seq
#
#                 # if this is not the key for a Seq object (i.e., not a sequence)...
#                 else:
#
#                     # format the combined information for each category based on which category it is;
#                     #  this format will directly affect how the information will appear in the fasta file headers
#
#                     # if it is the read length, then sum all of the values and confirm the concat length is correct
#                     if re.search(r'len$', k, re.I):
#
#                         # compare expected and observed read length of the concatenated read before adding info to header
#                         exp_rl = sum([int(x) for x in pulled_data_dict[k]])  # calc expected length based on sum of each region length
#                         obs_rl = len(joined_data['seq'])  # calc observed length based on the length of the concatenated read
#
#                         # if the expected and observed read lengths for the concatenated read match...
#                         if exp_rl == obs_rl:
#
#                             # then create a header component for read length
#                             joined_data['id'].append(f'{concat_name} {obs_rl}bp')
#
#                         else:
#                             # if these values don't match, then print an error and return None (i.e., don't concat)
#                             msg = f'ERROR. The sum of the values in the headers for the post-ITSx region files does not '\
#                                   f'match the length of the concatenated read. Please check the ITSx output for the '\
#                                   f'following regions:\n'\
#                                   f'   read ID = {read}\n' \
#                                   f'   concat regions = {regions_to_concat}\n' \
#                                   f'   expected length (nt) = {exp_rl}\n' \
#                                   f'   concat length (nt) =   {obs_rl}\n'
#                             return print(msg)
#
#                     # if this key is the read count for full-length reads (number of identical full-length reads/sample)
#                     elif re.search(r'^full\-len|count', k, re.I):
#
#                         # if all read counts are identical among the regions (they should be)...
#                         if len(set(pulled_data_dict[k])) == 1:
#
#                             # then create the header component for read count
#                             joined_data['id'].append(f'{k}={int(pulled_data_dict[k][0])}')  # if all same, doesn't matter which you chose
#
#                         # there's likely some larger error at play if a given read doesn't match in read count among regions
#                         else:
#                             msg = f'ERROR. The read counts for the regions of this read are different lengths:\n'\
#                                   f'  {pulled_data_dict[k]}\n'
#                             return print(msg)
#
#                     # if this key is not recognized as a default header key, then make generic version to add
#                     else:
#                         joined_data['id'].append(f'{k} = {pulled_data_dict[k]}')
#
#             # if there is missing information for at least one of the regions to concatenate...
#             else:
#                 # print an error message, and return Non
#                 # technically, this function should never end up running if there are missing regions for a read
#                 #   so if this error message prints, there is likely an issue elsewhere
#                 msg = f'Did not successfully pull all {k} from the regions to concatenate ' \
#                       f'({len(pulled_data_dict[k])} out of {len(regions_to_concat)} were pulled). ' \
#                       f'See the .no_concat.txt file for more details. \n'
#                 return print(msg)
#
#         return joined_data
#
#     # if a region is missing for a read, determine whether ITSx couldn't find or if another unknown issue
#     def explain_missing_region(read_id, err_region, sample_dir):
#         '''
#         If a read is missing from its sample's subregion fasta file, determine whether the missing subregion
#         was not detected by ITSx (an ITSx error) or if an error independent of ITSx occurred.
#         :param read_id: the read ID of the sample that was missing a subregion, which is the first element
#         in the fasta file header
#         :param err_region: the region of the sequence that could not be located (e.g., ITS1, ITS2)
#         :param sample_dir: directory containing the output from ITSx for this sample
#         :return: if an ITSx error, returns 'itsx' and logs error to common .no_concat.txt file for
#         the bioinformatics run that this sample is processed with and does not concatenate any
#         regions for this read; if not an ITSx error, returns 'other' and prints error message; sys.exit() has
#         to be carried out independent of this function in order to interrupt concat_regions()
#         '''
#
#         # determine whether an ITSx error by checking the <sample-id>.problematic.txt output file from ITSx
#
#         # get the path to the <sample-id>.problematic.txt output from ITSx
#         itsx_detect_error_path = next(sample_dir.glob('*problematic*'))
#
#         # create generator to find line that contains this read in the problematic.txt file, yield its error message
#         def find_problematic_read(err_file, read_id):
#             for err in err_file.readlines():
#                 if err.startswith(read_id):
#                     yield err.split('\t')[-1].strip()
#                 else:
#                     continue
#
#         # get the error message created by ITSx for this read
#         with open(itsx_detect_error_path, 'rt') as err_records:
#             errors = next(find_problematic_read(err_file=err_records, read_id=read_id))
#
#         # if the missing region is located in the ITSx error, then ITSx couldn't detect the region
#
#         # get the regex used to search files, etc. produced by ITSx for this region
#         region_regex = region_search_dict[err_region]
#
#         # search for the region causing the error in the error message for this read produced by ITSx
#         if re.search(region_regex, errors, re.I):
#             return 'itsx'
#         # if the region is not detected in the error message from ITSx for this read, then its some other unknown issue
#         else:
#             return 'other'
#
#     # create a dictionary of regions to concat, where keys are formatted regions and
#     #   values are regex to use for file search
#     region_search_dict = create_region_dict(regions_to_concat)
#
#     ## PER SAMPLE ########################################################################################
#
#     # go through each sample's directory...
#     for sample_dir in input_dirs:
#
#         print(f'sample dir = {sample_dir}\n')
#
#         # log progress by printing the name of the sample that is being processed; will append SUCCESS to string
#         #   if process completes successfully
#         print(f'\n{sample_dir.name} running...', end='\r')
#
#         # create a list of the fasta files to concat based on input regions
#         fastas_to_concat = [search_path_with_regex(sample_dir, regex=r) for r in region_search_dict.values()]
#
#         # get the sample ID from any of these fasta files
#         sample_name = get_sample_id(file_path = fastas_to_concat[0])
#
#         # get the file prefix used for the ITSx output files, will be used as prefix for the concatenated file output
#         file_prefix = re.search(PREFIX_RE, fastas_to_concat[0].name, re.I).group(0)
#
#         # create a list of errors when encountered in reads; formatted to write this list to an error log file
#         #  as soon as this sample is completed
#         # PATCH - I keep getting the same read printing over and over, so I think once the value that checks whether
#         #  an error has occurred sets to True, it isn't reseting such that every read after the error read triggers
#         #  the write of the same read information to the file
#         reads_with_missing_regions = set()
#
#         # collect read IDs of reads in this sample that itsx flagged as chimeric; will skip over completely later on
#         chimera_fasta = search_path_with_regex(sample_dir, regex=r'\.chimeric\.')
#
#         # if no itsx output for chimeric reads are located for this sample, then create an empty list of chimeras
#         if chimera_fasta is None:
#             chimer_read_ids = []
#         # if the itsx chimeric output fasta is located, extract the read IDs from this file
#         else:
#             with open(chimera_fasta, 'rt') as chimer_in:
#                 chimer_headers = [line for line in chimer_in.readlines() if line.startswith('>')]
#                 chimer_read_ids = [read.split(';')[0].replace('>', '') for read in chimer_headers]
#
#         # summarize information on errors and missing regions for each sample
#         reads_with_nonitsx_errors = []
#         reads_without_issue = []
#         reads_with_itsx_errors = []
#         error_occurred = False
#         regions_missing_in_sample = False
#         full_read_count = 0
#
#         # auto-assign name of joined regions for common concats, otherwise prompt user for input
#         formatted_input_regions = list(region_search_dict.keys())
#         if formatted_input_regions == ['ITS1', '5.8S', 'ITS2']:
#             concat_name = 'full-ITS'
#         elif formatted_input_regions == ['ITS1', '5.8S', 'ITS2', 'LSU']:
#             concat_name = 'ITS-LSU'
#         else:
#             msg = f'Please provide a name to assign to this combination of subregions (use '\
#                   f'hyphens instead of spaces):\n'
#             default_name = '-'.join(formatted_input_regions)
#             concat_name = prompt_generic(message=msg, auto_respond=settings['automate']['auto_respond'],
#                                          auto_response=default_name)
#
#         # create a list of components to include in the new fasta files within read header
#         #  added these to the default settings, but not yet in the user settings
#         details_to_include = settings['formatting']['fasta_headers']
#
#         ## PER REGION ########################################################################################
#
#         # for this sample, go through each of the itsx output files for the regions to concat, collecting info
#         concat_dict = {}
#         for fasta in fastas_to_concat:
#
#             ## PER READ #######################################################################################
#
#             for record in SeqIO.parse(fasta, 'fasta'):
#
#                 read_id = re.search(READ_ID_RENAMED_RE, record.description).group(0)
#
#                 # if this read ID matches one in ITSx's chimeric read output file, then skip to the next read
#                 if read_id in chimer_read_ids:
#                     continue
#
#                 # get information from the post-itsx read header, which depends to an extent on a certain header format
#                 try:
#                     region = re.search(READ_REGION_RENAMED_RE, record.description, re.I).group(0)
#                     length = re.search(READ_LEN_RENAMED_RE, record.description, re.I).group(0)
#                     derep01_reads = re.search(READ_COUNT_RENAMED_RE, record.description, re.I).group(0)
#                     details = {k: v for k, v in zip(details_to_include, [record.seq, length, derep01_reads])}
#
#                     if read_id in concat_dict:  # if this sample is already in the dictionary...
#                         concat_dict[read_id].update({region: details})  # append the new region info to it
#                     else:  # if it isn't yet in the dictionary...
#                         concat_dict[read_id] = {region: details}  # add a key value pair for this sample
#
#                 # if the itsx output isn't formatted as expected, can't get important info, so print error and exit
#                 except AttributeError:
#                     print(f'{sample_dir.name} running... ERROR!')
#                     print(f'  There was an issue extracting information from the fasta header of the input file:\n '
#                           f'     {fasta} \n'
#                           f'  This issue occurred for read:\n'
#                           f'     {read_id} \n'
#                           f'  This will occur if the sequence files produced by ITSx have not yet \n'
#                           f'  been reformatted. Please run the rename_read_header function and retry.\n')
#                     return sys.exit(1)
#
#         ## PER READ ############################################################################################
#
#         # create BioSeq records for each read for each region to concatenate
#         concatenated_records = []
#
#         for read in concat_dict.keys():
#
#             # keep track of the number of reads in this sample
#             full_read_count += 1
#
#             regions_missing_in_read = False
#
#             # create key-value pairs where key is the info to include in BioSeq record, and empty value list to add to
#             record_details = {k:[] for k in details_to_include}
#
#             # missing region
#             missing_regions = []
#
#             ## PER REGION #######################################################################################
#
#             # pull info from each region for this read
#             for region in region_search_dict.keys():
#
#                 # first try to index the given region for this read...
#                 try:
#                     # get the sequence and read description details for this region
#                     region_dict = concat_dict[read][region]
#
#                     # add the read description to the record_details dictionary for this read
#                     for det in details_to_include:
#                         if det in record_details:
#                             record_details[det].append(region_dict[det])  # append before sum to ensure all regions there
#                         else:
#                             record_details[det] = region_dict[det]  # append before sum to ensure all regions there
#
#                 # if the region doesn't exist for this read, determine whether an error from ITSx or unknown source
#                 except KeyError:
#
#                     # determine the cause of the missing region for this read; will log accordingly within function
#                     cause = explain_missing_region(read_id=read, err_region=region, sample_dir=sample_dir)
#
#                     # if the region is missing because ITSx could not detect it within the full-length read...
#                     if cause == 'itsx':
#
#                         # change this variable to True, used for formatting the error log file
#                         regions_missing_in_read = True
#
#                         # add this region to a list of missing regions for this read; these regions will be written to
#                         #  a log file once all regions to concatenate for this read have been searched
#                         missing_regions.append(region)
#
#                     # if ITSx is not the reason that this region is missing, and the reason is therefore unknown..
#                     else:
#
#                         # if ITSx can't find the end of the 5.8S region, it will not produce the ITS2 either;
#                         #   however, it will not log this error in the <sample-id>.problematic.txt log file
#                         #   so, if an itsx issue with 5.8S has been found for this read, chances are ITSx also could
#                         #   not find the ITS2 and it isn't some extraneous error in locating the ITS2 seq for this read
#                         if regions_missing_in_read:
#                             missing_regions.append(region)
#                         else:
#                             error_occurred = True
#                             print(f'{sample_dir.name} running... ERROR!')
#                             error_msg = (f'ERROR. The {region} region could not be located and ITSx did not log an error '
#                                          f'with this read, so some unknown issue is preventing the detection of the '
#                                          f'{region} for this read. The regions were not concatenated. Please check the '
#                                          f'ITSx output file for this sample to troubleshoot.\n'
#                                          f'   sample: {sample_name}\n'
#                                          f'   read:   {read}\n')
#                             print(f'   {error_msg}')
#                             # break out of this loop, i.e., stop going through regions for this read and move to next read
#                             #  breaking here will prevent any record from being made for this read, concat won't occur
#                             # break
#
#             ## PER READ #######################################################################################
#
#             # if there were any regions not detected by ITSx, then write this information out to a log file
#             if len(missing_regions) > 0:
#
#                 reads_with_itsx_errors.append(read)
#
#                 # format the sample ID, read ID, and non-detected regions to add to list that will append to log file
#                 reads_with_missing_regions.add(f'   {read} = {missing_regions}\n')
#
#             # if a non-ITSx error occurred...
#             elif error_occurred:
#                 reads_with_nonitsx_errors.append(read)
#
#             # if there weren't any issues detecting the regions required to concatenate...
#             else:
#
#                 combined_data = join_if_all_present(record_details)
#
#                 if (len(combined_data['id']) == 0) or (combined_data is None):
#                     reads_with_nonitsx_errors.append(read)
#                 else:
#                     reads_without_issue.append(read)
#                     concat_record = SeqRecord(combined_data['seq'],
#                                               name=read,
#                                               id=header_delim.join([read] + combined_data['id']),
#                                               description='')
#                     concatenated_records.append(concat_record)
#
#         error_reads = reads_with_nonitsx_errors + reads_with_itsx_errors
#         total_errs = len(error_reads)
#
#         if verbose:
#             print(f'Total reads with errors: {total_errs} ({(total_errs / full_read_count)*100:.2f}%)\n'
#                   f'   ITSx errors:     {len(reads_with_itsx_errors)}\n'
#                   f'   non-ITSx errors: {len(reads_with_nonitsx_errors)}\n'
#                   f'Total reads without errors: {len(reads_without_issue)}\n')
#
#         ## PER SAMPLE #######################################################################################
#
#         # write out the concatenated read records to a new concatenated reads fasta file
#         fasta_out = sample_dir / f'{file_prefix}_{sample_name}.{concat_name}.fasta'
#
#         SeqIO.write(concatenated_records, fasta_out, 'fasta')
#
#         # write out the information of the reads that had regions missing and weren't concatenated
#         if len(reads_with_missing_regions) > 0:
#             with open(err_log_path, 'at') as fout:
#
#                 # create a sample header to better distinguish between samples in the error log
#                 fout.write(f'\n{sample_name}  -----------------------------\n')
#
#                 # write the name of the read and the list of missing sequence regions for this read
#                 for read_error in reads_with_missing_regions:
#                     fout.write(read_error)
#
#             print(f'{sample_dir.name} running... WARNING!')
#             print(f'  discarded reads in which ITSx could not detect subregions.')
#         else:
#             if fasta_out.is_file():
#                 print(f'{sample_dir.name} running... success!')
#             else:
#                 print(f'{sample_dir.name} running... ERROR!\n')
#                 print(f'  concatenated output fasta file not detected.')
#
#     return None

# def check_concat_output(itsx_dir, full_len_dir, num_bp_compare, file_map, write_to_log=True, same_threshold=99, header_delim=None):
#
#     # import configuration settings
#     settings = get_settings(file_map)
#     run_name = settings['run_details']['run_name']
#
#
#     # check which value to use for the read header delimiter
#     if header_delim is None:
#         header_delim = settings['formatting']['header_delim']
#     else:
#         pass
#
#     bp = int(num_bp_compare)
#
#     # import records as dictionary, where the read ID is the key
#     def name_as_key(record):
#         header = record.id
#         sample_read_id = header.split(header_delim)[0]
#         return sample_read_id
#
#     # import full ITS records and LSU (for revcomp comparison)
#     full_its_records = SeqIO.to_dict(SeqIO.parse(next(itsx_dir.glob('*full_ITS*')), 'fasta'), key_function=name_as_key)
#     lsu_records = SeqIO.to_dict(SeqIO.parse(next(itsx_dir.glob('*LSU*')), 'fasta'), key_function=name_as_key)
#
#     # collect the full length _ bp sequences at start, full ITS _ bp sequences at start, and LSU _ bp seqs at end
#     post_itsx_diffs = {'read_id':[],
#                        'full_start':[],
#                        'its_start':[],
#                        'lsu_end':[]}
#     for record in SeqIO.parse(full_len_dir, 'fasta'):
#
#         sample_id = get_sample_id(file_path = full_len_dir.name)
#         read_num = re.search(READ_ID_OG_RE, record.id).group(0)
#         read_id = '_'.join([sample_id, read_num])
#
#         renamed_record = SeqRecord(record.seq, id=read_id, description='')
#
#         full_len_seq = renamed_record.seq
#         full_its_seq = full_its_records[renamed_record.id].seq
#         try:
#             lsu_seq = lsu_records[renamed_record.id].seq
#         except:
#             lsu_seq = 'N'*bp
#
#
#         post_itsx_diffs['read_id'].append(renamed_record.id)
#         post_itsx_diffs['full_start'].append(str(full_len_seq[:bp]))
#         post_itsx_diffs['its_start'].append(str(full_its_seq[:bp]))
#         post_itsx_diffs['lsu_end'].append(str(lsu_seq[-bp:]))
#
#     # create df from this dictionary
#     diff_df = pd.DataFrame.from_dict(post_itsx_diffs)
#
#     # compare the read slices first as full-ITS fwd, then full-LSU revcomp
#     comparison_result = []  # True if same
#     equal_as = []  # str describing nature of match: fwd, revcomp, no LSU, unresolved (error)
#     for i in range(diff_df.shape[0]):
#         full_fwd = Seq(diff_df['full_start'][i])
#         its_fwd = Seq(diff_df['its_start'][i])
#         lsu_end = Seq(diff_df['lsu_end'][i])
#         if full_fwd == its_fwd:
#             comparison_result.append(True)
#             equal_as.append('forward_strand')
#             continue
#         else:
#             if 'N' in lsu_end:
#                 comparison_result.append(False)
#                 equal_as.append('no_lsu')
#             else:
#                 its_revcomp = lsu_end.reverse_complement()
#                 if full_fwd == its_revcomp:
#                     comparison_result.append(True)
#                     equal_as.append('reverse_comp')
#                     continue
#                 else:
#                     comparison_result.append(False)
#                     equal_as.append('unresolved')
#
#     diff_df['comparison'] = comparison_result
#     diff_df['equal_as'] = equal_as
#
#     percent_same = (sum(diff_df['comparison']) / diff_df.shape[0])*100
#     revcomp_detect = ((diff_df[diff_df['equal_as'] == 'reverse_comp'].shape[0]) /sum(diff_df['comparison']))*100
#     forward_detect = ((diff_df[diff_df['equal_as'] == 'forward_strand'].shape[0]) /sum(diff_df['comparison']))*100
#     missing_lsu = ((diff_df[diff_df['equal_as'] == 'no_lsu'].shape[0]) /sum(diff_df['comparison']))*100
#
#     if write_to_log:
#         log_dir = itsx_dir / 'itsx_inspection_logs'
#         log_dir.mkdir(exist_ok=True)
#         diff_df.to_csv(log_dir / f'itsx_inspection_comparison.csv', index=False)
#         with open((log_dir / 'itsx_inspection_summary.log'), 'wt') as fout:
#             fout.write(f'sample: {get_sample_id(file_path = full_len_dir.name)}\n')
#             fout.write(f'date: {datetime.today().strftime("%Y-%M-%d")}\n\n')
#             fout.write(f'All full-length reads were compared to their ITS region after running ITSx. The first'
#                        f' {bp} bp of the full-length reads were aligned to the first {bp} bp of their corresponding '
#                        f'full ITS read that was concatenated after running ITSx. If the first {bp} bp did not '
#                        f'match, then the reverse complement of the last {bp} bp of their LSU region was compared '
#                        f'to the first {bp} bp of the full-length read. This is because reads are not reoriented '
#                        f'until done so by ITSx. \n\n')
#             fout.write(f'For {percent_same:.2f}% of these comparisons (out of {diff_df.shape[0]} total comparisons), '
#                        f'the first {bp} bp of the full length read matched the first {bp} bp of either the forward '
#                        f'ITS or reverse complement of the end of the LSU for that read.\n\n')
#             fout.write(f'Out of the {percent_same:.2f}% that matched, {forward_detect:.2f}% were matched by comparing '
#                        f'the full-length read to the start of the ITS, while {revcomp_detect:.2f}% were matched by '
#                        f'comparing the full-length read to the reverse complement of the end of the LSU region.\n\n')
#             if missing_lsu > 0:
#                 fout.write(f'{int(missing_lsu*diff_df.shape[0])} reads ({missing_lsu:.2f}% of total) did not have '
#                            f'an LSU region to compare to, since ITSx could not locate the LSU for the read(s). This '
#                            f'accounts for the proportion of reads that could not be matched.\n')
#     else:
#         if percent_same < same_threshold:
#             print(f'WARNING. The percent of reads that matched to their pre-ITSx full-length read was above '
#                   f'{same_threshold}%, at {percent_same}.\n')
#
#     return None

def create_query_fasta(input_dir, output_path, reference_dir, file_fmt='fasta', keep_unsorted=False):

    # process input_dir to allow for either a directory or list of files; either way, must be Path object(s)

    # if the input directory is a list...
    if isinstance(input_dir, list):
        # go through each item in the list and confirm that they are path objects
        for i in input_dir:
            if is_pathclass(i, exit_if_false=False):
                pass
    elif is_pathclass(input_dir, exit_if_false=False):
        if input_dir.is_dir():
            pass
    else:
        err_msg = (f'The input for the function {create_query_fasta.__name__} can be either a PosixPath to a directory, '
                   f'or a list of PosixPaths. However, the input provided:\n'
                   f'   {input_dir}\n'
                   f'is neither of these. Please provide either a PosixPath to a directory, or a list of PosixPath file'
                   f'paths.\n')
        exit_process(messages=err_msg)

    # import settings from the configuration .toml file
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    # rename the read headers, post-chimera detection, so that the read ID contains the sample ID with sample= prefix
    rename_read_header(input_dir=input_dir, run_name=run_name, append_sample_str=True)

    # create an empty list to add all records (reads) from all input fasta files in input_dir
    total_record_list = []

    # create a set of input sample IDs, to be used when checking the output fasta file that is created
    samples_in_input = set()

    # count number of reads read in from input file
    input_read_count = 0

    # go through each sequence file in the input path...
    for file in input_dir.glob(SEQ_FILE_GLOB):

        # get the sample ID from the file name, and add it to the input sample set
        samples_in_input.add(get_sample_id(file))

        # parse each record (read) from each sample file...
        for record in SeqIO.parse(file, file_fmt):

            # add to input read counter
            input_read_count += 1

            # add it to the total record list
            total_record_list.append(record)

    # write out the list of all records across all samples to a query fasta file
    output_path.mkdir(exist_ok=True)  # make directory of input path
    seq_platform = get_seq_platform(file)  # get the sequence platform from the last file iterated through
    combined_unsorted_out = (output_path / f'{QUERY_PREFIX}_{seq_platform}_{run_name}_unsorted').with_suffix('.' + file_fmt)
    SeqIO.write(total_record_list, combined_unsorted_out, file_fmt)


    ## SORT OUTPUT BY DECREASING READ ABUNDANCE

    # create a file path to write the sorted sequence file to
    combined_sorted_out = (output_path / f'{QUERY_PREFIX}_{seq_platform}_{run_name}_sorted').with_suffix('.' + file_fmt)

    # assemble the vsearch command that will take the unsorted query file and arrange all seqs by decreasing read abundance
    vsearch_sort_cmd = ['vsearch', '--sortbysize', combined_unsorted_out,
                        '--output', combined_sorted_out, '--sizeout']

    # execute the sorting vsearch command
    run_subprocess(vsearch_sort_cmd, dest_dir=output_path, run_name=run_name,
                   program='vsearch-sort', separate_sample_output=True,
                   auto_respond=settings['automate']['auto_respond'])


    ## CHECK OUTPUT FILE

    # create an empty set of read IDs (which are the sample IDs, after renaming) to add IDs in output to
    samples_in_output = set()

    # check read counts in output
    output_read_count = 0

    # read in the output file to confirm that all samples that should be represented, are represented
    with open(combined_sorted_out) as fasta_in:

        # go through each record in the output file...
        for record in SeqIO.parse(fasta_in, file_fmt):

            # add to the output file read count counter
            output_read_count += 1

            # get the sample ID portion of the read header (i.e., the first item,
            #   removing the read abundance annotation)
            read_id = record.id.split(settings['formatting']['header_delim'])[0]

            # add the read_id to the set of samples confirmed to be in the output
            samples_in_output.add(read_id)

    # compare the samples that are in the input to those in the output
    dropped_samples = list(samples_in_input.difference(samples_in_output))

    # if there are samples from the input that aren't in the output file...
    if len(dropped_samples) > 0:

        # if check_chimeras() has run and produced its summary file, compare to this file
        no_nonchim_path = input_dir.parent / f'{run_name}_no-nonchim-reads.txt'
        if no_nonchim_path.is_file():
            # check if the samples are missing because the non-chimera fasta files are empty (no reads)

            # locate the samples with empty non-chimera fasta files by checking the
            #    no-nonchim summary from check_chimeras()
            empty_samples = pd.read_table(no_nonchim_path, delimiter='\t')['sample id'].to_list()

            # if the list of no-nonchim samples matches the list of dropped samples, then everything worked as expected
            if empty_samples == dropped_samples:
                pass
            else:
                dropped_samples_unaccounted = list(set(empty_samples).difference(set(dropped_samples)))
                print(f'ERROR. {len(dropped_samples_unaccounted)} samples are missing from the query fasta file output:\n'
                      f'  {combined_sorted_out.name}\n'
                      f'This fasta file should contain reads from every sample that has reads in the non-chimera output '
                      f'fasta files. After accounting for samples with no non-chimeric reads, the following samples '
                      f'were expected to be present in the query fasta file but are not:\n')
                print_indented_list(dropped_samples_unaccounted)
                exit_process(message='')  # leave out message, I print above so I can use my fnc to print the indented list

        # if check_chimeras() has not run, or is not currently the script being run, then write out missing samples
        #   to a summary file in the input directory to this function
        else:

            dropped_samples_path = input_dir.parent / f'{run_name}_missing-from-pooled-seqs.txt'

            print(f'WARNING. {len(dropped_samples)} samples are missing from the pooled sequencing .fasta file. '
                  f'See a list of these missing samples in the summary file:\n'
                  f'   {dropped_samples_path}')

            with open(dropped_samples_path, 'wt') as fout:
                fout.write(f'The following samples are missing from the pooled sequencing file:\n')
                for sample in dropped_samples:
                    fout.write(f'   {sample}\n')

    # if there are no dropped samples detected, then the everything worked as expected, and none of the samples were
    #   devoid of non-chimera reads after chimera detection
    else:

        # if you want to keep the unsorted version of the query output file, do nothing
        if keep_unsorted:
            pass
        # if you do not want to keep the unsorted version of the query output file...
        else:
            # remove the unsorted version of the query fasta file
            combined_unsorted_out.unlink()
            # rename the sorted version by dropping its prefix
            combined_sorted_out_nosuffix = (combined_sorted_out.parent / combined_sorted_out.stem.split('_')[:-1]).with_suffix(combined_sorted_out.suffix)
            combined_sorted_out.rename(combined_sorted_out_nosuffix)


    print(f'Number of input reads = {input_read_count}\n'
          f'Number of output reads = {output_read_count}\n')

    return combined_sorted_out


def denoise(input_files, reference_dir, alpha, minsize, clust_threshold, pool_samples=True):

    ## GET SETTINGS ###################################

    # read in the configuration settings, get bioinformatics run name from settings
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']


    ## CREATE OUTPUT DIRECTORIES ######################

    # make the main output directory for denoised output
    denoise_parent = mkdir_exist_ok(new_dir=file_map['pipeline-output']['denoised'])

    # within main denoise output directory, create a subdirectory for the centroid sequence outputs
    denoise_path = mkdir_exist_ok(new_dir=f'{DENOISE_PREFIX}_{run_name}', parent_dir=denoise_parent)


    ## FILTER OUT NON-ILLUMINA FILES ###################

    illumina_files = []
    non_illumina_files = []

    for file in input_files:

        # infer the sequencing platform from the file name
        file_platform = get_seq_platform(file)

        # sort file paths based on whether illumina, or generally not illumina
        if file_platform == 'illumina':
            illumina_files.append(file)
        else:
            non_illumina_files.append(file)


    ## DENOISE #########################################

    # if you want to pool all reads across samples (default)...
    if pool_samples:

        # use create_query_fasta to rename read headers, combine into one fasta, and return the path to this new file
        input_parent_dir = illumina_files[0].parent  # create_query_fasta requires a directory, not a list
        combined_fasta = create_query_fasta(input_dir=input_parent_dir, file_map=file_map, output_path=denoise_path)

        # create an output file path to write denoised centroid sequences to
        denoise_output = add_prefix(file_path=combined_fasta, prefix=DENOISE_PREFIX, dest_dir=denoise_path, action=None)
        denoise_otutab_out = (denoise_output.parent / (denoise_output.stem + '_otu-table')).with_suffix('.txt')
        denoise_clust_out = (denoise_output.parent / (denoise_output.stem + '_clusters')).with_suffix('.uc')

        # assemble the vsearch UNOISE command for the command line
        vsearch_unoise_cmd = ['vsearch', '--cluster_unoise', combined_fasta,
                              '--centroids', denoise_output,
                              '--otutabout', denoise_otutab_out,
                              '--uc', denoise_clust_out,
                              '--minsize', str(minsize),
                              '--unoise_alpha', str(alpha),
                              '--id', str(clust_threshold),
                              '--relabel', 'ASV']

        # execute the vsearch UNOISE command
        run_subprocess(vsearch_unoise_cmd, dest_dir=denoise_parent, run_name=run_name, program='unoise',
                       auto_respond=settings['automate']['auto_respond'])

    else:

        # create a dictionary of the input and output files for each sample, to use further down when comparing read
        #   counts before and after denoising
        input_output_dict = {}

        # go through each Illumina sample...
        for illumina_file in illumina_files:

            # create an output file path to write denoised centroid sequences to
            denoise_output = add_prefix(file_path=illumina_file, prefix=DENOISE_PREFIX,
                                     dest_dir=denoise_path, action=None)

            # add this sample's input and output paths to the input/output dictionary
            input_output_dict.update({illumina_file: denoise_out})

            # assemble the vsearch UNOISE command for the command line
            vsearch_unoise_cmd = ['vsearch', '--cluster_unoise', illumina_file,
                                  '--centroids', denoise_output,
                                  '--minsize', str(minsize),
                                  '--unoise_alpha', str(alpha),
                                  '--id', str(clust_threshold)]

            # execute the vsearch UNOISE command
            run_subprocess(vsearch_unoise_cmd, dest_dir=denoise_parent, run_name=run_name, program='unoise',
                           auto_respond=settings['automate']['auto_respond'])


    ## LOG FILES THAT WERE NOT DENOISED ###############

    # if any of the input files were not Illumina sequences
    if len(non_illumina_files) > 0:

        # print a warning, and log these files to a .log file

        non_illumina_log = (denoise_parent / f'{run_name}_not-denoised').with_suffix('.log')

        print(f'WARNING. {len(non_illumina_files)} sequencing files that were provided to the function '
              f'{func.__name__} were not detected as Illumina sequences. UNOISE is only suitable for denoising '
              f'Illumina sequences, so these files were not denoised. See the output file\n'
              f'  {non_illumina_log.name}\n'
              f'for a list of the files that were not denoised.\n')

        with open(non_illumina_log, 'wt') as log_out:

            # write a header for the file
            log_out.write(f'bioinformatics run = {run_name}\n')
            log_out.write(f'The following sequences were determined to be non-Illumina sequencing files, and therefore '
                          f'were not denoised by UNOISE:\n')

            for err_file in non_illumina_files:
                logout.write(f'  {err_file}\n')

    # if all of the input reads are Illumina sequences, do nothing
    else:
        pass

    ## CONFIRM READ COUNTS DECREASED ##################

    # create a dictionary of unwanted read changes, either increase or no change of read counts after denoising
    read_count_change_issues = {'increase': {'num samples': 0,
                                             'samples': []},
                                'no change': {'num samples': 0,
                                              'sample': []}}

    # go through each pair of input and output files (i.e., for every sample)
    for file_in, denoise_out in input_output_dict.items():

        ## BEFORE DENOISING
        # instantiate counter to track the number of sequences before denoising
        read_count_before = 0

        # go through each record in this sample, and add to counter for each distinct read that is parsed
        with open(file_in) as input_fasta:
            for record in SeqIO.parse(input_fasta, 'fasta'):
                read_count_before += 1

        ## AFTER DENOISING
        # instantiate counter to track the number of sequences after denoising
        read_count_after = 0

        # go through each record in this sample, and add to counter for each distinct read that is parsed
        with open(denoise_out) as denoise_fasta:
            for record in SeqIO.parse(denoise_fasta, 'fasta'):
                read_count_after += 1

        ## COMPARISON
        # if there was no change in read count, this might be an issue
        if read_count_after == read_count_before:
            result = 'no change'

        # if there was an *increase* in read count, this is almost surely an issue
        elif read_count_after > read_count_before:
            result = 'increase'

        # if there was a decrease in read count, this would make sense, do nothing and move to next sample
        else:
            continue

        # if either read change error/warning conditions are met, then add to corresponding location in dictionary
        sample_id = get_sample_id(file_in)
        read_count_change_issues[result]['sample'].append(sample_id)
        read_count_change_issues[result]['num samples'] += 1

    ## ASSESS READ COUNT CHANGES

    # pull the number of samples for each situation, used several times below
    increase_count = read_count_change_issues['increase']['num samples']  # how many samples increased in read count?
    nochange_count = read_count_change_issues['no change']['num samples']  # how many samples didn't change read count?

    # if either issue contains samples...
    if (increase_count > 0) or (nochange_count > 0):

        # create a path for an output log file, with .json format
        change_error_path = (denoise_parent / f'{run_name}_read-count-issues').with_suffix('.json')

        # write the dictionary out to a .json file; write out both increase and no change, even if one has no samples
        with open(change_error_path, 'wt') as json_out:
            json.dump(read_count_change_issues, json_out)

        # calc percent of total input samples that had undesirable read count changes
        perc_issue = ((increase_count + nochange_count) / len(illumina_files))*100

        # print a warning before completing function
        print(f'WARNING. {increase_count + nochange_count} out of {len(illumina_files)} samples ({perc_issue:.2f}%) '
              f'had an unusual change in read count after denoising with vsearch:\n'
              f'  increase in read count  = {increase_count}\n'
              f'  no change in read count = {nochange_count}\n'
              f'Please see the samples belonging to each of these unusual changes in read counts in the summary '
              f'file:\n'
              f'  {change_error_path}')

    return None



def create_otu_fasta(query_fasta, reference_dir, otu_label='OTU', keep_abund=True, clust_threshold=None):

    # confirm that the input directory is a Path object
    is_pathclass(input_dir, exit_if_false=False)

    # import settings from the configuration .toml file
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    # if no clustering threshold is provided directly to the function, then look in the settings config file to
    #  find the clustering threshold to use
    if clust_threshold is None:
        clust_threshold = settings['otu_clustering']['min_threshold']
    else:
        pass

    # confirm that the clustering threshold is a float between 0 and 1; if not, try to scale to value between 0 and 1

    # if it between 0 and 1, continue without doing anything
    if 0 <= clust_threshold <= 1:
        pass

    # if it is not between 0 and 1...
    else:

        # see if creating a fraction percentage from the provided value will force between 0 and 1
        if 0 <= (clust_threshold/100) <= 1:
            clust_threshold = clust_threshold / 100

        # if it can't be automatically scaled to the correct float range, then return an error and exit
        else:
            err_msg = (f'ERROR. The provided clustering threshold of {clust_threshold} must be between 0 and 1, and '
                       f'could not automatically be scaled to a float. Please adjust this setting in the '
                       f'configuration to a value between 0 and 1, then proceed (e.g., a 97% clustering threshold '
                       f'should be entered as 0.97, not 97.\n')
            exit_process(message=err_msg)

    # create an output name for the otu fasta file, in same location as the query fasta file
    otus_out = query_fasta.parent / f'{run_name}_otus.fasta'
    otu_tab_out = query_fasta.parent / f'{run_name}_otu-tab.txt'

    # determine whether to keep the read abundance part of the header in the OTU .fasta output
    if keep_abund:
        abund_cmd = '--sizeout'
    else:
        abund_cmd = ''

    # cluster the reads in the query fasta, which contain all denoised reads in this bioinformatics run, across all
    #   samples
    vsearch_clust_cmd = ['vsearch', '--cluster_size', query_fasta, '--centroids', otus_out,
                         '--otutabout', otu_tab_out,'--id', clust_threshold, '', otu_label, abund_cmd]

    run_subprocess(vsearch_clust_cmd, dest_dir=query_fasta.parent, run_name=run_name,
                   program='uclust', separate_sample_output=True, auto_respond=settings['automate']['auto_respond'])

    return query_fasta.parent




def choose_representative(input_files, file_map):
    '''
    Decide which PacBio read to use as the representative read.

    :param input_files:
    :param file_map:
    :return:
    '''

    # create a BioSeq object to store the representative reads into
    # combine all representative reads into a single file

# def create_blast_db(config_dict, file_map, taxa_list=None):
#     '''
#     Create a reference dataset for BLAST+ search.
#
#     :param config_dict: configuration file dictionary; required to get the list
#     of reference databases to use, as defined by user in the configuration file.
#     :param file_map: file mapping, from mapping.py
#     :param taxa_list: option; list of taxa to include if wanting to limit search
#     to a specific taxonomic group
#     :return: returns path to the newly created db for blastn search; creates a custom
#     database using BLAST+ makeblastdb from the command line
#     '''
#     tax_settings = config_dict['taxonomy']
#     db_dir = file_map['config']['reference-db']
#
#     include_genbank = tax_settings['refdb']['genbank']['include']
#     include_unite = tax_settings['refdb']['unite']['include']
#     include_custom = tax_settings['refdb']['custom']['include']
#     include_maarjam = tax_settings['refdb']['maarjam']['include']
#
#     included_tag = ''
#     db_include_list = []
#     if include_genbank:
#         genbank_fasta = db_dir.glob('*genbank*')
#         db_include_list.append(genbank_fasta)
#         included_tag += 'G'
#
#     if include_unite:
#         unite_fasta = db_dir.glob('*unite*')
#         db_include_list.append(unite_fasta)
#         included_tag += 'U'
#
#     if include_custom:
#         custom_fasta = db_dir.glob('*custom*')
#         db_include_list.append(custom_fasta)
#         included_tag += 'C'
#
#     if include_maarjam:
#         maarjam_fasta = db_dir.glob('*maarjam*')
#         db_include_list.append(maarjam_fasta)
#         included_tag += 'M'
#
#     custom_search_records = []
#     for db in db_include_list:
#         for record in SeqIO.parse(db, 'fasta'):
#             if taxa_list is None:
#                 custom_search_records.append(record)
#             else:
#                 pass  # NEED TO SEE HOW TO GET ONLY SPECIFIC TAXONOMY, NOT SURE HOW FORMATTED
#
#     tax_out_dir = file_map['pipeline-output']['taxonomy']
#
#     search_date = datetime.today.strftime('%Y-%M-%d')
#     output_basename = f'custom-ref-{included_tag}_{search_date}'
#     output_path = (tax_out_dir / output_basename).with_suffix('.fasta')
#
#     SeqIO.write(custom_search_records, output_path, 'fasta')
#
#     output_db = tax_out_dir / f'{output_basename}_blast'
#     blast_cmd = ['makeblastdb', '-in', output_path, '-title', output_db, '-dbtype', 'nucl',
#                  output_db]
#
#     run_subprocess(blast_cmd, dest_dir = output_db)
#
#     return output_db

def assign_taxonomy(otu_fasta, query_fasta, reference_dir, method=None):

    # import settings from the configuration file
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    # define list of available methods
    available_methods = ['amptk', 'rdp', 'blastn']

    # if a method to assign taxonomy isn't provided to the function, then look in the configuration file
    if method is None:
        method = settings['taxonomy']['method']

    # based on the specified method, compile the command for the CLI
    if method == 'blastn':
        pass
        # ref_db = create_blast_db(config_dict, file_map, taxa_list=None)
        #
        # blast_out = (tax_output / f'{run_name}').with_suffix('.txt')
        # blast_cmd = ['blastn', '-query', ref_db, '-out', blast_out]
        #
        # run_subprocess(blast_cmd, dest_dir = tax_output)

    elif method == 'rdp':
        pass


    ## AMPTK

    elif method == 'amptk':

        # import the amptk databases based on type of sequences to assign taxonomy to
        # hm easier in theory; could be a mix, filename wouldn't indicate, only read headers would
        amptk_db_install_cmd = ['amptk', 'install', '-i', 'ITS']
        run_subprocess(amptk_db_install_cmd, dest_dir=otu_fasta.parent, run_name=run_name, program='amptk-db',
                       separate_sample_output=True, auto_respond=False)

        amptk_cmd = ['amptk', 'taxonomy', '-i', otu_fasta, '-f', SOMETHING, '-m',
                     SOMETHING, '-d', 'ITS1']

    else:
        error_msg = (f'ERROR. The provided method of assigning taxonomy, {method}, is not currently '
                     f'recognized among the list of available methods. Please choose one of the following '
                     f'accepted methods:\n')
        print_indented_list(available_methods)
        exit_process(error_msg)

    error_msg = (f'ERROR. The provided method of assigning taxonomy, {method}, is not currently '
                 f'available for use in the pipeline.\n')
    exit_process(error_msg)

    return None

def cluster_reads(input_files, output_dir, reference_dir, clust_threshold, clust_method, group=True):

    ## IMPORT SETTINGS FOR PIPELINE ###################

    # import configuration settings
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']


    ## CREATE OUTPUT FILE PATHS #######################

    # create a directory for all clustering output
    clust_parent = mkdir_exist_ok(new_dir=output_dir)

    # create a directory within the main clustering output for this particular pipeline run
    clust_output = mkdir_exist_ok(new_dir=f'./{CLUSTER_PREFIX}_{run_name}', parent_dir=clust_parent)


    ## FLATTEN INPUT FILE LIST ########################

    # if input files are directories (pacbio), create a list of all sequence files within all input directories
    updated_file_list = []
    for file_in in input_files:
        if file_in.is_dir():
            for fasta_file in file_in.glob(SEQ_FILE_GLOB):
                updated_file_list.append(fasta_file)
        else:
            updated_file_list.append(file_in)


    ## COMBINE SAMPLE SEQUENCES PRE-CLUSTERING ########

    # if group=True, gather all like sequences together into a single .fasta before clustering
    if group:

        ## SORT BY DNA REGION ##

        # create a dictionary where the key is the DNA region that will match the file tags and values are empty list
        clust_regions = [region_tag.lower() for region_tag in POST_ITSX_SUFFIXES.values() if
                            not (region_tag in ['5_8S', 'LSU'])]
        input_files_sorted = {region: [] for region in clust_regions}

        # create a regex that will search for any of the regions in the input .fasta files
        clust_regions_re = '|'.join(clust_regions)

        # go through the input files and sort files by region
        for file_in in updated_file_list:
            wanted_region_found = re.search(clust_regions_re, file_in.name, re.I)
            if wanted_region_found:
                input_files_sorted[wanted_region_found.group(0).lower()].append(file_in)
            else:
                continue

        ## COMBINE BY DNA REGION ##

        # create an output file path for these grouped sequences in the bioinfo run clustering directory
        combined_seq_output = mkdir_exist_ok(
            new_dir=f'pre-{CLUSTER_PREFIX}_grouped-sequences',
            parent_dir=clust_output,
        )

        # create a dict of the combined sequence files created below, will be input for vsearch clustering
        seq_files_to_cluster = {}

        # go through the list of sorted input files by DNA region...
        for dna_region, region_files in input_files_sorted.items():

            # gather all sequence records into a list of sequence records for this DNA region
            region_seq_records = []

            # create an output file path for these grouped sequences in the bioinfo run clustering directory
            region_seq_output = combined_seq_output / f'pre-{CLUSTER_PREFIX}_{run_name}_{dna_region}.fasta'

            # add this output file path to the list of combined sequence files to cluster
            seq_files_to_cluster.update({dna_region: region_seq_output})

            # go through each sample's sequence file for this DNA region...
            for sample_fasta in region_files:

                # get the sample ID from the input sequence file name
                sample_id = get_sample_id(sample_fasta, platform='itslsu').replace(f'{run_name}_', '')

                # go through each record (read) in this sample's sequence file...
                for record in SeqIO.parse(sample_fasta, 'fasta'):

                    # update the read header to include sample=<sample-id>; as first item, which is required
                    #  in order for vsearch clustering to recognize the sample ID for the OTU table
                    # otu=<read-id> is also required in order for the sequence read IDs to be used as the OTU name
                    updated_record_header = f'sample={sample_id};otu={record.id}'

                    # create a new sequence record using this updated header
                    seq_record_updated = SeqRecord(
                        id = updated_record_header,
                        name = updated_record_header,
                        description = updated_record_header,
                        seq = record.seq,
                    )

                    # add this updated sequence record to the list of sequence records for this DNA region
                    region_seq_records.append(seq_record_updated)

            # once new sequence records have been created for all samples for this DNA region, write out .fasta
            SeqIO.write(region_seq_records, region_seq_output, 'fasta')

    # if group=False, cluster reads only within a sequence file, not among input sequence files
    else:

        # COME BACK AND ADD THIS #

        seq_files_to_cluster = {'NA': updated_file_list}


    ## CLUSTER SEQUENCE FILE READS W/ VSEARCH #########

    ## CLUSTER BY DNA REGION ##

    for dna_region, clust_input in seq_files_to_cluster.items():

        # create output file paths for the centroid sequences and OTU table
        centroid_seqs_out = clust_output / f'{CLUSTER_PREFIX}_{run_name}_{dna_region}_centroids-{str(clust_threshold)}.fasta'
        otu_table_out = clust_output / f'{CLUSTER_PREFIX}_{run_name}_{dna_region}_otu-table-{str(clust_threshold)}.txt'

        # compile the clustering command for vsearch
        vsearch_clust_cmd = ['vsearch', '--cluster_fast', clust_input,
                             '--id', str(clust_threshold), '--iddef', str(clust_method),
                             '--centroids', centroid_seqs_out, '--clusterout_id',
                             '--otutabout', otu_table_out]

        # execute the compiled command for vsearch clustering
        run_subprocess(
            cli_command_list = vsearch_clust_cmd,
            dest_dir = clust_output,
            run_name = run_name,
            program = 'vsearch-clust',
            separate_sample_output = True,
            auto_respond = settings['automate']['auto_respond'],
        )

    # return the output path for clustered sequences from this bioinformatics run
    return clust_output
