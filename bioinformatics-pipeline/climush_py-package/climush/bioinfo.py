from pathlib import Path
from Bio.Seq import Seq
from datetime import datetime
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from climush.constants import NOPHIX_PREFIX, SEQ_FILE_GLOB, TRIMMED_PREFIX, DEREP_PREFIX, QUALFILT_PREFIX
from climush.utilities import *

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

            fwd_read_out = add_prefix(file_path=file, prefix=TRIMMED_PREFIX, dest_dir=trim_path, action='None')
            rev_read_out = add_prefix(file_path=paired_dict[file], prefix=TRIMMED_PREFIX, dest_dir=trim_path, action='None')

            cutadapt_cmd = ['cutadapt', '--report=minimal',
                            '-a', fwd_primer, '-A', rev_primer,
                            '-n', '2', '-e', str(max_err),
                            '-o', fwd_read_out, '-p', rev_read_out,
                            file, paired_dict[file]]

            run_subprocess(cutadapt_cmd, dest_dir=trim_primers_parent)

    else:

        for file in input_dir:

            trim_out = add_prefix(file_path=file, prefix=TRIMMED_PREFIX, dest_dir=trim_path, action='None')

            if cutadapt_settings['keep_untrimmed']:
                notrim_out = add_prefix(file_path=file, prefix=flip_prefix(TRIMMED_PREFIX),
                                        dest_dir=notrim_path, action='None')
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
    with open((trim_primers_parent / 'cutadapt.out'), 'r') as fin:
        cutadapt_df = pd.read_table(fin)
        sum_in = cutadapt_df['in_reads'].sum(0)
        sum_out = cutadapt_df['out_reads'].sum(0)
        percent_lost = ((sum_in - sum_out) / (sum_in)) * 100
        if percent_lost > max_untrimmed:
            msg = f'After primer trimming, {percent_lost:.2f}% of the input reads were lost, which is ' \
                  f'above the user-defined maximum threshold of {max_untrimmed}%.\n'
            exit_process(message=msg)
        else:
            print(f'{percent_lost:.2f}% of input reads were lost to primer trimming. This is above the user-provided '
                  f'input of {max_untrimmed}%, so proceeding to next step...\n')

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
                                           dest_dir=qfilt_path, action='None')
            nofilt_out = add_prefix(file_path=file, prefix=flip_prefix(QUALFILT_PREFIX),
                                        dest_dir=nofilt_path, action='None')

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
                                       dest_dir=qfilt_path, action='None')
            nofilt_fwd_out = add_prefix(file_path=file, prefix=flip_prefix(QUALFILT_PREFIX),
                                        dest_dir=nofilt_path, action='None')

            qfilt_rev_out = add_prefix(file_path=paired_dict[file], prefix=QUALFILT_PREFIX,
                                       dest_dir=qfilt_path, action='None')
            nofilt_rev_out = add_prefix(file_path=paired_dict[file], prefix=flip_prefix(QUALFILT_PREFIX),
                                        dest_dir=nofilt_path, action='None')

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

        itsx_output_basename = add_prefix(file_path=file, prefix=ITSX_PREFIX, dest_dir=itsx_path).with_suffix('')

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