from mapping import filepath_map as fpm
import matplotlib.pyplot as plt
from climush.utilities import *
import math
import numpy as np

# CREATE A LIST OF SAMPLES FROM ORD IN 2024-04 SPOROCARPS

# fl_its1_path = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/pipeline-output/derep-subregions/derep02_bioinfo-test')


base_path = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/pipeline-output/separated-subregions/itsx_psf2404-2407')
ord127_test = base_path / 'itsx_pacbio_sporocarp-f_2024-04_ORD127'

ord127_fullits = next(ord127_test.glob('*full_ITS*'))

read_counts = {}
with open(ord127_fullits, 'rt') as fin:
    for line in fin.readlines():
        if line.startswith('>'):
            read_id = line.split(';')[0]
            read_count = int(line.split(';')[-1].replace('full-len_copies=', ''))
            read_counts[read_id] = read_count
        else:
            continue

read_count_arr = np.array(list(read_counts.values()))
read_count_diffs = np.diff(read_count_arr)

fig, (ax1, ax2) = plt.subplots(1,2)

ax1.bar(x=[(f'ASV{x}-{x + 1}') for x in np.arange(1, len(read_count_arr) + 1)], height=read_count_diffs)
ax1.yaxis.set_inverted(True)  # invert y-axis
[ax.set_visible(False) for ax in ax1.get_xticklabels()]  # remove x-axis labels
ax1.axhline(y=np.mean(read_count_diffs), linestyle='dashed', color='red')  # calc mean diff, include as horizontal line
ax1.set_xlabel(f'mean diff: {np.mean(read_count_diffs)} nt\nASV count: {len(read_count_arr)}')

ax2.bar(x=[f'ASV{x}' for x in np.arange(len(read_count_arr))], height=read_count_arr)

plt.show()










fl_samples = {}
total_sample_count = 0
for file in ord127_test.glob('*ORD*full_ITS*'):

    total_sample_count += 1
    sample_file_empty = False

    # get the sample ID; create empty lists for each ITS1 read length and the number of exact copies of full-length
    #  reads for the sequence
    sample_id = re.search(r'(?<=\_)ORD\d{1,3}(?=\.)', file.name, re.I).group(0)
    read_copies_its1 = []
    read_lengths = []

    # open the ITS1 file from the ITSx output
    with open(file, 'rt') as fin:
        for line in fin.readlines():
            # look at the header line to gather information about the sequence
            if line.startswith('>'):
                # get the number of full-length read copies that this ITS1 read was extracted from
                try:
                    identical_its1 = re.search(r'(?<=size=)\d{1,4}$', line, re.I).group(0)
                    read_copies_its1.append(int(identical_its1))
                except AttributeError:
                    sample_file_empty = True
                    print(f'Sample does not contain read copy info: {sample_id}')

            else:
                # get the read length for this ITS1 sequence
                try:
                    itsx_length = len(line)
                    read_lengths.append(int(itsx_length))  # add as int
                except AttributeError:
                    continue

    # add this info for the sample to the fl_samples dictionary
    if sample_file_empty:
        continue
    else:
        fl_samples.update({sample_id: {'read_counts_its1':read_copies_its1,
                                       'read_lengths': read_lengths}})

print(f'There are {len(fl_samples)} sporocarp samples out of {total_sample_count} from the Ordway-Swisher site in '
      f'Florida in the sequencing run submitted in 2024-04.\n')

threshold_range = [40,60]
one_top_read = []
two_top_reads = []
many_top_reads = []
borderline = []
for sample in fl_samples:
    sample_read_counts = fl_samples[sample]['read_counts_its1']
    total_read_count = sum(sample_read_counts)
    sample_read_counts.sort(reverse=True)
    for c,count in enumerate(sample_read_counts):
        seq_read_count_per = math.ceil((count / total_read_count) * 100)  # round up
        if seq_read_count_per >= 70:
            one_top_read.append(sample)
            break
        else:
            next_read_count_per = math.ceil((sample_read_counts[c+1] / total_read_count) * 100)  # round up

            # if both this current read and the next comprise 40-60% of total reads individually...
            if (threshold_range[0] <= seq_read_count_per <= threshold_range[1]) and \
                    (threshold_range[0] <= next_read_count_per <= threshold_range[1]):
                two_top_reads.append(sample)
                break

            # if at least this current read or the next comprise 40-60% of the total reads individually...
            elif (threshold_range[0] <= seq_read_count_per <= threshold_range[1]) or \
                    (threshold_range[0] <= next_read_count_per <= threshold_range[1]):
                borderline.append(sample)
                break

            # if neither this current read or the next comprise 40-60% of the total reads individually...
            else:
                many_top_reads.append(sample)
                break

many_top_reads.sort()

print(f'Samples with one top read: {len(one_top_read)}\n'
      f'Samples with two top reads: {len(two_top_reads)}\n'
      f'Samples with no top reads: {len(many_top_reads)}\n'
      f'Samples with borderline number of reads: {len(borderline)}\n')

for s in many_top_reads:
    counts = fl_samples[s]['read_counts_its1']
    print(f'{s}: {counts}\n')

# calc diff between consecutive read sizes of no stand-out top read samples:
# set the number of rows and columns for the grid-arranged subplots
fig_ncol = 4
fig_nrow = math.ceil(len(many_top_reads) / fig_ncol)  # round up to nearest whole num
fig = plt.figure(figsize=(8.5, 11), layout='constrained')
gsp = fig.add_gridspec(fig_nrow, fig_ncol)
axs = gsp.subplots(sharex=False, sharey=False)
# fig, axes = plt.subplots(fig_nrow, fig_ncol, figsize=(8.5,11))
fig.suptitle('Difference in Consecutive Read Counts for Samples with Multiple Top Reads', fontsize=15)
# fig.ylabel('Difference in Read Count (nt)')
# fig.xlabel('ASV Read Count Comparison')

# iterate through samples with too many top reads, use index for plotting
current_row = 0
for s,sample in enumerate(many_top_reads):

    current_col = s % fig_ncol

    if (current_col == 0) and (s > 0):
        current_row += 1

    current_axis = axs[current_row, current_col]

    # create an array of read counts for this sample; already in desc order
    count_arr = np.array(fl_samples[sample]['read_counts_its1'])
    # calculate consecutive difference in read counts
    count_diffs = np.diff(count_arr)
    count_diffs_mean = np.mean(count_diffs).round(2)

    # plot a figure of the difference in consec read counts for this sample in its grid
    current_axis.bar(x=[(f'ASV{x}-{x+1}') for x in np.arange(1,len(count_diffs)+1)], height=count_diffs)
    current_axis.set_title(sample)  # set subplot title to sample ID
    current_axis.yaxis.set_inverted(True)  # invert y-axis
    [ax.set_visible(False) for ax in current_axis.get_xticklabels()]  # remove x-axis labels
    current_axis.axhline(y=count_diffs_mean, linestyle='dashed', color='red')  # calc mean diff, include as horizontal line
    current_axis.set_xlabel(f'mean diff: {count_diffs_mean} nt\nASV count: {len(count_arr)}')

# remove plots in gridspace that are empty
for ax1 in range(axs.shape[0]):
    for ax2 in range(axs.shape[1]):
        # not a great fix, but checks if the axes obj has a title, if not its probably empty so remove grid
        if axs[ax1,ax2].get_title() == '':
            axs[ax1,ax2].set_visible(False)
        else:
            continue

# once all plotted, show fig and close out
# plt.show()
plt.savefig('../ord_test/ord_multi-top-reads.pdf')
plt.close()

# def get_representative_read(input_files, file_map, threshold_range=[40,60]):
#     '''
#     Determine the representative read for sporocarp PacBio sequences.
#
#     :param threshold_range: range of values to use as the minimum and maximum for the percent of reads
#     represented by two sequences
#     :return: write out representative read to a single fasta file for all sequences in this bioinformatics run
#     '''
#
#     input_files = list(fl_its1_path.glob('*ORD*ITS1*'))
#     file_map = fpm
#     threshold_range = [40,60]
#
#     all_reads = {}
#     for file in input_files:
#         # sample_id = get_sample_id(file_path=file)
