from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import re

OUTLIER_QUANTILE = 0.7
pd.set_option('display.max_columns', 10)
pd.set_option('display.max_rows', 10)

def get_study_site(sample_id):
    '''
    Get the formatted string of the study cite from a sample ID.

    :param sample_id: sample ID from the ITS variant dataframe
    :return: string of the spelled out study site
    '''

    # if the sample ID starts with a number, it is from oregon
    if re.search(r'^\d', sample_id):
        return 'Oregon'

    # if the sample ID does not start with a number...
    else:

        # try to extract a three-letter climush study site prefix
        climush_site_match = re.search(r'^[A-Z]{3,4}', sample_id)
        if climush_site_match:

            # if the match is to HAD, it is Heather Dawson's project (not climush)
            if climush_site_match.group() == 'HAD':
                or_heather = {
                    'collector': 'Heather Dawson',
                    'state': 'Oregon',
                    'country': 'USA',
                }
                return or_heather

            # if the match is to KMS, it is undergrad Kyla Schmitt's project (not climush)
            elif climush_site_match.group() == 'KMS':
                or_kyla = {
                    'collector': 'Kyla Schmitt',
                    'state': 'Oregon',
                    'country': 'USA',
                }
                return or_kyla

            # otherwise, it should be a climush sample
            else:
                # get just the string for the study site that matches the regex
                climush_site_str = climush_site_match.group()

                # create a dictionary of each of the top ITS variant number sample ID climush sites
                climush_site_convert = {
                    'CDR': {
                        'name': 'Cedar Creek Biological Station',
                        'state': 'Minnesota',
                        'country': 'USA',
                    },
                    'KON': {
                        'name': 'Konza Prairie Biological Station',
                        'state': 'Kansas',
                        'country': 'USA',
                    },
                    'NWT': {
                        'name': 'Niwot Ridge LTER',
                        'state': 'Colorado',
                        'country': 'USA',
                    },
                }

                # select the climush site info from this dictionary that matches the climush site string
                return climush_site_convert[climush_site_str]



## FILE PATHS ##########################################################################################################

# PROJECT PATHS #

# path to the main directory for the ITS variant analysis
itsvar_main_path = Path('/Users/carolyndelevich/main/github_repos/climush/data-analysis/its-variants/')

# path to the data/ subdirectory in the ITS variant analysis directory
itsvar_data_path = itsvar_main_path / 'data'

# path to the figures/ subdirectory in the ITS variant analysis directory
itsvar_figs_path = itsvar_main_path / 'figures'
itsvar_figs_path.mkdir(exist_ok=True)

# INPUT PATHS #

# path to the merged and cleaned dataframe produced by the script its-variants_data-cleaning.py
# THIS IS AN ISSUE WITH AUTO-DETECTING WHICH TO USE; USE LATEST, BUT DIDNT' WANT TO CODE AUTO-DETECT LATEST
itsvar_counts_path = itsvar_data_path / f'pacbio_its-variants_2025-06-24.csv'

# OUTPUT PATHS #

# output path for the histogram produced here
date_suffix = datetime.now().strftime('%Y-%m-%d')
hist_output_path = itsvar_figs_path / f'its-variants_histogram_{date_suffix}.png'

########################################################################################################################


## IMPORT DATA #########################################################################################################

# import the cleaned dataframe produced in the previous script, its-variants_data-cleaning.py
itsvar_counts_df = pd.read_csv(itsvar_counts_path)

########################################################################################################################


## IDENTIFY OUTLIERS ###################################################################################################

# GET CUTOFF #

# determine what will be considered an outlier so I know where to create the y-axis break, below
# because I'm plotting a histogram, outliers are based on counts of each value in the dataset, not the value itself

# create an array of value counts in the ITS variant dataset
itsvar_valcounts = np.unique(itsvar_counts_df['its_variant_count'], return_counts=True)

# subset this two-dim array to get only the counts of each value (i.e., the y-axis in the histogram)
itsvar_counts = itsvar_valcounts[1]

# use the constant OUTLIER_QUANTILE (defined at top of this script) as the cut-off for what is considered an outlier
#    from the counts of each ITS variant value
outlier_cutoff = np.quantile(itsvar_counts, OUTLIER_QUANTILE)
# print(f'The {OUTLIER_QUANTILE} quantile used as the outlier cutoff results in any ITS variant count '
#       f'greater than {}')

# FIND OUTLIER IDS #

# what is the largest number of ITS variants?
itsvar_count_max = itsvar_counts_df['its_variant_count'].max()
itsvar_count_min = itsvar_counts_df['its_variant_count'].min()
itsvar_count_mean = itsvar_counts_df['its_variant_count'].mean()
itsvar_count_std = itsvar_counts_df['its_variant_count'].std()
itsvar_count_median = itsvar_counts_df['its_variant_count'].median()
itsvar_count_mode = itsvar_counts_df['its_variant_count'].mode().iloc[0]
itsvar_count_q95 = np.quantile(itsvar_counts_df['its_variant_count'], 0.95)
print(f'Summary of the number of ITS variants per sequenced sporocarp:\n\n'
      f'   max number of variants = {itsvar_count_max}\n'
      f'   min number of variants = {itsvar_count_min}\n\n'
      f'   mean number of variants (+/- stdev) = {itsvar_count_mean:.2f} +/- {itsvar_count_std:.2f}\n'
      f'   median number of variants = {itsvar_count_median:.2f}\n\n'
      f'   most common number of variants = {itsvar_count_mode}\n')

# what IDs / species have the most ITS variants?
itsvar_max_sporocarp = itsvar_counts_df[itsvar_counts_df['its_variant_count'] == itsvar_count_max]
print(f'The sporocarp with the highest number of ITS variants ({itsvar_count_max}):\n'
      f'   {itsvar_max_sporocarp["identification"]}\n')


## SUMMARIZE ONLY TOP ITS VARIANT COUNTS ##

# subset the original dataframe
# itsvar_counts_top = itsvar_counts_df[itsvar_counts_df['Nr dominant variants'] >= (itsvar_count_mean*2)][['identification','sample_id','Nr dominant variants']].sort_values('Nr dominant variants', ascending=False).reset_index(drop=True)
itsvar_counts_top = itsvar_counts_df[itsvar_counts_df['its_variant_count'] >= (itsvar_count_mean*2)].sort_values('its_variant_count', ascending=False).reset_index(drop=True)

# list of top sample IDs
top_sampids = itsvar_counts_top['sample_id'].to_list()
top_sampids.sort()

# create an empty list to which the formatted info string of each sporocarp with a top number of variants is added below
itsvar_counts_top_fmt = []

# pull out the sporocarp ID, the sample ID, and the number of dominant variants for each top sporocarp to
#   format them together for a well-formatted print-out of these values
for i in np.arange(itsvar_counts_top.shape[0]):

    # pull out the strings from the filtered top variant table
    sp_id = itsvar_counts_top['identification'].iloc[i]         # species ID
    samp_id = itsvar_counts_top['sample_id'].iloc[i]            # sample ID
    var_ct = itsvar_counts_top['its_variant_count'].iloc[i]  # variant count

    # put together this information into a single string for this sporocarp then add to the list
    fmt_str = f'{samp_id} | {sp_id} [n = {var_ct}]'
    itsvar_counts_top_fmt.append(fmt_str)

# join the strings of each individually formatted sporocarp info into a single string, where each sporocarp
#   info line is indented and one their own lines
itsvar_counts_top_fmt = '\n   '.join(itsvar_counts_top_fmt)  # replace the individual strings with single string

# print out a summary of these sporocarps that have 2x the mean number of ITS variant counts
print(f'The sporocarps that had twice the number of the average ITS variant count of '
      f'{itsvar_count_mean:.2f} (> {(itsvar_count_mean*2):.2f} variants):\n'
      f'   {itsvar_counts_top_fmt}')

########################################################################################################################

## BREAK Y-AXIS ########################################################################################################

# this code was generated from the following matplotlib tutorial:
#   https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html


# # create a function so parameter changes are easier to re-run
# def ybreak_histogram(dataframe=itsvar_counts_df, x_val='its_variant_count', xlims=[0,13], majority_ylims=[448,544], outlier_ylims=[0,22], save_fig=False, **kwargs):
#     '''
#
#     :param dataframe:
#     :param save_fig: True/False; if True, the resulting figure is saved to a file; if
#     False, the figure is displayed but not saved to a file
#     :param kwargs:
#     :return:
#     '''


dataframe=itsvar_counts_df
x_val='its_variant_count'
xlims=[0,13]
majority_ylims=[448,544]
outlier_ylims=[0,22]
save_fig=False


## SET UP SUBPLOTS ##

# create two subplots, where the upper subplot shows the majority of the data and the
#  lower subplot shows a 'zoomed-in' display of the outliers / tail of the distribution
fig, (ax_upper, ax_lower) = plt.subplots(
    nrows = 2,
    ncols = 1,
    sharex=True,
)

# increase the space between the upper and lower subplots
fig.subplots_adjust(hspace=0.1)

## PLOT DATA ##

# plot the same histogram in both subplot figures
hist_upper = sns.histplot(data=dataframe, x=x_val, ax=ax_upper, discrete=True)
hist_lower = sns.histplot(data=dataframe, x=x_val, ax=ax_lower, discrete=True)

## SET AXES LIMITS ##

# change the y-limit of the upper subplot to 'zoom-in' on the majority of the data (largest bins)
ax_upper.set_ylim(
    majority_ylims[0],
    majority_ylims[1],
)

# change the y-limit of the lower subplot to 'zoom-in' on the outliers of the data (tail of distribution)
ax_lower.set_ylim(
    outlier_ylims[0],
    outlier_ylims[1],  # make this a little larger than what your last tick will be so there's room for the break line
)

# change the x-limits for aesthetics (must do both as x-axis is shared, i.e., scaled same)
ax_upper.set_xlim(xlims)
ax_lower.set_xlim(xlims)

## FORMAT SPINES ##

# hide the spines between the upper and lower subplots
ax_upper.spines.bottom.set_visible(False)  # bottom for subplot on top (majority)
ax_lower.spines.top.set_visible(False)     # top for subplot on bottom (outliers)

# hide the absolute top (top of the top subplot) spine and absolute right (right of both subplots) spine
ax_upper.spines.top.set_visible(False)    # top spine of top subplot
ax_upper.spines.right.set_visible(False)  # right spine of top
ax_lower.spines.right.set_visible(False)  # right spine of bottom

## FORMAT TICKS ##

# only include ticks on the very bottom of the plot, i.e., the bottom ticks of the bottom subplot
ax_upper.tick_params(
    axis='x',
    which='both',
    top=False,
    bottom=False,
    labelbottom=False
)

# change the tick labels and ticks so each value is shown (1-12) and not just even numbers
ax_lower.set_xticks(
    ticks=np.arange(xlims[0],xlims[1]),
    labels=np.arange(xlims[0],xlims[1]),
)

## MARK BREAK IN Y-AXIS ##

# create cut-out slanted lines that represent the y-axis break

# proportion of vertical to horizontal extent of the slanted line
d = 0.5

# settings to use to draw the break line
kwargs = dict(
    marker=[(-1,-d), (1,d)],  # marker style string; x/y pairs used for Path vertices w/ center at 0
    markersize=12,            # size of the marker in points (float)
    linestyle='none',         # linestyle
    color='k',                # color of the line; k = black (shorthand notation)
    markeredgecolor='k',      # marker edge color; k = black (shorthand notation)
    markeredgewidth=1,        # marker edge width in points (float)
    clip_on=False,            # when False, artists will be visible outside the Axes
)

# transform the ax1 and ax2 axes using the cut-out slanted line break settings defined above
# ax1.plot([0,1], [0,0], transform=ax1.transAxes, **kwargs)
# ax2.plot([0,1], [1,1], transform=ax2.transAxes, **kwargs)

ax_upper.plot(0, 0, transform=ax_upper.transAxes, **kwargs)
ax_lower.plot(0, 1, transform=ax_lower.transAxes, **kwargs)

## AXES LABELS ##

# remove all subplot labels
ax_upper.set_ylabel('')
ax_lower.set_ylabel('')
ax_upper.set_xlabel('')
ax_lower.set_xlabel('')

# add whole-figure labels for y-axis, x-axis, and figure title
fig.supxlabel('Number of ITS Variants')                                    # x-axis label
fig.supylabel('Number of Sequenced Sporocarps')                            # y-axis label
fig.suptitle('Distribution of the Number of\nITS Variants per Sporocarp')  # figure title

## BAR LABELS ##

# add bin/bar heights to the upper subplot; don't need to exclude zeros due to change in y-axis lims (starts at 450)
ax_upper.bar_label(hist_upper.containers[0], padding=3)

# create a list of bar labels (bin height) for the lower subplot, removing any zero values so they aren't displayed
lower_barlabels = [h2 if (h2 := v2.get_height()) != 0 else '' for v2 in hist_lower.containers[0]]
ax_lower.bar_label(hist_lower.containers[0], labels=lower_barlabels, padding=3)

## ANNOTATIONS ##

# transform the x-axis and y-axis into standardized units for locating arrows and boxes
## NEED TRANSFORMATION HERE

# set style for the annotation text boxes
bbox_settings = {
    'boxstyle': 'round',
    'fc': '0.8',
}

# set style for the arrows
arrow_settings = {
    'arrowstyle': '->',
    'connectionstyle': 'angle,angleA=0,angleB=90,rad=10',
}

# create label for the sporocarp with the most ITS variants
most_itsvars_row = dataframe[dataframe[x_val] == itsvar_count_max]
most_itsvars_spid = most_itsvars_row['identification'].iloc[0]
most_itsvars_sitedict = get_study_site(most_itsvars_row['sample_id'].iloc[0])
most_itsvars_annotation = (f'sporocarp ID: {most_itsvars_spid}\n'
                           f'origin: {most_itsvars_sitedict["name"]}\n'
                           f'{most_itsvars_sitedict["state"]}, {most_itsvars_sitedict["country"]}')

# label the sporocarp with the most ITS variants
offset = 72
ax_lower.annotate(
    text=most_itsvars_annotation,                         # the text for the annotation
    xy=(itsvar_valcounts[0][-1],itsvar_valcounts[1][-1]), # the x-y coordinate of the annotation
    xytext=(-2*offset, offset),                           # position to place the text
    # textcoords='',                                  #
    bbox=bbox_settings,
    arrowprops=arrow_settings,
)



## SAVE OR DISPLAY FIG ##

if save_fig:
    plt.savefig(hist_output_path, dpi=300)
else:
    plt.show()

plt.close()

# return None

########################################################################################################################

# ybreak_histogram()
