import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
# choose color palette based on the quality score folder
palette_qualscore = sns.color_palette("rocket", 3)

if quality_score == '20':
    qualscore_color = palette_qualscore[0]
elif quality_score == '30':
    qualscore_color = palette_qualscore[1]
else:
    qualscore_color = palette_qualscore[2]

def single_panel_histogram(variable, caption, labels, output_tag, annotate=True):
    '''
    Plot + save a single panel histogram.

    Plots a histogram with summary statistics and a caption for a single grouping
    of samples (i.e., a single histogram panel). Use multi_panel_histogram() if
    you want multiple histograms plotted together.
    :param variable: a list, Series, or 1-d array of values to plot
    :param caption: caption to use for the figure, as a string (f-string included)
    :param labels: dictionary of labels/values to use for the x-axis, y-axis, title,
    and sample_size, where the key is the label category (x-axis, y-axis, title, sample_size)
    and the value is a string to use for the label; title will be added to a pre-defined title
    that includes the quality score and sequencing run, followed by a colon (e.g., 'Q40 Reads (sp23):');
    sample_size can be an int or string, will not be used mathematically
    :param output_tag: suffix used at the end of the plot file name that
    describes the contents of the figure; should be short, will individual words or abbreviated words
    separated by hyphens; should not include filetype
    :param annotate: Boolean; whether to annotate the plot with summary statistics, which
    includes the mean, mode, min, max, number of values (length of the input variable), and
    sample size (number of sporocarp samples described by the variable; depending on the variable,
    number of values and sample size may be equivalent)
    :return: saves plot to the summary folder of this sequencing run
    '''

    fig_width = 7  # define first, as it's used in estimate_caption_adj
    caption_text_size = 10  # also need caption text size to estimate caption adjustment

    if caption is None:
        caption_adj = 0
    else:
        def estimate_caption_adj():
            '''
            Estimates the adjustment to the size of the figure required
            to accommodate a caption of a given length.
            :return: float; size of adjustment needed in normalized
            figure coordinates
            '''
            num_char = len(caption)

            # calculate line height based on font size using golden ratio, convert to inch
            line_height = (caption_text_size*1.61803399) / 72

            # calculate characters per line (CPL) for this font size and figure width
            px_width = fig_width * 72  # convert inches to points
            char_constant = 1.76  # estimated for matplotlib default font, DejaVu Sans Serif
            cpl = px_width * (char_constant/caption_text_size)

            num_lines = np.ceil(num_char/cpl)  # round up
            print(f"The caption is {num_char} characters long and had a cpl of {cpl}; therefore"
                  f"the number of required lines is {num_lines}")
            return num_lines * line_height
        caption_adj = estimate_caption_adj()
        print(f"Adjusting the figure by {caption_adj} inches.")

    fig_height = 4.5 + caption_adj

    fig = plt.figure(figsize=(fig_width, fig_height))  # use subplots to make room for the caption below

    ax1_stop = 4  # set up ax1 to occupy 3/4 of figure space; will switch to 4/4 (all) if annotations are turned off

    # annotate first; sets subplot spacing for histogram based on T/F for annotations
    if annotate:  # annotate unless set to False
        ax2 = fig.add_subplot(1,5,5)  # create subplot to add text annotations to, will span final column
        ax2.set_axis_off()  # remove plot elements (e.g., tics, spines)
        ax2.annotate(f"{labels['x-axis']} Summary\n", xy=(0, 0.95), ha='left', va='top', fontsize=10, linespacing=1.75)
        summary_stats_str = f"mean: {np.mean(variable):.1f}\n" \
                            f"mode: {stats.mode(variable)[0][0]} (n = {stats.mode(variable)[1][0]})\n" \
                            f"min:  {np.min(variable)}\n" \
                            f"max:  {np.max(variable)}\n" \
                            f"n = {len(variable)} reads from {labels['sample_size']} samples"
        ax2.annotate(summary_stats_str, xy=(0.1, 0.89), ha='left', va='top', fontsize=9, linespacing=1.75, wrap=True)
    else:
        ax1_stop = 5

    # plot histogram
    ax1 = plt.subplot(1,5,(1,ax1_stop))  # plot overall will have 4 equal size columns, hist will span first three
    sns.histplot(variable, color=qualscore_color, kde=True, ax=ax1)  # plot hist and label
    ax1.set_ylabel(labels['y-axis'])
    ax1.set_xlabel(labels['x-axis'])
    ax1.spines['top'].set_visible(False)  # remove top plot border
    ax1.spines['right'].set_visible(False)  # remove bottom plot border

    fig.tight_layout(w_pad=0, rect=(0,(caption_adj/fig_height),1,0.95))
    fig.text(0.5, (caption_adj/fig_height), caption, wrap=True, ha='center', va='top', fontsize=caption_text_size)  # add caption
    fig.suptitle(f"Q{quality_score} Reads ({sequencing_run}): {labels['title']}", fontsize=12, fontweight='demibold')  # add title

    output_path = f"{summary_path}{run_name}{output_tag}.png"
    fig.savefig(output_path)  # save figure
    plt.close()

    if os.path.isfile(output_path):
        print(f"The histogram '{labels['title']}' has been saved: {output_path}")  # return print confirmation that figure is saved
    else:
        print(f"The histogram '{labels['title']}' could not be saved.")

    return None
