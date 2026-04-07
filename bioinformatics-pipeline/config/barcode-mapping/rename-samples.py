from pathlib import Path
import shutil
import re
import json

import pandas as pd

pd.set_option('display.max_columns', 15)

SAMPLEID_DELIM = '_'

# path to the main climush project directory
CLIMUSH_MAIN = Path('/Users/carolyndelevich/main/github_repos/climush')

# path to the bioinformatics pipeline subdirectory of the climush project
climush_bioinfo_main = CLIMUSH_MAIN / 'bioinformatics-pipeline'

# path to the climush python package; NEON site domain codes here
climush_bioinfo_pypackage = climush_bioinfo_main / 'climush_py-package'
climush_bioinfo_pypackage_src = climush_bioinfo_pypackage / 'climush'

# configuration directory for bioinformatics
climush_bioinfo_config_main = climush_bioinfo_main / 'config'

# directory within the bioinformatics configuration directory containing the barcode mapping
#   information for multiplexed samples
climush_bioinfo_config_bcmap = climush_bioinfo_config_main / 'barcode-mapping'

# subdirectory in barcode mapping into which the original versions of the barcode mapping
#   files, with the incorrect sample ID formats, will be copied
climush_bioinfo_config_bcmap_originals = climush_bioinfo_config_bcmap / 'original-mappings'
climush_bioinfo_config_bcmap_originals.mkdir(exist_ok=True)


# input file -- barcode pairs and sample IDs for the pacbio soil-litter samples
psl_bcmap_pathin = climush_bioinfo_config_bcmap / 'pacbio_soil-litter_2023-10_barcode-mapping.xlsx'

# reference file -- NEON domain codes and associated site aliases
neon_domains_pathin = climush_bioinfo_pypackage_src / 'site-codes.json'


## MAKE COPY ##

# before importing the barcode mapping dataframe to update, make a copy of it in the original mapping
#   subdirectory of the barcode mapping folder

## SOLUTION REQUIRES PYTHON / PATHLIB V. 3.14.0
# psl_bcmap_pathin.copy_into(
#     target = climush_bioinfo_config_bcmap_originals,
#     preserve_metadata = True,
# )
psl_bcmap_copy_pathout = climush_bioinfo_config_bcmap_originals / psl_bcmap_pathin.name
shutil.copy2(
    src = psl_bcmap_pathin,
    dst = psl_bcmap_copy_pathout,
)



## IMPORT ##

# import the NEON domains site alias dictionary
with open(neon_domains_pathin, 'r') as json_in:
    neon_domains = json.load(json_in)

# import each tab of the dataframe; returns dictionary where keys are the tab names and values are the dataframes in the tab
psl_bcmap_input = pd.read_excel(
    psl_bcmap_pathin,
    sheet_name = None,  # import all tabs of Excel sheet as dictionary
)

# check that the barcode map has the expected number of tabs (one per pool, two total pools)
tab_count_exp = 2
tab_count_obs = len(psl_bcmap_input)
if tab_count_obs == tab_count_exp:
    pass
else:
    input_tab_names = '\n\t'.join(list(psl_bcmap_input.keys()))
    err_msg = (f'Expected {tab_count_exp} tabs in the input Excel file (one for each pool), '
               f'got {tab_count_obs}.:\n'
               f'\t{input_tab_names}\n')
    raise ValueError(err_msg)


# functions to relabel the label components of the input sample IDs to match the output format
def relabel_dnasource(label_torename: str)->str:

    if label_torename == 'S':
        return 'soil'
    elif label_torename == 'L':
        return 'litter'
    else:
        err_msg = f'Unrecognized {label_type}: {label_torename}.\n'
        raise ValueError(err_msg)

def relabel_site(label_torename: str)->str:

    for neon_domain, site_aliases in neon_domains.items():

        for alias_type, alias_value in site_aliases.items():

            if isinstance(alias_value, str):
                if label_torename == alias_value:
                    return neon_domain
                else:
                    continue

            # the D16 Oregon domain has two site aliases (PIS, HJA) so search list
            else:
                if label_torename in alias_value:
                    return neon_domain
                else:
                    continue

    err_msg = (f'The site label {label_torename} is not a recognized alias for any of the '
               f'NEON sites.\n')
    raise ValueError(err_msg)

def relabel_treatment(label_torename: str, treatment_delim='-')->str:

    if treatment_delim in label_torename:
        treatment_sublabels = label_torename.split(treatment_delim)
        sublabel_count_exp = 2
        sublabel_count_obs = len(treatment_sublabels)
        if sublabel_count_obs == sublabel_count_exp:
            trtmnt01 = treatment_sublabels[0].lower()
            trtmnt02 = treatment_sublabels[1].lower()
        else:
            trtmnt_fmt = '\n\t\t'.join(treatment_sublabels)
            err_msg = (f'The number of sublabels in the treatment label {label_torename} after '
                       f'separating on the \'{treatment_delim}\' delimiter is greater than '
                       f'the expected number of sublabels:\n'
                       f'\texpected: {sublabel_count_exp}\n'
                       f'\tobserved: {sublabel_count_obs}\n'
                       f'\t\t{trtmnt_fmt}\n')
            raise ValueError(err_msg)
    else:
        err_msg = (f'The treatment label delimiter, {treatment_delim}, was not '
                   f'detected in the input string: {label_torename}\n')
        raise ValueError(err_msg)

    # what characters tend to be describing the habitat (oak, conifer, grassland) of the sample?
    habitat_chars = ['o', 'c', 'g']

    burnhist_sublabel = None
    habitat_sublabel = None

    # check for any of the habitat characters in the treatment substrings
    habitat_multiple_detects = []
    for h in habitat_chars:
        if (h in trtmnt01) and (h in trtmnt02):
            habitat_multiple_detects.append(h)
        elif (h not in trtmnt01) and (h not in trtmnt02):
            continue
        else:
            habitat_sublabel = h.upper()
            if h in trtmnt01:
                burnhist_sublabel = trtmnt02[0].upper()
            else:
                burnhist_sublabel = trtmnt01[0].upper()

    if (burnhist_sublabel is None) or (habitat_sublabel is None):
        if burnhist_sublabel is None:
            err_msg = f'Burn history sublabel not found in: {label_torename}'
        elif habitat_sublabel is None:
            err_msg = f'Habitat sublabel not found in: {label_torename}'
        else:
            err_msg = f'Neither the burn history sublabel nor the habitat sublabel found in: {label_torename}\n'
        raise ValueError(err_msg)

    else:
        return burnhist_sublabel + habitat_sublabel

def relabel_subplot(label_torename: str)->str:

    # search for subplot number, accounting for possible A/B/C etc. collections from a subplot
    number_found = re.search(r'\d{1,2}[A-Z]?', label_torename)
    if number_found:
        subplot_num = number_found.group(0).zfill(2)
        return subplot_num
    else:
        err_msg = f'No number was found in the {label_type} string: {label_torename}\n'
        raise ValueError(err_msg)


# how many labels does the updated sample ID have?
updated_label_count = 6

# create a sequence run prefix for these samples; will be the prefix for non-sample controls,
#   like NTC1 / NTC2 samples and MockCommunity1
sequence_run_prefix = 'pacbio_soil-litter_2023-10'

# create a dictionary for relabeling, where the key is the index position of the label in the
#   updated sample name and the value is either (1) a string, if standardized for all samples, or
#   (2) a function that will take the original label and recode it to the correct format
relabel_methods = {
    0: 'pacbio',
    1: relabel_dnasource,
    2: '2023-10',
    3: relabel_site,
    4: relabel_treatment,
    5: relabel_subplot,
}

# create a dictionary that remaps the position of the original input labels and the corrected
#   output label; if an output label isn't included in the original input sample ID, or the
#   input sample ID has a label that is not included in the renamed output sample ID, then omit
#   it entirely from this dictionary
sample_position_map = {
    0: 3,  # site
    1: 4,  # treatment
    3: 5,  # subplot
    4: 1,  # eDNA source
}


# iterate through each tab of the input dataframe
for tab_name, tab_bcmap in psl_bcmap_input.items():

    # find name of sample ID column
    sample_column_matches = [col for col in tab_bcmap.columns if re.search('samp', col, re.I)]
    if len(sample_column_matches) == 0:
        err_msg = ('Could not locate the name of the column containing the sample ID by '
                   'search for the substring \'samp\'\n')
        raise KeyError(err_msg)
    else:
        if len(sample_column_matches) == 1:
            sample_colname = sample_column_matches[0]
        else:
            sample_column_rematches = [col for col in tab_bcmap.columns if re.search('sample', col, re.I)]
            if len(sample_column_rematches) == 1:
                sample_colname = sample_column_rematches[0]
            else:
                err_msg = (f'Could not locate a distinct sample ID column using the following substrings:\n'
                           f'\t\'samp\'   = {len(sample_column_matches)} matches\n'
                           f'\t\'sample\' = {len(sample_column_rematches)} matches\n')
                raise KeyError(err_msg)

    # create a dictionary of original and updated sample IDs (updated empty values for now)
    sample_rename_map = {sampid_original: '' for sampid_original in tab_bcmap[sample_colname]}


    # iterate through each sample ID and create a correctly formatted sample ID based on the original
    for sampid_original in sample_rename_map:

        # split the label into its components by the delimiter
        sampid_original_labels = sampid_original.split(SAMPLEID_DELIM)

        # if the sample ID is composed of a single label, it is like a control sample
        if len(sampid_original_labels) == 1:

            # check if control sample is numbered
            ctrl_num_found = re.search(r'\d+', sampid_original)
            if ctrl_num_found:
                ctrl_num = ctrl_num_found.group(0).zfill(2)
            else:
                ctrl_num = None

            # mock community
            if re.search(r'mock', sampid_original, re.I):
                ctrl_type = 'mockcomm'

            # negative control
            elif re.search(r'ntc', sampid_original, re.I):
                ctrl_type = 'negctrl'

            # unrecognized
            else:
                err_msg = f'Unrecognized control sample: {sampid_original}\n'
                raise ValueError(err_msg)

            # assemble components of the control sample ID
            if ctrl_num is None:
                sampid_updated = SAMPLEID_DELIM.join([sequence_run_prefix, ctrl_type])
            else:
                ctrl_type_numbered = ctrl_type + '-' + ctrl_num
                sampid_updated = SAMPLEID_DELIM.join([sequence_run_prefix, ctrl_type_numbered])

        # if the sample ID has at least two labels, it is likely NOT a control sample
        else:

            # create a dictionary in which the recoded labels will be stored, by index they appear in output sample ID
            sampid_updated_labels = {i:'' for i in range(updated_label_count)}

            # iterate through the position and label value of each label
            for label_original_pos, label_original in enumerate(sampid_original_labels):

                # check if this label is used in the renamed sample ID
                if label_original_pos in sample_position_map:

                    # if the label should be recoded, get the position that the updated
                    #   label will appear in in the updated sample ID
                    label_updated_pos = sample_position_map[label_original_pos]

                    # get the method to use to relabel this original label
                    relabel_method = relabel_methods[label_updated_pos]

                    # if the 'method' is a string, then it is a constant; simply add to updated labels
                    if isinstance(relabel_method, str):
                        sampid_updated_labels.update({label_updated_pos: relabel_method})

                    # if the method is in fact a callable function, run the function on the original
                    #   label then add the output to the updated labels dictionary
                    else:
                        sampid_updated_labels.update({label_updated_pos: relabel_method(label_original)})

        # go through the updated labels and fill in any labels for the output relabel sample ID that
        #   aren't in the original input sample ID
        for i in range(updated_label_count):
            if sampid_updated_labels[i] == '':
                sampid_updated_labels.update({i:relabel_methods[i]})
            else:
                continue

        # after recoding all original labels that are to be included in the output relabeled sample ID,
        #   join together the labels into a single updated sample ID
        sampid_updated = SAMPLEID_DELIM.join(sampid_updated_labels.values())

        # add the updated sample ID to the renaming map
        sample_rename_map.update({sampid_original: sampid_updated})


    # use the sample rename map to replace the original input sample IDs in the input dataframe
    tab_bcmap.replace(
        {sample_colname: sample_rename_map},
        inplace=True,
    )

# export the dataframe tabs with the updated sample IDs to the location of the original input Excel file
with pd.ExcelWriter(
    psl_bcmap_pathin,
    mode = 'a',
    engine = 'openpyxl',
    if_sheet_exists = 'replace',
) as writer:
    for tab_name, tab_bcmap in psl_bcmap_input.items():
        tab_bcmap.to_excel(
            writer,
            sheet_name=tab_name
        )