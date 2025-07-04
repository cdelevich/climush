from multiprocessing.managers import Value
from pathlib import Path
import pandas as pd
import numpy as np
import re, pickle

# path to the IDP table of degeneracy codes and their associated bases
degencode_df_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/primer-detection/degeneracy-codes/idp_degeneracy-codes_2025-07-04.csv')

# path to the climush python package, to which the pickle of the degeneracy codes will be written
climush_package_path = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/climush_py-package/climush')

# path to the pickle output of the conversion dictionary
degencode_pickle_path = climush_package_path / 'degeneracy_codes.pickle'

# read in the degeneracy code table as a pandas df
degencode_df = pd.read_csv(degencode_df_path)

# get the name for the degeneracy code column and the mixed base column
degeneracy_code_col = 'Mixed-base code for DNA'
mixed_based_col = 'Mixed bases required'

# create a dictionary called degeneracy_codes, first including the degeneracy codes as the keys
degeneracy_codes = { code:[] for code in degencode_df[degeneracy_code_col] }

# split the values in the mixed base column first by the comma separator
degencode_df['mixed_bases_split'] = degencode_df[mixed_based_col].apply(lambda x: x.split(','))

# create an empty list to add any base reading errors to
base_errors = {
    'multi-letter':[],
    'unrecognized punctuation':[],
}

# create a quick function that can be used on a single base or a list of bases (in case of split on /)
def check_base_updated_dict(bases, dict_key, error_dict=base_errors,dict_name=degeneracy_codes):

    # create a list if input is single string so I can loop either way
    if isinstance(bases, str):
        bases = [bases]
    elif isinstance(bases, list):
        pass
    else:
        err_msg = f'Input type {type(bases)} is not valid; must be type str or list.'
        raise ValueError(err_msg)

    for b in bases:

        # confirm that this base has only one letter
        if len(b) == 1:

            # add this base to the list of bases corresponding to this degen code
            dict_name[dict_key].append(b)

        # if this base isn't one letter...
        else:

            # add it to the list of errors
            error_dict['multi-letter'].append(b)

    return None

# reformat and add bases to value list in degeneracy_codes dictionary
for i in np.arange(degencode_df.shape[0]):

    # get the list of values in the split mixed bases column, first stripping any whitespace
    base_list = [ base.strip() for base in degencode_df['mixed_bases_split'].iloc[i] ]

    # get the name of the degeneracy code for this list of bases
    degen_code = degencode_df[degeneracy_code_col].iloc[i]

    # iterate through this list and check the formatting of the base
    for b in base_list:

        # search for any punctuation in the base, which will then need to be reformatted
        punctuation_list = re.findall(r'\W', b)

        # if the base has any punctuation included, reformat
        if len(punctuation_list) > 0:

            # go through each punctuation detected...
            for punct in punctuation_list:

                # if a forward slash detected
                if punct == '/':

                    # split into the two component bases
                    comp_bases = b.split('/')

                    # add each base to the dictionary
                    check_base_updated_dict(
                        bases=comp_bases,
                        dict_key=degen_code,
                    )

                # if a parenthesis (either direction) is in the punctuations detected...
                elif (punct == '(') or (punct == ')'):

                    # search for only the valid letters in this base string, returning as separate values
                    comp_bases = re.findall(r'[A-Z]', b)

                    # add each base to the dictionary
                    check_base_updated_dict(
                        bases=comp_bases,
                        dict_key=degen_code,
                    )

                # if other punctuation, I don't know what to do so add to errors
                else:
                    base_errors['unrecognized punctuation'].append(f'{degen_code} = {b}')

        # if no punctuation was found...
        else:

            # add the base to the dictionary
            check_base_updated_dict(
                bases=b,
                dict_key=degen_code,
            )


# check if any errors were detected
for error_categ, errors in base_errors.items():
    if len(errors) > 0:
        errrors_fmt = '\n   '.join(errors)
        err_msg = (f'{len(errors)} {error_categ} errors were encountered when looking through the mixed base lists:\n'
                   f'   {errrors_fmt}')
        raise KeyboardInterrupt(err_msg)
    else:
        pass

# export the dictionary as a pickle for import to constants.py
with open(degencode_pickle_path, 'wb') as pickle_out:
    pickle.dump(degeneracy_codes, pickle_out)

# try to read in the degeneracy codes constant from constants.py
from climush.constants import DEGENERACY_CODES

msg_base = (f'The addition of DEGENERACY_CODES as a constant in the climush package '
            f'module constants.py was ')

if isinstance(DEGENERACY_CODES, dict):
    success_msg = msg_base + 'successful.\n'
    print(success_msg)
else:
    err_msg = msg_base + 'unsuccessful.\n'
    raise ImportError(err_msg)