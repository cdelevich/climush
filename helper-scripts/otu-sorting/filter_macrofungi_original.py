########################################################################################################################
## FROM WITHIN: /Volumes/wd_18TB/back-ups/dropbox/dropbox_all_2024-03-31/projects/molecular-work/pacbio_bioinformatics/older-pipelines/modified-pipeline_2310/09_assign-taxonomy/genbank/05_clean_genbank-info.py
########################################################################################################################

approved_macrofungi = pd.read_csv("../../miscell/climush_macrofungi_231009.csv")  # table of certified macrofungi

# lastly, check if the GenBank ID is a macrofungus
final_df['is_macrofungus'] = None

# ugh, rename taxonomy cols to run next script

# get family exceptions
# get a list of families to exclude
fam_step01 = approved_macrofungi[~approved_macrofungi['family_exceptions'].isnull()]['family_exceptions'].values
fam_step02 = [ i.split(", ") for i in fam_step01 ]
family_exceptions = list(np.concatenate(fam_step02))

for row in range(final_df.shape[0]):

    # print(f"Checking for row {row}...")

    test_row = final_df.iloc[row]

    # get name of tax column from genbank df
    tax_cols = [ col for col in final_df.columns if col.startswith("tax") ]


    # get list of taxonomy from vsearch df
    tax_list = [ test_row[col] for col in tax_cols ]

    # print(f"     taxonomy:{tax_list}")

    # can check right away, if family in exceptions list, exclude
    family = test_row['tax_family']

    if test_row['tax_family'] in family_exceptions:
        # print(f"OTU excluded as macrofungus because of family exception.")
        final_df['is_macrofungus'][row] = False
        continue  # continue to next row to check

    for i in range(len(tax_list)):

        # does the current taxonomic level have an id?
        current_lacks_taxon = pd.isnull(tax_list[i])
        if current_lacks_taxon:  # if it doesn't, then it isn't a macrofungus
            final_df['is_macrofungus'][row] = False
            # print(f"     Not a macrofungus based on condition 01, at {tax_cols[i]}.")
            break
        else:  # if it does, is it macro?
            current_is_macro = approved_macrofungi[tax_cols[i]].str.contains(tax_list[i]).any()

        # does the next taxonomic level have an id?
        macro_to_check = (sum(~approved_macrofungi[tax_cols[i+1]].isna()) > 0)
        if macro_to_check:  # are there *any* taxa to check at this taxonomic level? (esp. for species after Elaphomyces)
            next_lacks_taxon = pd.isnull(tax_list[i+1])
            no_further_checks = (sum(pd.isnull(approved_macrofungi[approved_macrofungi[tax_cols[i]] == tax_list[i]][tax_cols[i+1]])) > 0)
            if next_lacks_taxon and current_is_macro:  # if it doesn't, but current is macro, depends on tax level...
                if i < 2:  # if its above order...
                    final_df['is_macrofungus'][row] = False
                    # print(f"     Is not a macrofungus based on condition 02-01A, at {tax_cols[i]}.")
                    break
                else:  # else, if we made it this far, we can guess it is
                    final_df['is_macrofungus'][row] = True
                    # print(f"     Is a macrofungus based on condition 02-01B, at {tax_cols[i]}.")
                    break
            elif no_further_checks:  # if the next taxon lvl for this taxa is 'NA', then it's macro
                final_df['is_macrofungus'][row] = True
                # print(f"     Is a macrofungus based on condition 02-02, at {tax_cols[i]}.")
                break
            else:  # if it does, is it macro?
                next_is_macro = approved_macrofungi[tax_cols[i + 1]].str.contains(tax_list[i + 1]).any()
        else:  # if there aren't further restrictions to check, its macro
            final_df['is_macrofungus'][row] = True
            # print(f"     Is a macrofungus based on condition 02-03, at {tax_cols[i]}.")
            break

        # based on current tax and next, is it a macrofungus?
        if current_is_macro and next_is_macro:  # if they are both in macro table...
            continue  # keep looking
        elif current_is_macro and not next_is_macro:  # if current is, but next is not...
            # print(f"     no_further_checks = {no_further_checks}")
            final_df['is_macrofungus'][row] = False
            # print(f"     Is not a macrofungus based on condition 03-01, at {tax_cols[i]}.")
            break
        else:  # current is not macro...
            final_df['is_macrofungus'][row] = False
            # print(f"     Is not a macrofungus based on condition 03-02, at {tax_cols[i]}.")
            break

    # print(f"     Is is macro: **{tax_table['is_macrofungus'][row]}**")






########################################################################################################################
## FROM WITHIN: /Volumes/wd_18TB/back-ups/dropbox/dropbox_all_2024-03-31/projects/molecular-work/pacbio_bioinformatics/older-pipelines/modified-pipeline_2310/09_assign-taxonomy/vsearch/02_tabulate_vsearch-tax_output.py
########################################################################################################################

approved_macrofungi = pd.read_csv("../../miscell/climush_macrofungi_231009.csv")  # table of certified macrofungi

# get a list of families to exclude
fam_step01 = approved_macrofungi[~approved_macrofungi['family_exceptions'].isnull()]['family_exceptions'].values
fam_step02 = [ i.split(", ") for i in fam_step01 ]
family_exceptions = list(np.concatenate(fam_step02))

tax_table['is_macrofungus'] = None

for row in range(tax_table.shape[0]):

    print(f"Checking for row {row}...")

    test_row = tax_table.iloc[row]

    # get name of tax column from vsearch df
    tax_cols = [ col for col in tax_table.columns if col.startswith("tax") ]

    # get list of taxonomy from vsearch df
    tax_list = [ test_row[col] for col in tax_cols ]

    print(f"     taxonomy:{tax_list}")

    # can check right away, if family in exceptions list, exclude
    family = test_row['tax_family']

    if test_row['tax_family'] in family_exceptions:
        print(f"OTU excluded as macrofungus because of family exception.")
        tax_table['is_macrofungus'][row] = False
        continue  # continue to next row to check

    for i in range(len(tax_list)):

        # does the current taxonomic level have an id?
        current_lacks_taxon = pd.isnull(tax_list[i])
        if current_lacks_taxon:  # if it doesn't, then it isn't a macrofungus
            tax_table['is_macrofungus'][row] = False
            print(f"     Not a macrofungus based on condition 01, at {tax_cols[i]}.")
            break
        else:  # if it does, is it macro?
            current_is_macro = approved_macrofungi[tax_cols[i]].str.contains(tax_list[i]).any()

        # does the next taxonomic level have an id?
        macro_to_check = (sum(~approved_macrofungi[tax_cols[i+1]].isna()) > 0)
        if macro_to_check:  # are there *any* taxa to check at this taxonomic level? (esp. for species after Elaphomyces)
            next_lacks_taxon = pd.isnull(tax_list[i+1])
            no_further_checks = (sum(pd.isnull(approved_macrofungi[approved_macrofungi[tax_cols[i]] == tax_list[i]][tax_cols[i+1]])) > 0)
            if next_lacks_taxon and current_is_macro:  # if it doesn't, but current is macro, depends on tax level...
                if i < 2:  # if its above order...
                    tax_table['is_macrofungus'][row] = False
                    print(f"     Is not a macrofungus based on condition 02-01A, at {tax_cols[i]}.")
                    break
                else:  # else, if we made it this far, we can guess it is
                    tax_table['is_macrofungus'][row] = True
                    print(f"     Is a macrofungus based on condition 02-01B, at {tax_cols[i]}.")
                    break
            elif no_further_checks:  # if the next taxon lvl for this taxa is 'NA', then it's macro
                tax_table['is_macrofungus'][row] = True
                print(f"     Is a macrofungus based on condition 02-02, at {tax_cols[i]}.")
                break
            else:  # if it does, is it macro?
                next_is_macro = approved_macrofungi[tax_cols[i + 1]].str.contains(tax_list[i + 1]).any()
        else:  # if there aren't further restrictions to check, its macro
            tax_table['is_macrofungus'][row] = True
            print(f"     Is a macrofungus based on condition 02-03, at {tax_cols[i]}.")
            break

        # based on current tax and next, is it a macrofungus?
        if current_is_macro and next_is_macro:  # if they are both in macro table...
            continue  # keep looking
        elif current_is_macro and not next_is_macro:  # if current is, but next is not...
            print(f"     no_further_checks = {no_further_checks}")
            tax_table['is_macrofungus'][row] = False
            print(f"     Is not a macrofungus based on condition 03-01, at {tax_cols[i]}.")
            break
        else:  # current is not macro...
            tax_table['is_macrofungus'][row] = False
            print(f"     Is not a macrofungus based on condition 03-02, at {tax_cols[i]}.")
            break

    # print(f"     Is is macro: **{tax_table['is_macrofungus'][row]}**")
