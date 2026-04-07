from pathlib import Path
from Bio import SeqIO
import re

# path to output directory for taxonomic troubleshooting
tax_troubleshoot_output = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/illumina-2505/scripts/troubleshooting/taxonomic-assignment/output/')

# path to the Amanita-genus file
amanita_genus_in = tax_troubleshoot_output / 'ITS_search01_2025-06-30.fasta'

# path to the output Amanita exitialis file
amanita_exitialis_out = (tax_troubleshoot_output / f'{amanita_genus_in.stem}_amanita-exitialis').with_suffix('.fasta')

# open file and only take records that match the species Amanita exitialis
target_species = 'amanita exitialis'

# create an empty list to add target sequence records to
target_species_seqrecords = []

max_record_display = 20
record_count = 0
with open(amanita_genus_in, 'r') as fasta_in:

    for amanita_record in SeqIO.parse(fasta_in, format='fasta'):


        # get the taxonomy string from this record
        tax_str = amanita_record.description.split(';')[-1]

        # split the tax levels; take only the last in the list
        last_tax_lvl = tax_str.split(',')[-1]

        ## REMOVE AFTER TESTING ##
        record_count += 1
        if max_record_display > record_count:
            print_info = True
        else:
            print_info = False

        if print_info:
            print(f'Checking the taxonomy string: \nfull header: {amanita_record.id}\nlast tax str:{last_tax_lvl}\n')



        # if the last part of the tax levels list starts with an s, it is a species
        if last_tax_lvl.startswith('s'):

            # compare the species string of this record to the target species string (case insensitive)
            matches_target = re.search(f'(?<=s:){target_species}$', last_tax_lvl, re.I)

            # if a match is found...
            if matches_target:

                if print_info:
                    print(f'   has species\nmatches target\n')

                # add this Amanita record to the list of target species sequence records
                target_species_seqrecords.append(amanita_record)

            # if a match is not found...
            else:

                if print_info:
                    print(f'   has species\ndoes not match target\n')

                # move on to the next sequence record in the file
                continue

        # if the last part of the tax levels list does not start with an s...
        else:

            if print_info:
                print(f'   does not have species\n')

            # this reference sequence doesn't have a species assigned, so skip over check for a species match
            continue


# write the matching species records out to a file
if len(target_species_seqrecords) > 0:
    print(f'Writing {len(target_species_seqrecords)} matches to the {target_species} target species to the '
          f'sequence file:\n'
          f'   {amanita_exitialis_out}\n')
    with open(amanita_exitialis_out, 'w') as fasta_out:
        SeqIO.write(target_species_seqrecords, fasta_out, format='fasta')
else:
    print(f'No sequence records matching the target species {target_species} was found in the input '
          f'sequence file:\n'
          f'   {amanita_genus_in}')