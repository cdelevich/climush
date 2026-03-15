from pathlib import Path
import json

neon_domains_jsonout = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/climush_py-package/climush/site-codes.json')

# create a dictionary of the NEON ecoregion domain numbers and all associated aliases for each domain
neon_domains = {
    'D01': {
        'ecoregion': 'Northeast',
        'ecoregion abbr': 'NE',
        'site': 'Harvard Forest',
        'site abbr': 'HFMA',
        'state': 'Massachussettes',
        'state abbr': 'MA',
    },
    'D03': {
        'ecoregion': 'Southeast',
        'ecoregion abbr': 'SE',
        'site': 'Ordway-Swisher',
        'site abbr': 'ORD',
        'state': 'Florida',
        'state abbr': 'FL',
    },
    'D05': {
        'ecoregion': 'Great Lakes',
        'ecoregion abbr': 'GRTL',
        'site': 'Cedar Creek',
        'site abbr': 'CDR',
        'state': 'Minnesota',
        'state abbr': 'MN',
    },
    'D06': {
        'ecoregion': 'Prairie Peninsula',
        'ecoregion abbr': 'PRAI',
        'site': 'Konza Prairie',
        'site abbr': 'KON',
        'state': 'Kansas',
        'state abbr': 'KS',
    },
    'D13': {
        'ecoregion': 'Southern Rockies / Colorado Plateau',
        'ecoregion abbr': 'SROC',
        'site': 'Niwot Ridge',
        'site abbr': 'NWT',
        'state': 'Colorado',
        'state abbr': 'CO',
    },
    'D14': {
        'ecoregion': 'Desert Southwest',
        'ecoregion abbr': 'DSW',
        'site': 'Santa Rita Experimental Range',
        'site abbr': 'SRE',
        'state': 'Arizona',
        'state abbr': 'AZ',
    },
    'D16': {
        'ecoregion': 'Pacific Northwest',
        'ecoregion abbr': 'PNW',
        'site': ['Mt. Pisgah', 'H.J. Andrews'],
        'site abbr': ['PIS', 'HJA'],
        'state': 'Oregon',
        'state abbr': 'OR',
    },
    'D19': {
        'ecoregion': 'Taiga',
        'ecoregion abbr': 'TAIG',
        'site': 'Bonanza Creek',
        'site abbr': 'BNZ',
        'state': 'Alaska',
        'state abbr': 'AK',
    },
}

with open(neon_domains_jsonout, 'w') as json_out:
    json.dump(neon_domains, json_out)