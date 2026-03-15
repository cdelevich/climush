from pathlib import Path
import json
from typing import Literal, get_args

# import NEON domain codes and site names from climush package data
climush_data_main = Path(__file__ / 'data')
neon_domain_codes_pathin = climush_data_main / 'site-codes.json'
with open(neon_domain_codes_pathin, 'r') as codes_in:
    neon_domain_codes = json.load(codes_in)

class SampleNameError(Exception):
    """
    Exception raised when sample names are not formatted correctly.
    """

    def __init__(self, message):

        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message}'


# class Collection:
#
#     def __init__(self):
#
#         self._site

# class Sample(Collection):

class Sample:

    # argument options for on_name_error parameter
    _NAME_ERR_OPTS = Literal['stop', 'rename', 'ignore']

    neon_domain_codes = neon_domain_codes

    label_delim = '_'

    sample_labels = {
        0: 'seq_platform',
        1: 'dna_source',
        2: 'collection_date',
        3: 'domain_code',
        4: 'treatment',
        5: 'subplot',
    }

    def __init__(self, sample_name, on_name_error: _NAME_ERR_OPTS='rename'):

        valid_input = get_args(self._NAME_ERR_OPTS)
        assert on_name_error in valid_input, f'{on_name_error} not in {valid_input}'

        self._on_name_error = on_name_error
        self.sample_name = self.verify_sample_name(sample_name)

        self._fieldsite   = ''
        self._dnasource   = ''
        self._collectdate = ''

    @classmethod
    def from_filename(cls, filepath: Path):
        """
        Initialize a Sample object from a file name.

        :param filepath:
        :return: a class instance of Sample
        """

        file_suffixes = '.'.join(filepath.suffixes)
        sample_name = filepath.name.replace(file_suffixes, '')

        return Sample(sample_name)

    def fix_sample_name(self, sample_name):
        """
        Reformat a sample name that is not correctly formatted as a CliMush
        sample name, according to the standardized CliMush sample naming
        convention.
        :param sample_name: incorrectly formatted sample name
        :return: a correctly formatted sample name
        """

    def verify_sample_name(self, sample_name):

        sample_name_labels = sample_name.split(self.label_delim)

        if len(sample_name_labels) == len(self.sample_labels):
            ...

        else:
            if self._on_name_error == 'rename':
                return self.fix_sample_name(sample_name)
            elif self.on_name_error == 'ignore':
                return sample_name
            else:
                raise SampleNameError('Sample name not formatted correctly.')

    @property
    def fieldsite(self):

        sample_name_labels = self.sample_name.split(self.label_delim)
        domaincode_label_pos = [pos for pos, label in self.sample_labels.items() if label == 'domain_code'][0]
        domaincode = sample_name_labels[domaincode_label_pos]

        self.neon_domain_codes[domaincode]['site']

        return self._fieldsite

    @fieldsite.setter
    def fieldsite(self, fieldsite):
        self._fieldsite = fieldsite

    @property
    def dnasource(self):
        return self._dnasource

    @dnasource.setter
    def dnasource(self, dnasource):
        self._dnasource = dnasource

    @property
    def collectdate(self):
        return self._collectdate

    @collectdate.setter
    def collectdate(self, collectdate):
        self._collectdate = collectdate

class envSample(Sample):

    def __init__(self):

        self._subplot = ''

class Endophyte(envSample):

    def __init__(self):

        self._tissuetype   = ''
        self._plantspecies = ''

    @property
    def tissuetype(self):
        return self._tissuetype

    @tissuetype.setter
    def tissuetype(self, tissuetype):
        return self._tissuetype

    @property
    def plantspecies(self):
        return self._plantspecies

    @plantspecies.setter
    def plantspecies(self, plantspecies):
        self._plantspecies = plantspecies

class Sporocarp(Sample):

    def __init__(self):

        self._fieldname = ''

    @property
    def fieldname(self):
        return self._fieldname

    @fieldname.setter
    def fieldname(self, fieldname):
        self._fieldname = fieldname