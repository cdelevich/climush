from pathlib import Path
import sys
from climush.utilities import mkdir_exist_ok

# path to the pipeline scripts
try:
    PIPELINE = Path(__file__).parent.resolve()
except NameError:
    PIPELINE = Path('mapping.py').parent.resolve()
# WHY DO I HAVE THIS?
PIPE_DEMUX = next(PIPELINE.glob('*demu*'))

# path to the main pipeline directory
ROOT = PIPELINE.parent.resolve()

# path to helper scripts
HELPERS = ROOT.parent.resolve() / 'helper-scripts'

# THIS ONE WAS RENAMED FROM PIPE_RENAME
HELPERS_RENAME = next(HELPERS.glob('*nam*'))  # dir containing renaming scripts/files
HELPERS_SORT = next(HELPERS.glob('*sort*'))

# CONFIGURATION
CONFIG = ROOT / 'config'
BC_MAP = next(CONFIG.glob('*barcode*'))
REF_DB = CONFIG / 'reference-db'
CHIM_DB = REF_DB / 'chimera'
TAX_DB = REF_DB / 'taxonomy'

# SEQUENCES
SEQUENCES = ROOT / 'sequences'
SORT_PACBIO = SEQUENCES / 'pacbio'
SORT_ILLUMINA = SEQUENCES / 'illumina'
SORT_SANGER = SEQUENCES / 'sanger'

SORT_RENAME = SEQUENCES / 'needs_rename'
SORT_DEMUX = SEQUENCES / 'needs_demux'
SORT_UNCLEAR = SEQUENCES / 'needs_clarification'

# PIPELINE-OUTPUT
PIPELINE_OUT = ROOT / 'pipeline-output'
mkdir_exist_ok(PIPELINE_OUT)  # make the main output dir; will execute if this script is run in any pipeline script

# don't make these directories unless they produce output (so mkdir in their pipeline script, otherwise don't create)
DEMUX = PIPELINE_OUT / 'demultiplexed'
PREFILTER = PIPELINE_OUT / 'prefiltered'
PREFILTER_01 = 'prefilt01_no-phix'
PREFILTER_02 = 'prefilt02_no-ambig'
PRIMER_TRIM = PIPELINE_OUT / 'primers-trimmed'
DEREP01 = PIPELINE_OUT / 'derep-full-length'
QFILT = PIPELINE_OUT / 'quality-filtered'
MERGE = PIPELINE_OUT / 'merged'
ITSX = PIPELINE_OUT / 'separated-subregions'
CHIMERA = PIPELINE_OUT / 'chimera-checked'
DEREP02 = PIPELINE_OUT / 'derep-subregions'
CLUSTER = PIPELINE_OUT / 'otus-clustered'
TAX = PIPELINE_OUT / 'taxonomy-assigned'
OTU_TAB = PIPELINE_OUT / 'otu-table'
SUMMARY = PIPELINE_OUT / 'summary'

filepath_map = {'root': ROOT,
                'config': {'main':CONFIG,
                           'bc_mapping':BC_MAP,
                           'reference-db': {'main':REF_DB,
                                            'chimera':CHIM_DB,
                                            'taxonomy':TAX_DB}},
                'pipeline': {'main':PIPELINE,
                             'demux':PIPE_DEMUX},
                'sequences': {'main': SEQUENCES,
                              'pacbio':SORT_PACBIO,
                              'illumina':SORT_ILLUMINA,
                              'sanger':SORT_SANGER,
                              'rename':SORT_RENAME,
                              'demux':SORT_DEMUX,
                              'unclear':SORT_UNCLEAR},
                'pipeline-output': {'main': PIPELINE_OUT,
                                    'demultiplexed': DEMUX,
                                    'prefiltered': {'main': PREFILTER,
                                                    PREFILTER_01: PREFILTER / PREFILTER_01,
                                                    PREFILTER_02: PREFILTER / PREFILTER_02},
                                    'primers-trimmed': PRIMER_TRIM,
                                    'quality-filtered': QFILT,
                                    'derep-full-length': DEREP01,
                                    'merged': MERGE,
                                    'separated-subregions': ITSX,
                                    'chimera-checked': CHIMERA,
                                    'derep-subregions': DEREP02,
                                    'otus-clustered': CLUSTER,
                                    'taxonomy': TAX,
                                    'otu-table': OTU_TAB,
                                    'summary': SUMMARY}
                }