from pathlib import Path

# PIPELINE
try:
    PIPELINE = Path(__file__).parent.resolve()
except NameError:
    PIPELINE = Path('mapping.py').parent.resolve()

PIPE_RENAME = next(PIPELINE.glob('*rename*'))
PIPE_DEMUX = next(PIPELINE.glob('*demu*'))

ROOT = PIPELINE.parent.resolve()

# CONFIGURATION
CONFIG = ROOT / 'config'
BC_MAP = next(CONFIG.glob('*barcode*'))
REF_DB = CONFIG / 'reference-db'
CHIM_DB = REF_DB / 'chimera'
TAX_DB = REF_DB / 'taxonomy'

# PYTHON PACKAGE
PACKAGE = ROOT / 'climush_py-package'

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
DEMUX = PIPELINE_OUT / 'demultiplexed'
PREFILTER = PIPELINE_OUT / 'prefiltered'
PRIMER_TRIM = PIPELINE_OUT / 'primers-trimmed'
DEREP01 = PIPELINE_OUT / 'derep-full-length'
QFILT = PIPELINE_OUT / 'quality-filtered'
MERGE = PIPELINE_OUT / 'merged'
ITSX = PIPELINE_OUT / 'separate-subregions'
CHIMERA = PIPELINE_OUT / 'chimera-checked'
DEREP02 = PIPELINE_OUT / 'derep-subregions'
CLUSTER = PIPELINE_OUT / 'otus-clustered'
TAX = PIPELINE_OUT / 'taxonomy-assigned'
OTU_TAB = PIPELINE_OUT / 'otu-table'

filepath_map = {'root': ROOT,
                'config': {'main':CONFIG,
                           'bc_mapping':BC_MAP,
                           'reference-db': {'main':REF_DB,
                                            'chimera':CHIM_DB,
                                            'taxonomy':TAX_DB}},
                'pipeline': {'main':PIPELINE,
                             'rename':PIPE_RENAME,
                             'demux':PIPE_DEMUX},
                'package': {'main': PACKAGE},
                'sequences': {'main': SEQUENCES,
                              'pacbio':SORT_PACBIO,
                              'illumina':SORT_ILLUMINA,
                              'sanger':SORT_SANGER,
                              'rename':SORT_RENAME,
                              'demux':SORT_DEMUX,
                              'unclear':SORT_UNCLEAR},
                'pipeline-output': {'main': PIPELINE_OUT,
                                    'demultiplexed': DEMUX,
                                    'prefiltered': PREFILTER,
                                    'primers-trimmed': PRIMER_TRIM,
                                    'derep-full-length': DEREP01,
                                    'quality-filtered': QFILT,
                                    'merged': MERGE,
                                    'separate-subregions': ITSX,
                                    'chimera-checked': CHIMERA,
                                    'derep-subregions': DEREP02,
                                    'otus-clustered': CLUSTER,
                                    'taxonomy': TAX,
                                    'otu-table': OTU_TAB}}