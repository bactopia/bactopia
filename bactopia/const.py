"""
Constant values used throughout the Bactopia package.
"""

RESULT_TYPES = {
    "annotation": "annotation",
    "antimicrobial-resistance": "amr",
    "ariba": "ariba",
    "assembly": "assembly",
    "blast": "blast",
    "error": "error",
    "generic": "generic",
    "kmers": "kmers",
    # "logs": "logs",
    "mapping": "mapping",
    "minmers": "minmers",
    "mlst": "mlst",
    "quality-control": "qc",
    "variants": "variants",
}

IGNORE_LIST = [
    '.nextflow',
    '.nextflow.log',
    'bactopia',
    'bactopia-info',
    'bactopia-tools',
    'logs',
    'nf-reports',
    'software_versions.yml',
    'software_versions_mqc.yml',
    'staphopia',
    'work'
]
