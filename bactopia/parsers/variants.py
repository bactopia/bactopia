"""
Parsers for Variant related results.
"""
from .generic import get_file_type
RESULT_TYPE = 'variants'
ACCEPTED_FILES = ["txt"]


def parse(filename: str) -> dict:
    """
    Check input file is an accepted file, then select the appropriate parsing method.

    Args:
        filename (str): input file to be parsed

    Returns:
        dict: parsed results
    """
    filetype = get_file_type(ACCEPTED_FILES, filename)
    if filetype == "txt":
        return _parse_variants(filename)


def _parse_variants(filename: str) -> dict:
    """
    Parse Snippy summary text file.

    Args:
        filename (str): input file to be parsed

    Returns:
        dict: parsed Snippy summary file
    """
    from os.path import basename
    results = {}
    with open(filename, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            if not line.startswith("ReadFiles"):
                key, val = line.split("\t")
                if key == "Reference":
                    results[key] = basename(val).split('-')[-1].split(".")[0]
                else:
                    results[key] = val.lstrip()
    return results


def get_parsable_list(path: str, name: str) -> list:
    """
    Generate a list of parsable files.

    Args:
        path (str): a path to expected Bactopia results
        name (str): the name of sample to test

    Returns:
        list: information about the status of parsable files
    """
    import os
    parsable_results = []
    for variant_source in ['auto', 'user']:
        variant_dir = f"{path}/{name}/{RESULT_TYPE}/{variant_source}"
        if os.path.exists(variant_dir):
            with os.scandir(variant_dir) as dirs:
                for reference in dirs:
                    reference_dir = f"{variant_dir}/{reference}"
                    for result in ACCEPTED_FILES:
                        result_name = None
                        optional = True
                        filename = None

                        if result.endswith('txt'):
                            result_name = 'stats'
                            filename = f"{reference_dir}/{name}.{result}"

                        parsable_results.append({
                            'result_name': result_name,
                            'files': [filename],
                            'optional': optional,
                            'missing': False if os.path.exists(filename) else True
                        })

    return parsable_results
