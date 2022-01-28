"""
Parsers for Annotation related results.
"""
from .generic import get_file_type
RESULT_TYPE = 'annotation'
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
        return _parse_annotation(filename)


def _parse_annotation(filename: str) -> dict:
    """
    Parse Prokka summary text file.

    Args:
        filename (str): input file to be parsed

    Returns:
        dict: the parsed Prokka summary
    """
    results = {}
    with open(filename, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            key, val = line.split(":")
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
    for result in ACCEPTED_FILES:
        result_name = None
        optional = False
        filename = None

        if result.endswith('txt'):
            result_name = 'stats'
            filename = f"{path}/{name}/{RESULT_TYPE}/{name}.{result}"

        parsable_results.append({
            'result_name': result_name,
            'files': [filename],
            'optional': optional,
            'missing': False if os.path.exists(filename) else True
        })

    return parsable_results
