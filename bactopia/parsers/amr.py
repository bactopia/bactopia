"""
Parsers for Antimicrobial Resistance related results.
"""
from .generic import get_file_type, parse_table
RESULT_TYPE = 'antimicrobial-resistance'
ACCEPTED_FILES = ["gene-report.txt", "protein-report.txt"]


def parse(filename: str) -> dict:
    """
    Check input file is an accepted file, then select the appropriate parsing method.

    Args:
        filename (str): input file to be parsed

    Returns:
        dict: parsed results
    """
    filetype = get_file_type(ACCEPTED_FILES, filename)
    if filetype.endswith("report.txt"):
        return _parse_amrfinder_report(filename)


def _parse_amrfinder_report(filename: str) -> dict:
    """
    Parse the AMRFinder report file.

    Args:
        filename (str): input file to be parsed

    Returns:
        dict: the parsed AMRFinder+ results
    """
    return parse_table(filename)


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

        if result.endswith('gene-report.txt'):
            result_name = 'gene-report'
            filename = f"{path}/{name}/{RESULT_TYPE}/{name}-{result}"
        elif result.endswith('protein-report.txt'):
            result_name = 'protein-report'
            filename = f"{path}/{name}/{RESULT_TYPE}/{name}-{result}"

        parsable_results.append({
            'result_name': result_name,
            'files': [filename],
            'optional': optional,
            'missing': False if os.path.exists(filename) else True
        })

    return parsable_results
