"""
Parsers for Assembly related results.
"""
from .generic import get_file_type, parse_json, parse_table
RESULT_TYPE = 'assembly'
ACCEPTED_FILES = ["fna.json", "checkm-results.txt", "transposed_report.tsv"]


def parse(filename: str) -> dict:
    """
    Check input file is an accepted file, then select the appropriate parsing method.

    Args:
        filename (str): input file to be parsed

    Returns:
        dict: parsed results
    """
    filetype = get_file_type(ACCEPTED_FILES, filename)
    if filetype == "fna.json":
        return parse_json(filename)
    elif filetype.endswith("checkm-results.txt") or filetype.endswith("transposed_report.tsv"):
        return parse_table(filename)[0]


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

        if result.endswith('.json'):
            result_name = 'stats'
            filename = f"{path}/{name}/{RESULT_TYPE}/{name}.{result}"
        elif result.endswith("checkm-results.txt"):
            result_name = 'checkm'
            filename = f"{path}/{name}/{RESULT_TYPE}/checkm/{result}"
        elif result.endswith("transposed_report.tsv"):
            result_name = 'quast'
            filename = f"{path}/{name}/{RESULT_TYPE}/quast/{result}"

        parsable_results.append({
            'result_name': result_name,
            'files': [filename],
            'optional': optional,
            'missing': False if os.path.exists(filename) else True
        })

    return parsable_results
