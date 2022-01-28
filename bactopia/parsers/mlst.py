"""
Parsers for MLST related results.
"""
from .generic import get_file_type, parse_json, parse_table
RESULT_TYPE = 'mlst'
ACCEPTED_FILES = ["blast.json", "mlst_report.tsv"]


def parse(filename: str) -> dict:
    """
    Check input file is an accepted file, then select the appropriate parsing method.

    Args:
        filename (str): input file to be parsed

    Returns:
        dict: parsed results
    """
    filetype = get_file_type(ACCEPTED_FILES, filename)
    if filetype == "blast.json":
        return parse_json(filename)
    elif filetype == "mlst_report.tsv":
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
    mlst_dir = f"{path}/{name}/{RESULT_TYPE}"
    if os.path.exists(mlst_dir):
        with os.scandir(mlst_dir) as dirs:
            for schema_dir in dirs:
                schema = schema_dir.name
                missing = True
                blast = f"{mlst_dir}/{schema}/blast/{name}-blast.json"
                parsable_results.append({
                    'result_name': f"{schema}-blast",
                    'files': [blast],
                    'optional': True,
                    'missing': False if os.path.exists(blast) else True
                })

                ariba = f"{mlst_dir}/{schema}/ariba/mlst_report.tsv"
                parsable_results.append({
                    'result_name': f"{schema}-ariba",
                    'files': [ariba],
                    'optional': True,
                    'missing': False if os.path.exists(ariba) else True
                })

    return parsable_results
