"""
Parsers for Error related results.
"""
from .generic import get_file_type, parse_json, parse_table
ACCEPTED_FILES = ["error.txt"]
ERROR_TYPES = [
    "assembly",
    "different-read-count",
    "genome-size",
    "low-read-count",
    "low-sequence-depth",
    "low-basepair-proportion",
    "paired-end",
]


def parse(filename: str) -> dict:
    """
    Check input file is an accepted file, then select the appropriate parsing method.

    Args:
        filename (str): input file to be parsed

    Returns:
        list: observed error and a brief description
    """
    filetype = get_file_type(ACCEPTED_FILES, filename)
    error = filename

    if error.endswith("genome-size-error.txt"):
        return _format_error(["genome-size-error", "Poor estimate of genome size"])
    elif error.endswith("low-read-count-error.txt"):
        return _format_error(["low-read-count-error", "Low number of reads"])
    elif error.endswith("low-sequence-depth-error.txt"):
        return _format_error(["low-sequence-depth-error", "Low depth of sequencing"])
    elif error.endswith("paired-end-error.txt"):
        return _format_error(["paired-end-error", "Paired-end reads were not in acceptable format"])
    elif error.endswith("different-read-count-error.txt"):
        return _format_error(["different-read-count-error", "Paired-end read count mismatch"])
    elif error.endswith("low-basepair-proportion-error.txt"):
        return _format_error(["low-basepair-proportion-error", "Paired-end basepair counts are out of accesptable proportions"])
    elif error.endswith("assembly-error.txt"):
        return _format_error(["assembly-error", "Assembled size was not withing an acceptable range"])
    return _format_error(["unknown-error", "Unknown Error"])


def _format_error(error: list) -> dict:
    """
    Convert the error type list to a dict.
    Args:
        error (list): a two element list with the error type and description

    Returns:
        dict: explicit names for the list elements
    """
    return {'error_type': error[0], 'description': error[1]}
