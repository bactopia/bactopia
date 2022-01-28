"""
Bactopia's parser entry-point.

Example: bactopia.parse(result_type, filename)
"""
import errno
import os
from collections import OrderedDict
from typing import Union
from . import parsers
from .const import RESULT_TYPES, IGNORE_LIST


def parse(result_type: str, *files: str) -> Union[list, dict]:
    """
    Use the result type to automatically select the appropriate parsing method for an input.

    Args:
        result_type (str): the type of results (e.g. assembly, mlst, qc, etc...)
        *files (str): one or more input files to be parsed

    Raises:
        FileNotFoundError: the input file could not be found
        ValueError: the result type is not an accepted type

    Returns:
        Union[list, dict]: The results parsed for a given input.
    """
    if result_type in RESULT_TYPES:
        for f in files:
            if not os.path.exists(f):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f)
        
        return getattr(parsers, result_type).parse(*files)
    else:
        raise ValueError(f"'{result_type}' is not an accepted result type. Accepted types: {', '.join(RESULT_TYPES)}")


def parse_genome_size(gs_file: str) -> int:
    """
    Parse genome size from input file

    Args:
        gs_file (str): File containing the genome size of the sample

    Returns:
        int: genome size
    """
    with open(gs_file, 'rt') as gs_fh:
        return int(gs_fh.readline().rstrip())


def _is_bactopia_dir(path: str, name: str) -> list:
    """
    Check if a directory contains Bactopia output and any errors.

    Args:
        path (str): a path to expected Bactopia results
        name (str): the name of sample to test

    Returns:
        list: 0 (bool): path looks like Bactopia, 1 (list): any errors found
    """
    from .parsers.error import ERROR_TYPES
    errors = []
    is_bactopia = os.path.exists(f"{path}/{name}/{name}-genome-size.txt")

    for error_type in ERROR_TYPES:
        filename = f"{path}/{name}/{name}-{error_type}-error.txt"
        if os.path.exists(filename):
            is_bactopia = True
            errors.append(parsers.error.parse(filename))

    return [is_bactopia, errors]


def get_bactopia_files(path: str, name: str) -> dict:
    """
    Build a list of all parsable Bactopia files.

    Args:
        path (str): a path to expected Bactopia results
        name (str): the name of sample to test

    Returns:
        dict: path and info on all parsable Bactopia files
    """
    path = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    is_bactopia, errors = _is_bactopia_dir(path, name)
    bactopia_files = OrderedDict()
    bactopia_files['has_errors'] = True if errors else False
    bactopia_files['errors'] = errors
    bactopia_files['ignored'] = False
    bactopia_files['message'] = ""
    bactopia_files['files'] = OrderedDict()

    if is_bactopia:
        if not errors:
            bactopia_files['genome_size'] = parse_genome_size(f"{path}/{name}/{name}-genome-size.txt")
            for result_type in RESULT_TYPES:
                result_key = result_type
                if result_type == "amr":
                    result_key = "antimicrobial-resistance"
                elif result_type == "qc":
                    result_key = "quality-control"

                if result_type not in ['error', 'generic', 'kmers']:
                    bactopia_files['files'][result_key] = getattr(parsers, result_type).get_parsable_list(path, name)
    else:
        bactopia_files['ignored'] = True
        if name not in IGNORE_LIST:
            bactopia_files['message'] = f"'{path}/{name}' is does not look like Bactopia directory."
        else:
            bactopia_files['message'] = f"'{path}/{name}' is on the Bactopia ignore list."

    return bactopia_files


def parse_bactopia_files(path: str, name: str) -> dict:
    """
    Parse all results associated with an input sample.

    Args:
        path (str): a path to expected Bactopia results
        name (str): the name of sample to test

    Returns:
        dict: The parsed set of results associated with the input sample
    """
    from bactopia.parsers.qc import is_paired
    bactopia_files = get_bactopia_files(path, name)
    bactopia_results = OrderedDict((
        ('sample', name),
        ('genome_size', None),
        ('is_paired', None),
        ('has_errors', bactopia_files['has_errors']),
        ('errors', bactopia_files['errors']),
        ('has_missing', False),
        ('missing', []),
        ('ignored', bactopia_files['ignored']),
        ('message', bactopia_files['message']),
        ('results', OrderedDict())
    ))

    if not bactopia_results['has_errors'] and not bactopia_results['ignored']:
        bactopia_results['genome_size'] = bactopia_files['genome_size']
        for result_type, results in bactopia_files['files'].items():
            bactopia_results['is_paired'] = is_paired(path, name)
            bactopia_results['results'][result_type] = OrderedDict()
            result_key = result_type
            if result_type == "-_resistance":
                result_key = "amr"
            elif result_type == "quality-control":
                result_key = "qc"
            
            for result in results:
                if result['missing']:
                    if not result['optional']:
                        bactopia_results['has_missing'] = True
                        bactopia_results['missing'].append([result_type, result["files"]])
                    else:
                        bactopia_results['results'][result_type][result['result_name']] = {}
                else:
                    bactopia_results['results'][result_type][result['result_name']] = parse(result_key, *result['files'])
            
    return bactopia_results


def parse_bactopia_directory(path: str) -> list:
    """
    Scan a Bactopia directory and return parsed results.

    Args:
        path (str):  a path to expected Bactopia results

    Returns:
        list: Parsed results for all samples in a Bactopia directory
    """
    results = []

    with os.scandir(path) as dirs:
        for directory in dirs:
            if directory.name not in IGNORE_LIST:
                sample = directory.name
                results.append(parse_bactopia_files(path, sample))

    return results
