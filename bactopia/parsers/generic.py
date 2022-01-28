"""
Shared functions used by parsers.
"""
from typing import Union


def get_file_type(extensions: list, filename: str) -> str:
    """
    Checks if the extension of a given file matches a set of expected extensions.

    Args:
        extensions (list): a set of expected file extensions
        filename (str): a file name to test is extension is expected

    Raises:
        ValueError: given file does not match an expected extension

    Returns:
        str: the matched extension
    """
    for ext in extensions:
        if filename.endswith(ext):
            return ext

    raise ValueError(f"'{filename}' is not an accepted result file. Accepted extensions: {', '.join(extensions)}")


def parse_table(csvfile: str, delimiter: str = '\t', has_header: bool = True) -> Union[list, dict]:
    """
    Parse a delimited file.

    Args:
        csvfile (str): input delimited file to be parsed
        delimiter (str, optional): delimter used to separate column values. Defaults to '\t'.
        has_header (bool, optional): the first line should be treated as a header. Defaults to True.

    Returns:
        Union[list, dict]: A dict is returned if a header is present, otherwise a list is returned
    """
    import csv
    data = []
    with open(csvfile, 'rt') as fh:
        for row in csv.DictReader(fh, delimiter=delimiter) if has_header else csv.reader(fh, delimiter=delimiter):
            data.append(row)
    return data


def parse_json(jsonfile: str) -> Union[list, dict]:
    """
    Parse a JSON file.

    Args:
        jsonfile (str): input JSON file to be read

    Returns:
        Union[list, dict]: the values oarsed from the JSON file
    """
    import json
    with open(jsonfile, 'rt') as fh:
        return json.load(fh)
