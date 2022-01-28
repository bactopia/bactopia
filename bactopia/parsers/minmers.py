"""
Parsers for Minmer related results.
"""
from .generic import get_file_type, parse_table
RESULT_TYPE = 'minmers'
ACCEPTED_FILES = ["refseq-k21.txt", "plsdb-k21.txt", "genbank-k21.txt", "genbank-k31.txt", "genbank-k51.txt"]


def parse(filename: str) -> dict:
    """
    Check input file is an accepted file, then select the appropriate parsing method.

    Args:
        filename (str): input file to be parsed

    Returns:
        dict: parsed results
    """
    filetype = get_file_type(ACCEPTED_FILES, filename)
    if filetype.startswith("genbank"):
        return _parse_sourmash(filename)
    elif filetype == "refseq-k21.txt" or filetype == "plsdb-k21.txt":
        return parse_table(filename)


def _parse_sourmash(filename: str) -> dict:
    """
    Parse Sourmash output.
    
    Example Format:
        overlap     p_query p_match
        ---------   ------- --------
        2.7 Mbp       7.3%   99.3%      Staphylococcus aureus (** 2 equal matches)
        430.0 kbp     1.1%    0.5%      Tetrahymena thermophila
        90.0 kbp      0.2%    0.1%      Paramecium tetraurelia
        80.0 kbp      0.2%    2.7%      Staphylococcus aureus (** 8 equal matches)
        80.0 kbp      0.2%    1.2%      Staphylococcus haemolyticus
        10.0 kbp      0.0%    0.2%      Alcanivorax xenomutans

        74.6% (28.0 Mbp) of hashes have no assignment.

    Args:
        filename (str): input file to be parsed

    Returns:
        dict: the parsed Sourmash results
    """
    import re
    re_sourmash = re.compile(
        r'(?P<overlap>[0-9]+.[0-9]+ [A-Za-z]+)\s+(?P<p_query>[0-9]+.[0-9]+%)\s+(?P<p_match>[0-9]+.[0-9]+%)\s+(?P<match>.*)'
    )
    count = 0
    data = {"matches": [], "no_assignment": ""}
    with open(filename, "rt") as fh:
        parse_row = False
        parse_no_assignment = False
        for line in fh:
            line = line.rstrip()
            if parse_no_assignment:
                data["no_assignment"] = line
            elif parse_row:
                if line:
                    re_match = re_sourmash.match(line)
                    data['matches'].append({
                        'overlap': re_match.group('overlap'),
                        'p_query': re_match.group('p_query'),
                        'p_match': re_match.group('p_match'),
                        'match': re_match.group('match')
                    })
                    count += 1
                else:
                    parse_no_assignment = True
            elif line.startswith("----"):
                parse_row = True
    return data


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

        if result.endswith('.txt'):
            result_name = result.replace('.txt', '')
            filename = f"{path}/{name}/{RESULT_TYPE}/{name}-{result}"

        parsable_results.append({
            'result_name': result_name,
            'files': [filename],
            'optional': True,
            'missing': False if os.path.exists(filename) else True
        })

    return parsable_results
