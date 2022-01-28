"""
Parsers for QC (FASTQ) related results.
"""
import os
from .generic import get_file_type, parse_json
RESULT_TYPE = 'quality-control'
ACCEPTED_FILES = ["final.json", "original.json"]


def parse(r1: str, r2: str = None) -> dict:
    """
    Check input file is an accepted file, then select the appropriate parsing method.

    Args:
        r1 (str): input file associated with R1 or SE FASTQ
        r2 (str, optional): input file associated with R2 FASTQ. Defaults to None.

    Raises:
        ValueError: summary results to not have a matching origin (e.g. original vs final FASTQ)

    Returns:
        dict: parsed results
    """
    filetype = get_file_type(ACCEPTED_FILES, r1)
    filetype2 = filetype
    if r2:
        filetype2 = get_file_type(ACCEPTED_FILES, r2)

    if r1.endswith(".json"):
        if r2 and filetype != filetype2:
            raise ValueError(f"Original and Final QC files were mixed. R1: {filetype}, R2: {filetype2}")
        return _merge_qc_stats(parse_json(r1), parse_json(r2)) if r2 else parse_json(r1)


def _merge_qc_stats(r1: dict, r2: dict) -> dict:
    """
    Merge appropriate metrics (e.g. coverage) for R1 and R2 FASTQs.

    Args:
        r1 (dict): parsed metrics associated with R1 FASTQ
        r2 (dict): parsed metrics associated with R2 FASTQ

    Returns:
        dict: the merged FASTQ metrics
    """
    from statistics import mean
    merged = {
        'qc_stats': {},
        'r1_per_base_quality': r1['per_base_quality'],
        'r2_per_base_quality': r2['per_base_quality'],
        'r1_read_lengths': r1['read_lengths'],
        'r2_read_lengths': r2['read_lengths']
    }
    for key in r1['qc_stats']:
        if key in ['total_bp', 'coverage', 'read_total']:
            merged['qc_stats'][key] = r1['qc_stats'][key] + r2['qc_stats'][key] if r2 else r1['qc_stats'][key]
        else:
            val = mean([r1['qc_stats'][key], r2['qc_stats'][key]]) if r2 else r1['qc_stats'][key]
            merged['qc_stats'][key] = f'{val:.4f}' if isinstance(val, float) else val

    return merged


def is_paired(path: str, name: str) -> bool:
    """
    Check if in input sample had paired-end or single-end reads

    Args:
        path (str): a path to expected Bactopia results
        name (str): the name of sample to test

    Raises:
        ValueError: Processed FASTQ(s) could not be found.

    Returns:
        bool: True: reads are paired, False: reads are single-end
    """
    r1 = f"{path}/{name}/quality-control/{name}_R1.fastq.gz"
    r2 = f"{path}/{name}/quality-control/{name}_R2.fastq.gz"
    se = f"{path}/{name}/quality-control/{name}.fastq.gz"
    if os.path.exists(r1) and os.path.exists(r2):
        return True
    elif os.path.exists(se):
        return False
    else:
        raise ValueError(f"Processed FASTQs not found in {path}/{name}/quality-control/")
 

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
        filename = None
        r1 = None
        r2 = None
        se = None

        if result.endswith('original.json'):
            result_name = 'original'
            r1 = f"{path}/{name}/{RESULT_TYPE}/summary-original/{name}_R1-{result}"
            r2 = f"{path}/{name}/{RESULT_TYPE}/summary-original/{name}_R2-{result}"
            se = f"{path}/{name}/{RESULT_TYPE}/summary-original/{name}-{result}"
        elif result.endswith('final.json'):
            result_name = 'final'
            r1 = f"{path}/{name}/{RESULT_TYPE}/summary-final/{name}_R1-{result}"
            r2 = f"{path}/{name}/{RESULT_TYPE}/summary-final/{name}_R2-{result}"
            se = f"{path}/{name}/{RESULT_TYPE}/summary-final/{name}-{result}"

        if (se):
            if os.path.exists(se):
                parsable_results.append({
                    'result_name': result_name,
                    'files': [se],
                    'optional': False,
                    'missing': False
                })
            else:
                missing = True
                if os.path.exists(r1) and os.path.exists(r2):
                    missing = False
                parsable_results.append({
                    'result_name': result_name,
                    'files': [r1, r2],
                    'optional': False,
                    'missing': missing
                })

    return parsable_results
