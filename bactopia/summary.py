from collections import OrderedDict
from bactopia.parse import parse_bactopia_directory


def summarize(path: str) -> dict:
    """
    Creates summary reports for a Bactopia directory.

    Args:
        path (str): Path to a directory containing Bactopia results

    Returns:
        dict: [description]
    """
    results = parse_bactopia_directory(path)


def get_rank(cutoff: dict, coverage: float, quality: float, length: int, contigs: int, genome_size: int, is_paired: bool) -> list:
    """
    Determine the rank (gold, silver, bronze, fail) based on user cutoffs.

    Args:
        cutoff (dict): Cutoffs set by users to determine rank
        coverage (float): Estimated coverage of the sample
        quality (float): Per-read average quality
        length (int): Median length of reads
        contigs (int): Total number of contigs
        genome_size (int): Genome size of sample used in analysis
        is_paired (bool): Sample used paired-end reads

    Returns:
        list: the rank and reason for the ranking
    """
    rank = None
    reason = []
    coverage = float(f'{float(coverage):.2f}')
    quality = float(f'{float(quality):.2f}')
    length = round(float(f'{float(length):.2f}'))
    contigs = int(contigs)
    genome_size = int(genome_size)
    gold = cutoff['gold']
    silver = cutoff['silver']
    bronze = cutoff['bronze']

    if coverage >= gold['coverage'] and quality >= gold['quality'] and length >= gold['length'] and contigs <= gold['contigs'] and is_paired:
        reason.append('passed all cutoffs')
        rank = 'gold'
    elif coverage >= silver['coverage'] and quality >= silver['quality'] and length >= silver['length'] and contigs <= silver['contigs'] and is_paired:
        if coverage < gold['coverage']:
            reason.append(f"Low coverage ({coverage:.2f}x, expect >= {gold['coverage']}x)")
        if quality < gold['quality']:
            reason.append(f"Poor read quality (Q{quality:.2f}, expect >= Q{gold['quality']})")
        if length < gold['length']:
            reason.append(f"Short read length ({length}bp, expect >= {gold['length']} bp)")
        if contigs > gold['contigs']:
            reason.append(f"Too many contigs ({contigs}, expect <= {gold['contigs']})")
        rank = 'silver'
    elif coverage >= bronze['coverage'] and quality >= bronze['quality'] and length >= bronze['length'] and contigs <= bronze['contigs']:
        if coverage < silver['coverage']:
            reason.append(f"Low coverage ({coverage:.2f}x, expect >= {silver['coverage']}x)")
        if quality < silver['quality']:
            reason.append(f"Poor read quality (Q{quality:.2f}, expect >= Q{silver['quality']})")
        if length < silver['length']:
            reason.append(f"Short read length ({length}bp, expect >= {silver['length']} bp)")
        if contigs > silver['contigs']:
            reason.append(f"Too many contigs ({contigs}, expect <= {silver['contigs']})")
        if not is_paired:
            reason.append(f"Single-end reads")
        rank = 'bronze'

    if not rank:
        rank = 'exclude'

    if coverage < bronze['coverage']:
        reason.append(f"Low coverage ({coverage:.2f}x, expect >= {bronze['coverage']}x)")
    if quality < bronze['quality']:
        reason.append(f"Poor read quality (Q{quality:.2f}, expect >= Q{bronze['quality']})")
    if length < bronze['length']:
        reason.append(f"Short read length ({length:.2f}bp, expect >= {bronze['length']} bp)")
    if contigs > bronze['contigs']:
        reason.append(f"Too many contigs ({contigs}, expect <= {bronze['contigs']})")

    if cutoff['min-assembled-size']:
        if genome_size < cutoff['min-assembled-size']:
            reason.append(f"Assembled size is too small ({genome_size} bp, expect <= {cutoff['min-assembled-size']})")

    if cutoff['max-assembled-size']:
        if genome_size < cutoff['max-assembled-size']:
            reason.append(f"Assembled size is too large ({genome_size} bp, expect <= {cutoff['max-assembled-size']})")

    reason = ";".join(sorted(reason))
    return [rank, reason]


def print_failed(failed: list, spaces: int = 8) -> str:
    """
    Format the strings of samples that failed

    Args:
        failed (list): A list of samples that failed for a particular reason
        spaces (int, optional): Total number of spaces to indent. Defaults to 8.

    Returns:
        str: The set of formatted strings
    """
    lines = []
    for key, val in sorted(failed.items()):
        if key != 'failed-cutoff':
            lines.append(f'{spaces * " "}{key.replace("-", " ").title()}: {len(val)}')
    return "\n".join(lines) if lines else ""


def print_cutoffs(cutoffs: list, spaces: int = 8) -> str:
    """
    Format strings for samples that failed a cutoff.

    Args:
        cutoffs (list): A list of samples that failed for a cutoff
        spaces (int, optional): Total number of spaces to indent. Defaults to 8.

    Returns:
        str:  The set of formatted strings
    """
    lines = []
    for key, val in sorted(cutoffs.items()):
        lines.append(f'{spaces * " "}{key}: {val}')
    return "\n".join(lines) if lines else ""


def gather_results(sample: dict, rank: str, reason: str) -> dict:
    """
    Aggregate results into an unnested dictionary.

    Args:
        sample (dict): The results associated with a sample
        rank (str): The rank of a sample
        reason (str): The reason for the given rank

    Returns:
        dict: An unnested dictionary of results
    """
    results = OrderedDict((
        ('sample', sample['sample']),
        ('is_paired', sample['is_paired']),
        ('rank', rank),
        ('reason', reason),
        ('estimated_genome_size', sample['genome_size'])
    ))
    results.update(sample['results']['assembly']['stats'])
    results['checkm_lineage'] = sample['results']['assembly']["checkm"]["Marker lineage"]
    results['checkm_completeness'] = sample['results']['assembly']["checkm"]["Completeness"]
    results['checkm_contamination'] = sample['results']['assembly']["checkm"]["Contamination"]
    results['checkm_heterogeneity'] = sample['results']['assembly']["checkm"]["Strain heterogeneity"]
    results.update(_prefix_keys(sample['results']['quality-control']['original']['qc_stats'], 'original'))
    results.update(_prefix_keys(sample['results']['quality-control']['final']['qc_stats'], 'final'))
    results.update(_remove_keys(sample['results']['annotation']['stats'], ['organism', 'contigs', 'bases']))
    results.update(_add_minmers(sample['results']['minmers']))
    results.update(_add_mlst(sample['results']['mlst']))
    return results


def _add_mlst(mlsts: dict) -> dict:
    """
    Read through MLST results and create column each schema.

    Args:
        mlsts (dict): The MLST results associated with a sample

    Returns:
        dict: Per schema MLST hits
    """
    results = OrderedDict()
    for key, vals in mlsts.items():
        schema, tool = key.split('-')
        prefix = f"mlst_{tool}" if schema == "default" else f"mlst_{schema}_{tool}"

        if tool == "blast":
            results[f"{prefix}_st"] = vals['ST']['st']
            results[f"{prefix}_loci"] = len(vals) - 1
            results[f"{prefix}_perfect_matches"] = vals['ST']['perfect_matches']
        else:
            for k, v in vals.items():
                results[f"{prefix}_{k.lower()}"] = v
    return results


def _add_minmers(minmers: dict) -> dict:
    """
    Read through minmer results and create column for top hit.

    Args:
        minmers (dict): Mash and Sourmash results against RefSeq and GenBank

    Returns:
        dict: Top hit description for each set of databases
    """
    results = OrderedDict()
    for key in ['refseq-k21', 'genbank-k21', 'genbank-k31', 'genbank-k51']:
        if key in minmers:
            prefix = key.replace('-', '_')
            if len(minmers[key]):
                if key.startswith('genbank'):
                    # Sourmash keys: "overlap", "p_query", "p_match", "match"
                    if minmers[key]['matches']:
                        results[f'{prefix}_match'] = minmers[key]['matches'][0]['match'].split("(")[0].rstrip()
                        results[f'{prefix}_overlap'] = minmers[key]['matches'][0]['overlap']
                        results[f'{prefix}_p_query'] = minmers[key]['matches'][0]['p_query']
                        results[f'{prefix}_p_match'] = minmers[key]['matches'][0]['p_match']
                    else:
                        results[f'{prefix}_match'] = None
                        results[f'{prefix}_overlap'] = None
                        results[f'{prefix}_p_query'] = None
                        results[f'{prefix}_p_match'] = None

                    results[f'{prefix}_no_assignment'] = minmers[key]['no_assignment']
                    results[f'{prefix}_total'] = len(minmers[key]['matches'])
                    
                else:
                    # Mash keys: "identity", "shared-hashes", "median-multiplicity", "p-value", "query-ID", "query-comment"
                    results[f'{prefix}_id'] = minmers[key][0]['query-ID']
                    results[f'{prefix}_identity'] = minmers[key][0]['identity']
                    results[f'{prefix}_shared_hashes'] = minmers[key][0]['shared-hashes']
                    results[f'{prefix}_median_multiplicity'] = minmers[key][0]['median-multiplicity']
                    results[f'{prefix}_p_value'] = minmers[key][0]['p-value']
                    results[f'{prefix}_comment'] = minmers[key][0]['query-comment']
                    results[f'{prefix}_total'] = len(minmers[key])
    return results


def _remove_keys(results: dict, remove: list) -> dict:
    """
    Remove a set of keys from a dictionary.

    Args:
        results (dict): The dictionary of results
        remove (list): the keys to remove

    Returns:
        dict: The altered dictionary
    """
    removed = {}
    for key, val in results.items():
        if key not in remove:
            removed[key] = val
    return removed


def _prefix_keys(results: dict, prefix: str) -> dict:
    """
    Add a prefix to existing keys

    Args:
        results (dict): The dictionary of results
        prefix (str): A string to prefix each key with

    Returns:
        dict: The result dictionary with prefixed keys.
    """
    prefixed = {}
    for key, val in results.items():
        prefixed[f'{prefix}_{key}'] = val
    return prefixed
