/**
 * Calculate pairwise SNP distances from sequence alignments.
 *
 * This subworkflow uses [snp-dists](https://github.com/tseemann/snp-dists) to compute
 * pairwise SNP distance matrices from multiple sequence alignments. It reads an
 * alignment file (typically a core-genome alignment) and calculates the number
 * of SNP differences between each pair of sequences, producing a distance matrix
 * useful for phylogenetic and epidemiological analyses.
 *
 * @status stable
 * @keywords snp, distance, alignment, phylogeny, core-genome
 * @tags complexity:simple input-type:single output-type:multiple
 * @citation snpdists
 *
 * @modules snpdists
 *
 * @input tuple(meta, alignment)
 * - `meta`: Groovy Map containing sample information
 * - `alignment`: Multiple sequence alignment in FASTA format
 *
 * @output tsv             Pairwise SNP distance matrix in TSV format
 * @output results         Aggregated results channel containing all output files
 * @output logs            Aggregated logs channel containing all execution logs
 * @output nf_logs         Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions        Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SNPDISTS as SNPDISTS_MODULE } from '../../modules/snpdists/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SNPDISTS {
    take:
    alignment: Channel<Tuple<Map, Set<Path>>>

    main:
    SNPDISTS_MODULE(alignment)

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([SNPDISTS_MODULE.out.tsv])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([SNPDISTS_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([SNPDISTS_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([SNPDISTS_MODULE.out.versions])
}
