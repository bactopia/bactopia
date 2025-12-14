/**
 * Gene Allele Mutation Microbial Assessment.
 *
 * This subworkflow performs rapid identification, classification, and annotation of
 * translated gene matches from sequencing data using [GAMMA](https://github.com/rastanton/GAMMA).
 * The tool screens input sequences against a protein database to identify gene
 * variants, mutations, and allele types, providing detailed annotation and classification.
 *
 * @status stable
 * @keywords gene, allele, mutation, variant, antimicrobial resistance
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation gamma
 *
 * @modules csvtk_concat, gamma
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for gene allele identification
 *
 * @input db
 * Protein database file for sequence comparison (required)
 *
 * @output gamma        GAMMA allele identification and annotation results
 * @output merged_gamma Combined TSV file containing GAMMA results from all samples
 * @output psl          PSL (Pattern Space Layout) output format
 * @output fasta        Matched protein sequences in FASTA format
 * @output gff          Annotated features in GFF3 format
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { GAMMA as GAMMA_MODULE } from '../../modules/gamma/main'
include { CSVTK_CONCAT          } from '../../modules/csvtk/concat/main'
include { flattenPaths          } from 'plugin/nf-bactopia'
include { gather                } from 'plugin/nf-bactopia'

workflow GAMMA {
    take:
    assembly: Channel<Tuple<Map, Path>>
    db: Path

    main:
    GAMMA_MODULE(assembly, db)
    CSVTK_CONCAT(gather(GAMMA_MODULE.out.gamma, 'gamma'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    gamma: Channel<Tuple<Map, Set<Path>>> = GAMMA_MODULE.out.gamma
    merged_gamma: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    psl: Channel<Tuple<Map, Set<Path>>> = GAMMA_MODULE.out.psl
    fasta: Channel<Tuple<Map, Set<Path>>> = GAMMA_MODULE.out.fasta
    gff: Channel<Tuple<Map, Set<Path>>> = GAMMA_MODULE.out.gff

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.gamma,
        CSVTK_CONCAT.out.csv,
        GAMMA_MODULE.out.psl,
        GAMMA_MODULE.out.fasta,
        GAMMA_MODULE.out.gff
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        GAMMA_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
