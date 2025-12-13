/**
 * Systematically search for anti-phage defense systems.
 *
 * This subworkflow uses [DefenseFinder](https://github.com/mdmparis/defense-finder) to identify and classify
 * anti-phage defense systems in bacterial genomes. It detects defense genes, HMM hits, and complete
 * defense systems, providing comprehensive analysis of bacterial antiviral mechanisms.
 *
 * @status stable
 * @keywords bacteria, assembly, anti-phage, defense systems, immunity
 * @tags complexity:complex input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation defensefinder
 *
 * @modules defensefinder_run, defensefinder_update, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format for defense system detection
 *
 * @output genes_tsv          Individual anti-phage genes detected in TSV format
 * @output merged_genes_tsv   Combined anti-phage genes from all samples in a single TSV file
 * @output hmmer_tsv          HMM search results for defense system proteins in TSV format
 * @output merged_hmmer_tsv   Combined HMM results from all samples in a single TSV file
 * @output systems_tsv        Complete anti-phage defense systems identified in TSV format
 * @output merged_systems_tsv Combined defense systems from all samples in a single TSV file
 * @output proteins           Defense system protein sequences in FASTA format
 * @output proteins_index     Diamond protein index for rapid similarity searches
 * @output macsydata_raw      Raw MacSyFinder output files for detailed analysis
 * @output results            Aggregated results channel containing all output files
 * @output logs               Aggregated logs channel containing all execution logs
 * @output nf_logs            Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions           Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { DEFENSEFINDER_UPDATE           } from '../../modules/defensefinder/update/main'
include { DEFENSEFINDER_RUN              } from '../../modules/defensefinder/run/main'
include { CSVTK_CONCAT as GENES_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as HMMER_CONCAT   } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as SYSTEMS_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths                   } from 'plugin/nf-bactopia'
include { gather                         } from 'plugin/nf-bactopia'

workflow DEFENSEFINDER {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    DEFENSEFINDER_UPDATE()
    DEFENSEFINDER_RUN(assembly, DEFENSEFINDER_UPDATE.out.db)

    // Merge results
    GENES_CONCAT(gather(DEFENSEFINDER_RUN.out.genes_tsv, 'defensefinder-genes'), 'tsv', 'tsv')
    HMMER_CONCAT(gather(DEFENSEFINDER_RUN.out.hmmer_tsv, 'defensefinder-hmmer'), 'tsv', 'tsv')
    SYSTEMS_CONCAT(gather(DEFENSEFINDER_RUN.out.systems_tsv, 'defensefinder-systems'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    genes_tsv: Channel<Tuple<Map, Set<Path>>> = DEFENSEFINDER_RUN.out.genes_tsv
    merged_genes_tsv: Channel<Tuple<Map, Set<Path>>> = GENES_CONCAT.out.csv
    hmmer_tsv: Channel<Tuple<Map, Set<Path>>> = DEFENSEFINDER_RUN.out.hmmer_tsv
    merged_hmmer_tsv: Channel<Tuple<Map, Set<Path>>> = HMMER_CONCAT.out.csv
    systems_tsv: Channel<Tuple<Map, Set<Path>>> = DEFENSEFINDER_RUN.out.systems_tsv
    merged_systems_tsv: Channel<Tuple<Map, Set<Path>>> = SYSTEMS_CONCAT.out.csv
    proteins: Channel<Tuple<Map, Set<Path>>> = DEFENSEFINDER_RUN.out.proteins
    proteins_index: Channel<Tuple<Map, Set<Path>>> = DEFENSEFINDER_RUN.out.proteins_index
    macsydata_raw: Channel<Tuple<Map, Set<Path>>> = DEFENSEFINDER_RUN.out.macsydata_raw

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        DEFENSEFINDER_RUN.out.genes_tsv,
        GENES_CONCAT.out.csv,
        DEFENSEFINDER_RUN.out.hmmer_tsv,
        HMMER_CONCAT.out.csv,
        DEFENSEFINDER_RUN.out.systems_tsv,
        SYSTEMS_CONCAT.out.csv,
        DEFENSEFINDER_RUN.out.proteins,
        DEFENSEFINDER_RUN.out.proteins_index,
        DEFENSEFINDER_RUN.out.macsydata_raw
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        DEFENSEFINDER_RUN.out.logs,
        GENES_CONCAT.out.logs,
        HMMER_CONCAT.out.logs,
        SYSTEMS_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        DEFENSEFINDER_RUN.out.nf_logs,
        GENES_CONCAT.out.nf_logs,
        HMMER_CONCAT.out.nf_logs,
        SYSTEMS_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        DEFENSEFINDER_RUN.out.versions,
        GENES_CONCAT.out.versions,
        HMMER_CONCAT.out.versions,
        SYSTEMS_CONCAT.out.versions
    ])
}
