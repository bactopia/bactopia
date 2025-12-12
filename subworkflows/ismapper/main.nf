/**
 * Identify transposase insertion sites in bacterial genomes.
 *
 * This subworkflow maps insertion sequence (IS) positions in bacterial genomes
 * using [ISMapper](https://github.com/jhawkey/IS_mapper). The tool identifies
 * transposase insertion sites from short read sequence data by mapping reads
 * to reference sequences and detecting insertion sites with high precision.
 *
 * @status stable
 * @keywords insertion, sequence, transposase, mobile genetic elements
 * @tags complexity:simple input-type:multiple output-type:single
 * @citation ismapper
 *
 * @modules ismapper
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Paired-end sequencing reads for IS element mapping
 *
 * @input reference
 * Reference genome in FASTA format for mapping
 *
 * @input insertions
 * Insertion sequence reference file containing IS elements to map
 *
 * @output results     ISMapper results including insertion site coordinates and supporting read information
 * @output logs        Execution logs from the ISMapper analysis
 * @output nf_logs     Nextflow execution scripts and logs for debugging
 * @output versions    Software version information
 */


nextflow.preview.types = true

include { ISMAPPER as ISMAPPER_MODULE } from '../../modules/ismapper/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow ISMAPPER {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    reference: Path
    insertions: Path

    main:
    ISMAPPER_MODULE(reads, reference, insertions)

    emit:
    results: Channel<Tuple<Map, Path>> = flattenPaths([ISMAPPER_MODULE.out.supplemental])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([ISMAPPER_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([ISMAPPER_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([ISMAPPER_MODULE.out.versions])
}
