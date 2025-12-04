//
// assembler - Assembly of Illumina and ONT reads
//
nextflow.preview.types = true

include { ASSEMBLER as ASSEMBLER_MODULE } from '../../../modules/bactopia/assembler/main'
include { CSVTK_CONCAT                  } from '../../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow ASSEMBLER {
    take:
    reads: Channel<Tuple<Map, Set<Path>, Set<Path>>>

    main:
    ASSEMBLER_MODULE(reads)
    CSVTK_CONCAT(gather(ASSEMBLER_MODULE.out.tsv, 'assembly-scan'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    fna: Channel<Tuple<Map, Set<Path>>> = ASSEMBLER_MODULE.out.fna
    fna_fq: Channel<Tuple<Map, Set<Path>, Set<Path>>> = ASSEMBLER_MODULE.out.fna_fq
    tsv: Channel<Tuple<Map, Path>> = ASSEMBLER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ASSEMBLER_MODULE.out.fna,
        ASSEMBLER_MODULE.out.tsv,
        ASSEMBLER_MODULE.out.error,
        ASSEMBLER_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ASSEMBLER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ASSEMBLER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        ASSEMBLER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
