//
// staphopiasccmec - Primer based SCCmec typing of Staphylococcus aureus genomes
//
nextflow.preview.types = true

include { STAPHOPIASCCMEC as STAPHOPIASCCMEC_MODULE } from '../../modules/staphopiasccmec/main'
include { CSVTK_CONCAT                              } from '../../modules/csvtk/concat/main'
include { flattenPaths                              } from 'plugin/nf-bactopia'
include { gather                                    } from 'plugin/nf-bactopia'

workflow STAPHOPIASCCMEC {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    STAPHOPIASCCMEC_MODULE(fasta)
    CSVTK_CONCAT(gather(STAPHOPIASCCMEC_MODULE.out.tsv, 'staphopiasccmec'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = STAPHOPIASCCMEC_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        STAPHOPIASCCMEC_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        STAPHOPIASCCMEC_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        STAPHOPIASCCMEC_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        STAPHOPIASCCMEC_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
