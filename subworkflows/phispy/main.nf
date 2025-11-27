//
// phispy - Predict prophages in bacterial genomes
//
nextflow.preview.types = true

include { PHISPY as PHISPY_MODULE } from '../../modules/phispy/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow PHISPY {
    take:
    gbk: Channel<Tuple<Map, Set<Path>>>

    main:
    PHISPY_MODULE(gbk)
    CSVTK_CONCAT(gather(PHISPY_MODULE.out.tsv, 'phispy'), 'tsv', 'tsv')

    emit:
    tsv: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PHISPY_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PHISPY_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        PHISPY_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
