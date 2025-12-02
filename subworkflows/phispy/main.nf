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
    information: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.information
    bacteria_fasta: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.bacteria_fasta
    bacteria_gbk: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.bacteria_gbk
    phage_fasta: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.phage_fasta
    phage_gbk: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.phage_gbk
    prophage_gff: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.prophage_gff
    prophage_tbl: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.prophage_tbl
    prophage_tsv: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.prophage_tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PHISPY_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        PHISPY_MODULE.out.information,
        PHISPY_MODULE.out.bacteria_fasta,
        PHISPY_MODULE.out.bacteria_gbk,
        PHISPY_MODULE.out.phage_fasta,
        PHISPY_MODULE.out.phage_gbk,
        PHISPY_MODULE.out.prophage_gff,
        PHISPY_MODULE.out.prophage_tbl,
        PHISPY_MODULE.out.prophage_tsv
    ])
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
