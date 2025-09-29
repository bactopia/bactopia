//
// eggnog - Functional annotation of proteins using orthologous groups and phylogenies
//
include { EGGNOG_DOWNLOAD } from '../../modules/eggnog/download/main'
include { EGGNOG_MAPPER } from '../../modules/eggnog/mapper/main'

workflow EGGNOG {
    take:
    faa // channel: [ val(meta), [ fasta ] ]
    database
    download_eggnog

    main:
    if (download_eggnog) {
        // Force EGGNOG_MAPPER to wait
        EGGNOG_DOWNLOAD()
        EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
    } else {
        EGGNOG_MAPPER(faa, database)
    }

    emit:
    // Individual outputs
    hits = EGGNOG_MAPPER.out.hits
    seed_orthologs = EGGNOG_MAPPER.out.seed_orthologs
    annotations = EGGNOG_MAPPER.out.annotations
    xlsx = EGGNOG_MAPPER.out.xlsx
    orthologs = EGGNOG_MAPPER.out.orthologs
    genepred = EGGNOG_MAPPER.out.genepred
    gff = EGGNOG_MAPPER.out.gff
    no_anno = EGGNOG_MAPPER.out.no_anno
    pfam = EGGNOG_MAPPER.out.pfam

    // Generic aggregate outputs
    results = EGGNOG_MAPPER.out.hits.mix(
        EGGNOG_MAPPER.out.seed_orthologs,
        EGGNOG_MAPPER.out.annotations,
        EGGNOG_MAPPER.out.xlsx,
        EGGNOG_MAPPER.out.orthologs,
        EGGNOG_MAPPER.out.genepred,
        EGGNOG_MAPPER.out.gff,
        EGGNOG_MAPPER.out.no_anno,
        EGGNOG_MAPPER.out.pfam
    )
    logs = EGGNOG_MAPPER.out.logs
    nf_logs = EGGNOG_MAPPER.out.nf_begin.mix(
        EGGNOG_MAPPER.out.nf_err,
        EGGNOG_MAPPER.out.nf_log,
        EGGNOG_MAPPER.out.nf_out,
        EGGNOG_MAPPER.out.nf_run,
        EGGNOG_MAPPER.out.nf_sh,
        EGGNOG_MAPPER.out.nf_trace
    )
    versions = EGGNOG_MAPPER.out.versions
}
