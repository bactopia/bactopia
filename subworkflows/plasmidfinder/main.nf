//
// plasmidfinder - Plasmid identification from assemblies
//
nextflow.preview.types = true

include { PLASMIDFINDER as PLASMIDFINDER_MODULE } from '../../modules/plasmidfinder/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow PLASMIDFINDER {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    PLASMIDFINDER_MODULE(fasta)

    // Merge results
    ch_merge_plasmidfinder = PLASMIDFINDER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'plasmidfinder'], tsv]}
    CSVTK_CONCAT(ch_merge_plasmidfinder, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = PLASMIDFINDER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    json = PLASMIDFINDER_MODULE.out.json
    txt = PLASMIDFINDER_MODULE.out.txt
    genome_seq = PLASMIDFINDER_MODULE.out.genome_seq
    plasmid_seq = PLASMIDFINDER_MODULE.out.plasmid_seq

    // Generic aggregate outputs
    results = PLASMIDFINDER_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv,
        PLASMIDFINDER_MODULE.out.json,
        PLASMIDFINDER_MODULE.out.txt,
        PLASMIDFINDER_MODULE.out.genome_seq,
        PLASMIDFINDER_MODULE.out.plasmid_seq
    )
    logs = PLASMIDFINDER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = PLASMIDFINDER_MODULE.out.nf_begin.mix(
        PLASMIDFINDER_MODULE.out.nf_err,
        PLASMIDFINDER_MODULE.out.nf_log,
        PLASMIDFINDER_MODULE.out.nf_out,
        PLASMIDFINDER_MODULE.out.nf_run,
        PLASMIDFINDER_MODULE.out.nf_sh,
        PLASMIDFINDER_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = PLASMIDFINDER_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
