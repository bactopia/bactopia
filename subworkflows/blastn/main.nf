//
// blastn - Search against nucleotide BLAST databases using nucleotide queries
//
include { BLAST_BLASTN as BLASTN_MODULE } from '../../modules/blast/blastn/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow BLASTN {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    // Set up query file
    QUERY = params.blastn_query ? file(params.blastn_query) : []

    // Run BLASTN
    BLASTN_MODULE(reads, QUERY)
    ch_versions = ch_versions.mix(BLASTN_MODULE.out.versions)
    ch_logs = ch_logs.mix(BLASTN_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BLASTN_MODULE.out.nf_logs)

    // Merge results
    BLASTN_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'blastn'], tsv]}.set{ ch_merge_blastn }
    CSVTK_CONCAT(ch_merge_blastn, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CSVTK_CONCAT.out.nf_logs)

    emit:
    tsv = BLASTN_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions // channel: [ versions.yml ]
}
