//
// blastp - Search against protein BLAST databases using protein queries
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'blastp')
options.args = [
    "-outfmt '6 ${params.blastp_outfmt}'",
    "-qcov_hsp_perc ${params.blastp_qcov_hsp_perc}",
    "-max_target_seqs ${params.blastp_max_target_seqs}",
    params.blastp_opts
].join(' ').replaceAll("\\s{2,}", " ").trim()
QUERY = params.blastp_query ? file(params.blastp_query) : []

include { BLAST_BLASTP as BLASTP_MODULE } from '../../../modules/nf-core/blast/blastp/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'blastp-concat', process_name: params.merge_folder] )

workflow BLASTP {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    // Run BLASTP
    BLASTP_MODULE(reads, QUERY)
    ch_versions = ch_versions.mix(BLASTP_MODULE.out.versions)

    // Merge results
    BLASTP_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'blastp'], tsv]}.set{ ch_merge_blastp }
    CSVTK_CONCAT(ch_merge_blastp, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = BLASTP_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
