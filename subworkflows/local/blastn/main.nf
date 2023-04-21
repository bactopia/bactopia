//
// blastn - Search against nucleotide BLAST databases using nucleotide queries
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'blastn')
options.args = [
    "-outfmt '6 ${params.blastn_outfmt}'",
    "-perc_identity ${params.blastn_perc_identity}",
    "-qcov_hsp_perc ${params.blastn_qcov_hsp_perc}",
    "-max_target_seqs ${params.blastn_max_target_seqs}",
    params.blastn_opts
].join(' ').replaceAll("\\s{2,}", " ").trim()
QUERY = params.blastn_query ? file(params.blastn_query) : []

include { BLAST_BLASTN as BLASTN_MODULE } from '../../../modules/nf-core/blast/blastn/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'blastn-concat', process_name: params.merge_folder] )

workflow BLASTN {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    // Run BLASTN
    BLASTN_MODULE(reads, QUERY)
    ch_versions = ch_versions.mix(BLASTN_MODULE.out.versions)

    // Merge results
    BLASTN_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'blastn'], tsv]}.set{ ch_merge_blastn }
    CSVTK_CONCAT(ch_merge_blastn, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = BLASTN_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
