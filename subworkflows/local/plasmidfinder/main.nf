//
// plasmidfinder - Plasmid identification from assemblies
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'plasmidfinder')
options.args = [
    "-l ${params.pf_mincov}",
    "-t ${params.pf_threshold}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { PLASMIDFINDER as PLASMIDFINDER_MODULE } from '../../../modules/nf-core/plasmidfinder/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'plasmidfinder-concat', process_name: params.merge_folder] )

workflow PLASMIDFINDER {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    PLASMIDFINDER_MODULE(fasta)
    ch_versions = ch_versions.mix(PLASMIDFINDER_MODULE.out.versions.first())

    // Merge results
    PLASMIDFINDER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'plasmidfinder'], tsv]}.set{ ch_merge_plasmidfinder }
    CSVTK_CONCAT(ch_merge_plasmidfinder, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    json = PLASMIDFINDER_MODULE.out.json
    txt = PLASMIDFINDER_MODULE.out.txt
    tsv = PLASMIDFINDER_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    genome_seq = PLASMIDFINDER_MODULE.out.genome_seq
    plasmid_seq = PLASMIDFINDER_MODULE.out.plasmid_seq
    versions = ch_versions
}
