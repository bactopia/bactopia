//
// tbprofiler - Detect resistance and lineages of Mycobacterium  tuberculosis genomes
//

tbprofiler_args = [
    params.skip_resistance ? "" : "--resistance",
    params.skip_kaptive ? "" : "--kaptive",
    params.force_index ? "--force_index" : "",
    "--min_identity ${params.min_identity}",
    "--min_coverage ${params.min_coverage}",
    "--min_spurious_identity ${params.min_spurious_identity}",
    "--min_spurious_coverage ${params.min_spurious_coverage}",
    "--min_kaptive_confidence ${params.min_kaptive_confidence}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { TBPROFILER as TBPROFILER_MODULE }  from '../../../modules/nf-core/modules/tbprofiler/profile/main' addParams( options: [args: "${tbprofiler_args}", is_module: true] )
if (params.is_subworkflow) {
    include { CSVTK_CONCAT } from '../../../modules/nf-core/modules/csvtk/concat/main' addParams( options: [publish_to_base: true] )
}

workflow TBPROFILER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    TBPROFILER_MODULE(reads)
    if (params.is_subworkflow) {
        TBPROFILER_MODULE.out.csv.collect{meta, csv -> csv}.map{ csv -> [[id:'tbprofiler'], csv]}.set{ ch_merge_tbprofiler }
        CSVTK_CONCAT(ch_merge_tbprofiler, 'csv', 'csv')
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions.first())
    }
    ch_versions = ch_versions.mix(TBPROFILER_MODULE.out.versions.first())

    emit:
    bam = TBPROFILER_MODULE.out.bam
    csv = TBPROFILER_MODULE.out.csv
    merged_csv = CSVTK_CONCAT.out.csv
    json = TBPROFILER_MODULE.out.json
    vcf = TBPROFILER_MODULE.out.vcf
    versions = ch_versions // channel: [ versions.yml ]
}
