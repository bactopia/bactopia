//
// tbprofiler - Detect resistance and lineages of Mycobacterium  tuberculosis genomes
//
include { initOptions } from '../../../lib/nf/functions'

// tbprofiler profile
profiler_opts = initOptions(params.containsKey("options") ? params.options : [:], 'tbprofiler')
profiler_opts.args = [
    params.call_whole_genome ? "--call_whole_genome" : "",
    params.calling_params ? "--calling_params ${params.calling_params}" : "",
    params.suspect ? "--suspect" : "",
    params.no_flagstat ? "--no_flagstat" : "",
    params.no_delly ? "--no_delly" : "",
    "--mapper ${params.mapper}",
    "--caller ${params.caller}",
    params.tbprofiler_opts ? "${params.calling_params}" : "",
].join(' ').replaceAll("\\s{2,}", " ").trim()

// tb-profiler collate options
collate_opts = initOptions(params.containsKey("options") ? params.options : [:], 'tbprofiler-collate')
collate_opts.is_module = false
collate_opts.args = [
    params.itol ? "--itol" : "",
    params.full ? "--full" : "",
    params.all_variants ? "--all_variants" : "",
    params.mark_missing ? "--mark_missing" : "",
].join(' ').replaceAll("\\s{2,}", " ").trim()
collate_opts.logs_subdir = 'tbprofiler-collate'
collate_opts.process_name = params.merge_folder

include { TBPROFILER_PROFILE }  from '../../../modules/nf-core/tbprofiler/profile/main' addParams( options: profiler_opts )
include { TBPROFILER_COLLATE }  from '../../../modules/nf-core/tbprofiler/collate/main' addParams( options: collate_opts )

workflow TBPROFILER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    TBPROFILER_PROFILE(reads)
    ch_versions = ch_versions.mix(TBPROFILER_PROFILE.out.versions)

    // Merge results
    TBPROFILER_PROFILE.out.json.collect{meta, json -> json}.map{ json -> [[id:'tbprofiler'], json]}.set{ ch_merge_tbprofiler }
    TBPROFILER_COLLATE(ch_merge_tbprofiler)
    ch_versions = ch_versions.mix(TBPROFILER_COLLATE.out.versions)

    emit:
    bam = TBPROFILER_PROFILE.out.bam
    csv = TBPROFILER_PROFILE.out.csv
    merged_csv = TBPROFILER_COLLATE.out.csv
    json = TBPROFILER_PROFILE.out.json
    txt = TBPROFILER_PROFILE.out.txt
    vcf = TBPROFILER_PROFILE.out.vcf
    versions = ch_versions // channel: [ versions.yml ]
}
