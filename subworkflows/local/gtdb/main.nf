//
// gtdb - Identify marker genes and assign taxonomic classifications
//
classify_args = [
    params.gtdb_use_scratch ? "--scratch_dir ${params.gtdb_tmp}" : "",
    params.gtdb_debug ? "--debug" : "",
    params.force_gtdb ? "--force" : "",
    "--tmpdir ${params.gtdb_tmp}",
    "--min_perc_aa ${params.min_perc_aa}",
    "--min_af ${params.min_af}",
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { GTDBTK_SETUPDB as SETUPDB } from '../../../modules/nf-core/modules/gtdbtk/setupdb/main' addParams( options: [publish_to_base: true, is_module: true] )
include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../../../modules/nf-core/modules/gtdbtk/classifywf/main' addParams( options: [args: "${classify_args}", is_module: true] )

workflow GTDB {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    if (params.download_gtdb) {
        // Force CLASSIFY to wait
        SETUPDB()
        CLASSIFY(fasta, SETUPDB.out.db)
    } else {
        CLASSIFY(fasta, file("${params.gtdb}/*"))
    }
    ch_versions = ch_versions.mix(CLASSIFY.out.versions.first())

    emit:
    results = CLASSIFY.out.results
    versions = ch_versions // channel: [ versions.yml ]
}
