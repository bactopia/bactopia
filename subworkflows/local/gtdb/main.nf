//
// gtdb - Identify marker genes and assign taxonomic classifications
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'gtdb')
options.is_module = params.wf == 'gtdb' ? true : false

classify_args = [
    params.gtdb_use_scratch ? "--scratch_dir ${params.gtdb_tmp}" : "",
    params.gtdb_debug ? "--debug" : "",
    params.force_gtdb ? "--force" : "",
    "--tmpdir ${params.gtdb_tmp}",
    "--min_perc_aa ${params.min_perc_aa}",
    "--min_af ${params.min_af}",
].join(' ').replaceAll("\\s{2,}", " ").trim()
DATABASE_DIR = ! params.download_gtdb ? file(params.gtdb) : []
include { GTDBTK_SETUPDB as SETUPDB } from '../../../modules/nf-core/gtdbtk/setupdb/main' addParams( options: options + [publish_to_base: true] )
include { GTDBTK_CLASSIFYWF as CLASSIFY } from '../../../modules/nf-core/gtdbtk/classifywf/main' addParams( options: options + [args: "${classify_args}"] )

workflow GTDB {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()

    if (params.download_gtdb) {
        // Force CLASSIFY to wait
        SETUPDB()

        if (params.gtdb_save_as_tarball) {
            CLASSIFY(fasta, SETUPDB.out.db_tarball)
        } else {
            CLASSIFY(fasta, SETUPDB.out.db)
        }
    } else {
        CLASSIFY(fasta, DATABASE_DIR)
    }
    ch_versions = ch_versions.mix(CLASSIFY.out.versions.first())

    emit:
    results = CLASSIFY.out.results
    versions = ch_versions // channel: [ versions.yml ]
}
