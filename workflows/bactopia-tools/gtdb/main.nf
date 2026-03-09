#!/usr/bin/env nextflow
/**
 * Identify marker genes and assign taxonomic classifications using GTDB.
 *
 * This Bactopia Tool uses [GTDB-Tk's](https://github.com/Ecogenomics/GTDBTk) classify
 * workflow to assign taxonomic classifications to samples using the
 * [Genome Taxonomy Database](https://gtdb.ecogenomic.org/).
 *
 * @status stable
 * @keywords taxonomy, classification, marker genes, phylogeny, gtdb, bactopia-tool
 * @tags complexity:complex input-type:parameter output-type:multiple features:bactopia-tool,taxonomy,database-dependent
 * @citation gtdb
 *
 * @subworkflows bactopiatool_init, gtdb
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input gtdb
 * Path to GTDB database directory
 *
 * @input download_gtdb
 * Download GTDB database if not found
 *
 * @input gtdb_save_as_tarball
 * Save GTDB database as tarball after download
 *
 * @section Taxonomic Classification
 * @publish *.summary.tsv        Taxonomic classification summary
 * @publish *.gtdbtk.tsv        Detailed GTDB-Tk results
 *
 * @section Merged Results
 * @publish gtdb.tsv            Merged TSV file containing GTDB results from all samples
 *
 * @section Execution Logs
 * @publish logs/**             Tool execution logs
 * @publish logs/nf-*           Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml        Software version information
 */
nextflow.preview.types = true

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    gtdb                 : Path
    download_gtdb        : Boolean
    gtdb_save_as_tarball : Boolean
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { GTDB              } from '../../../subworkflows/gtdb/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    GTDB(
        BACTOPIATOOL_INIT.out.assembly,
        params.gtdb,
        params.download_gtdb,
        params.gtdb_save_as_tarball
    )

    ch_sample_nf_logs = GTDB.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = GTDB.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }

    publish:
    sample_outputs = GTDB.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = GTDB.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    sample_outputs {
        path { r ->
            r.results      >> "${r.meta.output_dir}/"
            r.supplemental >> "${r.meta.output_dir}/"
            r.logs         >> "${r.meta.logs_dir}/"
            r.versions     >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }
    run_outputs {
        path { r ->
            r.results  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
