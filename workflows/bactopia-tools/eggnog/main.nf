#!/usr/bin/env nextflow
/**
 * Functional annotation of proteins using orthologous groups and phylogenies.
 *
 * This Bactopia Tool uses [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) to assign
 * functional annotation to protein sequences. eggNOG-mapper uses orthologous groups and phylogenies
 * from the eggNOG database to more precisely functionally annotate than traditional homology methods.
 *
 * @status stable
 * @keywords functional annotation, orthology, proteins, eggnog, bactopia-tool
 * @tags complexity:complex input-type:parameter output-type:multiple features:bactopia-tool,functional-annotation,orthology,database-dependent
 * @citation eggnog
 *
 * @subworkflows bactopiatool_init, eggnog
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input eggnog_db
 * Path to eggNOG database for functional annotation
 *
 * @input download_eggnog
 * Download eggNOG database if not found locally
 *
 * @section Annotation
 * @publish *.emapper.annotations      Results from the annotation phase
 * @publish *.emapper.hits             Results from the search phase (HMMER, Diamond or MMseqs2)
 * @publish *.emapper.seed_orthologs   Results from parsing the hits
 * @publish *.emapper.annotations.xlsx Annotations in Excel format
 * @publish *.emapper.orthologs        List of orthologs found for each query
 * @publish *.emapper.genepred.fasta   Sequences of predicted CDS
 * @publish *.emapper.gff              GFF of predicted CDS
 * @publish *.emapper.no_annotations.fasta Sequences without annotation
 * @publish *.emapper.pfam             Positions of PFAM domains identified
 *
 * @section Execution Logs
 * @publish logs/**                    Tool execution logs
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    eggnog_db       : Path
    download_eggnog : Boolean
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { EGGNOG            } from '../../../subworkflows/eggnog/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    EGGNOG(
        BACTOPIATOOL_INIT.out.proteins,
        params.eggnog_db,
        params.download_eggnog,
        params.eggnog_save_as_tarball
    )

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = EGGNOG.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = EGGNOG.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = EGGNOG.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = EGGNOG.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
