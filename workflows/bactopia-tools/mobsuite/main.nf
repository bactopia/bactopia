#!/usr/bin/env nextflow
/**
 * Reconstruction and annotation of plasmids from bacterial genome assemblies.
 *
 * This Bactopia Tool uses [MOB-suite](https://github.com/phac-nml/mob-suite) to identify, reconstruct,
 * and annotate plasmid sequences from draft genome assemblies. It separates plasmid from chromosomal
 * contigs, groups plasmid sequences into clusters, and provides comprehensive plasmid typing
 * and mobility information.
 *
 * @status stable
 * @keywords plasmid, reconstruction, annotation, mobile genetic elements, typing
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, mobsuite, mobsuite_db
 *
 * @subworkflows bactopiatool_init, mobsuite
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish chromosome.fasta           FASTA file containing all contigs identified as chromosomal
 * @publish contig_report.txt          Report assigning each contig to chromosome or plasmid group
 * @publish plasmid_*.fasta            Individual FASTA files for each reconstructed plasmid
 * @publish *-mobtyper.txt            MOB-typer report with plasmid typing and mobility information
 *
 * @section Merged Results
 * @publish mobsuite.tsv              Merged TSV file containing MOB-suite results from all samples
 *
 * @section Execution Logs
 * @publish logs/mobsuite/*           Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                 Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml              Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { MOBSUITE          } from '../../../subworkflows/mobsuite/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    MOBSUITE(BACTOPIATOOL_INIT.out.assembly)
    ch_sample_nf_logs = MOBSUITE.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = MOBSUITE.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    publish:
    sample_outputs = MOBSUITE.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = MOBSUITE.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}
output {
    sample_outputs {
        path { r ->
            r.results  >> "${r.meta.output_dir}/"
            r.logs     >> "${r.meta.logs_dir}/"
            r.versions >> "${r.meta.logs_dir}/"
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
