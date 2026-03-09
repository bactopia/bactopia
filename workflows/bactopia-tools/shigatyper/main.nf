#!/usr/bin/env nextflow
/**
 * Rapid determination of Shigella serotypes from sequencing reads.
 *
 * This Bactopia Tool uses [ShigaTyper](https://github.com/CFSAN-Biostatistics/shigatyper) to rapidly
 * determine *Shigella* serotypes from both Illumina (single or paired-end) and Oxford Nanopore
 * reads. ShigaTyper performs k-mer based analysis targeting specific antigenic determinants
 * and marker genes to predict serotypes according to the White-Kauffmann-Le Minor classification
 * scheme. The tool supports multiple sequencing platforms and provides detailed hit statistics
 * for each target gene, enabling rapid serotype identification for epidemiological investigations
 * and outbreak response.
 *
 * @status stable
 * @keywords Shigella, serotyping, k-mer, Illumina, Nanopore, epidemiology, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, shigatyper
 *
 * @subworkflows bactopiatool_init, shigatyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                              Tab-delimited file containing predicted Shigella serotype for each sample
 * @publish *-hits.tsv                         Detailed statistics about each individual gene hit used for serotype prediction
 *
 * @section Merged Results
 * @publish shigatyper.tsv                     Merged TSV file containing ShigaTyper results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                            Tool execution logs including ShigaTyper output
 * @publish logs/nf-*                          Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                       Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { SHIGATYPER        } from '../../../subworkflows/shigatyper/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    SHIGATYPER(BACTOPIATOOL_INIT.out.reads)
    ch_sample_nf_logs = SHIGATYPER.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = SHIGATYPER.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    sample_outputs = SHIGATYPER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = SHIGATYPER.out.run_outputs
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
