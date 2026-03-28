#!/usr/bin/env nextflow
/**
 * Fast alignment-free computation of whole-genome Average Nucleotide Identity.
 *
 * This Bactopia Tool uses [FastANI](https://github.com/ParBLiSS/FastANI) to calculate the average
 * nucleotide identity (ANI) between samples. It can also calculate ANI against reference genomes
 * by downloading RefSeq assemblies using NCBI genome download.
 *
 * @status stable
 * @keywords ani, average nucleotide identity, similarity, comparative genomics, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,comparative
 * @citation fastani
 *
 * @subworkflows bactopiatool_init, fastani, ncbigenomedownload
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input fastani_reference
 * Path to reference FASTA file for ANI comparison
 *
 * @input fastani_pairwise
 * Perform pairwise ANI calculation between all samples
 *
 * @input species
 * Species name to download all RefSeq genomes for comparison
 *
 * @input accession
 * Specific NCBI Assembly RefSeq accession to download
 *
 * @input accessions
 * Path to file containing list of NCBI accessions to download
 *
 * @section Per-Sample Results
 * @publish *.tsv            FastANI results of samples against reference
 *
 * @section Merged Results
 * @publish fastani.tsv       Merged TSV file containing ANI results from all samples
 *
 * @section Execution Logs
 * @publish logs/**   Tool execution logs
 * @publish logs/nf-* Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml Software version information
 *
 * @citation Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S [High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries.](https://dx.doi.org/10.1038/s41467-018-07641-9) _Nat. Commun._ 9, 5114 (2018)
   */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    fastani_reference : Path?
    fastani_pairwise  : Boolean
    species           : String?
    accession         : String?
    accessions        : Path?
}

include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools/main'
include { FASTANI            } from '../../../subworkflows/fastani/main'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'

workflow {
    main:
    BACTOPIATOOL_INIT()

    // Reference if applicable
    ch_reference = channel.empty() as Channel<Record>
    if (params.fastani_reference) {
        ch_reference = ch_reference.mix(
            channel.of(record(_meta: [id: params.fastani_reference.getSimpleName()], fna: params.fastani_reference))
        )
    }

    // Download if applicable
    if (params.species || params.accession || params.accessions) {
        NCBIGENOMEDOWNLOAD(params.accessions)
        ch_reference = ch_reference.mix(
            NCBIGENOMEDOWNLOAD.out.bactopia_tools.map { meta, path ->
                record(_meta: meta, fna: path)
            }
        )
    }

    // Add query if pairwise
    ch_query = BACTOPIATOOL_INIT.out.assembly
    if (params.fastani_pairwise) {
        ch_reference = ch_reference.mix(ch_query)
        ch_query = ch_reference
    }

    // Run FastANI
    FASTANI(ch_query, ch_reference)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = FASTANI.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = FASTANI.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = FASTANI.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = FASTANI.out.run_outputs
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
