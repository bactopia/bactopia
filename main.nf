#!/usr/bin/env nextflow
/**
 * Comprehensive bacterial analysis pipeline for complete genomic characterization.
 *
 * This workflow performs end-to-end analysis including quality control, assembly,
 * annotation, antimicrobial resistance detection, MLST typing, and optional
 * pathogen-specific analysis through Merlin. It processes raw sequencing reads
 * and produces a complete genomic characterization suitable for downstream analysis.
 *
 * @status stable
 * @keywords bacteria, assembly, annotation, AMR, MLST, genomics, pipeline
 * @tags complexity:complex input-type:parameter output-type:multiple features:aggregation,conditional-logic,database-dependent
 *
 * @subworkflows utils_bactopia, amrfinderplus, bactopia_assembler, bactopia_datasets,
 *               bactopia_gather, bactopia_qc, bactopia_sketcher, bakta, merlin, mlst, prokka
 *
 * @input rundir
 * Directory containing raw sequencing reads
 *
 * @input adapters
 * Path to adapter sequences file for removal during QC
 *
 * @input phix
 * Path to PhiX sequences for contamination removal during QC
 *
 * @input use_bakta
 * Use Bakta for genome annotation instead of Prokka
 *
 * @input bakta_db
 * Path to Bakta database for annotation
 *
 * @input download_bakta
 * Download Bakta database if not provided
 *
 * @input bakta_save_as_tarball
 * Save Bakta database as tarball for reuse
 *
 * @input bakta_proteins
 * Path to trusted protein sequences for Bakta annotation
 *
 * @input bakta_prodigal_tf
 * Path to Prodigal training file for Bakta
 *
 * @input bakta_replicons
 * Path to replicon sequences for Bakta
 *
 * @input prokka_proteins
 * Path to protein sequences for Prokka annotation
 *
 * @input prokka_prodigal_tf
 * Path to Prodigal training file for Prokka
 *
 * @input emmtyper_blastdb
 * Path to emmtyper BLAST database for Merlin
 *
 * @input hicap_database_dir
 * Path to HiCap database directory for Merlin
 *
 * @input hicap_model_fp
 * Path to HiCap model file for Merlin
 *
 * @input ask_merlin
 * Enable Merlin pathogen-specific analysis
 *
 * @input spatyper_repeats
 * Path to Spatyper repeats database for Merlin
 *
 * @input spatyper_repeat_order
 * Path to Spatyper repeat order file for Merlin
 *
 * @section Quality Control
 * @publish supplemental/*_fastqc.*          FastQC quality control reports for raw and cleaned reads
 * @publish supplemental/*-NanoPlot.*       NanoPlot reports for Nanopore reads
 * @publish supplemental/*.fastp.*          Fastp quality reports (when applicable)
 * @publish supplemental/*_original.json    Quality metrics for original reads
 * @publish supplemental/*_final.json       Quality metrics for final reads
 *
 * @section Assembly
 * @publish *.fasta                           Assembled genome sequences
 * @publish assembly-stats.tsv              Assembly quality metrics
 * @publish merged-assembly-stats.tsv        Consolidated assembly statistics
 *
 * @section Annotation
 * @note Output format depends on chosen annotation tool (Bakta or Prokka)
 * @publish *.gff.gz                         Genome annotation in GFF3 format (compressed)
 * @publish *.gbk.gz                         Genome annotation in GenBank format (compressed)
 * @publish *.faa.gz                         Protein sequences (compressed)
 * @publish *.fna.gz                         Nucleotide sequences from annotation (compressed)
 * @publish *.ffn.gz                         Feature nucleotide sequences (compressed)
 * @publish annotation.tsv                   Annotation summary tables
 * @publish blastdb.*                        BLAST database created from annotation
 *
 * @section Typing
 * @publish mlst.tsv                          MLST sequence type results
 * @publish merged-mlst.tsv                   Consolidated MLST results
 *
 * @section Antimicrobial Resistance
 * @publish amrfinderplus.tsv                AMR gene detection results
 * @publish amrfinderplus.mutation.tsv       AMR point mutation results
 * @publish merged-amrfinderplus.tsv         Consolidated AMR results
 *
 * @section Comparative Analysis
 * @publish *-k21.msh                        Mash sketch files (k=21)
 * @publish *-k31.msh                        Mash sketch files (k=31)
 * @publish *-mash-refseq88-*.txt            Mash screening results against RefSeq
 * @publish *.sig                            Sourmash signatures
 * @publish sourmash-*.txt                   Sourmash classification results
 *
 * @section Pathogen-Specific Analysis
 * @note Only created if --ask_merlin is enabled
 * @publish merlin/clermontyping/*           E. coli phylogroup typing
 * @publish merlin/ectyper/*                 Enterotoxigenic E. coli typing
 * @publish merlin/shigatyper/*              Shigella serotype prediction
 * @publish merlin/shigapass/*               Shigella passive surveillance
 * @publish merlin/shigeifinder/*            Shigella and EIEC detection
 * @publish merlin/stecfinder/*              STEC detection and typing
 * @publish merlin/emmtyper/*                S. pyogenes emm typing
 * @publish merlin/hicap/*                   H. influenzae capsular typing
 * @publish merlin/hpsuissero/*              H. parasuis serotyping
 * @publish merlin/kleborate/*               Klebsiella species typing
 * @publish merlin/staphtyper/*              S. aureus spa typing
 * @publish merlin/agrvate/*                 S. aureus agr typing
 * @publish merlin/sccmec/*                  S. aureus SCCmec typing
 *
 * @section Merged Results
 * @note Run-level aggregated results from all samples
 * @publish samplesheet.tsv                  Sample metadata and quality metrics
 *
 * @section Execution Logs
 * @publish logs/**                          Tool execution logs
 * @publish logs/nf-*                        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                     Software version information
 */
nextflow.enable.types = true

params {
    rundir : String

    // Tool-specific parameters
    // QC
    adapters              : Path?
    phix                  : Path?

    // Annotation
    use_bakta             : Boolean
    bakta_db              : Path?
    download_bakta        : Boolean
    bakta_save_as_tarball : Boolean
    bakta_proteins        : Path?
    bakta_prodigal_tf     : Path?
    bakta_replicons       : Path?
    prokka_proteins       : Path?
    prokka_prodigal_tf    : Path?

    // Merlin
    emmtyper_blastdb      : Path?
    hicap_database_dir    : Path?
    hicap_model_fp        : Path?
    ask_merlin            : Boolean
    spatyper_repeats      : Path?
    spatyper_repeat_order : Path?
    staphscan_db_mlst     : Path?
}

// Core
include { BACTOPIA_INIT       } from './subworkflows/utils/bactopia/main'
include { AMRFINDERPLUS       } from './subworkflows/amrfinderplus/main'
include { ASSEMBLER           } from './subworkflows/bactopia/assembler/main'
include { DATASETS            } from './subworkflows/bactopia/datasets/main'
include { GATHER              } from './subworkflows/bactopia/gather/main'
include { SKETCHER            } from './subworkflows/bactopia/sketcher/main'
include { MLST                } from './subworkflows/mlst/main'
include { QC                  } from './subworkflows/bactopia/qc/main'

// Annotation wih Bakta or Prokka
include { BAKTA               } from './subworkflows/bakta/main'
include { PROKKA              } from './subworkflows/prokka/main'

// Merlin
include { MERLIN              } from './subworkflows/merlin/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_samples = BACTOPIA_INIT()

    // Core steps
    ch_datasets = DATASETS()

    // Gather samples in one place
    ch_gather = GATHER(ch_samples)

    // QC samples
    ch_qc = QC(ch_gather.reads, params.adapters, params.phix)

    // Assemble genomes
    ch_assembler = ASSEMBLER(ch_qc.reads)

    // Sketch and query
    ch_sketcher = SKETCHER(ch_assembler.assembly, ch_datasets.mash_db, ch_datasets.sourmash_db)

    // Annotate samples
    ch_annotations = channel.empty()
    ch_annotation_outputs = null
    if (params.use_bakta) {
        ch_bakta = BAKTA(
            ch_assembler.assembly,
            params.bakta_db,
            params.download_bakta,
            params.bakta_save_as_tarball,
            params.bakta_proteins,
            params.bakta_prodigal_tf,
            params.bakta_replicons
        )
        ch_annotation_outputs = ch_bakta
        ch_annotations = ch_bakta.annotations
    } else {
        ch_prokka = PROKKA(
            ch_assembler.assembly,
            params.prokka_proteins,
            params.prokka_prodigal_tf
        )
        ch_annotation_outputs = ch_prokka
        ch_annotations = ch_prokka.annotations
    }

    // AMR
    ch_amrfinderplus = AMRFINDERPLUS(ch_annotations, ch_datasets.amrfinderplus_db)

    // MLST
    ch_mlst = MLST(ch_assembler.assembly, ch_datasets.mlst_db)

    // Collect all outputs
    ch_sample_outputs = ch_gather.sample_outputs
        .mix(ch_qc.sample_outputs)
        .mix(ch_assembler.sample_outputs)
        .mix(ch_sketcher.sample_outputs)
        .mix(ch_annotation_outputs.sample_outputs)
        .mix(ch_amrfinderplus.sample_outputs)
        .mix(ch_mlst.sample_outputs)

    ch_run_outputs = ch_gather.run_outputs
        .mix(ch_qc.run_outputs)
        .mix(ch_assembler.run_outputs)
        .mix(ch_sketcher.run_outputs)
        .mix(ch_annotation_outputs.run_outputs)
        .mix(ch_amrfinderplus.run_outputs)
        .mix(ch_mlst.run_outputs)

    // Merlin
    if (params.ask_merlin) {
        ch_merlin = MERLIN(
            ch_assembler.assembly_reads,
            ch_datasets.mash_db,
            // emmtyper
            params.emmtyper_blastdb,
            // hicap
            params.hicap_database_dir,
            params.hicap_model_fp,
            // staphtyper
            params.spatyper_repeats,
            params.spatyper_repeat_order,
            params.staphscan_db_mlst
        )
        ch_sample_outputs = ch_sample_outputs.mix(ch_merlin.sample_outputs)
        ch_run_outputs = ch_run_outputs.mix(ch_merlin.run_outputs)
    }

    publish:
    // Per-sample
    sample_outputs = ch_sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_sample_outputs)
    // Run-level
    run_outputs = ch_run_outputs
    run_nf_logs = collectNextflowLogs(ch_run_outputs)
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
