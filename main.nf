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
 * @subworkflows bactopia_init, amrfinderplus, assembler, datasets, gather, sketcher,
 * @subworkflows mlst, qc, bakta, prokka, merlin
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
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    // QC
    adapters              : Path? = "${projectDir}/data/empty/EMPTY_ADAPTERS"    // TODO: Remove when Path? is fixed
    phix                  : Path? = "${projectDir}/data/empty/EMPTY_PHIX"        // TODO: Remove when Path? is fixed

    // Annotation
    use_bakta             : Boolean
    bakta_db              : Path? = "${projectDir}/data/empty/EMPTY_DB"          // TODO: Remove when Path? is fixed
    download_bakta        : Boolean
    bakta_save_as_tarball : Boolean
    bakta_proteins        : Path? = "${projectDir}/data/empty/EMPTY_PROTEINS"    // TODO: Remove when Path? is fixed
    bakta_prodigal_tf     : Path? = "${projectDir}/data/empty/EMPTY_PRODIGAL_TF" // TODO: Remove when Path? is fixed
    bakta_replicons       : Path? = "${projectDir}/data/empty/EMPTY_REPLICONS"   // TODO: Remove when Path? is fixed
    prokka_proteins       : Path
    prokka_prodigal_tf    : Path? = "${projectDir}/data/empty/EMPTY_PRODIGAL_TF" // TODO: Remove when Path? is fixed

    // Merlin
    emmtyper_blastdb      : Path? = "${projectDir}/data/empty/EMPTY_DB"          // TODO: Remove when Path? is fixed
    hicap_database_dir    : Path? = "${projectDir}/data/empty/EMPTY_DB"          // TODO: Remove when Path? is fixed
    hicap_model_fp        : Path? = "${projectDir}/data/empty/EMPTY_PROTEINS"    // TODO: Remove when Path? is fixed
    ask_merlin            : Boolean
    spatyper_repeats      : Path? = "${projectDir}/data/empty/EMPTY_DB"          // TODO: Remove when Path? is fixed
    spatyper_repeat_order : Path? = "${projectDir}/data/empty/EMPTY_PROTEINS"    // TODO: Remove when Path? is fixed
}

// Core
include { BACTOPIA_INIT   } from './subworkflows/utils/bactopia/main'
include { AMRFINDERPLUS   } from './subworkflows/amrfinderplus/main'
include { ASSEMBLER       } from './subworkflows/bactopia/assembler/main'
include { DATASETS        } from './subworkflows/bactopia/datasets/main'
include { GATHER          } from './subworkflows/bactopia/gather/main'
include { SKETCHER        } from './subworkflows/bactopia/sketcher/main'
include { MLST            } from './subworkflows/mlst/main'
include { QC              } from './subworkflows/bactopia/qc/main'

// Annotation wih Bakta or Prokka
include { BAKTA           } from './subworkflows/bakta/main'
include { PROKKA          } from './subworkflows/prokka/main'

// Merlin
include { MERLIN          } from './subworkflows/merlin/main'

workflow {
    main:
    BACTOPIA_INIT()

    // Core steps
    DATASETS()

    // Gather samples in one place
    GATHER(BACTOPIA_INIT.out.samples)

    // QC samples
    QC(GATHER.out.reads, params.adapters, params.phix)

    // Assemble genomes
    ASSEMBLER(QC.out.reads)

    // Sketch and query
    SKETCHER(ASSEMBLER.out.assembly, DATASETS.out.mash_db, DATASETS.out.sourmash_db)

    // Annotate samples
    ch_annotations = channel.empty()
    ch_annotation_sample_outputs = channel.empty()
    if (params.use_bakta) {
        BAKTA(
            ASSEMBLER.out.assembly,
            params.bakta_db,
            params.download_bakta,
            params.bakta_save_as_tarball,
            params.bakta_proteins,
            params.bakta_prodigal_tf,
            params.bakta_replicons
        )
        ch_annotation_sample_outputs = BAKTA.out.sample_outputs
        ch_annotations = BAKTA.out.annotations
    } else {
        PROKKA(
            ASSEMBLER.out.assembly,
            params.prokka_proteins,
            params.prokka_prodigal_tf
        )
        ch_annotation_sample_outputs = PROKKA.out.sample_outputs
        ch_annotations = PROKKA.out.annotations
    }

    // AMR
    AMRFINDERPLUS(ch_annotations, DATASETS.out.amrfinderplus_db)

    // MLST
    MLST(ASSEMBLER.out.assembly, DATASETS.out.mlst_db)

    // Collect all sample-level outputs
    ch_sample_outputs = GATHER.out.sample_outputs
        .mix(QC.out.sample_outputs)
        .mix(ASSEMBLER.out.sample_outputs)
        .mix(SKETCHER.out.sample_outputs)
        .mix(ch_annotation_sample_outputs)
        .mix(AMRFINDERPLUS.out.sample_outputs)
        .mix(MLST.out.sample_outputs)

    // Collect all run-level outputs (only subworkflows that have them)
    ch_run_outputs = GATHER.out.run_outputs
        .mix(ASSEMBLER.out.run_outputs)
        .mix(AMRFINDERPLUS.out.run_outputs)
        .mix(MLST.out.run_outputs)

    // Merlin
    if (params.ask_merlin) {
        MERLIN(
            ASSEMBLER.out.assembly_reads,
            DATASETS.out.mash_db,
            // emmtyper
            params.emmtyper_blastdb,
            // hicap
            params.hicap_database_dir,
            params.hicap_model_fp,
            // staphtyper
            params.spatyper_repeats,
            params.spatyper_repeat_order
        )
        ch_sample_outputs = ch_sample_outputs.mix(MERLIN.out.sample_outputs)
        ch_run_outputs = ch_run_outputs.mix(MERLIN.out.run_outputs)
    }

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = ch_sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = ch_run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = ch_sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = ch_run_outputs
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
