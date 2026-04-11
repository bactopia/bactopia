#!/usr/bin/env nextflow
/**
 * Comprehensive analysis pipeline for Staphylococcus aureus isolates.
 *
 * This workflow performs complete bacterial analysis including quality control,
 * assembly, annotation, antimicrobial resistance detection, MLST typing,
 * and Staphylococcus-specific analysis using [Spatyper](https://github.com/HCGB-IGTP/spaTyper),
 * [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE), and [SCCmecFinder](https://github.com/rpetit3/sccmec).
 * It processes raw sequencing reads and produces a comprehensive genomic characterization for S. aureus isolates.
 *
 * @status stable
 * @keywords Staphylococcus aureus, assembly, annotation, AMR, MLST, spa typing, agr typing, sccmec
 * @tags complexity:complex input-type:parameter output-type:multiple features:aggregation,conditional-logic,database-dependent
 * @citation staphopia
 *
 * @subworkflows utils_bactopia, amrfinderplus, bactopia_assembler, bactopia_datasets,
 *               bactopia_gather, bactopia_qc, bactopia_sketcher, bakta, mlst, prokka, staphtyper
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
 * Path to trusted protein sequences for Bakta
 *
 * @input bakta_prodigal_tf
 * Path to Prodigal training file for Bakta
 *
 * @input bakta_replicons
 * Path to replicon sequences for Bakta
 *
 * @input prokka_proteins
 * Path to protein sequences for Prokka
 *
 * @input prokka_prodigal_tf
 * Path to Prodigal training file for Prokka
 *
 * @input spatyper_repeats
 * Path to repeats database for Spatyper
 *
 * @input spatyper_repeat_order
 * Path to repeat order file for Spatyper
 *
 * @section Quality Control
 * @publish supplemental/*_fastqc.*         FastQC quality control reports for raw and cleaned reads
 * @publish supplemental/*-NanoPlot.*      NanoPlot reports for Nanopore reads
 * @publish supplemental/*.fastp.*         Fastp quality reports (when applicable)
 *
 * @section Assembly
 * @publish *.fna                         Assembled genome sequences in FASTA format
 * @publish assembly-stats.tsv            Assembly quality metrics per sample
 *
 * @section Annotation
 * @note Output format depends on chosen annotation tool (Bakta or Prokka)
 * @publish *.gff.gz                      Genome annotation in GFF3 format (compressed)
 * @publish *.gbk.gz                      Genome annotation in GenBank format (compressed)
 * @publish *.faa.gz                      Protein sequences (compressed)
 * @publish *.fna.gz                      Nucleotide sequences from annotation (compressed)
 * @publish annotation.tsv                Annotation summary tables
 *
 * @section Typing
 * @publish mlst.tsv                      MLST sequence type results
 * @publish agrvate-*                     Agr locus typing results
 * @publish spatyper-*                    spa typing results
 * @publish sccmec-*                      SCCmec typing results (targets, regions, details)
 *
 * @section Antimicrobial Resistance
 * @publish amrfinderplus.tsv             AMR gene detection results
 * @publish amrfinderplus.mutation.tsv    AMR point mutation results
 *
 * @section Comparative Analysis
 * @publish *-k21.msh                     Mash sketch files (k=21)
 * @publish *-k31.msh                     Mash sketch files (k=31)
 * @publish *-mash-refseq88-*.txt         Mash screening results against RefSeq
 * @publish *.sig                         Sourmash signatures
 * @publish sourmash-*.txt                Sourmash classification results
 *
 * @section Merged Results
 * @note Run-level aggregated results from all samples
 * @publish merged-assembly-stats.tsv     Consolidated assembly statistics
 * @publish merged-mlst.tsv               Consolidated MLST results
 * @publish staphtyper.tsv                Consolidated Staphylococcus typing summary
 *
 * @section Execution Logs
 * @publish logs/**                       Tool execution logs
 * @publish logs/nf-*                     Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                  Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    adapters              : Path?
    phix                  : Path?
    use_bakta             : Boolean
    bakta_db              : Path?
    download_bakta        : Boolean
    bakta_save_as_tarball : Boolean
    bakta_proteins        : Path?
    bakta_prodigal_tf     : Path?
    bakta_replicons       : Path?
    prokka_proteins       : Path?
    prokka_prodigal_tf    : Path?
    spatyper_repeats      : Path?
    spatyper_repeat_order : Path?
}

// Core
include { BACTOPIA_INIT       } from '../../subworkflows/utils/bactopia'
include { AMRFINDERPLUS       } from '../../subworkflows/amrfinderplus/main'
include { ASSEMBLER           } from '../../subworkflows/bactopia/assembler/main'
include { DATASETS            } from '../../subworkflows/bactopia/datasets/main'
include { GATHER              } from '../../subworkflows/bactopia/gather/main'
include { SKETCHER            } from '../../subworkflows/bactopia/sketcher/main'
include { MLST                } from '../../subworkflows/mlst/main'
include { QC                  } from '../../subworkflows/bactopia/qc/main'

// Annotation wih Bakta or Prokka
include { BAKTA               } from '../../subworkflows/bakta/main'
include { PROKKA              } from '../../subworkflows/prokka/main'

// Merlin
include { STAPHTYPER          } from '../../subworkflows/staphtyper/main'
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

    // Staphtyper
    ch_staphtyper = STAPHTYPER(
        ch_assembler.assembly,
        params.spatyper_repeats,
        params.spatyper_repeat_order
    )

    // Collect all outputs
    ch_sample_outputs = ch_gather.sample_outputs
        .mix(ch_qc.sample_outputs)
        .mix(ch_assembler.sample_outputs)
        .mix(ch_sketcher.sample_outputs)
        .mix(ch_annotation_outputs.sample_outputs)
        .mix(ch_amrfinderplus.sample_outputs)
        .mix(ch_mlst.sample_outputs)
        .mix(ch_staphtyper.sample_outputs)

    ch_run_outputs = ch_gather.run_outputs
        .mix(ch_qc.run_outputs)
        .mix(ch_assembler.run_outputs)
        .mix(ch_sketcher.run_outputs)
        .mix(ch_annotation_outputs.run_outputs)
        .mix(ch_amrfinderplus.run_outputs)
        .mix(ch_mlst.run_outputs)
        .mix(ch_staphtyper.run_outputs)

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
