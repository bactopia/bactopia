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
 *
 * @subworkflows bactopia_init, amrfinderplus, assembler, datasets, gather, sketcher,
 * @subworkflows mlst, qc, bakta, prokka, staphtyper
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
    rundir   : String

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
include { BACTOPIA_INIT   } from '../../subworkflows/utils/bactopia'
include { AMRFINDERPLUS   } from '../../subworkflows/amrfinderplus/main'
include { ASSEMBLER       } from '../../subworkflows/bactopia/assembler/main'
include { DATASETS        } from '../../modules/bactopia/datasets/main'
include { GATHER          } from '../../subworkflows/bactopia/gather/main'
include { SKETCHER        } from '../../subworkflows/bactopia/sketcher/main'
include { MLST            } from '../../subworkflows/mlst/main'
include { QC              } from '../../subworkflows/bactopia/qc/main'

// Annotation wih Bakta or Prokka
include { BAKTA           } from '../../subworkflows/bakta/main'
include { PROKKA          } from '../../subworkflows/prokka/main'

// Merlin
include { STAPHTYPER      } from '../../subworkflows/staphtyper/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIA_INIT()

    // Core steps
    DATASETS()

    // Gather samples in one place
    GATHER(BACTOPIA_INIT.out.samples)
    ch_results = ch_results.mix(GATHER.out.results)
    ch_logs = ch_logs.mix(GATHER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GATHER.out.nf_logs)
    ch_versions = ch_versions.mix(GATHER.out.versions)

    // QC samples
    QC(
        GATHER.out.raw_fastq,
        params.adapters,
        params.phix
    )
    ch_results = ch_results.mix(QC.out.results)
    ch_logs = ch_logs.mix(QC.out.logs)
    ch_nf_logs = ch_nf_logs.mix(QC.out.nf_logs)
    ch_versions = ch_versions.mix(QC.out.versions)

    // Assemble genomes
    ASSEMBLER(QC.out.fastq)
    ch_results = ch_results.mix(ASSEMBLER.out.results)
    ch_logs = ch_logs.mix(ASSEMBLER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(ASSEMBLER.out.nf_logs)
    ch_versions = ch_versions.mix(ASSEMBLER.out.versions)

    // Sketch and query
    SKETCHER(ASSEMBLER.out.fna, DATASETS.out.mash_db, DATASETS.out.sourmash_db)
    ch_results = ch_results.mix(SKETCHER.out.results)
    ch_logs = ch_logs.mix(SKETCHER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SKETCHER.out.nf_logs)
    ch_versions = ch_versions.mix(SKETCHER.out.versions)

    // Annotate samples
    ch_annotations = channel.empty()
    if (params.use_bakta) {
        BAKTA(
            ASSEMBLER.out.fna,
            params.bakta_db,
            params.download_bakta,
            params.bakta_save_as_tarball,
            params.bakta_proteins,
            params.bakta_prodigal_tf,
            params.bakta_replicons
        )
        ch_results = ch_results.mix(BAKTA.out.results)
        ch_logs = ch_logs.mix(BAKTA.out.logs)
        ch_nf_logs = ch_nf_logs.mix(BAKTA.out.nf_logs)
        ch_versions = ch_versions.mix(BAKTA.out.versions)
        ch_annotations = ch_annotations.mix(BAKTA.out.annotations)
    } else {
        PROKKA(
            ASSEMBLER.out.fna,
            params.prokka_proteins,
            params.prokka_prodigal_tf
        )
        ch_results = ch_results.mix(PROKKA.out.results)
        ch_logs = ch_logs.mix(PROKKA.out.logs)
        ch_nf_logs = ch_nf_logs.mix(PROKKA.out.nf_logs)
        ch_versions = ch_versions.mix(PROKKA.out.versions)
        ch_annotations = ch_annotations.mix(PROKKA.out.annotations)
    }

    // AMR
    AMRFINDERPLUS(ch_annotations, DATASETS.out.amrfinderplus_db)
    ch_results = ch_results.mix(AMRFINDERPLUS.out.results)
    ch_logs = ch_logs.mix(AMRFINDERPLUS.out.logs)
    ch_nf_logs = ch_nf_logs.mix(AMRFINDERPLUS.out.nf_logs)
    ch_versions = ch_versions.mix(AMRFINDERPLUS.out.versions)

    // MLST
    MLST(ASSEMBLER.out.fna, DATASETS.out.mlst_db)
    ch_results = ch_results.mix(MLST.out.results)
    ch_logs = ch_logs.mix(MLST.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MLST.out.nf_logs)
    ch_versions = ch_versions.mix(MLST.out.versions)

    // Staphtyper
    STAPHTYPER(
        ASSEMBLER.out.fna,
        params.spatyper_repeats,
        params.spatyper_repeat_order
    )
    ch_results = ch_results.mix(STAPHTYPER.out.results)
    ch_logs = ch_logs.mix(STAPHTYPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(STAPHTYPER.out.nf_logs)
    ch_versions = ch_versions.mix(STAPHTYPER.out.versions)

    // Branch the based on scope (sample or run)
    ch_final_results = ch_results.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_logs = ch_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_nf_logs = ch_nf_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_versions = ch_versions.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    publish:
    run_results = ch_final_results.run
    run_logs = ch_final_logs.run
    run_nf_logs = ch_final_nf_logs.run
    run_versions = ch_final_versions.run
    sample_results = ch_final_results.sample
    sample_logs = ch_final_logs.sample
    sample_nf_logs = ch_final_nf_logs.sample
    sample_versions = ch_final_versions.sample
}

output {
    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        }
    }
    run_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${meta.logs_dir}/nf${file.name}"
        }
    }
    sample_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}
