/**
 * Assemble bacterial genomes using automated assembler selection.
 *
 * This subworkflow automatically selects the optimal assembly strategy based on input read types:
 * - **Short Paired-End Reads:** Uses [Shovill](https://github.com/tseemann/shovill) (SKESA/SPAdes wrapper)
 * - **Short Single-End Reads:** Uses [Shovill](https://github.com/rpetit3/shovill) (SKESA/SPAdes wrapper)
 * - **Long Reads:** Uses [Dragonflye](https://github.com/rpetit3/dragonflye) (Flye/Miniasm wrapper)
 * - **Hybrid Assembly:** Uses [Unicycler](https://github.com/rrwick/Unicycler) or Dragonflye with short-read polishing
 *
 * The workflow performs individual assemblies per sample and aggregates assembly statistics
 * across all samples using [assembly-scan](https://github.com/rpetit3/assembly-scan) for comprehensive quality assessment.
 *
 * @status stable
 * @keywords bacteria, assembly, hybrid, shovill, dragonflye, unicycler, illumina, nanopore
 * @tags complexity:high input-type:single output-type:multiple features:aggregation, conditional-logic, alternative-execution
 * @citation any2fasta, assembly-scan, bwa, dragonflye, flash, flye, medaka, megahit, miniasm, minimap2, nanoq, pigz, pilon, racon, rasusa, raven, samclip, samtools, shovill, shovill-se, skesa, spades, unicycler, velvet
 *
 * @modules bactopia_assembler, csvtk_concat
 *
 * @input tuple(meta, fq, extra)
 * - `meta`: Groovy Map containing sample information
 * - `fq`: Primary reads (Illumina paired-end or Nanopore)
 * - `extra`: Secondary reads for hybrid assembly or polishing (Optional)
 *
 * @output fna        Assembled contigs in FASTA format
 * @output fna_fq     Tuple containing assembly and primary reads for downstream analysis
 * @output tsv        Per-sample tab-delimited assembly statistics (N50, length, coverage)
 * @output merged_tsv Consolidated assembly statistics report across all samples
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { ASSEMBLER as ASSEMBLER_MODULE } from '../../../modules/bactopia/assembler/main'
include { CSVTK_CONCAT                  } from '../../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow ASSEMBLER {
    take:
    reads: Channel<Tuple<Map, Set<Path>, Set<Path>>>

    main:
    ASSEMBLER_MODULE(reads)
    CSVTK_CONCAT(gather(ASSEMBLER_MODULE.out.tsv, 'assembly-scan'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    fna: Channel<Tuple<Map, Set<Path>>> = ASSEMBLER_MODULE.out.fna
    fna_fq: Channel<Tuple<Map, Set<Path>, Set<Path>>> = ASSEMBLER_MODULE.out.fna_fq
    tsv: Channel<Tuple<Map, Path>> = ASSEMBLER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ASSEMBLER_MODULE.out.fna,
        ASSEMBLER_MODULE.out.tsv,
        ASSEMBLER_MODULE.out.error,
        ASSEMBLER_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ASSEMBLER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ASSEMBLER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        ASSEMBLER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
