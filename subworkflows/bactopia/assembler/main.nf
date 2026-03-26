/**
 * Assemble bacterial genomes using automated assembler selection.
 *
 * This subworkflow automatically selects the optimal assembly strategy based on input read types:
 * - **Short Paired-End Reads:** Uses [Shovill](https://github.com/tseemann/shovill) (SKESA/SPAdes wrapper)
 * - **Short Single-End Reads:** Uses [Shovill-SE](https://github.com/rpetit3/shovill) (SKESA/SPAdes wrapper)
 * - **Long Reads:** Uses [Dragonflye](https://github.com/rpetit3/dragonflye) (Flye/Miniasm wrapper)
 * - **Hybrid Assembly:** Uses [Unicycler](https://github.com/rrwick/Unicycler) or Dragonflye with short-read polishing
 *
 * The workflow performs individual assemblies per sample and aggregates assembly statistics
 * across all samples using [assembly-scan](https://github.com/rpetit3/assembly-scan) for
 * comprehensive quality assessment.
 *
 * @status stable
 * @keywords bacteria, assembly, hybrid, shovill, dragonflye, unicycler, illumina, nanopore
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,alternative-execution
 * @citation any2fasta, assembly_scan, bwa, dragonflye, flash, flye, medaka, megahit, miniasm, minimap2, nanoq, pigz, pilon, racon, rasusa, raven, samclip, samtools, shovill, shovill_se, skesa, spades, unicycler, velvet
 *
 * @modules bactopia/assembler, csvtk/concat
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`  : Illumina R1 reads (paired-end forward)
 * - `r2`  : Illumina R2 reads (paired-end reverse)
 * - `se`  : Single-end Illumina reads
 * - `lr`  : Long reads (ONT/PacBio) for long-read or hybrid assembly
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited report of assembly statistics (N50, length, coverage)
 * - `supplemental`: Supplemental files including assembly graphs and tool-specific logs
 * - `error`: Captured error messages if assembly fails
 *
 * @output run_outputs
 * - `csv`: Aggregated assembly statistics from all samples
 */
nextflow.preview.types = true

include { ASSEMBLER as ASSEMBLER_MODULE } from '../../../modules/bactopia/assembler/main'
include { CSVTK_CONCAT                  } from '../../../modules/csvtk/concat/main'
include { filterWithData                } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow ASSEMBLER {
    take:
    samples: Channel<Record>

    main:
    ASSEMBLER_MODULE(samples)
    CSVTK_CONCAT(gather(ASSEMBLER_MODULE.out, 'tsv', [name: 'assembly-scan']), 'tsv', 'tsv')

    emit:
    // Downstream inputs
    assembly = filterWithData(ASSEMBLER_MODULE.out, ['fna'])
    assembly_reads = filterWithData(ASSEMBLER_MODULE.out, ['fna', 'r1', 'r2', 'se', 'lr'])

    // Published outputs
    sample_outputs = ASSEMBLER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
