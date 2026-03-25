/**
 * Perform capsular typing of Streptococcus pneumoniae from NGS data.
 *
 * This subworkflow uses [PneumoCaT](https://github.com/ukhsa-collaboration/PneumoCaT) to
 * identify serotype-specific capsular loci and determine serotypes from next-generation
 * sequencing data. It provides comprehensive serotype determination including coverage
 * statistics and confidence scores for each sample.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords streptococcus pneumoniae, serotype, capsular typing, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent
 * @citation pneumocat
 *
 * @modules pneumocat
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads (not supported by PneumoCaT)
 * - `lr`: Long reads (not supported by PneumoCaT)
 *
 * @output sample_outputs
 * - `xml`: The PneumoCaT result files in XML format
 * - `txt`: A file containing the coverage information across the genes
 */
nextflow.preview.types = true

include { PNEUMOCAT as PNEUMOCAT_MODULE } from '../../modules/pneumocat/main'

workflow PNEUMOCAT {
    take:
    reads: Channel<Record>

    main:
    PNEUMOCAT_MODULE(reads)

    emit:
    sample_outputs = PNEUMOCAT_MODULE.out
}
