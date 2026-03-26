/**
 * Determine the agr, spa and SCCmec types for _Staphylococcus aureus_ genomes.
 *
 * This subworkflow performs comprehensive typing of *Staphylococcus aureus* genomes by
 * determining the agr locus type using [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE),
 * spa repeat type using [spaTyper](https://github.com/HCGB-IGTP/spaTyper), and SCCmec element
 * type using SCCmec typing. It combines results from multiple typing methods to provide
 * a complete characterization of *S. aureus* strains.
 *
 * @status stable
 * @keywords staphylococcus aureus, agr typing, spa typing, sccmec, strain characterization
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation, database-dependent
 * @citation agrvate, spatyper, sccmec
 *
 * @subworkflows agrvate, spatyper, sccmec
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input repeats
 * Optional spa repeats database for improved spa typing
 *
 * @input repeat_order
 * Optional spa repeat order file for improved spa typing
 *
 * @output sample_outputs
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { AGRVATE  } from '../agrvate/main'
include { SPATYPER } from '../spatyper/main'
include { SCCMEC   } from '../sccmec/main'

workflow STAPHTYPER {
    take:
    assembly: Channel<Record>
    repeats: Path?
    repeat_order: Path?

    main:
    // agrvate - agr locus type and agr operon variants
    AGRVATE(assembly)

    // spatyper - spa typing
    SPATYPER(assembly, repeats, repeat_order)

    // sccmec - SCCmec type based on targets and full cassettes
    SCCMEC(assembly)

    emit:
    // Published outputs
    sample_outputs = AGRVATE.out.sample_outputs.mix(SPATYPER.out.sample_outputs, SCCMEC.out.sample_outputs)
    run_outputs = AGRVATE.out.run_outputs.mix(SPATYPER.out.run_outputs, SCCMEC.out.run_outputs)
}
