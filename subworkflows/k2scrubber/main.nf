/**
 * Remove human reads from metagenomic data using Kraken2.
 *
 * This subworkflow identifies and removes human sequences from metagenomic reads using [Kraken2](https://github.com/DerrickWood/kraken2)
 * with a specialized human genome database. It downloads the k2_HPRC human reference database,
 * classifies reads, and separates human (contaminant) sequences from microbial reads for downstream analysis.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, human decontamination, read filtering, kraken2
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent, resource-download
 * @citation kraken2
 *
 * @modules wget, kraken2
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @output sample_outputs
 *   - `kraken2_report`: Standard Kraken2 report containing taxonomic abundance counts
 *   - `scrub_report`: Summary report of reads removed during host scrubbing (optional)
 *   - `special_meta`: A simplified metadata map for internal use
 *   - `classified`: Reads assigned to a taxon in the database (FASTQ)
 *   - `unclassified`: Reads NOT assigned to any taxon (FASTQ)
 *   - `classified_extra`: Duplicate classified channel with placeholder for pipeline routing
 *   - `unclassified_extra`: Duplicate unclassified channel with placeholder for pipeline routing
 */
nextflow.preview.types = true

include { WGET    } from '../../modules/wget/main'
include { KRAKEN2 } from '../../modules/kraken2/main'
include { gather  } from 'plugin/nf-bactopia'

workflow K2SCRUBBER {
    take:
    reads: Channel<Record>

    main:
    WGET([
        "name": "k2scrubber",
        "save_as": "k2_HPRC_20230810.tar.gz",
        "url": "https://zenodo.org/records/8339732/files/k2_HPRC_20230810.tar.gz?download=1"
    ])
    KRAKEN2(reads, WGET.out.download)

    emit:
    sample_outputs = KRAKEN2.out
}
