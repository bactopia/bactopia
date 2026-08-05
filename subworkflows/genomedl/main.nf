/**
 * Download genome assemblies from NCBI Datasets.
 *
 * This subworkflow downloads genome assemblies using
 * [genome-dl](https://github.com/rpetit3/genome-dl), which resolves accessions to their latest
 * version and subsamples species queries before fetching files from the NCBI FTP site. The
 * downloaded assemblies are fanned out into per-genome records for downstream analysis, and the
 * first assembly is exposed separately for use as a reference genome.
 *
 * @status stable
 * @keywords download, ncbi, datasets, genome, assembly, refseq
 * @tags complexity:moderate input-type:single output-type:multiple features:resource-download,internet-access
 * @citation genome_dl
 *
 * @modules genomedl as genomedl_module
 *
 * @input accessions
 * A file containing NCBI Assembly accessions, one per line. May be combined with the `--accession` and `--species` parameters.
 *
 * @output sample_outputs
 * - `fna`: Genomic nucleotide sequences in FASTA format
 * - `gbff`: GenBank format genome sequences
 * - `wgs_gbk`: WGS master records in GenBank format
 * - `gff`: Genome annotations in GFF3 format
 * - `gtf`: Genome annotations in GTF format
 * - `faa`: Protein sequences in FASTA format
 * - `gpff`: Protein sequences in GenPept format
 * - `cds`: CDS nucleotide sequences in FASTA format
 * - `translated_cds`: CDS protein sequences in FASTA format
 * - `rna`: RNA feature nucleotide sequences in FASTA format
 * - `features`: Feature table with locations and attributes
 * - `report`: Assembly report with unit and sequence relationships
 * - `stats`: Assembly statistics
 * - `metadata`: NCBI Datasets metadata for each downloaded assembly
 * - `summary`: Human-readable run summary of the version, parameters, and results
 * - `json`: Machine-readable run report of the parameters, results, and per-assembly metadata
 *
 * @output run_outputs
 *
 * @output assemblies
 * - `fna`: Individual downloaded assembly in FASTA format
 *
 * @output reference
 * First downloaded genome for use as a reference. Prefers GenBank (`--format genbank`) over
 * FASTA, since annotation-aware consumers such as Snippy require a GenBank reference.
 */
nextflow.enable.types = true

include { GENOMEDL as GENOMEDL_MODULE } from '../../modules/genomedl/main'

workflow GENOMEDL {

    take:
    accessions: Path?

    main:
    ch_genomedl = GENOMEDL_MODULE(accessions)
    // Fan out on the named `fna` field rather than `results`, which also carries the metadata
    // TSV, run summary, and JSON report that genome-dl always writes
    ch_assemblies = ch_genomedl.map { r -> r.fna }.flatten().map { path ->
        def sample_name = file(path).getSimpleName()
        record(meta: record(id: sample_name, name: sample_name), fna: path)
    }
    // Prefer GenBank over FASTA: consumers of `reference` (Snippy) need the annotations, and
    // `--format` decides which of the two genome-dl actually wrote
    ch_reference = ch_genomedl.map { r -> r.gbff ? r.gbff : r.fna }.flatten().first()

    emit:
    // Downstream inputs
    assemblies = ch_assemblies
    reference = ch_reference
    // Published outputs
    sample_outputs = ch_genomedl
    run_outputs = channel.empty()
}
