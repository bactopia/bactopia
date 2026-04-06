/**
 * Download bacterial genomes from NCBI's RefSeq database.
 *
 * This subworkflow downloads complete and draft bacterial genomes using the
 * [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) tool. It fetches
 * genome assemblies in various formats including GenBank, GFF, and FASTA files
 * along with associated annotation files and statistics.
 *
 * @status stable
 * @keywords download, ncbi, refseq, genome, assembly, database
 * @tags complexity:moderate input-type:single output-type:multiple features:resource-download,database-dependent
 * @citation ncbigenomedownload
 *
 * @modules ncbigenomedownload as ncbigenomedownload_module
 *
 * @input accessions
 * A file containing NCBI accession numbers, one per line. If empty, will download all genomes matching the specified criteria.
 *
 * @output sample_outputs
 *   - `gbff`: GenBank format genome sequences
 *   - `fna`: Genomic nucleotide sequences in FASTA format
 *   - `gff`: Genome annotations in GFF3 format
 *   - `faa`: Protein sequences in FASTA format
 *   - `gpff`: Protein sequences in GenPept format
 *   - `wgs_gbk`: WGS master records in GenBank format
 *   - `cds`: CDS nucleotide sequences in FASTA format
 *   - `rna`: RNA product sequences in FASTA format
 *   - `rna_fna`: RNA feature nucleotide sequences in FASTA format
 *   - `features`: Feature table with locations and attributes
 *   - `rm`: RepeatMasker output (optional)
 *   - `report`: Assembly report with unit and sequence relationships
 *   - `stats`: Assembly statistics
 *   - `accessions`: Generated accession list files
 * @output bactopia_tools  Downloaded files formatted for Bactopia Tools workflows
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { NCBIGENOMEDOWNLOAD as NCBIGENOMEDOWNLOAD_MODULE } from '../../modules/ncbigenomedownload/main'

workflow NCBIGENOMEDOWNLOAD {

    take:
    accessions: Value<Path?>

    main:
    NCBIGENOMEDOWNLOAD_MODULE(accessions)
    ch_assemblies = NCBIGENOMEDOWNLOAD_MODULE.out.map { r -> r.results }.flatten().map { path ->
        record(meta: [id: file(path).getSimpleName()], fna: path)
    }
    ch_reference = NCBIGENOMEDOWNLOAD_MODULE.out.map { r -> r.results }.flatten().first()

    emit: // bactopia-lint: ignore S005, S010
    // Downstream inputs
    assemblies = ch_assemblies
    reference = ch_reference
    // Published outputs
    sample_outputs = NCBIGENOMEDOWNLOAD_MODULE.out
    run_outputs = channel.empty()
}
