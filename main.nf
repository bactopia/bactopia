#! /usr/bin/env nextflow
import groovy.json.JsonSlurper
import groovy.text.SimpleTemplateEngine
import groovy.util.FileNameByRegexFinder
import java.io.RandomAccessFile;
import java.util.zip.GZIPInputStream;
import java.nio.file.Path
import java.nio.file.Paths
import nextflow.util.SysHelper
PROGRAM_NAME = workflow.manifest.name
VERSION = workflow.manifest.version

// Adjust memory/cpu requests for standard profile only
MAX_MEMORY = ['standard', 'docker', 'singularity'].contains(workflow.profile) ? get_max_memory(params.max_memory).GB : (params.max_memory).GB
MAX_MEMORY_INT = MAX_MEMORY.toString().split(" ")[0]
MAX_CPUS = ['standard', 'docker', 'singularity'].contains(workflow.profile) ? get_max_cpus(params.cpus.toInteger()) : params.cpus.toInteger()
MAX_CPUS_75 = Math.round(MAX_CPUS * 0.75)
MAX_CPUS_50 = Math.round(MAX_CPUS * 0.50)
USING_MERGE = false

// Validate parameters
if (params.help || params.help_all || params.conda_help) print_usage();
if (params.nfdir) print_basedir();
if (params.nfdir) print_basedir();
if (params.available_datasets && params.datasets) print_available_datasets(params.datasets);
if (workflow.commandLine.trim().endsWith(workflow.scriptName)) print_usage();
if (params.example_fastqs) print_example_fastqs();
if (params.version) print_version();
run_type = check_input_params()
check_input_fastqs(run_type)
if (params.check_fastqs) print_check_fastqs(run_type);

// Setup output directories
outdir = params.outdir ? params.outdir : './'

// Setup some defaults
log.info "${PROGRAM_NAME} - ${VERSION}"

AMR_DATABASES = []
ARIBA_DATABASES = []
MINMER_DATABASES = []
MLST_DATABASES = []
REFERENCES = []
PRIMERS = []
BLAST_GENE_FASTAS = []
BLAST_PRIMER_FASTAS = []
BLAST_PROTEIN_FASTAS = []
MAPPING_FASTAS = []
PROKKA_PROTEINS = file(params.empty_proteins)
PRODIGAL_TF = file(params.empty_tf)
REFSEQ_SKETCH = []
REFSEQ_SKETCH_FOUND = false
SPECIES = format_species(params.species)
SPECIES_GENOME_SIZE = null
print_efficiency()
setup_datasets()


process gather_fastqs {
    /* Gather up input FASTQs for analysis. */
    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "bactopia.versions"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    tag "${sample}"

    input:
    set val(sample), val(sample_type), val(single_end), file(r1: '*???-r1'), file(r2: '*???-r2'), file(extra) from create_input_channel(run_type)

    output:
    file "*-error.txt" optional true
    set val(sample), val(final_sample_type), val(single_end),
        file("fastqs/${sample}*.fastq.gz"), file("extra/*.gz") optional true into FASTQ_PE_STATUS
    file "${task.process}/*" optional true
    file "bactopia.versions" optional true
    file "multiple-read-sets-merged.txt" optional true

    shell:
    bactopia_version = VERSION
    nextflow_version = nextflow.version
    is_assembly = sample_type.startsWith('assembly') ? true : false
    is_compressed = false
    no_cache = params.no_cache ? '-N' : ''
    use_ena = params.use_ena
    if (task.attempt >= 4) {
        if (use_ena) {
            // Try SRA
            use_ena = false
        } else {
            // Try ENA
            use_ena = true
        }
    }
    if (extra) {
        is_compressed = extra.getName().endsWith('gz') ? true : false
    }
    section = null
    if (sample_type == 'assembly_accession') {
        section = sample.startsWith('GCF') ? 'refseq' : 'genbank'
    }
    fcov = params.coverage.toInteger() == 0 ? 150 : Math.round(params.coverage.toInteger() * 1.5)
    final_sample_type = sample_type
    if (sample_type == 'hybrid-merge-pe') {
        final_sample_type = 'hybrid'
    } else if (sample_type == 'merge-pe') {
        final_sample_type = 'paired-end'
    } else if (sample_type == 'merge-se') {
        final_sample_type = 'single-end'
    }
    template(task.ext.template)
}

process fastq_status {
    /* Determine if FASTQs are PE or SE, and if they meet minimum basepair/read counts. */
    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    input:
    set val(sample), val(sample_type), val(single_end), file(fq), file(extra) from FASTQ_PE_STATUS

    output:
    file "*-error.txt" optional true
    set val(sample), val(sample_type), val(single_end),
        file("fastqs/${sample}*.fastq.gz"), file(extra) optional true into ESTIMATE_GENOME_SIZE
    file "${task.process}/*" optional true

    shell:
    single_end = fq[1] == null ? true : false
    qin = sample_type.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    template(task.ext.template)
}

process estimate_genome_size {
    /* Estimate the input genome size if not given. */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    input:
    set val(sample), val(sample_type), val(single_end), file(fq), file(extra) from ESTIMATE_GENOME_SIZE

    output:
    file "${sample}-genome-size-error.txt" optional true
    file("${sample}-genome-size.txt") optional true
    set val(sample), val(sample_type), val(single_end),
        file("fastqs/${sample}*.fastq.gz"), file(extra), file("${sample}-genome-size.txt") optional true into QC_READS, QC_ORIGINAL_SUMMARY
    file "${task.process}/*" optional true

    shell:
    genome_size = SPECIES_GENOME_SIZE
    template(task.ext.template)
}

process qc_reads {
    /* Cleanup the reads using Illumina-Cleanup */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*error.txt"

    input:
    set val(sample), val(sample_type), val(single_end), file(fq), file(extra), file(genome_size) from QC_READS

    output:
    file "*-error.txt" optional true
    file "quality-control/*"
    set val(sample), val(single_end),
        file("quality-control/${sample}*.fastq.gz") optional true into COUNT_31MERS, ARIBA_ANALYSIS,
                                                                       MINMER_SKETCH, CALL_VARIANTS,
                                                                       MAPPING_QUERY
    set val(sample), val(sample_type), val(single_end),
        file("quality-control/${sample}*.fastq.gz"), file(extra),
        file(genome_size) optional true into ASSEMBLY

    set val(sample), val(single_end),
        file("quality-control/${sample}*.{fastq,error-fq}.gz"),
        file(genome_size) optional true into QC_FINAL_SUMMARY
    file "${task.process}/*" optional true

    shell:
    qc_ram = task.memory.toString().split(' ')[0]
    is_assembly = sample_type.startsWith('assembly') ? true : false
    qin = sample_type.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    adapters = params.adapters ? file(params.adapters) : 'adapters'
    phix = params.phix ? file(params.phix) : 'phix'
    template(task.ext.template)
}

process qc_original_summary {
    /* Run FASTQC on the input FASTQ files. */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"

    input:
    set val(sample), val(sample_type), val(single_end), file(fq), file(extra), file(genome_size) from QC_ORIGINAL_SUMMARY

    output:
    file "quality-control/*"
    file "${task.process}/*" optional true

    shell:
    template(task.ext.template)
}

process qc_final_summary {
    /* Run FASTQC on the cleaned up FASTQ files. */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "quality-control/*"

    input:
    set val(sample), val(single_end), file(fq), file(genome_size) from QC_FINAL_SUMMARY

    output:
    file "quality-control/*"
    file "${task.process}/*" optional true

    shell:
    template(task.ext.template)
}


process assemble_genome {
    /* Assemble the genome using Shovill, SKESA is used by default */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "assembly/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${sample}-assembly-error.txt"

    input:
    set val(sample), val(sample_type), val(single_end), file(fq), file(extra), file(genome_size) from ASSEMBLY

    output:
    file "assembly/*"
    file "${sample}-assembly-error.txt" optional true
    set val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("assembly/${sample}.{fna,fna.gz}") optional true into SEQUENCE_TYPE
    set val(sample), val(single_end), file("assembly/${sample}.{fna,fna.gz}") optional true into MAKE_BLASTDB
    set val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("assembly/${sample}.{fna,fna.gz}"), file("total_contigs_*") optional true into ANNOTATION
    set val(sample), file("assembly/${sample}.{fna,fna.gz}"), file(genome_size) optional true into ASSEMBLY_QC
    file "${task.process}/*" optional true

    shell:
    contig_namefmt = params.contig_namefmt ? params.contig_namefmt : "${sample}_%05d"
    shovill_ram = task.memory.toString().split(' ')[0]
    opts = params.shovill_opts ? "--opts '${params.shovill_opts}'" : ""
    kmers = params.shovill_kmers ? "--kmers '${params.shovill_kmers}'" : ""
    nostitch = params.nostitch ? "--nostitch" : ""
    nocorr = params.nocorr ? "--nocorr" : ""
    no_miniasm = params.no_miniasm ? "--no_miniasm" : ""
    no_rotate = params.no_rotate ? "--no_rotate" : ""
    no_pilon = params.no_pilon ? "--no_pilon" : ""
    keep = params.keep_all_files ? "--keep 3" : "--keep 1"
    use_original_assembly = null
    if (sample_type.startsWith('assembly')) {
        use_original_assembly = params.reassemble ? false : true
    }
    template(task.ext.template)
}

process assembly_qc {
    /* Assess the quality of the assembly using QUAST and CheckM */
    tag "${sample} - ${method}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/assembly", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${method}/*"

    input:
    set val(sample), file(fasta), file(genome_size) from ASSEMBLY_QC
    each method from Channel.fromList(['checkm', 'quast'])

    output:
    file "${method}/*"
    file "${task.process}/*" optional true

    shell:
    //CheckM Related
    full_tree = params.full_tree ? '' : '--reduced_tree'
    checkm_ali = params.checkm_ali ? '--ali' : ''
    checkm_nt = params.checkm_nt ? '--nt' : ''
    force_domain = params.force_domain ? '--force_domain' : ''
    no_refinement = params.no_refinement ? '--no_refinement' : ''
    individual_markers = params.individual_markers ? '--individual_markers' : ''
    skip_adj_correction = params.skip_adj_correction ? '--skip_adj_correction' : ''
    skip_pseudogene_correction = params.skip_pseudogene_correction ? '--skip_pseudogene_correction' : ''
    ignore_thresholds = params.ignore_thresholds ? '--ignore_thresholds' : ''
    template(task.ext.template)
}


process make_blastdb {
    /* Create a BLAST database of the assembly using BLAST */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "blastdb/*"

    input:
    set val(sample), val(single_end), file(fasta) from MAKE_BLASTDB

    output:
    file("blastdb/*")
    set val(sample), file("blastdb/*") into BLAST_GENES, BLAST_PRIMERS, BLAST_PROTEINS
    file "${task.process}/*" optional true

    shell:
    template(task.ext.template)
}


process annotate_genome {
    /* Annotate the assembly using Prokka, use a proteins FASTA if available */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "annotation/${sample}*"

    input:
    set val(sample), val(single_end), file(fq), file(fasta), file(total_contigs) from ANNOTATION
    file prokka_proteins from PROKKA_PROTEINS
    file prodigal_tf from PRODIGAL_TF

    output:
    file "annotation/${sample}*"
    set val(sample), file("annotation/${sample}.{ffn,ffn.gz}") optional true into PLASMID_BLAST
    set val(sample),
        file("annotation/${sample}.{ffn,ffn.gz}"),
        file("annotation/${sample}.{faa,faa.gz}") optional true into ANTIMICROBIAL_RESISTANCE
    file "${task.process}/*" optional true

    shell:
    gunzip_fasta = fasta.getName().replace('.gz', '')
    contig_count = total_contigs.getName().replace('total_contigs_', '')
    genus = "Genus"
    species = "species"
    proteins = ""
    if (prokka_proteins.getName() != 'EMPTY_PROTEINS') {
        proteins = "--proteins ${prokka_proteins}"
        if (SPECIES.contains("-")) {
            genus = SPECIES.split('-')[0].capitalize()
            species = SPECIES.split('-')[1]
        } else {
            genus = SPECIES.capitalize()
            species = "spp."
        }
    }

    prodigal = ""
    if (prodigal_tf.getName() != 'EMPTY_TF' && !params.skip_prodigal_tf) {
        prodigal = "--prodigaltf ${prodigal_tf}"
    }

    compliant = params.compliant ? "--compliant" : ""
    locustag = "--locustag ${sample}"
    renamed = false
    // Contig ID must <= 37 characters
    if ("gnl|${params.centre}|${sample}_${contig_count}".length() > 37) {
        locustag = ""
        compliant = "--compliant"
        renamed = true
    }
    addgenes = params.nogenes ? "" : "--addgenes"
    addmrna = params.addmrna ? "--addmrna" : ""
    rawproduct = params.rawproduct ? "--rawproduct" : ""
    cdsrnaolap = params.cdsrnaolap ? "--cdsrnaolap" : ""
    norrna = params.norrna ? "--norrna" : ""
    notrna = params.notrna ? "--notrna" : ""
    rnammer = params.rnammer ? "--rnammer" : ""
    rfam = params.rnammer ? "--rfam" : ""
    template(task.ext.template)
}


process count_31mers {
    /* Count 31mers in the reads using McCortex */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/kmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.ctx"

    input:
    set val(sample), val(single_end), file(fq) from COUNT_31MERS

    output:
    file "${sample}.ctx"
    file "${task.process}/*" optional true

    shell:
    m = task.memory.toString().split(' ')[0].toInteger() * 1000 - 500
    template(task.ext.template)
}


process sequence_type {
    /* Determine MLST types using ARIBA and BLAST */
    tag "${sample} - ${schema} - ${method}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/mlst/${schema}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${method}/*"

    input:
    set val(sample), val(single_end), file(fq), file(assembly) from SEQUENCE_TYPE
    each file(dataset) from MLST_DATABASES

    output:
    file "${method}/*"
    file "${task.process}/*" optional true

    when:
    MLST_DATABASES.isEmpty() == false

    shell:
    method = dataset =~ /.*blastdb.*/ ? 'blast' : 'ariba'
    dataset_tarball = file(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '').split('-')[1]
    schema = dataset_tarball.split('-')[0]
    noclean = params.ariba_no_clean ? "--noclean" : ""
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    template(task.ext.template)
}


process ariba_analysis {
    /* Run reads against all available (if any) ARIBA datasets */
    tag "${sample} - ${dataset_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/ariba", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${dataset_name}/*"

    input:
    set val(sample), val(single_end), file(fq) from ARIBA_ANALYSIS
    each file(dataset) from ARIBA_DATABASES

    output:
    file "${dataset_name}/*"
    file "${task.process}/*" optional true

    when:
    single_end == false && ARIBA_DATABASES.isEmpty() == false

    shell:
    dataset_tarball = file(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    noclean = params.ariba_no_clean ? "--noclean" : ""
    template(task.ext.template)
}


process minmer_sketch {
    /*
    Create minmer sketches of the input FASTQs using Mash (k=21,31) and
    Sourmash (k=21,31,51)
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/minmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.{msh,sig}"

    input:
    set val(sample), val(single_end), file(fq) from MINMER_SKETCH

    output:
    file("${sample}*.{msh,sig}")
    set val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("${sample}.sig") into MINMER_QUERY
    set val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("${sample}-k31.msh") into DOWNLOAD_REFERENCES
    file "${task.process}/*" optional true

    shell:
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template(task.ext.template)
}


process minmer_query {
    /*
    Query minmer sketches against pre-computed RefSeq (Mash, k=21) and
    GenBank (Sourmash, k=21,31,51)
    */
    tag "${sample} - ${dataset_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/minmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.txt"

    input:
    set val(sample), val(single_end), file(fq), file(sourmash) from MINMER_QUERY
    each file(dataset) from MINMER_DATABASES

    output:
    file "*.txt"
    file "${task.process}/*" optional true

    when:
    MINMER_DATABASES.isEmpty() == false

    shell:
    dataset_name = dataset.getName()
    mash_w = params.screen_w ? "-w" : ""
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template(task.ext.template)
}


process call_variants {
    /*
    Identify variants (SNPs/InDels) against a set of reference genomes
    using Snippy.
    */
    tag "${sample} - ${reference_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/variants/user", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${reference_name}/*"

    input:
    set val(sample), val(single_end), file(fq) from CALL_VARIANTS
    each file(reference) from REFERENCES

    output:
    file "${reference_name}/*"
    file "${task.process}/*" optional true

    when:
    REFERENCES.isEmpty() == false

    shell:
    snippy_ram = task.memory.toString().split(' ')[0]
    reference_name = reference.getSimpleName()
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    bwaopt = params.bwaopt ? "--bwaopt 'params.bwaopt'" : ""
    fbopt = params.fbopt ? "--fbopt 'params.fbopt'" : ""
    template(task.ext.template)
}


process download_references {
    /*
    Download the nearest RefSeq genomes (based on Mash) to have variants called against.

    Exitcode 75 is due to being unable to download from NCBI (e.g. FTP down at the time)
    Downloads will be attempted 300 times total before giving up. On failure to download
    variants will not be called against the nearest completed genome.
    */
    tag "${sample} - ${params.max_references} reference(s)"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/variants/auto", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: 'mash-dist.txt'

    input:
    set val(sample), val(single_end), file(fq), file(sample_sketch) from DOWNLOAD_REFERENCES
    file(refseq_sketch) from REFSEQ_SKETCH

    output:
    set val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz"), file("genbank/*.gbk") optional true into CALL_VARIANTS_AUTO
    file("mash-dist.txt")
    file "${task.process}/*" optional true

    when:
    REFSEQ_SKETCH_FOUND == true

    shell:
    no_cache = params.no_cache ? '-N' : ''
    tie_break = params.random_tie_break ? "--random_tie_break" : ""
    total = params.max_references
    template(task.ext.template)
}


process call_variants_auto {
    /*
    Identify variants (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${sample} - ${reference_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/variants/auto", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${reference_name}/*"

    input:
    set val(sample), val(single_end), file(fq), file(reference) from create_reference_channel(CALL_VARIANTS_AUTO)

    output:
    file "${reference_name}/*"
    file "${task.process}/*" optional true

    shell:
    snippy_ram = task.memory.toString().split(' ')[0]
    reference_name = reference.getBaseName().split("${sample}-")[1].split(/\./)[0]
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    bwaopt = params.bwaopt ? "--bwaopt 'params.bwaopt'" : ""
    fbopt = params.fbopt ? "--fbopt 'params.fbopt'" : ""
    template(task.ext.template)
}


process antimicrobial_resistance {
    /*
    Query nucleotides and proteins (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "logs/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${amrdir}/*"

    input:
    set val(sample), file(genes), file(proteins) from ANTIMICROBIAL_RESISTANCE
    each file(amrdb) from AMR_DATABASES

    output:
    file "${amrdir}/*"
    file "logs/*" optional true

    shell:
    amrdir = "antimicrobial-resistance"
    plus = params.amr_plus ? "--plus" : ""
    report_common = params.amr_report_common ? "--report_common" : ""
    organism_gene = ""
    organism_protein = ""
    if (params.amr_organism) {
        organism_gene = "-O ${params.amr_organism} --point_mut_all ${amrdir}/${sample}-gene-point-mutations.txt"
        organism_protein = "-O ${params.amr_organism} --point_mut_all ${amrdir}/${sample}-protein-point-mutations.txt"
    }
    template(task.ext.template)
}

process blast_primers {
    /*
    Query primer FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "primers/*.{json,json.gz}"

    input:
    set val(sample), file(blastdb) from BLAST_PRIMERS
    file(query) from Channel.from(BLAST_PRIMER_FASTAS).collect()

    output:
    file("primers/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_PRIMER_FASTAS.isEmpty() == false

    shell:
    template(task.ext.template)
}

process blast_proteins {
    /*
    Query protein FASTA files against annotated assembly using BLAST
    */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "proteins/*.{json,json.gz}"

    input:
    set val(sample), file(blastdb) from BLAST_PROTEINS
    file(query) from Channel.from(BLAST_PROTEIN_FASTAS).collect()

    output:
    file("proteins/*.{json,json.gz}")
    file "${task.process}/*" optional true

    when:
    BLAST_PROTEIN_FASTAS.isEmpty() == false

    shell:
    template(task.ext.template)
}


process mapping_query {
    /*
    Map FASTQ reads against a given set of FASTA files using BWA.
    */
    tag "${sample}}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "mapping/*"

    input:
    set val(sample), val(single_end), file(fq) from MAPPING_QUERY
    file(query) from Channel.from(MAPPING_FASTAS).collect()

    output:
    file "mapping/*"
    file "${task.process}/*" optional true

    when:
    MAPPING_FASTAS.isEmpty() == false

    shell:
    bwa_mem_opts = params.bwa_mem_opts ? params.bwa_mem_opts : ""
    bwa_aln_opts = params.bwa_aln_opts ? params.bwa_aln_opts : ""
    bwa_samse_opts = params.bwa_samse_opts ? params.bwa_samse_opts : ""
    bwa_sampe_opts = params.bwa_sampe_opts ? params.bwa_sampe_opts : ""
    template(task.ext.template)
}


workflow.onComplete {
    workDir = new File("${workflow.workDir}")
    workDirSize = toHumanString(workDir.directorySize())

    println """

    Bactopia Execution Summary
    ---------------------------
    Command Line    : ${workflow.commandLine}
    Resumed         : ${workflow.resume}

    Completed At    : ${workflow.complete}
    Duration        : ${workflow.duration}
    Success         : ${workflow.success}
    Exit Code       : ${workflow.exitStatus}
    Error Report    : ${workflow.errorReport ?: '-'}

    Launch Dir      : ${workflow.launchDir}
    Working Dir     : ${workflow.workDir} (Total Size: ${workDirSize})
    Working Dir Size: ${workDirSize}

    """
}


// Utility functions
def toHumanString(bytes) {
    // Thanks Niklaus
    // https://gist.github.com/nikbucher/9687112
    base = 1024L
    decimals = 3
    prefix = ['', 'K', 'M', 'G', 'T']
    int i = Math.log(bytes)/Math.log(base) as Integer
    i = (i >= prefix.size() ? prefix.size()-1 : i)
    return Math.round((bytes / base**i) * 10**decimals) / 10**decimals + prefix[i]
}

def print_version() {
    println(PROGRAM_NAME + ' ' + VERSION)
    exit 0
}


def print_basedir() {
    println("BACTOPIA_DIR = ${baseDir}")
    exit 0
}

def print_example_fastqs() {
    log.info """
        Printing example input for "--fastqs"

        sample	runtype	r1	r2	extra
        SA103113	assembly			/example/SA103113.fna.gz
        SA110685	hybrid	/example/SA110685_R1.fastq.gz	/example/SA110685_R2.fastq.gz	/example/SA110685.fastq.gz
        SA123186	paired-end	/example/SA123186_R1.fastq.gz	/example/SA123186_R2.fastq.gz
        SA123456	single-end	/example/SA12345.fastq.gz
    """.stripIndent()
    exit 0
}


def print_check_fastqs(run_type) {
    if (run_type == "fastqs") {
        log.info 'Printing what would have been processed. Each line consists of an array of'
        log.info 'five elements: [SAMPLE_NAME, RUNTYPE, IS_SINGLE_END, [FASTQ_1], [FASTQ_2], EXTRA]'
        log.info ''
        log.info 'Found:'
        create_input_channel(run_type).view()
        log.info ''
        exit 0
    } else {
        log.error '"--check_fastqs" requires "--fastqs" to be set.'
        exit 1
    }
}


def print_available_datasets(dataset) {
    exit_code = 0
    if (dataset) {
        if (file("${dataset}/summary.json").exists()) {
            available_datasets = read_json("${dataset}/summary.json")
            log.info 'Printing the available pre-configured dataset.'
            log.info "Database Location (--datasets): ${dataset}"
            log.info ''
            if (available_datasets.size() > 0) {
                IGNORE = ['species-specific']
                GENERAL = ['ariba', 'minmer', 'plasmid']
                available_datasets.each { key, value ->
                    if (GENERAL.contains(key) == 'ariba') {
                        log.info "${key.capitalize()}"
                        value.each {
                            log.info "\tFound ${it}"
                        }
                    } else if (key == 'species-specific') {
                        value.each { species, sets ->
                            log.info "${species.capitalize().replace('-', ' ')} (use --species \"${species}\")"
                            sets.each { set_name, set_path ->
                                if (set_name == 'optional') {
                                    log.info "\tOptional:"
                                    set_path.each {
                                        log.info "\t\tFound ${it}"
                                    }
                                } else {
                                    log.info "\tFound ${set_name}=${set_path}"
                                }
                            }
                            log.info ''
                        }
                    } else {
                        log.info "${key.capitalize()}"
                        value.each {
                            log.info "\tFound ${it}"
                        }
                    }
                    log.info ''
                }
            }
        } else {
            log.error "Please verify the PATH is correct and ${dataset}/summary.json" +
                     " exists, if not try rerunning 'bactopia datasets'."
            exit_code = 1
        }
    } else {
        log.error "Please use '--datasets' to specify the path to pre-built datasets."
        exit_code = 1
    }
    exit exit_code
}


def read_json(json_file) {
    slurp = new JsonSlurper()
    return slurp.parseText(file(json_file).text)
}


def format_species(species) {
    /* Format species name to accepted format. */
    name = false
    if (species) {
        name = species.toLowerCase()
        if (species.contains(" ")) {
            name = name.replace(" ", "-")
        }
    }

    return name
}


def get_max_memory(requested) {
    available = Math.floor(Double.parseDouble(SysHelper.getAvailMemory().toString().split(" ")[0])).toInteger()

    if (available < requested) {
        log.warn "Maximum memory (${requested}) was adjusted to fit your system (${available})"
        return available
    }

    return requested
}


def get_max_cpus(requested) {
    available = SysHelper.getAvailCpus()
    if (available < requested) {
        log.warn "Maximum CPUs (${requested}) was adjusted to fit your system (${available})"
        return available
    }

    return requested
}

def dataset_exists(dataset_path) {
    if (file(dataset_path).exists()) {
        return true
    } else {
        log.warn "Warning: ${dataset_path} does not exist and will not be used for analysis."
    }
}

def setup_datasets() {
    genome_size = ['min': 0, 'median': 0, 'mean': 0, 'max': 0]
    if (params.datasets) {
        dataset_path = params.datasets
        available_datasets = read_json("${dataset_path}/summary.json")

        // Antimicrobial Resistance Datasets
        if (params.skip_amr) {
            log.info "Found '--skip_amr', datasets for AMRFinder+ will not be used for analysis."
        } else {
            if (available_datasets.containsKey('antimicrobial-resistance')) {
                available_datasets['antimicrobial-resistance'].each {
                    if (dataset_exists("${dataset_path}/antimicrobial-resistance/${it.name}")) {
                        AMR_DATABASES << file("${dataset_path}/antimicrobial-resistance/${it.name}")
                    }
                }
                print_dataset_info(AMR_DATABASES, "Antimicrobial resistance datasets")
            }
        }


        // Ariba Datasets
        if (available_datasets.containsKey('ariba')) {
            available_datasets['ariba'].each {
                if (dataset_exists("${dataset_path}/ariba/${it.name}")) {
                    ARIBA_DATABASES << file("${dataset_path}/ariba/${it.name}")
                }
            }
            print_dataset_info(ARIBA_DATABASES, "ARIBA datasets")
        }

        // RefSeq/GenBank Check
        if (available_datasets.containsKey('minmer')) {
            if (available_datasets['minmer'].containsKey('sketches')) {
                available_datasets['minmer']['sketches'].each {
                    if (dataset_exists("${dataset_path}/minmer/${it}")) {
                        MINMER_DATABASES << file("${dataset_path}/minmer/${it}")
                    }
                }
            }
        }

        // PLSDB Check
        if (available_datasets.containsKey('plasmid')) {
            if (available_datasets['plasmid'].containsKey('sketches')) {
                if (dataset_exists("${dataset_path}/plasmid/${available_datasets['plasmid']['sketches']}")) {
                    MINMER_DATABASES << file("${dataset_path}/plasmid/${available_datasets['plasmid']['sketches']}")
                }
            }
        }
        print_dataset_info(MINMER_DATABASES, "minmer sketches/signatures")

        if (SPECIES) {
            if (available_datasets.containsKey('species-specific')) {
                if (available_datasets['species-specific'].containsKey(SPECIES)) {
                    species_db = available_datasets['species-specific'][SPECIES]
                    if (species_db.containsKey('genome_size')) {
                        genome_size = species_db['genome_size']
                    }

                    if (params.genome_size) {
                        if (['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
                            SPECIES_GENOME_SIZE = genome_size[params.genome_size]
                        } else {
                            SPECIES_GENOME_SIZE = params.genome_size
                        }

                        if (SPECIES_GENOME_SIZE > 0) {
                            log.info "Will use ${SPECIES_GENOME_SIZE} bp for genome size"
                        } else if (SPECIES_GENOME_SIZE == 0) {
                            log.info "Found ${SPECIES_GENOME_SIZE} bp for genome size, it will be estimated."
                        }
                    } else {
                        SPECIES_GENOME_SIZE = null
                        log.info "Genome size will be estimated."
                    }

                    if (species_db.containsKey('annotation')) {
                        if (species_db['annotation'].containsKey('proteins')) {
                            prokka = "${dataset_path}/${species_db['annotation']['proteins']}"
                            if (dataset_exists(prokka)) {
                                PROKKA_PROTEINS = file(prokka)
                                log.info "Found Prokka proteins file"
                                log.info "\t${PROKKA_PROTEINS}"
                            }
                        }

                        if (species_db['annotation'].containsKey('training_set')) {
                            prodigal_tf = "${dataset_path}/${species_db['annotation']['training_set']}"
                            if (dataset_exists(prodigal_tf)) {
                                PRODIGAL_TF = file(prodigal_tf)
                                log.info "Found Prodigal training file"
                                log.info "\t${PRODIGAL_TF}"
                            }
                        }
                    }

                    if (species_db.containsKey('minmer')) {
                        if (species_db['minmer'].containsKey('mash')) {
                            refseq_minmer = "${dataset_path}/${species_db['minmer']['mash']}"
                            if (dataset_exists(refseq_minmer)) {
                                REFSEQ_SKETCH = file(refseq_minmer)
                                REFSEQ_SKETCH_FOUND = true
                                log.info "Found Mash Sketch of RefSeq genomes"
                                log.info "\t${REFSEQ_SKETCH}"
                            }
                        }
                    }

                    if (species_db.containsKey('mlst')) {
                        species_db['mlst'].each { schema, vals ->
                            vals.each { key, val ->
                                if (key != "last_updated") {
                                    if (dataset_exists("${dataset_path}/${val}")) {
                                        MLST_DATABASES << file("${dataset_path}/${val}")
                                    }
                                }
                            }
                        }
                        print_dataset_info(MLST_DATABASES, "MLST datasets")
                    }

                    if (species_db.containsKey('optional')) {
                        if (species_db['optional'].containsKey('reference-genomes')) {
                            file("${dataset_path}/${species_db['optional']['reference-genomes']}").list().each() {
                                if (dataset_exists("${dataset_path}/${species_db['optional']['reference-genomes']}/${it}")) {
                                    REFERENCES << file("${dataset_path}/${species_db['optional']['reference-genomes']}/${it}")
                                }
                            }
                            print_dataset_info(REFERENCES, "reference genomes")
                        }

                        if (species_db['optional'].containsKey('mapping-sequences')) {
                            file("${dataset_path}/${species_db['optional']['mapping-sequences']}").list().each() {
                                if (dataset_exists("${dataset_path}/${species_db['optional']['mapping-sequences']}/${it}")) {
                                    MAPPING_FASTAS << file("${dataset_path}/${species_db['optional']['mapping-sequences']}/${it}")
                                }
                            }
                            print_dataset_info(MAPPING_FASTAS, "FASTAs to align reads against")
                        }

                        if (species_db['optional'].containsKey('blast')) {
                            species_db['optional']['blast'].each() {
                                blast_type = it
                                temp_path = "${dataset_path}/${it}"
                                file(temp_path).list().each() {
                                    if (dataset_exists("${temp_path}/${it}")) {
                                        if (blast_type.contains('blast/genes')) {
                                            BLAST_GENE_FASTAS << file("${temp_path}/${it}")
                                        } else if (blast_type.contains('blast/primers')) {
                                            BLAST_PRIMER_FASTAS << file("${temp_path}/${it}")
                                        } else {
                                            BLAST_PROTEIN_FASTAS << file("${temp_path}/${it}")
                                        }
                                    }
                                }
                            }
                            print_dataset_info(BLAST_GENE_FASTAS, "gene FASTAs to query with BLAST")
                            print_dataset_info(BLAST_PRIMER_FASTAS, "primer FASTAs to query with BLAST")
                            print_dataset_info(BLAST_PROTEIN_FASTAS, "protein FASTAs to query with BLAST")
                        }
                    }
                } else {
                    log.error "Species '${params.species}' not available, please check spelling " +
                              "or use '--available_datasets' to verify the dataset has been set " +
                              "up. Exiting"
                    exit 1
                }
            } else {
                log.error "Species '${params.species}' not available, please check spelling or " +
                          "use '--available_datasets' to verify the dataset has been set up. Exiting"
                exit 1
            }
        } else {
            log.info "--species not given, skipping the following processes (analyses):"
            log.info "\tsequence_type"
            log.info "\tcall_variants"
            log.info "\tprimer_query"
            if (['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
                log.error "Asked for genome size '${params.genome_size}' which requires a " +
                          "species to be given. Please give a species or specify " +
                          "a valid genome size. Exiting"
                exit 1
            }
        }
    } else {
        log.info "--datasets not given, skipping the following processes (analyses):"
        log.info "\tsequence_type"
        log.info "\tariba_analysis"
        log.info "\tminmer_query"
        log.info "\tcall_variants"
        log.info "\tprimer_query"
    }

    if (params.disable_auto_variants) {
        REFSEQ_SKETCH_FOUND = false
    }

    if (USING_MERGE) {
        log.info "\n"
        log.warn "One or more samples consists of multiple read sets and will have their reads merged. "
        log.warn "This is an experimental feature, please use with caution."
        log.info "\n"
    }

    log.info "\nIf something looks wrong, now's your chance to back out (CTRL+C 3 times). "
    log.info "Sleeping for ${params.sleep_time} seconds..."
    sleep(params.sleep_time * 1000)
}


def is_positive_integer(value, name) {
    error = 0
    if (value.getClass() == Integer) {
        if (value < 0) {
            log.error('Invalid input (--'+ name +'), "' + value + '"" is not a positive integer.')
            error = 1
        }
    } else {
        if (!value.isInteger()) {
            log.error('Invalid input (--'+ name +'), "' + value + '"" is not numeric.')
            error = 1
        } else if (value.toInteger() < 0) {
            log.error('Invalid input (--'+ name +'), "' + value + '"" is not a positive integer.')
            error = 1
        }
    }
    return error
}


def file_exists(file_name, parameter) {
    if (!file(file_name).exists()) {
        log.error('Invalid input ('+ parameter +'), please verify "' + file_name + '" exists.')
        return 1
    } else if (['--R1', '--R2', '--SE', '--assembly'].contains(parameter)) {
        if (!isGzipped(file_name)) {
            log.error('Invalid input ('+ parameter +'), please verify "' + file_name + '" is GZIP compressed.')
            return 1
        }
    }
    return 0
}


def check_unknown_params() {
    valid_params = []
    error = 0
    new File("${baseDir}/conf/params.config").eachLine { line ->
        if (line.contains("=")) {
            valid_params << line.trim().split(" ")[0]
        }
    }

    params.each { k,v ->
        if (!valid_params.contains(k)) {
            if (k != "container-path") {
                log.error("'--${k}' is not a known parameter")
                error = 1
            }
        }
    }

    return error
}


def check_input_params() {
    error = 0
    run_type = null

    // Check for unexpected paramaters
    error += check_unknown_params()

    if (params.fastqs) {
        error += file_exists(params.fastqs, '--fastqs')
        run_type = "fastqs"
    } else if  (params.R1 && params.R2 && params.SE && params.hybrid && params.sample) {
        error += file_exists(params.R1, '--R1')
        error += file_exists(params.R2, '--R2')
        error += file_exists(params.SE, '--SE')
        run_type = "hybrid"
    } else if  (params.R1 && params.R2 && params.SE) {
        log.error "Cannot use --R1, --R2, and --SE together, unless --hybrid is used."
    } else if  (params.R1 && params.R2 && params.sample) {
        error += file_exists(params.R1, '--R1')
        error += file_exists(params.R2, '--R2')
        run_type = "paired-end"
    } else if (params.SE && params.sample) {
        error += file_exists(params.SE, '--SE')
        run_type = "single-end"
    } else if (params.assembly && params.sample) {
        error += file_exists(params.assembly, '--assembly')
        run_type = "assembly"
    } else if (params.accessions) {
        error += file_exists(params.accessions, '--accessions')
        run_type = "is_accessions"
    } else if (params.accession) {
        run_type = "is_accession"
    } else {
        log.error """
        One or more required parameters are missing, please check and try again.

        Required Parameters:
            ### For Procesessing Multiple Samples
            --fastqs STR            An input file containing the sample name and
                                        absolute paths to FASTQ/FASTAs to process

            ### For Processing A Single Sample
            --R1 STR                First set of reads for paired end in compressed (gzip)
                                        FASTQ format

            --R2 STR                Second set of reads for paired end in compressed (gzip)
                                        FASTQ format

            --SE STR                Single end set of reads in compressed (gzip) FASTQ format

            --hybrid                The SE should be treated as long reads for hybrid assembly.

            --sample STR            The name of the input sequences

            ### For Downloading from SRA/ENA or NCBI Assembly
            **Note: Assemblies will have error free Illumina reads simulated for processing.**
            --accessions            An input file containing ENA/SRA Experiment accessions or
                                        NCBI Assembly accessions to be processed

            --accession             A single ENA/SRA Experiment accession or NCBI Assembly accession
                                        to be processed

            ### For Processing an Assembly
            **Note: The assembly will have error free Illumina reads simulated for processing.**
            --assembly STR          A assembled genome in compressed FASTA format.
        """.stripIndent()
        error += 1
    }

    error += is_positive_integer(params.cpus, 'cpus')
    error += is_positive_integer(params.max_time, 'max_time')
    error += is_positive_integer(params.max_memory, 'max_memory')
    error += is_positive_integer(params.max_downloads, 'max_downloads')
    error += is_positive_integer(params.max_retry, 'max_retry')
    error += is_positive_integer(params.min_time, 'min_time')
    error += is_positive_integer(params.shovill_ram, 'shovill_ram')
    error += is_positive_integer(params.snippy_ram, 'snippy_ram')
    error += is_positive_integer(params.cortex_ram, 'cortex_ram')
    error += is_positive_integer(params.qc_ram, 'qc_ram')
    error += is_positive_integer(params.minmer_ram, 'minmer_ram')

    if (params.max_downloads >= 10) {
        log.warn "Please be aware the value you have set for --max_downloads (${params.max_downloads}) may cause NCBI " +
                 "to temporarily block your IP address due to too many queries at once."
    }

    if (params.genome_size) {
        if (!['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
            error += is_positive_integer(params.genome_size, 'genome_size')
        }
    }

    if (!['megahit', 'velvet', 'skesa', 'spades', 'unicycler'].contains(params.assembler)) {
        log.error "Invalid assembler (--assembler ${params.assembler}), must be 'megahit', " +
                  "'skesa', 'velvet', 'spades', or 'unicycler'. Please correct to continue."
        error += 1
    }

    if (!['dockerhub', 'github', 'quay'].contains(params.registry)) {
        log.error "Invalid registry (--registry ${params.registry}), must be 'dockerhub', " +
                    "'github' or 'quay'. Please correct to continue."
        error += 1
    }

    if (params.min_time > params.max_time) {
        log.error "The value for min_time (${params.min_time}) exceeds max_time (${params.max_time}), Please correct to continue."
        error += 1
    }

    if (params.adapters) {
        error += file_exists(params.adapters, '--adapters')
    }

    if (params.phix) {
        error += file_exists(params.phix, '--phix')
    }

    if (params.datasets) {
        if (!file("${params.datasets}/summary.json").exists()) {
            log.error "Please verify the PATH is correct for '--datasets'. Unable " +
                      "to open ${params.datasets}/summary.json"
            error += 1
        }
    }

    // Check for existing output directory
    if (!workflow.resume) {
        if (!file(params.outdir).exists() && !params.force) {
            log.error("Output directory (${params.outdir}) exists, Bactopia will not continue unless '--force' is used.")
            error += 1
        }
    }

    if (error > 0) {
        log.error('Cannot continue, please see --help for more information')
        exit 1
    }

    return run_type
}

def handle_multiple_fqs(read_set) {
    def fqs = []
    def String[] reads = read_set.split(",");
    reads.each { fq ->
        fqs << file(fq)
    }
    return fqs
}



def process_fastqs(line) {
    /* Parse line and determine if single end or paired reads*/
    if (line.runtype == 'single-end') {
        return tuple(line.sample, line.runtype, true, [file(line.r1)], [null], null)
    } else if (line.runtype == 'paired-end') {
        return tuple(line.sample, line.runtype, false, [file(line.r1)], [file(line.r2)], null)
    } else if (line.runtype == 'hybrid') {
        return tuple(line.sample, line.runtype, false, [file(line.r1)], [file(line.r2)], file(line.extra))
    } else if (line.runtype == 'assembly') {
        return tuple(line.sample, line.runtype, false, [null], [null], file(line.extra))
    } else if (line.runtype == 'merge-pe') {
        return tuple(line.sample, line.runtype, false, handle_multiple_fqs(line.r1), handle_multiple_fqs(line.r2), null)
    } else if (line.runtype == 'hybrid-merge-pe') {
        return tuple(line.sample, line.runtype, false, handle_multiple_fqs(line.r1), handle_multiple_fqs(line.r2), file(line.extra))
    } else if (line.runtype == 'merge-se') {
        return tuple(line.sample, line.runtype, false, handle_multiple_fqs(line.r1), [null], null)
    } else {
        log.error("Invalid run_type ${line.runtype} found, please correct to continue. Expected: single-end, paired-end, hybrid, merge-pe, hybrid-merge-pe, merge-se, or assembly")
        exit 1
    }
}

def process_accessions(accession) {
    /* Parse line and determine if single end or paired reads*/
    if (accession.length() > 0) {
        if (accession.startsWith('GCF') || accession.startsWith('GCA')) {
            return tuple(accession.split(/\./)[0], "assembly_accession", false, [null], [null], null)
        } else if (accession.startsWith('DRX') || accession.startsWith('ERX') || accession.startsWith('SRX')) {
            return tuple(accession, "sra_accession", false, [null], [null], null)
        } else {
            log.error("Invalid accession: ${accession} is not an accepted accession type. Accessions must be Assembly (GCF_*, GCA*) or Exeriment (DRX*, ERX*, SRX*) accessions. Please correct to continue.\n\nYou can use 'bactopia search' to convert BioProject, BioSample, or Run accessions into an Experiment accession.")
            exit 1
        }
    }
}


def create_input_channel(run_type) {
    if (run_type == "fastqs") {
        return Channel.fromPath( file(params.fastqs) )
            .splitCsv(header: true, sep: '\t')
            .map { row -> process_fastqs(row) }
    } else if (run_type == "is_accessions") {
        return Channel.fromPath( file(params.accessions) )
            .splitText()
            .map { line -> process_accessions(line.trim()) }
    } else if (run_type == "is_accession") {
        return [process_accessions(params.accession)]
    } else if (run_type == "paired-end") {
        return [tuple(params.sample, run_type, false, [file(params.R1)], [file(params.R2)], null)]
    } else if (run_type == "hybrid") {
        return [tuple(params.sample, run_type, false, [file(params.R1)], [file(params.R2)], file(params.SE))]
    } else if (run_type == "assembly") {
        return [tuple(params.sample, run_type, false, [null], [null], file(params.assembly))]
    } else {
        return [tuple(params.sample, run_type, true, [file(params.SE)], [null], null)]
    }
}


def create_reference_channel(channel_input) {
    return channel_input.map { it ->
        def split_references = []
        def sample = it[0]
        def single_end = it[1]
        def fastq = it[2]
        if (it[3] instanceof List) {
            it[3].each {
                split_references.add(tuple(sample, single_end, fastq, it))
            }
        } else {
            split_references.add(tuple(sample, single_end, fastq, it[3]))
        }
        return split_references
     }.flatMap { it -> L:{ it.collect { it } } }
}


def isGzipped(f) {
    /* https://github.com/ConnectedPlacesCatapult/TomboloDigitalConnector/blob/master/src/main/java/uk/org/tombolo/importer/ZipUtils.java */
    int magic = 0;
    try {
        RandomAccessFile raf = new RandomAccessFile(f, "r");
        magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
        raf.close();
    } catch (Throwable e) {
        log.error("Failed to check if gzipped " + e.getMessage());
        return false;
    }

    return magic == GZIPInputStream.GZIP_MAGIC;
}


def check_input_fastqs(run_type) {
    /* Read through --fastqs and verify each input exists. */

    if (run_type == "fastqs") {
        samples = [:]
        error = false
        has_valid_header = false
        line = 1
        file(params.fastqs).splitEachLine('\t') { cols ->
            if (line == 1) {
                if (cols[0] == 'sample' && cols[1] == 'runtype' && cols[2] == 'r1' && cols[3] == 'r2' && cols[4] == 'extra') {
                    has_valid_header = true
                }
            } else {
                if (samples.containsKey(cols[0])) {
                    samples[cols[0]] = samples[cols[0]] + 1
                } else {
                    samples[cols[0]] = 1
                }
                if (cols[2]) {
                    count = 0
                    cols[2].split(',').each{ fq ->
                        if (!file(fq).exists()) {
                            log.error "LINE " + line + ':ERROR: Please verify ' + fq + ' exists, and try again'
                            error = true
                        }
                        count = count + 1
                    }
                    if (count > 1) {
                        USING_MERGE = true
                    }
                }
                if (cols[3]) {
                    cols[3].split(',').each{ fq ->
                        if (!file(fq).exists()) {
                            log.error "LINE " + line + ':ERROR: Please verify ' + fq + ' exists, and try again'
                            error = true
                        }
                    }
                }
                if (cols[4]) {
                    if (!file(cols[4]).exists()) {
                        log.error "LINE " + line + ':ERROR: Please verify ' + cols[4]+ ' exists, and try again'
                        error = true
                    }
                }
            }

            line = line + 1
        }

        samples.each{ sample, count ->
            if (count > 1) {
                error = true
                log.error 'Sample name "'+ sample +'" is not unique, please revise sample names'
            }
        }

        if (!has_valid_header) {
            error = true
            log.error 'The header line (line 1) does not follow expected structure. (sample, runtype, r1, r2, extra)'
        }

        if (error) {
            log.error 'Verify sample names are unique and/or FASTA/FASTQ paths are correct'
            log.error 'See "--example_fastqs" for an example'
            log.error 'Exiting'
            exit 1
        }
    }
}


def get_absolute_path(file_path) {
    // Thanks Fabian Steeg
    // https://stackoverflow.com/questions/3204955/converting-relative-paths-to-absolute-paths
    File file_obj = new File("${workflow.launchDir}/${file_path}");
    String absolute_path = file_obj.getAbsolutePath(); // may throw IOException

    return absolute_path
}


def get_canonical_path(file_path) {
    // Thanks Fabian Steeg
    // https://stackoverflow.com/questions/3204955/converting-relative-paths-to-absolute-paths
    File file_obj = new File("${workflow.launchDir}/${file_path}");
    String canonical_path = file_obj.getCanonicalPath(); // may throw IOException

    return canonical_path
}


def print_dataset_info(dataset_list, dataset_info) {
    if (dataset_list.size() > 0) {
        log.info "Found ${dataset_list.size()} ${dataset_info}"
        count = 0
        try {
            dataset_list.each {
                if (count < 5) {
                    log.info "\t${it}"
                } else {
                    log.info "\t...More than 5, truncating..."
                    throw new Exception("break")
                }
                count++
            }
        } catch (Exception e) {}
    }
}

def print_efficiency() {
    if (['standard', 'docker', 'singularity'].contains(workflow.profile)) {
        // This is a local run on a single machine
        total_cpus = Runtime.runtime.availableProcessors()
        tasks = total_cpus / MAX_CPUS
        log.info ""
        log.info """
            Each task will use ${MAX_CPUS} CPUs out of the available ${total_cpus} CPUs. At most ${tasks} task(s) will be run at
            a time, this can affect the efficiency of Bactopia.
        """.stripIndent()
        log.info ""
    }

}


def print_usage() {
    usage_text = params.help_all ? full_help() : basic_help()
    log.info"""
    ${PROGRAM_NAME} v${VERSION}
    ${basic_help()}
    ${params.help_all ? full_help() : ""}
    """.stripIndent()

    if (params.conda_help) {
        // Cleanup up the directory
        // This is only meant to be used with tests for conda build
        file("./work/").deleteDir()
        file("./.nextflow/").deleteDir()
        def files = new FileNameByRegexFinder().getFileNames('./', '.nextflow.log*')
        files.each { new File(it).delete()}
    }
    exit 0
}


def basic_help() {
    genome_size = params.genome_size ? params.genome_size : "Mash Estimate"
    return """
    Required Parameters:
        ### For Procesessing Multiple Samples
        --fastqs STR            An input file containing the sample name and
                                    absolute paths to FASTQ/FASTAs to process

        ### For Processing A Single Sample
        --R1 STR                First set of reads for paired end in compressed (gzip)
                                    FASTQ format

        --R2 STR                Second set of reads for paired end in compressed (gzip)
                                    FASTQ format

        --SE STR                Single end set of reads in compressed (gzip) FASTQ format

        --hybrid                The SE should be treated as long reads for hybrid assembly.

        --sample STR            The name of the input sequences

        ### For Downloading from SRA/ENA or NCBI Assembly
        **Note: Assemblies will have error free Illumina reads simulated for processing.**
        --accessions            An input file containing ENA/SRA Experiment accessions or
                                    NCBI Assembly accessions to be processed

        --accession             A single ENA/SRA Experiment accession or NCBI Assembly accession
                                    to be processed

        ### For Processing an Assembly
        **Note: The assembly will have error free Illumina reads simulated for processing.**
        --assembly STR          A assembled genome in compressed (gzip) FASTA format.

        --reassemble            The simulated reads will be used to create a new assembly.
                                    Default: Use the original assembly, do not reassemble

    Dataset Parameters:
        --datasets DIR          The path to available datasets that have
                                    already been set up

        --species STR           Determines which species-specific dataset to
                                    use for the input sequencing

    Optional Parameters:
        --coverage INT          Reduce samples to a given coverage
                                    Default: ${params.coverage}x

        --genome_size INT       Expected genome size (bp) for all samples, a value of '0'
                                    will disable read error correction and read subsampling.
                                    Special values (requires --species):
                                        'min': uses minimum completed genome size of species
                                        'median': uses median completed genome size of species
                                        'mean': uses mean completed genome size of species
                                        'max': uses max completed genome size of species
                                    Default: Mash estimate

        --outdir DIR            Directory to write results to
                                    Default: ${params.outdir}

    Nextflow Queue Parameters:
        At execution, Nextflow creates a queue and the number of slots in the queue is determined by the total number
        of cores on the system. When a task is submitted to the queue, the total number of slots it occupies is
        determined by the value set by "--cpus".

        This can have a significant effect on the efficiency of the Nextflow's queue system. If "--cpus" is set to a
        value that is equal to the number of cores availabe, in most cases only a single task will be able to run
        because its occupying all available slots.

        When in doubt, "--cpus 4" is a safe bet, it is also the default value if you don't use "--cpus".

        --max_retry INT         Maximum times to retry a process before allowing it to fail.
                                    Default: ${params.max_retry}

        --min_time INT          The minimum number of minutes a single task should run before being halted.
                                    Default: ${params.min_time} minutes

        --max_time INT          The maximum number of minutes a single task should run before being halted.
                                    Default: ${params.max_time} minutes

        --max_memory INT        The maximum amount of memory (Gb) allowed to a single task.
                                    Default: ${params.max_memory} Gb

        --cpus INT              Number of processors made available to a single task.
                                    Default: ${params.cpus}

        -qs INT                 Nextflow queue size. This parameter is very useful to limit the total number of
                                    processors used on desktops, laptops or shared resources.
                                    Default: Nextflow defaults to the total number of processors on your system.


    Nextflow Related Parameters:
        --infodir DIR           Directory to write Nextflow summary files to
                                    Default: ${params.infodir}

        --condadir DIR          Directory to Nextflow should use for Conda environments
                                    Default: Bactopia's Nextflow directory

        --registry STR          Docker registry to pull containers from.
                                    Available options: dockerhub, quay, or github
                                    Default: dockerhub

        --singularity_cache STR Directory where remote Singularity images are stored. If using a cluster, it must
                                    be accessible from all compute nodes.
                                    Default: NXF_SINGULARITY_CACHEDIR evironment variable, otherwise ${params.singularity_cache}

        --queue STR             The name of the queue(s) to be used by a job scheduler (e.g. AWS Batch or SLURM).
                                    If using multiple queues, please seperate queues by a comma without spaces.
                                    Default: ${params.queue}

        --disable_scratch       All intermediate files created on worker nodes of will be transferred to the head node.
                                    Default: Only result files are transferred back

        --nfconfig STR          A Nextflow compatible config file for custom profiles. This allows
                                    you to create profiles specific to your environment (e.g. SGE,
                                    AWS, SLURM, etc...). This config file is loaded last and will
                                    overwrite existing variables if set.
                                    Default: Bactopia's default configs

        --nfdir                 Print directory Nextflow has pulled Bactopia to

        --overwrite             Nextflow will overwrite existing output files.
                                    Default: ${params.overwrite}

        --sleep_time            After reading datases, the amount of time (seconds) Nextflow
                                    will wait before execution.
                                    Default: ${params.sleep_time} seconds

        --publish_mode          Set Nextflow's method for publishing output files. Allowed methods are:
                                    'copy' (default)    Copies the output files into the published directory.

                                    'copyNoFollow' Copies the output files into the published directory
                                                   without following symlinks ie. copies the links themselves.

                                    'link'    Creates a hard link in the published directory for each
                                              process output file.

                                    'rellink' Creates a relative symbolic link in the published directory
                                              for each process output file.

                                    'symlink' Creates an absolute symbolic link in the published directory
                                              for each process output file.

                                    Default: ${params.publish_mode}

        --force                 Nextflow will overwrite existing output files.
                                    Default: ${params.force}

        -resume                 Nextflow will attempt to resume a previous run. Please notice it is
                                    only a single '-'

        --cleanup_workdir       After Bactopia is successfully executed, the work directory will be deleted.
                                    Warning: by doing this you lose the ability to resume workflows.

    Useful Parameters:
        --skip_logs             Logs for each process per sample will not be kept.

        --available_datasets    Print a list of available datasets found based
                                    on location given by "--datasets"

        --example_fastqs        Print example of expected input for FASTQs file

        --check_fastqs          Verify "--fastqs" produces the expected inputs

        --compress              Compress (gzip) select outputs (e.g. annotation, variant calls)
                                    to reduce overall storage footprint.

        --keep_all_files        Keeps all analysis files created. By default, intermediate
                                    files are removed. This will not affect the ability
                                    to resume Nextflow runs, and only occurs at the end
                                    of the process.

        --version               Print workflow version information

        --help                  Show this message and exit

        --help_all              Show a complete list of adjustable parameters
    """
}


def full_help() {
    return """
    Additional Parameters:
    The description of the following parameters were taken from the program for
    which they apply to.

    Many of the default values were also taken from the program for which they
    apply to.

    AWS Batch Profile  (-profile awsbatch) Parameters:
        --aws_region STR        AWS Region to be used by Nextflow
                                    Default: ${params.aws_region}

        --aws_volumes STR       Volumes to be mounted from the EC2 instance to the Docker container
                                    Default: ${params.aws_volumes}

        --aws_cli_path STR       Path to the AWS CLI for Nextflow to use.
                                    Default: ${params.aws_cli_path}

        --aws_upload_storage_class STR
                                The S3 storage slass to use for storing files on S3
                                    Default: ${params.aws_upload_storage_class}

        --aws_max_parallel_transfers INT
                                The number of parallele transfers between EC2 and S3
                                    Default: ${params.aws_max_parallel_transfers}

        --aws_delay_between_attempts INT
                                The duration of sleep (in seconds) between each transfer between EC2 and S3
                                    Default: ${params.aws_delay_between_attempts}

        --aws_max_transfer_attempts INT
                                The maximum number of times to retry transferring a file between EC2 and S3
                                    Default: ${params.aws_max_transfer_attempts}

        --aws_max_retry INT     The maximum number of times to retry a process on AWS Batch
                                    Default: ${params.aws_max_retry}

        --aws_ecr_registry STR  The ECR registry containing Bactopia related containers.
                                    Default: Use the registry given by --registry

    ENA Download Parameters:
        --max_downloads INT     Maximum number of FASTQs to download at once.
                                    Warning: Setting this value too high can lead to NCBI temporarily
                                             blocking your IP addess. 3-5 is reasonable, >10 is likely
                                             to be excessive.
                                    Default: ${params.max_downloads}

        --use_ena               Download FASTQs from ENA with Aspera Connect.
                                    Default: Download from SRA

        --no_cache              Don't cache the assembly summary file from ncbi-genome-download

    FASTQ Minimum Requirements Parameters:
        --min_basepairs INT     The minimum amount of input sequenced basepairs required
                                    to continue downstream analyses.
                                    Default: ${params.min_basepairs}

        --min_coverage INT      The minimum coverage of input sequences required
                                    to continue downstream analyses.
                                    Default: ${params.min_coverage}

        --min_reads INT         The minimum amount of input sequenced reads required
                                    to continue downstream analyses.
                                    Default: ${params.min_reads}

        --min_proportion FLOAT  The minimum proportion of basepairs for paired-end reads to continue
                                    downstream analyses. Example: If set to 0.75 the R1 and R2 must
                                    have > 75% proportion of reads (e.g. R1 100bp, R2 75bp, not
                                    R1 100bp, R2 50bp)
                                    Default: ${params.min_proportion}

        --skip_fastq_check      The input FASTQs will not be check to verify they meet the
                                    minimum requirements to be processed. This parameter
                                    is useful if you are confident your sequences will
                                    pass the minimum requirements.

    Estimate Genome Size Parameters:
        Only applied if the genome size is estimated.

        --min_genome_size INT   The minimum estimated genome size allowed for the input sequence
                                to continue downstream analyses.
                                Default: ${params.min_genome_size}

        --max_genome_size INT   The maximum estimated genome size allowed for the input sequence
                                to continue downstream analyses.
                                Default: ${params.max_genome_size}

    QC Reads Parameters:
        --skip_qc               The QC step qill be skipped and it will be assumed the inputs
                                    sequences have already been QCed.

        --skip_error_correction FLASH error correction of paired-end reads will be skipped.

        --qc_ram INT            Try to keep RAM usage below this many GB
                                    Default: ${params.qc_ram} GB

        --adapters FASTA        Illumina adapters to remove
                                    Default: BBmap adapters

        --adapter_k INT         Kmer length used for finding adapters. Adapters
                                    shorter than k will not be found
                                    Default: ${params.adapter_k}

        --phix FASTA            phiX174 reference genome to remove
                                    Default: NC_001422

        --phix_k INT            Kmer length used for finding phiX174.
                                    Contaminants shorter than k will not be
                                    found
                                    Default: ${params.phix_k}

        --ktrim STR             Trim reads to remove bases matching reference
                                    kmers
                                    Values:
                                        f (do not trim)
                                        r (trim to the right, Default)
                                        l (trim to the left)

        --mink INT              Look for shorter kmers at read tips down to this
                                    length, when k-trimming or masking. 0 means
                                    disabled. Enabling this will disable
                                    maskmiddle
                                    Default: ${params.mink}

        --hdist INT             Maximum Hamming distance for ref kmers (subs only)
                                    Memory use is proportional to (3*K)^hdist
                                    Default: ${params.hdist}

        --tpe BOOL              When kmer right-trimming, trim both reads to the
                                    minimum length of either
                                    Values:
                                        f (do not equally trim)
                                        t (equally trim to the right, Default)

        --tbo BOOL              Trim adapters based on where paired reads overlap
                                    Values:
                                        f (do not trim by overlap)
                                        t (trim by overlap, Default)

        --qtrim STR             Trim read ends to remove bases with quality
                                    below trimq. Performed AFTER looking for
                                    kmers
                                    Values:
                                        rl (trim both ends, Default)
                                        f (neither end)
                                        r (right end only)
                                        l (left end only)
                                        w (sliding window)

        --trimq FLOAT           Regions with average quality BELOW this will be
                                    trimmed if qtrim is set to something other
                                    than f
                                    Default: ${params.trimq}

        --maq INT               Reads with average quality (after trimming)
                                    below this will be discarded
                                    Default: ${params.maq}

        --minlength INT         Reads shorter than this after trimming will be
                                    discarded. Pairs will be discarded if both
                                    are shorter
                                    Default: ${params.minlength}

        --ftm INT               If positive, right-trim length to be equal to
                                    zero, modulo this number
                                    Default: ${params.ftm}

        --tossjunk              Discard reads with invalid characters as bases
                                    Values:
                                        f (keep all reads)
                                        t (toss reads with ambiguous bases, Default)

        --qout STR              Output quality offset
                                    Values:
                                        33 (PHRED33 offset quality scores, Default)
                                        64 (PHRED64 offset quality scores)
                                        auto (keeps the current input offset)

        --maxcor INT            Max number of corrections within a 20bp window
                                    Default: ${params.maxcor}

        --sampleseed INT        Set to a positive number to use as the rng seed
                                    for sampling
                                    Default: ${params.sampleseed}


    Assembly Parameters:
        Standard Assembly:
        --shovill_ram INT       Try to keep RAM usage below this many GB
                                    Default: ${params.shovill_ram} GB

        --assembler STR         Assembler: megahit velvet skesa spades unicycler
                                    Default: ${params.assembler}

        --min_contig_len INT    Minimum contig length <0=AUTO>
                                    Default: ${params.min_contig_len}

        --min_contig_cov INT    Minimum contig coverage <0=AUTO>
                                    Default: ${params.min_contig_cov}

        --contig_namefmt STR    Format of contig FASTA IDs in 'printf' style
                                    Default: "SAMPLE_NAME_%05d"

        --shovill_opts STR      Extra assembler options in quotes eg.
                                    spades: "--untrusted-contigs locus.fna" ...
                                    Default: ''

        --shovill_kmers STR     K-mers to use <blank=AUTO>
                                    Default: ''

        --trim                  Enable adaptor trimming

        --nostitch              Disable read stitching

        --nocorr                Disable post-assembly correction

        Hybrid Assembly:
        --unicycler_ram INT       Try to keep RAM usage below this many GB
                                    Default: ${params.unicycler_ram} GB

        --unicycler_mode STR    Bridging mode used by Unicycler, choices are:
                                    conservative = smaller contigs, lowest
                                                   misassembly rate
                                    normal = moderate contig size and
                                             misassembly rate (Default)
                                    bold = longest contigs, higher misassembly
                                           rate

        --min_polish_size INT   Contigs shorter than this value (bp) will not be
                                    polished using Pilon
                                    Default: ${params.min_polish_size}

        --min_component_size INT
                                Graph components smaller than this size (bp) will
                                    be removed from the final graph
                                    Default: ${params.min_component_size}

        --min_dead_end_size INT
                                Graph dead ends smaller than this size (bp) will
                                    be removed from the final graph
                                    Default: ${params.min_dead_end_size}

        --no_miniasm            Skip miniasm+Racon bridging
                                    Default: Produce long-read bridges

        --no_rotate             Do not rotate completed replicons to start at a
                                    standard gene

        --no_pilon              Do not use Pilon to polish the final assembly

    Assembly Quality Control Parameters:
        --skip_checkm           CheckM analysis will be skipped. This is useful for systems
                                    with less than 8GB of memory.

        --checkm_unique INT     Minimum number of unique phylogenetic markers required
                                    to use lineage-specific marker set.
                                    Default: ${params.checkm_unique}

        --checkm_multi INT      Maximum number of multi-copy phylogenetic markers before
                                    defaulting to domain-level marker set.
                                    Default: ${params.checkm_multi}

        --aai_strain FLOAT      AAI threshold used to identify strain heterogeneity
                                    Default: ${params.aai_strain}

        --checkm_length FLOAT   Percent overlap between target and query
                                    Default: ${params.checkm_length}

        --full_tree             Use the full tree (requires ~40GB of memory) for determining
                                    lineage of each bin.
                                    Default: Use reduced tree (<16gb memory)

        --skip_pseudogene_correction
                                Skip identification and filtering of pseudogene

        --ignore_thresholds     Ignore model-specific score thresholds

        --checkm_ali            Generate HMMER alignment file for each bin

        --checkm_nt             Generate nucleotide gene sequences for each bin

        --force_domain          Use domain-level sets for all bins

        --no_refinement         Do not perform lineage-specific marker set refinement

        --individual_markers    Treat marker as independent (i.e., ignore co-located
                                    set structure.

        --skip_adj_correction   Do not exclude adjacent marker genes when estimating
                                    contamination

        --contig_thresholds STR Comma-separated list of contig length thresholds
                                    Default: ${params.contig_thresholds}

        --plots_format STR      Save plots in specified format.
                                    Supported formats: emf, eps, pdf, png, ps, raw,
                                                        rgba, svg, svgz
                                    Default: ${params.plots_format}

    Count 31mers Parameters:
        --cortex_ram INT        Try to keep RAM usage below this many GB
                                    Default: ${params.cortex_ram} GB

        --keep_singletons       Keep all counted 31-mers
                                    Default: Filter out singletons

    Annotation Parameters:
        --compliant             Force Genbank/ENA/DDJB compliance: --genes --mincontiglen ${params.min_contig_len} --centre '${params.centre}'
                                    Default: ${params.compliant}

        --centre STR            Sequencing centre ID
                                    Default: '${params.centre}'

        --addmrna               Add 'mRNA' features for each 'CDS' feature

        --rawproduct            Do not clean up /product annotation

        --cdsrnaolap            Allow [tr]RNA to overlap CDS

        --prokka_evalue STR     Similarity e-value cut-off
                                    Default: ${params.prokka_evalue}

        --prokka_coverage INT   Minimum coverage on query protein
                                     Default: ${params.prokka_coverage}

        --nogenes               Do not add 'gene' features for each 'CDS' feature

        --norrna                Don't run rRNA search

        --notrna                Don't run tRNA search

        --rnammer               Prefer RNAmmer over Barrnap for rRNA prediction

        --rfam                  Enable searching for ncRNAs with Infernal+Rfam

        --skip_prodigal_tf      If a Prodigal training file was found, it will not be used

    Minmer Sketch Parameters:
        --mash_sketch INT       Sketch size. Each sketch will have at most this
                                    many non-redundant min-hashes.
                                    Default: ${params.mash_sketch}

        --sourmash_scale INT    Choose number of hashes as 1 in FRACTION of
                                    input k-mers
                                    Default: ${params.sourmash_scale}


    Minmer Query Parameters:
        --minmer_ram INT        Try to keep RAM usage below this many GB
                                    Default: ${params.minmer_ram} GB

        --screen_w              Winner-takes-all strategy for identity estimates.
                                    After counting hashes for each query, hashes
                                    that appear in multiple queries will be
                                    removed from all except the one with the best
                                    identity (ties broken by larger query), and
                                    other identities will be reduced. This
                                    removes output redundancy, providing a rough
                                    compositional outline.
                                    Default: True

        --screen_i FLOAT        Minimum identity to report. Inclusive unless set
                                    to zero, in which case only identities greater
                                    than zero (i.e. with at least one shared hash)
                                    will be reported. Set to -1 to output
                                    everything.
                                    Default: ${params.screen_i}

    Ariba Parameters:
        --nucmer_min_id INT     Minimum alignment identity (delta-filter -i)
                                    Default: ${params.nucmer_min_id}

        --nucmer_min_len INT    Minimum alignment length (delta-filter -i)
                                    Default: ${params.nucmer_min_len}

        --nucmer_breaklen INT   Value to use for -breaklen when running nucmer
                                    Default: ${params.nucmer_breaklen}

        --assembly_cov INT      Target read coverage when sampling reads for
                                    assembly
                                    Default: ${params.assembly_cov}

        --min_scaff_depth INT   Minimum number of read pairs needed as evidence
                                    for scaffold link between two contigs
                                    Default: ${params.min_scaff_depth}

        --spades_options STR    Extra options to pass to Spades assembler
                                    Default: ${params.spades_options}

        --assembled_threshold FLOAT (between 0 and 1)
                                If proportion of gene assembled (regardless of
                                    into how many contigs) is at least this
                                    value then the flag gene_assembled is set
                                    Default: ${params.assembled_threshold}

        --gene_nt_extend INT    Max number of nucleotides to extend ends of gene
                                    matches to look for start/stop codons
                                    Default: ${params.gene_nt_extend}

        --unique_threshold FLOAT (between 0 and 1)
                                If proportion of bases in gene assembled more
                                    than once is <= this value, then the flag
                                    unique_contig is set
                                    Default: ${params.unique_threshold}

        --ariba_no_clean        Do not clean up intermediate files created by
                                    Ariba. By default, the local assemblies are
                                    deleted.

    Call Variant Parameters:
        --snippy_ram INT        Try and keep RAM under this many GB
                                    Default: ${params.snippy_ram} GB

        --mapqual INT           Minimum read mapping quality to consider
                                    Default: ${params.mapqual}

        --basequal INT          Minimum base quality to consider
                                    Default: ${params.basequal}

        --mincov INT            Minimum site depth to for calling alleles
                                    Default: ${params.mincov}

        --minfrac FLOAT         Minimum proportion for variant evidence (0=AUTO)
                                    Default: ${params.minfrac}

        --minqual INT           Minimum QUALITY in VCF column 6
                                    Default: ${params.minqual}

        --maxsoft INT           Maximum soft clipping to allow
                                    Default: ${params.maxsoft}

        --bwaopt STR            Extra BWA MEM options, eg. -x pacbio
                                    Default: ''

        --fbopt STR             Extra Freebayes options,
                                    eg. --theta 1E-6 --read-snp-limit 2
                                    Default: ''

    Nearest Neighbor Reference Genomes:
        --max_references INT    Maximum number of nearest neighbor reference genomes to
                                    download for variant calling.
                                    Default: 1

        --random_tie_break      On references with matching distances, randomly select one.
                                    Default: Pick earliest accession number

        --disable_auto_variants Disable automatic selection of reference genome based on
                                    Mash distances. This will not skip reference genomes
                                    you may have provided in the optional folders.

    BLAST Parameters:
        --perc_identity INT     Percent identity
                                    Default: ${params.perc_identity}

        --qcov_hsp_perc INT     Percent query coverage per hsp
                                    Default: ${params.qcov_hsp_perc}

        --max_target_seqs INT   Maximum number of aligned sequences to
                                    keep
                                    Default: ${params.max_target_seqs}

    Mapping Parameters:
        --keep_unmapped_reads   Keep unmapped reads, this does not affect variant
                                    calling.

        --bwa_mem_opts STR      Extra BWA MEM options
                                    Default: ''

        --bwa_aln_opts STR      Extra BWA ALN options
                                    Default: ''

        --bwa_samse_opts STR    Extra BWA SAMSE options
                                    Default: ''

        --bwa_sampe_opts STR    Extra BWA SAMPE options
                                    Default: ''

        --bwa_n INT             Maximum number of alignments to output in the XA
                                    tag for reads paired properly. If a read has
                                    more than INT hits, the XA tag will not be
                                    written.
                                    Default: ${params.bwa_n}

    Antimicrobial Resistance Parameters:
        --skip_amr              AMRFinder+ analysis will be skipped. This is useful
                                    if the AMRFinder+ software and database versions are
                                    no longer compatible.

        --amr_ident_min         Minimum identity for nucleotide hit (0..1). -1
                                    means use a curated threshold if it exists and
                                    0.9 otherwise
                                    Default: ${params.amr_ident_min}

        --amr_coverage_min      Minimum coverage of the reference protein (0..1)
                                    Default: ${params.amr_coverage_min}

        --amr_organism          Taxonomy group: Campylobacter, Escherichia, Klebsiella
                                    Salmonella, Staphylococcus, Vibrio
                                    Default: ''

        --amr_translation_table NCBI genetic code for translated BLAST
                                    Default: ${params.amr_translation_table}

        --amr_plus              Add the plus genes to the report

        --amr_report_common     Suppress proteins common to a taxonomy group

    """
}
