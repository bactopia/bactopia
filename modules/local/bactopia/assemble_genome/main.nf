nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "assemble_genome"

process ASSEMBLE_GENOME {
    /* Assemble the genome using Shovill, SKESA is used by default */
    tag "${meta.id}"
    label "max_cpu_75"
    label PROCESS_NAME

    publishDir "${params.outdir}/${meta.id}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME, ignore: ["-genome-size.txt", ".fastq.gz"]) }

    input:
    tuple val(meta), path(fq), path(extra), path(genome_size)

    output:
    tuple val(meta), path(genome_size), path("results/${meta.id}.{fna,fna.gz}"), file("total_contigs_*"), emit: fna, optional: true
    tuple val(meta), path(fq), path("results/${meta.id}.{fna,fna.gz}"), emit: fna_fastq, optional: true
    tuple val(meta), path("blastdb/*"), emit: blastdb, optional: true
    path "results/*"
    path "${meta.id}-assembly-error.txt", optional: true
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    shell:
    // Unicycler
    no_miniasm = params.no_miniasm ? "--no_miniasm" : ""
    no_rotate = params.no_rotate ? "--no_rotate" : ""
    no_pilon = params.no_pilon ? "--no_pilon" : ""
    keep = params.keep_all_files ? "--keep 3" : "--keep 1"
    is_hybrid = meta.runtype == "hybrid" ? "-l ${extra}" : ""
    unicycler_opts = "${no_miniasm} ${no_rotate} ${no_pilon} ${is_hybrid} ${keep}"

    // Shovill
    contig_namefmt = params.contig_namefmt ? params.contig_namefmt : "${meta.id}_%05d"
    shovill_ram = task.memory.toString().split(' ')[0]
    opts = params.shovill_opts ? "--opts '${params.shovill_opts}'" : ""
    kmers = params.shovill_kmers ? "--kmers '${params.shovill_kmers}'" : ""
    nostitch = params.nostitch ? "--nostitch" : ""
    nocorr = params.nocorr ? "--nocorr" : ""
    shovill_mode = single_end == false ? "shovill --R1 ${fq[0]} --R2 ${fq[1]} --noreadcorr ${nostitch}" : "shovill-se --SE ${fq[0]}"
    shovill_opts = "${opts} ${kmers} ${nocorr}"

    // Assembly inputs
    use_original_assembly = null
    if (meta.runtype.startsWith('assembly')) {
        use_original_assembly = params.reassemble ? false : true
    }
    '''
    OUTDIR=results
    GENOME_SIZE=`head -n 1 !{genome_size}`
    if [ "!{use_original_assembly}" == "true" ]; then
        mkdir ${OUTDIR}
        gzip -cd !{extra} > ${OUTDIR}/!{meta.id}.fna
    elif [[ "!{meta.runtype}" == "hybrid"  || "!{params.assembler}" == "unicycler" ]]; then
        unicycler -1 !{fq[0]} -2 !{fq[1]} -o ${OUTDIR}/ --no_correct \
            --min_fasta_length !{params.min_contig_len} \
            --threads !{task.cpus} \
            --mode !{params.unicycler_mode} \
            --min_polish_size !{params.min_polish_size} \
            --min_component_size !{params.min_component_size} \
            --min_dead_end_size !{params.min_dead_end_size} !{unicycler_opts} > unicycler.stdout.txt 2> unicycler.stderr.txt
        sed -r 's/^>([0-9]+)(.*)/>!{meta.id}_\\1\\2/' ${OUTDIR}/assembly.fasta > ${OUTDIR}/!{meta.id}.fna

        # Capture Version
        unicycler --version > unicycler.version.txt 2>&1
    else
        !{shovill_mode} --depth 0 --gsize ${GENOME_SIZE} \
            --outdir ${OUTDIR} \
            --force \
            --minlen !{params.min_contig_len} \
            --mincov !{params.min_contig_cov} \
            --namefmt "!{contig_namefmt}" \
            --keepfiles \
            --cpus !{task.cpus} \
            --ram !{shovill_ram} \
            --assembler !{params.assembler} !{shovill_opts} > shovill.stdout.txt 2> shovill.stderr.txt
        mv ${OUTDIR}/contigs.fa ${OUTDIR}/!{meta.id}.fna

        # Capture version
        shovill --version > shovill.version.txt
        shovill --check > shovil-depends.version.txt 2>&1
        if [ "!{params.assembler}" == "spades" ]; then
            spades.py --version > spades.version.txt 2>&1
        elif [ "!{params.assembler}" == "skesa" ]; then
            skesa --version 2>&1 | tail -n 1 > skesa.version.txt 2>&1
        elif [ "!{params.assembler}" == "velvet" ]; then
            velvetg | grep "^Version" > velvet.version.txt 2>&1
        else
            megahit --version > megahit.version.txt 2>&1
        fi
    fi

    # Check quality of assembly
    TOTAL_CONTIGS=`grep -c "^>" ${OUTDIR}/!{meta.id}.fna || true`
    touch "total_contigs_${TOTAL_CONTIGS}"
    if [ "${TOTAL_CONTIGS}" -gt "0" ]; then
        assembly-scan ${OUTDIR}/!{meta.id}.fna > ${OUTDIR}/!{meta.id}.fna.json 2> assembly-scan.stderr.txt
        assembly-scan --version > assembly-scan.version.txt 2>&1
        TOTAL_CONTIG_SIZE=`grep "total_contig_length" ${OUTDIR}/!{meta.id}.fna.json | sed -r 's/.*: ([0-9]+)/\1/'`
        if [ ${TOTAL_CONTIG_SIZE} -lt "!{params.min_genome_size}" ]; then
            mv ${OUTDIR}/!{meta.id}.fna ${OUTDIR}/!{meta.id}-error.fna
            mv ${OUTDIR}/!{meta.id}.fna.json ${OUTDIR}/!{meta.id}-error.fna.json
            echo "!{meta.id} assembled size (${TOTAL_CONTIG_SIZE} bp) is less than the minimum allowed genome
                    size (!{params.min_genome_size} bp). If this is unexpected, please investigate !{meta.id} to
                    determine a cause (e.g. metagenomic, contaminants, etc...) for the poor assembly.
                    Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                    based analysis of !{meta.id} will be discontinued." | \
            sed 's/^\\s*//' > !{meta.id}-assembly-error.txt
        else
            # Make BLASTDB
            mkdir blastdb
            cat ${OUTDIR}/!{meta.id}.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for !{meta.id}" -out blastdb/!{meta.id}
            makeblastdb -version > makeblastdb.version.txt 2>&1
        fi
    else
        echo "!{meta.id} assembled successfully, but 0 contigs were formed. Please investigate
                !{meta.id} to determine a cause (e.g. metagenomic, contaminants, etc...) for this
                outcome. Further assembly-based analysis of !{meta.id} will be discontinued." | \
        sed 's/^\\s*//' > !{meta.id}-assembly-error.txt
    fi

    # Cleanup and compress
    if [ "!{params.keep_all_files}" == "false" ]; then
        # Remove intermediate files
        rm -rfv ${OUTDIR}/shovill.bam* ${OUTDIR}/flash.extendedFrags* ${OUTDIR}/flash.notCombined* \
                ${OUTDIR}/skesa.fasta.* ${OUTDIR}/*.fq.gz ${OUTDIR}/00*.gfa ${OUTDIR}/pilon_polish*
    fi

    if [[ !{params.skip_compression} == "false" ]]; then
        # Compress based on matched extensions
        find ${OUTDIR}/ -type f | \
            grep -E "\\.fna$|\\.fasta$|\\.fa$|\\.gfa$|\\.fastg$|\\.LastGraph$" | \
            xargs -I {} pigz -n --best -p !{task.cpus} {}
    fi
    '''

    stub:
    """
    mkdir assembly
    mkdir fastqs
    mkdir ${PROCESS_NAME}
    touch total_contigs_${meta.id}
    touch ${meta.id}-assembly-error.txt
    touch fastqs/${meta.id}.fastq.gz
    touch assembly/${meta.id}
    touch assembly/${meta.id}.fna
    touch assembly/${meta.id}.fna.gz
    touch ${PROCESS_NAME}/${sampmeta.idle}
    """
}
