nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "sequence_type"

process SEQUENCE_TYPE {
    /* Determine MLST types using ARIBA and BLAST */
    tag "${meta.id} - ${schema}"
    label PROCESS_NAME

    publishDir "${params.outdir}/${meta.id}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, process_name:PROCESS_NAME, logs_subdir: schema) }

    input:
    tuple val(meta), path(fq), path(assembly)
    each path(dataset)

    output:
    path "results/*", emit: results
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    shell:
    dataset_tarball = dataset.getName()
    schema = dataset_tarball.replace('.tar.gz', '')
    noclean = params.ariba_no_clean ? "--noclean" : ""
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    blast_opts = params.skip_compression ? "" : "--compressed"
    '''
    tar -xzvf !{dataset_tarball}

    # Run BLAST
    OUTDIR="results/!{schema}"
    mkdir -p ${OUTDIR}/blast
    mlst-blast.py !{assembly} !{schema}/blastdb ${OUTDIR}/blast/!{meta.id}-blast.json --cpu !{task.cpus} !{blast_opts} 1> mlst-blast.stdout.txt 2> mlst-blast.stderr.txt

    # Run Ariba
    if [ "!{meta.single_end}" == "false" ]; then
        mv !{schema}/ariba/ref_db ./
        ariba run ref_db !{fq[0]} !{fq[1]} ariba \
            --nucmer_min_id !{params.nucmer_min_id} \
            --nucmer_min_len !{params.nucmer_min_len} \
            --nucmer_breaklen !{params.nucmer_breaklen} \
            --assembly_cov !{params.assembly_cov} \
            --min_scaff_depth !{params.min_scaff_depth} \
            --assembled_threshold !{params.assembled_threshold} \
            --gene_nt_extend !{params.gene_nt_extend} \
            --unique_threshold !{params.unique_threshold} \
            --threads 1 \
            --force \
            --verbose !{noclean} !{spades_options} 1> ariba.stdout.txt 2> ariba.stderr.txt
        mv ariba ${OUTDIR}/
    else
        mkdir -p ${OUTDIR}/ariba
        echo "Ariba cannot be run on single end reads" >${OUTDIR}/ariba/ariba-not-run.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    sequence_type:
        ariba:  $(echo $(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*$//')
        blastn: $(echo $(blastn -version 2>&1) | sed 's/^.*blastn: //;s/ .*$//')
    END_VERSIONS
    '''

    stub:
    dataset_tarball = path(dataset).getName()
    schema = dataset_tarball.replace('.tar.gz', '')
    """
    mkdir -p results/${schema}/ariba
    touch results/${schema}/ariba/ariba-not-run.txt
    """
}
