nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options = initOptions(params.containsKey('options') ? params.options : [:], 'sequence_type')

process SEQUENCE_TYPE {
    /* Determine MLST types using ARIBA and BLAST */
    tag "${meta.id} - ${schema}"
    label "max_cpus_1"
    label "sequence_type"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options, logs_subdir: schema) }

    input:
    tuple val(meta), path(assembly), path(fq)
    each path(dataset)

    output:
    path "results/*", emit: results
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    shell:
    dataset_tarball = dataset.getName()
    schema = dataset_tarball.replace('.tar.gz', '')
    noclean = params.mlst_ariba_no_clean ? "--noclean" : ""
    spades_options = params.mlst_spades_options ? "--spades_options '${params.mlst_spades_options}'" : ""
    blast_opts = params.skip_compression ? "" : "--compressed"
    '''
    tar -xzvf !{dataset_tarball}

    # Run BLAST
    OUTDIR="results/!{schema}"
    mkdir -p ${OUTDIR}/blast
    mlst-blast.py !{assembly} !{schema}/blastdb ${OUTDIR}/blast/!{meta.id}-blast.json --cpu !{task.cpus} !{blast_opts}

    # Run Ariba
    if [ "!{meta.single_end}" == "false" ]; then
        mv !{schema}/ariba/ref_db ./
        ariba run ref_db !{fq[0]} !{fq[1]} ariba \
            --nucmer_min_id !{params.mlst_nucmer_min_id} \
            --nucmer_min_len !{params.mlst_nucmer_min_len} \
            --nucmer_breaklen !{params.mlst_nucmer_breaklen} \
            --assembly_cov !{params.mlst_assembly_cov} \
            --min_scaff_depth !{params.mlst_min_scaff_depth} \
            --assembled_threshold !{params.mlst_assembled_threshold} \
            --gene_nt_extend !{params.mlst_gene_nt_extend} \
            --unique_threshold !{params.mlst_unique_threshold} \
            --threads 1 \
            --force \
            --verbose !{noclean} !{spades_options}
        mv ariba ${OUTDIR}/
    else
        mkdir -p ${OUTDIR}/ariba
        echo "Ariba cannot be run on single end reads" >${OUTDIR}/ariba/ariba-not-run.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ariba:  $(echo $(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*$//')
        blastn: $(echo $(blastn -version 2>&1) | sed 's/^.*blastn: //;s/ .*$//')
    END_VERSIONS
    '''
}
