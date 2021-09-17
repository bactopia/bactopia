nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "minmer_sketch"

process MINMER_SKETCH {
    /*
    Create minmer sketches of the input FASTQs using Mash (k=21,31),
    Sourmash (k=21,31,51), and McCortex (k=31)
    */
    tag "${meta.id}"
    label PROCESS_NAME

    publishDir "${params.outdir}/${meta.id}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME, ignore: [".fastq.gz"]) }

    input:
    tuple val(meta), path(fq)

    output:
    tuple val(meta), path(fq), path("${meta.id}.sig"), emit: sketch
    path("${meta.id}*.{msh,sig}")
    path("${meta.id}.ctx"), optional: true
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    shell:
    fastq = meta.single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    mccortex_fq = meta.single_end ? "-1 ${fq[0]}" : "-2 ${fq[0]}:${fq[1]}"
    m = task.memory.toString().split(' ')[0].toInteger() * 1000 - 500
    '''
    gzip -cd !{fastq} | mash sketch -o !{meta.id}-k21 -k 21 -s !{params.mash_sketch} -r -I !{meta.id} - > mash.stdout.txt 2> mash.stderr.txt
    gzip -cd !{fastq} | mash sketch -o !{meta.id}-k31 -k 31 -s !{params.mash_sketch} -r -I !{meta.id} - >> mash.stdout.txt 2>> mash.stderr.txt
    sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=!{params.sourmash_scale} --merge !{meta.id} -o !{meta.id}.sig !{fastq} > sourmash.stdout.txt 2> sourmash.stderr.txt

    if [[ "!{params.count_31mers}" == "true" ]]; then
        mccortex31 build -f -k 31 -s !{meta.id} !{mccortex_fq} -t !{task.cpus} -m !{m}mb -q temp_counts > mccortex.stdout.txt 2> mccortex.stderr.txt

        if [ "!{params.keep_singletons}" == "false" ]; then
            # Clean up Cortex file (mostly remove singletons)
            mccortex31 clean -q -B 2 -U2 -T2 -m !{m}mb -o !{meta.id}.ctx temp_counts >> mccortex.stdout.txt 2>> mccortex.stderr.txt
            rm temp_counts
        else
            mv temp_counts !{meta.id}.ctx
        fi
        mccortex31 2>&1 | grep "version" > mccortex31.version.txt 2>&1
    fi

    # Capture versions
    mash --version > mash.version.txt 2>&1
    sourmash --version > sourmash.version.txt 2>&1
    '''

    stub:
    """
    touch ${sample}.sig
    touch ${sample}-k31.msh
    """
}
