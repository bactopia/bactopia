
nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "qc_reads"

process QC_READS {
    /* Clean up Illumina reads */
    tag "${meta.id}"
    label "max_cpus"
    label PROCESS_NAME

    publishDir "${params.outdir}/${meta.id}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME, ignore:[ '-genome-size.txt', extra]) }

    input:
    tuple val(meta), path(fq), path(extra), path(genome_size)

    output:
    tuple val(meta), path("results/${meta.id}*.fastq.gz"), emit: fastq, optional: true
    tuple val(meta), path("results/${meta.id}*.fastq.gz"), path(extra), path(genome_size), emit: fastq_assembly, optional: true
    path "results/*"
    path "*.std{out,err}.txt", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version
    path "*-error.txt", optional: true


    shell:
    meta.single_end = fq[1] == null ? true : false
    qc_ram = task.memory.toString().split(' ')[0]
    is_assembly = meta.runtype.startsWith('assembly') ? true : false
    qin = meta.runtype.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    adapters = params.adapters ? path(params.adapters) : 'adapters'
    phix = params.phix ? path(params.phix) : 'phix'
    adapter_opts = meta.single_end ? "" : "in2=${fq[1]} out2=adapter-r2.fq"
    phix_opts = meta.single_end ? "" : "in2=adapter-r2.fq out2=phix-r2.fq"
    lighter_opts = meta.single_end ? "" : "-r phix-r2.fq"
    reformat_opts = meta.single_end ? "" : "in2=phix-r2.cor.fq out2=subsample-r2.fq"
    '''
    mkdir -p results
    ERROR=0
    GENOME_SIZE=`head -n 1 !{genome_size}`
    MIN_BASEPAIRS=$(( !{params.min_coverage}*${GENOME_SIZE} ))
    TOTAL_BP=$(( !{params.coverage}*${GENOME_SIZE} ))

    if [ "!{params.skip_qc}" == "true" ]; then
        echo "Sequence QC was skipped for !{meta.id}" > results/!{meta.id}-qc-skipped.txt
        if [[ -L "!{fq[0]}" ]]; then
            if [ "!{meta.single_end}" == "false" ]; then
                # Paired-End Reads
                ln -s `readlink !{fq[0]}` results/!{meta.id}_R1.fastq.gz
                ln -s `readlink !{fq[1]}` results/!{meta.id}_R2.fastq.gz
            else
                # Single-End Reads
                ln -s `readlink !{fq[0]}` results/!{meta.id}.fastq.gz
            fi
        else
            if [ "!{meta.single_end}" == "false" ]; then
                # Paired-End Reads
                cp !{fq[0]} results/!{meta.id}_R1.fastq.gz
                cp !{fq[1]} results/!{meta.id}_R2.fastq.gz
            else
                # Single-End Reads
                cp !{fq[0]} results/!{meta.id}.fastq.gz
            fi
        fi
    else
        bbduk.sh --version 2>&1 | grep " version" > bbduk.version.txt 2>&1

        # Remove Adapters
        bbduk.sh -Xmx!{qc_ram}g \
            in=!{fq[0]} out=adapter-r1.fq !{adapter_opts} \
            ref=!{adapters} \
            k=!{params.adapter_k} \
            ktrim=!{params.ktrim} \
            mink=!{params.mink} \
            hdist=!{params.hdist} \
            tpe=!{params.tpe} \
            tbo=!{params.tbo} \
            threads=!{task.cpus} \
            ftm=!{params.ftm} \
            !{qin} ordered=t \
            stats=bbduk-adapter.stdout.txt 2> bbduk-adapter.stderr.txt

        # Remove PhiX
        bbduk.sh -Xmx!{qc_ram}g \
            in=adapter-r1.fq out=phix-r1.fq !{phix_opts} \
            ref=!{phix} \
            k=!{params.phix_k} \
            hdist=!{params.hdist} \
            tpe=!{params.tpe} \
            tbo=!{params.tbo} \
            qtrim=!{params.qtrim} \
            trimq=!{params.trimq} \
            minlength=!{params.minlength} \
            minavgquality=!{params.maq} \
            !{qin} qout=!{params.qout} \
            tossjunk=!{params.tossjunk} \
            threads=!{task.cpus} \
            ordered=t stats=bbduk-phix.stdout.txt 2> bbduk-phix.stderr.txt

        # Error Correction
        if [ "!{params.skip_error_correction}" == "false" ]; then
            lighter -v > lighter.version.txt 2>&1
            lighter -od . -r phix-r1.fq !{lighter_opts} -K 31 ${GENOME_SIZE} -maxcor 1 -zlib 0 -t !{task.cpus} 1> lighter.stdout.txt 2> lighter.stderr.txt
        else
            echo "Skipping error correction"
            ln -s phix-r1.fq phix-r1.cor.fq
            if [ "!{meta.single_end}" == "false" ]; then
                ln -s phix-r2.fq phix-r2.cor.fq
            fi
        fi

        # Reduce Coverage
        if (( ${TOTAL_BP} > 0 )); then
            reformat.sh -Xmx!{qc_ram}g \
                in=phix-r1.cor.fq out=subsample-r1.fq !{reformat_opts} \
                samplebasestarget=${TOTAL_BP} \
                sampleseed=!{params.sampleseed} \
                overwrite=t 1> reformat.stdout.txt 2> reformat.stderr.txt
        else
            echo "Skipping coverage reduction"
            ln -s phix-r1.cor.fq subsample-r1.fq
            if [ "!{meta.single_end}" == "false" ]; then
                ln -s phix-r2.cor.fq subsample-r2.fq
            fi
        fi

        # Compress
        if [ "!{meta.single_end}" == "false" ]; then
            pigz -p !{task.cpus} -c -n subsample-r1.fq > results/!{meta.id}_R1.fastq.gz
            pigz -p !{task.cpus} -c -n subsample-r2.fq > results/!{meta.id}_R2.fastq.gz
        else
            pigz -p !{task.cpus} -c -n subsample-r1.fq > results/!{meta.id}.fastq.gz
        fi

        if [ "!{params.keep_all_files}" == "false" ]; then
            # Remove intermediate FASTQ files
            rm *.fq
        fi
    fi

    # Quality stats before and after QC
    fastqc -version > fastqc.version.txt 2>&1
    fastq-scan -v fastq-scan.version.txt 2>&1
    mkdir results/summary/
    if [ "!{meta.single_end}" == "false" ]; then
        # Paired-End Reads
        # fastq-scan
        gzip -cd !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > results/summary/!{meta.id}_R1-original.json
        gzip -cd !{fq[1]} | fastq-scan -g ${GENOME_SIZE} > results/summary/!{meta.id}_R2-original.json
        gzip -cd results/!{meta.id}_R1.fastq.gz | fastq-scan -g ${GENOME_SIZE} > results/summary/!{meta.id}_R1-final.json
        gzip -cd results/!{meta.id}_R2.fastq.gz | fastq-scan -g ${GENOME_SIZE} > results/summary/!{meta.id}_R2-final.json

        # FastQC
        ln -s !{fq[0]} !{meta.id}_R1-original.fastq.gz
        ln -s !{fq[1]} !{meta.id}_R2-original.fastq.gz
        ln -s results/!{meta.id}_R1.fastq.gz !{meta.id}_R1-final.fastq.gz
        ln -s results/!{meta.id}_R2.fastq.gz !{meta.id}_R2-final.fastq.gz
        fastqc --noextract -f fastq -t !{task.cpus} !{meta.id}_R1-original.fastq.gz !{meta.id}_R2-original.fastq.gz !{meta.id}_R1-final.fastq.gz !{meta.id}_R2-final.fastq.gz
    else
        # Single-End Reads
        # fastq-scan
        gzip -cd !{fq[0]} | fastq-scan -g ${GENOME_SIZE} > results/summary/!{meta.id}-original.json
        gzip -cd results/!{meta.id}.fastq.gz | fastq-scan -g ${GENOME_SIZE} > results/summary/!{meta.id}-final.json

        # FastQC 
        ln -s !{fq[0]} !{meta.id}-original.fastq.gz
        ln -s results/!{meta.id}.fastq.gz !{meta.id}-final.fastq.gz
        fastqc --noextract -f fastq -t !{task.cpus} !{meta.id}-original.fastq.gz !{meta.id}-final.fastq.gz
    fi
    mv *_fastqc.html *_fastqc.zip results/summary/

    # Final QC check
    gzip -cd results/*.fastq.gz | fastq-scan > temp.json
    FINAL_BP=$(grep "total_bp" temp.json | sed -r 's/.*:[ ]*([0-9]+),/\\1/')
    FINAL_READS=$(grep "read_total" temp.json | sed -r 's/.*:[ ]*([0-9]+),/\\1/')
    # rm temp.json

    if [ ${FINAL_BP} -lt ${MIN_BASEPAIRS} ]; then
        ERROR=1
        echo "After QC, !{meta.id} FASTQ(s) contain ${FINAL_BP} total basepairs. This does
                not exceed the required minimum ${MIN_BASEPAIRS} bp (!{params.min_coverage}x coverage). Further analysis 
                is discontinued." | \
        sed 's/^\\s*//' > !{meta.id}-low-sequence-depth-error.txt
    fi

    if [ ${FINAL_BP} -lt "!{params.min_basepairs}" ]; then
        ERROR=1
        echo "After QC, !{meta.id} FASTQ(s) contain ${FINAL_BP} total basepairs. This does
                not exceed the required minimum !{params.min_basepairs} bp. Further analysis
                is discontinued." | \
        sed 's/^\\s*//' >> !{meta.id}-low-sequence-depth-error.txt
    fi

    if [ ${FINAL_READS} -lt "!{params.min_reads}" ]; then
        ERROR=1
        echo "After QC, !{meta.id} FASTQ(s) contain ${FINAL_READS} total reads. This does
                not exceed the required minimum !{params.min_reads} reads count. Further analysis
                is discontinued." | \
        sed 's/^\\s*//' > !{meta.id}-low-read-count-error.txt
    fi

    if [ "!{is_assembly}" == "true" ]; then
        touch results/reads-simulated-from-assembly.txt
    fi

    if [ "${ERROR}" -eq "1" ]; then
        if [ "!{meta.single_end}" == "false" ]; then
            mv results/!{meta.id}_R1.fastq.gz results/!{meta.id}_R1.error-fastq.gz
            mv results/!{meta.id}_R2.fastq.gz results/!{meta.id}_R2.error-fastq.gz
        else
            mv results/!{meta.id}.fastq.gz results/!{meta.id}.error-fastq.gz
        fi
    fi
    '''

    stub:
    """
    mkdir results
    touch ${meta.id}-error.txt
    touch results/${meta.id}.fastq.gz
    touch results/${meta.id}.error-fastq.gz
    """
}
