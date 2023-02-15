// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES      = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options        = initOptions(params.options ? params.options : [:], 'qc')
options.ignore = ['.fna.gz']
options.btype  = options.btype ?: "main"
conda_tools    = "bioconda::bactopia-qc=1.0.0"
conda_name     = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env      = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process QC {
    tag "${meta.id}"
    label "base_mem_4gb"
    label "qc_reads"

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bactopia-qc:1.0.0--hdfd78af_0' :
        'quay.io/biocontainers/bactopia-qc:1.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fq), path(extra)
    path adapters
    path phix

    output:
    tuple val(meta), path("results/${prefix}*.fastq.gz"), emit: fastq, optional: true
    tuple val(meta), path("results/${prefix}*.fastq.gz"), path(extra), emit: fastq_assembly, optional: true
    path "results/*"
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions
    path "*-error.txt", optional: true

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    meta.single_end = fq[1] == null ? true : false
    is_assembly = meta.runtype.startsWith('assembly') ? true : false
    qin = meta.runtype.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    adapter_file = adapters.getName() == 'EMPTY_ADAPTERS' ? 'adapters' : adapters
    phix_file = phix.getName() == 'EMPTY_PHIX' ? 'phix' : phix
    adapter_opts = meta.single_end ? "" : "in2=repair-r2.fq out2=adapter-r2.fq"
    phix_opts = meta.single_end ? "" : "in2=adapter-r2.fq out2=phix-r2.fq"
    lighter_opts = meta.single_end ? "" : "-r phix-r2.fq"
    reformat_opts = meta.single_end ? "" : "in2=filt-r2.fq out2=subsample-r2.fq"
    fastp_fqs = meta.single_end ? "" : "--in2 ${fq[1]} --out2 filt-r2.fq --detect_adapter_for_pe"

    // set Xmx to 95% of what was allocated, to avoid going over
    xmx = Math.round(task.memory.toBytes()*0.95)
    """
    mkdir -p results
    touch results/.${meta.runtype}
    ERROR=0
    MIN_COVERAGE=\$(( ${params.min_coverage}*${meta.genome_size} ))
    TOTAL_BP=\$(( ${params.coverage}*${meta.genome_size} ))

    if [[ "${params.skip_qc}" == "true" ]]; then
        echo "Sequence QC was skipped for ${prefix}" > results/${prefix}-qc-skipped.txt
        if [ "${meta.single_end}" == "false" ]; then
            # Paired-End Reads
            cp ${fq[0]} results/${prefix}_R1.fastq.gz
            cp ${fq[1]} results/${prefix}_R2.fastq.gz
        else
            # Single-End Reads
            cp ${fq[0]} results/${prefix}.fastq.gz
        fi
    else
        if [[ "${meta.runtype}" == "ont" ]]; then
            # Remove Adapters
            porechop --input ${fq[0]} ${params.porechop_opts} \
                --format fastq \
                --threads ${task.cpus} > adapter-r1.fq

            # Quality filter
            nanoq --min-len ${params.ont_minlength} \
                  --min-qual ${params.ont_minqual} \
                  --input adapter-r1.fq 1> filt-r1.fq
        elif [[ "${params.use_fastp}" == "true" ]]; then
            # QC with fastp
            mkdir -p results/summary/
            fastp \
                --in1 ${fq[0]} --out1 filt-r1.fq ${fastp_fqs} \
                --thread ${task.cpus} \
                --json results/summary/${prefix}.fastp.json \
                --html results/summary/${prefix}.fastp.html ${params.fastp_opts} 2> ${prefix}-fastp.log
        else
            # Illumina Reads

            # Validate paired-end reads if necessary
            if [[ "${meta.single_end}" == "false" ]]; then
                # Make sure paired-end reads have matching IDs
                repair.sh \
                    in=${fq[0]} \
                    in2=${fq[1]} \
                    out=repair-r1.fq \
                    out2=repair-r2.fq \
                    outs=repair-singles.fq \
                    ain=${params.ain}

                if [ ! -s repair-r1.fq ]; then
                    ERROR=1
                    echo "After validating read pairs, ${prefix} FASTQs are empty. Please check the input FASTQs.
                        Further analysis is discontinued." | \
                    sed 's/^\\s*//' >> ${prefix}-paired-match-error.txt
                fi
            else
                gunzip -c ${fq[0]} > repair-r1.fq 
            fi

            if [ "\${ERROR}" -eq "0" ]; then
                # Remove Adapters
                bbduk.sh -Xmx${xmx} \
                    in=repair-r1.fq out=adapter-r1.fq ${adapter_opts} \
                    ref=${adapter_file} \
                    k=${params.adapter_k} \
                    ktrim=${params.ktrim} \
                    mink=${params.mink} \
                    hdist=${params.hdist} \
                    tpe=${params.tpe} \
                    tbo=${params.tbo} \
                    threads=${task.cpus} \
                    ftm=${params.ftm} \
                    ${qin} ordered=t ${params.bbduk_opts}

                if [ ! -s adapter-r1.fq ]; then
                    ERROR=1
                    echo "After adapter removal, ${prefix} FASTQs are empty. Please check the input FASTQs.
                        Further analysis is discontinued." | \
                    sed 's/^\\s*//' >> ${prefix}-adapter-qc-error.txt
                fi
            fi

            if [ "\${ERROR}" -eq "0" ]; then
                # Remove PhiX
                bbduk.sh -Xmx${xmx} \
                    in=adapter-r1.fq out=phix-r1.fq ${phix_opts} \
                    ref=${phix_file} \
                    k=${params.phix_k} \
                    hdist=${params.hdist} \
                    tpe=${params.tpe} \
                    tbo=${params.tbo} \
                    qtrim=${params.qtrim} \
                    trimq=${params.trimq} \
                    minlength=${params.minlength} \
                    minavgquality=${params.maq} \
                    ${qin} qout=${params.qout} \
                    tossjunk=${params.tossjunk} \
                    threads=${task.cpus} \
                    ordered=t ${params.bbduk_opts}

                if [ ! -s phix-r1.fq ]; then
                    ERROR=1
                    echo "After PhiX removal, ${prefix} FASTQs are empty. Please check the input FASTQs.
                        Further analysis is discontinued." | \
                    sed 's/^\\s*//' >> ${prefix}-phix-qc-error.txt
                fi
            fi
        fi

        # Error Correction
        if [ "\${ERROR}" -eq "0" ]; then
            if [[ "${params.use_fastp}" == "true" ]]; then
                echo "Skipping error correction since fastp was used"
            elif [[ "${meta.runtype}" == "ont" ]]; then
                echo "Skipping error correction. Have a recommended ONT error corrector? Let me know!"
            else
                if [ "${params.skip_error_correction}" == "false" ]; then
                    lighter -od . -r phix-r1.fq ${lighter_opts} -K 31 ${meta.genome_size} -maxcor 1 -zlib 0 -t ${task.cpus}
                    mv phix-r1.cor.fq filt-r1.fq
                    if [[ "${meta.single_end}" == "false" ]]; then
                        mv phix-r2.cor.fq filt-r2.fq
                    fi
                else
                    echo "Skipping error correction"
                    ln -s phix-r1.fq filt-r1.fq
                    if [ "${meta.single_end}" == "false" ]; then
                        ln -s phix-r2.fq filt-r2.fq
                    fi
                fi
            fi
        fi

        # Reduce Coverage
        if [ "\${ERROR}" -eq "0" ]; then
            if (( \${TOTAL_BP} > 0 )); then
                if [[ "${meta.runtype}" == "ont" ]]; then
                    rasusa -i filt-r1.fq \
                        -c ${params.coverage} \
                        -g ${meta.genome_size} \
                        -s ${params.sampleseed} 1> subsample-r1.fq
                else
                    if [ -f filt-r1.fq ]; then
                        reformat.sh -Xmx${xmx} \
                            in=filt-r1.fq out=subsample-r1.fq ${reformat_opts} \
                            samplebasestarget=\${TOTAL_BP} \
                            sampleseed=${params.sampleseed} \
                            overwrite=t
                    fi
                fi
            else
                echo "Skipping coverage reduction"
                ln -s filt-r1.fq subsample-r1.fq
                if [ "${meta.single_end}" == "false" ]; then
                    ln -s filt-r2.fq subsample-r2.fq
                fi
            fi
        fi

        # Compress
        if [ "\${ERROR}" -eq "0" ]; then
            if [ "${meta.single_end}" == "false" ]; then
                pigz -p ${task.cpus} -c -n subsample-r1.fq > results/${prefix}_R1.fastq.gz
                pigz -p ${task.cpus} -c -n subsample-r2.fq > results/${prefix}_R2.fastq.gz
            else
                pigz -p ${task.cpus} -c -n subsample-r1.fq > results/${prefix}.fastq.gz
            fi

            if [ "${params.keep_all_files}" == "false" ]; then
                # Remove remaining intermediate FASTQ files
                rm *.fq
            fi
        fi
    fi

    # Quality stats before and after QC
    if [ "\${ERROR}" -eq "0" ]; then
        mkdir -p results/summary/
        # fastq-scan
        if [[ "${meta.single_end}" == "false" ]]; then
            # Paired-End Reads
            gzip -cd ${fq[0]} | fastq-scan -g ${meta.genome_size} > results/summary/${prefix}_R1-original.json
            gzip -cd ${fq[1]} | fastq-scan -g ${meta.genome_size} > results/summary/${prefix}_R2-original.json
            gzip -cd results/${prefix}_R1.fastq.gz | fastq-scan -g ${meta.genome_size} > results/summary/${prefix}_R1-final.json
            gzip -cd results/${prefix}_R2.fastq.gz | fastq-scan -g ${meta.genome_size} > results/summary/${prefix}_R2-final.json
        else
            # Single-End Reads
            gzip -cd ${fq[0]} | fastq-scan -g ${meta.genome_size} > results/summary/${prefix}-original.json
            gzip -cd results/${prefix}.fastq.gz | fastq-scan -g ${meta.genome_size} > results/summary/${prefix}-final.json
        fi

        # FastQC and NanoPlot
        if [[ "${params.skip_qc_plots}" == "false" ]]; then
            if [[ "${meta.runtype}" == "ont" ]]; then
                mkdir results/summary/${prefix}-original results/summary/${prefix}-final
                NanoPlot ${params.nanoplot_opts} \
                    --threads ${task.cpus} \
                    --fastq ${fq[0]} \
                    --outdir results/summary/${prefix}-original/ \
                    --prefix ${prefix}-original_
                cp results/summary/${prefix}-original/${prefix}-original_NanoPlot-report.html results/summary/${prefix}-original_NanoPlot-report.html
                tar -cvf - results/summary/${prefix}-original/ | pigz --best -p ${task.cpus} > results/summary/${prefix}-original_NanoPlot.tar.gz

                NanoPlot ${params.nanoplot_opts} \
                    --threads ${task.cpus} \
                    --fastq results/${prefix}.fastq.gz \
                    --outdir results/summary/${prefix}-final/ \
                    --prefix ${prefix}-final_
                cp results/summary/${prefix}-final/${prefix}-final_NanoPlot-report.html results/summary/${prefix}-final_NanoPlot-report.html
                tar -cvf - results/summary/${prefix}-final/ | pigz --best -p ${task.cpus} > results/summary/${prefix}-final_NanoPlot.tar.gz
                rm -rf results/summary/${prefix}-original/ results/summary/${prefix}-final/
            else
                if [ "${meta.single_end}" == "false" ]; then
                    # Paired-End Reads
                    ln -s ${fq[0]} ${prefix}_R1-original.fastq.gz
                    ln -s ${fq[1]} ${prefix}_R2-original.fastq.gz
                    ln -s results/${prefix}_R1.fastq.gz ${prefix}_R1-final.fastq.gz
                    ln -s results/${prefix}_R2.fastq.gz ${prefix}_R2-final.fastq.gz
                    fastqc --noextract -f fastq -t ${task.cpus} ${prefix}_R1-original.fastq.gz ${prefix}_R2-original.fastq.gz ${prefix}_R1-final.fastq.gz ${prefix}_R2-final.fastq.gz
                else
                    # Single-End Reads
                    ln -s ${fq[0]} ${prefix}-original.fastq.gz
                    ln -s results/${prefix}.fastq.gz ${prefix}-final.fastq.gz
                    fastqc --noextract -f fastq -t ${task.cpus} ${prefix}-original.fastq.gz ${prefix}-final.fastq.gz
                fi
                mv *_fastqc.html *_fastqc.zip results/summary/
            fi
        fi
    fi

    # Final QC check
    if [ "${params.skip_fastq_check}" == "false" ]; then
        # Only check for errors if we haven't already found them
        if [ "\${ERROR}" -eq "0" ]; then
            gzip -cd results/*.fastq.gz | fastq-scan -g ${meta.genome_size} > temp.json
            FINAL_BP=\$(grep "total_bp" temp.json | sed -r 's/.*:[ ]*([0-9]+),/\\1/')
            rm temp.json

            if [ \${FINAL_BP} -lt \${MIN_COVERAGE} ]; then
                ERROR=1
                echo "After QC, ${prefix} FASTQ(s) contain \${FINAL_BP} total basepairs. This does
                        not exceed the required minimum \${MIN_COVERAGE} bp (${params.min_coverage}x coverage). Further analysis 
                        is discontinued." | \
                sed 's/^\\s*//' > ${prefix}-low-sequence-coverage-error.txt
                ERROR=2
            fi

            # Check paired-end reads have same read counts
            OPTS="--sample ${prefix} --min_basepairs ${params.min_basepairs} --min_reads ${params.min_reads} --min_proportion ${params.min_proportion} --runtype ${meta.runtype}"
            if [ -f  "results/${prefix}_R2.fastq.gz" ]; then
                # Paired-end
                gzip -cd results/${prefix}_R1.fastq.gz | fastq-scan > r1.json
                gzip -cd results/${prefix}_R2.fastq.gz | fastq-scan > r2.json
                if ! reformat.sh in1=results/${prefix}_R1.fastq.gz in2=results/${prefix}_R2.fastq.gz ${qin} out=/dev/null 2> ${prefix}-paired-end-error.txt; then
                    ERROR=2
                    echo "${prefix} FASTQs contains an error. Please check the input FASTQs.
                        Further analysis is discontinued." | \
                    sed 's/^\\s*//' >> ${prefix}-paired-end-error.txt
                else
                    rm -f ${prefix}-paired-end-error.txt
                fi
                if ! check-fastqs.py --fq1 r1.json --fq2 r2.json \${OPTS}; then
                    ERROR=2
                fi
                rm r1.json r2.json
            else
                # Single-end
                gzip -cd results/${prefix}.fastq.gz | fastq-scan > r1.json
                if ! check-fastqs.py --fq1 r1.json \${OPTS}; then
                    ERROR=2
                fi
                rm r1.json
            fi
        fi
    fi

    if [ "${is_assembly}" == "true" ]; then
        touch results/reads-simulated-from-assembly.txt
    fi

    if [ "\${ERROR}" -eq "1" ]; then
        if [ "${meta.single_end}" == "false" ]; then
            cp ${fq[0]} results/${prefix}_R1.error-fastq.gz
            cp ${fq[1]} results/${prefix}_R2.error-fastq.gz
            if [ ! -s repair-singles.fq ]; then
                pigz -p ${task.cpus} -c -n repair-singles.fq > results/${prefix}.error-fastq.gz
            fi
        else
            cp ${fq[0]} results/${prefix}.error-fastq.gz
        fi
    elif [ "\${ERROR}" -eq "2" ]; then
        if [ "${meta.single_end}" == "false" ]; then
            if [ -f results/${prefix}_R1.fastq.gz ]; then
                mv results/${prefix}_R1.fastq.gz results/${prefix}_R1.error-fastq.gz
                mv results/${prefix}_R2.fastq.gz results/${prefix}_R2.error-fastq.gz

                if [ -s repair-singles.fq ]; then
                    pigz -p ${task.cpus} -c -n repair-singles.fq > results/${prefix}.error-fastq.gz
                fi
            fi
        else
            if [ -f results/${prefix}.fastq.gz ]; then
                mv results/${prefix}.fastq.gz results/${prefix}.error-fastq.gz
            fi
        fi
    fi

    # Capture versions
    if [[ "${params.skip_qc_plots}" == "false" ]]; then
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbduk: \$(echo \$(bbduk.sh --version 2>&1) | sed 's/^.*BBMap version //;s/ .*\$//')
        fastp: \$(echo \$(fastp --version 2>&1) | sed -e "s/fastp //g")
        fastqc: \$(echo \$(fastqc --version 2>&1) | sed 's/^.*FastQC v//')
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        lighter: \$(echo \$(lighter -v 2>&1) | sed 's/Lighter v//')
        nanoplot: \$(echo \$(NanoPlot -v 2>&1) | sed 's/NanoPlot //')
        nanoq: \$(echo \$(nanoq --version 2>&1) | sed 's/nanoq //')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        porechop: \$(echo \$(porechop --version 2>&1))
        rasusa: \$(echo \$(rasusa --version 2>&1) | sed 's/rasusa //')
    END_VERSIONS
    else
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbduk: \$(echo \$(bbduk.sh --version 2>&1) | sed 's/^.*BBMap version //;s/ .*\$//')
        fastp: \$(echo \$(fastp --version 2>&1) | sed -e "s/fastp //g")
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        lighter: \$(echo \$(lighter -v 2>&1) | sed 's/Lighter v//')
        nanoq: \$(echo \$(nanoq --version 2>&1) | sed 's/nanoq //')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        porechop: \$(echo \$(porechop --version 2>&1))
        rasusa: \$(echo \$(rasusa --version 2>&1) | sed 's/rasusa //')
    END_VERSIONS
    fi
    """
}
