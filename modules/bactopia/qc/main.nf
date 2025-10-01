process QC {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fq), path(extra)
    path adapters
    path phix

    output:
    tuple val(meta), path("results/${prefix}*.fastq.gz"), path("extra/*"), emit: fastq, optional: true
    tuple val(meta), path("results/${prefix}*.fastq.gz")                 , emit: fastq_only, optional: true
    tuple val(meta), path("*-error.txt")   , emit: error, optional: true
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    meta.single_end = fq[1] == null ? true : false
    meta.runtype = _meta.runtype
    is_assembly = meta.runtype.startsWith('assembly') ? true : false
    qin = meta.runtype.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    adapter_file = adapters.getName() == 'EMPTY_ADAPTERS' ? 'adapters' : adapters
    phix_file = phix.getName() == 'EMPTY_PHIX' ? 'phix' : phix
    adapter_opts = meta.single_end ? "" : "in2=repair-r2.fq out2=adapter-r2.fq"
    phix_opts = meta.single_end ? "" : "in2=adapter-r2.fq out2=phix-r2.fq"
    lighter_opts = meta.single_end ? "" : "-r phix-r2.fq"
    reformat_opts = meta.single_end ? "" : "in2=filt-r2.fq out2=subsample-r2.fq"
    fastp_fqs = meta.single_end ? "" : "--in2 ${fq[1]} --out2 filt-r2.fq --detect_adapter_for_pe"
    ont_fq = meta.runtype == 'ont' ? fq[0] : extra[0]
    meta.single_end = meta.runtype == 'short_polish' ? true : meta.single_end

    // set Xmx to 95% of what was allocated, to avoid going over
    xmx = Math.round(task.memory.toBytes()*0.95)
    """
    mkdir -p results
    touch results/.${meta.runtype}
    ERROR=0
    MIN_COVERAGE=\$(( ${task.ext.min_coverage}*${meta.genome_size} ))
    TOTAL_BP=\$(( ${task.ext.coverage}*${meta.genome_size} ))
    ENABLE_ONT=0
    ENABLE_ILLUMINA=0
    IS_HYBRID=0

    if [ "${meta.runtype}" == "hybrid" ] || [ "${meta.runtype}" == "short_polish" ]; then
        ENABLE_ONT=1
        ENABLE_ILLUMINA=1
        IS_HYBRID=1
    elif [ "${meta.runtype}" == "ont" ]; then
        ENABLE_ONT=1
    else
        ENABLE_ILLUMINA=1
    fi

    mkdir extra/
    if [[ "${is_assembly}" == "true" ]]; then
        # Copy the assembly over to extra
        cp ${extra[0]} extra/
    else
        touch extra/EMPTY_EXTRA
    fi

    if [[ "${task.ext.skip_qc}" == "true" ]]; then
        echo "Sequence QC was skipped for ${prefix}" > results/${prefix}-qc-skipped.txt
        if [ "\${IS_HYBRID}" -eq "1" ]; then
            # Illumina Paired-End Reads and Nanopore Reads
            cp ${fq[0]} results/${prefix}_R1.fastq.gz
            cp ${fq[1]} results/${prefix}_R2.fastq.gz
            cp ${extra[0]} results/${prefix}.fastq.gz
        elif [ "${meta.single_end}" == "false" ]; then
            # Paired-End Reads
            cp ${fq[0]} results/${prefix}_R1.fastq.gz
            cp ${fq[1]} results/${prefix}_R2.fastq.gz
        else
            # Single-End Reads
            cp ${fq[0]} results/${prefix}.fastq.gz
        fi
    else
        if [ "\${ENABLE_ONT}" -eq "1" ]; then
            # QC the Nanopore reads
            if [[ "${task.ext.use_porechop}" == "true" ]]; then
                # Remove Adapters
                porechop --input ${ont_fq} ${task.ext.porechop_opts} \
                    --format fastq \
                    --threads ${task.cpus} > adapter-ont.fq
            else
                cp ${ont_fq} adapter-ont.fq
            fi

            # Filter reads based on length and quality
            nanoq --input adapter-ont.fq \
                --min-qual ${task.ext.ont_min_qual} \
                --min-len ${task.ext.ont_min_len} \
                --max-len ${task.ext.ont_max_len} \
                --report nanoq.json | \
            rasusa --input /dev/stdin \
                --depth \${TOTAL_BP} \
                --genome-size ${meta.genome_size} \
                --seed ${task.ext.sampleseed} > results/${prefix}.fastq

            # Compress
            pigz -c -p ${task.cpus} -n results/${prefix}.fastq > results/${prefix}.fastq.gz
            rm results/${prefix}.fastq

            # Create a symlink to results with _SE extension for bactopia-tools
            cd results
            ln -s ${prefix}.fastq.gz ${prefix}_SE.fastq.gz
            cd ../

            TOTAL_BP_ONT=\$(grep -H "total_bp" nanoq.json | sed -r 's/.*:([0-9]+),?/\\1/')
            if [ -s results/${prefix}.fastq.gz ]; then
                FILTERED_BP_ONT=\$(gzip -cd results/${prefix}.fastq.gz | fastq-scan | grep -H "total_bp" | sed -r 's/.*:([0-9]+),?/\\1/')
            else
                FILTERED_BP_ONT=0
            fi
            rm adapter-ont.fq

            # Generate Reports
            if [[ "${task.ext.skip_qc_plots}" == "false" ]]; then
                NanoPlot --threads ${task.cpus} \
                    --fastq ${ont_fq} \
                    --N50 \
                    --title "${meta.id} ${meta.runtype} Raw QC" \
                    --color ${task.ext.nanoplot_color} \
                    -f ${task.ext.nanoplot_format} \
                    --plots ${task.ext.nanoplot_plots} \
                    -o ont-original

                NanoPlot --threads ${task.cpus} \
                    --fastq results/${prefix}.fastq.gz \
                    --N50 \
                    --title "${meta.id} ${meta.runtype} Filtered QC" \
                    --color ${task.ext.nanoplot_color} \
                    -f ${task.ext.nanoplot_format} \
                    --plots ${task.ext.nanoplot_plots} \
                    -o ont-final

                # Cleanup
                rm *.fastq.gz 2> /dev/null || true
                mv ont-original results/
                mv ont-final results/
            fi

            if [ "\${IS_HYBRID}" -eq "0" ]; then
                if [ "\${FILTERED_BP_ONT}" -lt "\${MIN_COVERAGE}" ]; then
                    ERROR=1
                    echo "${prefix} ONT estimated coverage (\${FILTERED_BP_ONT} bp) does not exceed the minimum 
                        coverage required (\${MIN_COVERAGE} bp). If this is unexpected, please
                        investigate ${prefix} to determine a cause (e.g. metagenomic, contaminants, etc...) 
                        for the lack of coverage." | \
                    sed 's/^\\s*//' > ${prefix}-low-ont-coverage-error.txt
                fi
            fi
        fi

        if [ "\${ENABLE_ILLUMINA}" -eq "1" ]; then
            # Setup BBDuk preprocessing for fastp and/or Lighter
            if [ "${task.ext.use_bbmap}" == "true" ]; then
                # Randomly Subsample to 2x the expected coverage (prevent large fastq manipulation)
                reformat.sh -Xmx${xmx} \
                    in=${fq[0]} out=subsample-r1.fq samplebasestarget=\$((\${TOTAL_BP}*2)) \
                    ${qin} interleaved=f overwrite=t ${reformat_opts} idmodulo=1

                # Correct reads to be in the proper order
                repair.sh -Xmx${xmx} \
                    in1=subsample-r1.fq out1=repair-r1.fq outs=repair-singles.fq ${adapter_opts} overwrite=t 

                # Remove Adapters
                bbduk.sh -Xmx${xmx} \
                    in=repair-r1.fq out=adapter-r1.fq ref=${adapter_file} \
                    stats=results/bbduk-adapter.txt threads=${task.cpus} overwrite=t \
                    ktrim=r mink=${task.ext.mink} hdist=${task.ext.hdist} k=${task.ext.adapter_k} tpe=${task.ext.tpe} tbo=${task.ext.tbo} ${adapter_opts}

                # Remove PhiX
                bbduk.sh -Xmx${xmx} \
                    in=adapter-r1.fq out=phix-r1.fq ref=${phix_file} \
                    stats=results/bbduk-phix.txt threads=${task.cpus} overwrite=t \
                    k=${task.ext.phix_k} hdist=${task.ext.hdist} ${phix_opts}

                # Lighter Error Correction
                lighter -t ${task.cpus} -r phix-r1.fq -k 31 ${meta.genome_size} 0.1 ${lighter_opts}
                LIGHTER_COMPLETE=\$?
                rm -rf subsample-r1.fq subsample-r2.fq repair-r1.fq repair-r2.fq
                rm -rf adapter-r1.fq adapter-r2.fq phix-r1.fq phix-r2.fq

                if [ "\${LIGHTER_COMPLETE}" -eq 0 ]; then
                    if [ "${meta.single_end}" == "false" ]; then
                        # Paired-End Reads
                        mv phix-r1.cor.fq filt-r1.fq
                        mv phix-r2.cor.fq filt-r2.fq
                    else
                        # Single-End Reads
                        mv phix-r1.cor.fq filt-r1.fq
                    fi
                else
                    # Lighter returned an error
                    ERROR=1
                    echo "${prefix} had too few reads for Lighter error correction" | \
                    sed 's/^\\s*//' > ${prefix}-low-read-count-error.txt
                fi
            else
                # Skip QC preprocessing
                if [ "${meta.single_end}" == "false" ]; then
                    # Paired-End Reads
                    gzip -cd ${fq[0]} > filt-r1.fq
                    gzip -cd ${fq[1]} > filt-r2.fq
                else
                    # Single-End Reads
                    gzip -cd ${fq[0]} > filt-r1.fq
                fi
            fi

            if [ "\${ERROR}" -eq "0" ]; then
                # fastp adapter trimming and quality filtering
                fastp --in1 filt-r1.fq ${fastp_fqs} \
                    --out1 final-r1.fq \
                    --unpaired1 final-unpaired.fq --unpaired2 final-unpaired.fq \
                    --failed_out ${prefix}-error-reads.fq \
                    ${task.ext.args} \
                    --thread ${task.cpus} \
                    --json results/${prefix}-fastp.json \
                    --html results/${prefix}-fastp.html

                rm filt-r1.fq

                # Failed reads
                FAILED_READS=\$(head -n 1 ${prefix}-error-reads.fq | wc -l)
                if [ "\${FAILED_READS}" -gt 0 ]; then
                    ERROR=2
                    echo "${prefix} has reads that failed fastp. Please review the read set to
                        determine a cause (e.g. adapters, below --qualified_quality_phred, etc...) 
                        for the failure. The failed reads have been written to the logs/ folder." | \
                    sed 's/^\\s*//' > ${prefix}-fastp-error.txt
                    pigz -p ${task.cpus} -c -n ${prefix}-error-reads.fq > ${prefix}-error-reads.fq.gz
                fi
                rm ${prefix}-error-reads.fq
                
                # Sub-sampling 
                reformat.sh -Xmx${xmx} \
                    in=final-r1.fq out=subsample-r1.fq samplebasestarget=\${TOTAL_BP} seed=${task.ext.sampleseed} \
                    ${qin} interleaved=f overwrite=t ${reformat_opts}
                rm final-r1.fq

                if [ "${meta.single_end}" == "false" ]; then
                    rm filt-r2.fq final-r2.fq
                    mv subsample-r1.fq results/${prefix}_R1.fastq
                    mv subsample-r2.fq results/${prefix}_R2.fastq
                    pigz -p ${task.cpus} -c -n results/${prefix}_R1.fastq > results/${prefix}_R1.fastq.gz
                    pigz -p ${task.cpus} -c -n results/${prefix}_R2.fastq > results/${prefix}_R2.fastq.gz
                    rm results/${prefix}_R1.fastq results/${prefix}_R2.fastq

                    # Total up reads
                    gzip -cd results/${prefix}_R1.fastq.gz | fastq-scan > r1.json
                    gzip -cd results/${prefix}_R2.fastq.gz | fastq-scan > r2.json
                    FINAL_BP=\$(jq '.qc_stats.read_total * .qc_stats.read_mean' r1.json r2.json | awk '{s+=$1} END {printf "%.0f", s}')
                    FINAL_READS=\$(jq '.qc_stats.read_total' r1.json r2.json | awk '{s+=$1} END {printf "%d", s}')
                    rm r1.json r2.json
                else
                    mv subsample-r1.fq results/${prefix}.fastq
                    pigz -p ${task.cpus} -c -n results/${prefix}.fastq > results/${prefix}.fastq.gz
                    rm results/${prefix}.fastq

                    # Total up reads
                    gzip -cd results/${prefix}.fastq.gz | fastq-scan > r1.json
                    FINAL_BP=\$(jq '.qc_stats.read_total * .qc_stats.read_mean' r1.json | awk '{printf "%.0f", $1}')
                    FINAL_READS=\$(jq '.qc_stats.read_total' r1.json | awk '{printf "%d", $1}')
                    rm r1.json
                fi

                if [ "\${FINAL_BP}" -lt "\${MIN_COVERAGE}" ]; then
                    ERROR=1
                    echo "${prefix} Illumina estimated coverage (\${FINAL_BP} bp) does not exceed the minimum 
                        coverage required (\${MIN_COVERAGE} bp). If this is unexpected, please
                        investigate ${prefix} to determine a cause (e.g. metagenomic, contaminants, etc...) 
                        for the lack of coverage." | \
                    sed 's/^\\s*//' > ${prefix}-low-illumina-coverage-error.txt
                fi

                # Generate Reports
                if [[ "${task.ext.skip_qc_plots}" == "false" ]]; then
                    if [ "${meta.single_end}" == "false" ]; then
                        # Paired-End Reads
                        ln -s ${fq[0]} ${prefix}_R1.fastq.gz
                        ln -s ${fq[1]} ${prefix}_R2.fastq.gz
                        fastqc --noextract -f fastq -t ${task.cpus} ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz
                        mv ${prefix}_R1_fastqc.html results/original-r1_fastqc.html
                        mv ${prefix}_R1_fastqc.zip results/original-r1_fastqc.zip
                        mv ${prefix}_R2_fastqc.html results/original-r2_fastqc.html
                        mv ${prefix}_R2_fastqc.zip results/original-r2_fastqc.zip

                        # Post QC
                        fastqc --noextract -f fastq -t ${task.cpus} results/${prefix}_R1.fastq.gz results/${prefix}_R2.fastq.gz
                        mv results/${prefix}_R1_fastqc.html results/final-r1_fastqc.html
                        mv results/${prefix}_R1_fastqc.zip results/final-r1_fastqc.zip
                        mv results/${prefix}_R2_fastqc.html results/final-r2_fastqc.html
                        mv results/${prefix}_R2_fastqc.zip results/final-r2_fastqc.zip
                    else
                        # Single-End Reads
                        ln -s ${fq[0]} ${prefix}.fastq.gz
                        fastqc --noextract -f fastq -t ${task.cpus} ${prefix}.fastq.gz

                        # Pre QC
                        mv ${prefix}_fastqc.html results/original_fastqc.html
                        mv ${prefix}_fastqc.zip results/original_fastqc.zip

                        # Post QC
                        fastqc --noextract -f fastq -t ${task.cpus} results/${prefix}.fastq.gz
                        mv results/${prefix}_fastqc.html results/final_fastqc.html
                        mv results/${prefix}_fastqc.zip results/final_fastqc.zip
                    fi

                    # Cleanup
                    rm *.fastq.gz
                fi

                # Only Paired-End Reads have a separate unpaired FASTQ
                if [ -s final-unpaired.fq ]; then
                    pigz -p ${task.cpus} -c -n final-unpaired.fq > results/${prefix}-unpaired.fastq.gz
                fi
                rm final-unpaired.fq

                # Store extra reads (e.g. repair.sh singletons, final.unpaired.fq)
                EXTRA_READS=\$(head -n 1 repair-singles.fq | wc -l)
                if [[ "\${EXTRA_READS}" -gt 0 ]]; then
                    cat repair-singles.fq >> extra-reads.fq
                    pigz -p ${task.cpus} -c -n extra-reads.fq > results/${prefix}-extra.fastq.gz
                fi
                rm -f repair-singles.fq extra-reads.fq
            fi
        fi
    fi

    if [ "\${ERROR}" -eq "1" ]; then
        if [ "\${IS_HYBRID}" -eq "1" ]; then
            cp ${fq[0]} results/${prefix}_R1.error-fastq.gz
            cp ${fq[1]} results/${prefix}_R2.error-fastq.gz
            cp ${extra[0]} results/${prefix}.error-fastq.gz
            if [ ! -s repair-singles.fq ]; then
                pigz -p ${task.cpus} -c -n repair-singles.fq > results/${prefix}.error-fastq.gz
            fi
        elif [ "${meta.single_end}" == "false" ]; then
            cp ${fq[0]} results/${prefix}_R1.error-fastq.gz
            cp ${fq[1]} results/${prefix}_R2.error-fastq.gz
            if [ ! -s repair-singles.fq ]; then
                pigz -p ${task.cpus} -c -n repair-singles.fq > results/${prefix}.singles-fastq.gz
            fi
        else
            cp ${fq[0]} results/${prefix}.error-fastq.gz
        fi
    elif [ "\${ERROR}" -eq "2" ]; then
        if [ -f results/${prefix}_R1.fastq.gz ]; then
            mv results/${prefix}_R1.fastq.gz results/${prefix}_R1.error-fastq.gz
            mv results/${prefix}_R2.fastq.gz results/${prefix}_R2.error-fastq.gz

            if [ -s repair-singles.fq ]; then
                pigz -p ${task.cpus} -c -n repair-singles.fq > results/${prefix}.singles-fastq.gz
            fi
        fi

        if [ -f results/${prefix}.fastq.gz ]; then
            mv results/${prefix}.fastq.gz results/${prefix}.error-fastq.gz
        fi
    fi

    # Capture versions
    if [[ "${task.ext.skip_qc_plots}" == "false" ]]; then
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbduk: \$(echo \$(bbduk.sh --version 2>&1) | sed 's/^.*BBTools version //;s/ .*\$//')
        fastp: \$(echo \$(fastp --version 2>&1) | sed -e "s/fastp //g")
        fastqc: \$(echo \$(fastqc --version 2>&1) | sed 's/^.*FastQC v//')
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        lighter: \$(echo \$(lighter -v 2>&1) | sed 's/Lighter v//')
        nanoplot: \$(echo \$(NanoPlot -v 2>&1) | sed 's/.*NanoPlot //')
        nanoq: \$(echo \$(nanoq --version 2>&1) | sed 's/nanoq //')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        porechop: \$(echo \$(porechop --version 2>&1))
        rasusa: \$(echo \$(rasusa --version 2>&1) | sed 's/rasusa //')
    END_VERSIONS
    else
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbduk: \$(echo \$(bbduk.sh --version 2>&1) | sed 's/^.*BBTools version //;s/ .*\$//')
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
