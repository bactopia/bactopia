/**
 * Search, validate, gather, or simulate input samples.
 *
 * This process is the entry point for data ingestion. It handles:
 * - **Validation:** Verifies FASTQ formatting and gzip integrity.
 * - **Merging:** Combines multiple runs (lanes) into a single sample.
 * - **Downloading:** Fetches reads (SRA/ENA) or assemblies (NCBI) from accessions.
 * - **Simulation:** Generates synthetic reads from assemblies using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) to enable read-based analysis.
 *
 * Uses explicit positional tuple slots for input and output reads:
 * - Input accepts Set<Path> for each slot (pre-merge, supports multiple files)
 * - Output emits Path? for each slot (post-merge, single consolidated file or null)
 *
 * @status stable
 * @keywords fastq, validation, sra, ena, download, merging, simulation, art, ncbi
 * @tags complexity:complex input-type:multiple output-type:multiple features:internet-access,resource-download,conditional-logic
 * @citation bactopia, art, fastq_dl, fastq_scan, ncbigenomedownload, pigz
 *
 * @input tuple(meta, r1_files, r2_files, se_files, lr_files)
 * - `meta`: Groovy Map containing sample information
 * - `r1_files`: Illumina R1 read files (Set, for merging multiple runs)
 * - `r2_files`: Illumina R2 read files (Set, for merging multiple runs)
 * - `se_files`: Single-end read files (Set, for merging multiple runs)
 * - `lr_files`: Long read files (ONT/PacBio) or assembly for simulation
 *
 * @output reads       A tuple with explicit read slots: (meta, r1, r2, se, lr) where each is Path?
 * @output tsv         A tab-delimited metadata file describing the valid samples
 * @output error       Captured error messages for validation or download failures
 * @output logs        Optional software execution logs containing warnings/errors
 * @output nf_logs     Nextflow execution scripts and logs for debugging
 * @output versions    A YAML formatted file with software versions
 */
nextflow.preview.types = true

process GATHER {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, r1_files, r2_files, se_files, lr_files) : Tuple<Map, Set<Path>, Set<Path>, Set<Path>, Set<Path>>

    stage:
    stageAs '*???-r1', r1_files
    stageAs '*???-r2', r2_files
    stageAs '*???-se', se_files
    stageAs '*???-lr', lr_files

    output:
    reads    = tuple(meta, file("fastqs/${prefix}_R1.fastq.gz", optional: true), file("fastqs/${prefix}_R2.fastq.gz", optional: true), file("fastqs/${prefix}.fastq.gz", optional: true), file("extra/${prefix}.fastq.gz", optional: true))
    tsv      = tuple(meta, files("${prefix}-meta.tsv"))
    error    = tuple(meta, files("*-{error,merged}.txt", optional: true))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    if (task.ext.wf == "teton") {
        meta.output_dir = "${prefix}/teton/main/${task.ext.process_name}/${task.ext.subdir}"
        meta.logs_dir = "${prefix}/teton/main/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    }
    else {
        meta.output_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}"
        meta.logs_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    }
    meta.process_name = task.ext.process_name
    meta.original_runtype = _meta.runtype
    meta.genome_size = _meta.genome_size
    meta.species = _meta.species

    // WF specific parameters
    runtype = meta.original_runtype
    is_assembly = runtype.startsWith('assembly') ? true : false

    // Get first file from lr_files for assembly/hybrid cases
    def lr_file = lr_files.size() > 0 ? lr_files.toList()[0] : null
    is_compressed = lr_file ? (lr_file.getName().endsWith('gz') ? true : false) : false

    no_cache = task.ext.no_cache ? '-N' : ''
    archive = task.ext.use_ena ? (task.attempt >= 4 ? "SRA" : "ENA") : "SRA"
    section = runtype == 'assembly_accession' ? (prefix.startsWith('GCF') ? 'refseq' : 'genbank') : null
    fcov = task.ext.coverage.toInteger() == 0 ? 150 : Math.round(task.ext.coverage.toInteger() * 1.5)

    // Normalize runtype for downstream
    if (runtype == 'hybrid-merge-pe') {
        meta.runtype = 'hybrid'
    }
    else if (runtype == 'short_polish-merge-pe') {
        meta.runtype = 'short_polish'
    }
    else if (runtype == 'merge-pe') {
        meta.runtype = 'paired-end'
    }
    else if (runtype == 'merge-se') {
        meta.runtype = 'single-end'
    }
    else if (runtype == 'sra_accession_ont') {
        meta.runtype = 'ont'
    }
    else {
        meta.runtype = runtype
    }
    meta.is_compressed = task.ext.skip_compression ? false : true

    // Determine what reads we have based on the explicit slots
    has_r1 = r1_files.size() > 0
    has_r2 = r2_files.size() > 0
    has_se = se_files.size() > 0
    has_lr = lr_files.size() > 0

    // For paired-end, get the first file from each set (single file case) or handle merge
    r1_first = has_r1 ? r1_files.toList()[0] : null
    r2_first = has_r2 ? r2_files.toList()[0] : null
    se_first = has_se ? se_files.toList()[0] : null
    lr_first = has_lr ? lr_files.toList()[0] : null

    qin = is_assembly ? 'qin=33' : 'qin=auto'
    """
    MERGED="multiple-read-sets-merged.txt"
    mkdir -p fastqs
    mkdir -p extra

    if [ "${runtype}" == "paired-end" ]; then
        # Paired-End Reads
        cp -L ${r1_first} fastqs/${prefix}_R1.fastq.gz
        cp -L ${r2_first} fastqs/${prefix}_R2.fastq.gz
    elif [ "${runtype}" == "single-end" ]; then
        # Single-End Reads
        cp -L ${se_first} fastqs/${prefix}.fastq.gz
    elif [ "${runtype}" == "ont" ] || [ "${runtype}" == "pacbio" ]; then
        # Long reads (Nanopore or PacBio)
        cp -L ${lr_first} extra/${prefix}.fastq.gz
    elif  [ "${runtype}" == "hybrid" ] || [ "${runtype}" == "short_polish" ]; then
        # Paired-End Reads + Long Reads
        cp -L ${r1_first} fastqs/${prefix}_R1.fastq.gz
        cp -L ${r2_first} fastqs/${prefix}_R2.fastq.gz
        cp -L ${lr_first} extra/${prefix}.fastq.gz
    elif [ "${runtype}" == "merge-pe" ] || [ "${runtype}" == "hybrid-merge-pe" ] || [ "${runtype}" == "short_polish-merge-pe" ]; then
        # Merge Paired-End Reads
        echo "This sample had reads merged." > \${MERGED}
        echo "R1:" >> \${MERGED}
        find -name "*-r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print \$5"\t"\$9}' >> \${MERGED}
        find -name "*-r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/${prefix}_R1.fastq.gz
        echo "Merged R1:" >> \${MERGED}
        ls -l fastqs/${prefix}_R1.fastq.gz | awk '{print \$5"\t"\$9}' >> \${MERGED}

        echo "R2:" >> \${MERGED}
        find -name "*-r2" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print \$5"\t"\$9}' >> \${MERGED}
        find -name "*-r2" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/${prefix}_R2.fastq.gz
        echo "Merged R2:" >> \${MERGED}
        ls -l fastqs/${prefix}_R2.fastq.gz | awk '{print \$5"\t"\$9}' >> \${MERGED}

        if [ "${runtype}" == "hybrid-merge-pe" ] || [ "${runtype}" == "short_polish-merge-pe" ]; then
            cp -L ${lr_first} extra/${prefix}.fastq.gz
        fi
    elif [ "${runtype}" == "merge-se" ]; then
        # Merge Single-End Reads
        echo "This sample had reads merged." > \${MERGED}
        echo "SE:" >> \${MERGED}
        find -name "*-se" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print \$5"\t"\$9}' >> \${MERGED}
        find -name "*-se" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/${prefix}.fastq.gz
        echo "Merged SE:" >> \${MERGED}
        ls -l fastqs/${prefix}.fastq.gz | awk '{print \$5"\t"\$9}' >> \${MERGED}
    elif [ "${runtype}" == "sra_accession" ] || [ "${runtype}" == "sra_accession_ont" ]; then
        if [ "${task.attempt}" == "${task.ext.max_retry}" ]; then
            echo "Unable to download ${prefix} from both SRA and ENA ${task.ext.max_retry} times. This may or may
                not be a temporary connection issue. Rather than stop the whole Bactopia run,
                further analysis of ${prefix} will be discontinued." | \\
            sed 's/^\\s*//' > ${prefix}-fastq-download-error.txt
            exit
        else
            # Download accession from ENA/SRA
            fastq-dl \\
                --accession ${prefix} \\
                --provider ${archive} \\
                --cpus ${task.cpus} \\
                --outdir fastqs/ \\
                --prefix ${prefix} \\
                --group-by-experiment

            # Move ONT reads to the extra folder (lr slot)
            if [ "${runtype}" == "sra_accession_ont" ]; then
                if [ -f "fastqs/${prefix}.fastq.gz" ]; then
                    mv fastqs/${prefix}.fastq.gz extra/${prefix}.fastq.gz
                fi
            fi
        fi
    elif [ "${is_assembly}" == "true" ]; then
        if [ "${runtype}" == "assembly_accession" ]; then
            if [ "${task.attempt}" == "${task.ext.max_retry}" ]; then
                echo "Unable to download ${prefix} from NCBI Assembly ${task.ext.max_retry} times. This may or may
                    not be a temporary connection issue. Rather than stop the whole Bactopia run,
                    further analysis of ${prefix} will be discontinued." | \\
                sed 's/^\\s*//' > ${prefix}-assembly-download-error.txt
                exit
            else
                # Verify Assembly accession
                check-assembly-accession.py ${prefix} > accession.txt 2> check-assembly-accession.txt

                if [ -s "accession.txt" ]; then
                    # Download from NCBI assembly and simulate reads
                    mkdir fasta/
                    ncbi-genome-download bacteria -o ./ -F fasta -p ${task.cpus} \\
                                                -u "https://ftp.ncbi.nlm.nih.gov/genomes" \\
                                                -s ${section} -A accession.txt -r 50 ${no_cache}
                    find . -name "*${prefix}*.fna.gz" | xargs -I {} mv {} fasta/
                    rename 's/(GC[AF]_\\d+).*/\$1.fna.gz/' fasta/*
                    gzip -cd fasta/${prefix}.fna.gz > ${prefix}-art.fna
                    rm check-assembly-accession.txt
                else
                    mv check-assembly-accession.txt ${prefix}-assembly-accession-error.txt
                    exit
                fi
            fi
        elif [ "${runtype}" == "assembly" ]; then
            if [ "${is_compressed}" == "true" ]; then
                gzip -cd ${lr_first} > ${prefix}-art.fna
            else
                cat ${lr_first} > ${prefix}-art.fna
            fi
        fi

        # Simulate reads from assembly, reads are 250bp without errors
        art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov ${fcov} -ir 0 -ir2 0 -dr 0 -dr2 0 -rs ${task.ext.sampleseed}\\
                     -na -qL 33 -qU 40 -o ${prefix}_R --id ${prefix} -i ${prefix}-art.fna

        mv ${prefix}_R1.fq fastqs/${prefix}_R1.fastq
        mv ${prefix}_R2.fq fastqs/${prefix}_R2.fastq
        pigz -p ${task.cpus} --fast fastqs/*.fastq
        # Keep the original assembly in extra/lr slot
        cp ${prefix}-art.fna extra/${prefix}.fna
        pigz -p ${task.cpus} --best extra/${prefix}.fna
        # Move to lr slot format
        mv extra/${prefix}.fna.gz extra/${prefix}.fastq.gz || true
    fi

    # Validate input FASTQs
    if [ "${task.ext.skip_fastq_check}" == "false" ]; then
        ERROR=0
        # Check paired-end reads have same read counts
        OPTS="--sample ${prefix} --min_basepairs ${task.ext.min_basepairs} --min_reads ${task.ext.min_reads} --min_proportion ${task.ext.min_proportion} --runtype ${runtype}"
        if [ -f  "fastqs/${prefix}_R2.fastq.gz" ]; then
            # Paired-end
            IS_PAIRED="true"
            gzip -cd fastqs/${prefix}_R1.fastq.gz | fastq-scan > r1.json
            gzip -cd fastqs/${prefix}_R2.fastq.gz | fastq-scan > r2.json
            if ! reformat.sh in1=fastqs/${prefix}_R1.fastq.gz in2=fastqs/${prefix}_R2.fastq.gz ${qin} out=/dev/null 2> ${prefix}-paired-end-error.txt; then
                ERROR=1
                echo "${prefix} FASTQs contains an error. Please check the input FASTQs.
                    Further analysis is discontinued." | \\
                sed 's/^\\s*//' >> ${prefix}-paired-end-error.txt
            else
                rm -f ${prefix}-paired-end-error.txt
            fi

            if [[ -s r1.json ]] && [[ -s r2.json ]]; then
                if ! check-fastqs.py --fq1 r1.json --fq2 r2.json \${OPTS}; then
                    ERROR=1
                fi
            else
                NOT_GZIP=0
                if ! gzip -t fastqs/${prefix}_R1.fastq.gz; then
                    NOT_GZIP=1
                elif ! gzip -t fastqs/${prefix}_R2.fastq.gz; then
                    NOT_GZIP=1
                fi

                if [ "\${NOT_GZIP}" -eq "0" ]; then
                    echo "${prefix} FASTQs are empty. Please check the input FASTQs.
                        Further analysis is discontinued." | \\
                    sed 's/^\\s*//' > ${prefix}-empty-error.txt
                    ERROR=1
                else
                    echo "${prefix} FASTQs failed Gzip tests. Please check the input FASTQs.
                        Further analysis is discontinued." | \\
                    sed 's/^\\s*//' > ${prefix}-gzip-error.txt
                    ERROR=1
                fi
            fi
            rm r1.json r2.json
        else
            # Single-end
            IS_PAIRED="false"
            gzip -cd fastqs/${prefix}.fastq.gz | fastq-scan > r1.json

            if [[ -s r1.json ]]; then
                if ! check-fastqs.py --fq1 r1.json \${OPTS}; then
                    ERROR=1
                fi
            else
                if ! gzip -t fastqs/${prefix}.fastq.gz; then
                    echo "${prefix} FASTQs failed Gzip tests. Please check the input FASTQs.
                        Further analysis is discontinued." | \\
                    sed 's/^\\s*//' > ${prefix}-gzip-error.txt
                    ERROR=1
                elif ! gzip -t fastqs/${prefix}_R2.fastq.gz; then
                    echo "${prefix} FASTQs are empty. Please check the input FASTQs.
                        Further analysis is discontinued." | \\
                    sed 's/^\\s*//' > ${prefix}-empty-error.txt
                    ERROR=1
                fi
            fi
            rm r1.json
        fi

        # Short polish should not be considered paired-end
        if [ "${runtype}" == "short_polish" ]; then
            IS_PAIRED="false"
        fi

        # Failed validations so, let's keep them from continuing
        if [ "\${ERROR}" -eq "1" ]; then
            mv fastqs/ failed-tests-fastqs/
        fi
    fi

    # Determine paired status
    IS_PAIRED="unknown"
    if [ -f  "fastqs/${prefix}_R2.fastq.gz" ]; then
        # Paired-end
        IS_PAIRED="true"
    else
        # Single-end
        IS_PAIRED="false"
    fi

    # Short polish should not be considered paired-end
    if [ "${runtype}" == "short_polish" ]; then
        IS_PAIRED="false"
    fi

    # Dump meta values to a TSV
    echo "sample<TAB>runtype<TAB>original_runtype<TAB>is_paired<TAB>is_compressed<TAB>species<TAB>genome_size" | sed 's/<TAB>/\t/g' > ${prefix}-meta.tsv
    echo "${meta.name}<TAB>${meta.runtype}<TAB>${meta.original_runtype}<TAB>\$IS_PAIRED<TAB>${meta.is_compressed}<TAB>${meta.species}<TAB>${meta.genome_size}" | sed 's/<TAB>/\t/g' >> ${prefix}-meta.tsv

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        art: \$(echo \$(art_illumina --help 2>&1) | sed 's/^.*Version //;s/ .*\$//')
        fastq-dl: \$(echo \$(fastq-dl --version 2>&1) | sed 's/fastq-dl, version //')
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        ncbi-genome-download: \$(echo \$(ncbi-genome-download --version 2>&1))
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
    END_VERSIONS
    """
}
