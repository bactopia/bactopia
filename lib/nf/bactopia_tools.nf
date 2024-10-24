/*
========================================================================================
    Nextflow Functions specific to Bactopia Tools
========================================================================================
/*

/*
========================================================================================
    Public: collect_samples(bactopia_dir, extension, exclude_list, include_list)

    bactopia_dir -> the path to bactopia outputs
    extension -> the file extension to collect (e.g. 'fastq', 'fna', gff', ...)
    exclude_list -> a text file with a subset of samples to exclude from analysis
    include_list -> a text file with a subset of samples to only include from analysis
========================================================================================
*/
def collect_samples(bactopia_dir, extension, include_list, exclude_list) {
    include_all = true
    inclusions = []
    exclusions = []
    IGNORE_LIST = ['.nextflow', 'bactopia-info', 'bactopia-tools', 'work', 'bactopia-runs']
    if (include_list) {
        file(include_list).eachLine { line ->
            inclusions << line.trim().split('\t')[0]
        }
        include_all = false
        log.info "Including ${inclusions.size} samples for analysis"
    }
    else if (exclude_list) {
        file(exclude_list).eachLine { line ->
            exclusions << line.trim().split('\t')[0]
        }
        log.info "Excluding ${exclusions.size} samples from the analysis"
    }
    
    sample_list = []
    missing = []
    file("${bactopia_dir}/").eachFile { item ->
        if( item.isDirectory() ) {
            sample = item.getName()
            if (!IGNORE_LIST.contains(sample)) {
                if (inclusions.contains(sample) || include_all) {
                    if (!exclusions.contains(sample)) {
                        if (_is_sample_dir(sample, bactopia_dir)) {
                            sample = _collect_inputs(sample, bactopia_dir, extension)
                            if (sample instanceof String) {
                                missing << sample
                            } else {
                                sample_list << sample
                            }
                        } else {
                            log.info "${sample} does not appear to be a Bactopia sample, skipping..."
                        }
                    }
                }
            }
        }
    }

    log.info "Found ${sample_list.size} samples to process"
    if (missing.size() > 0) {
        log.info "${missing.size} samples were excluded due to missing files. They are:"
        for (sample in missing) {
            log.info "    ${sample}"
        }
    }
    log.info "\nIf this looks wrong, now's your chance to back out (CTRL+C 3 times)."
    log.info "Sleeping for 5 seconds..."
    log.info "--------------------------------------------------------------------"
    sleep(5000)
    return sample_list
}

/*
========================================================================================
    Public: collect_samples(file_path)

    key -> parameter that represents the filr
    local_file -> path to the files to collect
    local_pattern -> Pattern to match incase `local_file` is a directory
========================================================================================
*/
def collect_local_files(local_file, local_pattern) {
    local_list = []
    if (local_file) {
        if (file(local_file).exists()) {
            if (file(local_file).isDirectory()) {
                file("${local_file}/${local_pattern}" ).each{ 
                    local_list << tuple([[id: it.getSimpleName()], it])                    
                }
            } else {
                local_list << tuple([id: file(local_file).getSimpleName()], file(local_file))
            }
        }
    }
    return local_list
}

/*
========================================================================================
    Private: _is_sample_dir(sample, dir)

    This function does a quick check to see if the *-genome-size.txt file exists, which
    is a required output and all samples processed by Bactopia should have it.

    sample -> The name of the sample
    dir -> The base directory of sample results
========================================================================================
*/
def _is_sample_dir(sample, dir) {
    return file("${dir}/${sample}").exists()
}

/*
========================================================================================
    Private: _collect_inputs(sample, dir, extension)

    This function collects outputs based on file extension. As more outputs are needed 
    the PATHS map will need to be updated.

    sample -> The name of the sample
    dir -> The base directory of sample results
    extension -> The file extension to match.
========================================================================================
*/
def _collect_inputs(sample, dir, extension) {
    PATHS = [:]
    PATHS.blastdb = "annotator"
    PATHS.fastq = "qc"
    PATHS.fna = "assembler"
    PATHS.faa = "annotator"
    PATHS.gbk = "annotator"
    PATHS.gff = "annotator"
    PATHS.meta = "gather"

    base_dir = "${dir}/${sample}/main/"
    se = "${base_dir}/${PATHS['fastq']}/${sample}.fastq.gz"
    ont = "${base_dir}/${PATHS['fastq']}/${sample}-final_NanoPlot-report.html"
    pe1 = "${base_dir}/${PATHS['fastq']}/${sample}_R1.fastq.gz"
    pe2 = "${base_dir}/${PATHS['fastq']}/${sample}_R2.fastq.gz"
    fna = "${base_dir}/${PATHS['fna']}/${sample}.fna"
    meta = "${base_dir}/${PATHS['meta']}/${sample}-meta.tsv"

    ont = false
    if (file("${base_dir}/${PATHS['fastq']}/summary/${sample}-final_NanoPlot-report.html").exists()) {
        // the se read is ONT data
        ont = true
    }

    if (extension == "illumina_fastq") {
        // Prioritize PE reads first
        if (file(pe1).exists() && file(pe2).exists()) {
            return tuple([id:sample, single_end:false, runtype:'illumina'], [file(pe1), file(pe2)])
        } else if (file(se).exists() && !ont) {
            return tuple([id:sample, single_end:true, runtype:'illumina'], [file(se)])
        }
    } else if (extension == 'fastq') {
        if (file(se).exists()) {
            if (ont) {
                return tuple([id:sample, single_end:true, runtype:'ont'], [file(se)])
            } else {
                return tuple([id:sample, single_end:true, runtype:'illumina'], [file(se)])
            }
        } else if (file(pe1).exists() && file(pe2).exists()) {
            return tuple([id:sample, single_end:false, runtype:'illumina'], [file(pe1), file(pe2)])
        }
    } else if (extension == 'fna_fastq') {
        if (file(se).exists()) {
            runtype = "illumina"
            if (ont) {
                runtype = "ont"
            }

            if (file("${fna}.gz").exists()) {
                return tuple([id:sample, single_end:true, is_compressed:true, runtype:runtype], [file("${fna}.gz")], [file(se)])
            } else if (file(fna).exists()) {
                return tuple([id:sample, single_end:true, is_compressed:false, runtype:runtype], [file("${fna}")], [file(se)])
            }
        } else if (file(pe1).exists() && file(pe2).exists()) {
            if (file("${fna}.gz").exists()) {
                return tuple([id:sample, single_end:false, is_compressed:true, runtype:'illumina'], [file("${fna}.gz")], [file(pe1), file(pe2)])
            } else if (file(fna).exists()) {
                return tuple([id:sample, single_end:false, is_compressed:false, runtype:'illumina'], [file("${fna}")], [file(pe1), file(pe2)])
            }
        }
    } else if (extension == 'fna_faa_gff') {
        // Default to Bakta faa
        fna = "${base_dir}/${PATHS['faa']}/bakta/${sample}.fna"
        faa = "${base_dir}/${PATHS['faa']}/bakta/${sample}.faa"
        gff = "${base_dir}/${PATHS['faa']}/bakta/${sample}.gff3"
        if (!file("${faa}").exists() && !file("${faa}.gz").exists()) {
            // Fall back on Prokka
            fna = "${base_dir}/${PATHS['faa']}/prokka/${sample}.fna"
            faa = "${base_dir}/${PATHS['faa']}/prokka/${sample}.faa"
            gff = "${base_dir}/${PATHS['faa']}/prokka/${sample}.gff"
        }

        if (file("${fna}.gz").exists() && file("${faa}.gz").exists() && file("${gff}.gz").exists()) {
            return tuple([id:sample, is_compressed:true], [file("${fna}.gz")], [file("${faa}.gz")], [file("${gff}.gz")])
        } else if (file(fna).exists() && file(faa).exists() && file(gff).exists()) {
            return tuple([id:sample, is_compressed:false], [file("${fna}")], [file("${faa}")], [file("${gff}")])
        }
    } else if (extension == 'fna_faa') {
        // Default to Bakta faa
        faa = "${base_dir}/${PATHS['faa']}/bakta/${sample}.faa"
        if (!file("${faa}").exists() && !file("${faa}.gz").exists()) {
            // Fall back on Prokka
            faa = "${base_dir}/${PATHS['faa']}/prokka/${sample}.faa"
        }

        if (file("${fna}.gz").exists() && file("${faa}.gz").exists()) {
            return tuple([id:sample, is_compressed:true], [file("${fna}.gz")], [file("${faa}.gz")])
        } else if (file(fna).exists() && file(faa).exists()) {
            return tuple([id:sample, is_compressed:false], [file("${fna}")], [file("${faa}")])
        }
    } else if (extension == 'fna_meta') {
        // include the meta file
        if (file("${fna}.gz").exists() && file(meta).exists()) {
            return tuple([id:sample, is_compressed:true], [file("${fna}.gz")], [file(meta)])
        } else if (file(fna).exists() && file(meta).exists()) {
            return tuple([id:sample, is_compressed:false], [file("${fna}")], [file(meta)])
        }
    } else if (extension == 'blastdb') {
        // Default to Bakta blastdb
        input = "${base_dir}/${PATHS[extension]}/bakta/${sample}-${extension}.tar.gz"
        if (!file("${input}").exists()) {
            // Fall back on Prokka
            input = "${base_dir}/${PATHS[extension]}/prokka/${sample}-${extension}.tar.gz"
        }

        if (file("${input}").exists()) {
            return tuple([id:sample], [file("${input}")])
        }
    } else {
        input = "${base_dir}/${PATHS[extension]}/${sample}.${extension}"
        if (extension == "gbk") {
            // Default to Bakta (gbff)
            input = "${base_dir}/${PATHS[extension]}/bakta/${sample}.gbff"
            if (!file("${input}").exists() && !file("${input}.gz").exists()) {
                // Fall back on Prokka (gbk)
                input = "${base_dir}/${PATHS[extension]}/prokka/${sample}.${extension}"
            }
        } else if (extension == "gff") {
            // Default to Bakta (gff3)
            input = "${base_dir}/${PATHS[extension]}/bakta/${sample}.gff3"
            if (!file("${input}").exists() && !file("${input}.gz").exists()) {
                // Fall back on Prokka (gff)
                input = "${base_dir}/${PATHS[extension]}/prokka/${sample}.${extension}"
            }
        } else if (extension == "faa") {
            // Default to Bakta faa
            input = "${base_dir}/${PATHS[extension]}/bakta/${sample}.${extension}"
            if (!file("${input}").exists() && !file("${input}.gz").exists()) {
                // Fall back on Prokka
                input = "${base_dir}/${PATHS[extension]}/prokka/${sample}.${extension}"
            }
        }

        if (file("${input}.gz").exists()) {
            return tuple([id:sample, is_compressed:true], [file("${input}.gz")])
        } else if (file(input).exists()) {
            return tuple([id:sample, is_compressed:false], [file("${input}")])
        }
    }

    return sample
}
