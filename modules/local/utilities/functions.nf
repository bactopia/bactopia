/*
  Utility functions used in Bactopia DSL2 module files
*/
import nextflow.util.SysHelper

def get_resources(profile, max_memory, max_cpus) {
        // Adjust memory/cpu requests for standard profile only
    def Map resources = [:]
    resources.MAX_MEMORY = ['standard', 'docker', 'singularity'].contains(profile) ? _get_max_memory(max_memory).GB : (max_memory).GB
    resources.MAX_MEMORY_INT = resources.MAX_MEMORY.toString().split(" ")[0]
    resources.MAX_CPUS = ['standard', 'docker', 'singularity'].contains(profile) ? _get_max_cpus(max_cpus.toInteger()) : max_cpus.toInteger()
    resources.MAX_CPUS_75 = Math.round(resources.MAX_CPUS * 0.75)
    resources.MAX_CPUS_50 = Math.round(resources.MAX_CPUS * 0.50)
    return resources
}

/* Get the maximum available memory for the given system */
def _get_max_memory(requested) {
    available = Math.floor(Double.parseDouble(SysHelper.getAvailMemory().toString().split(" ")[0])).toInteger()
    if (available < requested) {
        log.warn "Maximum memory (${requested}) was adjusted to fit your system (${available})"
        return available
    }

    return requested
}

/* Get the maximum available cpus for the given system */
def _get_max_cpus(requested) {
    available = SysHelper.getAvailCpus()
    if (available < requested) {
        log.warn "Maximum CPUs (${requested}) was adjusted to fit your system (${available})"
        return available
    }
    
    return requested
}

def save_files(Map args) {
    /* Modeled after nf-core/modules saveFiles function */
    final_output = ""
    found_ignore = false
    logs_subdir = args.containsKey('logs_subdir') ? args.logs_subdir : ""
    if (args.filename.endsWith('.version.txt') || args.filename.endsWith('.stderr.txt') || args.filename.endsWith('.stdout.txt')) {
        // Its a version file or  program specific log files
        final_output = "logs/${args.process_name}/${logs_subdir}/${args.filename}"
    } else if (args.filename.startsWith('.command')) {
        // Its a Nextflow process file, rename to "nf-<PROCESS_NAME>.*"
        ext = args.filename.replace(".command.", "")
        final_output = "logs/${args.process_name}/${logs_subdir}/nf-${args.process_name}.${ext}"
    } else {
        // Its a program output
        filename = args.filename
        if (filename.startsWith("results/")) {
            filename = args.filename.replace("results/","")
        }

        // *-error.txt should be at the base dir and 'blastdb' should go in blast folder
        final_output = null
        if (filename.endsWith("-error.txt")) {
            final_output = filename
        } else if (filename.startsWith("blastdb/")) {
            final_output = "blast/${filename}"
        } else {
            final_output = "${params.publish_dir[args.process_name]}/${filename}"
        }
        
        if (args.containsKey('ignore')) {
            args.ignore.each {
                if (filename.endsWith("${it}")) {
                    final_output = null
                }
            }
        }
    }

    return final_output ? final_output.replace("//", "/") : final_output
}
