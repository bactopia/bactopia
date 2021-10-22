/*
========================================================================================
    Nextflow Functions specific to Bactopia
========================================================================================
*/
import nextflow.util.SysHelper

def get_schemas() {
    def schemas = []
    def is_workflow = params.workflows[params.wf].containsKey('is_workflow')

    if (is_workflow == true) {
        // Named workflow based of Bactopia
        schemas << 'conf/schema/bactopia.json'
    } else {
        // Bactopia Tool
        schemas << 'conf/schema/bactopia-tools.json'
    }

    if (params.workflows[params.wf].containsKey('includes')) {
        // Wrapper around multiple workflows
        schemas += _get_include_schemas(params.workflows[params.wf]["includes"])
    } else if (params.workflows[params.wf].containsKey('modules')) {
        // Workflow or Subworkflow
        schemas += _get_module_schemas(params.workflows[params.wf]["modules"])
    } else if (params.workflows[params.wf].containsKey('path')) {
        // Module
        schemas << params.workflows[params.wf].path
    }

    schemas << 'conf/schema/generic.json'
    return schemas
}

def _get_include_schemas(includes) {
    def include_schemas = []
    includes.each { it ->
        if (params.workflows[it].containsKey('modules')) {
            include_schemas += _get_module_schemas(params.workflows[it]['modules'])
        }
    }
    return include_schemas
}

def _get_module_schemas(modules) {
    def module_schemas = []
    modules.each { it ->
        module_schemas << "${params.workflows[it].path}/params.json"
    }
    return module_schemas
}

def get_resources(profile, max_memory, max_cpus) {
    /* Adjust memory/cpu requests for standard profile only */
    def Map resources = [:]
    resources.MAX_MEMORY = ['standard', 'docker', 'singularity'].contains(profile) ? _get_max_memory(max_memory).GB : (max_memory).GB
    resources.MAX_MEMORY_INT = resources.MAX_MEMORY.toString().split(" ")[0]
    resources.MAX_CPUS = ['standard', 'docker', 'singularity'].contains(profile) ? _get_max_cpus(max_cpus.toInteger()) : max_cpus.toInteger()
    resources.MAX_CPUS_75 = Math.round(resources.MAX_CPUS * 0.75)
    resources.MAX_CPUS_50 = Math.round(resources.MAX_CPUS * 0.50)
    return resources
}

def _get_max_memory(requested) {
    /* Get the maximum available memory for the given system */
    available = Math.floor(Double.parseDouble(SysHelper.getAvailMemory().toGiga().toString().split(" ")[0])).toInteger()
    if (available < requested) {
        log.warn "Maximum memory (${requested}) was adjusted to fit your system (${available})"
        return available
    }

    return requested
}

def _get_max_cpus(requested) {
    /* Get the maximum available cpus for the given system */
    available = SysHelper.getAvailCpus()
    if (available < requested) {
        log.warn "Maximum CPUs (${requested}) was adjusted to fit your system (${available})"
        return available
    }
    
    return requested
}

def print_efficiency(cpus) {
    /* Inform user how local bactiopia run will use resources */
    if (['standard', 'docker', 'singularity'].contains(workflow.profile)) {
        // This is a local run on a single machine
        available = SysHelper.getAvailCpus()
        tasks = available / cpus
        log.info """
            Each task will use ${cpus} CPUs out of the available ${available} CPUs. At most 
            ${tasks} task(s) will be run at a time, this can affect the efficiency 
            of Bactopia. You can use the '-qs' parameter to alter the number of 
            tasks to run at a time (e.g. '-qs 2', means only 2 tasks or a maximum 
            of ${2 * cpus} CPUs will be used at once)
        """.stripIndent()
        log.info ""
    }
}

def saveFiles(Map args) {
    /* Modeled after nf-core/modules saveFiles function */
    final_output = ""
    found_ignore = false
    logs_subdir = args.containsKey('logs_subdir') ? args.logs_subdir : ""
    subworkflow = args.containsKey('subworkflow') ? args.subworkflow : ""
    if (args.filename) {
        if (args.filename.equals('versions.yml') && !System.getenv("BACTOPIA_TEST")) {
            // Do not publish versions.yml unless running from pytest workflow
            // Adapted from nf-core/modules
            return null
        } else if (args.filename.startsWith('.command')) {
            // Its a Nextflow process file, rename to "nf-<PROCESS_NAME>.*"
            ext = args.filename.replace(".command.", "")
            final_output = "logs/${subworkflow}/${args.process_name}/${logs_subdir}/nf-${args.process_name}.${ext}"
        } else if (args.filename.endsWith('.stderr.txt') || args.filename.endsWith('.stdout.txt') || args.filename.endsWith('.log')  || args.filename.endsWith('.err') || args.filename.equals('versions.yml')) {
            // Its a version file or  program specific log files
            final_output = "logs/${subworkflow}/${args.process_name}/${logs_subdir}/${args.filename}"
        } else {
            // Its a program output
            filename = args.filename
            if (filename.startsWith("results/")) {
                filename = args.filename.replace("results/","")
            }

            // *-error.txt should be at the base dir and 'blastdb' should go in blast folder
            final_output = null
            if (filename.endsWith("-error.txt") || args.publish_to_base == true) {
                final_output = filename
            } else if (filename.startsWith("blastdb/")) {
                final_output = "blast/${filename}"
            } else if (filename.startsWith("total_contigs_")) {
                final_output = null
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
}

/*
========================================================================================
    Nextflow Functions specific to nf-core/modules
    Taken from nf-vore/modules functions.nf
========================================================================================
*/

def getSoftwareName(task_process, full_name) {
    /* Extract name of software tool from process name using $task.process */
    if (full_name == true) {
        return task_process.tokenize(':')[-1].toLowerCase()
    } else {
        return task_process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()
    }
}

def getProcessName(task_process) {
    /* Extract name of module from process name using $task.process */
    return task_process.tokenize(':')[-1].toLowerCase()
}

def initOptions(Map args) {
    /* Function to initialise default values and to generate a Groovy Map of available options for nf-core modules */
    def Map options = [:]
    options.args            = args.args ?: ''
    options.args2           = args.args2 ?: ''
    options.args3           = args.args3 ?: ''
    options.publish_by_meta = args.publish_by_meta ?: []
    options.publish_dir     = args.publish_dir ?: ''
    options.publish_files   = args.publish_files
    options.suffix          = args.suffix ?: ''
    options.subworkflow     = args.subworkflow ?: ''
    options.publish_to_base = args.publish_to_base ?: false
    options.full_software_name = args.full_software_name ?: false
    return options
}
