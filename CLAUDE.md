# Bactopia Pipeline Reference for Claude

This document contains essential information about the Bactopia pipeline architecture, patterns, and conventions to help Claude understand the codebase structure and make informed decisions.

## Quick Reference

### Common Patterns
- **Entry workflow template**: See [Entry Workflow Structure](#entry-workflow-structure)
- **Module template**: See [Module Template](#module-template)
- **Subworkflow template**: See [Subworkflow Interface](#subworkflow-interface)

### Frequently Asked Questions
- **Why `file()` vs `files()`?**: See [File Output Types](#file-output-types)
- **What are EMPTY_* files?**: See [Path? Optional Parameters](#path-optional-parameters)
- **How to add a new tool?**: See [Adding a New Tool](#adding-a-new-tool)

### Common Gotchas
- Don't "fix" Path? workarounds - they're intentional
- Always emit 4 channels from subworkflows (results, logs, nf_logs, versions)
- Module inputs use `Tuple<Map, Set<Path>>`, not `Tuple<Map, Path>`

### Decision Trees

#### I Need To
- **Add a new tool**: Follow [Adding a New Tool](#adding-a-new-tool)
- **Debug a type error**: Check [Common Misconceptions](#common-misconceptions) first
- **Modify input handling**: See [Module Input/Output Patterns](#module-inputoutput-patterns)
- **Update configuration**: See [Configuration System](#configuration-system)
- **Write tests**: See [Testing Framework](#testing-framework)

## Troubleshooting

### Common Error Messages
- **"Cannot cast X to Y"**: Usually means using `file()` when you need `files()`
- **"Missing required parameter"**: Check schema.json and params.config
- **"Channel not found"**: Verify all 4 channels are emitted from subworkflows

### Debugging Tips
1. Check the Important Notes section first
2. Look for TODO comments - they indicate known limitations
3. Verify meta map fields are correctly set
4. Ensure consistent typing across connected components

## Table of Contents
1. [Repository Structure](#repository-structure)
2. [Architecture Overview](#architecture-overview)
3. [Static Typing Conventions](#static-typing-conventions)
4. [Design Patterns](#design-patterns)
5. [Component Tag Taxonomy](#component-tag-taxonomy)
6. [Configuration System](#configuration-system)
7. [Testing Framework](#testing-framework)
8. [Development Workflow](#development-workflow)
9. [Known Limitations](#known-limitations)
10. [Common Examples](#common-examples)
11. [Important Notes](#important-notes)
12. [Version Requirements](#version-requirements)
13. [Groovydoc Documentation Standards](#groovydoc-documentation-standards)
14. [Glossary](#glossary)

## Repository Structure

```
bactopia/
├── [.github/](#github)           # Github Actions workflows and issue templates
├── [.vscode/](#vscode)           # Visual Studio Code settings and configurations
├── [bin/](#bin)                   # Helper scripts and utilities
├── [conf/](#configuration-system) # Configuration files
├── [data/](#data)                 # Static data and resources
│   └── [empty/](#empty-files)     # Empty placeholder files
├── [modules/](#tier-3-modules-modules) # Individual process implementations
│   ├── [bactopia/](#core-modulesmodulesbactopia) # Core pipeline processes
│   └── {tool}/                    # External tool implementations
├── [subworkflows/](#tier-2-subworkflows-subworkflows) # Reusable workflow components
│   ├── [bactopia/](#core-subworkflowssubworkflowsbactopia) # Essential pipeline components
│   ├── [utils/](#utility-subworkflowssubworkflowsutils) # Helper workflows
│   └── {tool}/                    # Tool-specific processing logic
├── [tests/](#testing-framework)   # Test suite
├── [workflows/](#tier-1-workflows-workflows) # Entry point workflows
│   ├── [bactopia-tools/](#bactopia-tools-workflowsbactopia-tools) # Workflows requiring Bactopia outputs
│   └── {named_workflow}/         # Standalone workflows
├── main.nf                       # [Main Bactopia workflow](#main-workflow-mainnf)
├── nextflow.config               # [Global configuration](#nextflow-config-root)
├── nextflow_schema.json          # [Parameter validation schema](#schema-validation)
└── README.md                     # Project documentation
```

### Directory Details

#### `/.github/`
GitHub Actions and community management:
- `workflows/` - CI/CD pipeline automation
    - `all-bactopia-tests.yml` - Run all Bactopia tests using Docker, Singularity and Conda
    - `conda-build.manual.yml` - Create a development Conda package made available on rpetit3 channel
- `ISSUE_TEMPLATE/` - Standardized issue templates

#### `/.vscode/`
Visual Studio Code workspace settings

#### `/bin/`
Collection of utility scripts used by the pipeline:
- `bactopia/` - Main Bactopia command-line interface
- `bactopia-prokka` - Prokka wrapper script
- `check-fastqs.py` - FASTQ file validation
- `kraken-bracken-summary.py` - Kraken/Bracken result aggregation
- `mlst-blast.py` - MLST BLAST analysis
- `teton-prepare.py` - Teton workflow preparation
- Various helper scripts for data processing

#### `/conf/`
Pipeline configuration files:
- `base.config` - Base configuration settings
- `params.config` - Parameter definitions
- `profiles.config` - Execution profile configurations
- `test.config` - Test-specific configuration
- `workflows.yaml` - Workflow metadata and definitions
- `params/` - Tool-specific parameter files
- `schema/` - JSON schema definitions

#### `/data/`
Static data and resources:
- `empty/` - Empty placeholder files for Path? workarounds ([see explanation](#empty_-files))
- `conda/` - Conda environment specifications for the development environment
- `proteins.faa` - Reference protein database used with Prokka
- `*.txt` - Accession lists and reference data
- `*.png,*.jpg` - Documentation images and logos

#### `/modules/`
Nextflow Process implementations:
- The singular building blocks of the pipeline
- Each module contains the actual tool execution logic
- The sub-folders are named by the tool they implement, and may contain multiple sub-folders for different uses of the tool
    - Example Single Subfolder: `bracken/`
    - Example Multiple Subfolders: `blast/` -> `blastn/`, `blastp/`, etc.
- Each module folder contains the following files:
    - `main.nf` - Process definitions
    - `meta.yaml` - Module metadata
    - `params.config` - Module-specific parameters
    - `process.config` - Process resource configurations
    - `schema.json` - Module-specific parameter schema in JSON format (used in creating workflow-specific schemas)

#### `/subworkflows/`
Reusable Nextflow workflow components:
- Subworkflows are collections of other subworkflows and/or modules
- Each subworkflow contains a set of related processes that can be reused across multiple workflows
- The sub-folders are named by the tool they implement, and may contain multiple sub-folders for different uses of the tool
    - Example Single Subfolder: `bracken/`
    - Example Multiple Subfolders: `bactopia/` -> `assembler/`, `qc/`, etc.
- `utils/` - Utility subworkflows for gathering inputs in Bactopia, Named workflows, and Bactopia Tools
- Each subworkflow folder contains the following files:
    - `main.nf` - Subworkflow definitions
    - `meta.yaml` - Subworkflow metadata

#### `/tests/`
Bactopia Pipeline test suite:
- `main.nf.test` - Main workflow tests
- `main.nf.test.snap` - Test snapshots
- `nf-test.config` - Test configuration
- Uses nf-test framework for pipeline testing

#### `/workflows/`
Entry point workflows built using subworkflows and modules:
- Every workflow has its own directory with the following files:
    - `tests/` - Workflow-specific tests previously described in `/tests/`
    - `main.nf` - Entry workflow script
    - `nextflow_schema.json` - Workflow-specific parameter schema in JSON format
    - `nextflow.config` - Workflow-specific configuration
- Workflow types:
    - `bactopia-tools/` - Individual Bactopia Tool workflows
        - Examples: abricate, prokka, blastn, etc.
        - Requires outputs from Bactopia main workflow, or Named Workflows
    - Named workflows are standalone workflows that do not require Bactopia outputs
        - Examples: clean-yer-reads, staphopia, teton, etc.

## Architecture Overview

### Three-Tier Structure
The Bactopia pipeline follows a clear three-tier architecture:

1. **Workflows** (`/main.nf` and `/workflows/`)
    - Entry points for pipeline execution
    - Main Bactopia workflow: `main.nf`
    - Tool-specific workflows: `/workflows/bactopia-tools/*/main.nf`
    - Named workflows: `/workflows/{named_workflow}/main.nf` (excluding `bactopia-tools`)
    - Handle user-facing parameters and orchestration

2. **Subworkflows** (`/subworkflows/`)
    - Reusable workflow components
    - Core functionality: `/subworkflows/bactopia/`
    - Tool-specific: `/subworkflows/{toolname}/`
    - Utilities: `/subworkflows/utils/`

3. **Modules** (`/modules/`)
    - Individual process implementations
    - Tool execution logic
    - Basic building blocks of the pipeline

### Tier 1: Workflows (`/workflows/`)

Entry point workflows that users directly interact with. Each workflow directory contains:

#### Required Files
- **`main.nf`** - Main workflow script
    - Must start with `#!/usr/bin/env nextflow`
    - Until Nextflow v26.04, must include `nextflow.preview.types = true`
    - Must include a `params` block defining workflow parameters used within the entry workflow
        - `params` should only be available in the entry workflow

        ```nextflow
        params {
            rundir   : String

            // Tool-specific parameters
            // QC
            adapters              : Path?
            phix                  : Path?

            // Annotation
            use_bakta             : Boolean

            // Merlin
            ask_merlin            : Boolean
        }
        ```

    - Includes an initialization workflow
        - `BACTOPIA_INIT` for Bactopia and Named Workflows

        ```nextflow
        include { BACTOPIA_INIT } from '../subworkflows/bactopia/main.nf'
        ```

        - `BACTOPIATOOL_INIT` for Bactopia Tools

        ```nextflow
        include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main.nf'
        ```

    - Includes tool-specific subworkflows

        ```nextflow
        include { AMRFINDERPLUS   } from './subworkflows/amrfinderplus/main'
        include { ASSEMBLER       } from './subworkflows/bactopia/assembler/main'
        include { DATASETS        } from './modules/bactopia/datasets/main'
        ```

    - Workflow block
        - Defines main execution logic
        - Follows 4-channel output pattern (results, logs, nf_logs, versions)
        - Implements branching by scope (run/sample) for outputs
        - Defines a publish block for run/sample outputs
            - Uses consistent naming for output channels:
                - `run_results`, `sample_results`
                - `run_logs`, `sample_logs`
                - `run_nf_logs`, `sample_nf_logs`
                - `run_versions`, `sample_versions`

        ```nextflow
        workflow {
            main:
            // Initialize output channels
            ch_results = channel.empty() as Channel<Tuple<Map, Path>>
            ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
            ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
            ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

            // Execute subworkflows
            BACTOPIA_INIT()

            // Core steps
            DATASETS()

            // Gather samples in one place
            GATHER(BACTOPIA_INIT.out.samples)
            ch_results = ch_results.mix(GATHER.out.results)
            ch_logs = ch_logs.mix(GATHER.out.logs)
            ch_nf_logs = ch_nf_logs.mix(GATHER.out.nf_logs)
            ch_versions = ch_versions.mix(GATHER.out.versions)

            // ... CODE OMITTED FOR BREVITY ...

            // Branch the based on scope (sample or run)
            ch_final_results = ch_results.branch{ meta, _file ->
                run: meta.scope == 'run'
                sample: meta.scope == 'sample'
            }

            ch_final_logs = ch_logs.branch{ meta, _file ->
                run: meta.scope == 'run'
                sample: meta.scope == 'sample'
            }

            ch_final_nf_logs = ch_nf_logs.branch{ meta, _file ->
                run: meta.scope == 'run'
                sample: meta.scope == 'sample'
            }

            ch_final_versions = ch_versions.branch{ meta, _file ->
                run: meta.scope == 'run'
                sample: meta.scope == 'sample'
            }

            publish:
            run_results = ch_final_results.run
            run_logs = ch_final_logs.run
            run_nf_logs = ch_final_nf_logs.run
            run_versions = ch_final_versions.run
            sample_results = ch_final_results.sample
            sample_logs = ch_final_logs.sample
            sample_nf_logs = ch_final_nf_logs.sample
            sample_versions = ch_final_versions.sample
        }
        ```

    - Output block
        - Defines run-level and sample-level outputs for all four channels
        - Uses `meta.output_dir` and `meta.logs_dir` from modules for output paths
        - Follows consistent pattern in all workflows

        ```nextflow
        output {
            // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
            run_results: Channel<Tuple<Map, Path>> {
                path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
            }
            run_logs: Channel<Tuple<Map, Path>> {
                path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
            }
            run_nf_logs: Channel<Tuple<Map, Path>> {
                path { meta, file ->
                    file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
                }
            }
            run_versions: Channel<Tuple<Map, Path>> {
                path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
            }

            // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
            sample_results: Channel<Tuple<Map, Path>> {
                path { meta, _file -> "${meta.output_dir}/" }
            }
            sample_logs: Channel<Tuple<Map, Path>> {
                path { meta, _file -> "${meta.logs_dir}/" }
            }
            sample_nf_logs: Channel<Tuple<Map, Path>> {
                path { meta, file ->
                    file >> "${meta.logs_dir}/nf${file.name}"
                }
            }
            sample_versions: Channel<Tuple<Map, Path>> {
                path { meta, _file -> "${meta.logs_dir}/" }
            }
        }
        ```

- **`nextflow_schema.json`** - Parameter validation schema
    - JSON schema defining all workflow parameters
    - Generated by merging module schemas (`modules/{tool_name}/schema.json`) included in the workflow as defined in `conf/workflows.yaml`
        - Merging done with `bactopia-merge-schemas` from [bactopia-py](https://github.com/bactopia/bactopia-py)
    - Used to validate user inputs at runtime

- **`nextflow.config`** - Workflow-specific configuration
    - Workflow-specific settings, all workflows will follow similar patterns
    - All `nextflow.config` should include:
        - `manifest` block
        - `params.workflow` block
        - includes `params.config` and `process.config` from included modules
        - include `conf/base.config`, `conf/params.conf`, `conf/profiles.config`
        - include params config based on workflow
            - Main Bactopia: `conf/params/bactopia.config`
            - Bactopia Tools: `conf/params/bactopia-tools.config`
            - Named Workflows: `conf/params/{named_workflow}.config`

- **`tests/`** directory (optional)
    - Workflow-specific tests
    - Uses nf-test framework
    - Test config (`conf/test.config`) includes global test settings
    - Contains:
        - `.nftignore` - Files to ignore during testing
        - `main.nf.test` - Workflow test definitions
        - `main.nf.test.snap` - Snapshot of expected outputs
        - `nf-test.config` - test-specific configuration
    - Test data is hosted at [bactopia-tests](https://github.com/bactopia/bactopia-tests)
        - Environment variable `BACTOPIA_TESTS` must be set and point to the local path of the cloned repository

#### Workflow Types

1. **Main Workflow** (`/main.nf`)
    - Primary Bactopia pipeline entry point
    - Orchestrates the "Bactopia" analysis
    - Includes `BACTOPIA_INIT` subworkflow

2. **Bactopia Tools** (`/workflows/bactopia-tools/`)
    - Additional workflows (e.g., abricate, prokka, blastn) which utilize Bactopia outputs
    - Require outputs from main Bactopia workflow or Named Workflows
    - Example structure
        - Simple: `/workflows/bactopia-tools/abricate/main.nf`
        - Complex: `/workflows/bactopia-tools/pangenome/main.nf` (may include multiple subworkflows/modules)
    - All follow the same patterns as Bactopia workflow
        - Includes `BACTOPIATOOL_INIT` subworkflow

3. **Named Workflows**
    - Standalone workflows that don't require Bactopia outputs
    - Restructure steps from main Bactopia workflow
    - Examples:
        - clean-yer-reads - only runs the `GATHER` and `QC` steps from Bactopia
        - staphopia - Includes all Bactopia steps, in addition to `STAPHTYPER`
        - teton - Includes `GATHER`, `SCRUBBER` and `BRACKEN` steps
    - Have the same input requirements as Bactopia
    - Follow the same patterns as Bactopia workflow
        - Includes `BACTOPIA_INIT` subworkflow

### Tier 2: Subworkflows (`/subworkflows/`)

Reusable workflow components that combine multiple modules or other subworkflows.

#### Required Files
- **`main.nf`** - Subworkflow definition
    - No `#!/usr/bin/env nextflow` shebang (not an entry point)
    - Must include `nextflow.preview.types = true` in Nextflow versions prior to v26.04
    - Includes necessary modules and/or subworkflows
    - Must always include `flattenPaths` and `gather` from `plugin/nf-bactopia` for output management
        - `flattenPaths` - Converts all output channels into `Channel<Tuple<Map, Path>>`
        - `gather` - Merges multiple channels into a single channel for aggregate outputs (e.g., summaries)
    - Defines workflow with sections:
        - `take` - Statically typed input channels
        - `main` - Core processing logic
        - `emit` - Output channels
            - output channels must include:
                - individual outputs from each module/subworkflow
                - generic aggregate outputs using `flattenPaths` and `gather`
                    - `results`, `logs`, `nf_logs`, `versions`
    - Simple example - `/subworkflows/abricate/main.nf`

    ```nextflow
    //
    // abricate - Mass screening of contigs for antimicrobial and virulence genes
    //
    nextflow.preview.types = true

    include { ABRICATE_RUN     } from '../../modules/abricate/run/main'
    include { ABRICATE_SUMMARY } from '../../modules/abricate/summary/main'
    include { flattenPaths     } from 'plugin/nf-bactopia'
    include { gather           } from 'plugin/nf-bactopia'

    workflow ABRICATE {
        take:
        fasta: Channel<Tuple<Map, Set<Path>>>

        main:
        ABRICATE_RUN(fasta)
        ABRICATE_SUMMARY(gather(ABRICATE_RUN.out.report, 'abricate'))

        emit:
        // Individual outputs
        tsv: Channel<Tuple<Map, Path>> = ABRICATE_RUN.out.report
        merged_tsv: Channel<Tuple<Map, Path>> = ABRICATE_SUMMARY.out.report

        // Generic aggregate outputs
        results: Channel<Tuple<Map, Path>> = flattenPaths([
            ABRICATE_RUN.out.report,
            ABRICATE_SUMMARY.out.report
        ])
        logs: Channel<Tuple<Map, Path>> = flattenPaths([
            ABRICATE_RUN.out.logs,
            ABRICATE_SUMMARY.out.logs
        ])
        nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
            ABRICATE_RUN.out.nf_logs,
            ABRICATE_SUMMARY.out.nf_logs
        ])
        versions: Channel<Tuple<Map, Path>> = flattenPaths([
            ABRICATE_RUN.out.versions,
            ABRICATE_SUMMARY.out.versions
        ])
    }
    ```

    - Complex example - `/subworkflows/pangenome/main.nf`

    ```nextflow
    //
    // pangenome - Pangenome analysis with optional core-genome phylogeny
    //
    nextflow.preview.types = true

    include { PIRATE       } from '../pirate/main'
    include { ROARY        } from '../roary/main'
    include { PANAROO      } from '../panaroo/main'
    include { SNPDISTS     } from '../snpdists/main'
    include { flattenPaths } from 'plugin/nf-bactopia'
    include { gather       } from 'plugin/nf-bactopia'

    workflow PANGENOME {
        take:
        gff        : Channel<Tuple<Map, Set<Path>>>
        use_pirate : Boolean
        use_roary  : Boolean

        main:

        // Initialize channels
        ch_aln = channel.empty() as Channel<Tuple<Map, Path>>
        ch_csv = channel.empty() as Channel<Tuple<Map, Path>>
        ch_results = channel.empty() as Channel<Tuple<Map, Path>>
        ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
        ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
        ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

        // Execute subworkflows
        // Choose pangenome tool based on params
        if (use_pirate) {
            PIRATE(gff)
            ch_aln = PIRATE.out.aln
            ch_csv = PIRATE.out.csv
            ch_results = PIRATE.out.results
            ch_logs = PIRATE.out.logs
            ch_nf_logs = PIRATE.out.nf_logs
            ch_versions = PIRATE.out.versions
        } else if (use_roary) {
            ROARY(gff)
            ch_aln = ROARY.out.aln
            ch_csv = ROARY.out.csv
            ch_results = ROARY.out.results
            ch_logs = ROARY.out.logs
            ch_nf_logs = ROARY.out.nf_logs
            ch_versions = ROARY.out.versions
        } else {
            PANAROO(gff)
            ch_aln = PANAROO.out.filtered_aln
            ch_csv = PANAROO.out.csv
            ch_results = PANAROO.out.results
            ch_logs = PANAROO.out.logs
            ch_nf_logs = PANAROO.out.nf_logs
            ch_versions = PANAROO.out.versions
        }

        // Per-sample SNP distances
        ch_unmasked_aln = ch_aln.map({ _meta, aln -> 
            tuple([name: "core-genome.distance", process_name: "snpdists"], aln)
        })
        SNPDISTS(ch_unmasked_aln)
        ch_results = ch_results.mix(SNPDISTS.out.results)
        ch_logs = ch_logs.mix(SNPDISTS.out.logs)
        ch_nf_logs = ch_nf_logs.mix(SNPDISTS.out.nf_logs)
        ch_versions = ch_versions.mix(SNPDISTS.out.versions)

        emit:
        // Individual outputs
        aln: Channel<Tuple<Map, Path>> = ch_aln
        csv: Channel<Tuple<Map, Path>> = ch_csv

        // Generic aggregate outputs
        results: Channel<Tuple<Map, Path>> = flattenPaths([ch_results])
        logs: Channel<Tuple<Map, Path>> = flattenPaths([ch_logs])
        nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([ch_nf_logs])
        versions: Channel<Tuple<Map, Path>> = flattenPaths([ch_versions])
    }
    ```

- **`meta.yaml`** - Subworkflow metadata
    - Defines subworkflow name and description
    - Specifies input/output contracts
    - Used by workflow documentation

    ```yaml
    name: SUBWORKFLOW_NAME
    description: Brief description of the subworkflow purpose
    input:
      - name: fasta
        type: Tuple<Map, Path>
        description: FASTA file to process
    output:
      - name: results
        type: Tuple<Map, Path>
        description: Analysis results
      - name: logs
        type: Tuple<Map, Set<Path>>
        description: Execution logs
    ```

#### Subworkflow Organization

1. **Core Subworkflows** (`/subworkflows/bactopia/`)
    - Essential pipeline components
    - Examples: assembler/, qc/, sketcher/
    - Handle fundamental processing steps
    - These can still be used in other workflows outside of Bactopia

2. **Tool Subworkflows** (`/subworkflows/{tool}/`)
    - Tool-specific processing logic
    - May contain multiple variants (e.g., blast/blastn/, blast/blastp/)
    - Example: `/subworkflows/abricate/main.nf`

3. **Utility Subworkflows** (`/subworkflows/utils/`)
    - House the `BACTOPIA_INIT` and `BACTOPIATOOL_INIT` subworkflows
        - `subworkflows/utils/bactopia/main.nf` - Initializes inputs for Bactopia and Named Workflows
        - `subworkflows/utils/bactopia-tools/main.nf` - Initializes inputs for Bactopia Tools workflows

### Tier 3: Modules (`/modules/`)

Individual process implementations that execute specific tools or commands.

#### Required Files
- **`main.nf`** - Process definition
    - modules can be invoked by subworkflows or directly by workflows, not by other modules
    - module testing is handled by workflow tests
    - Must start with `nextflow.preview.types = true` in Nextflow versions prior to v26.04
    - Defines single process with input/output specifications
    - the `tag` provides a label for monitoring progress
    - the `label` defines resource requirements available from `conf/base.config`
    - A new `meta` map is created in each process to define standard fields:
        - `id` - Unique identifier for the process run
        - `name` - Sample or prefix name
        - `scope` - 'sample' or 'run'
        - `output_dir` - Output directory path
        - `logs_dir` - Logs directory path
        - `process_name` - Name of the process/tool
        - Can be defined from existing variables such as `_meta` and `task.ext` values
        - The `meta` map is used to capture generic metadata for outputs
    - `file` is used for single files, `files` for multiple files (typically identified by wildcards)
        - `file` are `Tuple<Map, Path>`
        - `files` are `Tuple<Map, Set<Path>>`
        - Optional files should use `optional: true` in `file()` and `files()`
    - Modules may contain many outputs, but the standard outputs which should always be included are:
        - `logs` - Optional tool logs
        - `nf_logs` - Nextflow command logs (`.command.*` files)
        - `versions` - Software version information
    - the output block may include any number of individual result files/channels
    - all outputs should be emitted as `Tuple<Map, Path>` or `Tuple<Map, Set<Path>>` with the `meta` map
    - versions are captured in a `versions.yml` file using a heredoc for easy parsing
        - The version block should capture the process name and tool version
        - the expected output is yaml formatted to include process name and each tool name/version

        ```nextflow
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ${task.ext.toolName}: \$(echo \$(MODULE --version 2>&1) | sed 's/^.*MODULE //' )
        END_VERSIONS
        ```

        - The output format should follow:

        ```yaml
        "PROCESS_NAME":
            TOOL_NAME: VERSION_STRING
        ```

        - Simple example - `/modules/abricate/run/main.nf`

        ```yaml
        "ABRICATE_RUN":
            abricate: 1.0.1
        ```

    - Example module:

    ```nextflow
    nextflow.preview.types = true

    process MODULE {
        tag "${prefix}"
        label 'process_single'

        conda "${task.ext.condaDir}/${task.ext.toolName}"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker }"

        input:
        (_meta, input): Tuple<Map, Set<Path>>

        output:
        report   = tuple(meta, file("${prefix}.txt"))
        logs     = tuple(meta, files("*.{log,err}", optional: true))
        nf_logs  = tuple(meta, files(".command.*"))
        versions = tuple(meta, file("versions.yml"))

        script:
        prefix = task.ext.prefix ?: "${_meta.name}"

        // Create a new meta variable
        meta = [:]
        meta.id = "${prefix}-${task.process}"
        meta.name = prefix
        meta.scope = task.ext.scope
        meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
        meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
        meta.process_name = task.ext.process_name
        """
        MODULE \\
            $task.ext.args \\
            --threads $task.cpus > ${prefix}.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            MODULE: \$(echo \$(MODULE --version 2>&1) | sed 's/^.*MODULE //' )
        END_VERSIONS
        """
    }
    ```

- **`meta.yaml`** - Module metadata
    - Defines module name and description
    - Specifies input/output contracts

    ```yaml
    name: MODULE_NAME
    description: Brief description of the module purpose
    input:
      - name: fasta
        type: Tuple<Map, Path>
        description: FASTA file to process
    output:
      - name: results
        type: Tuple<Map, Path>
        description: Analysis results
      - name: logs
        type: Tuple<Map, Set<Path>>
        description: Execution logs
    ```

- **`params.config`** - Module parameters
    - Default parameter values
    - Parameter descriptions and types
    - all `params.config` files should define a `run_name`
    - Unlike workflow-level `params` blocks, module-level `params` blocks do not require type annotations
    - Module-level `params` can be used in config files, but within Nextflow scripts they can only be used in entry workflows

    ```groovy
    /*
    This file includes default parameter values for abricate.
    */

    params {
        // Abricate
        abricate_db = "ncbi"
        abricate_minid = 80
        abricate_mincov = 80
        run_name = "${params.workflow.name}-${params.abricate_db}"
    }
    ```

- **`process.config`** - Process configuration
    - Define `ext` variables that can be accessed within the module process as `task.ext`
    - `ext.args` - CLI options for the module process
        - Puts options into a single string for easier management
        - Maybe multiple `args` denoted with a number (e.f., `ext.args1`, `ext.args2`, etc.)
    - Environment specific information
        - `ext.toolName` - Conda environment name
        - `ext.docker` - Docker image name
        - `ext.image` - Singularity image URL
        - `ext.condaDir` - Conda environments base directory
    - Other useful metadata
        - `ext.wf` - Workflow name
        - `ext.scope` - 'sample' or 'run'
        - `ext.subdir` - Subdirectory for outputs
        - `ext.logs_subdir` - Subdirectory for logs
        - `ext.prefix` - Output file prefix
        - `ext.process_name` - Name of the process/tool, when it differs from module name
    - Required for all modules:
        - `args`, `toolName`, `docker`, `image`, `condaDir`, `wf`, `scope`, `process_name`

    ```groovy
    /*
    This file includes default process values for abricate.
    */

    process {
        withName: 'ABRICATE_RUN|ABRICATE_SUMMARY' {
            // Optional arguments
            ext.args = [
                "--db ${params.abricate_db}",
                "--minid ${params.abricate_minid}",
                "--mincov ${params.abricate_mincov}"
            ].join(' ').replaceAll("\\s{2,}", " ").trim()

            // Environment information
            ext.toolName = "bioconda::abricate=1.0.1".replace("=", "-").replace(":", "-").replace(" ", "-")
            ext.docker = "biocontainers/abricate:1.0.1--ha8f3691_1"
            ext.image = "https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1"
            ext.condaDir = "${params.condadir}"

            ext.wf = params.wf
            ext.scope = "sample"
        }

        withName: 'ABRICATE_RUN' {
            ext.subdir = params.abricate_db
            ext.logs_subdir = ""
            ext.process_name = "abricate"
        }

        withName: 'ABRICATE_SUMMARY' {
            ext.logs_subdir = "abricate-concat"
            ext.prefix = "abricate-${params.abricate_db}"
            ext.process_name = params.merge_folder
            ext.subdir = params.abricate_db
            ext.scope = "run"
        }
    }
    ```

- **`schema.json`** - Parameter schema
    - JSON schema specific to a module
    - Used to generate workflow schemas (`nextflow_schema.json`) by merging all included module schemas
    - Validates parameter types and values

    ```json
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/bactopia/bactopia/master/modules/abricate/run/schema.json",
        "title": "Abricate Module",
        "description": "A module for mass screening of contigs for antimicrobial and virulence genes",
        "type": "object",
        "$defs": {
            "abricate_parameters": {
                "title": "Abricate Parameters",
                "type": "object",
                "description": "",
                "default": "",
                "fa_icon": "fas fa-exclamation-circle",
                "properties": {
                    "abricate_db": {
                        "type": "string",
                        "default": "ncbi",
                        "description": "Database to use",
                        "fa_icon": "fas fa-expand-arrows-alt"
                    },
                    "abricate_minid": {
                        "type": "integer",
                        "default": 80,
                        "description": "Minimum DNA percent identity",
                        "fa_icon": "fas fa-expand-arrows-alt"
                    },
                    "abricate_mincov": {
                        "type": "integer",
                        "default": 80,
                        "description": "Minimum DNA percent coverage",
                        "fa_icon": "fas fa-expand-arrows-alt"
                    }
                }
            }
        },
        "allOf": [
            {
                "$ref": "#/$defs/abricate_parameters"
            }
        ]
    }
    ```

#### Module Input/Output Patterns

**Standard Single Input:**

```nextflow
    input:
    (_meta, assembly): Tuple<Map, Set<Path>>
```

**Multiple Files Input:**

```nextflow
    input:
    (_meta, assembly): Tuple<Map, Set<Path>>
    db: Path
```

**Optional Input with Workaround:**

```nextflow
    input:
    (_meta, fasta) : Tuple<Map, Set<Path>>
    proteins       : Path?
    prodigal_tf    : Path?

    script:
    def proteins_opt = proteins.getName() != "EMPTY_PROTEINS" ? "--proteins ${proteins.getName()}" : "" 
    def prodigal_opt = prodigal_tf.getName() != "EMPTY_PRODIGAL_TF" ? "--prodigaltf ${prodigal_tf.getName()}" : "" 
```

**Valued Parameters:**

```nextflow
    (_meta, assembly): Tuple<Map, Set<Path>>
    species: String
```

**Standard Output (3-Channel Pattern):**

```nextflow
    output:
    logs        = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs     = tuple(meta, files(".command.*"))
    versions    = tuple(meta, file("versions.yml"))
```

**Multiple Output Sets:**

```nextflow
    output:
    annotations = tuple(meta, files("${prefix}.{fna,fna.gz}"), files("${prefix}.{faa,faa.gz}"), files("${prefix}.{gff,gff.gz}"))
    gff         = tuple(meta, file("${prefix}.{gff,gff.gz}"))
    gbk         = tuple(meta, file("${prefix}.{gbk,gbk.gz}"))
    fna         = tuple(meta, file("${prefix}.{fna,fna.gz}"))
    faa         = tuple(meta, file("${prefix}.{faa,faa.gz}"))
    ffn         = tuple(meta, file("${prefix}.{ffn,ffn.gz}"))
    sqn         = tuple(meta, file("${prefix}.{sqn,sqn.gz}"))
    fsa         = tuple(meta, file("${prefix}.{fsa,fsa.gz}"))
    tbl         = tuple(meta, file("${prefix}.{tbl,tbl.gz}"))
    txt         = tuple(meta, file("${prefix}.txt"))
    tsv         = tuple(meta, file("${prefix}.tsv"))
    blastdb     = tuple(meta, files("${prefix}-blastdb.tar.gz"))
    logs        = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs     = tuple(meta, files(".command.*"))
    versions    = tuple(meta, file("versions.yml"))
```

#### Module Naming Conventions
- Process names: `UPPER_SNAKE_CASE` (e.g., `ABRICATE`, `PROKKA`)
- File names: `main.nf`
- Output channels: `logs`, `nf_logs`, `versions` (standard)
- Meta fields: `id`, `name`, `scope`, `output_dir`, `logs_dir`, `process_name`

### Data Flow Architecture
```
User Input
    ↓
Workflow (Tier 1) - Parameter handling, orchestration
    ↓
Subworkflow (Tier 2) - Combines related functionality
    ↓
Module (Tier 3) - Executes individual tools
    ↓
Results (4 channels: results, logs, nf_logs, versions)
```

### Component Relationships
- **Workflows** include **Subworkflows**
- **Subworkflows** include **Modules**
- **Modules** are never directly included by **Workflows**
- This hierarchical approach promotes code reuse and maintainability

## Static Typing Conventions

### Universal Typing
- ALL Nextflow files have `nextflow.preview.types = true` at the top
- Static typing is fully implemented across all 251 files
- No files are missing type annotations

### File Output Types

#### Single Files: `file()`
- Returns: `Path`
- Used in: `Tuple<Map, Path>`
- Example: `versions = tuple(meta, file("versions.yml"))`

#### Multiple Files: `files()`
- Returns: `Set<Path>`
- Used in: `Tuple<Map, Set<Path>>`
- Example: `logs = tuple(meta, files("*.{log,err}", optional: true))`

**Key Rule**: Use `file()` for single, known files. Use `files()` for wildcards, optional files, or multiple files.

### Channel Type Declarations
Standard pattern for workflow-level channels:

```nextflow
ch_results = channel.empty() as Channel<Tuple<Map, Path>>
ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
ch_versions = channel.empty() as Channel<Tuple<Map, Path>>
```

### Module Input Types
- Single assemblies: `Tuple<Map, Path>`
- Read pairs (R1/R2): `Tuple<Map, List<Path>>`
- File collections: `Tuple<Map, Set<Path>>`

## Design Patterns

### Entry Workflow Structure

All entry workflows follow this consistent pattern:

```nextflow
#!/usr/bin/env nextflow
nextflow.preview.types = true

params {
    rundir : String
    // Tool-specific parameters with types
}

include { BACTOPIATOOL } from '../../../subworkflows/utils/bactopia/main'
// Bactopia Tools workflows include BACTOPIATOOL_INIT
// include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { TOOL_NAME } from '../../../subworkflows/{tool}/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    TOOL_NAME(BACTOPIATOOL_INIT.out.samples, ...)

    // Collect outputs
    ch_results = ch_results.mix(TOOL_NAME.out.results)
    ch_logs = ch_logs.mix(TOOL_NAME.out.logs)
    ch_nf_logs = ch_nf_logs.mix(TOOL_NAME.out.nf_logs)
    ch_versions = ch_versions.mix(TOOL_NAME.out.versions)

    // Branch based on scope
    ch_final_results = ch_results.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }
    // ... similar for logs, nf_logs, versions

    publish:
    // Always includes: run_results, run_logs, run_nf_logs, run_versions, sample_results, sample_logs, sample_nf_logs, sample_versions
    run_results = ch_final_results.run
    sample_results = ch_final_results.sample
    // ... similar for other channels
}

output {
    // Run-level outputs
    run_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    // ... similar for other channels

    // Sample-level outputs
    sample_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    // ... similar for other channels
}
```

### Meta Map Structure

Standard meta fields used across all modules:

```groovy
meta.id = "${prefix}-${task.process}"
meta.name = prefix
meta.scope = task.ext.scope          // 'sample' or 'run'
meta.output_dir = "..."               // Used in output blocks
meta.logs_dir = "..."                 // Used in output blocks
meta.process_name = task.ext.process_name
```

### Channel Output Pattern for Subworkflows

All sub-workflows emit four standard "aggregate" output types:

1. **results**: Analysis results
2. **logs**: Tool execution logs
3. **nf_logs**: Nextflow execution logs
4. **versions**: Software version information

Each of these four must be converted from `Tuple<Map, Set<Path>>` to `Tuple<Map, Path>` using `flattenPaths` before being emitted.

Subworkflows must also emit the individual output channels from each included module/subworkflow.
These are emitted as-is without flattening.

### Channel Output Pattern for Modules

All modules emit three standard output types:

1. **results**: Analysis results
2. **logs**: Tool execution logs
3. **versions**: Software version information

Each of these three, makes use of `files` making them type `Tuple<Map, Set<Path>>`

Modules may emit additional output channels as needed. These will follow same `tuple` pattern,
but the type will depend on whether it is a single file (`file() -> Tuple<Map, Path>`) or
multiple files (files() -> `Tuple<Map, Set<Path>>`).

### Output Directory Pattern
Modules set output directories that are used in entry workflow output blocks:

```groovy
// In module with sample scope
meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"

// In module with run scope
meta.output_dir = "${task.ext.process_name}"
meta.logs_dir = "${task.ext.process_name}/logs/${task.ext.logs_subdir}/${task.ext.subdir}"

// In entry workflow output block
sample_results: Channel<Tuple<Map, Path>> {
    path { meta, _file -> "${meta.output_dir}/" }
}
```

## Component Tag Taxonomy

Bactopia uses a comprehensive tag system to classify components (modules, subworkflows, and workflows) based on their complexity, input/output patterns, and special features. These tags help quickly identify component characteristics and are documented in [`tag-taxonomy.md`](tag-taxonomy.md).

### Tag Categories

#### 1. Complexity Tags
Indicate the overall complexity of the component:
- `[complexity: simple]` - Basic tools with straightforward I/O (e.g., abricate)
- `[complexity: moderate]` - Tools with multiple options or moderate complexity (e.g., quast)
- `[complexity: complex]` - Tools with conditional logic, many options, or complex setup (e.g., prokka)

#### 2. Input Type Tags
Describe the input requirements:
- `[input-type: single]` - Single primary input
- `[input-type: multiple]` - Multiple separate inputs required

#### 3. Output Type Tags
Describe the output patterns:
- `[output-type: single]` - Single output file
- `[output-type: multiple]` - Multiple output files

#### 4. Feature Tags
Describe specific technical features or capabilities:

**Common Technical Features:**
- `path-workarounds` - Uses EMPTY_* files for optional parameters
- `conditional-input` - Accepts optional inputs
- `conditional-logic` - Contains if/else statements in script
- `custom-outputs` - Non-standard output channel patterns
- `compression` - Handles file compression/decompression
- `database-dependent` - Requires external databases

**Processing Features:**
- `aggregation` - Combines multiple results into summary
- `filtering` - Filters input data
- `validation` - Validates input format
- `directory-output` - Creates organized output directory structure
- `archive-output` - Creates compressed archives (tar/zip)
- `resource-download` - Downloads external databases, datasets, or files

### Implementation Format
Tags are placed in comment blocks at the top of main.nf files:
```nextflow
/* Tags: [complexity: simple] [input-type: single] [output-type: multiple] [features: database-dependent] */
```

### Usage Examples

**Simple Module (abricate):**
`[complexity: simple] [input-type: single] [output-type: multiple] [features: database-dependent]`

**Complex Module (prokka):**
`[complexity: complex] [input-type: single] [output-type: multiple] [features: path-workarounds, conditional-input, conditional-logic, compression, archive-output]`

**Aggregation Module (csvtk concat):**
`[complexity: simple] [input-type: multiple] [output-type: single] [features: aggregation]`

## Configuration System

### Configuration Hierarchy
The Bactopia pipeline uses a hierarchical configuration system that allows settings to be overridden at multiple levels:

    1. Global: nextflow.config
    2. Base: conf/base.config
    3. Profiles: conf/profiles.config
    4. Workflow: {workflow}/nextflow.config
    5. Module: {module}/process.config

### Main Configuration Files

#### `nextflow.config` (Root)
Global pipeline configuration including:
    - Process definitions and resource allocation
    - Profile configurations (conda, docker, singularity, aws, etc.)
    - Timeline and report settings
    - Manifest for Docker/Singularity images

#### `conf/base.config`
Base configuration settings applied to all executions:
    - Default process labels and container configurations
    - Basic parameter definitions
    - Default resource limits

#### `conf/profiles.config`
Execution environment profiles:
    - **conda**: Conda environment management
    - **docker**: Docker container execution
    - **singularity**: Singularity container execution
    - **aws**: AWS batch execution
    - **local**: Local machine execution
    - **standard**: Standard resource allocation

### Parameter System

#### Parameter Sources
Parameters are defined in multiple locations:
    1. `params.config` - Global defaults
    2. `{tool}/params.config` - Tool-specific defaults
    3. `{workflow}/nextflow_schema.json` - Workflow validation
    4. Command line overrides

#### Schema Validation
- `nextflow_schema.json` - Main pipeline schema
- `{workflow}/nextflow_schema.json` - Workflow-specific schemas
- Generated from module `schema.json` files
- Validates parameters at runtime

### Configuration Inheritance
Configuration values follow this precedence:
    1. Command line arguments (highest)
    2. Workflow-specific config
    3. Tool-specific params
    4. Global params.config
    5. Default values (lowest)

## Testing Framework

### nf-test Framework
Bactopia uses the nf-test framework for pipeline testing with these key components:

#### Test Organization
- `/tests/main.nf.test` - Main workflow tests
- `{workflow}/tests/` - Workflow-specific tests
- `.nf-test.config` - Test configuration
- `*.test.snap` - Test snapshots for output validation

#### Test Structure
```groovy
test("test_name") {
    when {
        process {
            """
            Test setup code
            """
        }
        then {
            assert workflow.completed
            assert workflow.success
            // Additional assertions
        }
    }
}
```

#### Test Categories
    1. **Unit Tests** - Individual module testing
    2. **Integration Tests** - Subworkflow testing
    3. **End-to-End Tests** - Full workflow testing
    4. **Regression Tests** - Output validation

### Test Data Management
- Test inputs stored in `tests/data/`
- Expected outputs in test snapshots
- Temporary test files managed by nf-test
- Test isolation through separate work directories

### Running Tests
```bash
# All tests
nf-test run

# Specific test
nf-test run tests/main.nf.test

# With specific profile
nf-test run -profile conda
```

## Development Workflow

### Adding a New Tool

#### 1. Create Module (if needed)

1. Create directory: `modules/{tool}/`
2. Create required files:
    - `main.nf` - Process definition
    - `meta.yaml` - Module metadata
    - `params.config` - Tool parameters
    - `process.config` - Resource configuration
    - `schema.json` - Parameter schema

#### 2. Create Subworkflow (if needed)

1. Create directory: `subworkflows/{tool}/`
2. Create required files:
    - `main.nf` - Subworkflow definition
    - `meta.yaml` - Subworkflow metadata

#### 3. Create Entry Workflow

1. Create directory: `workflows/{tool}/` for Named Workflows, or `workflows/bactopia-tools/{tool}/` for Bactopia Tools
2. Create required files:
    - `main.nf` - Entry workflow
    - `nextflow_schema.json` - Parameter schema
    - `nextflow.config` - Workflow configuration
    - `tests/` - Test suite

#### 4. Update Configuration
- Add tool to `conf/workflows.yaml`
- Include in appropriate tool categories
- Update citation information

#### 5. Add Tests
- Create test in `tests/` or workflow directory
- Define test cases and expected outputs
- Update test snapshots

### Development Checklist
- [ ] Module follows typing conventions
- [ ] Uses consistent meta map structure
- [ ] Emits all four standard channels
- [ ] Includes version tracking
- [ ] Has comprehensive tests
- [ ] Documentation updated
- [ ] Tested across profiles

### Code Quality Standards
- Follow existing patterns and conventions
- Use meaningful variable and process names
- Include inline comments for complex logic
- Validate all inputs and outputs
- Handle edge cases gracefully

## Known Limitations

### Path? Optional Parameters
Nextflow has incomplete support for optional Path parameters (`Path?`). This is a known bug
and is being worked on by the Nextflow team. Until this is resolved, Bactopia uses a workaround
for all optional Path? parameters.

**Workaround**: Use EMPTY_* placeholder files

```groovy
// This pattern is used throughout the codebase
bakta_db : Path? = "${projectDir}/data/empty/EMPTY_DB"  // TODO: Remove when Path? is fixed
```

**Detection in scripts**:

```groovy
def proteins_opt = proteins.getName() != "EMPTY_PROTEINS" ? "--proteins ${proteins.getName()}" : ""
```

**Important**: These are TEMPORARY workarounds. Do NOT attempt to "fix" them. They will be removed
when Nextflow properly supports Path?.

### EMPTY_* Files

When a module input is optional (Path?), an EMPTY_* placeholder file is used to indicate the
absence of a file. While they are all "empty files", please use the appropriate EMPTY_* file
for the expected file type.

Empty placeholder files are stored in `/data/empty/`:
- `EMPTY_ADAPTERS` - For adapter files
- `EMPTY_DB` - For database files
- `EMPTY_EXTRA` - For extra files
- `EMPTY_PHIX` - For PhiX files
- `EMPTY_PRODIGAL_TF` - For Prodigal training files
- `EMPTY_PROTEINS` - For protein files
- `EMPTY_R1` - For read 1 files
- `EMPTY_R2` - For read 2 files
- `EMPTY_REPLICONS` - For replicon files
- `EMPTY_TF` - For training files

The default value for optional Path? parameters should point to these files would be:

```groovy
param_name : Path? = "${projectDir}/data/empty/EMPTY_FILETYPE"
```

## Common Examples

### Module Template

```nextflow
nextflow.preview.types = true

process TOOL_NAME {
    tag "${prefix}"
    label 'process_single'  // or process_low, process_medium, process_high

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly): Tuple<Map, Set<Path>>

    output:
    report   = tuple(meta, file("${prefix}.txt"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create meta object
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    """
    # Tool execution commands
    tool command $assembly > ${prefix}.txt

    # Version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version 2>&1)
    END_VERSIONS
    """
}
```

### Subworkflow Interface

```nextflow

nextflow.preview.types = true

include { TOOL_A       } from '../tool_a/main'
include { TOOL_B       } from '../tool_b/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow SUBWORKFLOW_NAME {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    // Process channels
    TOOL_A(fasta)
    TOOL_B(gather(TOOL_A.out.results))

    emit:
    // Individual outputs
    tool_a_results: Channel<Tuple<Map, Path>> = TOOL_A.out.results
    tool_b_results: Channel<Tuple<Map, Set<Path>>> = TOOL_B.out.results

    // Aggregate outputs
    results = flattenPaths([
        TOOL_A.out.results,
        TOOL_B.out.results
    ])
    logs = flattenPaths([
        TOOL_A.out.logs,
        TOOL_B.out.logs
    ])
    nf_logs = flattenPaths([
        TOOL_A.out.nf_logs,
        TOOL_B.out.nf_logs
    ])
    versions = flattenPaths([
        TOOL_A.out.versions,
        TOOL_B.out.versions
    ])
}
```

## Important Notes

### Known Bugs
1. **Nextflow Path? Bug** - Optional Path parameters are not fully supported
2. **Workarounds are Temporary** - EMPTY_* files are used until the bug is fixed as tracked by TODO statements

### What NOT to Do
1. **Don't "fix" Path? workarounds** - They're necessary until Nextflow improves ([see explanation](#path-optional-parameters))
2. **Don't change `file()` to `files()`** - The difference is intentional ([see File Output Types](#file-output-types))
3. **Don't modify typing without understanding** - Types are deliberately chosen ([see Static Typing Conventions](#static-typing-conventions))
4. **Don't alter the subworkflow 4-channel pattern** - It's fundamental to the architecture ([see Channel Output Pattern for Subworkflows](#channel-output-pattern-for-subworkflows))

### Best Practices
1. **Follow existing patterns** - Don't reinvent unless necessary
2. **Use consistent meta fields** - Maintain the standard structure
3. **Document deviations** - If you must break a pattern, explain why
4. **Test with `nextflow config`** - Validate type annotations

### When Adding New Components
1. **Entry workflows**: Follow the abricate/main.nf template for Bactiopia Tools, and the Bactopia/main.nf template for Named workflows
2. **Subworkflows**: Emit all four standard channels
3. **Subworkflows**: Use `flattenPaths` for aggregate outputs and `gather` for inputs to summarizing modules
4. **Modules**: Use the module template above
5. **Modules**: Always create `process.config`, `params.config` and `schema.json`
6. **Modules**: Always emit logs, nf_logs, and versions outputs

### Common Misconceptions
- `Tuple<Map, Set<Path>>` vs `Tuple<Map, Path>`: Not an inconsistency, but based on using `files()` vs `file()`
- Channel type variations: Actually consistent across workflows
- Complex tuple patterns: Appropriately complex for the use case

## Version Requirements

### Nextflow Versions
- **Minimum**: v25.10.x
- **Recommended**: v25.10.x or later (for better Path? support)
- **Current limitation**: Path? requires workarounds until v26.04+

### Breaking Changes to Watch
- Future versions may remove need for EMPTY_* files
- Static typing may become default (no preview flag needed)

## Groovydoc Documentation Standards

Bactopia uses standardized Groovydoc templates for consistent documentation across all components. Three distinct templates exist for different component types:

### Component Type Distinctions
- **Modules**: Internal tool execution with `@output` channels
- **Subworkflows**: Internal orchestration with `@output` channels
- **Entry Workflows**: Published results for end users with `@publish` files

### Standard Fields
- **Required**: @status, @keywords, @citation
- **Optional**: @name (workflows only), @subworkflows, @modules, @note
- **Modules/Subworkflows**: @tags for classification
- **Entry Workflows**: @section to group related outputs

### Formatting Rules
- **No type annotations** in documentation (visible in code)
- **Multi-line inputs** for clarity
- **Single-line outputs** for readability
- **No comment markers** (e.g., // --- OUTPUTS ---)

## Module Template

For individual processes/modules that execute tools:

```groovy
/**
 * <Short summary of what this module does>.
 *
 * <Detailed description of the process>. You can use [Markdown links](url)
 * and multiple lines to explain the tool's purpose.
 *
 * @status <stable|beta|deprecated>
 * @keywords <comma, separated, keywords>
 * @tags complexity:<simple|moderate|complex> input-type:<single|multiple> output-type:<single|multiple> features:<features>
 * @citation <comma, separated, bibtex_keys>
 *
 * @note <Optional: Additional note about behavior/requirements>
 *
 * @input <channel_name>
 * <Description of a simple input>.
 *
 * @input tuple(<name1>, <name2>)
 * - `<name1>`: <Description of the first element>
 * - `<name2>`: <Description of the second element>
 *
 * @output <channel_name>    <Description of the output channel>
 * @output logs              Optional tool execution logs
 * @output nf_logs           Nextflow execution logs
 * @output versions          Software version information (YAML format)
 */
```

### Module Tag Categories
- **complexity**: `simple`, `moderate`, or `complex`
- **input-type**: `single` or `multiple`
- **output-type**: `single` or `multiple`
- **features**: Comma-separated list of applicable features:
  - `database-dependent`: Requires external database
  - `path-workarounds`: Uses EMPTY_* files for optional parameters
  - `conditional-input`: Accepts optional inputs
  - `conditional-logic`: Contains if/else statements
  - `compression`: Handles file compression/decompression
  - `archive-output`: Creates compressed archives
  - `resource-download`: Downloads external resources
  - `filtering`: Filters input data
  - `custom-outputs`: Non-standard output patterns
  - `aggregation`: Combines multiple results

## Subworkflow Template

For workflow components that orchestrate modules/subworkflows:

```groovy
/**
 * <Short summary of what this workflow does>.
 *
 * <Detailed description of the orchestration>. Can include conditional logic,
 * aggregation steps, and the relationship between included components.
 *
 * <Explain the flow, e.g. "Creates a pangenome using optional annotation with
 * Prokka, followed by optional phylogeny construction with IQ-TREE">.
 *
 * @status <stable|beta|deprecated>
 * @keywords <comma, separated, keywords>
 * @tags complexity:<simple|moderate|complex> input-type:<single|multiple> output-type:<multiple> features:<features>
 * @citation <comma, separated, bibtex_keys>
 *
 * @subworkflows <optional> <comma, separated, list_of_subworkflows>
 * @modules <optional> <comma, separated, list_of_modules>
 *
 * @note <Optional: Additional notes about workflow behavior>
 *
 * @input <channel_or_param_name>
 * <Description of input>.
 *
 * @input tuple(<name1>, <name2>)
 * - `<name1>`: <Description>
 * - `<name2>`: <Description>
 *
 * @output <specific_output> <Description of a specific, useful output channel>
 * @output results           Aggregated results channel containing all output files
 * @output logs              Aggregated logs channel containing all execution logs
 * @output nf_logs           Aggregated Nextflow execution logs from all processes
 * @output versions          Aggregated version information from all executed tools
 */
```

### Subworkflow Tag Categories
- **features**: Can include `aggregation`, `conditional-logic`, `components` (when using @subworkflows/@modules)

## Entry Workflow Template

For user-facing entry workflows (Bactopia Tools and Named Workflows):

```groovy
/**
 * Bactopia Tool: <Tool Name>.
 *
 * <Short summary of the tool>.
 *
 * <Detailed description>. Uses [Tool Name](URL) to <task>.
 *
 * @status <stable|beta|deprecated>
 * @keywords <comma, separated, keywords>
 * @citation <comma, separated, bibtex_keys>
 *
 * @subworkflows <optional> <comma, separated, list_of_subworkflows>
 * @modules <optional> <comma, separated, list_of_modules>
 *
 * @note <Optional: General tool notes/requirements>
 *
 * @input <param_name>
 * <Description of a specific parameter input>.
 *
 * @section <Group Name>
 * @note <Optional: Section-specific note, e.g. "Only created if --save_raw is used">
 * @publish <file_pattern>    <Description of the published file>
 * @publish <file_pattern>    <Description of the published file>
 */
```

## Quick Reference

```groovy
// Module: Use @output for channels
@output report     Abricate screening results
@output logs       Tool execution logs

// Subworkflow: Use @output for channels, include 4 standard channels
@output results    Aggregate of all result files
@output logs       Aggregate of all execution logs
@output nf_logs    Nextflow execution logs
@output versions   Software version information

// Entry Workflow: Use @publish for published files
@section Analysis Results
@publish *.vcf     Variant calls in VCF format
@publish *.txt     Summary statistics
```

## Examples

### Simple Module Example (Abricate)

```groovy
/**
 * Screen assemblies for antimicrobial resistance genes.
 *
 * This process screens assembled contigs against resistance gene databases
 * using [Abricate](https://github.com/tseemann/abricate). It quickly identifies
 * antimicrobial resistance, virulence genes, and other relevant sequences in
 * bacterial genomes.
 *
 * @status stable
 * @keywords bacteria, assembly, antimicrobial resistance, screening
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent
 * @citation abricate
 *
 * @input tuple(meta, assembly)
 * Input assembly and metadata. The meta map contains sample information
 * such as sample ID, name, and other metadata. The assembly set contains
 * the assembled contigs in FASTA format to be screened.
 *
 * @output report    Abricate screening results (TSV format with resistance gene hits)
 * @output logs      Optional tool log files containing execution information
 * @output nf_logs   Nextflow execution logs including command, stdout, stderr, and exit code
 * @output versions  Software version information in YAML format
 */
```

### Complex Module Example (Prokka)

```groovy
/**
 * Annotate bacterial genomes with functional information.
 *
 * This process annotates bacterial contigs or complete genomes using [Prokka](https://github.com/tseemann/prokka).
 * It rapidly calls genes, translates them, and searches them against multiple protein databases
 * to produce comprehensive annotation in various standard formats.
 *
 * @status stable
 * @keywords bacteria, annotation, genome, prokaryote, functional annotation
 * @tags complexity:complex input-type:single output-type:multiple features:path-workarounds conditional-input conditional-logic compression archive-output
 * @citation prokka
 *
 * @note Optional: Proteins and training files improve annotation quality
 *
 * @input tuple(meta, assembly)
 * Input assembly for annotation. The meta map contains sample information,
 * and the assembly set contains the assembled contigs in FASTA format.
 *
 * @input proteins
 * Optional protein sequences for homology search. When provided, these trusted
 * protein sequences are used to improve annotation accuracy through homology.
 *
 * @input prodigal_tf
 * Optional Prodigal training file. Species-specific training data that improves
 * gene prediction accuracy.
 *
 * @output annotations Complete annotation package containing GFF, GBK, FASTA, and other formats
 * @output gff         Genome annotation in GFF3 format (standard for genome browsers)
 * @output gbk         Genome annotation in GenBank format (rich format suitable for NCBI submission)
 * @output logs        Optional tool log files containing execution information
 * @output nf_logs     Nextflow execution logs including command, stdout, stderr, and exit code
 * @output versions    Software version information in YAML format
 */
```

### Subworkflow Example (Pangenome)

```groovy
/**
 * Perform pangenome analysis with optional core-genome phylogeny.
 *
 * This subworkflow creates a pangenome from GFF3 annotation files using one of three
 * tools: PIRATE (default), Roary, or Panaroo. It generates core-genome alignments
 * and gene presence/absence matrices, followed by SNP distance calculations using
 * snp-dists. The workflow conditionally executes the selected pangenome tool based
 * on Boolean parameters.
 *
 * @status stable
 * @keywords alignment, core-genome, pan-genome, phylogeny, comparative genomics
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation pirate, panaroo, roary, snpdists
 *
 * @subworkflows pirate, panaroo, roary, snpdists
 *
 * @note Optional: SNP distance calculation is always performed on core-genome alignment
 *
 * @input gff
 * GFF3 annotation files from assembled genomes. Each tuple contains metadata
 * about the sample and a set of GFF3 files representing the annotation.
 *
 * @input use_pirate
 * Use PIRATE for pangenome analysis (default). When true, executes PIRATE
 * which is suitable for highly diverse datasets.
 *
 * @output aln          Core-genome alignment file containing genes present across all input genomes
 * @output csv          Gene presence/absence matrix showing which genes are present in each genome
 * @output results      Aggregate of all result files from pangenome analysis and SNP distances
 * @output logs         Aggregate of all log files from executed tools
 * @output nf_logs      Nextflow execution logs from all processes
 * @output versions     Software version information from all executed tools
 */
```

### Entry Workflow Example (Pangenome Bactopia Tool)

```groovy
/**
 * Bactopia Tool: Pangenome Analysis.
 *
 * Performs comprehensive pangenome analysis from bacterial genomes.
 * Creates gene presence/absence matrices and builds phylogenetic trees.
 *
 * @status stable
 * @keywords pangenome, comparative genomics, phylogeny
 * @citation pirate, panaroo, roary, iqtree, scoary, clonalframeml
 *
 * @subworkflows ncbigenomedownload, prokka, pangenome, clonalframeml, iqtree, scoary
 *
 * @note Optional: Requires trait file for GWAS analysis with SCOARY
 *
 * @input scoary_traits
 * Path to trait file for genome-wide association studies.
 *
 * @section Pangenome Analysis
 * @publish *.csv         Gene presence/absence matrix
 * @publish *.aln         Core-genome alignment
 * @publish *.tree        Phylogenetic tree file
 *
 * @section GWAS Analysis
 * @note Only created if --scoary_traits is specified
 * @publish scoary/*.csv  Association results between genes and traits
 * @publish scoary/*.txt  Statistical summary of GWAS results
 *
 * @section Recombination Analysis
 * @note Only created if --skip_recombination is false
 * @publish *.masked.aln Recombination-masked alignment
 * @publish clonalframe/*.json ClonalFrameML analysis results
 *
 * @section Versions
 * @publish versions.yml   Software version information
 */
```

## Glossary

- **Bactopia Tools**: Workflows that require Bactopia outputs (e.g., abricate, prokka)
- **Named Workflows**: Standalone workflows independent of Bactopia (e.g., clean-yer-reads, staphopia)
- **Scope**: Either 'sample' (per-sample) or 'run' (aggregate)
- **Meta map**: Standardized metadata structure with fields like id, name, scope, output_dir, logs_dir
- **4-channel pattern**: Standard output pattern (results, logs, nf_logs, versions) for subworkflows
- **Path? workarounds**: Temporary solutions using EMPTY_* placeholder files
- **flattenPaths**: Utility function to convert Set<Path> to Path for outputs
- **gather**: Utility function to merge channels for aggregate operations
- **nf-test**: Testing framework used by Bactopia

## Conclusion

The Bactopia pipeline has a well-designed, consistently implemented static typing system. What might appear as "inconsistencies" are often intentional design patterns that serve specific purposes.

When working with this codebase:
1. Trust the existing patterns
2. Document any deviations
3. Focus on maintaining consistency
4. Remember that workarounds are temporary and tracked

This reference should be consulted whenever making decisions about the codebase architecture or typing conventions.
