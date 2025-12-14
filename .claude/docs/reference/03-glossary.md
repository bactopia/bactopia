# Glossary

## B
**Bactopia Tools**
- Workflows that require outputs from the main Bactopia workflow or Named Workflows
- Examples: abricate, prokka, blastn
- Located in `workflows/bactopia-tools/`

**Bakta**
- A tool for rapid, standardized annotation of bacterial genomes
- Known for being resource-intensive (complex module)

## C
**Channel**
- Nextflow's primary data structure for connecting processes
- Types: `Channel<Tuple<Map, Path>>`, `Channel<Tuple<Map, Set<Path>>>`

**Conda**
- Package and environment management system
- One of the execution profiles available in Bactopia

## D
**Docker**
- Container platform for reproducible computing
- One of the execution profiles available in Bactopia

## E
**EMPTY_* Files**
- Placeholder files used as temporary workaround for Nextflow Path? limitations
- Stored in `/data/empty/`
- Examples: `EMPTY_PROTEINS`, `EMPTY_DB`, `EMPTY_ADAPTERS`

## F
**file() vs files()**
- `file()`: Returns single `Path`, used for known individual files
- `files()`: Returns `Set<Path>`, used for wildcards, optional files, or multiple files

**flattenPaths**
- Utility function from plugin/nf-bactopia
- Converts `Tuple<Map, Set<Path>>` to `Tuple<Map, Path>` for outputs
- Essential for subworkflow aggregate outputs

## G
**gather**
- Utility function from plugin/nf-bactopia
- Merges multiple channels for aggregate operations
- Used to combine inputs for summarizing modules

**GroovyDoc**
- Documentation format used across all Bactopia components
- Standardized templates for modules, subworkflows, and workflows

## M
**module.config**
- Combined configuration file for modules
- Contains both parameter defaults and process settings
- Replaces the previous separate `params.config` and `process.config` files

**Meta Map**
- Standardized metadata structure attached to all outputs
- Required fields: `id`, `name`, `scope`, `output_dir`, `logs_dir`, `process_name`
- Used for organizing and routing outputs

**Modules**
- Individual process implementations that execute specific tools
- Located in `/modules/`
- The basic building blocks of the pipeline
- Emit 3 standard channels: logs, nf_logs, versions

## N
**Named Workflows**
- Standalone workflows independent of Bactopia outputs
- Examples: clean-yer-reads, staphopia, teton
- Located in `/workflows/` (not in `bactopia-tools/`)

**Nextflow**
- Workflow management system for data-driven computational pipelines
- Requires v25.10.x or later for Bactopia

**nf-test**
- Testing framework used by Bactopia
- Provides comprehensive testing capabilities for workflows

**nf_logs**
- Standard output channel containing Nextflow execution information
- Includes `.command.*` files (command, stdout, stderr, exit code)

## P
**Path?**
- Optional Path parameter type in Nextflow
- Currently has incomplete support, requiring EMPTY_* workarounds
- Expected to be fixed in Nextflow v26.04+

**Process**
- Nextflow's abstraction for executing commands or scripts
- All tool execution happens in processes

## R
**results**
- Standard output channel for analysis results
- One of the 4 required channels for subworkflows

## S
**Scope**
- Either 'sample' (per-sample) or 'run' (aggregate)
- Determines output organization and processing level
- Set in meta map for each component

**Singularity**
- Container platform designed for HPC systems
- One of the execution profiles available in Bactopia

**Subworkflows**
- Reusable workflow components that combine modules
- Located in `/subworkflows/`
- Must emit exactly 4 standard channels

## T
**Tuple<Map, Path>**
- Single file output type
- Used with `file()` output declaration

**Tuple<Map, Set<Path>>**
- Multiple file output type
- Used with `files()` output declaration
- Used for module inputs when multiple files are expected

## V
**versions**
- Standard output channel for software version information
- Contains `versions.yml` file with YAML-formatted version data
- Required for all modules

## W
**Workflows**
- Entry point workflows that users directly interact with
- Located in `/workflows/`
- Handle user-facing parameters and orchestration

## Quick Reference

### Required Channels
- **Modules**: logs, nf_logs, versions (plus specific outputs)
- **Subworkflows**: results, logs, nf_logs, versions (always)
- **Workflows**: Various @publish files organized by @section

### Common Patterns
- **4-channel pattern**: Standard for subworkflows
- **Path? workarounds**: Use EMPTY_* files
- **flattenPaths**: Convert Set<Path> to Path
- **gather**: Merge channels for aggregation

### Component Types
- **Simple**: Linear execution, single tool
- **Moderate**: Branching logic, database dependency
- **Complex**: Orchestration, multiple aggregation steps

## See Also
- [Repository Structure](../project/01-repository-structure.md) - For directory organization
- [Technical Specifications](../standards/03-technical-specs.md) - For implementation details
- [Examples](../reference/01-examples.md) - For practical implementations