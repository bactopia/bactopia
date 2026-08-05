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
- Types: `Channel<Record>` for typed module/subworkflow data; `Channel<Path>` / `Channel<Map>` for utility flows

**Conda**
- Package and environment management system
- One of the execution profiles available in Bactopia

## D
**Docker**
- Container platform for reproducible computing
- One of the execution profiles available in Bactopia

## E
**EMPTY_* Files (removed)**
- Previously used as placeholder files to work around Nextflow `Path?` limitations
- Replaced with native `Path?` types and null checks in Nextflow 26.04+
- Optional parameters now use `Path?` in take blocks and null guards in scripts

## F
**file() vs files()**
- `file()`: Returns single `Path`, used for known individual files
- `files()`: Returns `Set<Path>`, used for wildcards, optional files, or multiple files

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
- Emit a single `record(...)` with standard fields (`meta`, `results`, `logs`, `nf_logs`, `versions`) plus tool-specific outputs

## N
**Named Workflows**
- Standalone workflows independent of Bactopia outputs
- Examples: clean-yer-reads, staphopia, teton
- Located in `/workflows/` (not in `bactopia-tools/`)

**Nextflow**
- Workflow management system for data-driven computational pipelines
- Requires Nextflow v26 or later (per `nextflowVersion` in `nextflow.config`)

**nf-test**
- Testing framework used by Bactopia
- Provides comprehensive testing capabilities for workflows

**nf_logs**
- Standard output channel containing Nextflow execution information
- Includes `.command.*` files (command, stdout, stderr, exit code)

## P
**Path?**
- Optional Path parameter type in Nextflow
- Natively supported in Nextflow 26.04+ with null checks
- Used in take blocks for optional file inputs (e.g., `proteins: Path?`)

**Process**
- Nextflow's abstraction for executing commands or scripts
- All tool execution happens in processes

## R
**Record**
- Statically-typed value class introduced in Nextflow 26.04 preview
- Primary container for module/subworkflow inputs and outputs (replaces the pre-typed tuple-of-map-and-paths pattern)
- Constructed once and immutable; workflow-dependent fields use inline ternaries inside the constructor

**results**
- Standard field on every module record, carrying publishable analysis output paths

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
- Emit two record channels: `sample_outputs` (per-sample record passthrough) and `run_outputs` (run-scoped aggregated record)

## V
**versions**
- Standard field on every module record, holding the YAML-formatted `versions.yml` content
- Required for all modules

## W
**Workflows**
- Entry point workflows that users directly interact with
- Located in `/workflows/`
- Handle user-facing parameters and orchestration

## Quick Reference

### Required Emissions
- **Modules**: single `record(...)` with `meta`, `results`, `logs`, `nf_logs`, `versions` plus tool-specific fields
- **Subworkflows**: `sample_outputs` (per-sample record passthrough) + `run_outputs` (aggregated record)
- **Workflows**: Various `@publish` files organized by `@section`

### Common Patterns
- **Record**: Statically-typed, immutable container for module/subworkflow I/O
- **Path?**: Native null support for optional inputs
- **gather**: Merge channels for aggregation

### Component Types
- **Simple**: Linear execution, single tool
- **Moderate**: Branching logic, database dependency
- **Complex**: Orchestration, multiple aggregation steps

## See Also
- [Repository Structure](../project/01-repository-structure.md) - For directory organization
- [Technical Specifications](../standards/03-technical-specs.md) - For implementation details
- [Examples](../reference/01-examples.md) - For practical implementations
