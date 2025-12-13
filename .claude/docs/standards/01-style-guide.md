# GroovyDoc Style Guide

## Overview
This guide defines the standardized GroovyDoc templates and formatting conventions used across all Bactopia components for consistent documentation.

## Component Type Distinctions
- **Modules**: Internal tool execution with `@output` channels
- **Subworkflows**: Internal orchestration with `@output` channels
- **Entry Workflows**: Published results for end users with `@publish` files

## Standard Fields
- **Required**: `@status`, `@keywords`, `@citation`
- **Optional**: `@name` (workflows only), `@subworkflows`, `@modules`, `@note`
- **Modules/Subworkflows**: `@tags` for classification
- **Entry Workflows**: `@section` to group related outputs

## Formatting Rules
- **No type annotations** in documentation (visible in code)
- **Multi-line inputs** for clarity
- **Single-line outputs** for readability
- **No comment markers** (e.g., `// --- OUTPUTS ---`)
- **Links**: Always HTTPS to source repository (GitHub), not documentation
- **Alignment**: Pad output names so descriptions start at the same column

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
 * @output logs              Optional software execution logs containing warnings/errors
 * @output nf_logs           Nextflow execution scripts and logs for debugging
 * @output versions          A YAML formatted file with software versions
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
 * @output nf_logs           Aggregated Nextflow execution scripts and logs for debugging from all processes
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
@output nf_logs    Nextflow execution scripts and logs for debugging
@output versions   Software version information

// Entry Workflow: Use @publish for published files
@section Analysis Results
@publish *.vcf     Variant calls in VCF format
@publish *.txt     Summary statistics
```

## See Also
- [Logic Rules](../standards/02-logic-rules.md) - For determining complexity and input/output types
- [Technical Specifications](../standards/03-technical-specs.md) - For variable naming conventions