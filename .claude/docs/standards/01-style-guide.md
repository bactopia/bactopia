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

### The @note Tag
Use `@note` for special requirements, caveats, or important information:
- Database requirements
- Optional features
- Conditional behavior
- Important caveats

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
 * @input record(meta, <input_name>)
 * - `meta`: Groovy Map containing sample information
 * - `<input_name>`: <Description of the input file>
 *
 * @input <param_name>
 * <Description of a simple input parameter>.
 *
 * @output record(meta, <field1>, <field2>, results, logs, nf_logs, versions)
 * - `<field1>`: <Description of tool-specific output field>
 * - `<field2>`: <Description of tool-specific output field>
 */
```

### Module Tag Categories
- **complexity**: `simple`, `moderate`, or `complex`
- **input-type**: `single` or `multiple`
- **output-type**: `single` or `multiple`
- **features**: Comma-separated list (NO SPACES after commas) of applicable features:
  - `database-dependent`: Requires external database
  - `path-workarounds`: Uses EMPTY_* files for optional parameters
  - `conditional-input`: Accepts optional inputs
  - `conditional-logic`: Contains if/else statements
  - `compression`: Handles file compression/decompression
  - `archive-output`: Creates compressed archives
  - `resource-download`: Downloads external resources
  - `internet-access`: Requires internet connection for downloads
  - `alternative-execution`: Multiple tool options (e.g., assembler can use shovill, dragonflye, unicycler)
  - `filtering`: Filters input data
  - `custom-outputs`: Non-standard output patterns
  - `aggregation`: Combines multiple results

**Example**: `features:database-dependent,conditional-logic,compression` (NOT `features:database-dependent, conditional-logic, compression`)

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
 * @input record(meta, <input_name>)
 * - `meta`: Groovy Map containing sample information
 * - `<input_name>`: <Description of the input>
 *
 * @input <param_name>
 * <Description of a simple input parameter>.
 *
 * @output sample_outputs
 * - `<field>`: <Description of tool-specific output field from the module record>
 *
 * @output run_outputs
 * - `<field>`: <Description of aggregated output from CSVTK_CONCAT>
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
// Module: Use @output record(...) listing all fields, describe only tool-specific ones
@output record(meta, report, results, logs, nf_logs, versions)
- `report`: Abricate screening results

// Subworkflow: Emit sample_outputs (module record) and run_outputs (aggregated record)
@output sample_outputs
- `report`: Per-sample Abricate screening results
@output run_outputs
- `csv`: A merged TSV file with results from all samples

// Entry Workflow: Use @publish for published files
@section Analysis Results
@publish *.vcf     Variant calls in VCF format
@publish *.txt     Summary statistics
```

## See Also
- [Logic Rules](../standards/02-logic-rules.md) - For determining complexity and input/output types
- [Technical Specifications](../standards/03-technical-specs.md) - For variable naming conventions