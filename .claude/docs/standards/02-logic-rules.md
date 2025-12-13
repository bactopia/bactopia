# Logic and Taxonomy

## Overview
This guide defines the decision-making logic and taxonomy used to classify Bactopia components. It provides clear rubrics for determining complexity, input/output types, and feature classifications.

## Complexity Rubric

### Simple
- **Definition**: Linear execution, single tool wrapper
- **Characteristics**:
  - Straightforward input processing
  - No conditional branching
  - Minimal parameter configuration
  - Single tool execution
- **Examples**: abricate, ariba, blast

### Moderate
- **Definition**: Branching logic, database dependency, or basic aggregation
- **Characteristics**:
  - Conditional tool selection based on parameters
  - Requires external database
  - Basic result aggregation
  - Multiple processing steps but linear flow
- **Examples**: kraken2, quast, mlst

### Complex
- **Definition**: Orchestration of subworkflows, massive databases, or multiple distinct aggregation steps
- **Characteristics**:
  - Coordinates multiple subworkflows/modules
  - Complex conditional logic
  - Multiple aggregation steps
  - Large database requirements (e.g., Bakta)
  - Intricate output management
- **Examples**: prokka, pangenome, bactopia (main workflow)

## Input-Type Logic

### Single Input
- **Definition**: The `take` block defines exactly **1 Channel**
- **Note**: Do not count `path` or `string` parameters
- **Examples**:
  ```nextflow
  workflow EXAMPLE {
      take:
      fasta: Channel<Tuple<Map, Set<Path>>>  // Only 1 channel

      // String/Path parameters don't count
      db: Path
      species: String
  }
  ```

### Multiple Inputs
- **Definition**: The `take` block defines **2 or more Channels**
- **Examples**:
  ```nextflow
  workflow EXAMPLE {
      take:
      fasta: Channel<Tuple<Map, Set<Path>>>     // Channel 1
      gff: Channel<Tuple<Map, Set<Path>>>       // Channel 2
      parameters: Channel<Map>                  // Channel 3
  }
  ```

## Output-Type Logic

### Single Output
- **Definition**: Component produces one primary output file/channel
- **Typical for**: Simple conversion tools
- **Example**: Single report file

### Multiple Outputs
- **Definition**: Component produces multiple distinct outputs
- **Typical for**::
  - Annotation tools (multiple file formats)
  - Analysis tools (results + logs + versions)
  - Subworkflows (aggregated results)

## Feature Classification

### Technical Features

#### database-dependent
- Component requires external database to function
- Examples: abricate, kraken2, mlst
- **Key indicators**: Database download/selection parameters

#### path-workarounds
- Uses EMPTY_* files for optional Path? parameters
- Temporary workaround for Nextflow Path? limitations
- **Key indicators**: Parameters with default EMPTY_* values
- Examples: prokka, bakta

#### conditional-input
- Accepts optional inputs that modify behavior
- **Key indicators**: Optional parameters that change processing
- Examples: Optional protein files, training data

#### conditional-logic
- Contains if/else statements in script
- Branching execution based on parameters
- **Key indicators**: Multiple execution paths
- Examples: pangenome tool selection, optional analysis steps

#### compression
- Handles file compression/decompression
- **Key indicators**: .gz file handling, compression parameters

#### archive-output
- Creates compressed archives (tar/zip)
- **Key indicators**: Archive creation, tar/zip outputs

#### resource-download
- Downloads external databases, datasets, or files
- **Key indicators**: Download steps, external resource fetching

#### filtering
- Filters input data based on criteria
- **Key indicators**: Filtering parameters, subset selection

#### custom-outputs
- Non-standard output channel patterns
- **Key indicators**: Complex output structure beyond standard 4-channel

### Processing Features

#### aggregation
- Combines multiple results into summary
- **Key indicators**: gather/flattenPaths usage, result merging
- Examples: csvtk concat, summary reports

## Component Classification Guide

### Modules
- Always have `logs`, `nf_logs`, `versions` outputs
- Additional outputs based on tool function
- Complexity determined by tool requirements and script logic

### Subworkflows
- Always emit 4 standard channels: `results`, `logs`, `nf_logs`, `versions`
- Must use `flattenPaths` for aggregate outputs
- Must use `gather` for inputs to summarizing modules
- Complexity based on orchestration complexity

### Entry Workflows
- Organize outputs by `@section` groups
- Use `@publish` for end-user files
- Handle both run-level and sample-level outputs

## Decision Trees

### Determining Complexity
1. Does it coordinate multiple components? → **Complex**
2. Does it have conditional branching? → **Moderate** to **Complex**
3. Does it require a database? → **Moderate**
4. Is it a simple tool wrapper? → **Simple**

### Determining Input-Type
1. Count channels in `take` block
2. Ignore `path` and `string` parameters
3. 1 channel → **Single**
4. 2+ channels → **Multiple**

### Determining Output-Type
1. Count distinct output files/channels
2. 1 primary output → **Single**
3. Multiple outputs → **Multiple**

## See Also
- [Style Guide](../standards/01-style-guide.md) - For template formats
- [Technical Specifications](../standards/03-technical-specs.md) - For implementation details