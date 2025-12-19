# Development Workflow

## Overview
This guide outlines the standard workflow for developing and adding components to the Bactopia pipeline, ensuring consistency and maintainability across the codebase.

## Adding a New Tool

### Step 1: Create Module (if needed)

The module is the basic building block that executes a specific tool.

1. **Create directory**:
   ```bash
   mkdir modules/{tool_name}
   ```

2. **Create required files**:
   - `main.nf` - Process definition with GroovyDoc documentation
   - `module.config` - Module parameters and process configuration
   - `schema.json` - Parameter validation schema

3. **Key requirements**:
   - Include `nextflow.preview.types = true` at the top
   - Define inputs appropriately (`Tuple<Map, Path>` for single files, `Tuple<Map, Set<Path>>` for multiple)
   - Emit `logs`, `nf_logs`, and `versions` channels
   - Include version tracking in `versions.yml`

### Step 2: Create Subworkflow (if needed)

The subworkflow orchestrates one or more modules to provide higher-level functionality.

1. **Create directory**:
   ```bash
   mkdir subworkflows/{tool_name}
   ```

2. **Create required files**:
   - `main.nf` - Subworkflow definition with GroovyDoc documentation

3. **Key requirements**:
   - Include necessary modules/subworkflows
   - Always include `flattenPaths` and `gather` from plugin
   - Emit exactly 4 channels: `results`, `logs`, `nf_logs`, `versions`
   - Use `flattenPaths` for aggregate outputs

### Step 3: Create Entry Workflow

The entry workflow is what users interact with directly.

1. **Determine workflow type**:
   - **Bactopia Tool**: Requires Bactopia outputs → `workflows/bactopia-tools/{tool}/`
   - **Named Workflow**: Standalone → `workflows/{tool}/`

2. **Create directory structure**:
   ```bash
   mkdir workflows/{bactopia-tools/}{tool_name}/
   mkdir workflows/{bactopia-tools/}{tool_name}/tests/
   ```

3. **Create required files**:
   - `main.nf` - Entry workflow script
   - `nextflow_schema.json` - Parameter validation schema
   - `nextflow.config` - Workflow-specific configuration
   - `tests/main.nf.test` - Workflow tests

4. **Key requirements**:
   - Start with `#!/usr/bin/env nextflow`
   - Include appropriate initialization (BACTOPIA_INIT or BACTOPIATOOL_INIT)
   - Define all 4 output channels and branch by scope
   - Include publish block with run/sample outputs

### Step 4: Update Configuration

1. **Add to workflow registry**:
   - Edit `conf/workflows.yaml`
   - Add tool metadata and categorization
   - Include citation information

2. **Verify configuration inheritance**:
   - Global: `nextflow.config`
   - Base: `conf/base.config`
   - Profiles: `conf/profiles.config`
   - Workflow: `{workflow}/nextflow.config`

### Step 5: Add Tests

1. **Create test structure**:
   ```groovy
   test("tool_test_name") {
       when {
           process {
               """
               # Test setup code
               """
           }
       }
       then {
           assert workflow.completed
           assert workflow.success
           # Add specific assertions
       }
   }
   ```

2. **Test data**:
   - Use test data from bactopia-tests repository
   - Set `BACTOPIA_TESTS` environment variable
   - Create comprehensive test cases

## Code Quality Standards

### Required Patterns

1. **Static Typing**:
   - All files must have type annotations
   - Use proper channel declarations
   - Consistent typing across components

2. **Meta Map Structure**:
   ```groovy
   meta.id = "${prefix}-${task.process}"
   meta.name = prefix
   meta.scope = task.ext.scope
   meta.output_dir = "..."
   meta.logs_dir = "..."
   meta.process_name = task.ext.process_name
   ```

3. **Channel Patterns**:
   - Modules: Use `files()` for outputs
   - Subworkflows: Always emit 4 standard channels
   - Workflows: Branch outputs by scope (run/sample)

4. **Documentation**:
   - Follow GroovyDoc standards
   - Include all required fields (@status, @keywords, @citation)
   - Add appropriate tags for classification

### Resource Labels

Modules must declare resource requirements using labels defined in `conf/base.config`. Labels scale with retry attempts.

#### Primary Labels

Choose ONE primary label based on typical resource needs:

| Label | CPUs | Memory | Time | Use For |
|-------|------|--------|------|---------|
| `process_single` | 1 | 4GB | 2h | Single-threaded tools that cannot parallelize |
| `process_low` | 4 | 8GB | 4h | Lightweight tools, quick operations |
| `process_medium` | 8 | 32GB | 12h | Moderate workloads, BLAST-based tools |
| `process_high` | 12 | 64GB | 24h | Memory-intensive, complex analyses |

#### Modifier Labels

Add these IN ADDITION to a primary label when needed:

| Label | Effect | Use For |
|-------|--------|---------|
| `process_long` | Time → 96h | Long-running jobs (pangenome, phylogeny) |
| `process_high_memory` | Memory → 128GB | Extreme memory needs (GTDB-Tk) |

#### Example Usage

```groovy
// In main.nf
process TOOL_NAME {
    label 'process_low'           // Primary label

    // ... rest of process
}

// For long-running + high memory (e.g., gtdbtk/classifywf)
process GTDBTK_CLASSIFYWF {
    label 'process_high'
    label 'process_high_memory'   // Adds to process_high

    // ... rest of process
}

// For pangenome analysis
process PANAROO_RUN {
    label 'process_high'
    label 'process_long'          // Extends time limit

    // ... rest of process
}
```

#### Modules by Label

**process_single** (single-threaded):
- `abricate/run`, `rgi/heatmap`

**process_low** (lightweight):
- `mlst`, `prokka`, `snippy/run`, `csvtk/concat`
- `mykrobe/predict`, `shigatyper`, `seroba/run`
- Most typing tools (pasty, spatyper, lissero, etc.)

**process_medium** (moderate):
- `bracken`, `blast/*`, `sistr`, `gubbins`
- `ismapper`, `checkm2/predict`, `snippy/core`
- `phispy`, `mobsuite/recon`

**process_high** (resource-intensive):
- `gtdbtk/classifywf`, `panaroo/run`

**process_long** (extended time):
- `panaroo/run`, `checkm2/download`, `gtdbtk/download`

#### Choosing the Right Label

1. **Start with `process_low`** for most tools
2. **Upgrade to `process_medium`** for BLAST-based, database searches, or multi-step tools
3. **Use `process_high`** only for known memory-intensive tools
4. **Add `process_long`** if jobs regularly exceed 12 hours
5. **Add `process_high_memory`** only if 64GB is insufficient

### Multi-Process Module Structure

When a tool requires multiple distinct operations, split it into subdirectories within the module folder. This maintains modularity and allows each process to have its own configuration.

#### When to Split

Split a module into subdirectories when:

1. **Separate execution stages**: Download database vs. run analysis
2. **Different scopes**: Per-sample processing vs. run-level aggregation
3. **Different resource requirements**: Lightweight summary vs. intensive computation
4. **Reusable components**: Processes that can be called independently

#### Directory Structure

```bash
modules/{tool_name}/
├── run/                    # Primary analysis process
│   ├── main.nf
│   ├── module.config
│   └── schema.json
└── summary/                # Aggregation process (if needed)
    ├── main.nf
    ├── module.config
    └── schema.json
```

#### Common Patterns

| Tool | Subdirectories | Purpose |
|------|----------------|---------|
| `abricate/` | `run/`, `summary/` | Per-sample screening → Run-level aggregation |
| `bakta/` | `download/`, `run/` | Database download → Annotation |
| `tbprofiler/` | `profile/`, `collate/` | Per-sample profiling → Run-level collation |
| `gtdbtk/` | `download/`, `classifywf/` | Database download → Classification |
| `checkm2/` | `download/`, `predict/` | Database download → Quality assessment |
| `snippy/` | `run/`, `core/` | Per-sample mapping → Core genome extraction |
| `rgi/` | `main/`, `heatmap/` | AMR detection → Visualization |

#### Configuration for Multi-Process Modules

Use shared and process-specific configuration blocks in `module.config`:

```groovy
process {
    // Shared configuration for all processes
    withName: 'TOOL_RUN|TOOL_SUMMARY' {
        ext.toolName = "bioconda-tool-1.0.0"
        ext.docker = "biocontainers/tool:1.0.0"
        ext.image = "https://depot.galaxyproject.org/singularity/tool:1.0.0"
        ext.condaDir = "${params.condadir}"
        ext.wf = params.wf
    }

    // Process-specific: per-sample
    withName: 'TOOL_RUN' {
        ext.scope = "sample"
        ext.process_name = "tool"
        ext.subdir = ""
    }

    // Process-specific: run-level
    withName: 'TOOL_SUMMARY' {
        ext.scope = "run"
        ext.process_name = params.merge_folder
        ext.prefix = "tool-summary"
    }
}
```

#### Naming Conventions

- **`run/`**: Primary per-sample analysis
- **`summary/` or `collate/`**: Run-level aggregation
- **`download/`**: Database or resource downloads
- **`main/`**: Primary process when other names don't fit
- Use lowercase, descriptive names that match the tool's terminology

### Best Practices

1. **Follow existing patterns** - Don't reinvent unless necessary
2. **Use consistent naming** - Maintain naming conventions
3. **Validate all inputs** - Include proper validation
4. **Handle edge cases** - Graceful error handling
5. **Include version tracking** - Always generate versions.yml

## Development Checklist

Before submitting a new tool:

- [ ] Module created with all required files
- [ ] Module follows typing conventions
- [ ] Module uses consistent meta map structure
- [ ] Module emits logs, nf_logs, and versions
- [ ] Subworkflow emits all 4 standard channels (if applicable)
- [ ] Entry workflow follows standard pattern
- [ ] Workflow schema validates correctly
- [ ] Tests created and passing
- [ ] Documentation updated (GroovyDoc)
- [ ] Tool added to workflow registry
- [ ] Citation information included
- [ ] Code follows style guidelines

## Common Anti-patterns to Avoid

1. **Hard-coding paths** - Use relative paths and meta.output_dir
2. **Skipping version tracking** - Always include versions.yml
3. **Inconsistent typing** - Use consistent types across connected components
4. **Missing standard channels** - Subworkflows must emit 4 channels
5. **Breaking the 3-tier architecture** - Don't include modules directly in workflows

## See Also
- [Testing Framework](../project/04-testing-framework.md) - For detailed testing guidelines
- [Technical Specifications](../standards/03-technical-specs.md) - For implementation details
- [Style Guide](../standards/01-style-guide.md) - For documentation standards