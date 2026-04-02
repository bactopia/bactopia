# Repository Structure

## Overview
The Bactopia repository follows a well-organized three-tier architecture that separates concerns between entry points (workflows), reusable components (subworkflows), and individual tool implementations (modules).

## Directory Structure

```
bactopia/
├── .claude/                   # AI context documentation (modular)
├── .github/                   # Github Actions workflows and issue templates
├── .vscode/                   # Visual Studio Code settings and configurations
├── bin/                       # Helper scripts and utilities
├── conf/                      # Configuration files
├── data/                      # Static data and resources
│   └── empty/                 # Empty placeholder files
├── modules/                   # Individual process implementations
│   └── bactopia/              # Core pipeline processes
├── subworkflows/              # Reusable workflow components
│   ├── bactopia/              # Essential pipeline components
│   └── utils/                 # Helper workflows
├── tests/                     # Test suite
├── workflows/                 # Entry point workflows
│   └── bactopia-tools/        # Workflows requiring Bactopia outputs
├── main.nf                    # Main Bactopia workflow
├── nextflow.config            # Global configuration
├── nextflow_schema.json       # Parameter validation schema
├── CLAUDE.md                  # AI Context Master Map
└── README.md                  # Project documentation
```

## Key Directories

### `/modules/` (Tier 3)
- **Purpose**: Individual process implementations that execute specific tools
- **Structure**: Each tool has its own subdirectory
- **Contents**:
  - `main.nf` - Process definition with GroovyDoc documentation
  - `module.config` - Module parameters and process configuration
  - `schema.json` - Parameter schema for validation
- **Count**: 96 modules
- **Examples**: `abricate/`, `prokka/`, `kraken2/`

### `/subworkflows/` (Tier 2)
- **Purpose**: Reusable workflow components that combine modules
- **Categories**:
  - `bactopia/` - Core pipeline functionality (contains nested subworkflows)
  - `utils/` - Helper workflows (contains nested subworkflows)
  - `{tool}/` - Tool-specific processing logic
- **Contents**:
  - `main.nf` - Subworkflow definition with GroovyDoc documentation
- **Count**: 87 subworkflows
- **Key Requirement**: Must always emit 4 channels (results, logs, nf_logs, versions)

### `/workflows/` (Tier 1)
- **Purpose**: Entry point workflows for users
- **Types**:
  - **Main Workflow**: `main.nf` - Primary Bactopia pipeline
  - **Bactopia Tools**: `workflows/bactopia-tools/` - Require Bactopia outputs
  - **Named Workflows**: Standalone workflows independent of Bactopia
- **Contents**:
  - `main.nf` - Entry workflow script
  - `nextflow_schema.json` - Workflow-specific parameter schema
  - `nextflow.config` - Workflow-specific configuration
  - `tests/` - Workflow-specific tests

### `/conf/`
- **Purpose**: Configuration files and parameter definitions
- **Key Files**:
    - `base.config` - Base configuration settings
    - `params.config` - Global parameter definitions
    - `profiles.config` - Execution profile configurations
    - `test.config` - Test configuration settings
- **Subdirectories**:
    - `params/` - Workflow-specific parameter configurations
        - `bactopia.config` - Main Bactopia parameters
        - `bactopia-tools.config` - Shared bactopia-tools parameters
        - `cleanyerreads.config` - CleanYerReads workflow parameters
        - `staphopia.config` - Staphopia workflow parameters
        - `teton.config` - Teton workflow parameters
    - `schema/` - JSON schemas for parameter validation
        - `bactopia.json`, `bactopia-tools.json`, `generic.json`, etc.

### `/data/`
- **Purpose**: Static data and resources
- **Contents**:
    - `empty/` - Placeholder files for Path? workarounds (EMPTY_*)
    - `conda/` - Development environment specifications
    - `citations.yml` - Tool citations and references
    - `workflows.yml` - Workflow metadata and definitions
    - `proteins.faa` - Protein reference file
    - Image assets (logos, banners)

### `/bin/`
- **Purpose**: CLI wrapper scripts for the bioconda `bactopia` package
- **Contents**:
  - `bactopia` - Main CLI wrapper (routes subcommands, runs Nextflow)
  - `clean-yer-reads`, `cyr`, `staphopia`, `teton` - Workflow alias wrappers
- **Note**: Pipeline utility scripts (FASTQ validation, coverage processing, etc.) have been migrated to bactopia-py as `bactopia-*` CLI entry points

## Three-Tier Architecture

### Data Flow
```
User Input
    ↓
Workflows (Tier 1) - Parameter handling, orchestration
    ↓
Subworkflows (Tier 2) - Combine related functionality
    ↓
Modules (Tier 3) - Execute individual tools
    ↓
Results (4 channels: results, logs, nf_logs, versions)
```

### Component Relationships
- **Workflows** include **Subworkflows**
- **Subworkflows** include **Modules**
- **Modules** are never directly included by **Workflows**
- This hierarchical approach promotes code reuse and maintainability

## File Naming Conventions

### Workflows
- Main entry: `main.nf`
- Config: `nextflow.config`
- Schema: `nextflow_schema.json`

### Modules
- Primary file: `main.nf` (includes GroovyDoc documentation)
- Configuration: `module.config`
- Schema: `schema.json`

### Subworkflows
- Primary file: `main.nf` (includes GroovyDoc documentation)

## Organization Patterns

### Module Organization
- Single tool per directory
- Multiple variants as subdirectories when needed
    - Example: `blast/blastn/`, `blast/blastp/`, `blast/blastx/`, `blast/tblastn/`, `blast/tblastx/`
- Some modules use a `/run/` subdirectory for the main process
    - Example: `abricate/run/`, `bakta/run/`, `prokka/` (direct)
    - Used when modules have multiple related processes (e.g., download + run)

### Subworkflow Organization
- By functionality or tool
- Core Bactopia components in `subworkflows/bactopia/`
- Utilities in `subworkflows/utils/`

### Workflow Organization
- Entry points in root directory
- Bactopia Tools in `workflows/bactopia-tools/` (66 workflows, 69 total across all tiers)
- Named workflows as separate directories under `workflows/`:
    - `workflows/cleanyerreads/` - Read cleaning workflow
    - `workflows/staphopia/` - Staphylococcus-focused analysis
    - `workflows/teton/` - Taxonomic classification workflow

## See Also
- [Development Workflow](../project/02-development-workflow.md) - For adding new components
- [Technical Specifications](../standards/03-technical-specs.md) - For implementation details
- [Configuration System](../project/03-configuration-system.md) - For configuration hierarchy