# Repository Structure

## Overview
The Bactopia repository follows a well-organized three-tier architecture that separates concerns between entry points (workflows), reusable components (subworkflows), and individual tool implementations (modules).

## Directory Structure

```
bactopia/
├── .ai-context/               # AI context documentation (modular)
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
  - `main.nf` - Process definition
  - `meta.yaml` - Module metadata
  - `params.config` - Module-specific parameters
  - `process.config` - Process resource configurations
  - `schema.json` - Parameter schema for validation
- **Examples**: `abricate/`, `prokka/`, `kraken2/`

### `/subworkflows/` (Tier 2)
- **Purpose**: Reusable workflow components that combine modules
- **Categories**:
  - `bactopia/` - Core pipeline functionality
  - `utils/` - Helper workflows (BACTOPIA_INIT, BACTOPIATOOL_INIT)
  - `{tool}/` - Tool-specific processing logic
- **Contents**:
  - `main.nf` - Subworkflow definition
  - `meta.yaml` - Subworkflow metadata
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
  - `workflows.yaml` - Workflow metadata

### `/data/`
- **Purpose**: Static data and resources
- **Contents**:
  - `empty/` - Placeholder files for Path? workarounds
  - `conda/` - Development environment specifications
  - Reference databases and files

### `/bin/`
- **Purpose**: Utility scripts and helper tools
- **Examples**:
  - `bactopia` - Main CLI tool
  - `check-fastqs.py` - FASTQ validation
  - `kraken-bracken-summary.py` - Result aggregation

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

### Modules & Subworkflows
- Primary file: `main.nf`
- Metadata: `meta.yaml`
- Configuration: `params.config`, `process.config`
- Schema: `schema.json` (modules only)

## Organization Patterns

### Module Organization
- Single tool per directory
- Multiple variants as subdirectories when needed
  - Example: `blast/blastn/`, `blast/blastp/`, `blast/blastx/`

### Subworkflow Organization
- By functionality or tool
- Core Bactopia components in `subworkflows/bactopia/`
- Utilities in `subworkflows/utils/`

### Workflow Organization
- Entry points in root directory
- Bactopia Tools in `workflows/bactopia-tools/`
- Named workflows as top-level directories

## See Also
- [Development Workflow](../project/02-development-workflow.md) - For adding new components
- [Technical Specifications](../standards/03-technical-specs.md) - For implementation details
- [Configuration System](../project/03-configuration-system.md) - For configuration hierarchy