# Bactopia Pipeline Reference for Claude

This document serves as the AI Context Master Map for the Bactopia pipeline. It provides entry points to modular documentation for understanding the codebase structure, patterns, and conventions.

## Project Overview

Bactopia is a flexible pipeline for bacterial genome analysis. It follows a three-tier architecture:
- **Workflows** (Tier 1): User-facing entry points
- **Subworkflows** (Tier 2): Reusable orchestration components
- **Modules** (Tier 3): Individual tool implementations

The pipeline uses standardized GroovyDoc documentation and static typing throughout all components.

## Documentation Index

### Standards and Conventions
- **[Style Guide & Templates](.claude/docs/standards/01-style-guide.md)**
    - *Read this for*: GroovyDoc templates, header format, and tag ordering
    - Visual formatting rules for all component types

- **[Logic & Taxonomy](.claude/docs/standards/02-logic-rules.md)**
    - *Read this for*: Determining complexity, input/output types
    - Decision-making logic for component classification

- **[Technical Specifications](.claude/docs/standards/03-technical-specs.md)**
    - *Read this for*: Variable naming, type conventions, Path? workarounds
    - Implementation details and "gotchas"

- **[Subworkflow Documentation](.claude/docs/standards/04-subworkflow-documentation.md)**
    - *Read this for*: Complete methodology for documenting subworkflows
    - Step-by-step process with examples and best practices

- **[Module Documentation](.claude/docs/standards/05-module-documentation.md)**
    - *Read this for*: Complete methodology for documenting modules
    - Detailed patterns and examples for individual tool implementations

- **[Workflow Documentation](.claude/docs/standards/06-workflow-documentation.md)**
    - *Read this for*: Complete methodology for documenting entry workflows
    - User-facing documentation patterns with @publish and @section organization

- **[Tier Architecture](.claude/docs/standards/07-tier-architecture.md)**
    - *Read this for*: Formalized rules for workflows, subworkflows, and modules
    - Tier responsibilities, allowed operations, plugin functions, ext system, catalog.json

### Project Documentation
- **[Repository Structure](.claude/docs/project/01-repository-structure.md)**
    - *Read this for*: Directory organization and three-tier architecture
    - Physical layout of the codebase

- **[Development Workflow](.claude/docs/project/02-development-workflow.md)**
    - *Read this for*: Adding new tools and components
    - Step-by-step development guide with checklist

- **[Configuration System](.claude/docs/project/03-configuration-system.md)**
    - *Read this for*: Understanding parameter hierarchy
    - Configuration inheritance and profile management

- **[Testing Framework](.claude/docs/project/04-testing-framework.md)**
    - *Read this for*: Writing and running tests
    - nf-test framework usage and patterns

### Reference Materials
- **[Examples](.claude/docs/reference/01-examples.md)**
    - *Read this for*: Concrete implementation examples
    - Annotated examples of modules, subworkflows, and workflows

- **[Troubleshooting](.claude/docs/reference/02-troubleshooting.md)**
    - *Read this for*: Common error solutions
    - Debugging tips and problem resolution

- **[Glossary](.claude/docs/reference/03-glossary.md)**
    - *Read this for*: Definitions of Bactopia-specific terms
    - Quick reference for terminology and concepts

- **[Plugin Functions](.claude/docs/reference/04-plugin-functions.md)**
    - *Read this for*: Understanding `gather()` and `flattenPaths()` functions
    - Channel manipulation utilities from nf-bactopia plugin

- **[task.ext Properties](.claude/docs/reference/05-task-ext-properties.md)**
    - *Read this for*: Configuring module behavior via task.ext
    - Complete reference for all task.ext properties used in module.config files

## AI Agent Instructions

When working with this codebase:

1. **Read this Master Map first** to understand the structure
2. **Load only modules relevant to your current task** to maintain context efficiency
3. **For documenting modules**: Read [.claude/docs/standards/05-module-documentation.md](.claude/docs/standards/05-module-documentation.md) for complete methodology and examples
4. **For documenting subworkflows**: Read [.claude/docs/standards/04-subworkflow-documentation.md](.claude/docs/standards/04-subworkflow-documentation.md) for complete methodology and examples
5. **Always check** [.claude/docs/standards/03-technical-specs.md](.claude/docs/standards/03-technical-specs.md) for variable naming and technical conventions

## Quick Reference

### Common Tasks

**Adding a new tool**:
1. Read [Development Workflow](.claude/docs/project/02-development-workflow.md)
2. Follow the step-by-step guide
3. Use templates from [Module Documentation](.claude/docs/standards/05-module-documentation.md)

**Debugging type errors**:
1. Check [Technical Specifications](.claude/docs/standards/03-technical-specs.md)
2. Review [Troubleshooting](.claude/docs/reference/02-troubleshooting.md)
3. Look for Path? workarounds

**Understanding architecture**:
1. Start with [Repository Structure](.claude/docs/project/01-repository-structure.md)
2. Review three-tier architecture
3. Study [Examples](.claude/docs/reference/01-examples.md)

### Key Patterns

**Module inputs**: Record-typed with named parameters (e.g., `(_meta: Map, assembly: Path): Record`)
**Module outputs**: Single `record()` with named fields (downstream) + generic fields (publishing)
**Subworkflow outputs**: Emit `sample_outputs` (module record passthrough) and `run_outputs` (aggregated)
**Optional parameters**: Use EMPTY_* files (temporary workaround)

### Important Reminders
- **Don't "fix" Path? workarounds** - They're intentional
- **Use `file()` for single files, `files()` for multiple**
- **Follow existing patterns** - Don't reinvent unless necessary
- **Always use 4 spaces for indentation** in all code blocks and lists, with the exception of YAML files which use 2 spaces
