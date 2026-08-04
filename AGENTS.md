# Bactopia Pipeline Reference for AI Agents

This document serves as the AI Context Master Map for the Bactopia pipeline, following the [agents.md](https://agents.md/) convention. It provides entry points to modular documentation for understanding the codebase structure, patterns, and conventions.

## Project Overview

Bactopia is a flexible pipeline for bacterial genome analysis. It follows a three-tier architecture:
- **Workflows** (Tier 1): User-facing entry points
- **Subworkflows** (Tier 2): Reusable orchestration components
- **Modules** (Tier 3): Individual tool implementations

The pipeline uses standardized GroovyDoc documentation and static typing throughout all components.

## Documentation Index

### Standards and Conventions
- **[Style Guide & Templates](.agents/docs/standards/01-style-guide.md)**
    - *Read this for*: GroovyDoc templates, header format, and tag ordering
    - Visual formatting rules for all component types

- **[Logic & Taxonomy](.agents/docs/standards/02-logic-rules.md)**
    - *Read this for*: Determining complexity, input/output types
    - Decision-making logic for component classification

- **[Technical Specifications](.agents/docs/standards/03-technical-specs.md)**
    - *Read this for*: Variable naming, type conventions, Path? optional inputs
    - Implementation details and conventions

- **[Subworkflow Documentation](.agents/docs/standards/04-subworkflow-documentation.md)**
    - *Read this for*: Complete methodology for documenting subworkflows
    - Step-by-step process with examples and best practices

- **[Module Documentation](.agents/docs/standards/05-module-documentation.md)**
    - *Read this for*: Complete methodology for documenting modules
    - Detailed patterns and examples for individual tool implementations

- **[Workflow Documentation](.agents/docs/standards/06-workflow-documentation.md)**
    - *Read this for*: Complete methodology for documenting entry workflows
    - User-facing documentation patterns with @publish and @section organization

- **[Tier Architecture](.agents/docs/standards/07-tier-architecture.md)**
    - *Read this for*: Formalized rules for workflows, subworkflows, and modules
    - Tier responsibilities, allowed operations, plugin functions, ext system, catalog.json

### Project Documentation
- **[Repository Structure](.agents/docs/project/01-repository-structure.md)**
    - *Read this for*: Directory organization and three-tier architecture
    - Physical layout of the codebase

- **[Development Workflow](.agents/docs/project/02-development-workflow.md)**
    - *Read this for*: Adding new tools and components
    - Step-by-step development guide with checklist

- **[Configuration System](.agents/docs/project/03-configuration-system.md)**
    - *Read this for*: Understanding parameter hierarchy
    - Configuration inheritance and profile management

- **[Testing Framework](.agents/docs/project/04-testing-framework.md)**
    - *Read this for*: Writing and running tests
    - nf-test framework usage and patterns

### Reference Materials
- **[Examples](.agents/docs/reference/01-examples.md)**
    - *Read this for*: Concrete implementation examples
    - Annotated examples of modules, subworkflows, and workflows

- **[Troubleshooting](.agents/docs/reference/02-troubleshooting.md)**
    - *Read this for*: Common error solutions
    - Debugging tips and problem resolution

- **[Glossary](.agents/docs/reference/03-glossary.md)**
    - *Read this for*: Definitions of Bactopia-specific terms
    - Quick reference for terminology and concepts

- **[Plugin Functions](.agents/docs/reference/04-plugin-functions.md)**
    - *Read this for*: Understanding `gather()` and `flattenPaths()` functions
    - Channel manipulation utilities from nf-bactopia plugin

- **[task.ext Properties](.agents/docs/reference/05-task-ext-properties.md)**
    - *Read this for*: Configuring module behavior via task.ext
    - Complete reference for all task.ext properties used in module.config files

- **[Skills](.agents/docs/reference/06-skills.md)**
    - *Read this for*: project-local skill inventory and when to invoke `skill-creator`
    - Catalog of AI tooling built on top of `bactopia-*` CLIs

## AI Agent Instructions

When working with this codebase:

1. **Read this Master Map first** to understand the structure
2. **Load only modules relevant to your current task** to maintain context efficiency
3. **For documenting modules**: Read [.agents/docs/standards/05-module-documentation.md](.agents/docs/standards/05-module-documentation.md) for complete methodology and examples
4. **For documenting subworkflows**: Read [.agents/docs/standards/04-subworkflow-documentation.md](.agents/docs/standards/04-subworkflow-documentation.md) for complete methodology and examples
5. **Always check** [.agents/docs/standards/03-technical-specs.md](.agents/docs/standards/03-technical-specs.md) for variable naming and technical conventions
6. **Always use the `bactopia-dev` conda env for all project tooling** — `ruff`, `bactopia-*` CLIs (`bactopia-lint`, `bactopia-test`, `bactopia-merge-schemas`, `bactopia-catalog`, `bactopia-citations`), and `nf-test`. Invoke via `conda run -n bactopia-dev <cmd>` (or activate the env first). Never report a check as SKIP because a tool is "not on PATH" without trying this env.

## Quick Reference

### Common Tasks

**Adding a new tool**:
1. Read [Development Workflow](.agents/docs/project/02-development-workflow.md)
2. Follow the step-by-step guide
3. Use templates from [Module Documentation](.agents/docs/standards/05-module-documentation.md)

**Debugging type errors**:
1. Check [Technical Specifications](.agents/docs/standards/03-technical-specs.md)
2. Review [Troubleshooting](.agents/docs/reference/02-troubleshooting.md)
3. Look for Path? optional input patterns

**Understanding architecture**:
1. Start with [Repository Structure](.agents/docs/project/01-repository-structure.md)
2. Review three-tier architecture
3. Study [Examples](.agents/docs/reference/01-examples.md)

**Creating or editing a skill**:
1. Use the `skill-creator` skill — do not hand-scaffold `SKILL.md` files
2. See [Skills](.agents/docs/reference/06-skills.md) for the project's skill conventions and inventory

### Key Patterns

**Module inputs**: Record-typed with named parameters (e.g., `record(meta: Record, fna: Path)`)
**Module outputs**: Single `record()` with named fields (downstream) + generic fields (publishing)
**Subworkflow outputs**: Emit `sample_outputs` (module record passthrough) and `run_outputs` (aggregated)
**Optional parameters**: Use `Path?` types with `?` suffix in GroovyDoc

### Important Reminders
- **Use `file()` for single files, `files()` for multiple**
- **Follow existing patterns** - Don't reinvent unless necessary
- **Always use 4 spaces for indentation** in all code blocks and lists, with the exception of YAML files which use 2 spaces
