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