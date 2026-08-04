# Troubleshooting

## Overview
This guide provides solutions to common issues, error messages, and debugging tips for working with the Bactopia pipeline.

## Common Error Messages

### Type Errors

#### "Cannot cast X to Y"
**Cause**: Usually means using `file()` when you need `files()`
**Solution**:
- Check if you're returning multiple files when a single file is expected
- Use `file()` for single known files, `files()` for multiple or wildcard patterns
- Example fix: `record(meta: meta, results: [ files("*.txt") ])` — use `files()` (plural) when the glob can match multiple files

#### "Missing required parameter"
**Cause**: Parameter not defined in configuration or schema
**Solution**:
- Check `schema.json` and `params.config` files
- Ensure parameter is defined with correct type
- Verify parameter is being passed correctly

#### "Channel not found"
**Cause**: Workflow trying to access a channel that doesn't exist
**Solution**:
- Verify the subworkflow emits both `sample_outputs` (per-sample record passthrough) and `run_outputs` (aggregated run-scope record)
- Check channel names match exactly
- Ensure proper channel mixing/branching

### Path and File Issues

#### "No such file or directory"
**Cause**: File path doesn't exist or is incorrect
**Solution**:
- Verify relative vs absolute paths
- Ensure input files are properly staged
- For optional inputs, confirm `Path?` parameters are null when not provided

### Configuration Issues

#### "Profile not found"
**Cause**: Trying to use undefined execution profile
**Solution**:
- Check available profiles in `conf/profiles.config`
- Use standard profiles: conda, docker, singularity, local
- Create custom profile if needed

#### "Schema validation failed"
**Cause**: Parameter doesn't match schema definition
**Solution**:
- Check parameter type and value against `nextflow_schema.json`
- Ensure required parameters are provided
- Validate parameter format

## Debugging Tips

### 1. Check the Basics First
- Verify Nextflow version meets requirements (v26+)
- Check that all required files exist
- Ensure proper permissions on input/output directories

### 2. Verify Meta Record Fields
Ensure the meta record contains required fields:

```groovy
meta.id
meta.name
meta.scope  // 'sample' or 'run'
meta.output_dir
meta.logs_dir
meta.process_name
```

### 3. Check Consistent Typing
Verify consistent typing across connected components:
- `Channel<Record>` for module and subworkflow take blocks
- `record(meta: Record, ...)` for module inputs, with `Path?` for optional file slots
- A single `record(...)` output per module; downstream fields (e.g., `tsv: Path`, `fna: Path`) plus the shared `results` / `logs` / `nf_logs` / `versions` fields

### 4. Validate Channel Patterns
Ensure proper channel patterns:
- Modules: emit a single `record(...)` output whose fields include the downstream payload plus `results`, `logs`, `nf_logs`, `versions`
- Subworkflows: emit `sample_outputs` (per-sample record passthrough) and `run_outputs` (aggregated record); both are `Channel<Record>`
- Workflows: branch outputs by scope (run/sample) before handing them to publishers

## Common Problems and Solutions

### Modules

#### Module Not Finding Input Files
**Symptoms**: Process fails with file not found error
**Solution**:
- Check input channel name matches
- Verify input tuple structure
- Ensure files are staged correctly

#### Module Producing No Output
**Symptoms**: Process completes but no output files
**Solution**:
- Check tool execution command
- Verify working directory
- Look for error messages in log files
- Check if input data is valid

#### Version File Not Created
**Symptoms**: Missing versions.yml file
**Solution**:
- Ensure heredoc syntax is correct
- Check command substitution
- Verify file permissions

### Subworkflows

#### Missing Standard Channels
**Symptoms**: Workflow fails with missing channel error
**Solution**:
- Ensure both `sample_outputs` and `run_outputs` are emitted
- Use `gather()` / `gatherCsvtk()` / `gatherFields()` from the `nf-bactopia` plugin to aggregate per-sample records into run-scope inputs (see [Plugin Functions](04-plugin-functions.md) for the full menu)
- Check channel names in emit block

#### Type Mismatch in Aggregation
**Symptoms**: Type error when mixing channels
**Solution**:
- Ensure consistent record types across channels being mixed
- Use proper type annotations on take blocks (`Channel<Record>`)
- Check that `gather()` / `gatherCsvtk()` projections name fields that actually exist on the source record

### Workflows

#### Output Path Issues
**Symptoms**: Files not appearing in expected location
**Solution**:
- Check `meta.output_dir` values
- Verify publish block syntax
- Ensure proper branching by scope

#### Parameter Not Passed Through
**Symptoms**: Module not receiving parameter
**Solution**:
- Check parameter inheritance chain
- Verify configuration includes
- Use `--` prefix for command line

## Container Entrypoint Workarounds

Some bioinformatics containers have custom entrypoints that initialize environments or set variables. Nextflow overrides container entrypoints with `/bin/bash`, which can break tools that depend on this initialization.

### The Problem

When Nextflow runs a container, it sets `--entrypoint /bin/bash` (Docker) or similar for Singularity/Apptainer. Containers that rely on custom entrypoints (e.g., `/usr/local/env-execute`) to:

- Activate conda environments
- Set environment variables
- Configure tool paths

...will fail because their initialization scripts never run.

### The Solution: Manual Activation

Check for and source the container's activation script at the start of your script block:

```bash
# Nextflow changes the container --entrypoint to /bin/bash
# (container default entrypoint: /usr/local/env-execute)
# Check for container variable initialisation script and source it.
if [ -f "/usr/local/env-activate.sh" ]; then
    set +u  # Otherwise, errors out because of various unbound variables
    . "/usr/local/env-activate.sh"
    set -u
fi
```

### Why `set +u` and `set -u`?

- Nextflow scripts run with `set -eu` by default (exit on error, error on unbound variables)
- Container activation scripts often reference unbound variables during initialization
- `set +u` temporarily disables the unbound variable check
- `set -u` re-enables it after sourcing completes

### Modules Using This Pattern

| Module | Container Issue |
|--------|-----------------|
| `busco` | BUSCO container requires environment activation for Augustus config |

### When to Use This Pattern

Use this workaround when:

1. **Tool fails with environment errors**: Missing variables or paths
2. **Container works outside Nextflow**: Same container works with `docker run` but fails in pipeline
3. **Container has custom entrypoint**: Check Dockerfile for non-standard `ENTRYPOINT`

### Debugging Container Issues

To check if a container has a custom entrypoint:

```bash
# Docker
docker inspect <image> | grep -A5 "Entrypoint"

# Singularity - check runscript
singularity inspect --runscript <image.sif>
```

If the entrypoint is not `/bin/bash` or `/bin/sh`, the container may need the activation workaround.

### Additional Container Workarounds

#### Augustus Config Directory (BUSCO-specific)

BUSCO requires a writable Augustus config directory. When running in containers, this directory may be read-only:

```bash
# If the augustus config directory is not writable, copy to writeable area
if [ ! -w "${AUGUSTUS_CONFIG_PATH}" ]; then
    AUG_CONF_DIR=$( mktemp -d -p $PWD )
    cp -r $AUGUSTUS_CONFIG_PATH/* $AUG_CONF_DIR
    export AUGUSTUS_CONFIG_PATH=$AUG_CONF_DIR
    echo "New AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH}"
fi
```

## Performance Issues

### Slow Execution
**Causes and Solutions**:
- **Insufficient resources**: Increase CPU/memory in process config
- **Large datasets**: Consider chunking or parallelization
- **Network latency**: Use local storage for temporary files
- **Container overhead**: Use native execution if possible

### Memory Issues
**Causes and Solutions**:
- **Large file processing**: Increase memory limit
- **Multiple parallel processes**: Limit concurrency
- **Memory leaks**: Restart pipeline, check tool versions

## Getting Help

### Information to Collect
When seeking help, gather:
1. Full error message with stack trace
2. Nextflow version and configuration
3. Command used to run pipeline
4. Relevant input file information
5. `.nextflow.log` file

### Where to Get Help
- **GitHub Issues**: For bugs and feature requests
- **Discord/Slack**: For community support
- **Documentation**: Check for existing solutions
- **Examples**: See if similar use case exists

## Best Practices to Avoid Issues

### Development
1. **Test incrementally**: Start small, add complexity gradually
2. **Use type annotations**: Catch errors early
3. **Follow patterns**: Use existing templates and conventions
4. **Document deviations**: Explain why patterns were broken

### Execution
1. **Use test data**: Verify with small datasets first
2. **Check resources**: Ensure adequate CPU/memory
3. **Monitor logs**: Watch for warning messages
4. **Validate outputs**: Check if results make sense

### Maintenance
1. **Keep updated**: Use latest Nextflow and tool versions
2. **Regular testing**: Run test suite before changes
3. **Backup configs**: Version control configuration files
4. **Document changes**: Keep changelog updated

## See Also
- [Technical Specifications](../standards/03-technical-specs.md) - For type-related issues
- [Testing Framework](../project/04-testing-framework.md) - For writing tests to prevent issues
- [Development Workflow](../project/02-development-workflow.md) - For development best practices
