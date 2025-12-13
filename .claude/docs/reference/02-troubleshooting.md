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
- Example fix: `tuple(meta, files("*.txt"))` instead of `tuple(meta, file("*.txt"))`

#### "Missing required parameter"
**Cause**: Parameter not defined in configuration or schema
**Solution**:
- Check `schema.json` and `params.config` files
- Ensure parameter is defined with correct type
- Verify parameter is being passed correctly

#### "Channel not found"
**Cause**: Workflow trying to access a channel that doesn't exist
**Solution**:
- Verify all 4 channels are emitted from subworkflows
- Check channel names match exactly
- Ensure proper channel mixing/branching

### Path and File Issues

#### "No such file or directory"
**Cause**: File path doesn't exist or is incorrect
**Solution**:
- Check if EMPTY_* files are being used correctly for optional inputs
- Verify relative vs absolute paths
- Ensure input files are properly staged

#### "EMPTY_PROTEINS" in output
**Cause**: Optional parameter not provided, workaround file name showing
**Solution**:
- This is expected behavior for optional inputs
- Check parameter name and provide file if needed
- See [Technical Specifications](../standards/03-technical-specs.md) for Path? workarounds

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
- Verify Nextflow version meets requirements (v25.10.x+)
- Check that all required files exist
- Ensure proper permissions on input/output directories

### 2. Look for TODO Comments
TODO comments in the code indicate known limitations or temporary workarounds
```groovy
// TODO: Remove when Path? is fixed
```

### 3. Verify Meta Map Fields
Ensure meta map contains required fields:
```groovy
meta.id
meta.name
meta.scope  // 'sample' or 'run'
meta.output_dir
meta.logs_dir
meta.process_name
```

### 4. Check Consistent Typing
Verify consistent typing across connected components:
- `Tuple<Map, Set<Path>>` for module inputs
- `Tuple<Map, Path>` for single file outputs
- `Tuple<Map, Set<Path>>` for multiple file outputs

### 5. Validate Channel Patterns
Ensure proper channel patterns:
- Modules: Use specific output channels
- Subworkflows: Always emit 4 standard channels
- Workflows: Branch outputs by scope (run/sample)

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
- Ensure all 4 channels are emitted
- Use `flattenPaths` for aggregate outputs
- Check channel names in emit block

#### Type Mismatch in Aggregation
**Symptoms**: Type error when mixing channels
**Solution**:
- Ensure consistent types across channels
- Use proper type annotations
- Check `flattenPaths` usage

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