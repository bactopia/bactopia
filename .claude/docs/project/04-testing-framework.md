# Testing Framework

## Overview
Bactopia uses the nf-test framework for comprehensive pipeline testing. This ensures reliability and correctness across all components, from individual modules to end-to-end workflows.

## nf-test Framework

### Key Components
- **Test files**: Files ending with `.nftest` or `.test`
- **Test configuration**: `.nf-test.config`
- **Test snapshots**: Files ending with `.test.snap`
- **Test data**: External test data repository

### Test Organization

#### Main Tests
- **`/tests/main.nf.test`** - Main workflow tests
- **`/tests/main.nf.test.snap`** - Expected outputs snapshot

#### Workflow-Specific Tests
- **`{workflow}/tests/`** - Tests specific to individual workflows
- **`{workflow}/tests/*.nftest`** - Test files for the workflow
- **`{workflow}/tests/*.test.snap`** - Workflow-specific snapshots

#### Test Configuration
- **`/tests/nf-test.config`** - Global test configuration
- **`{workflow}/tests/nf-test.config`** - Workflow-specific test config

## Test Structure

### Basic Test Template
```groovy
test("test_name") {
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
        // Additional assertions
    }
}
```

### Test Categories

#### 1. Unit Tests
- **Purpose**: Test individual modules in isolation
- **Scope**: Single process execution
- **Focus**: Input/output validation, parameter handling

#### 2. Integration Tests
- **Purpose**: Test subworkflow orchestration
- **Scope**: Multiple connected processes
- **Focus**: Data flow between components

#### 3. End-to-End Tests
- **Purpose**: Test complete workflows
- **Scope**: Full pipeline execution
- **Focus**: Real-world usage scenarios

#### 4. Regression Tests
- **Purpose**: Validate output consistency
- **Scope**: Compare with expected outputs
- **Focus**: Preventing unintended changes

## Test Data Management

### Test Data Source
- **Repository**: [bactopia-tests](https://github.com/bactopia/bactopia-tests)
- **Environment Variable**: `BACTOPIA_TESTS` must point to cloned repository
- **Structure**: Organized by tool/workflow

### Test Data Types
- **Input files**: Sample data for testing
- **Reference files**: Expected outputs
- **Configuration**: Test-specific configs

### Test Isolation
- Each test runs in a separate work directory
- Temporary files managed by nf-test
- No interference between tests

## Running Tests

### All Tests
```bash
nf-test run
```

### Specific Test
```bash
nf-test run tests/main.nf.test
```

### With Profile
```bash
nf-test run -profile conda
```

### Debug Mode
```bash
nf-test run --debug
```

### Test Report
```bash
nf-test run --report tests/report.html
```

## Test Configuration

### Global Config (`tests/nf-test.config`)
```groovy
testConfig {
    nextflow {
        // Path to Nextflow binary
        executable = 'nextflow'

        // Base configuration
        configFile = 'nextflow.config'

        // Profiles to use
        profiles = ['conda', 'test']
    }

    // Test environment variables
    env {
        BACTOPIA_TESTS = '/path/to/bactopia-tests'
    }
}
```

### Workflow-Specific Config
```groovy
testConfig {
    // Override global config for specific workflow
    workflow {
        name = 'my-tool'
        profiles = ['conda', 'test']
    }
}
```

## Writing Tests

### Testing Parameters
```groovy
test("parameter_validation") {
    when {
        process {
            """
            // Test with specific parameters
            --tool_option value
            --threshold 90
            """
        }
    }
    then {
        assert workflow.success
        // Verify parameter effects
    }
}
```

### Testing Outputs
```groovy
test("output_validation") {
    when {
        // Setup code
    }
    then {
        assert workflow.success
        assert path("${workflow.workDir}/output.txt").exists()
        // Additional output validations
    }
}
```

### Testing Error Conditions
```groovy
test("error_handling") {
    when {
        process {
            """
            // Provide invalid input
            --invalid_option
            """
        }
    }
    then {
        assert !workflow.success
        assert workflow.error != null
    }
}
```

## Best Practices

### Test Design
1. **Test one thing** - Each test should validate one specific behavior
2. **Use descriptive names** - Test names should clearly indicate what's being tested
3. **Include edge cases** - Test boundary conditions and error cases
4. **Make tests independent** - Tests should not depend on each other

### Test Data
1. **Use consistent test data** - Same inputs across tests for comparability
2. **Keep test data small** - Use minimal datasets for faster tests
3. **Document test data** - Include README with test data descriptions

### Assertions
1. **Be specific** - Assert exact conditions, not just success
2. **Include helpful messages** - Explain what failed
3. **Check outputs** - Verify file existence and content

## Common Test Patterns

### Module Testing
```groovy
test("module_basic_functionality") {
    when {
        process {
            """
            // Provide minimal required input
            input_file.txt
            """
        }
    }
    then {
        assert workflow.success
        // Check for expected outputs
    }
}
```

### Workflow Testing
```groovy
test("workflow_complete") {
    when {
        process {
            """
            // Run full workflow with sample data
            --sample_id test_sample
            input_reads/
            """
        }
    }
    then {
        assert workflow.success
        // Validate workflow outputs
    }
}
```

## Troubleshooting Tests

### Common Issues
1. **Missing test data** - Ensure BACTOPIA_TESTS is set correctly
2. **Container issues** - Check container availability
3. **Permission errors** - Verify file permissions
4. **Timeout** - Increase test timeout for long-running tests

### Debug Tips
1. **Use --debug flag** - See detailed execution information
2. **Check work directory** - Inspect intermediate files
3. **Run with -profile test** - Use test-specific configuration
4. **Check snapshots** - Ensure expected outputs are current

## See Also
- [Development Workflow](../project/02-development-workflow.md) - For creating tests with new tools
- [Repository Structure](../project/01-repository-structure.md) - For test file locations