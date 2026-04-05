# Testing Framework

## Overview
Bactopia uses the [nf-test](https://www.nf-test.com/) framework (v0.9.4) with Nextflow 26.01.0 for comprehensive pipeline testing. Tests cover individual modules, subworkflows, and full pipeline runs using snapshot-based assertions.

## Test Organization

### Directory Structure

Tests live alongside their components:

```text
modules/{tool}/tests/
    main.nf.test          # Test specification
    main.nf.test.snap     # Snapshot file (expected outputs)
    nf-test.config        # nf-test configuration
    nextflow.config       # Nextflow configuration for test
    .nftignore            # File patterns to exclude from snapshots

subworkflows/{name}/tests/
    main.nf.test
    main.nf.test.snap
    nf-test.config
    nextflow.config
    .nftignore

workflows/{name}/tests/
    main.nf.test
    main.nf.test.snap
    nf-test.config
    .nftignore

tests/                    # Top-level pipeline tests
    main.nf.test
    main.nf.test.snap
    nf-test.config
    .nftignore
```

### Test Block Types

| Block Type | Scope | Used For |
|------------|-------|----------|
| `nextflow_process` | Module | Testing individual processes |
| `nextflow_workflow` | Subworkflow | Testing workflow orchestration |
| `nextflow_pipeline` | Pipeline | Testing full end-to-end runs |

## Test Configuration

### nf-test.config (Module level)

```groovy
config {
    testsDir "."
    workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"
    configFile "nextflow.config"
    profile "docker"
    options "--is_ci --max_memory 8.GB"

    plugins {
        load "nft-utils@0.0.5"
    }
}
```

### nextflow.config (Module level)

```groovy
// Minimal config for module-level testing
nextflow.preview.types = true
nextflow.enable.strict = true

params {
    workflow {
        name = "tool_name"
        logo_name = "bactopia-tools"
        description = "Tool description"
        ext = "fna"
    }

    bactopia_version = '4.0.0'
    bactopia_cache = System.getenv("BACTOPIA_CACHEDIR") ?: "${System.getenv('HOME')}/.bactopia"
    condadir = "${params.bactopia_cache}/conda"
    wf = params.workflow.name
    merge_folder = "merged-results"
    test_data_dir = System.getenv("BACTOPIA_TESTS") ?: ""
    is_ci = true

    // Max Job Request Parameters
    max_retry = 1
    max_time = 2.h
    max_memory = 8.GB
    max_cpus = 2

    // Nextflow Profile Parameters
    registry = "quay.io"
    singularity_cache = "${params.bactopia_cache}/singularity"
    singularity_pull_docker_container = false
    container_opts = ""
}

includeConfig "../module.config"
includeConfig "../../../conf/base.config"
includeConfig "../../../conf/profiles.config"
```

## Writing Tests

### Module Test (nextflow_process)

```groovy
nextflow_process {
    name "Test MLST"
    script "../main.nf"
    process "MLST"
    tag "modules"
    tag "mlst"

    test("mlst - module - GCF_000017085") {
        when {
            params {
                test_data_dir = System.getenv("BACTOPIA_TESTS") ?: ""
            }
            process {
                """
                input[0] = Channel.of(
                    record(
                        meta: [name: "GCF_000017085"],
                        assembly: file("${params.test_data_dir}/data/species/staphylococcus_aureus/genome/GCF_000017085.fna")
                    )
                )
                input[1] = file("${params.test_data_dir}/data/datasets/mlst/mlst.tar.gz")
                """
            }
        }

        then {
            def record = process.out[0][0]
            assertAll(
                { assert process.success },
                { assert snapshot(
                    record.meta,
                    record.tsv,
                    record.versions
                ).match() }
            )
        }
    }
}
```

**Key patterns**:
- Input uses `record()` syntax matching the module's input declaration
- `input[0]` for the primary record input, `input[1]` for additional inputs (db, etc.)
- Output accessed via `process.out[0][0]` to get the first record
- Record fields accessed via `record.fieldName`
- Snapshot matching captures meta, tool-specific outputs, and versions

### Subworkflow Test (nextflow_workflow)

```groovy
nextflow_workflow {
    name "Test MLST"
    script "../main.nf"
    workflow "MLST"
    tag "subworkflows"
    tag "mlst"

    test("mlst - subworkflow - GCF_000017085") {
        when {
            params {
                test_data_dir = System.getenv("BACTOPIA_TESTS") ?: ""
            }
            workflow {
                """
                input[0] = Channel.of(
                    record(
                        meta: [name: "GCF_000017085"],
                        assembly: file("${params.test_data_dir}/data/species/staphylococcus_aureus/genome/GCF_000017085.fna")
                    )
                )
                input[1] = file("${params.test_data_dir}/data/datasets/mlst/mlst.tar.gz")
                """
            }
        }

        then {
            def outputs = workflow.out.sample_outputs[0]
            assertAll(
                { assert workflow.success },
                { assert snapshot(...).match() }
            )
        }
    }
}
```

### Pipeline Test (nextflow_pipeline)

```groovy
nextflow_pipeline {
    name "Test Pipeline"
    script "../main.nf"
    tag "workflows"
    tag "pipeline_name"

    test("pipeline - description") {
        when {
            params {
                outdir = "$outputDir"
                // test parameters
            }
        }

        then {
            def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true)
            def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: '.nftignore')
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    workflow.trace.succeeded().size(),
                    stable_name,
                    stable_path
                ).match() }
            )
        }
    }
}
```

**Key patterns**:
- Uses `getAllFilesFromDir` from the `nft-utils` plugin
- `ignoreFile: '.nftignore'` excludes unstable files from content checks
- Captures task count, file structure, and file content in snapshot

## Test Data

### Location

Test data is stored externally and referenced via the `BACTOPIA_TESTS` environment variable:
- **CI/CD path**: `/data/storage/bactopia-ci/bactopia-tests/data`
- **Repository**: [bactopia/bactopia-tests](https://github.com/bactopia/bactopia-tests)

### Directory Layout

```text
$BACTOPIA_TESTS/data/
├── species/
│   ├── staphylococcus_aureus/
│   │   └── genome/
│   │       └── GCF_000017085.fna
│   ├── neisseria_gonorrhoeae/
│   │   └── genome/
│   │       ├── GCF_001047255.fna
│   │       └── GCF_001047255.fna.gz
│   ├── haemophilus_influenzae/
│   │   └── genome/
│   │       ├── GCF_900478275.fna
│   │       └── GCF_900478275.fna.gz
│   └── portiera/
│       ├── compressed/
│       │   └── SRR2838702/main/qc/
│       │       ├── SRR2838702_R1.fastq.gz
│       │       └── SRR2838702_R2.fastq.gz
│       ├── genome/
│       │   └── GCF_000292685.fna.gz
│       └── illumina/
│           ├── SRR2838702.fna
│           ├── SRR2838702.faa
│           └── SRR2838702.gff
├── datasets/
│   ├── mlst/
│   │   └── mlst.tar.gz
│   ├── bakta/
│   │   └── light/
│   │       └── bakta-light.tar.gz
│   ├── amrfinderplus/
│   │   └── amrfinderplus.tar.gz
│   └── {other tool databases}
```

### Optional Inputs

Optional inputs use `Path?` types and are passed as `null` when not needed in a test:

```groovy
input[2] = null
```

## Snapshot Files

### Format

Snapshots are JSON files (`.nf.test.snap`) containing expected outputs with MD5 checksums:

```json
{
    "mlst - module - GCF_000017085": {
        "content": [
            {
                "id": "GCF_000017085-MLST",
                "logs_dir": "GCF_000017085/tools/mlst//logs/",
                "name": "GCF_000017085",
                "output_dir": "GCF_000017085/tools/mlst/",
                "process_name": "mlst",
                "scope": "sample"
            },
            "GCF_000017085.tsv:md5,67e1b29068b46c7ccd845ec36d41da89",
            [
                "versions.yml:md5,2d78d6f807230df01693f94e1d484ca5"
            ]
        ],
        "timestamp": "2026-03-11T05:58:00.940736931",
        "meta": {
            "nf-test": "0.9.4",
            "nextflow": "26.01.0"
        }
    }
}
```

**Structure**:
1. First element: metadata map (meta fields from the record)
2. Subsequent elements: file checksums (`filename:md5,hash`)
3. Arrays indicate `Set<Path>` fields (like `versions`)
4. Metadata section tracks nf-test and Nextflow versions

### Updating Snapshots

```bash
nf-test run path/to/tests/ --update-snapshot
```

## .nftignore Files

Exclude unstable or large files from snapshot content comparisons:

```text
**/*.{err,gz,html,log,pdf,stderr,stdout}
**/nf.command.*
bactopia-runs/**/nf-reports/*.{dot,html}

**/*.{corrections,hist,histogram,json,msh,sig,txt,tsv,zip}
```

Files matching these patterns are excluded from `getAllFilesFromDir` content hashing when using `ignoreFile: '.nftignore'`. File names are still tracked, only content comparison is skipped.

## Running Tests

### Individual Module Test

```bash
cd modules/mlst/tests
nf-test run main.nf.test
```

### With Specific Profile

```bash
nf-test run main.nf.test --profile docker
```

### Update Snapshots

```bash
nf-test run main.nf.test --update-snapshot
```

### Debug Mode

```bash
nf-test run main.nf.test --debug
```

## CI/CD Integration

Tests run via GitHub Actions on a self-hosted runner:

- **Workflow**: `.github/workflows/all-bactopia-tests.yml`
- **Trigger**: `workflow_dispatch` (manual)
- **Profiles tested**: Singularity, Docker, Conda
- **Parallelism**: Up to 20 parallel test groups, each with 5 workers

### Environment Variables

| Variable | Value |
|----------|-------|
| `BACTOPIA_TESTS` | `/data/storage/bactopia-ci/bactopia-tests/data` |
| `BACTOPIA_CONDA` | `/data/storage/bactopia-ci/envs/conda` |
| `BACTOPIA_SINGULARITY` | `/data/storage/bactopia-ci/envs/singularity` |

## Best Practices

1. **Test both compressed and uncompressed inputs** when a module supports both
2. **Use null** for optional `Path?` inputs that aren't needed in a test case
3. **Snapshot tool-specific fields** -- always include `record.meta`, tool outputs, and `record.versions`
4. **Keep test data small** -- use minimal genome assemblies and datasets
5. **Tag tests** with component type (`modules`, `subworkflows`, `workflows`) and tool name

## See Also
- [Development Workflow](02-development-workflow.md) - For creating tests with new tools
- [Repository Structure](01-repository-structure.md) - For test file locations
