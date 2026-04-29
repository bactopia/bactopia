# Testing Framework

## Overview
Bactopia uses the [nf-test](https://www.nf-test.com/) framework (≥0.9.5) with Nextflow (≥26) for comprehensive pipeline testing. Tests cover individual modules, subworkflows, and full pipeline runs using snapshot-based assertions. Day-to-day runs go through the [`bactopia-test`](#running-tests) CLI; the `/run-tests` skill wraps it for repeatable invocations and pairs with `/review-tests` for post-run triage.

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
nextflow.enable.types = true
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

> `params.workflow.ext` is a **string** at module-test scope (a single extension for the module's primary output). Workflow-level configs (`workflows/{name}/nextflow.config`) use the **list** form — e.g. `ext = ['fna']` — because workflows aggregate publishing across multiple module outputs. Match the surrounding layer when editing.
> `bactopia_version` is a placeholder — keep it in sync with the repo's `manifest.version` in [nextflow.config](../../../nextflow.config) when it drifts.

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
                        fna: file("${params.test_data_dir}/species/staphylococcus_aureus/uncompressed/GCF_000017085/main/assembler/GCF_000017085.fna")
                    )
                )
                input[1] = file("${params.test_data_dir}/datasets/mlst/mlst.tar.gz")
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
                        fna: file("${params.test_data_dir}/species/staphylococcus_aureus/compressed/GCF_000017085/main/assembler/GCF_000017085.fna.gz")
                    )
                )
                input[1] = file("${params.test_data_dir}/datasets/mlst/mlst.tar.gz")
                """
            }
        }

        then {
            def sample = workflow.out.sample_outputs[0]
            def run = workflow.out.run_outputs[0]
            assertAll(
                { assert workflow.success },
                { assert workflow.out.sample_outputs != null },
                { assert workflow.out.run_outputs != null },
                { assert snapshot(
                    sample.meta,
                    sample.tsv,
                    sample.versions,
                    run.meta,
                    run.versions
                ).match() },
                { assert sample.results != null },
                { assert run.csv != null },
                { assert run.results != null }
            )
        }
    }
}
```

**Key patterns**:
- Subworkflow tests access both `workflow.out.sample_outputs` (per-sample record passthrough) and `workflow.out.run_outputs` (aggregated run-level record)
- Snapshot matching spans both emits; non-snapshotted assertions (`!= null`) verify the aggregated convenience fields like `results` and `csv`

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
            def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: [])
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

Test data is stored externally and referenced via the `BACTOPIA_TESTS` environment variable, which points at the `data/` directory of the bactopia-tests repo. All in-tree paths shown below are relative to `${params.test_data_dir}` (i.e. `$BACTOPIA_TESTS`) — no leading `data/` prefix inside test files.

- **CI/CD path**: `/data/storage/bactopia-ci/bactopia-tests/data`
- **Repository**: [bactopia/bactopia-tests](https://github.com/bactopia/bactopia-tests)

### Directory Layout

Assemblies sit under `uncompressed/` or `compressed/` subtrees keyed by accession; reads sit under `reads/`; tool databases sit under `datasets/`:

```text
$BACTOPIA_TESTS/
├── species/
│   └── <species_name>/
│       ├── uncompressed/
│       │   └── <accession>/main/assembler/
│       │       └── <accession>.fna
│       ├── compressed/
│       │   └── <accession>/main/assembler/
│       │       └── <accession>.fna.gz
│       └── reads/
│           ├── illumina/
│           │   ├── <accession>_R1.fastq.gz
│           │   └── <accession>_R2.fastq.gz
│           └── <accession>/main/qc/
│               ├── <accession>_R1.fastq.gz
│               └── <accession>_R2.fastq.gz
└── datasets/
    ├── mlst/mlst.tar.gz
    ├── bakta/light/bakta-light.tar.gz
    ├── amrfinderplus/amrfinderplus.tar.gz
    ├── kraken2/k2_viral_20251015.tar.gz
    ├── generic/                    # shared reference assemblies & GFFs
    └── {other tool databases}
```

**Concrete examples** (from the live tests):
- `species/staphylococcus_aureus/uncompressed/GCF_000017085/main/assembler/GCF_000017085.fna`
- `species/staphylococcus_aureus/compressed/GCF_000017085/main/assembler/GCF_000017085.fna.gz`
- `species/portiera/reads/illumina/SRR2838702_R1.fastq.gz`
- `species/streptococcus_pneumoniae/reads/ERR1438863/main/qc/ERR1438863_R1.fastq.gz`
- `datasets/mlst/mlst.tar.gz`

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
            "nf-test": "<nf-test-version>",
            "nextflow": "<nextflow-version>"
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

Use `bactopia-test --generate` (see [Running Tests](#running-tests)). Unlike `nf-test --update-snapshot`, which overwrites in-place, `--generate` deletes the existing snapshot then runs the test twice — the first pass writes a fresh snapshot, the second pass verifies it reproduces.

```bash
bactopia-test \
    --bactopia-path /path/to/bactopia \
    --test-data $BACTOPIA_TESTS \
    --tier modules --include mlst \
    --generate
```

Raw `nf-test --update-snapshot` still works for debugging a single file in isolation:

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

### Recommended: `bactopia-test`

`bactopia-test` is the canonical entry point. It discovers `*.nf.test` files across the tier you select, runs them with consistent env/cache wiring, classifies outcomes into a summary table, and writes per-test logs into `logs/{timestamp}/` for `/review-tests` to interpret.

```bash
# Run a single module
bactopia-test \
    --bactopia-path /path/to/bactopia \
    --test-data $BACTOPIA_TESTS \
    --tier modules --include mlst \
    --profile docker

# Run every subworkflow except two
bactopia-test \
    --bactopia-path /path/to/bactopia \
    --test-data $BACTOPIA_TESTS \
    --tier subworkflows --exclude checkm,ismapper

# Regenerate snapshots for a workflow
bactopia-test \
    --bactopia-path /path/to/bactopia \
    --test-data $BACTOPIA_TESTS \
    --tier workflows --include abricate \
    --generate
```

The project-local `/run-tests` skill wraps these invocations (and enforces guardrails like never passing `--generate` without explicit opt-in). Prefer the skill for interactive sessions; call the CLI directly for scripting.

### `bactopia-test` flag reference

Sourced from `bactopia-test --help`:

| Flag | Group | Purpose |
|---|---|---|
| `--bactopia-path PATH` | Required | Path to the bactopia repo |
| `--test-data PATH` | Required | Test-data directory (sets `BACTOPIA_TESTS`); skip only with `--cleanup` |
| `--tier {modules,subworkflows,workflows,all}` | Selection | Which tier to run |
| `--include NAMES` | Selection | Comma-separated component names to run |
| `--exclude NAMES` | Selection | Comma-separated component names to skip |
| `--profile {docker,singularity,conda}` | Execution | Nextflow profile |
| `--condadir PATH` | Execution | Conda env cache (`NXF_CONDA_CACHEDIR` takes precedence) |
| `--singularity_cache PATH` | Execution | Singularity image cache (`NXF_SINGULARITY_CACHEDIR` takes precedence) |
| `--generate` | Execution | Delete snapshots and run twice to verify reproducibility |
| `--jobs N` | Execution | Parallel workers |
| `--fail-fast` | Execution | Stop on first failure |
| `--timeout MINUTES` | Execution | Per-test timeout (0 disables) |
| `--outdir PATH` | Output | Where `logs/` is written |
| `--json` | Output | Emit structured JSON summary |
| `--keep` | Cleanup | Preserve `.nf-test/` dirs and logs on pass |
| `--cleanup` | Cleanup | Remove all `.nf-test/` temp dirs and exit (no tests run) |

### Raw `nf-test` (debugging)

For ad-hoc debugging of a single file, `nf-test` can still be driven directly:

```bash
cd modules/mlst/tests
nf-test run main.nf.test                    # single module test
nf-test run main.nf.test --profile docker   # specific profile
nf-test run main.nf.test --update-snapshot  # in-place snapshot refresh
nf-test run main.nf.test --debug            # verbose debug output
```

Prefer `bactopia-test` for anything touching more than one file or heading for CI-style verification.

## CI/CD Integration

- **Workflow**: `.github/workflows/all-bactopia-tests.yml`
- **Trigger**: `workflow_dispatch` (manual) on a self-hosted runner
- **Profiles tested**: Singularity, Docker, Conda
- **Parallelism**: Up to 20 parallel subworkflow groups, each with 5 workers

> The current CI workflow still drives the legacy `pytest-workflow` suite (`pytest --wt 5 --git-aware --tag <subworkflow>`), not nf-test. Migrating CI to `bactopia-test` is in progress — until that lands, nf-test coverage is exercised locally (via `bactopia-test` / the `/run-tests` skill) rather than on every dispatch.

### Environment Variables

Nextflow consumes the `NXF_*` cache variables directly; these are what nf-test and `bactopia-test` expect:

| Variable | CI value | Consumed by |
|----------|----------|-------------|
| `BACTOPIA_TESTS` | `/data/storage/bactopia-ci/bactopia-tests/data` | Test files via `params.test_data_dir` |
| `NXF_CONDA_CACHEDIR` | `/data/storage/bactopia-ci/envs/conda` | Nextflow (Conda profile) |
| `NXF_SINGULARITY_CACHEDIR` | `/data/storage/bactopia-ci/envs/singularity` | Nextflow (Singularity profile) |

> The CI workflow also exports `BACTOPIA_CONDA` and `BACTOPIA_SINGULARITY` as script-internal shadows (passed into `bactopia-test` via `--condadir` / `--singularity_cache`), but these are workflow wiring — Nextflow itself reads the `NXF_*` forms. When setting caches outside CI, use the `NXF_*` variables.

## Best Practices

1. **Test both compressed and uncompressed inputs** when a module supports both
2. **Use null** for optional `Path?` inputs that aren't needed in a test case
3. **Snapshot tool-specific fields** -- always include `record.meta`, tool outputs, and `record.versions`
4. **Keep test data small** -- use minimal genome assemblies and datasets
5. **Tag tests** with component type (`modules`, `subworkflows`, `workflows`) and tool name

## See Also
- [Development Workflow](02-development-workflow.md) - For creating tests with new tools
- [Repository Structure](01-repository-structure.md) - For test file locations
