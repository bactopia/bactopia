# Plugin Functions

## Overview

The `nf-bactopia` plugin provides utility functions for channel manipulation that are used throughout Bactopia subworkflows. The primary function is `gather`, which extracts fields from module record outputs and aggregates them for downstream processing.

## Import Statement

```groovy
include { gather } from 'plugin/nf-bactopia'
```

## Functions

### gather

Extracts a specific field from module record outputs and aggregates them into a single collection for downstream processing (e.g., concatenation with CSVTK_CONCAT).

#### Signature

```groovy
gather(moduleOut, toolName, field: 'fieldName')
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `moduleOut` | Channel | Yes | The module's full record output channel (e.g., `MODULE.out`) |
| `toolName` | String | Yes | Tool identifier used as the `id` in output meta |
| `field` | String (named) | Yes | The record field name to extract from each record |

#### Returns

Single tuple: `[[id: toolName], Set<extracted_field_values>]`

- Extracts the specified field from each record in the channel
- The original meta maps are discarded
- All extracted field values are collected into a `Set`
- New meta contains only `id`

#### Usage Example

```groovy
// Aggregate TSV outputs for concatenation
CSVTK_CONCAT(
    gather(MLST_MODULE.out, 'mlst', field: 'tsv'),
    'tsv',
    'tsv'
)

// Aggregate report outputs
CSVTK_CONCAT(
    gather(ABRICATE_RUN.out, 'abricate', field: 'report'),
    'tsv',
    'tsv'
)
```

#### How It Works

1. Takes the full record output channel from a module
2. Extracts the value of the named `field:` from each record
3. Collects all extracted values into a set
4. Wraps them with a new meta map containing only the tool name

---

### combineWith

Creates a cartesian product by combining a gathered channel (typically single-item, from `gatherFields`) with a multi-item channel, merging each item into the gathered map under a specified field name. Replaces the deprecated Nextflow `each` input qualifier.

#### Signature

```groovy
combineWith(gathered, items, field: 'fieldName')
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gathered` | Channel or List | Yes | A gathered channel (e.g. from `gatherFields`) |
| `items` | Channel or List | Yes | Multi-item channel (e.g. individual reference paths) |
| `field` | String | Yes | Field name to assign each item in the output map |

#### Returns

Channel or List of maps, one per item, with the item merged into the gathered map under the field name.

#### Usage Example

```groovy
// Replace deprecated `each` qualifier for FastANI pairwise mode
ch_ref = reference.map { r -> r.fna }
FASTANI_MODULE(
    combineWith(
        gatherFields(query, [fna: 'query'], [name: 'fastani']),
        ch_ref,
        'reference'
    )
)
```

#### How It Works

1. Takes a gathered channel (e.g. `[meta:[name:fastani], query:[a.fna, b.fna]]`) and a multi-item channel (e.g. `ref1.fna`, `ref2.fna`)
2. Uses `combine` to create the cartesian product: `[[meta:..., query:...], ref1.fna]`, `[[meta:..., query:...], ref2.fna]`
3. Merges each item into the gathered map under the field name: `[meta:..., query:..., reference:ref1.fna]`, `[meta:..., query:..., reference:ref2.fna]`
4. The resulting maps match process Record inputs with named fields

---

### flattenPaths (deprecated)

`flattenPaths` was previously used in subworkflows to convert `Tuple<Map, Set<Path>>` channels to `Tuple<Map, Path>` for the old 4-channel output pattern. It is **no longer used** in subworkflows.

Subworkflows now pass through module record outputs directly:

```groovy
emit:
sample_outputs = MODULE.out
run_outputs = CSVTK_CONCAT.out
```

The function still exists in the plugin for backwards compatibility but should not be used in new code.

---

### validateParameters

Validates pipeline parameters before execution.

#### Purpose

- Checks that required inputs are present
- Validates parameter dependencies (e.g., if using X, must provide Y)
- Ensures specified files exist on the filesystem

#### Usage

Called automatically during workflow initialization. Not typically invoked directly in user code.

---

### bactopiaInputs

Initializes input channels from sample sheets and input specifications.

**Status**: Pending refactoring - documentation will be added after updates are complete.

---

## Channel Flow Diagram

```text
Module Record Outputs              Subworkflow
─────────────────────              ───────────

MODULE(input, db)
    │
    ├─→ MODULE.out ──────────────→ emit: sample_outputs
    │       (full record with meta, tsv, results, logs, etc.)
    │
    └─→ gather(MODULE.out, 'tool', field: 'tsv')
            │
            └─→ CSVTK_CONCAT ──→ emit: run_outputs
                    (aggregated record)
```

## Common Patterns

### Subworkflow Standard Pattern

```groovy
nextflow.preview.types = true

include { TOOL_MODULE } from '../../modules/tool/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { gather } from 'plugin/nf-bactopia'

workflow TOOL {
    take:
    assembly: Channel<Record>
    db: Path

    main:
    TOOL_MODULE(assembly, db)
    CSVTK_CONCAT(gather(TOOL_MODULE.out, 'tool', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = TOOL_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
```

### Aggregation for Summary

```groovy
// Collect all TSV outputs from the module and concatenate into single file
CSVTK_CONCAT(
    gather(MLST_MODULE.out, 'mlst', field: 'tsv'),
    'tsv',
    'tsv'
)
```

## See Also

- [Technical Specifications](../standards/03-technical-specs.md) - For channel type conventions
- [Examples](./01-examples.md) - For usage in context
- [Subworkflow Documentation](../standards/04-subworkflow-documentation.md) - For subworkflow patterns
