# Plugin Functions

## Overview

The `nf-bactopia` plugin provides utility functions for channel manipulation that are used throughout Bactopia subworkflows. The primary user surface is a set of channel-manipulation helpers — `gather`, `gatherCsvtk`, `gatherFields`, `filterWithData`, `combineWith`, and `collectNextflowLogs` — plus a handful of workflow-initialization helpers (`bactopiaInputs`, `bactopiaToolInputs`, `validateParameters`) that `utils_bactopia*` subworkflows wire in.

All functions are exposed as Nextflow `@Function` extensions and are imported individually:

```groovy
include { gather         } from 'plugin/nf-bactopia'
include { gatherCsvtk    } from 'plugin/nf-bactopia'
include { gatherFields   } from 'plugin/nf-bactopia'
include { filterWithData } from 'plugin/nf-bactopia'
include { combineWith    } from 'plugin/nf-bactopia'
include { collectNextflowLogs } from 'plugin/nf-bactopia'
```

Implementations live in [nf-bactopia/src/main/groovy/bactopia/plugin/BactopiaExtension.groovy](https://github.com/bactopia/nf-bactopia/blob/main/src/main/groovy/bactopia/plugin/BactopiaExtension.groovy) and [ChannelUtils.groovy](https://github.com/bactopia/nf-bactopia/blob/main/src/main/groovy/bactopia/plugin/utils/ChannelUtils.groovy).

## Channel Manipulation Functions

### gather

Extracts one field from each record in a channel and collects the values into a single `Set`, wrapped as a record-like map keyed by `meta`.

#### Signature

```groovy
gather(chResults, field, meta)
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chResults` | Channel or List of records | Yes | Direct-assigned result from a module call (e.g., `ch_tool = TOOL_MODULE(...)`) |
| `field` | String | Yes | Record field name to extract (positional, not named) |
| `meta` | Map | Yes | Output meta map; must contain `name`. All other keys pass through |

#### Returns

A single-element channel (or List) emitting one record-like map:

```
[meta: [name: 'tool', ...], <field>: Set<values>]
```

The original per-sample meta maps are discarded; only the supplied `meta` is retained. The output field name matches the input `field` (to rename for CSVTK_CONCAT, use `gatherCsvtk` instead).

#### Usage Example

```groovy
// Collect all per-sample JSON outputs for a downstream heatmap step
ch_rgi_heatmap = RGI_HEATMAP(gather(ch_rgi_main, 'json', [name: 'rgi']))
```

Reference: [subworkflows/rgi/main.nf:41](../../../subworkflows/rgi/main.nf#L41)

---

### gatherCsvtk

Identical to `gather`, but renames the extracted field to `csv` in the output record — matching the input field name that `CSVTK_CONCAT` expects. This is the dominant aggregation idiom across bacterial-typing subworkflows.

#### Signature

```groovy
gatherCsvtk(chResults, field, meta)
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chResults` | Channel or List of records | Yes | Direct-assigned module result |
| `field` | String | Yes | Record field name to extract (typically `'tsv'` or `'report'`) |
| `meta` | Map | Yes | Output meta map; must contain `name` |

#### Returns

```
[meta: [name: 'tool', ...], csv: Set<values>]
```

#### Usage Example

```groovy
ch_sistr = SISTR_MODULE(assembly)
ch_concat = CSVTK_CONCAT(gatherCsvtk(ch_sistr, 'tsv', [name: 'sistr']), 'tsv', 'tsv')
```

Reference: [subworkflows/sistr/main.nf:41-42](../../../subworkflows/sistr/main.nf#L41)

---

### gatherFields

Gathers multiple fields from records with an explicit input-to-output rename map. Used when a downstream process needs several aggregated fields in one record.

#### Signature

```groovy
gatherFields(chResults, fieldMapping, meta)
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chResults` | Channel or List of records | Yes | Direct-assigned module result |
| `fieldMapping` | `Map<String, String>` | Yes | `[inputField: outputField]` — extracts `inputField` from each record, emits as `outputField` |
| `meta` | Map | Yes | Output meta map; must contain `name` |

#### Returns

```
[meta: [name: 'tool', ...], <outputField1>: Set<values>, <outputField2>: Set<values>, ...]
```

#### Usage Examples

Single-field rename (pairs a query channel for `combineWith`):

```groovy
gatherFields(query, [fna: 'query'], [name: 'fastani'])
```

Reference: [subworkflows/fastani/main.nf:47](../../../subworkflows/fastani/main.nf#L47)

Multi-field rename (collects VCFs and aligned FASTAs for core-SNP analysis):

```groovy
ch_core_input = gatherFields(
    ch_snippy.variants,
    [vcf: '_vcf', aligned_fa: '_aligned_fa'],
    [name: 'core-snp']
)
```

Reference: [workflows/bactopia-tools/snippy/main.nf:104-108](../../../workflows/bactopia-tools/snippy/main.nf#L104)

---

### filterWithData

Filters records where at least one of the specified fields is non-null (and non-empty for collections), then projects each surviving record down to just `meta` plus the requested fields. Prevents downstream processes from receiving records with extra fields that would trip type-checking.

#### Signature

```groovy
filterWithData(input, fields)
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `input` | Channel or List of records | Yes | Record channel to filter |
| `fields` | `List<String>` | Yes | Field names to check for presence AND include in the output record |

#### Returns

A channel (or List) of projected records, each containing only `meta` plus the listed `fields`. Records where all listed fields are null are dropped.

#### Usage Example

Drop samples that have no reads, and project down to just the read fields:

```groovy
ch_ariba_run = ARIBA_RUN(filterWithData(reads, ['r1', 'r2']), ch_ariba_getref.map { r -> r.db })
```

Reference: [subworkflows/ariba/main.nf:53](../../../subworkflows/ariba/main.nf#L53)

---

### combineWith

Creates a cartesian product by combining a gathered channel (single-item, typically from `gatherFields`) with a multi-item channel, merging each item into the gathered map under a specified field name. Replaces the deprecated Nextflow `each` input qualifier.

#### Signature

```groovy
combineWith(gathered, items, field)
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gathered` | Channel or List | Yes | Single-item gathered channel (e.g., from `gatherFields`) |
| `items` | Channel or List | Yes | Multi-item channel (e.g., individual reference paths) |
| `field` | String | Yes | Field name under which each item is merged into the gathered map |

#### Returns

A channel (or List) of records, one per `item`, each with the item merged into the gathered map under `field`.

If `gathered` emits `[meta: [name: 'fastani'], query: [a.fna, b.fna]]` and `items` emits `ref1.fna` and `ref2.fna`, the result emits:

- `[meta: [name: 'fastani'], query: [a.fna, b.fna], reference: ref1.fna]`
- `[meta: [name: 'fastani'], query: [a.fna, b.fna], reference: ref2.fna]`

#### Usage Example

Pairs a query channel with each reference for FastANI pairwise mode (replacing the deprecated `each` qualifier):

```groovy
ch_ref = reference.map { r -> r.fna }
ch_fastani = FASTANI_MODULE(
    combineWith(
        gatherFields(query, [fna: 'query'], [name: 'fastani']),
        ch_ref,
        'reference'
    )
)
```

Reference: [subworkflows/fastani/main.nf:47](../../../subworkflows/fastani/main.nf#L47)

---

### collectNextflowLogs

Flat-maps each record's `nf_logs` field into individual `[meta, file]` tuples suitable for the workflow-tier `nf_logs` publishing pattern.

#### Signature

```groovy
collectNextflowLogs(chResults)
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chResults` | Channel or List of records | Yes | Records containing `meta` and `nf_logs` fields |

#### Returns

A channel (or List) of `[meta, file]` tuples — one per log file per record.

#### Usage Example

```groovy
include { collectNextflowLogs } from 'plugin/nf-bactopia'

// In workflow-tier main.nf, aggregate nf_logs for publishing
ch_nf_logs = collectNextflowLogs(ch_sample_outputs)
```

References: [main.nf:182](../../../main.nf#L182), [workflows/teton/main.nf:69](../../../workflows/teton/main.nf#L69)

---

## Workflow Initialization Functions

These are wired into the `utils_bactopia` and `utils_bactopia-tools` subworkflows — user workflows generally don't call them directly.

### bactopiaInputs

```groovy
bactopiaInputs(runType)
```

Collects sample inputs for standard Bactopia runs from the configured sample sheet. Returns a Map with `hasErrors`, `error`, `logs`, and `samples` keys.

### bactopiaToolInputs

```groovy
bactopiaToolInputs()
```

Same as `bactopiaInputs` but for Bactopia Tools workflows. Returns the same Map shape.

### validateParameters

```groovy
validateParameters(isBactopiaTool)
```

Validates pipeline parameters against the schema before the workflow body runs. Returns a Map of validation results and captured logs. Invoked during workflow initialization, not in user code.

---

## Channel Flow Diagram

```text
Subworkflow
───────────

ch_tool = TOOL_MODULE(assembly, db)
    │
    ├─→ emit: sample_outputs = ch_tool.sample_outputs
    │
    └─→ ch_concat = CSVTK_CONCAT(
            gatherCsvtk(ch_tool, 'tsv', [name: 'tool']),
            'tsv', 'tsv'
        )
            │
            └─→ emit: run_outputs = ch_concat.run_outputs
```

## Common Patterns

### Subworkflow Standard Pattern

```groovy
nextflow.enable.types = true

include { TOOL_MODULE  } from '../../modules/tool/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { gatherCsvtk  } from 'plugin/nf-bactopia'

workflow TOOL {
    take:
    assembly: Channel<Record>
    db: Path

    main:
    ch_tool = TOOL_MODULE(assembly, db)
    ch_concat = CSVTK_CONCAT(gatherCsvtk(ch_tool, 'tsv', [name: 'tool']), 'tsv', 'tsv')

    emit:
    sample_outputs = ch_tool.sample_outputs
    run_outputs = ch_concat.run_outputs
}
```

### Aggregation for Summary

```groovy
// Collect per-sample TSVs from the module and concatenate into a single file
ch_concat = CSVTK_CONCAT(
    gatherCsvtk(ch_mlst, 'tsv', [name: 'mlst']),
    'tsv',
    'tsv'
)
```

## See Also

- [Technical Specifications](../standards/03-technical-specs.md) — channel type conventions
- [Examples](./01-examples.md) — full subworkflow examples in context
- [Subworkflow Documentation](../standards/04-subworkflow-documentation.md) — subworkflow patterns
