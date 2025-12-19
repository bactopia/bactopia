# Plugin Functions

## Overview

The `nf-bactopia` plugin provides utility functions for channel manipulation that are used throughout Bactopia subworkflows. These functions handle the aggregation and transformation of channels, enabling the standard 4-channel output pattern.

## Import Statement

```groovy
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather } from 'plugin/nf-bactopia'
```

## Functions

### gather

Aggregates outputs from a channel into a single collection for downstream processing (e.g., concatenation with CSVTK_CONCAT).

#### Signature

```groovy
gather(chResults, toolName, args = '')
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chResults` | Channel or List | Yes | Channel or list of `[meta, output]` tuples |
| `toolName` | String | Yes | Tool identifier used as the `id` in output meta |
| `args` | String | No | Optional arguments string (default: `''`) |

#### Returns

Single tuple: `[[id: toolName, args: args], Set<outputs>]`

- The original meta maps are discarded
- All outputs are collected into a `Set`
- New meta contains only `id` and `args`

#### Usage Example

```groovy
// Aggregate TSV outputs for concatenation
CSVTK_CONCAT(
    gather(MLST_MODULE.out.tsv, 'mlst'),
    'tsv',
    'tsv'
)

// With optional args
CSVTK_CONCAT(
    gather(ABRICATE_RUN.out.report, 'abricate', '--header'),
    'tsv',
    'tsv'
)
```

#### Implementation

```groovy
static Object gather(Object chResults, String toolName, String args = '') {
    // Input validation
    if (chResults == null) {
        throw new IllegalArgumentException("chResults cannot be null")
    }
    if (toolName == null || toolName.trim().isEmpty()) {
        throw new IllegalArgumentException("toolName cannot be null or empty")
    }

    // Detect if input is a channel
    if (chResults instanceof DataflowReadChannel || chResults instanceof DataflowWriteChannel) {
        // Apply built-in operators to the channel
        return chResults
            .collect { _meta, output -> output }
            .map { output -> [[id: toolName, args: args], output.toSet()] }
    } else {
        // Input is a list - process directly
        def outputs = chResults.collect { tuple -> tuple[1] }
        return [[id: toolName, args: args], outputs.toSet()]
    }
}
```

---

### flattenPaths

Consolidates multiple output channels into a single results channel by flattening file collections. Transforms `Tuple<Map, Set<Path>>` to `Tuple<Map, Path>`.

#### Signature

```groovy
flattenPaths(channels)
```

#### Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `channels` | List | Yes | List of channels, each containing `[meta, files]` tuples |

#### Returns

Channel of `[meta, file]` tuples - one tuple per individual file.

- If `files` is a `Set<Path>` or `List<Path>`, each file becomes its own tuple
- If `files` is a single `Path`, it remains as `[meta, file]`
- Multiple channels are mixed before flattening

#### Usage Example

```groovy
// Standard subworkflow output pattern
emit:
results  = flattenPaths([
    TOOL_MODULE.out.report,
    TOOL_MODULE.out.summary,
    TOOL_MODULE.out.json
])
logs     = flattenPaths([TOOL_MODULE.out.logs])
nf_logs  = flattenPaths([TOOL_MODULE.out.nf_logs])
versions = flattenPaths([TOOL_MODULE.out.versions])
```

#### Implementation

```groovy
static Object flattenPaths(List channels) {
    // Input validation
    if (channels == null) {
        throw new IllegalArgumentException("channels cannot be null")
    }
    if (channels.isEmpty()) {
        return []
    }
    if (channels.any { it == null }) {
        throw new IllegalArgumentException("channels list cannot contain null elements")
    }

    // Check if we're dealing with channels or lists
    def firstItem = channels[0]
    def isChannel = firstItem instanceof DataflowReadChannel || firstItem instanceof DataflowWriteChannel

    if (isChannel) {
        // Handle single channel case
        if (channels.size() == 1) {
            return channels[0].flatMap { tuple ->
                def meta, files
                if (tuple instanceof List && tuple.size() >= 2) {
                    meta = tuple[0]
                    files = tuple[1]
                } else {
                    return [tuple]  // Return as-is if structure is unexpected
                }

                // Check if files is a collection or single path
                if (files instanceof Collection) {
                    return files.collect { file -> [meta, file] }
                } else {
                    return [[meta, files]]
                }
            }
        }

        // Mix multiple channels and flatten
        def mixed = channels[0]
        channels[1..-1].each { ch -> mixed = mixed.mix(ch) }
        return mixed.flatMap { tuple ->
            def meta, files
            if (tuple.size() == 1 && tuple[0] instanceof List && tuple[0].size() == 2) {
                // Handle mixed channels where tuple is wrapped: [[meta, files]]
                meta = tuple[0][0]
                files = tuple[0][1]
            } else if (tuple.size() == 2) {
                // Handle normal case: [meta, files]
                meta = tuple[0]
                files = tuple[1]
            } else {
                return [tuple]  // Return as-is if structure is unexpected
            }

            if (files instanceof Collection) {
                return files.collect { file -> [meta, file] }
            } else {
                return [[meta, files]]
            }
        }
    } else {
        // Process lists directly
        def result = []
        channels.each { list ->
            list.each { tuple ->
                def meta = tuple[0]
                def files = tuple[1]

                if (files instanceof Collection) {
                    files.each { file -> result.add([meta, file]) }
                } else {
                    result.add([meta, files])
                }
            }
        }
        return result
    }
}
```

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

```
Module Outputs                    Subworkflow Aggregation
─────────────────                 ───────────────────────

TOOL.out.report ─┐
TOOL.out.json   ─┼─→ flattenPaths([...]) ─→ results (Tuple<Map, Path>)
TOOL.out.txt    ─┘

TOOL.out.logs   ───→ flattenPaths([...]) ─→ logs    (Tuple<Map, Path>)

TOOL.out.tsv ───────→ gather(..., 'tool') ─→ CSVTK_CONCAT ─→ merged.tsv
```

## Common Patterns

### Subworkflow Standard Output

```groovy
workflow TOOL_WORKFLOW {
    take:
    ch_input

    main:
    TOOL_MODULE(ch_input)

    emit:
    results  = flattenPaths([TOOL_MODULE.out.report, TOOL_MODULE.out.summary])
    logs     = flattenPaths([TOOL_MODULE.out.logs])
    nf_logs  = flattenPaths([TOOL_MODULE.out.nf_logs])
    versions = flattenPaths([TOOL_MODULE.out.versions])
}
```

### Aggregation for Summary

```groovy
// Collect all TSV outputs and concatenate into single file
CSVTK_CONCAT(
    gather(MLST.out.tsv, 'mlst'),
    'tsv',
    'tsv'
)
```

## See Also

- [Technical Specifications](../standards/03-technical-specs.md) - For channel type conventions
- [Examples](./01-examples.md) - For usage in context
- [Subworkflow Documentation](../standards/04-subworkflow-documentation.md) - For subworkflow patterns
