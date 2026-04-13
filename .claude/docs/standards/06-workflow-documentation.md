# Entry Workflow Documentation Methodology

## Overview
This guide provides a comprehensive methodology for creating consistent and accurate GroovyDoc documentation for Bactopia entry workflows. Entry workflows are user-facing entry points that orchestrate subworkflows to produce published results.

## 1. Workflow Architecture

### 1.1 Entry Workflow Characteristics
Entry workflows are fundamentally different from modules and subworkflows:

1. **User-Facing Entry Points**: Workflows are what end users actually run
2. **@publish instead of @output**: Workflows publish files, don't emit channels
3. **Parameter-Based Inputs**: Workflows accept params, not typed channels
4. **Section Organization**: Results are organized into @section groups
5. **Parameter-Typed Tags**: Workflows use `@tags` with `input-type:parameter` (since inputs are params, not channels)
6. **Narrative Documentation**: More descriptive, user-focused explanations

### 1.2 Standard Workflow Components
- **Header Block**: GroovyDoc with user-facing description
- **Parameters Block**: Workflow-level input parameters
- **Include Statements**: Subworkflow dependencies
- **Workflow Definition**: Main orchestration logic
- **Publish Block**: Published file outputs

## 2. GroovyDoc Standard

### 2.1 Required Structure
Entry workflows must follow this exact structure:

```groovy
#!/usr/bin/env nextflow
/**
 * <Concise description of what the workflow does>.
 *
 * <Detailed description explaining the pipeline, tools used, and outputs>.
 * Include [ToolName](URL) links with https protocol.
 *
 * @status stable
 * @keywords <comma, separated, relevant keywords>
 * @tags complexity:<level> input-type:<type> output-type:<type> features:<features>
 * @citation <tool_name>  // Only if specific tool citation needed
 *
 * @subworkflows <comma separated list of subworkflows>
 *
 * @input <param_name>
 * <Description of the parameter>
 *
 * @input <optional_param>
 * <Description of optional parameter>
 *
 * @section <Section Name>
 * @note <Optional note about conditions>
 * @publish <file_pattern>    <Description of published files>
 * @publish <file_pattern>    <Description of published files>
 *
 * @section <Another Section>
 * @publish <file_pattern>    <Description of published files>
 */
nextflow.preview.types = true
```

**Important**: The GroovyDoc must immediately follow the shebang line (`#!/usr/bin/env nextflow`). The `nextflow.preview.types = true` declaration can come after the GroovyDoc.

## 3. Header Format

### 3.1 Title Patterns

#### Bactopia Tools
```groovy
/**
 * Mass screening of contigs for antimicrobial resistance and virulence genes.
 *
 * Performs mass screening of contigs for antimicrobial resistance and virulence genes
 * using [Abricate](https://github.com/tseemann/abricate) to screen assemblies against
 * multiple resistance and virulence gene databases.
 * ...
 */
```

#### Named Workflows
```groovy
/**
 * Taxonomic classification and abundance profiling of metagenomic reads.
 *
 * This workflow performs taxonomic classification using [Kraken2](https://github.com/DerrickWood/kraken2)
 * and [Bracken](https://github.com/jenniferlu717/Bracken), with optional host read removal.
 * ...
 */
```

### 3.2 Required Fields
- **@status**: Always `stable`
- **@keywords**: Comma-separated list relevant to the analysis type
- **@citation**: Only if specific tool citation needed (not always required)
- **@subworkflows**: Comma-separated list of subworkflows used

## 4. Parameter Documentation

### 4.1 Parameter Format

```groovy
@input <param_name>
<Description of the parameter>

@input <optional_param>
<Description of optional parameter> (Optional)
```

### 4.2 Generic Bactopia-Tools Parameters

All bactopia-tools workflows inherit common parameters that are typically NOT documented in the individual workflow GroovyDoc (they are inherited from the framework):
- `bactopia` - Path to Bactopia analysis directory
- `includes` - Sample inclusion filter
- `excludes` - Sample exclusion filter
- `workflow` - Internal workflow identifier

These are handled automatically by the bactopia-tools initialization subworkflow.

### 4.3 Common Parameter Types

#### Required Parameters
```groovy
@input rundir
Directory containing raw sequencing reads

@input kraken2_db
Path to Kraken2 database for classification
```

#### Optional Parameters
```groovy
@input use_srascrubber
Remove host reads using SRA scrubber before classification

@input use_bakta
Use Bakta for genome annotation instead of Prokka
```

#### Database Parameters
```groovy
@input bakta_db
Path to Bakta database for annotation

@input adapters
Path to adapter sequences file for removal during QC
```

## 5. Section Organization

### 5.1 Standard Sections

#### Quality Control (when applicable)
```groovy
@section Quality Control
@publish fastqc/*   FastQC quality control reports
@publish multiqc/*   MultiQC aggregated quality reports
```

#### Assembly
```groovy
@section Assembly
@publish *.fasta   Assembled genome sequences
@publish assembly-stats.txt Assembly quality metrics
```

#### Annotation
```groovy
@section Annotation
@publish *.gff    Genome annotation in GFF3 format
@publish *.gbk    Genome annotation in GenBank format
@publish *.faa    Protein sequences
@publish *.tsv    Annotation summary tables
```

#### Typing
```groovy
@section Typing
@publish mlst.txt   MLST sequence type results
```

#### Antimicrobial Resistance
```groovy
@section Antimicrobial Resistance
@publish amrfinderplus.tsv AMR gene detection results
```

#### Comparative Analysis
```groovy
@section Comparative Analysis
@publish mash-dist.tsv Mash distance matrix
@publish sketch.msh Mash sketch files
```

#### Standard Sections (always included)
```groovy
@section Execution Logs
@publish logs/**   Tool execution logs
@publish logs/nf-* Nextflow execution scripts and logs for debugging

@section Versions
@publish versions.yml Software version information
```

### 5.2 Section Patterns
- **Logical Grouping**: Related outputs grouped together
- **User-Centric**: Sections make sense from user perspective
- **Progressive**: Follows typical analysis pipeline flow
- **Consistent Naming**: Standardized section titles across workflows

## 6. Publishing Files

### 6.1 @publish Format
```groovy
@publish <file_pattern>    <Description of published files>
```

### 6.2 File Pattern Examples

#### Single Files
```groovy
@publish versions.yml Software version information
@publish abricate.tsv    Merged TSV file with Abricate results from all samples
```

#### Pattern Matching
```groovy
@publish *.txt    Tab-delimited report of screening results
@publish *.fasta  Assembled genome sequences
@publish logs/**  Tool execution logs
```

#### Directories
```groovy
@publish fastqc/*   FastQC quality control reports
@publish multiqc/*   MultiQC aggregated quality reports
```

### 6.3 Conditional Outputs
```groovy
@section Pathogen-Specific Analysis
@note Only created if --ask_merlin is enabled
@publish merlin/*    Pathogen-specific analysis results
```

## 7. Examples by Complexity

### 7.1 Simple Bactopia Tool (abricate)
```groovy
/**
 * Mass screening of contigs for antimicrobial resistance and virulence genes.
 *
 * This Bactopia Tool uses [Abricate](https://github.com/tseemann/abricate) to screen
 * assemblies against multiple resistance and virulence gene databases, including
 * NCBI, CARD, RESFINDER, ARG-ANNOT, VFDB, PLASMIDFINDER, ECOLI_VF, and MEGARES.
 * It processes a Bactopia analysis directory, runs Abricate on each sample, and
 * creates a merged summary report.
 *
 * @status stable
 * @keywords bacteria, antimicrobial resistance, virulence, screening, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation abricate, arg_annot, card, csvtk, ecoh, megares2, ncbi_reference_gene_catalog, plasmidfinder, resfinder, vfdb
 *
 * @subworkflows abricate, utils_bactopia-tools
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv            Tab-delimited report of Abricate screening results for each sample
 *
 * @section Merged Results
 * @publish abricate.tsv     Merged TSV report containing Abricate results from all samples
 *
 * @section Execution Logs
 * @publish logs/abricate/*  Tool execution logs (stdout/stderr)
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml     Software version information
 */
```

### 7.2 Complex Workflow with Parameters (teton)
```groovy
/**
 * Taxonomic classification and abundance profiling of metagenomic reads.
 *
 * This workflow performs metagenomic classification using [Kraken2](https://github.com/DerrickWood/kraken2)
 * and [Bracken](https://github.com/jenniferlu717/Bracken), with optional host read removal
 * using SRA Scrubber. It processes metagenomic sequencing reads to estimate bacterial
 * genome sizes and separate bacterial from non-bacterial organisms.
 *
 * @status stable
 * @keywords metagenomics, classification, kraken2, bracken, abundance, profiling
 * @tags complexity:complex input-type:parameter output-type:multiple features:aggregation,conditional-logic,database-dependent
 *
 * @subworkflows bactopia_gather, teton, utils_bactopia
 *
 * @input rundir
 * Directory containing metagenomic sequencing reads
 *
 * @input kraken2_db
 * Path to Kraken2 database for classification
 *
 * @input use_srascrubber
 * Remove host reads using SRA scrubber before classification
 *
 * @section Classification Results
 * @publish bacteria.tsv            Per-sample TSV files containing bacterial organisms and their properties
 * @publish nonbacteria.tsv         Per-sample TSV files containing non-bacterial organisms
 * @publish sizemeup.tsv            Per-sample TSV files with genome size estimates
 *
 * @section Merged Results
 * @publish merged-bacteria.tsv     Consolidated TSV file of all bacterial organisms across samples
 * @publish merged-nonbacteria.tsv  Consolidated TSV file of all non-bacterial organisms across samples
 * @publish merged-sizemeup.tsv     Consolidated TSV file of genome size estimates across samples
 * @publish report.tsv              Joined TSV file combining scrubber and classification results
 *
 * @section Execution Logs
 * @publish logs/**                 Tool execution logs
 * @publish logs/nf-*               Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml            Software version information
 */
```

### 7.3 Main Bactopia Workflow (complex)
```groovy
/**
 * Comprehensive bacterial analysis pipeline for complete genomic characterization.
 *
 * This workflow performs end-to-end analysis including quality control, assembly,
 * annotation, antimicrobial resistance detection, MLST typing, and optional
 * pathogen-specific analysis through Merlin. It processes raw sequencing reads
 * and produces a complete genomic characterization suitable for downstream analysis.
 *
 * @status stable
 * @keywords bacteria, assembly, annotation, AMR, MLST, genomics, pipeline
 * @tags complexity:complex input-type:parameter output-type:multiple features:aggregation,conditional-logic
 *
 * @subworkflows amrfinderplus, bactopia_assembler, bactopia_datasets, bactopia_gather,
 *               bactopia_qc, bactopia_sketcher, bakta, merlin, mlst, prokka, utils_bactopia
 *
 * @input rundir
 * Directory containing raw sequencing reads
 *
 * @input adapters
 * Path to adapter sequences file for removal during QC
 *
 * @input phix
 * Path to PhiX sequences for contamination removal during QC
 *
 * @input use_bakta
 * Use Bakta for genome annotation instead of Prokka
 *
 * @input bakta_db
 * Path to Bakta database for annotation
 *
 * @input download_bakta
 * Download Bakta database if not provided
 *
 * @input bakta_save_as_tarball
 * Save Bakta database as tarball for reuse
 *
 * @input bakta_proteins
 * Path to trusted protein sequences for Bakta annotation
 *
 * @input bakta_prodigal_tf
 * Path to Prodigal training file for Bakta
 *
 * @input bakta_replicons
 * Path to replicon sequences for Bakta
 *
 * @input prokka_proteins
 * Path to protein sequences for Prokka annotation
 *
 * @input prokka_prodigal_tf
 * Path to Prodigal training file for Prokka
 *
 * @input emmtyper_blastdb
 * Path to emmtyper BLAST database for Merlin
 *
 * @input hicap_database_dir
 * Path to HiCap database directory for Merlin
 *
 * @input hicap_model_fp
 * Path to HiCap model file for Merlin
 *
 * @input ask_merlin
 * Enable Merlin pathogen-specific analysis
 *
 * @input spatyper_repeats
 * Path to Spatyper repeats database for Merlin
 *
 * @input spatyper_repeat_order
 * Path to Spatyper repeat order file for Merlin
 *
 * @section Quality Control
 * @publish fastqc/*   FastQC quality control reports
 * @publish multiqc/*  MultiQC aggregated quality reports
 *
 * @section Assembly
 * @publish *.fasta             Assembled genome sequences
 * @publish assembly-stats.txt  Assembly quality metrics
 * @publish quast.html          QUAST assembly quality report
 *
 * @section Annotation
 * @publish *.gff               Genome annotation in GFF3 format
 * @publish *.gbk               Genome annotation in GenBank format
 * @publish *.faa               Protein sequences
 * @publish *.fna               Nucleotide sequences
 * @publish *.tsv               Annotation summary tables
 *
 * @section Typing
 * @publish mlst.txt            MLST sequence type results
 *
 * @section Antimicrobial Resistance
 * @publish amrfinderplus.tsv           AMR gene detection results
 * @publish amrfinderplus.mutation.tsv  AMR mutation results
 *
 * @section Comparative Analysis
 * @publish mash-dist.tsv       Mash distance matrix
 * @publish sketch.msh          Mash sketch files
 * @publish sourmash.sig        Sourmash signatures
 *
 * @section Pathogen-Specific Analysis
 * @note Only created if --ask_merlin is enabled
 * @publish merlin/             Merlin pathogen-specific analysis results
 *
 * @section Execution Logs
 * @publish logs/**             Tool execution logs
 * @publish logs/nf-*           Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml        Software version information
 */
```

## 8. Implementation Patterns

### 8.1 Standard Publish and Output Block Pattern

All workflows follow the same record-based publish/output structure. Since every subworkflow emits both `sample_outputs` and `run_outputs`, every workflow uses the same template:

```groovy
include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { TOOL                } from '../../../subworkflows/<tool>/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_tool = TOOL(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_tool.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_tool.sample_outputs)
    // Run-level
    run_outputs = ch_tool.run_outputs
    run_nf_logs = collectNextflowLogs(ch_tool.run_outputs)
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
```

#### Key rules:
- Subworkflow results are captured via direct assignment (`ch_tool = TOOL(...)`) and read as `ch_tool.sample_outputs` / `ch_tool.run_outputs` — do not use `TOOL.out.X` access.
- `flatten()` is called in the `output {}` block because `results` is a list of `files()` outputs (list of lists). This is a **publishing concern** -- subworkflows should NOT flatten.
- Every workflow includes both sample and run output sections. When a subworkflow uses `channel.empty()` for `run_outputs`, the closures simply never fire.
- The nf_logs extraction is delegated to the `collectNextflowLogs()` plugin function from `nf-bactopia`. It does the `flatMap`+`collect` internally and returns `[meta, file]` tuples ready for the `path { meta, f -> ... }` closures — workflows just import and call it.

### 8.2 Parameter Declarations

```groovy
params {
    rundir : String

    // Tool-specific parameters
    kraken2_db : Path
    use_srascrubber : Boolean
}
```

## 9. Quality Checklist

Before completing workflow documentation, verify:

- [ ] Shebang line is `#!/usr/bin/env nextflow`
- [ ] GroovyDoc immediately follows the shebang line
- [ ] `nextflow.preview.types = true` is included (usually after GroovyDoc)
- [ ] Header format matches pattern (Bactopia Tool vs Workflow Name)
- [ ] All required fields are present (@status, @keywords)
- [ ] Parameters use simple @input format (not @input tuple)
- [ ] Sections are logically organized and consistently named
- [ ] All outputs use @publish with clear descriptions
- [ ] File patterns are accurate and match actual outputs
- [ ] Conditional outputs include @note explaining conditions
- [ ] Subworkflows are listed with @subworkflows
- [ ] Links use https protocol
- [ ] @tags and complexity ratings are correctly specified

## 10. Missing Parameters

### 10.1 Identifying Missing Parameters

When documenting workflows, verify that all parameters referenced in the workflow logic are properly declared in the `params` block.

1. **Review the workflow body** - Look for any `params.XYZ` references
2. **Compare with params block** - Ensure all referenced parameters are declared
3. **Add missing parameters** with appropriate types

### 10.2 Adding Missing Parameters

When you find missing parameters, add them using this template:

```groovy
params {
    // Existing parameters...
    rundir : String

    // Added parameters with descriptive comments
    // <CATEGORY> parameters
    param_name : Type        // Brief description of purpose
    another_param : Path?    // Optional parameter with nullable type
}
```

### 10.3 Parameter Categories

Organize parameters with descriptive comments:

- **Input parameters**: Data paths and input specifications
- **Tool selection**: Boolean flags for choosing between tools
- **Reference data**: Paths to databases or reference files
- **Analysis options**: Parameters controlling analysis behavior
- **Resource settings**: Thread counts, memory limits, etc.

### 10.4 Example

Before (missing parameters):
```groovy
params {
    rundir : String
}

workflow {
    // Later in workflow...
    if (params.use_feature) {
        process_file(params.input_file)
    }
}
```

After (parameters added):
```groovy
params {
    rundir : String

    // Tool-specific parameters
    // Feature selection parameters
    use_feature : Boolean

    // Input parameters
    input_file : Path?
}
```

**Important**: Never add parameters to the documentation that don't exist in the actual params block. Always verify the implementation matches the documentation.

## 11. What to Avoid

1. **Do not** place anything between the shebang and GroovyDoc
2. **Do not** use @output in workflow GroovyDoc (only @publish)
3. **Do not** omit @tags or complexity ratings
4. **Do not** reference internal implementation details
5. **Do not** include Nextflow channel types
6. **Do not** use tuple formats for inputs
7. **Do not** include parameter defaults in documentation
8. **Do not** reference module.config or schema.json
9. **Do not** document parameters that don't exist in the params block

## 12. Vertical Alignment of @publish Tags

### 12.1 Alignment Methodology

To ensure consistent formatting across all workflows, @publish tags must be vertically aligned. Follow this process:

1. **Identify all file patterns** in the workflow and count their characters
2. **Find the longest file pattern** - this will be the baseline
3. **Determine target alignment**: The longest pattern should have exactly 2 spaces after it
4. **Calculate required padding** for each pattern:

   ```
   padding_needed = target_column - (11 + length(pattern))
   where target_column = 11 + length(longest_pattern) + 2
   ```

5. **Apply consistent spacing** so all descriptions start at the same column

### 12.2 Example

If file patterns are:
- `*.txt` (5 chars)
- `*.summary.txt` (13 chars) ← Longest
- `report.tsv` (10 chars)

Target column = 11 + 13 + 2 = 26

Result:

```groovy
@publish *.txt           Per-sample analysis results
@publish *.summary.txt   Summary statistics
@publish report.tsv      Consolidated report
```

## 13. Common Section Templates

### 13.1 Basic Sections

```groovy
@section Per-Sample Results
@publish *.txt      Per-sample analysis results

@section Merged Results
@publish merged.tsv Consolidated results from all samples
```

### 13.2 Analysis Type Sections
```groovy
@section Quality Control
@publish qc/*          Quality control metrics and reports

@section Assembly
@publish *.fasta       Assembled genome sequences
@publish *.stats       Assembly statistics

@section Annotation
@publish *.gff         Genome annotations in GFF3 format
@publish *.gbk         Annotations in GenBank format

@section Typing
@publish typing.txt    MLST or other typing results

@section Comparative Analysis
@publish distance.tsv  Distance matrix
@publish matrix.tsv    Presence/absence matrix
```

## 14. Future Considerations

### 14.1 Documentation Updates
When updating workflow documentation:
- Maintain consistency with established patterns
- Update sections when new features are added
- Keep parameter descriptions aligned with actual usage
- Ensure file patterns match actual outputs

### 14.2 New Workflow Creation
When creating new workflows:
- Follow the standard structure
- Use appropriate section organization
- Include all standard sections (Execution Logs, Versions)
- Document all parameters clearly
- Test that @publish patterns match actual outputs

## See Also
- [Style Guide](01-style-guide.md) - For general GroovyDoc templates
- [Module Documentation](05-module-documentation.md) - For module-specific documentation standards
- [Subworkflow Documentation](04-subworkflow-documentation.md) - For subworkflow-specific documentation standards
