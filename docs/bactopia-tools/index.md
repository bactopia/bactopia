# Overview
After Bactopia has completed, there will be [a lot of output files](output-overview/) 
for each __individual__ sample. The next step is to take these results and compare 
samples between each other. To aid in these comparative analyses, a set of predefined 
workflows, called _Bactopia Tools_, have been created. These workflows use the 
predictable directory structure of Bactopia to automate analyses.

### Common Inputs
With the exceptions of the `summary` tool, each Bactopia Tool will use the following 
input parameters:
```
    --bactopia STR          Directory containing Bactopia analysis results for all samples.

    --exclude STR           A text file containing sample names to exclude from the
                                analysis. The expected format is a single sample per line.

    --include STR           A text file containing sample names to include in the
                                analysis. The expected format is a single sample per line.
```

#### `--bactiopia`
This parameter tells each tool where to find your Bactopia outputs from your project. 
Using this path, the tool will identify the required inputs and begin analysis. What 
this means is there is no need for you to wrangle up input files for compartive analyses.

#### `--exclude`
What `--exclude` allows is for you to give a text file with a list of samples that 
*should probably* be excluded from further analyses. While you can produce this list
yourself, the `summary` tool will produce a list of samples that do not pass certain 
thresholds. These thresholds are based on read lengths, sequence quality, coverage 
and assembly quality. You can adjust these thresholds to meet your needs.

#### `--include`
Similarly, `--include` allows you to give a text file with a list of samples to be 
included in the analysis. This allows you to target your anlyses on a specific subset
of samples. An example of this may be to use the `fastani` tool to determine samples
with >95% ANI to a reference, then create a pan-genome with the `roary` tool using 
only the subset of samples.

### Available Tools
Below is a list of Bactopia Tools currently available. To learn more about each, 
please follow the link.

__[eggnog](/bactopia-tools/eggnog/)__  
Functional annotation using orthologous groups

__[fastani](/bactopia-tools/fastani/)__  
Pairwise average nucleotide identity

__[gtdb](/bactopia-tools/gtdb/)__  
Identify marker genes and assign taxonomic classifications

__[ismapper](/bactopia-tools/ismapper/)__  
Identify positions of insertion sites 

__[mashtree](/bactopia-tools/mashtree/)__  
Tree based on Mash distances

__[phyloflash](/bactopia-tools/phyloflash/)__  
16s extraction, alignment, and tree  

__[pirate](/bactopia-tools/pirate/)__  
Pan-genome, core-genome alignment and tree 

__[roary](/bactopia-tools/roary/)__  
Pan-genome, core-genome alignment and tree 

__[summary](/bactopia-tools/summary/)__  
A report summarizing Bactopia project

### Suggest A Tool
Do you have an idea or suggestion for an analysis that should be added to the set 
of Bactopia Tools? If so, please feel free to submit it to 
[Bactopia GitHub Issues](https://github.com/bactopia/bactopia/issues)!
