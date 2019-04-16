# About Bactopia
We (Tim Read and myself) initially released [Staphopia](https://staphopia.emory.edu/) in the summer 
of 2018, for the analysis of *Staphylococcus aureus* genomes. At conferences,
or via email, the most common question about Staphopia that we received was:

*"How can I use Staphopia for my bacteria of interest?"*.

We too, have always been interested in this, because we tend to work with other bacteria and would like a standard pipeline for incoming sequencing projects. Given both the community and self interest, we developed Bactopia, *a generic pipeline for the analysis of short-read bacterial sequences*! 

Side-note... we kept the *"-opia"* in the name to pay homage to Staphopia!

While the philosophy behind Staphopia and Bactopia is mostly the same. The
development of Bactopia from scratch has really allowed us to take what we
learned from Staphopia and use it to our advantage. Bactopia has been developed
with usability, portability, and speed in mind from the start.

Bactopia uses [Nextflow](https://www.nextflow.io/) to manage the workflow. We
have also prioritized software packages available from
[Bioconda](https://bioconda.github.io/) (or other
[Anaconda channels](https://anaconda.org/)) to make installation
as simple as possible for *all* users. This has also given us the opportunity to update out workflow with the latest methods.

## Workflow Overview

**INSERT IMAGE OF WORKFLOW**


