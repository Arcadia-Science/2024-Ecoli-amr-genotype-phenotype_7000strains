# TODO: Creating a 7000 species genotype-phenotype dataset of *E. coli*  and antimicrobial resistance phenotype

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)
[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)

## Purpose

In this work, we aimed to produce a large-scale genotype-phenotype dataset that can be used as a resource to serve as a testbed for developing further strategies in genotype-phenotype investigations or evolutionary research.

We leveraged the genetic information and antimicrobial resistance (AMR) phenotype data available for the bacterium *Escherichia* coli to construct our dataset, and took advantage of the existing knowledge about genetic variations and AMR phenotypes to validate our approach and dataset. We performed variant calling and compiled a genotype-phenotype dataset for more than 7,000 *E. coli* strains. Briefly, variant calling consists of identifying all genetic variations and their associated genotypes in a population compared to a reference genome. This is performed by aligning sequencing reads for each sample of the population against a reference genome, then identifying polymorphic regions in the population, and finally characterizing variants and their genotypes at each of these polymorphic regions.

Here, we have generated a dataset that successfully revealed significant genetic diversity and identified 2.4 million variants. By focusing on non-silent variants within genes associated with AMR, we confirmed the dataset's accuracy and provided insights into specific mutations contributing to resistance to trimethoprim.


## Installation and Setup

This repository uses Snakemake, Python and R.

This repository uses Snakemake to run the pipeline and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
TODO: Replace <NAME> with the name of your environment
mamba env create -n <NAME> --file envs/dev.yml
conda activate <NAME>
```

Snakemake manages rule-specific environments via the `conda` directive and using environment files in the [envs/](./envs/) directory. Snakemake itself is installed in the main development conda environment as specified in the [dev.yml](./envs/dev.yml) file.

**Tips for Developers**

You can use the following command to export your current conda environment to a `yml` file.  
This command will only export the packages that you have installed directly, not the ones that were installed as dependencies. When you're ready to share, please delete this section.

```{bash}
conda env export --from-history --no-builds > envs/dev.yml
```

## Data

All the data related to this repository are available on Zenodo (ADD LINK), and organized into the same structure. When data size permits, we have also included relevant data (inputs or results) to this repository. Their names, locations, and the code they were used for or generated from are described in detail in the following sections.

## Overview

The primary objective of this study was to integrate available genotypic and phenotypic data of *E. coli* to map genetic variations to known antimicrobial resistance (AMR) phenotypes. We first identified an *E. coli* population for which AMR phenotypes and genomes were documented, and we further characterized the phenotype distribution within our cohort. Then, given the genetic diversity observed in *E. coli*, we decided to construct a pangenome as a comprehensive reference and used it to conduct variant calling and characterize the genetic diversity across our population. Finally, we correlated the identified genetic variations with the known AMR phenotypes.

### Repo organization

This repository is divided into two sections. 
The first section, **Dataset_generation**, includes the code and information necessary to build the genotype dataset and perform the variant calling when size permitted. It covers major steps like the generation of the reference pangenome used for variant calling, the variant calling pipeline applied to each of the 7,000 strains, the filtering of false positive variants, and the annotation of the variants. 
The second section, **Dataset_analysis**, includes the code and information used to process and analyze the dataset and generate figures for the Pub (https://doi.org/10.57844/arcadia-d2cf-ebe5). It includes the preliminary analysis of AMR phenotypes within the population, and the analysis of variants in regards to knowm AMR phenotypes.

### Approach

#### Population-wide variant calling

The first step to generate our dataset involved mapping out all the genetic variations within the population to a reference genome. This process led to the creation of a genotype matrix, which summarizes the genotype of each individual at each variant across the genome. The procedure included several key steps: selecting and constructing a reference genome against which genetic variants are identified, identifying genetic variants within each strain, integrating all the genetic variants and their corresponding genotypes from all strains to construct the genotype matrix, filtering false positives and annotating the variants.
The codes and some data associated with this part are all shared in the **Dataset_generation** folder of this repo. 
To obtain the all input and output data and be able to recapitulate this analysis, please download the Zenodo repository: LINK TO ZENODO

##### Generation of the reference genome

The selection of a good reference genome is important for genotype-phenotype analysis and the precise identification and annotation of genetic variants in the population. The reference genome must provide comprehensive coverage and accurately represent the genetic diversity of the population (https://doi.org/10.1099/mgen.0.001021). 
*E. coli* is a highly recombinogenic species, exhibiting high genetic diversity among its strains, so we need a reference genome that encompasses this global diversity. Therefore, we have generated a pangenome using the genomes of the ECOR collection, which consists of 72 *E. coli* strains.

From these strains, we constructed a pangenome containing both coding sequences and intergenic regions (IGRs). 

To construct the pangenome, we first downloaded the genome files for every 72 strains from (https://doi.org/10.1128/mra.01133-18) using Batch Entrez [link] and the assembly_accession numbers from this table: Dataset_generation/data/ECOR72_SRA_and_assembly_accessions.csv.

We then used a snakefile (Dataset_generation/scripts/Snakemake_ECOR72_annotation) to unzip, reannotate the genomes using Prokka (https://doi.org/10.1093/bioinformatics/btu153), and sort the generated gff files into a specific folder for further analysis with Roary (https://doi.org/10.1093/bioinformatics/btv421)
`snakemake -s Dataset_generation/scripts/Snakemake_ECOR72_annotation --cores 8`

Next, we used Roary to create the pangenome of coding sequences. We created the pangenome using a 90% identity threshold between the proteins (`i -90`) and included the flag -mafft to use the aligniment function and obtain the final fasta file containing all the pangenome sequences.
We ran the follwoing command line:
`roary -e --mafft -f Dataset_generation/results/pangenome_cds -i 90 Gff_files/*.gff -p 8`  
where Gff_files is the location of all the previously generated Prokka Gff files.
The two main outputs of Roary used in this work are Dataset_generation/results/pangenome_cds/gene_presence_absence.csv, providing the presence-absence information in the different ECOR strains for all the identified genes (or locus) identified in the pangenome, and Dataset_generation/results/pangenome_cds/cds_sequences.fa a multi-sequence fasta file containing all sequences of the coding sequences pangenome ; We also share the summary file Dataset_generation/results/pangenome_cds/summary_statistics.txt
The rest of the output is available on Zenodo: ADD LINK

We further utilized Piggy (https://doi.org/10.1093/gigascience/giy015) to generate the pangenome of IGRs. We installed Piggy following the directions provided on the GitHub repository: https://github.com/harry-thorpe/piggy (`git clone https://github.com/harry-thorpe/piggy.git`). Running Piggy relies on the previous Roary run to create the pangenome of IGRs. Make sure to include all the Roary outputs (available on Zenodo) and not only the 3 files provided in this repository
Then, we ran the following command line to obtain Piggy output 
`piggy/bin/piggy --in_dir Gff_files --out_dir Dataset_generation/results/pangenome_igr --roary_dir Dataset_generation/results/pangenome_cds -t 5`
Again, only the two main outputs of Piggy used in this project are shared on this reposistory: Dataset_generation/results/pangenome_igr/IGR_presence_absence.csv and Dataset_generation/results/pangenome_igr/pangenome_igr.fasta 
The rest of the output is available on Zenodo: ADD LINK

Finally, we concatenated the fasta outputs of Roary and Piggy to generate the whole pangenome, including both coding sequences and IGRs (Dataset_generation/results/pangenome_whole/whole_pangenome.fasta)
`cat Dataset_generation/results/pangenome_cds/cds_sequences.fa Dataset_generation/results/pangenome_igr/pangenome_igr.fasta > Dataset_generation/results/pangenome_whole/whole_pangenome.fasta`

Eventually, we indexed the pangenome using `bwa index`
`bwa index Dataset_generation/results/pangenome_whole/whole_pangenome.fasta`

##### Download of sequencing reads from SRA
To come

##### Variant Calling in 7055 samples and merging
To come

##### Filtering variants
To come

##### Variant annotation
To come

#### Phenotype distribution analysis
To come

#### Antimicrobial resistance analysis
To come


**Tips for Developers**

You should consider having a quickstart guide for users who want to run the pipeline, and/or a demo dataset that they can use to test the pipeline.  
When you're ready to share, please delete this section.

### Compute Specifications

TODO: Describe what compute resources were used to run the analysis. For example, you could list the operating system, number of cores, RAM, and storage space.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).

---
## For Developers

This section contains information for developers who are working off of this template. Please delete this section when you're ready to share your repository.

### GitHub templates
This template uses GitHub templates to provide checklists when making new pull requests as well as templates for issues, which could be used to request new features or report bugs. These templates are stored in the [.github/](./.github/) directory.

### VSCode
This template includes recommendations to VSCode users for extensions, particularly the `ruff` linter. These recommendations are stored in `.vscode/extensions.json`. When you open the repository in VSCode, you should see a prompt to install the recommended extensions. 

### `.gitignore`
This template uses a `.gitignore` file to prevent certain files from being committed to the repository.

### `pyproject.toml`
`pyproject.toml` is a configuration file to specify your project's metadata and to set the behavior of other tools such as linters, type checkers etc. You can learn more [here](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)

### Linting
This template automates linting and formatting using GitHub Actions and the `ruff` and `snakefmt` linters. When you push changes to your repository, GitHub will automatically run the linter and report any errors, blocking merges until they are resolved.

### Testing
This template uses GitHub Actions to automate a test dry run of the pipeline. When you push changes to your repository, GitHub will automatically run the tests and report any errors, blocking merges until they are resolved.
