# Creating a 7000 species genotype-phenotype dataset of *E. coli*  and antimicrobial resistance phenotype

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)
[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)

## Purpose

In this work, we aimed to produce a large-scale genotype-phenotype dataset that can be used as a resource to serve as a testbed for developing further strategies in genotype-phenotype investigations or evolutionary research.

We leveraged the genetic information and antimicrobial resistance (AMR) phenotype data available for the bacterium *Escherichia coli* to construct our dataset, and took advantage of the existing knowledge about genetic variations and AMR phenotypes to validate our approach and dataset. We performed variant calling and compiled a genotype-phenotype dataset for more than 7,000 *E. coli* strains. Briefly, variant calling consists of identifying all genetic variations and their associated genotypes in a population compared to a reference genome. This is performed by aligning sequencing reads for each sample of the population against a reference genome, then identifying polymorphic regions in the population, and finally characterizing variants and their genotypes at each of these polymorphic regions.

Here, we have generated a dataset that successfully revealed significant genetic diversity and identified 2.4 million variants. By focusing on non-silent variants within genes associated with AMR, we confirmed the dataset's accuracy and provided insights into specific mutations contributing to resistance to trimethoprim.


## Installation and Setup

This repository uses Snakemake, Python and R.

This repository uses Snakemake to run two different the pipelines that will construct the dataset, Python and R scripts that analyse the dataset, and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). 

After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

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

TODO: Add details about the description of input / output data and links to Zenodo depositions, if applicable.

If you want to reproduce this work or access the data, please download the corresponding folder from Zenodo (ADD LINK). The paths mentioned in this ReadMe refer to the organization of the Zenodo folder, which is mirrored in this GitHub repository as well, and we have included relevant data here when their size permits.

## Overview

The primary objective of this study was to integrate available genotypic and phenotypic data of *E. coli* to map genetic variations to known antimicrobial resistance (AMR) phenotypes. We first identified an *E. coli* population for which AMR phenotypes and genomes were documented, and we further characterized the phenotype distribution within our cohort. Then, given the genetic diversity observed in *E. coli*, we decided to construct a pangenome as a comprehensive reference and used it to conduct variant calling and characterize the genetic diversity across our population. Finally, we correlated the identified genetic variations with the known AMR phenotypes.

### Repo organization

This repository is divided into two sections. 
The first section, **dataset_generation**, includes the code and information necessary to build the genotype dataset and perform the variant calling when size permitted. It covers major steps like the generation of the reference pangenome used for variant calling, the variant calling pipeline applied to each of the 7,000 strains, the filtering of false positive variants, and the annotation of the variants. 
The second section, **dataset_analysis**, includes the code and information used to process and analyze the dataset and generate figures for the Pub (https://doi.org/10.57844/arcadia-d2cf-ebe5). It includes the preliminary analysis of AMR phenotypes within the population, and the analysis of variants in regards to known AMR phenotypes.

### Approach

#### Population-wide variant calling

The first step to generate our dataset involved mapping out all the genetic variations within the population to a reference genome. This process led to the creation of a genotype matrix, which summarizes the genotype of each individual at each variant across the genome. The procedure included several key steps: selecting and constructing a reference genome against which genetic variants are identified, identifying genetic variants within each strain, integrating all the genetic variants and their corresponding genotypes from all strains to construct the genotype matrix, filtering false positives and annotating the variants.
The codes and some data associated with this part are all shared in the **dataset_generation** folder of this repo. 
To obtain the all input and output data and be able to recapitulate this analysis, please download the Zenodo repository
### TODO: ADD ZENODO link

##### Generation of the reference genome

The selection of a good reference genome is important for genotype-phenotype analysis and the precise identification and annotation of genetic variants in the population. The reference genome must provide comprehensive coverage and accurately represent the genetic diversity of the population (https://doi.org/10.1099/mgen.0.001021). 
*E. coli* is a highly recombinogenic species, exhibiting high genetic diversity among its strains, so we need a reference genome that encompasses this global diversity. Therefore, we have generated a pangenome using the genomes of the ECOR collection, which consists of 72 *E. coli* strains isolated from a wide variety of hosts and geographical locations. This collection offers a broad representation of the natural diversity of the species (https://doi.org/10.1128/jb.157.2.690-693.1984).

From these strains, we constructed a pangenome containing both coding sequences and intergenic regions (IGRs). 

To construct the pangenome, we first downloaded the genome files for each of the 72 strains from (https://doi.org/10.1128/mra.01133-18) using Batch Entrez API [link](https://www.ncbi.nlm.nih.gov/sites/batchentrez) and the assembly_accession numbers from this table: dataset_generation/data/ECOR72_SRA_and_assembly_accessions.csv.

We then used a snakefile (data_generation/scripts/Snakemake_ECOR72_annotation) to unzip, reannotate the genomes using Prokka (https://doi.org/10.1093/bioinformatics/btu153), and sort the generated gff files into a specific folder for further analysis with Roary (https://doi.org/10.1093/bioinformatics/btv421)
`snakemake -s dataset_generation/scripts/Snakemake_ECOR72_annotation --cores 8`

Next, we used Roary to create the pangenome of coding sequences. We created the pangenome using a 90% identity threshold between the proteins (`i -90`) and included the flag `-mafft` to use the aligniment function and obtain the final fasta file containing all the pangenome sequences.
We ran the following command line:
`roary -e --mafft -f dataset_generation/results/pangenome_cds -i 90 Gff_files/*.gff -p 8`  
where Gff_files is the location of all the previously generated Prokka Gff files.
The two main outputs of Roary used in this work are dataset_generation/results/pangenome_cds/gene_presence_absence.csv, providing the presence-absence information in the different ECOR strains for all the identified genes (or locus) identified in the pangenome, and dataset_generation/results/pangenome_cds/cds_sequences.fa a multi-sequence fasta file containing all sequences of the coding sequences pangenome ; We also share the summary file dataset_generation/results/pangenome_cds/summary_statistics.txt
The rest of the output is available on Zenodo: ADD LINK

We further utilized Piggy (https://doi.org/10.1093/gigascience/giy015) to generate the pangenome of IGRs. We installed Piggy following the directions provided on the GitHub repository: https://github.com/harry-thorpe/piggy (`git clone https://github.com/harry-thorpe/piggy.git` - Commit Hash 68079ae1c310865d9d3a54221f8f3b3993329081). 
Running Piggy relies on the previous Roary run to create the pangenome of IGRs. Make sure to include all the Roary outputs (available on Zenodo) and not only the 3 files provided in this repository
Then, we ran the following command line to obtain Piggy output 
`piggy/bin/piggy --in_dir Gff_files --out_dir dataset_generation/results/pangenome_igr --roary_dir dataset_generation/results/pangenome_cds -t 5`
Again, only the two main outputs of Piggy used in this project are shared on this reposistory: dataset_generation/results/pangenome_igr/IGR_presence_absence.csv and dataset_generation/results/pangenome_igr/pangenome_igr.fasta 
The rest of the output is available on Zenodo: 
### TODO: ADD ZENODO link

Finally, we concatenated the fasta outputs of Roary and Piggy to generate the whole pangenome, including both coding sequences and IGRs (dataset_generation/results/pangenome_whole/whole_pangenome.fasta)
`cat dataset_generation/results/pangenome_cds/cds_sequences.fa dataset_generation/results/pangenome_igr/pangenome_igr.fasta > dataset_generation/results/pangenome_whole/whole_pangenome.fasta`

Eventually, we indexed the pangenome using `bwa index`
`bwa index dataset_generation/results/pangenome_whole/whole_pangenome.fasta`

##### Download of sequencing reads from SRA

We obtained the SRA accession numbers for the sequencing files of the 7,055 selected strains (6,983 strains from the population and 72 strains from the ECOR collection).
`dataset_generation/data/sample_list_SRA.csv`

We used the GNU Parallel shell tool and the fasterq-dump tool from the sra-toolkit (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit), and downloaded the FASTQ files from paired-end sequencing for each strain. 

`parallel -j 64 fasterq-dump {} ::: $(cat sample_list_SRA.csv` 

##### Variant Calling in 7,055 samples and merging of individual variant calling format (VCF) files

Variant calling aims to identify the differences between a strain and a reference genome by aligning sequencing reads of the strains against the reference and identifying where they differ. Variant calling algorithms then assess these differences to determine the likely genetic variations versus sequencing errors. The confirmed variants are then recorded in a Variant Call Format (VCF) file that lists the position of each variant in the genome, the nature of the genetic change, and the quality of the call. 

**Variant calling on individual strains**
We performed variant calling independently on eac strain againts the pangenome, and generated a VCF file for each of them.
The variant calling workflow is integrated into a Snakefile `data_generation/scripts/variant_calling_pipeline`, allowing for efficient parallel processing of multiple samples.

To run the snakefile, we used the following command:
`nohup snakemake -s dataset_generation/scripts/variant_calling_pipeline --cores 54 > variant_calling.txt 2>&1 &`
 
The different steps of the workflow include (for each strain):
	-> input file: whole genome paired fastq files R1 and R2
	- fastq quality control, trimming of common Illumina adapter, filtering short reads using FastP (https://doi.org/10.1002/imt2.107)
	- concatenation of filtered fastq1 and fastq2  -> intermediary file saved
	- alignment against the pangenome using bwa mem (
https://doi.org/10.48550/arXiv.1303.3997) -> intermediary file saved
	- sorting and marking duplicates in alignment files using Samtools (https://doi.org/10.1093/gigascience/giab008)
	- addition of RG tags using Picard (https://broadinstitute.github.io/picard/)
	- files indexing using Samtools
	- generating mpileup files using bcftools (https://doi.org/10.1093/gigascience/giab008) -> intermediary file saved
	- variant calling on mpileup file using bcftools and generation of VCF file
	- indexing VCF file
    -> output: VCF file 

TODO: mention where individual VCF files and intermediary files are available 

**Merging individual VCF.gz files**
We merged all VCF files into a single VCF (`data_generation/results/vcf/merged_output_all.vcf.gz`). To do so, we first merged batches of 1,000 files using the `merge` function from bcftools. Each batch was then re-indexed before a final comprehensive merge was performed.

To batch merge the VCF files, we first created 7 lists of 1000 strains as txt files, and for each list we used the following commands to merge and then index the VCF file (example with list1, but we ran the same command for the seven lists)

```{bash}
nohup bcftools merge -O z -o dataset_generation/results/vcf/merged_output_list7.vcf.gz -l dataset_generation/data/vcf_merging/List_7_merging.txt > 2>&1 &
bcftools index dataset_generation/results/vcf/merged_vcf/merged_output_list1.vcf.gz
```

To eventually merge everything we ran the following command:

```{bash}
bcftools merge dataset_generation/results/vcf/merged_output_list1.vcf.gz dataset_generation/results/vcf/merged_output_list2.vcf.gz dataset_generation/results/vcf/merged_output_list3.vcf.gz dataset_generation/results/vcf/merged_output_list4.vcf.gz dataset_generation/results/vcf/merged_output_list5.vcf.gz dataset_generation/results/vcf/merged_output_list6.vcf.gz dataset_generation/results/vcf/merged_output_list7.vcf.gz -O z -o dataset_generation/results/vcf/merged_output_all.vcf.gz

bcftools index dataset_generation/results/vcf/merged_output_all.vcf.gz
```


##### Filtering variants

To ensure the accuracy and relevance of genetic data, it's essential to filter variants and minimize false positives. This involves removing low-quality variants by applying specific thresholds for read depth and quality scores, for instance. In this study, we aimed to implement reasonable filtering criteria that control for false positives without being overly stringent, thus avoiding excluding too many true positives. We applied two filters: QUAL and DP. The QUAL filter corresponds to the quality score of the variant calling process, reflecting confidence in the detected variant at a given position. We set the QUAL threshold at 30, representing a 99.9% probability that the variant is correctly identified, thereby filtering out any variant calls with a lower score. The DP filter represents the total read depth at the variant's position, reflecting the volume of data supporting the variant. 

** Determining the DP threshold **

We established this threshold based on the coverage data from the 72 ECOR strains used to generate the pangenome and our understanding of the presence or absence of each locus of the pangenome in these strains.

We defined the DP threshold as the minimum number of reads required to confidently assert the presence of a nucleotide (and, by extension, the locus) in a strain. To calculate this threshold, we analyzed the presence-absence data for the 72 ECOR strains, incorporating the coverage depth observed for each nucleotide that mapped against the pangenome.

*1-Extracting DP information at each nucleotide*
First we needed to extract the Locus, Position, and depth information for each nucleotide from the mpileup files of each strains. 
To do so, we wrote a custom Python script `data_generation/scripts/extract_mpileup_info.py` that obtains this information for individual strains, and used it within a snakefile (`data_generation/scripts/sequencing_depth_info_extraction`) to process individual mpileups in parallel. 
Briefly, this Python script opens a mpileup file into a panda dataframe (skipping all the comments lines at the beginning of the mpileup), turns the DP values into numerics, and then extracts only the columns for CHROM (which corresponds to the locus name), POS (which corresponds to the position of the nucleotide in CHROM) and DP (which corresponds to the sequencing depth for this nucleotide). Finally, this is saved into a tab delimited file.

We ran the snakefile using the following command line:

`nohup snakemake -s sequencing_depth_info --cores 54 > extracting_DP.txt 2>&1 &`

*2-Calulating average read depth per locus in ECOR strains*

Then, using the custom  Python script `data_generation/scripts/numpy_merge_ecor.py`, we consolidated all the nucleotide-level DP information of the 72 ECOR strains into a single file `data_generation/results/ecor72_DP/ecor72_array.txt`.

`nohup python data_generation/scripts/numpy_merge_ecor.py > 2>&1 &`

We further used the custom R script `data_generation/scripts/Ecor72_averaging_locusDP.ipynb` to calculate the average read depth per locus for each strain as the sum of reads per nucleotide for a locus divided by the locus length.  

*3-Calulating the locus-level DP threshold associated with the presence/absence of a locus*
Finally, in the R notebook `data_generation/scripts/ECOR72_and_DP_threshold_analysis.Rmd` we assessed locus read depth patterns in relation to the loci's presence-absence status by integrating read depth and presence-absence data for each locus across all 72 ECOR strains. 
Based on the results we chose a DP threshold of 19.28 that indicates that any nucleotide or locus with a read depth exceeding 19.28 is confidently considered present.

** Filtering **

We filtered the variants using the chosen DP and QUAL thresholds and generated and indexed the VCF file with filtered variants (`data_generation/results/vcf/filtered_output.vcf.gz`)

`bcftools view -i 'QUAL >= 30 & DP > 21.81' data_generation/results/vcf/merged_output_all.vcf.gz’| bgzip -c > data_generation/results/vcf/filtered_output.vcf.gz`

and indexed the VCF file

`bcftools index data_generation/results/vcf/filtered_output.vcf.gz`

##### Variant annotation

We used SnpEff (https://doi.org/10.4161/fly.19695) to annotate variants within the pangenome's coding sequences. SnpEff analyzes input variants from a VCF file by annotating them based on a predefined database that includes genes, gene annotations, and gene sequence information. 

** Installing SnpEff **
First, we downloaded and installed SnpEff
`wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip`

`unzip snpEff_latest_core.zip`

** Annotating pangenome’s CDS sequences with Prokka **

To run SnpEff on our data, we first needed to create a custom database form the pangenome. To generate this databse, we needed a GFF file describing the genes or CDS in the pangenome. Thus, we used Prokka to annotate these CDS sequences (or genes), running:

`prokka --force --outdir data_generation/data/pangenome_cds  --prefix genes data_generation/data/pangenome_cds/pangenome_cds.fa`

This generated the gff file: `data_generation/results/pangenome_cds/genes.gff`

** Creating the custom database for the pangenome **

We generated our custom database from our pangenome following the directions detailed on the SnpEff website and edited the snpEff.config accordingly.

We created the custom database utilizing the GFF format annotation file and used the options `-noCheckCds -noCheckProtein` to bypass checks for transcript and protein sequence information.
To create the database, we ran:
` java -jar snpEff.jar build -gff3 -v ecor72 -noCheckCds -noCheckProtein`

And we checked if the database was created by running:
`java -jar snpEff.jar dump ecor72 | less`

** Annotating the variants **

Once the database has been created, we annotated the filtered VCF file by running:

` nohup sh -c 'java -Xmx32g -jar snpEff.jar eff -s annot_summary_filtered.html data_generation/results/vcf/filtered_output.vcf.gz | bgzip > data_generation/results/vcf/annotated_output.vcf.gz ' > nohup.out 2>&1 &`

** Filtering nonsilent variants **

We used SnpSift (https://doi.org/10.3389/fgene.2012.00035), a component of the SnpEff suite, to filter and refine annotated VCF files. We excluded silent mutations to retain only variants associated with missense, nonsense, and frameshift mutations, resulting in the curated file (`data_generation/results/vcf/output.non_silent.vcf.gz`).

We ran:
`java -jar SnpSift.jar filter "(ANN[*].EFFECT has '\''missense_variant'\'' | ANN[*].EFFECT has '\''nonsense_variant'\'' | ANN[*].EFFECT has '\''frameshift_variant'\'')" data_generation/results/vcf/annotated_output.vcf.gz | bgzip > data_generation/results/vcf/output.non_silent.vcf.gz`


#### Phenotype distribution analysis
To come

#### Variants distribution analysis
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
