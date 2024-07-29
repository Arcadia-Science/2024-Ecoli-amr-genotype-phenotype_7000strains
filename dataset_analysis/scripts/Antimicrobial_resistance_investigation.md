# Objectives

In this notebook, we focus on genes that have previously been shown to
be associated with antimicrobial resistance (AMR), and we evaluate the
variants identified within these genes in our cohort. We focus on five
genes known to confer resistance: -tetA\_1 (contig: LMHPMMMF\_04732),
tetA\_1 (contig: APHKLHJA\_00520), and tet\_3(contig: FCDKFLAE\_04147)
associated with tetracycline resistance -drfD (NGHFEPFE\_01999)
associated with trimethoprim resistance -catA1 (DHJNCGMO\_04398)
associated with chloramphenicol resistance

After identifying the variants within each genes and their different
genotypes in the strains, we correlate that information with the AMR
phenotype for the right antibiotic in the corresponding strains, when
AMR phenotype information was available.

------------------------------------------------------------------------

# Environment setup

    rm(list=ls())

    library(tidyverse)

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

    library(VariantAnnotation)

    ## Loading required package: BiocGenerics
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min
    ## 
    ## Loading required package: MatrixGenerics
    ## Loading required package: matrixStats
    ## 
    ## Attaching package: 'matrixStats'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count
    ## 
    ## 
    ## Attaching package: 'MatrixGenerics'
    ## 
    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars
    ## 
    ## Loading required package: GenomeInfoDb
    ## Loading required package: S4Vectors
    ## Loading required package: stats4
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## Loading required package: GenomicRanges
    ## Loading required package: SummarizedExperiment
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## 
    ## Attaching package: 'Biobase'
    ## 
    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians
    ## 
    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians
    ## 
    ## Loading required package: Rsamtools
    ## Loading required package: Biostrings
    ## Loading required package: XVector
    ## 
    ## Attaching package: 'XVector'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact
    ## 
    ## 
    ## Attaching package: 'Biostrings'
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit
    ## 
    ## 
    ## Attaching package: 'VariantAnnotation'
    ## 
    ## The following object is masked from 'package:stringr':
    ## 
    ##     fixed
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     tabulate

    library(ComplexHeatmap)

    ## Loading required package: grid
    ## 
    ## Attaching package: 'grid'
    ## 
    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     pattern
    ## 
    ## ========================================
    ## ComplexHeatmap version 2.14.0
    ## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    ## Github page: https://github.com/jokergoo/ComplexHeatmap
    ## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    ## 
    ## If you use it in published research, please cite either one:
    ## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
    ## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    ##     genomic data. Bioinformatics 2016.
    ## 
    ## 
    ## The new InteractiveComplexHeatmap package can directly export static 
    ## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(ComplexHeatmap))
    ## ========================================

------------------------------------------------------------------------

# Importing and formating the data necessary for the analysis

## Importing the variant genotype matrix

In this section, we import the variant calling format (VCF) file
containing variant information and associated genotypes for the five
genes we decided to investigate.

    # Importing the vcf file 

    vcf_file_ab=readVcf('dataset_analysis/data/antimicrobial_resistance_analysis/resistance_output.non_silent.vcf.gz', genome="ecoli")
    geno_matrix_ab <- as.data.frame(geno(vcf_file_ab)$GT) # turn the vcf into a matrix with variants as rows and strains as column and genotype as information (1, means that the strain possess the alternative allele , '.' means that there is no called genotype. Note a '.' genotype does not necessarily means that the strain posses the ancestral allele, it could also mean that this strain does not posess that gene)


    # Formating the data for downstream process
     
      ## Keeping only strains were there is a called variants (=strains with known alternative allele)
    ab_data <- geno_matrix_ab %>%
      dplyr::select(where(~ !all(. == '.')))

     ## Adding the contig name as the row.name
    ab_data$LOC=rownames(ab_data)

     ## Adding a column with row.name and formatting to add contigs and pos information 
    ab_data <- ab_data %>%
      dplyr::select(last_col(), everything())
    ab_data <- separate(ab_data, LOC, into = c("contigs", "PosB"), sep = ":") 
    ab_data <-  separate(ab_data, PosB, into = c("pos", "allele"), sep = "_")
    ab_data$pos=as.character(ab_data$pos)

    head(ab_data[,1:10],10)

    ##                              contigs pos allele ERR1197946 ERR1197966
    ## LMHECDEF_03475:10_T/G LMHECDEF_03475  10    T/G          .          .
    ## LMHECDEF_03475:11_C/T LMHECDEF_03475  11    C/T          .          .
    ## LMHECDEF_03475:20_C/T LMHECDEF_03475  20    C/T          .          .
    ## LMHECDEF_03475:22_C/A LMHECDEF_03475  22    C/A          .          .
    ## LMHECDEF_03475:23_A/G LMHECDEF_03475  23    A/G          .          .
    ## LMHECDEF_03475:28_T/G LMHECDEF_03475  28    T/G          .          .
    ## LMHECDEF_03475:29_C/A LMHECDEF_03475  29    C/A          .          .
    ## LMHECDEF_03475:34_T/A LMHECDEF_03475  34    T/A          .          .
    ## LMHECDEF_03475:35_C/T LMHECDEF_03475  35    C/T          .          .
    ## LMHECDEF_03475:46_G/A LMHECDEF_03475  46    G/A          .          .
    ##                       ERR1197968 ERR1218535 ERR1218536 ERR1218538 ERR1218539
    ## LMHECDEF_03475:10_T/G          .          .          .          .          .
    ## LMHECDEF_03475:11_C/T          .          .          .          .          .
    ## LMHECDEF_03475:20_C/T          .          .          .          .          .
    ## LMHECDEF_03475:22_C/A          .          .          .          .          .
    ## LMHECDEF_03475:23_A/G          .          .          .          .          .
    ## LMHECDEF_03475:28_T/G          .          .          .          .          .
    ## LMHECDEF_03475:29_C/A          .          .          .          .          .
    ## LMHECDEF_03475:34_T/A          .          .          .          .          .
    ## LMHECDEF_03475:35_C/T          .          .          .          .          .
    ## LMHECDEF_03475:46_G/A          .          .          .          .          .

## Importing AMR contig variants general information

Next, for each variant, we add information regarding - the number of
allele - the number of strains with the alternative allele - the generic
name of the AMR gene (catA1, tetA\_3, etc…) - the antibiotic the variant
is expected to confer resistance to

    # Importing information

    res_freq=read.csv('dataset_analysis/data/antimicrobial_resistance_analysis/antibiotic_resistance_freq.csv') # file with all required information
    res_freq$pos=as.character(res_freq$pos)
    colnames(res_freq)=c('contigs','pos','n_alleles','n_strains','Type','Resistance_gene','Antibiotic')

    # Adding information to the variant genotype matrix
    ab_data=left_join(ab_data, res_freq, by=c('contigs','pos'))

    ab_data <- ab_data %>%
      filter(!is.na(Resistance_gene))

    ## Reformatting data to obtain by strain information 

    ab_data_rfmt=ab_data %>% gather(Strain,Called,ERR1197946:SRR850839)
    ab_data_rfmt=subset(ab_data_rfmt, ab_data_rfmt$Called!='.') #only keep entries when a alrernative allel was called

    head(ab_data_rfmt,10)

    ##             contigs pos allele n_alleles n_strains  Type Resistance_gene
    ## 351  FCDKFLAE_04147 572    G/A         2       110 LOCUS        tet(A)_3
    ## 466  FCDKFLAE_04147 482    A/G         2        14 LOCUS        tet(A)_3
    ## 729  FCDKFLAE_04147 572    G/A         2       110 LOCUS        tet(A)_3
    ## 787  FCDKFLAE_04147 548    C/A         2         3 LOCUS        tet(A)_3
    ## 792  FCDKFLAE_04147 572    G/A         2       110 LOCUS        tet(A)_3
    ## 841  FCDKFLAE_04147 440    G/A         2         3 LOCUS        tet(A)_3
    ## 918  FCDKFLAE_04147 572    G/A         2       110 LOCUS        tet(A)_3
    ## 981  FCDKFLAE_04147 572    G/A         2       110 LOCUS        tet(A)_3
    ## 1079 LMHPMMMF_04732 800    G/T         2         1 LOCUS        tet(A)_1
    ## 1107 FCDKFLAE_04147 572    G/A         2       110 LOCUS        tet(A)_3
    ##        Antibiotic     Strain Called
    ## 351  Tetracycline ERR1218538      1
    ## 466  Tetracycline ERR1218542      1
    ## 729  Tetracycline ERR1218546      1
    ## 787  Tetracycline ERR1218549      1
    ## 792  Tetracycline ERR1218549      1
    ## 841  Tetracycline ERR1218551      1
    ## 918  Tetracycline ERR1218552      1
    ## 981  Tetracycline ERR1218553      1
    ## 1079 Tetracycline ERR1218558      1
    ## 1107 Tetracycline ERR1218558      1

## Importing the phenotype matrix

Next we use the original phenotype matrix, to extract AMR information
for strains that have been identified with variants in the 5 contigs of
interest

In addition, we import a file that gives the correspondence between the
strain genome identification numbers and the strain genome SRA accession
numbers to be able to connect information from the phenotype matrix and
the genotype matrix.

    #Importing AMR phenotype map
    amr_phen=read.csv("dataset_analysis/results/phenotype_matrix_08302024.csv")[,-1]

    #Importing correspondence SRA and genome name
    SRA_to_genome=read.csv('dataset_analysis/data/antimicrobial_resistance_analysis/SRA_to_genome_name.csv')[,-1]

     ##Changing Genome name to SRA accession in AMR phenotype matrix

    amr_phen_sra=left_join(SRA_to_genome, amr_phen, by="Genome.Name")

     ##Extracting AMR phenotypes for species with variants in AMR known genes

    amr_phen_geneamr=subset(amr_phen_sra, amr_phen_sra$SRA.Accession%in%ab_data_rfmt$Strain)

------------------------------------------------------------------------

# Phenotypes correlation

The goal of the ext section is to evaluate the correlation between AMR
phenotype for an antibiotic and the presence/absence of genetic variant
in associated resistance gene. Here, we focus on three antibiotics:
chloramphenicol, trimethoprim and tetracycline. For each antibiotics, we
have identified genetic variants and associated strains. Wgen the AMR
phenotype is available for these strain-antibiotic combinations, we can
evaluate whether these variations are associated with antimicrobial
resistance.

## Chloramphenicol

    # Subsetting variants within chloramphenicol resistance gene (catA1)
    chlor_var=subset(ab_data_rfmt, ab_data_rfmt$Antibiotic=='Chloramphenicol')


    # Selection of the strains with available Chloramphenicol AMR phenotypes

    chlor_amr_all=amr_phen_sra %>% dplyr::select(Genome.Name,SRA.Accession, chloramphenicol) %>%
      na.omit()

    chlor_amr_res=subset(chlor_amr_all, chlor_amr_all$chloramphenicol=='Resistant')

    # Subset the phenotype information for the strains that have been identified with variant within the chloramphenicol resistance gene 
    chlor_amr=subset(chlor_amr_all, chlor_amr_all$SRA.Accession%in%chlor_var$Strain) 

    chlor_amr

    ##                         Genome.Name SRA.Accession chloramphenicol
    ## 2724 Escherichia coli strain PH92-1    SRR1610061       Resistant
    ## 2726   Escherichia coli strain PH93    SRR1610062       Resistant

    # Testing if the set of variants with documented AMR phenotype is enriched in Resistant phenotypes

    N=dim(chlor_amr_all)[1] #number of strains with known AMR for Chloramphenicol
    K=dim(chlor_amr_res)[1] #number of strains with known AMR for chloramphenicol and that are associated with a Resistant phenotype
    n=dim(chlor_amr)[1] #number of strains with knwon AMR for chloramphenicol and that have a variant in catA1
    k=dim(subset(chlor_amr,chlor_amr$chloramphenicol=='Resistant'))[1] #number of strains with known AMR and variant in catA1 and associated with Resistant phenotype

    p_value <- dhyper(k, K, N - K, n)

    p_value

    ## [1] 0.07892543

## Trimethoprim

    # Subseting variants within the Trimethoprim resistance gene drfD

    trim_var=subset(ab_data_rfmt, ab_data_rfmt$Antibiotic=='Trimethoprim') 


    # Selection of strains with known Trimethoprim AMR phenotypes
    trim_amr_all=amr_phen_sra %>% dplyr::select(Genome.Name,SRA.Accession, trimethoprim) %>%
      na.omit()

    trim_amr_res=subset(trim_amr_all, trim_amr_all$trimethoprim=='Resistant')

    # Subset the phenotype information for the strains that have been identified with variant within the trimethoprim resistance gene
    trim_amr=subset(trim_amr_all, trim_amr_all$SRA.Accession%in%trim_var$Strain)

    trim_amr

    ##                                    Genome.Name SRA.Accession trimethoprim
    ## 953     Escherichia coli O25:H4 strain ECO0056     ERR435030    Resistant
    ## 972      Escherichia coli NA:H6 strain ECO0076     ERR435033  Susceptible
    ## 982     Escherichia coli O25:H4 strain ECO0087     ERR434649  Susceptible
    ## 1003    Escherichia coli O25:H4 strain ECO0111     ERR434764    Resistant
    ## 1061    Escherichia coli O25:H4 strain ECO0177     ERR434779    Resistant
    ## 1065   Escherichia coli O83:H31 strain ECO0181     ERR435045    Resistant
    ## 1079    Escherichia coli O25:H4 strain ECO0207     ERR439636  Susceptible
    ## 1110     Escherichia coli O9:H4 strain ECO0240     ERR434460  Susceptible
    ## 1137    Escherichia coli O25:H4 strain ECO0270     ERR435107    Resistant
    ## 1170   Escherichia coli O83:H33 strain ECO0305     ERR435047  Susceptible
    ## 1231    Escherichia coli O18:H7 strain ECO0390     ERR435492  Susceptible
    ## 1246  Escherichia coli O107:H27 strain ECO0408     ERR434560  Susceptible
    ## 1309 Escherichia coli not known strain ECO0321     ERR434970  Susceptible
    ## 2273    Escherichia coli NA:H18 strain ECO0379     ERR434636    Resistant
    ## 2576             Escherichia coli strain PB339    ERR1218617    Resistant
    ## 2603       Escherichia coli strain W2_10_ERB10    ERR1218681  Susceptible
    ## 2614        Escherichia coli strain W2_12_ERB1    ERR1218685  Susceptible
    ## 2619         Escherichia coli strain W2_4_ERB4    ERR1218661    Resistant
    ## 2622         Escherichia coli strain W2_3_ERB7    ERR1218658    Resistant
    ## 2624        Escherichia coli strain W2_3_ERG14    ERR1218660    Resistant
    ## 2632         Escherichia coli strain W2_8_ERB3    ERR1218561  Susceptible
    ## 2756            Escherichia coli strain ENV737    ERR2580172  Susceptible
    ## 3723      Escherichia coli strain 503028_aEPEC     ERR134544    Resistant
    ## 3726      Escherichia coli strain 504005_aEPEC     ERR137803    Resistant
    ## 3735      Escherichia coli strain 202453_aEPEC     ERR134525    Resistant
    ## 3737        Escherichia coli strain G100788-1A     ERR175731  Susceptible
    ## 3751           Escherichia coli strain G503854     ERR178224  Susceptible
    ## 3767            Escherichia coli strain 202374     ERR178152  Susceptible
    ## 3770            Escherichia coli strain 403728     ERR178161    Resistant
    ## 3776            Escherichia coli strain 201350     ERR178216  Susceptible
    ## 3777      Escherichia coli strain 202443_aEPEC     ERR134524    Resistant
    ## 3809            Escherichia coli strain 504528     ERR178168    Resistant
    ## 3813            Escherichia coli strain 702328     ERR178174  Susceptible
    ## 3819      Escherichia coli strain 202423_aEPEC     ERR134523    Resistant
    ## 3823            Escherichia coli strain 503256     ERR178197    Resistant
    ## 3851      Escherichia coli strain 102328_aEPEC     ERR134517  Susceptible
    ## 3855            Escherichia coli strain 202474     ERR178153  Susceptible
    ## 3863           Escherichia coli strain G603423     ERR178228    Resistant
    ## 3877            Escherichia coli strain 500193     ERR178213    Resistant
    ## 3880           Escherichia coli strain G303212     ERR175730    Resistant

    # Testing if the set of variants with documented AMR phenotype is enriched in Resistant phenotypes

    N=dim(trim_amr_all)[1] #number of strains with known AMR for trimethoprim
    K=dim(trim_amr_res)[1] #number of strains with known AMR for trimethoprim and that are associated with a Resistant phenotype
    n=dim(trim_amr)[1] #number of strains with knwon AMR for trimethoprim and that have a variant in catA1
    k=dim(subset(trim_amr,trim_amr$trimethoprim=='Resistant'))[1] #number of strains with known AMR and variant in catA1 and associated with Resistant phenotype

    p_value <- dhyper(k, K, N - K, n)

    p_value

    ## [1] 0.1123763

Because the 40 strains for whcich we have variant and phenotype
information are split between Susceptible and Resistant phenotypes, we
investigated if specific positions within *dfrD* are associated with
Resistance

    trim_amr_phen=na.omit(trim_amr)
    trim_var_phen=subset(trim_var, trim_var$Strain%in%trim_amr_phen$SRA.Accession)
    trim_var_phen= trim_var_phen %>% dplyr::select(pos, Strain, Called) %>%
      spread(Strain, Called)

    trim_var_phen$pos=as.numeric(trim_var_phen$pos)

    df_hm=data.frame(pos=trim_var_phen$pos)

    ERR=c(colnames(trim_var_phen[,-1]))
    phen_loc=colnames(trim_var_phen)

    for (i in 1:length(ERR)){
      ERR_tmp=ERR[i]
      phen=subset(trim_amr_phen, trim_amr_phen$SRA.Accession==ERR_tmp)
      phen=phen[,3]
      
      vec_name=c(trim_var_phen$pos)
      loc=which(phen_loc==ERR_tmp)
      vec_var=trim_var_phen[,loc]
      names(vec_var)= vec_name
      vec_var[vec_var == 1] <- phen
      
      df_tmp=data.frame('Sample'=vec_var)
      rownames(df_tmp)=vec_name
      colnames(df_tmp)=ERR_tmp
      
      df_hm=cbind(df_hm,df_tmp)
      
    }

    df_hm <- df_hm[order(df_hm$pos), ]
    mat=as.matrix(t(df_hm[,-1]))
    mat[is.na(mat)] <- "NA"

    colors <- c("Susceptible" = "#7A77AB", "Resistant" = "#F7B846", "No_phenotype_available"="#EDE6DA")
    names(colors) <- c("Susceptible", "Resistant", "NA")

    Heatmap(mat, col = colors, name = "Categories",rect_gp = gpar(col = "white", lwd = 1),
            row_names_gp = gpar(fontsize = 7)) 

![](Antimicrobial_resistance_investigation_files/figure-markdown_strict/unnamed-chunk-7-1.png)

# Tetracycline

    # Subseting variants within Trimethoprim resistance genes (tet genes)
    tetra_var=subset(ab_data_rfmt, ab_data_rfmt$Antibiotic=='Tetracycline')

    # Selection of strains with known Tetracycline AMR phenotypes
    tetra_amr_all=amr_phen_sra %>% dplyr::select(Genome.Name,SRA.Accession, tetracycline)%>%
      na.omit()

    tetra_amr_res=subset(tetra_amr_all, tetra_amr_all$tetracycline=='Resistant')

    # Subset the phenotype information for the strains that have been identified with variant within the tetracycline resistance genes

    tetra_amr=subset(tetra_amr_all, tetra_amr_all$SRA.Accession%in%tetra_var$Strain)

    tetra_amr

    ##                               Genome.Name SRA.Accession tetracycline
    ## 2376          Escherichia coli strain 197    SRR4245481    Resistant
    ## 2528    Escherichia coli strain MUGSI_162    SRR4065703    Resistant
    ## 2536    Escherichia coli strain MUGSI_217    SRR4065676    Resistant
    ## 2546    Escherichia coli strain MUGSI_250    SRR4017901    Resistant
    ## 2733        Escherichia coli strain ENV49    ERR2580149    Resistant
    ## 2734        Escherichia coli strain ENV66    ERR2580150    Resistant
    ## 2736        Escherichia coli strain ENV76    ERR2580152    Resistant
    ## 2737       Escherichia coli strain ENV103    ERR2580153    Resistant
    ## 2741       Escherichia coli strain ENV225    ERR2580157    Resistant
    ## 2747       Escherichia coli strain ENV323    ERR2580163    Resistant
    ## 2748       Escherichia coli strain ENV326    ERR2580164    Resistant
    ## 2754       Escherichia coli strain ENV727    ERR2580171    Resistant
    ## 3678      Escherichia coli strain CFS0351    ERR3449882    Resistant
    ## 3722 Escherichia coli strain 401886_aEPEC     ERR137792    Resistant
    ## 3727 Escherichia coli strain 200708_aEPEC     ERR137782    Resistant
    ## 3734 Escherichia coli strain 401938_aEPEC     ERR137793    Resistant
    ## 3746 Escherichia coli strain 201589_aEPEC     ERR137816    Resistant
    ## 3749 Escherichia coli strain 201488_aEPEC     ERR137815    Resistant
    ## 3756 Escherichia coli strain 204302_aEPEC     ERR134528    Resistant
    ## 3771       Escherichia coli strain 402794     ERR178204    Resistant
    ## 3776       Escherichia coli strain 201350     ERR178216    Resistant
    ## 3815 Escherichia coli strain 200135_aEPEC     ERR134519    Resistant
    ## 3827 Escherichia coli strain 100383_aEPEC     ERR137807    Resistant
    ## 3838 Escherichia coli strain 402099_aEPEC     ERR134533    Resistant
    ## 3858 Escherichia coli strain 102366_aEPEC     ERR137810    Resistant
    ## 3865 Escherichia coli strain 402248_aEPEC     ERR134534    Resistant

    # Testing if the set of variants with documented AMR phenotype is enriched in Resistant phenotypes

    N=dim(tetra_amr_all)[1] #number of strains with known AMR for tetracycline
    K=dim(tetra_amr_res)[1] #number of strains with known AMR for tetracycline and that are associated with a Resistant phenotype
    n=dim(tetra_amr)[1] #number of strains with knwon AMR for tetracycline and that have a variant in catA1
    k=dim(subset(tetra_amr,tetra_amr$tetracycline=='Resistant'))[1] #number of strains with known AMR and variant in catA1 and associated with Resistant phenotype

    p_value <- dhyper(k, K, N - K, n)
    p_value

    ## [1] 1.092171e-06

------------------------------------------------------------------------

    sessionInfo()

    ## R version 4.2.3 (2023-03-15)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Monterey 12.5
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ComplexHeatmap_2.14.0       VariantAnnotation_1.44.1   
    ##  [3] Rsamtools_2.14.0            Biostrings_2.66.0          
    ##  [5] XVector_0.38.0              SummarizedExperiment_1.28.0
    ##  [7] Biobase_2.58.0              GenomicRanges_1.50.2       
    ##  [9] GenomeInfoDb_1.34.9         IRanges_2.32.0             
    ## [11] S4Vectors_0.36.2            MatrixGenerics_1.10.0      
    ## [13] matrixStats_1.2.0           BiocGenerics_0.44.0        
    ## [15] lubridate_1.9.3             forcats_1.0.0              
    ## [17] stringr_1.5.1               dplyr_1.1.4                
    ## [19] purrr_1.0.2                 readr_2.1.5                
    ## [21] tidyr_1.3.1                 tibble_3.2.1               
    ## [23] ggplot2_3.5.1               tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7             bit64_4.0.5              RColorBrewer_1.1-3      
    ##  [4] doParallel_1.0.17        filelock_1.0.3           progress_1.2.3          
    ##  [7] httr_1.4.7               tools_4.2.3              utf8_1.2.4              
    ## [10] R6_2.5.1                 DBI_1.2.2                colorspace_2.1-0        
    ## [13] GetoptLong_1.0.5         withr_3.0.0              tidyselect_1.2.1        
    ## [16] prettyunits_1.2.0        bit_4.0.5                curl_5.2.1              
    ## [19] compiler_4.2.3           cli_3.6.3                xml2_1.3.6              
    ## [22] DelayedArray_0.24.0      rtracklayer_1.58.0       scales_1.3.0            
    ## [25] rappdirs_0.3.3           digest_0.6.35            rmarkdown_2.26          
    ## [28] pkgconfig_2.0.3          htmltools_0.5.8.1        highr_0.10              
    ## [31] BSgenome_1.66.3          dbplyr_2.5.0             fastmap_1.1.1           
    ## [34] GlobalOptions_0.1.2      rlang_1.1.4              rstudioapi_0.16.0       
    ## [37] RSQLite_2.3.6            shape_1.4.6.1            BiocIO_1.8.0            
    ## [40] generics_0.1.3           BiocParallel_1.32.6      RCurl_1.98-1.14         
    ## [43] magrittr_2.0.3           GenomeInfoDbData_1.2.9   Matrix_1.6-5            
    ## [46] Rcpp_1.0.12              munsell_0.5.1            fansi_1.0.6             
    ## [49] lifecycle_1.0.4          stringi_1.8.3            yaml_2.3.8              
    ## [52] zlibbioc_1.44.0          BiocFileCache_2.6.1      blob_1.2.4              
    ## [55] parallel_4.2.3           crayon_1.5.2             lattice_0.22-6          
    ## [58] GenomicFeatures_1.50.4   circlize_0.4.16          hms_1.1.3               
    ## [61] KEGGREST_1.38.0          knitr_1.46               pillar_1.9.0            
    ## [64] rjson_0.2.21             codetools_0.2-20         biomaRt_2.54.1          
    ## [67] XML_3.99-0.16.1          glue_1.7.0               evaluate_0.23           
    ## [70] foreach_1.5.2            png_0.1-8                vctrs_0.6.5             
    ## [73] tzdb_0.4.0               gtable_0.3.5             clue_0.3-65             
    ## [76] cachem_1.0.8             xfun_0.43                restfulr_0.0.15         
    ## [79] iterators_1.0.14         GenomicAlignments_1.34.1 AnnotationDbi_1.60.2    
    ## [82] memoise_2.0.1            cluster_2.1.6            timechange_0.3.0
