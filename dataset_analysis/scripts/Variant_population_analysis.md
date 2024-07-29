# Objectives

In this notebook, we describe the output of the variant calling
investigation in the 7000+ strains cohort.

More specifically, we investigate: - how many variants are found in
coding sequences (CDS) and how many variants are found in intergenic
regions (IGR) - how many variants are observed per contig and IGR
(variant rate) - how many alleles are observed per variants - for each
variant: how many strains have the alternative variant(s)

Then, we focus on variants in CDS that are annotated as non-silent
variants and investigate: - number of strains with alternative allele
for non silent variants - the proportion of non-silent variant per
contig - functional annotation of the CDS associated with high and low
non-silent varaiations rates

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

    library(seqinr)

    ## 
    ## Attaching package: 'seqinr'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    library(ggridges)
    library(ape)

    ## 
    ## Attaching package: 'ape'
    ## 
    ## The following objects are masked from 'package:seqinr':
    ## 
    ##     as.alignment, consensus
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

    library(arcadiathemeR)

------------------------------------------------------------------------

# Investigation of the variant distribution

We start by analyzing the distribution of identified variants in the
whole cohort. At this stage, no variant annotation has happened and
there is no distinction between silent and non-silent variants. However,
variants have been filtered based on their reliability using DP and QUAL
filters (see Pub Approch section
(<https://doi.org/10.57844/arcadia-d2cf-ebe5>) or Github repo ReadMe).

## Data import

For this part of the analysis, we focus on whole variants identified in
the cohort regardless of the strains. The per strain genotype is not
necessary. The only information we need to characterize the variant
population is, for each identified variant, it contig name and its
position in this contig.

To perform the analysis, we need three input files: - the list of
variants characterized by their contig and position - the pangenome
fasta file to obtain contigs length - the variant ‘frequency’ file that
indicates for each variant the number of strains the associated
alternative allele(s) is(are) found in

    # Variant information (this corresponds to contig and position within contig information of variants) 

    variants=read.delim('dataset_analysis/data/variant_analysis/variants_pos.tsv',header=T, sep="\t")
    variants=rbind(names(variants),variants) # reformatting the headers
    variants[1,2]=12
    colnames(variants)=c('contigs','pos')

    # Pangenome fasta file

    pan_sequences <- read.fasta(file = "dataset_generation/results/pangenome_whole/whole_pangenome.fasta", seqtype = "DNA")

    #Variant 'frequency' information (provide information of the of alleles per variant and the number of strains where alternative allele(s) was(were) found)

    freq_info=read.delim('dataset_analysis/data/variant_analysis/allele_freqs.txt', header = TRUE, sep = "\t", row.names = NULL)
    freq_info_cleaned= freq_info %>% select(row.names,CHROM, POS,N_ALLELES) # clean up the data
    colnames(freq_info_cleaned)=c('contigs','pos','n_alleles','n_strains')

## Identifying CDS versus IGR

When we constructed the pangenome, we included both CDS and intergenic
regions (IGR), thus there are two types of contigs in the pangenome: CDS
and IGR. We first want to investigate how variants are distributed
within each contig type

IGR are identified by their name that starts with ‘Cluster\_’. Relying
on this specificity, we add a new column to tag CDS and IGR

    variants <- variants %>%
      mutate(Type = ifelse(grepl("^Cluster", contigs), "IGR", "CDS"))


    # We count how many entries for each Type
    table_loc_vs_igr=data.frame(table(variants$Type))

    table_loc_vs_igr

    ##   Var1    Freq
    ## 1  CDS 2086374
    ## 2  IGR  364810

## Calculation of variant rate per CDS and IGR

For each CDS or IGR we identify how many different positions are
associated with genetic variations. This can be simply done by looking
at how many time each cds name or igr name is seen in the table.

    # Calculation of the number of variants per contig

    table_variant_per_contig=data.frame(table(variants$contigs))

    colnames(table_variant_per_contig)=c('contigs','numb_var')

    table_variant_per_contig <- table_variant_per_contig %>%
      mutate(Type = ifelse(grepl("^Cluster", contigs), "IGR", "CDS"))

    head(table_variant_per_contig,10)

    ##           contigs numb_var Type
    ## 1  ABCDEHGJ_00151       36  CDS
    ## 2  ABCDEHGJ_00444       22  CDS
    ## 3  ABCDEHGJ_00553       31  CDS
    ## 4  ABCDEHGJ_00559        5  CDS
    ## 5  ABCDEHGJ_00562        7  CDS
    ## 6  ABCDEHGJ_00564       35  CDS
    ## 7  ABCDEHGJ_00844       50  CDS
    ## 8  ABCDEHGJ_01604       24  CDS
    ## 9  ABCDEHGJ_01669       50  CDS
    ## 10 ABCDEHGJ_01671       21  CDS

To compare the distribution of the number of variants per cds and per
igr, we calculate variant rate that corresponds to the contig length (in
bp) divided by the number of variants in this contig.

    # We use the pangenome fasta file to obtain the size of each contig

     ## Extract sequence names and their lengths
    seq_names <- sapply(pan_sequences, attr, "name")
    seq_lengths <- sapply(pan_sequences, length)

     ## Create a dataframe
    contig_length = data.frame(contigs = seq_names, len = seq_lengths)

     ## Combine contig length information with the variant information
    contig_length_var=subset(contig_length, contig_length$contigs%in%table_variant_per_contig$contigs)

    # Generate the dataframe with variant rate information

    table_nbvariant_rate=left_join(table_variant_per_contig, contig_length_var, by='contigs')

    table_nbvariant_rate$rate_var=table_nbvariant_rate$len/table_nbvariant_rate$numb_var

    # Plot the distribution of variant rate per cds or igr

    viop_nbvar=ggplot(table_nbvariant_rate, aes(x=Type, y=rate_var, fill=Type))+
      geom_violin(alpha=0.5)+theme_arcadia(x_axis_type = "categorical")+ ylab('variant rate') +
    scale_fill_arcadia(palette_name = "primary") + ylim(0,1000) +
     geom_boxplot(width=0.1, outlier.shape = NA) 

    viop_nbvar

    ## Warning: Removed 63 rows containing non-finite outside the scale range
    ## (`stat_ydensity()`).

    ## Warning: Removed 63 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](Variant_population_analysis_files/figure-markdown_strict/unnamed-chunk-5-1.png)

**Main metrics of varianr rate for CDS contigs**

    CDS_rate=subset(table_nbvariant_rate, table_nbvariant_rate$Type=='CDS')
    summary(CDS_rate$rate_var)

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##    1.889    4.226    7.925   37.493   18.500 2250.000

**Main metrics of varianr rate for IGR contigs**

    IGR_rate=subset(table_nbvariant_rate, table_nbvariant_rate$Type=='IGR')
    summary(IGR_rate$rate_var)

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   1.636   5.444   9.833  23.591  21.500 832.000

## Variant ‘frequency’

For each variant we evaluated the number of strains in which the
alternative allele(s) was(were) found. Although this isn’t the precise
measurement of allele frequency, which compares the presence of
alternative alleles to the ancestral allele, it still provides insights
into the commonality of variations.

    # Calculation of the number of variants per contig

      ## Reformating the dataframe

    freq_info_cleaned=subset(freq_info_cleaned,freq_info_cleaned$contigs%in%contig_length_var$contigs)

    freq_info_cleaned <- freq_info_cleaned %>%
      mutate(Type = ifelse(grepl("^Cluster", contigs), "IGR", "CDS"))

    freq_info_cleaned$n_alleles=as.numeric(freq_info_cleaned$n_alleles)
    freq_info_cleaned$n_strains=as.numeric(freq_info_cleaned$n_strains)

     ## Summarizing the number of variants and strains with alternative allele 
    table_strain=data.frame(table(freq_info_cleaned$n_strains))
    table_strain$Var1=as.numeric(table_strain$Var1)
    colnames(table_strain)=c('nb_strain','nb_variants_with_that_nb_strains_with_alt_allele')

    head(table_strain,10)

    ##    nb_strain nb_variants_with_that_nb_strains_with_alt_allele
    ## 1          1                                           527686
    ## 2          2                                           233987
    ## 3          3                                           147625
    ## 4          4                                           170274
    ## 5          5                                           104962
    ## 6          6                                            77557
    ## 7          7                                            56115
    ## 8          8                                            46140
    ## 9          9                                            40417
    ## 10        10                                            35547

------------------------------------------------------------------------

# Investigation of the non-silent variants

In this section we focus on the variants that have been annotated as
non-silent in CDS contigs. After annotating the variants using SnpEff
(<https://doi.org/10.4161/fly.19695>), we kept non-silent variants
defined as variants associated with missense, nonsense, and frame-shift
mutations.

## Data import

As previously, we only need the variant contig and position information,
and the per strain genotype is not necessary. The data we import for
this part are: - the non-silent variant contig and position
information - the variant ‘frequency’ information defined as, for each
variant, the number of strains where the alternative allele(s) is(are)
found.

    #Variant information

    variants_ns=read.delim('dataset_analysis/data/variant_analysis/variants_non_silent_pos.tsv',header=T, sep="\t")
      
    variants_ns=rbind(names(variants_ns),variants_ns) #change the header as the first row
    variants_ns[1,2]=14
    colnames(variants_ns)=c('contigs','pos')


    #Variant 'frequency' information (provide information of # of allele and number of strains were alternative allele(s) was(were) found)

    freq_info_ns=read.delim('dataset_analysis/data/variant_analysis/allele_non_silent_freqs.txt', header = TRUE, sep = "\t", row.names = NULL)
    freq_info_ns_cleaned= freq_info_ns %>% select(row.names,CHROM, POS,N_ALLELES) # clean up the data
    colnames(freq_info_ns_cleaned)=c('contigs','pos','n_alleles','n_strains')

## Non-silent variant ‘frequency’

We start by investigating the non-silent variant frequency

    # Reformatting frequency information
    freq_info_ns_cleaned=subset(freq_info_ns_cleaned,freq_info_ns_cleaned$contigs%in%contig_length_var$contigs)
    freq_info_ns_cleaned$n_alleles=as.numeric(freq_info_ns_cleaned$n_alleles)
    freq_info_ns_cleaned$n_strains=as.numeric(freq_info_ns_cleaned$n_strains)

    # Summarizing the number of variants and strains with alternative allele 
    table_ns_strain=data.frame(table(freq_info_ns_cleaned$n_strains))
    table_ns_strain$Var1=as.numeric(table_ns_strain$Var1)
    colnames(table_ns_strain)=c('nb_strain','nb_variants_with_that_nb_strains_with_alt_allele')

    head(table_ns_strain,10)

    ##    nb_strain nb_variants_with_that_nb_strains_with_alt_allele
    ## 1          1                                           250832
    ## 2          2                                            86968
    ## 3          3                                            47619
    ## 4          4                                            41377
    ## 5          5                                            26869
    ## 6          6                                            21202
    ## 7          7                                            14470
    ## 8          8                                            12064
    ## 9          9                                            11034
    ## 10        10                                             9711

## Investigation of non-silent variants and fraction of non-silent variants per contig

We next investigate what fraction of variant in each contig is actually
associated with non-silent variations. To do so, for each CDS, we
calculate the ration of the number of non-silent variant in a contig to
the total number of variants in this contig.  
This information will highlight cds that still have high or low rates of
non-silent variations. CDS dominated by non-silent variants may suggest
positive selection for variants conferring advantageous traits. In
contrast, cds predominantly containing silent variants likely highlight
essential, conserved cellular functions where non-silent mutations would
be deleterious and selected against.

    # Calculating the number of non-silent variant per annotated CDS
    numb_cds=as.data.frame(table(variants_ns$contigs))
    colnames(numb_cds)=c('contigs','numb_var_ns')

    # Subset original variant list with CDS were non silent variants were found
    CDS_rate_sub=subset(CDS_rate, CDS_rate$contigs%in%numb_cds$contigs)
    colnames(CDS_rate_sub)=c('contigs','numb_var_all','Type','len','rate_var_all')

    # Combine into a single dataframe with non-silent variant/per contig information and all variant information/contig

    CDS_rate_comp=left_join(numb_cds,CDS_rate_sub, by='contigs')

    # Calculate the fraction of non-silent variant to all variants at each contig 
    CDS_rate_comp$ratio=CDS_rate_comp$numb_var_ns/CDS_rate_comp$numb_var_all  

    dens_ratio=ggplot(CDS_rate_comp, aes(x=ratio))+
      geom_density(fill='#BABEE0')+theme_arcadia(x_axis_type = 'numerical')

    dens_ratio

![](Variant_population_analysis_files/figure-markdown_strict/unnamed-chunk-11-1.png)

## Functional analysis of CDS associated with high and low non-silent mutation rates

The functional analysis of the cds associated with high or low
non-silent mutation rates can add some interesting insight about what
functions are prone to modifications in the cohort.

Here we defined ‘High’ non-silent mutation rate has the CDS for which
non-silent variants represented 90% or more of the variants in a CDS.
Conversly, we defined low non-silent mutation rate as the CDS for which
non-silent variants represented 10% or less of the variants of a contig.

To perform the functional analysis, we obtained the COG annotations of
the CDS in the pangenome.

    # Extracting the contigs with high non silent mutation rates (ratio >90%) as the ones with high silent mutation rate (ratio <10%)
    CDS_rate_comp_high=subset(CDS_rate_comp, CDS_rate_comp$ratio>=0.90)
    CDS_rate_comp_low=subset(CDS_rate_comp, CDS_rate_comp$ratio<=0.1)

    # Importing eggNog-mapper information for the pangenome CDS

    eggnog=read.delim('dataset_analysis/data/variant_analysis/cds_eggNog.tsv', sep="\t")
    eggnog_short=eggnog %>% select(query, COG_category, Description) # select columns of interest in the whole eggNog output

    # Extracting COG information for the two CDS categories (high rate, low rate)

    eggnog_high=subset(eggnog_short, eggnog_short$query%in%CDS_rate_comp_high$contigs)
    colnames(eggnog_high)=c('contigs','COG','Description')
    eggnog_high_rate=left_join(CDS_rate_comp_high, eggnog_high, by='contigs')

    head(eggnog_high_rate,5)

    ##          contigs numb_var_ns numb_var_all Type  len rate_var_all ratio  COG
    ## 1 ABCDEHGJ_03231           2            2  CDS  591        295.5     1    S
    ## 2 ABCDEHGJ_03370           1            1  CDS  264        264.0     1 <NA>
    ## 3 ABCDEHGJ_03623           5            5  CDS  480         96.0     1 <NA>
    ## 4 ABCDEHGJ_03939           1            1  CDS  903        903.0     1    M
    ## 5 AIHHOGBN_03552           2            2  CDS 1164        582.0     1 <NA>
    ##                             Description
    ## 1 Type VI secretion system protein DotU
    ## 2                                  <NA>
    ## 3                                  <NA>
    ## 4                             chaperone
    ## 5                                  <NA>

    eggnog_low=subset(eggnog_short, eggnog_short$query%in%CDS_rate_comp_low$contigs)
    colnames(eggnog_low)=c('contigs','COG','Description')
    eggnog_low_rate=left_join(CDS_rate_comp_low, eggnog_low, by='contigs')

    head(eggnog_low_rate,5)

    ##          contigs numb_var_ns numb_var_all Type  len rate_var_all      ratio
    ## 1 ABCDEHGJ_03598           1           19  CDS 1299     68.36842 0.05263158
    ## 2 AKJGLPIO_00079           1           12  CDS 1023     85.25000 0.08333333
    ## 3 BCDBBOOO_01566           1           12  CDS 1227    102.25000 0.08333333
    ## 4 BCDBBOOO_02541           1           10  CDS  267     26.70000 0.10000000
    ## 5 BCDBBOOO_02935           1           16  CDS  363     22.68750 0.06250000
    ##    COG                                                 Description
    ## 1    K helix_turn_helix gluconate operon transcriptional repressor
    ## 2 <NA>                                                        <NA>
    ## 3 <NA>                                                        <NA>
    ## 4    K                   Prophage CP4-57 regulatory protein (AlpA)
    ## 5    H                        6-pyruvoyl tetrahydropterin synthase

To investigate potential biological functions associated with high and
low non-silent mutation, we investigated within each dataset (high and
low rates) the distribution of the different COG functional categories.
More specifically, we calculated for each COG category the percentage of
CDS it was found to be associated with.

       # Create a function that sorts the seen COG categories and how many times each category is annotated in the set of contig - The function also calculates the percentage of cds of the set, a given category is associated with

    egg_cat_inv=function(df_egg_rate) {
      df_egg_rate[is.na(df_egg_rate)]<-"-"
      df=as.data.frame(table(df_egg_rate$COG))
      
      df_count <- df %>%
      mutate(COG_cat = strsplit(as.character(Var1), "")) %>%
      unnest(COG_cat) %>%
      group_by(COG_cat) %>%
      summarise(Count = sum(Freq))
      
      df_count$prop=(df_count$Count/sum(df_count$Count))*100
      
      df_count
      
    }

     # Use the function to study the function distribution of the contig associated with high and low non-silent mutation rates

    egg_cat_high=egg_cat_inv(eggnog_high_rate)
    egg_cat_high$group='high'

    head(egg_cat_high,10)

    ## # A tibble: 10 × 4
    ##    COG_cat Count   prop group
    ##    <chr>   <int>  <dbl> <chr>
    ##  1 -         104 40.3   high 
    ##  2 C           3  1.16  high 
    ##  3 E           4  1.55  high 
    ##  4 F           2  0.775 high 
    ##  5 G           8  3.10  high 
    ##  6 H           1  0.388 high 
    ##  7 I           1  0.388 high 
    ##  8 J           2  0.775 high 
    ##  9 K          11  4.26  high 
    ## 10 L          19  7.36  high

    egg_cat_low=egg_cat_inv(eggnog_low_rate)
    egg_cat_low$group='low'

    head(egg_cat_low,10)

    ## # A tibble: 10 × 4
    ##    COG_cat Count  prop group
    ##    <chr>   <int> <dbl> <chr>
    ##  1 -          11 9.32  low  
    ##  2 A           1 0.847 low  
    ##  3 C           3 2.54  low  
    ##  4 D           6 5.08  low  
    ##  5 E           2 1.69  low  
    ##  6 F           3 2.54  low  
    ##  7 G           7 5.93  low  
    ##  8 H           4 3.39  low  
    ##  9 J           1 0.847 low  
    ## 10 K           8 6.78  low

    cog_fun_cat=read.csv('dataset_analysis/data/variant_analysis/COG_functional_categories.csv') # COG functional categories information

    # Data reformatting before plotting COG category distribution for high and low non-silent variant rates
    data_plot=bind_rows(egg_cat_high, egg_cat_low)
    data_plot=left_join(data_plot, cog_fun_cat, by='COG_cat')
    data_plot$COG_family[is.na(data_plot$COG_family)]<-"Poorly characterized"


    # Plot
    plot_comp= ggplot(data_plot, aes(x=COG_cat, y=group, col=COG_cat)) +
      geom_point(aes(size=prop, col=COG_family))+
      theme_arcadia(x_axis_type = 'categorical', y_axis_type = 'categorical') +
      xlab('COG functional categories') +
      scale_color_arcadia(palette_name = 'primary',reverse=FALSE) 

    plot_comp

![](Variant_population_analysis_files/figure-markdown_strict/unnamed-chunk-14-1.png)

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] arcadiathemeR_0.1.0 ape_5.7-1           ggridges_0.5.6     
    ##  [4] seqinr_4.2-36       lubridate_1.9.3     forcats_1.0.0      
    ##  [7] stringr_1.5.1       dplyr_1.1.4         purrr_1.0.2        
    ## [10] readr_2.1.5         tidyr_1.3.1         tibble_3.2.1       
    ## [13] ggplot2_3.5.1       tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1  xfun_0.43         lattice_0.22-6    colorspace_2.1-0 
    ##  [5] vctrs_0.6.5       generics_0.1.3    htmltools_0.5.8.1 yaml_2.3.8       
    ##  [9] utf8_1.2.4        rlang_1.1.4       pillar_1.9.0      glue_1.7.0       
    ## [13] withr_3.0.0       lifecycle_1.0.4   munsell_0.5.1     gtable_0.3.5     
    ## [17] evaluate_0.23     labeling_0.4.3    knitr_1.46        tzdb_0.4.0       
    ## [21] fastmap_1.1.1     parallel_4.2.3    fansi_1.0.6       highr_0.10       
    ## [25] Rcpp_1.0.12       scales_1.3.0      showtext_0.9-7    sysfonts_0.8.9   
    ## [29] farver_2.1.2      systemfonts_1.0.6 hms_1.1.3         digest_0.6.35    
    ## [33] stringi_1.8.3     showtextdb_3.0    grid_4.2.3        ade4_1.7-22      
    ## [37] cli_3.6.3         tools_4.2.3       magrittr_2.0.3    pkgconfig_2.0.3  
    ## [41] MASS_7.3-60.0.1   timechange_0.3.0  rmarkdown_2.26    rstudioapi_0.16.0
    ## [45] R6_2.5.1          nlme_3.1-164      compiler_4.2.3
