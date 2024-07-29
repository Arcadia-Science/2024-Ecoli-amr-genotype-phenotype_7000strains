### Objective

This R notebook contains the code used to analyse: - some meta
information associated with the different strains used in this study
(all strains but the ECOR72) - the antimicrobial resistance (AMR)
phenotypes in the dataset

------------------------------------------------------------------------

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

    library(ggwordcloud)
    library(arcadiathemeR)
    library(maps)

    ## 
    ## Attaching package: 'maps'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     map

------------------------------------------------------------------------

# Meta-information analysis of the cohort

## Data import

     ## Samples list for which we obtained genome sequencing files from SRA accession number
    sample_list=read.csv('dataset_generation/data/sample_list_SRA.csv')

     ## Genome information from BVBRC
    genome_data=read.csv('dataset_analysis/data/dataset_analysis/BVBRC_genome_May31.csv')

## Selection of the strains used in this study

    # 1- Extracting the Genome Names in the BV-BRC dataset associated with the SRA accession number we used to obtain the sequencing files

    samples=c(unique(sample_list$samples))
    genome_samples=subset(genome_data, genome_data$SRA.Accession%in%samples)

       ## Some Genome Names are associated with multiple SRA accession number and we need to treat them separately 

    out=setdiff(samples, genome_samples$SRA.Accession)
    pattern <- paste(out, collapse="|")
    matching_rows <- grep(pattern, genome_data$SRA.Accession)
    subset_df <- genome_data[matching_rows, ]

       ## And combine them 

    genome_samples=rbind(genome_samples,subset_df) # dataframe with all meta-information of the strains used in this study (excluding the ECOR72 strains)

## Selection of metadata of interest

The metadata we want to investigate are: - genome size (number of bp) -
number of CDS - collection year - isolation country - host

    metadata=genome_samples %>%
      select(Genome.Name, Size, CDS, Collection.Year, Isolation.Country, Host.Name)

    empty_counts <- colSums(metadata == "")

    empty_counts # indicates for each metadata type the number of strains for which the information is missing

    ##       Genome.Name              Size               CDS   Collection.Year 
    ##                 0                 0                 0                NA 
    ## Isolation.Country         Host.Name 
    ##                 3              3061

### Genome size

    min_size=min(metadata$Size)
    max_size=max(metadata$Size)
    mean_size=mean(metadata$Size)
    median_size=median(metadata$Size)

    size_plot <- ggplot(metadata, aes(x = Size)) +
      geom_histogram(binwidth = 50000, fill='#7A77AB',col='#484B50') +
      xlim(4000000, 6000000)+
      theme_arcadia(x_axis_type = "numerical") +
      xlab('Genome size (bp)')+
      ylab('Number of strains')

    min_size

    ## [1] 4069292

    max_size

    ## [1] 5980243

    mean_size

    ## [1] 5071975

    median_size

    ## [1] 5088889

    size_plot

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](Dataset_metainfo_AMR_analysis_files/figure-markdown_strict/unnamed-chunk-5-1.png)
\### CDS number

    min_cds=min(metadata$CDS)
    max_cds=max(metadata$CDS)
    mean_cds=mean(metadata$CDS)

    min_cds

    ## [1] 3992

    max_cds

    ## [1] 6135

    mean_cds

    ## [1] 5068.812

### Collection year

    year_counts <- table(metadata$Collection.Year)
    year_counts

    ## 
    ## 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 
    ##   74  184  217  226  223  224  286  347  431  361  380  301  243  320  485  342 
    ## 2017 
    ##  346

### Host

    hosts_counts <- as.data.frame(table(metadata$Host.Name))
    hosts_counts

    ##                          Var1 Freq
    ## 1                             3061
    ## 2                  Bos taurus   18
    ## 3            Cat, Felis catus    2
    ## 4      Chicken, Gallus gallus    1
    ## 5 Dog, Canis lupus familiaris    2
    ## 6         Human, Homo sapiens 3898
    ## 7             Pig, Sus scrofa    1

### Country

    country_counts <- as.data.frame(table(metadata$Isolation.Country))
    colnames(country_counts)=c('country','count')

    #Reformatting the data to use them with the map package and represent the data on a world map. This includes removing the strains without country information and consolidating strains from the United Kingdom and England into a single count.

    country_for_map=country_counts[-1,]
    country_for_map <- data.frame(
      country = c("United Kingdom", setdiff(country_for_map$country, c("England", "United Kingdom"))),
      count = c(sum(country_for_map$count[country_for_map$country %in% c("England", "United Kingdom")]), 
                country_for_map$count[!country_for_map$country %in% c("England", "United Kingdom")])
    )

    # Changing United Kingdom to UK to match the country name in the map package
    country_for_map$country <- gsub("United Kingdom", "UK", country_for_map$country)

    # Mapping the data onto a world map
    ## Generating the data
    world_map <- map_data("world")
    merged_data <- merge(world_map, country_for_map, by.x = "region", by.y = "country", all.x = TRUE)

    # Plot the world map with ggplot2
    ggplot() +
      geom_map(data = merged_data, map = world_map,
               aes(x = long, y = lat, map_id = region, fill = count),
               color = "white", size = 0.1) +
     scale_fill_gradient(low = "#282A49", high = "#97CD78", name = "Value") +
      theme_arcadia() 

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning in geom_map(data = merged_data, map = world_map, aes(x = long, y = lat,
    ## : Ignoring unknown aesthetics: x and y

![](Dataset_metainfo_AMR_analysis_files/figure-markdown_strict/unnamed-chunk-9-1.png)

------------------------------------------------------------------------

# AMR phenotype analysis

In this part we investigate the AMR information available for the
strains studied in the Pub.

## Data import

    amr_data=read.csv('dataset_analysis/data/dataset_analysis/BVBRC_genome_amr_May31.csv')

    # Extracting antibiotic info for the strains using their genomes names
    amr_data_sub=subset(amr_data, amr_data$Genome.Name%in%genome_samples$Genome.Name)

    #Correspondence antibiotic - antibiotic class
    antibiotic_class=read.csv('dataset_analysis/data/dataset_analysis/Antibiotic_class.csv')

## Reformatting the data, removing inconsistent information

For some strains, multiple AMR phenotypes were reported for the same
antibiotic, when these AMR phenotypes are inconsistent we remove them.

    # Keeping consistent repetitive AMR phenotype but keeping only a single row entry

    amr_data_sub_short=amr_data_sub%>% select(Genome.Name, Antibiotic, Resistant.Phenotype)%>%
      filter(Resistant.Phenotype != "")  %>% distinct()

    # Removing any inconsistent AMR phenotype for the same strain and antibiotic

    amr_data_filtered <- amr_data_sub_short %>%
      group_by(Genome.Name, Antibiotic) %>%
      filter(n_distinct(Resistant.Phenotype) == 1) %>%
      ungroup()

## Generating the phenotype matrix

The phenotype matrix provides AMR information for all strains (rows) and
antibiotics (columns) NAs indicate that no AMR information available

      # Preparing phenotype matrix

    amr_data_matrix=amr_data_filtered%>%spread(Antibiotic, Resistant.Phenotype)

    #write.csv(amr_data_matrix,"dataset_analysis/results/phenotype_matrix_08302024.csv")

## Data exploration

In this part, we explore and visualize different information related to
AMR phenotypes

### Number of Antibiotics or Combinations of antibiotics

    antibiotics=c(unique(amr_data_filtered$Antibiotic)) # 50 antibiotics with at least 1 known phenotype
    length(antibiotics)

    ## [1] 50

## Number of known AMR phenotypes per strain

    # For each strains counts how many AMR phenotypes are available
    count_sample=amr_data_filtered%>%group_by(Genome.Name)%>%
      summarize(count=n())

    head(count_sample,5)

    ## # A tibble: 5 × 2
    ##   Genome.Name                                           count
    ##   <chr>                                                 <int>
    ## 1 Escherichia coli 00116690-7bb9-11e9-a8d3-68b59976a384    10
    ## 2 Escherichia coli 002b98bc-7bb9-11e9-a8d3-68b59976a384    10
    ## 3 Escherichia coli 00455662-7bb9-11e9-a8d3-68b59976a384     9
    ## 4 Escherichia coli 006393f2-7bb9-11e9-a8d3-68b59976a384     9
    ## 5 Escherichia coli 007f095c-7bb9-11e9-a8d3-68b59976a384     9

    # Visualization of the distribution of known AMR phenotype per strain
    amrphen_per_strain=ggplot(count_sample, aes(x=count))+
      geom_histogram(bins=max(count_sample$count), fill='#97CD78', col='#484B50') + 
      xlab('# of knwon AMR phenotypes in a strain') +
      ylab ('# of strains in the cohort')+
      theme_arcadia(x_axis_type = "numerical")

    amrphen_per_strain

![](Dataset_metainfo_AMR_analysis_files/figure-markdown_strict/unnamed-chunk-14-1.png)

## Number of strains with known AMR phenotype per antibiotic

    # For each antibiotic we look at how many strains have available AMR phenotype information
    count_antibio=amr_data_filtered%>%group_by(Antibiotic)%>%
      summarize(count=n())

    head(count_antibio,10)

    ## # A tibble: 10 × 2
    ##    Antibiotic                  count
    ##    <chr>                       <int>
    ##  1 amikacin                     1827
    ##  2 amoxicillin                  1131
    ##  3 amoxicillin/clavulanic acid  3779
    ##  4 ampicillin                   4874
    ##  5 ampicillin/sulbactam           85
    ##  6 azithromycin                  163
    ##  7 aztreonam                     886
    ##  8 cefalexin                       1
    ##  9 cefalotin                     377
    ## 10 cefazolin                     192

    amrphen_per_antibiotic=ggplot(count_antibio, aes(x=count))+
      geom_histogram(binwidth =100, fill='#C85152', col='#484B50') + 
      xlab('# of knwown AMR phenotypes') +
      ylab ('# of antibiotics')+
      theme_arcadia(x_axis_type = "numerical")

    amrphen_per_antibiotic

![](Dataset_metainfo_AMR_analysis_files/figure-markdown_strict/unnamed-chunk-15-1.png)

## Antibiotics with 500 or more AMR phenotypes available: analysis of AMR distribution

    # Selecting antibiotics with more than 500 AMR phenotypes 
    list_antibio_high <- count_antibio %>%
      filter(count >= 500) %>%
      pull(Antibiotic) %>%
      unique()

    df_antibio_high=subset(amr_data_filtered, amr_data_filtered$Antibiotic%in%list_antibio_high)

    # For each antibiotic counts how many 'Resistant, 'Susceptible' and 'Intermediary' phenotype
    phenotypes_counts <- df_antibio_high %>%
      group_by(Antibiotic, Resistant.Phenotype) %>%
      summarise(Count = n(), .groups = 'drop')

    phen_distribution=ggplot(phenotypes_counts, aes(x=Antibiotic, y=Count, fill=Resistant.Phenotype))+
      geom_bar(stat='identity') + 
      xlab('Antibiotic') + ylab('# of strain')+
      theme_arcadia(x_axis_type = "categorical") +
      scale_fill_arcadia(palette_name = "secondary", reverse = TRUE)+
      theme(axis.text.x = element_text(angle=90,hjust=1))

    phen_distribution

![](Dataset_metainfo_AMR_analysis_files/figure-markdown_strict/unnamed-chunk-16-1.png)

## Resistant strains

    Resistant=subset(amr_data_filtered, amr_data_filtered$Resistant.Phenotype=='Resistant')

    list_antibio_resist=c(unique(Resistant$Antibiotic))

    #Number of resistant AMR per strain

    res_count=Resistant%>%group_by(Genome.Name)%>%
      summarize(count=n())

    head(res_count)

    ## # A tibble: 6 × 2
    ##   Genome.Name                                           count
    ##   <chr>                                                 <int>
    ## 1 Escherichia coli 00116690-7bb9-11e9-a8d3-68b59976a384     2
    ## 2 Escherichia coli 002b98bc-7bb9-11e9-a8d3-68b59976a384     3
    ## 3 Escherichia coli 00455662-7bb9-11e9-a8d3-68b59976a384     1
    ## 4 Escherichia coli 006393f2-7bb9-11e9-a8d3-68b59976a384     2
    ## 5 Escherichia coli 007f095c-7bb9-11e9-a8d3-68b59976a384     2
    ## 6 Escherichia coli 00b33114-7bb9-11e9-a8d3-68b59976a384     4

    #Resistance per antibiotic class

    res_class=Resistant %>% left_join(antibiotic_class, by="Antibiotic")%>%
      group_by(Class)%>%
      summarize(count=n())

    res_per_class=ggplot(res_class, aes(label=Class,size=count, color=Class))+
      geom_text_wordcloud() +
      scale_size_area(max_size = 10)+
      theme_arcadia() +
      scale_color_arcadia(palette_name = "primary")

    res_per_class

![](Dataset_metainfo_AMR_analysis_files/figure-markdown_strict/unnamed-chunk-17-1.png)

## Multi-drug resistant strains

Here we define Multi-drug resistant bacteria if they are resistant to 3
or more different class of antibiotic

    MDR_candidate=res_count %>%
      filter(count >= 3) %>%
      pull(Genome.Name) %>%
      unique()

    MDR_candidate_DF=subset(Resistant, Resistant$Genome.Name%in%MDR_candidate)
    MDR_candidate_DF=left_join(MDR_candidate_DF, antibiotic_class, by="Antibiotic")

    strains=c(unique(MDR_candidate_DF$Genome.Name))

    length(strains)

    ## [1] 1425

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
    ##  [1] maps_3.4.2          arcadiathemeR_0.1.0 ggwordcloud_0.6.2  
    ##  [4] lubridate_1.9.3     forcats_1.0.0       stringr_1.5.1      
    ##  [7] dplyr_1.1.4         purrr_1.0.2         readr_2.1.5        
    ## [10] tidyr_1.3.1         tibble_3.2.1        ggplot2_3.5.1      
    ## [13] tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1  xfun_0.43         colorspace_2.1-0  vctrs_0.6.5      
    ##  [5] generics_0.1.3    htmltools_0.5.8.1 yaml_2.3.8        utf8_1.2.4       
    ##  [9] rlang_1.1.4       gridtext_0.1.5    pillar_1.9.0      glue_1.7.0       
    ## [13] withr_3.0.0       lifecycle_1.0.4   commonmark_1.9.1  munsell_0.5.1    
    ## [17] gtable_0.3.5      evaluate_0.23     labeling_0.4.3    knitr_1.46       
    ## [21] tzdb_0.4.0        fastmap_1.1.1     markdown_1.12     fansi_1.0.6      
    ## [25] highr_0.10        Rcpp_1.0.12       scales_1.3.0      showtext_0.9-7   
    ## [29] sysfonts_0.8.9    farver_2.1.2      systemfonts_1.0.6 hms_1.1.3        
    ## [33] png_0.1-8         digest_0.6.35     stringi_1.8.3     showtextdb_3.0   
    ## [37] grid_4.2.3        cli_3.6.3         tools_4.2.3       magrittr_2.0.3   
    ## [41] crayon_1.5.2      pkgconfig_2.0.3   xml2_1.3.6        timechange_0.3.0 
    ## [45] rmarkdown_2.26    rstudioapi_0.16.0 R6_2.5.1          compiler_4.2.3
