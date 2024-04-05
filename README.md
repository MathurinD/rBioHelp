rBioHelp
================
Mathurin Dorel
2024-04-05

A set of functions that implement common operations used when analysing
human biological datasets. Aggregates several packages I like and
versions of datasets I use.

``` r
    library(rBioHelp)
```

    ## Loading required package: tidyverse

    ## ── Attaching core tidyverse packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.3     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
    ## Loading required package: ComplexHeatmap
    ## 
    ## Loading required package: grid
    ## 
    ## ========================================
    ## ComplexHeatmap version 2.15.4
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
    ## 
    ## 
    ## Loading required package: GO.db
    ## 
    ## Loading required package: AnnotationDbi
    ## 
    ## Loading required package: stats4
    ## 
    ## Loading required package: BiocGenerics
    ## 
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union
    ## 
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     Filter, Find, Map, Position, Reduce, anyDuplicated, append, as.data.frame, basename, cbind, colnames, dirname, do.call,
    ##     duplicated, eval, evalq, get, grep, grepl, intersect, is.unsorted, lapply, mapply, match, mget, order, paste, pmax, pmax.int,
    ##     pmin, pmin.int, rank, rbind, rownames, sapply, setdiff, sort, table, tapply, union, unique, unsplit, which.max, which.min
    ## 
    ## 
    ## Loading required package: Biobase
    ## 
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with 'browseVignettes()'. To cite Bioconductor, see 'citation("Biobase")', and for
    ##     packages 'citation("pkgname")'.
    ## 
    ## 
    ## Loading required package: IRanges
    ## 
    ## Loading required package: S4Vectors
    ## 
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-
    ## 
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     I, expand.grid, unname
    ## 
    ## 
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%
    ## 
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## 
    ## 
    ## Attaching package: 'AnnotationDbi'
    ## 
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select
    ## 
    ## 
    ## 
    ## 
    ## Loading required package: Rgraphviz
    ## 
    ## Loading required package: graph
    ## 
    ## 
    ## Attaching package: 'graph'
    ## 
    ## 
    ## The following object is masked from 'package:stringr':
    ## 
    ##     boundary
    ## 
    ## 
    ## 
    ## Attaching package: 'Rgraphviz'
    ## 
    ## 
    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     from, to
    ## 
    ## 
    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     from, to
    ## 
    ## 
    ## Loading required package: clusterProfiler
    ## 
    ## clusterProfiler v4.5.1  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
    ## 
    ## 
    ## Attaching package: 'clusterProfiler'
    ## 
    ## 
    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select
    ## 
    ## 
    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice
    ## 
    ## 
    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename
    ## 
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify
    ## 
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter
    ## 
    ## 
    ## Loading required package: eulerr
    ## 
    ## Loading required package: org.Hs.eg.db
    ## 
    ## 
    ## 
    ## Loading required package: wikiprofiler

``` r
    library(tidyverse)
    library(ComplexHeatmap)
    library(org.Hs.eg.db)
    # Toy data for illustration
    geneslist = c('MAP2K2','MAP2K1','MAPK1','MAPK3','PIK3CA','PIK3CG', 'MTOR','TSC1','TSC2','EGFR','GRB2','RAF1','BRAF','KRAS','PTPN11') # Example signalling
    geneslist2 = c('MCM5','PCNA','TYMS','FEN1','MCM7','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1','UHRF1','CENPU','HELLS','RFC2','POLR1B','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1','CHAF1B','MRPL36','E2F8')
    clusters = tibble(Cluster = c(rep(1, length(geneslist)), rep(2, length(geneslist2))), Genes=c(geneslist, geneslist2) ) %>% mutate(Entrez=mapIds(org.Hs.eg.db, Genes, 'ENTREZID','SYMBOL')) %>% mutate(Ensembl=mapIds(org.Hs.eg.db, Genes, 'ENSEMBL','SYMBOL')) %>% mutate(Uniprot=mapIds(org.Hs.eg.db, Genes, 'UNIPROT','SYMBOL'))
```

    ## 'select()' returned 1:1 mapping between keys and columns
    ## 'select()' returned 1:many mapping between keys and columns
    ## 'select()' returned 1:many mapping between keys and columns

``` r
    zscores = matrix(c(rnorm(5*58, -2), rnorm(5*58, 2)), nrow=nrow(clusters), dimnames=list(clusters$Ensembl, paste0('S', 1:10))) %>% as_tibble(rownames='g') # Pseudo z-score with two groups
    patients_info = tibble(Name=paste0('S', 1:10), Status=c(rep('Sick', 5), rep('Healthy', 5))) %>% column_to_rownames('Name')
```

The package provides several helper functions to compare gene clusters
(e.g DEGs after different treatments, coregulated modules, etc).

``` r
    # Simply plot the expression/zscores with sample annotation
    compareClusters(zscores, clusters, 1:2) # TODO have patients_info as argument
```

    ## Warning: The input is a data frame-like object, convert it to a matrix.

![](README_files/figure-gfm/genesets-1.png)<!-- -->

``` r
    # Alternative to clusterProfiler::dotplot
    wpenrich = enrichWP(clusters$Entrez, 'Homo sapiens')
    plotEnrichment(wpenrich)
```

![](README_files/figure-gfm/genesets-2.png)<!-- -->

``` r
    # Compute enrichment using multiple gene set sources for each gene group
    clusters_enrichment = suppressMessages(genesetsOverview(clusters))
    clusters_enrichment %>% tail %>% pretty_rich()
```

| ID                                                | Description                                       |    pvalue | genes                          |
|:--------------------------------------------------|:--------------------------------------------------|----------:|:-------------------------------|
| Macroautophagy                                    | Macroautophagy                                    | 0.0146905 | MTOR, TSC2                     |
| Regulation of HSF1-mediated heat shock response   | Regulation of HSF1-mediated heat shock response   | 0.0176852 | MAPK1, MAPK3                   |
| Senescence-Associated Secretory Phenotype (SASP)  | Senescence-Associated Secretory Phenotype (SASP)  | 0.0199697 | MAPK1, MAPK3                   |
| TP53 Regulates Metabolic Genes                    | TP53 Regulates Metabolic Genes                    | 0.0238671 | MTOR, TSC2                     |
| Constitutive Signaling by Aberrant PI3K in Cancer | Constitutive Signaling by Aberrant PI3K in Cancer | 0.0243747 | PIK3CA, PTPN11                 |
| unattributed                                      | Genes not in an enriched gene set                 | 1.0000000 | CDCA7 , CENPU , WDR76 , MRPL36 |

A graph between gene sets helps to better visualise the coverage by the
various gene sets enriched.

``` r
    clusters %>% group_by(Cluster) %>% group_map(function(xx,yy){genesetsOverview(xx) %>% mutate(Cluster=yy$Cluster) %>% suppressMessages}) -> c_riches # Clusters enrichment
    clusterSetsGraph(clusters, c_riches, top_sets=5, max_sets=50, graphic_scaling=0.8)
```

![](README_files/figure-gfm/genesets_graph-1.png)<!-- -->

``` r
    clusterSetsGraph(clusters, c_riches %>% lapply(filter, grepl('HALLMARK', ID)), top_sets=10, max_sets=50, graphic_scaling=0.8) # Restrict the plot to hallmark gene sets
```

![](README_files/figure-gfm/genesets_graph-2.png)<!-- -->

Improve wikiprofiler plots, allowing control over the color gradient and
plotting multiple datasets on the same pathway TODO list: sample name
legend

``` r
#wpplot.R:WPplot <- function(ftable, pathwayid='WP4172', high="red", low="blue", legend = TRUE, legend_x = 0.001, legend_y = 0.94, strip_out=0.01, fixed_range=F, lim=5) {
#wpplot.R:WPmultiplot <- function(ftable, pathwayid='WP4172', multi='Cell_line', high="red", low="blue", legend = TRUE, legend_x = 0.001, legend_y = 0.94, strip_out=0.01, fixed_range=F, lim=5) {
#wpplot.R:#' Helper function to change the background for a gene
#wpplot.R:replace_bg <- function(svg, gene, positions, colors) {
#wpplot.R:WPbgfill <- function(pp, value, high="red", low="blue", legend = TRUE, legend_x = 0.001, legend_y = 0.94, strip_out=0.01, fixed_range=F, lim=5) {

#Eventually will automatically update datasets on install (and keep old versions for reproducibility) once I find the time to code that.
```

Developped by Mathurin Dorel
