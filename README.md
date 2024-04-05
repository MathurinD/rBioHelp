rBioHelp
================
Mathurin Dorel
2024-04-05

A set of functions that implement common operations used when analysing
human biological datasets. Aggregates several packages I like and
versions of datasets I use.

``` r
    knitr::opts_chunk$set(warning=FALSE, cache=TRUE, message=FALSE)
    suppressPackageStartupMessages(library(rBioHelp))
    library(tidyverse)
    library(ComplexHeatmap)
    library(org.Hs.eg.db)
    # Toy data for illustration
    geneslist = c('MAP2K2','MAP2K1','MAPK1','MAPK3','PIK3CA','PIK3CG', 'MTOR','TSC1','TSC2','EGFR','GRB2','RAF1','BRAF','KRAS','PTPN11') # Example signalling
    geneslist2 = c('MCM5','PCNA','TYMS','FEN1','MCM7','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1','UHRF1','CENPU','HELLS','RFC2','POLR1B','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1','CHAF1B','MRPL36','E2F8')
    clusters = tibble(Cluster = c(rep(1, length(geneslist)), rep(2, length(geneslist2))), Genes=c(geneslist, geneslist2) ) %>% mutate(Entrez=mapIds(org.Hs.eg.db, Genes, 'ENTREZID','SYMBOL')) %>% mutate(Ensembl=mapIds(org.Hs.eg.db, Genes, 'ENSEMBL','SYMBOL')) %>% mutate(Uniprot=mapIds(org.Hs.eg.db, Genes, 'UNIPROT','SYMBOL'))

    # Pseudo z-score with two groups
    zscores = matrix(c(rnorm(5*58, -2), rnorm(5*58, 2)), nrow=nrow(clusters),
                     dimnames=list(clusters$Ensembl, paste0('S', 1:10))) %>% as_tibble(rownames='g')
    patients_info = tibble(
                           Name=paste0('S', 1:10),
                           Status=c(rep('Sick', 5), rep('Healthy', 5))
                    ) %>% column_to_rownames('Name')
```

The package provides several helper functions to compare gene clusters
(e.g DEGs after different treatments, coregulated modules, etc).

The `genesetsOverview` function performs enrichment on KEGG,
Wikipathway, GO, MsigDB hallmarks and Reactome to give a full overview
of the roles of the genes in the list. A speficic geneset source can
easily be extracted by filtering for ID patterns (hsa*, WP*, GO:*,
HALLMARK* and non matching respectivel).<br/> The last entry contains
the genes from the lists that are not in any of the enriched gene sets.
Those singletons highlight a gap in the gene sets databases (and
sometimes in knowledge).

``` r
    # Simply plot the expression/zscores with sample annotation
    compareClusters(zscores, clusters, 1:2) # TODO have patients_info as argument
```

![](README_files/figure-gfm/genesets-1.png)<!-- -->

``` r
    customVenn(list(MAPK=geneslist, S_Phase=geneslist2)) # Compare gene sets
```

![](README_files/figure-gfm/genesets-2.png)<!-- -->

``` r
    # Compute enrichment using multiple gene set sources for each gene group
    wpenrich = suppressMessages(genesetsOverview(clusters))
    plotEnrichment(wpenrich) # Alternative to clusterProfiler::dotplot
```

![](README_files/figure-gfm/genesets-3.png)<!-- -->

``` r
    # List of genes in the genesets
    wpenrich %>% {bind_rows(head(.), tail(.))} %>% pretty_rich()
```

| ID                                                | Description                                                                      |    pvalue | genes                                                                                          |
|:--------------------------------------------------|:---------------------------------------------------------------------------------|----------:|:-----------------------------------------------------------------------------------------------|
| WP466                                             | DNA replication                                                                  | 0.0000000 | MCM5 , PCNA , MCM7 , MCM4 , MCM6 , PRIM1, RFC2 , GMNN , POLD3, CDC45, CDC6 , POLA1             |
| WP2261                                            | Glioblastoma signaling pathways                                                  | 0.0000000 | MAP2K2, MAP2K1, MAPK1 , MAPK3 , PIK3CA, PIK3CG, TSC1 , TSC2 , EGFR , GRB2 , RAF1 , BRAF , KRAS |
| WP2446                                            | Retinoblastoma gene in cancer                                                    | 0.0000000 | RAF1 , PCNA , TYMS , MCM7 , MCM4 , RRM1 , MCM6 , PRIM1, CCNE2, POLD3, RRM2 , CDC45, POLA1      |
| WP3879                                            | 4 hydroxytamoxifen dexamethasone and retinoic acids regulation of p27 expression | 0.0000000 | MAP2K2, MAP2K1, MAPK1 , MAPK3 , PIK3CA, MTOR , TSC1 , TSC2 , RAF1                              |
| WP2263                                            | Androgen receptor network in prostate cancer                                     | 0.0000000 | MAP2K2, MAP2K1, MAPK1 , MAPK3 , PIK3CA, MTOR , TSC1 , TSC2 , GRB2 , RAF1 , PTPN11, MSH2 , BLM  |
| WP3676                                            | BDNF TrkB signaling                                                              | 0.0000000 | MAP2K1, MAPK1 , PIK3CG, MTOR , TSC1 , TSC2 , GRB2 , BRAF , KRAS                                |
| Macroautophagy                                    | Macroautophagy                                                                   | 0.0146905 | MTOR, TSC2                                                                                     |
| Regulation of HSF1-mediated heat shock response   | Regulation of HSF1-mediated heat shock response                                  | 0.0176852 | MAPK1, MAPK3                                                                                   |
| Senescence-Associated Secretory Phenotype (SASP)  | Senescence-Associated Secretory Phenotype (SASP)                                 | 0.0199697 | MAPK1, MAPK3                                                                                   |
| TP53 Regulates Metabolic Genes                    | TP53 Regulates Metabolic Genes                                                   | 0.0238671 | MTOR, TSC2                                                                                     |
| Constitutive Signaling by Aberrant PI3K in Cancer | Constitutive Signaling by Aberrant PI3K in Cancer                                | 0.0243747 | PIK3CA, PTPN11                                                                                 |
| unattributed                                      | Genes not in an enriched gene set                                                | 1.0000000 | CDCA7 , CENPU , WDR76 , MRPL36                                                                 |

A graph between gene sets helps to better visualise the coverage by the
various gene sets enriched.<br/> NOTE: The implementation does not show
the overlaps between the gene sets, which should always be considered
when analysing gene set enrichment results,

``` r
    clusters_enrichment = clusters %>% group_by(Cluster) %>%
        group_map(function(xx,yy){genesetsOverview(xx) %>% mutate(Cluster=yy$Cluster) %>% suppressMessages})
    clusterSetsGraph(clusters, clusters_enrichment, top_sets=5, max_sets=50, graphic_scaling=0.8)
```

![](README_files/figure-gfm/genesets_graph-1.png)<!-- -->

``` r
    clusterSetsGraph(clusters, clusters_enrichment %>% lapply(filter, grepl('HALLMARK', ID)),
                     top_sets=10, max_sets=50, graphic_scaling=0.8) # Restrict the plot to hallmark gene sets
```

![](README_files/figure-gfm/genesets_graph-2.png)<!-- -->

Improve wikiprofiler plots, allowing control over the color gradient and
plotting multiple datasets on the same pathway TODO list: - sample name
legend - fix little bug

``` r
    WPmultiplot(zscores %>% column_to_rownames('g'), pathwayid='WP2261')
```

![](README_files/figure-gfm/wikipathway-1.png)<!-- -->

``` r
    WPmultiplot(zscores %>% column_to_rownames('g'), pathwayid='WP3879') # Little bug
```

![](README_files/figure-gfm/wikipathway-2.png)<!-- -->

``` r
    WPmultiplot(zscores %>% column_to_rownames('g'), pathwayid='WP466')
```

![](README_files/figure-gfm/wikipathway-3.png)<!-- -->
