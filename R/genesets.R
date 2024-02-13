#' @import clusterProfiler
#' @import tidyverse

library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
reactome = read_tsv('/project/pe_data/reference/genesets/UniProt2Reactome.txt', col_names=c('Uniprot','RPID','URL','Name','Evidence','Species')) %>% filter(Species=='Homo sapiens') %>% mutate(Gene=mapIds(org.Hs.eg.db, keys=Uniprot, column='SYMBOL', keytype='UNIPROT'))
msigdb = read.gmt('/project/pe_data/reference/genesets/h.all.v7.4.symbols.gmt') %>% as_tibble%>% mutate(Entrez=mapIds(org.Hs.eg.db, keys=gene, column='ENTREZID', keytype='SYMBOL'))
clusterProfiler:::prepare_KEGG('hsa', "KEGG", 'kegg') -> KEGG_DATA

#' 
#' @param clusters CLICK cluster output with colums Cluster, Entrez and Uniprot
#' @rdname genesetsOverview
genesetsOverview <- function(clusters, ...)  {
    kegg_rich = clusters %>% pull(Entrez) %>% enrichKEGG('hsa', ...)
    wp_rich = clusters %>% pull(Entrez) %>% enrichWP('Homo sapiens', ...)
    reactome_rich = clusters %>% pull(Uniprot) %>% enricher(TERM2GENE=reactome %>% distinct(Name, Uniprot), ...)
    msig_rich = clusters %>% pull(Entrez) %>% enricher(TERM2GENE=msigdb%>% distinct(term, Entrez), ...)
    go_rich =clusters %>% pull('Entrez') %>% enrichGO(org.Hs.eg.db, ...)

    enriched_sets = NULL
    if (!is.null(wp_rich))
        enriched_sets = bind_rows(enriched_sets, wp_rich %>% as_tibble %>% mutate(genes=str_split(geneID, '/')))
    if (!is.null(kegg_rich))
        enriched_sets = bind_rows(enriched_sets, kegg_rich %>% as_tibble %>% mutate(genes=str_split(geneID, '/')) )
    if (!is.null(go_rich))
        enriched_sets = bind_rows(enriched_sets, go_rich %>% as_tibble %>% mutate(genes=str_split(geneID, '/')) )
    if (!is.null(msig_rich))
        enriched_sets = bind_rows(enriched_sets, msig_rich %>% as_tibble %>% mutate(genes=str_split(geneID, '/')) )
    if (!is.null(reactome_rich))
        enriched_sets = bind_rows(enriched_sets, reactome_rich %>% as_tibble %>% mutate(genes=str_split(geneID, '/')) %>% mutate(genes=lapply(genes, function(gg){mapIds(org.Hs.eg.db, gg, 'ENTREZID','UNIPROT')})))

    if (is.null(enriched_sets)) { stop('No gene sets enriched') }
    
    enriched_sets %>% pull(genes) %>% unlist %>% unique -> attributed_genes
    enriched_sets %>% bind_rows( tibble(ID='unattributed', Description='Genes not in an enriched gene set', pvalue=1, genes=list(clusters$Entrez[!clusters$Entrez %in% attributed_genes])) ) -> enriched_sets

    return(enriched_sets)
}
#' Wrapper for CLICK cluster enrichment.
#' @param clusters CLICK cluster output with colums Cluster, Entrez and Uniprot
#' @param selection Clusters on which gene set enrichment should be performed
#' @rdname genesetsOverview
genesetsFromCluster <- function(clusters, selection, ...) { genesetsOverview(clusters %>% filter(Cluster %in% selection), ...) }

#' Nicely display 
pretty_rich <- function(rich){
    knitr::kable(rich %>% mutate(genes=lapply(genes, function(gg){suppressMessages(mapIds(org.Hs.eg.db, gg, 'SYMBOL', 'ENTREZID'))})) %>% dplyr::select(ID, Description, pvalue, genes))
}
#' Plot heatmap for multiple clusters
#' @param zscores Z-normalised expression data with genes in rownames and samples in colnames.
#' @param clusters CLICK cluster output with colums Cluster and Ensembl
#' @param cids Ids of the clusters to plot
#' @rdname genesetsOverview
compareClusters <- function(zscores, clusters, cids) {
    fcl = clusters %>% filter(Cluster %in% cids) %>% mutate(Cluster=as.factor(Cluster)) %>% select(Cluster, Ensembl)
    print( zscores %>% filter(g %in% fcl$Ensembl) %>% column_to_rownames('g') %>% Heatmap(name='Expression z-score', bottom_annotation=HeatmapAnnotation(df=as.data.frame(patients_info)), column_title=paste('Clusters ', paste0(cids, collapse=',')), show_row_names=FALSE, show_row_dend=FALSE, right_annotation=rowAnnotation(df=fcl %>% column_to_rownames('Ensembl'))) )
}

#' Convert a two columns tibble to a named vector
#'
#' Wrapper for setNames
t2v <- function(tt) {
    nn = names(tt)
    if (length(nn) < 2) { stop("Wrong dimensions for t2v, the tible should have 2 columns") }
    setNames(tt[[nn[2]]], tt[[nn[1]]])
}
