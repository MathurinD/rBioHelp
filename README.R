#'---
#'title: "rBioHelp"
#'author: "Mathurin Dorel"
#'output:
#'  github_document:
#'      toc: true
#'---

#' A set of functions that implement common operations used when analysing human biological datasets.
#' Aggregates several packages I like and versions of datasets I use.
#'
#+ setup
    library(rBioHelp)
    library(tidyverse)
    library(ComplexHeatmap)
    library(org.Hs.eg.db)
    # Toy data for illustration
    geneslist = c('MAP2K2','MAP2K1','MAPK1','MAPK3','PIK3CA','PIK3CG', 'MTOR','TSC1','TSC2','EGFR','GRB2','RAF1','BRAF','KRAS','PTPN11') # Example signalling
    geneslist2 = c('MCM5','PCNA','TYMS','FEN1','MCM7','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1','UHRF1','CENPU','HELLS','RFC2','POLR1B','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1','CHAF1B','MRPL36','E2F8')
    clusters = tibble(Cluster = c(rep(1, length(geneslist)), rep(2, length(geneslist2))), Genes=c(geneslist, geneslist2) ) %>% mutate(Entrez=mapIds(org.Hs.eg.db, Genes, 'ENTREZID','SYMBOL')) %>% mutate(Ensembl=mapIds(org.Hs.eg.db, Genes, 'ENSEMBL','SYMBOL')) %>% mutate(Uniprot=mapIds(org.Hs.eg.db, Genes, 'UNIPROT','SYMBOL'))

    zscores = matrix(c(rnorm(5*58, -2), rnorm(5*58, 2)), nrow=nrow(clusters), dimnames=list(clusters$Ensembl, paste0('S', 1:10))) %>% as_tibble(rownames='g') # Pseudo z-score with two groups
    patients_info = tibble(Name=paste0('S', 1:10), Status=c(rep('Sick', 5), rep('Healthy', 5))) %>% column_to_rownames('Name')

#' The package provides several helper functions to compare gene clusters (e.g DEGs after different treatments, coregulated modules, etc).
#+ genesets
    # Simply plot the expression/zscores with sample annotation
    compareClusters(zscores, clusters, 1:2) # TODO have patients_info as argument

    # Alternative to clusterProfiler::dotplot
    wpenrich = enrichWP(clusters$Entrez, 'Homo sapiens')
    plotEnrichment(wpenrich)

    # Compute enrichment using multiple gene set sources for each gene group
    clusters_enrichment = suppressMessages(genesetsOverview(clusters))
    clusters_enrichment %>% tail %>% pretty_rich()

#' A graph between gene sets helps to better visualise the coverage by the various gene sets enriched.
#+ genesets_graph
    clusters %>% group_by(Cluster) %>% group_map(function(xx,yy){genesetsOverview(xx) %>% mutate(Cluster=yy$Cluster) %>% suppressMessages}) -> c_riches # Clusters enrichment
    clusterSetsGraph(clusters, c_riches, top_sets=5, max_sets=50, graphic_scaling=0.8)
    clusterSetsGraph(clusters, c_riches %>% lapply(filter, grepl('HALLMARK', ID)), top_sets=10, max_sets=50, graphic_scaling=0.8) # Restrict the plot to hallmark gene sets


#+ hidden, echo=FALSE
#genesets.R:genesetsOverview <- function(clusters, ...)  {
#genesets.R:genesetsFromCluster <- function(clusters, selection, ...) { genesetsOverview(clusters %>% filter(Cluster %in% selection), ...) }
#genesets.R:pretty_rich <- function(rich){
#genesets.R:compareClusters <- function(zscores, clusters, cids) {
#genesets.R:t2v <- function(tt) {
#genesets.R:nodeElement <- function(partial, graph, default_value) {
#genesets.R:clusterSetsGraph <- function(clusters, clusters_enrichment, top_sets=5, max_sets=50, lwd_scaling=10, node_label='Description', node_size='fixed') {
#go.R:findGenericGO <- function(goset, generic_set = GOBP1) {
#go.R:goToTerm <- function(goid) {
#heatmaps.R:starHeatmap2 <- function(hm_values, star_values, encoding="size", return_list=FALSE, range=-5:-2, hm_type='Value', star_type='log10(p-value)', threshold=2, ...) {
#kegg.R:plotKEGG <- function(kegg_result, max_terms=10, ylab='KEGG pathway', title=paste('Enrichment score of', ylab), ref_enrichment='', ref_rep='one') {
#ksea.R:readGMT <- function(gmt_file, format='hugo.') {
#ksea.R:ksea <- function(sig, signame='default') { # Implement from the methods of Wiredja 2017
#pca.R:setMethod(plotPCA, signature(object='prcomp'), function(object, components=1:2) {
#venn.R:combiSet <- function(vec, low_order=TRUE) {
#venn.R:customVenn <- function(vennlist, main='') {

#' Improve wikiprofiler plots, allowing control over the color gradient and plotting multiple datasets on the same pathway
#' TODO list: sample name legend
#+ wikipathway
    
#wpplot.R:WPplot <- function(ftable, pathwayid='WP4172', high="red", low="blue", legend = TRUE, legend_x = 0.001, legend_y = 0.94, strip_out=0.01, fixed_range=F, lim=5) {
#wpplot.R:WPmultiplot <- function(ftable, pathwayid='WP4172', multi='Cell_line', high="red", low="blue", legend = TRUE, legend_x = 0.001, legend_y = 0.94, strip_out=0.01, fixed_range=F, lim=5) {
#wpplot.R:#' Helper function to change the background for a gene
#wpplot.R:replace_bg <- function(svg, gene, positions, colors) {
#wpplot.R:WPbgfill <- function(pp, value, high="red", low="blue", legend = TRUE, legend_x = 0.001, legend_y = 0.94, strip_out=0.01, fixed_range=F, lim=5) {

#Eventually will automatically update datasets on install (and keep old versions for reproducibility) once I find the time to code that.

#' Developped by Mathurin Dorel
