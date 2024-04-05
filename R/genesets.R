#' @import clusterProfiler
#' @import tidyverse
#' @import org.Hs.eg.db
#' @import Rgraphviz

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

#' Build a named vector with all node names, so that it is valid for Rgraphviz::makeNodeAttrs
#'
#' @param partial A named vector with node names as names and a value to set as value
#' @param graph The graph for which the attribute will be set
#' @param default_value, the value to set for the other nodes, uses the node name if default_value=='names' (useful for setting label).
#' @examples
#' \dontrun{ makeNodeAttrs(graph, fill=nodeElement(setNames(1='blue'), graph, 'red'), label=nodeElement(setNames(1='text'), graph, "names")) }
nodeElement <- function(partial, graph, default_value) {
    nn = names(nodeData(graph))
    if (default_value=="names") {
        nn = setNames(nn, nn)
    } else {
        nn = setNames(rep(default_value, length(nn)), nn)
    }
    partial = partial[names(partial)[names(partial)%in%names(nn)]]
    nn[names(partial)]=partial
    return(nn)
}

#' Compute genesets enrichment for click clusters and plot the graph of gene sets relationships
#'
#' @param clusters Gene cluster assignement data.frame with columns Ensembl
#' @param clusters_enrichment A list with the enriched gene sets for each cluster
#' @param top_sets Number of top genesets to use for each CLICK cluster
#' @param max_sets Maximum number of genesets to plot
#' @param lwd_scaling Max line width, effectively the factor used to scale the line width in relation to the proportion of overlap in the gene sets. Use a high value if low overlaps are expected
#' @param cluster_downscaling CLICK clusters are way larger than genesets
#' @examples
#' \donotrun{
#' clusters %>% group_by(Cluster) %>% group_map(function(xx,yy){genesetsOverview(xx) %>% mutate(Cluster=yy$Cluster) %>% suppressMessages}) -> c_riches # Clusters enrichment
#' clusterSetsGraph(clusters, c_riches)
#'}
clusterSetsGraph <- function(clusters, clusters_enrichment, top_sets=5, max_sets=50, lwd_scaling=10, node_label='Description', node_size='fixed') {
    clusters_enrichment %>% lapply(arrange, pvalue) %>% lapply(head, top_sets) %>% lapply(function(xx){xx %>% separate(GeneRatio, c('d','n'), '/', FALSE, TRUE) %>% mutate(prop=d/n) }) %>% lapply(mutate, Description=gsub('HALLMARK_','', Description)) -> clusters_enrichment
    clusters_enrichment %>% bind_rows %>% filter(ID!='unattributed') %>% group_by(ID) %>% summarise(pvalue=min(pvalue)) %>% ungroup %>% arrange(pvalue) %>% head(n=max_sets) %>% distinct(ID) %>% pull -> genesets # Extract the most representative gene sets for the graph
    #missing_gs = genesets[!genesets %in% names(PATH2TERM)]
    #if (length(missing_gs) > 0) { warning(paste("Gene sets", paste(missing_gs, collapse=','), "not present in the PATH2TERM internal database")) }

    point2inch=72 # 1 inch (idth unit) is 72 pt (font unit)
    # TODO fill nodes with piecharts the number of genes attributed to a geneset
    # Plot cluster-geneset graph from enrichment results (edges are proportion of the gene set covered by the cluster list)
    clusters_enrichment %>% bind_rows %>% filter(ID %in% genesets) %>% select(`node_label`, Cluster, prop) %>% rename(from=`node_label`, to='Cluster', weight='prop') -> graph_edges
    #cgsedges %>% mutate(name=gsub('\\|','~',ename)) %>% group_by(name) %>% {setNames(group_split(.), group_keys(.)$name)} %>% {.[names(.) %in% edgeNames(cgsgraph)]}
    cgsgraph = graph_edges %>% select(from, to) %>% as.matrix %>% ftM2graphNEL(edgemode='undirected')
    clusters_enrichment %>% bind_rows %>% distinct(ID, Description) %>% filter(ID %in% names(nodeData(cgsgraph))) %>% {setNames(.[[node_label]], .$ID)} %>% nodeElement(cgsgraph, "names") %>% sapply(function(tt){gsub('\\n', '\\\\n', str_wrap(tt, 20))}) -> bnames
    if (node_size=='scaled') {
        gssizes = clusters_enrichment %>% bind_rows %>% group_by(ID) %>% summarise(n=mean(n)) %>% t2v %>% c(clusters %>% group_by(Cluster) %>% mutate(Cluster=gsub('^0', ' ', sprintf('%02d', Cluster))) %>% tally  %>% t2v) %>% {log10(.)/3} %>% round(1) %>% nodeElement(cgsgraph, 0.2)
    } else {
        gssizes = nodeElement(c(), cgsgraph, 1) # Fixed node size is more stable, 1 works well for pdfs
    }
    fontsize = min(gssizes)*72/2
    gspensizes = clusters_enrichment %>% bind_rows %>% group_by(ID) %>% summarise(n=mean(n)) %>% t2v %>% c(clusters %>% group_by(Cluster) %>% mutate(Cluster=gsub('^0', ' ', sprintf('%02d', Cluster))) %>% tally  %>% t2v) %>% nodeElement(cgsgraph, 0.2) %>% {round(./min(.)*1+0.2, 1)}
    cgsedges = bind_rows(graph_edges %>% mutate(ename=paste(from, gsub('^0', ' ', sprintf('%02d', to)), sep='~')), graph_edges %>% mutate(ename=paste(gsub('^0', ' ', sprintf('%02d', to)), from, sep='~'))) %>% filter(ename %in% edgeNames(cgsgraph)) %>% {list(
        lwd=setNames(ceiling(lwd_scaling*.$weight), .$ename)[edgeNames(cgsgraph)],
        fontsize=setNames(rep(fontsize*0.75, length(.$ename)), .$ename),
        #len=setNames(rep(max(gssizes), length(.$ename)), .$ename),
        minlen=setNames(rep(0.2, length(.$ename)), .$ename),
        label=setNames(sprintf("% 2d%%", round(100*.$weight)), .$ename)
    )} # Bug in Rgraphviz, the order of the edge weights must match the order by the graph, the naming is ignored (but used properly for labels)
    cgsnodes = makeNodeAttrs(cgsgraph,
        fillcolor=c("#8080DD","#DD8888")[1+sapply(names(nodeData(cgsgraph)), function(xx){grepl('^[ 0-9]',xx)})], # Different colors for genesets and clusters
        fontsize=fontsize,
        width=gssizes, height=gssizes/3,
        label=bnames,
        fixedsize=T )
    gsorc = list(list(graph=subGraph(names(nodeData(cgsgraph))[sapply(names(nodeData(cgsgraph)), function(xx){grepl('^[ 0-9]',xx)})], cgsgraph), cluster=FALSE, attrs=c(rank="sink")) ) # Bipartition geneset vs cluster (separated by the fact that cluster names are just numbers)
    plot(cgsgraph, 'neato', attrs=list(graph=list(rankdir="LR", rank="")), nodeAttrs=cgsnodes, subGList=gsorc, edgeAttrs=cgsedges)
}

