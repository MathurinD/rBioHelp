# Functions to display enrichment results from clusterProfiler
# Can be GO or KEGG output

# @import forecast

#' @title Plot KEGG result
#'
#' Compute and plots the rich factor for the enriched pathways or terms.
#' @param kegg_result Output of enrichKEGG or enrichGO
#' @param max_terms The maximum number of terms to display on the plot
#' @param ref_enrichment Output of enrichKEGG or enrichGO for the universe set
#' @param ref_rep Representation of the reference enrichment, either the richFactor is computed using the reference set count as the denominator ('one'), or the richFactor of the reference set is displayed on the plot ('plot')
#' @value The input with an extra column richFactor
#' @export
plotKEGG <- function(kegg_result, max_terms=10, ylab='KEGG pathway', title=paste('Enrichment score of', ylab), ref_enrichment='', ref_rep='one') {
    kegg_result %>% mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) -> kegg_tibble
    if (nrow(kegg_tibble)==0) {
        message('No enriched terms, skipping...')
    } else {
        gg = kegg_tibble %>% as_tibble %>% arrange(desc(richFactor)) %>% head(max_terms) %>%
              ggplot(aes(richFactor, fct_reorder(Description, richFactor))) +
              geom_segment(aes(xend=0, yend = Description)) +
              geom_point(aes(color=p.adjust, size = Count))+
              ggtitle(title) + ylab(ylab) +
              scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),guide=guide_colorbar(reverse=TRUE, order=1))
        if (class(ref_enrichment) %in% c('enrichResult','tibble')) {
            if (ref_rep == 'one') {
                gg = kegg_tibble %>% as_tibble %>%
                    left_join(select(as_tibble(ref_enrichment), Description, Count), by='Description') %>% mutate(richFactor=Count.x/Count.y) %>% # Compute enrichment relative to the reference set
                    rename(Count=Count.x) %>%
                    arrange(desc(richFactor)) %>% head(max_terms) %>% # Select most enriched terms
                    ggplot(aes(richFactor, fct_reorder(Description, richFactor))) +
                    geom_segment(aes(xend=0, yend = Description)) +
                    geom_point(aes(color=p.adjust, size = Count))+
                    ggtitle(title) + ylab(ylab) + xlab('richFactor relative to the reference') +
                    #scale_size_continuous(breaks=ceiling(exp(seq(0, 4.5, 0.5))), limits=c(0,1000)) +
                    scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),guide=guide_colorbar(reverse=TRUE, order=1))
            } else if (ref_rep == 'plot') {
                ref_enrichment = ref_enrichment %>% mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
                gg = gg + geom_point(aes(color=p.adjust, size = Count), data=filter(ref_enrichment, ID %in% gg$data$ID) %>% as_tibble)
            } else {
                stop('Invalid "ref_rep" argument in plotKEGG')
            }
        }
        print(gg)
    }
    invisible(kegg_tibble)
}

