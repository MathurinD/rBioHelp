# Functions to display enrichment results from clusterProfiler
# Can be GO or KEGG output

#' @import forecats

#' Plot KEGG result
#'
#' Compute and plots the rich factor for the enriched pathways or terms.
#' @param kegg_result Output of enrichKEGG or enrichGO
#' @param max_terms The maximum number of terms to display on the plot
#' @value The input with an extra column richFactor
plotKEGG <- function(kegg_result, max_terms=10) {
    kegg_result %>% mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) -> kegg_tibble
    if (nrow(kegg_tibble)==0) {
        message('No enriched terms, skipping...')
    } else {
        print( kegg_tibble %>% as_tibble %>% arrange(desc(richFactor)) %>% head(max_terms) %>%
              ggplot(aes(richFactor, fct_reorder(Description, richFactor))) +
              geom_segment(aes(xend=0, yend = Description)) +
              geom_point(aes(color=p.adjust, size = Count))+
              scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),guide=guide_colorbar(reverse=TRUE, order=1))
        )
    }
    invisible(kegg_tibble)
}

