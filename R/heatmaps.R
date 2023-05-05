#' Heatmap with stars for extra information
#'
#' Plot heatmap with stars in the tile to encode for another information layer (significance, effect size, etc)
#' The range and threshold defaults consider the star values are log10-pvalues,
#' @param hm_values Values for the colors of the heatmap tiles
#' @param star_values Values for the text (or number of stars) of the heatmap tiles. Should be all >=0 or all <=0.
#' @param encoding One of c('size', 'count', 'multi_count') for whether the value in star_values should be encoded by the size of one star or a number of stars on one or multiple lines
#' @param return_list If FALSE the heatmap is plotted with a legend, if TRUE returns a list with fields 'object' for the Heatmap object and 'annotation_legend_list' for the annotation corresponding to the encoding.
#' @param range The range of values that should be covered by the stars.
#' @param hm_type Name of the data type represented in the heatmap, used to name the legend.
#' @param star_type Name of the data type represented by the stars, used to name the legend.
#' @param threshold Threshold under which the stars are not plotted.
#' @export
#' @rdname starHeatmap
starHeatmap2 <- function(hm_values, star_values, encoding="size", return_list=FALSE, range=-5:-2, hm_type='Value', star_type='log10(p-value)', threshold=2, ...) {
    star_values = star_values %>% .[nrow(.):1,,drop=FALSE] %>% as.matrix # Ensure star_values is a matrix
    min_hm = quantile(as.matrix(hm_values[!is.infinite(as.matrix(hm_values))]), 0.05)
    max_hm = quantile(as.matrix(hm_values[!is.infinite(as.matrix(hm_values))]), 0.95)
    if (min_hm >= 0) { hm_colors = circlize::colorRamp2(c(0, max_hm), c("white","red")) }
    else if (max_hm <= 0) { hm_colors = circlize::colorRamp2(c(min_hm, 0), c("blue","white")) }
    else { hm_colors = circlize::colorRamp2(c(min_hm, 0, max_hm), c("blue","white","red")) }
    if (sign(range[1]) != sign(range[length(range)])) { stop('The range of starHeatmap2 must be of consistent sign') }
    star_values[is.infinite(star_values)] = sign(range[1]) * max(abs(range))
    #star_values[star_values > max(abs(range))] = sign(range[1]) * max(abs(range))
    if (encoding == 'size') {
        out = list(
            object = hm_values %>%
                Heatmap(name=hm_type, cell_fun=function(j, i, x, y, width, height, fill) {
                    if(!is.na(star_values[i,j]) && abs(star_values[i,j]) > threshold) { grid.text("*", x, y-unit(0.2, "lines"), gp=gpar(fontsize=30*abs(star_values[i,j])/max(abs(range)))) }
                }, show_row_dend=FALSE, show_column_dend=FALSE, row_names_side="left", col=hm_colors, border=TRUE, column_title_side="bottom", ...),
               annotation_legend_list = list(Legend(title=star_type, at=range, labels=range, pch="*", legend_gp=gpar(fontsize=30/max(abs(range))*abs(range)), type="points", background="white", grid_height=unit(2, "lines")))
        )
    } else if (encoding == 'count') {
        out = list(
            object = hm_values %>%
                Heatmap(name=hm_type, cell_fun=function(j, i, x, y, width, height, fill) {
                    if(!is.na(star_values[i,j]) && abs(star_values[i,j]) > threshold) { grid.text(paste0(rep("*", floor(abs(star_values[i,j])/max(abs(range))*length(range))), collapse=""), x, y-unit(0.2, "lines"), gp=gpar(fontsize=20)) }
                }, show_row_dend=FALSE, show_column_dend=FALSE, row_names_side="left", col=hm_colors, border=TRUE, column_title_side="bottom", ...),
                annotation_legend_list = list(Legend(title=star_type, at=range, labels=range, pch=floor(abs(range)/max(abs(range))*length(range)) %>% sapply(function(nn){paste0(rep("*", nn), collapse="")}), type="points"))
        )
    } else if (encoding == 'multi_count') {
        out = list(
            object = hm_values %>%
                Heatmap(name=hm_type, cell_fun=function(j, i, x, y, width, height, fill) {
                    if(!is.na(star_values[i,j]) && abs(star_values[i,j]) > threshold) { grid.text(gsub("(.{3})", "\\1\n", paste0(rep("*", floor(abs(star_values[i,j])/max(abs(range))*length(range))), collapse="")), x, y-unit(0.2, "lines"), gp=gpar(fontsize=20)) }
                }, show_row_dend=FALSE, show_column_dend=FALSE, row_names_side="left", col=hm_colors, border=TRUE, column_title_side="bottom", ...),
                annotation_legend_list = list(Legend(title=star_type, at=range, labels=range, pch=floor(abs(range)/max(abs(range))*5) %>% sapply(function(nn){paste0(rep("*", nn), collapse="")}), type="points"))
            )
    }
    if (return_list) {
        return(out)
    } else {
        #return(out$object+out$annotation_legend_list[[1]])
        return(do.call(draw, out))
    }
}
