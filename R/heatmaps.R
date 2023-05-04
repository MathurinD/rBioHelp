#' Heatmap with stars for extra data
#'
#' Plot heatmap with stars in the tile to encode for another data layer (significance, effect size, etc)
#' @param hm_values Values for the colors of the heatmap tiles
#' @param data Values for the text (or number of stars) of the heatmap tiles. Should be all >=0 or all <=0.
#' @param encoding One of c('size', 'count', 'multi_count') for whether the value in data should be encoded by the size of one star or a number of stars on one or multiple lines
#' @param return_list If FALSE the heatmap is plotted with a legend, if TRUE returns a list with fields 'object' for the Heatmap object and 'annotation_legend_list' for the annotation corresponding to the encoding.
#' @param range The range of values that should be covered by the stars.
#' @param hm_type Name of the data type represented in the heatmap, used to name the legend.
#' @param star_type Name of the data type represented by the stars, used to name the legend.
#' @param threshold Threshold under which the stars are not plotted.
#' @export
#' @rdname starHeatmap
starHeatmap2 <- function(hm_values, data, encoding="size", return_list=FALSE, range=-5:-1, hm_type='Value', star_type='log p-value', threshold=0, ...) {
    data = data %>% .[nrow(.):1,,drop=FALSE] %>% as.matrix
    min_hm = quantile(as.matrix(hm_values[!is.infinite(as.matrix(hm_values))]), 0.05)
    max_hm = quantile(as.matrix(hm_values[!is.infinite(as.matrix(hm_values))]), 0.95)
    if (min_hm >= 0) { hm_colors = circlize::colorRamp2(c(0, max_hm), c("white","red")) }
    else if (max_hm <= 0) { hm_colors = circlize::colorRamp2(c(min_hm, 0), c("blue","white")) }
    else { hm_colors = circlize::colorRamp2(c(min_hm, 0, max_hm), c("blue","white","red")) }
    data[is.infinite(data)] = max(abs(range))
    if (encoding == 'size') {
        out = list(
            object = hm_values %>%
                Heatmap(name=hm_type, cell_fun=function(j, i, x, y, width, height, fill) {
                    if(!is.na(data[i,j]) && abs(data[i,j]) > threshold) { grid.text("*", x, y-unit(0.2, "lines"), gp=gpar(fontsize=30*data[i,j]/max(abs(range)))) }
                }, show_row_dend=FALSE, show_column_dend=FALSE, row_names_side="left", col=hm_colors, border=TRUE, column_title_side="bottom", ...),
               annotation_legend_list = list(Legend(title=star_type, at=range, labels=range, pch="*", legend_gp=gpar(fontsize=30/max(abs(range))*abs(range)), type="points", background="white", grid_height=unit(2, "lines")))
        )
    } else if (encoding == 'multi_count') {
        out = list(
            object = hm_values %>%
                Heatmap(name=hm_type, cell_fun=function(j, i, x, y, width, height, fill) {
                    if(!is.na(data[i,j]) && abs(data[i,j]) > threshold) { if(is.infinite(data[i,j])){data[i,j]=max(abs(range))}; grid.text(paste0(rep("*", floor(abs(data[i,j])/max(abs(range))*5)), collapse=""), x, y-unit(0.2, "lines"), gp=gpar(fontsize=20)) }
                }, show_row_dend=FALSE, show_column_dend=FALSE, row_names_side="left", col=hm_colors, border=TRUE, column_title_side="bottom", ...),
                annotation_legend_list = list(Legend(title=star_type, at=range, labels=range, pch=floor(abs(range)/max(abs(range))*5) %>% sapply(function(nn){paste0(rep("*", nn), collapse="")}), type="points"))
        )
    } else if (encoding == 'count') {
        out = list(
            object = hm_values %>%
                Heatmap(name=hm_type, cell_fun=function(j, i, x, y, width, height, fill) {
                    if(!is.na(data[i,j]) && abs(data[i,j]) > threshold) { if(is.infinite(data[i,j])){data[i,j]=max(abs(range))};  grid.text(gsub("(.{3})", "\\1\n", paste0(rep("*", floor(abs(data[i,j]/10))), collapse="")), x, y-unit(0.2, "lines"), gp=gpar(fontsize=20)) }
                }, show_row_dend=FALSE, show_column_dend=FALSE, row_names_side="left", col=hm_colors, border=TRUE, column_title_side="bottom", ...),
                annotation_legend_list = list(Legend(title=star_type, at=range, labels=range, pch=floor(abs(range)/max(abs(range))*5) %>% sapply(function(nn){paste0(rep("*", nn), collapse="")}), type="points"))
            )
    }
    if (return_list) {
        return(out)
    } else {
        return(do.call(draw, out))
    }
}

