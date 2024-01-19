#' @title Wrapper for wikiprofiler wpplot and wp_bgfill
#' @import wikiprofiler
#' @import org.Hs.eg.db
#' @param ftable A data table with columns 'Gene' for gene names and 'beta' for the metric to plot
#' @export
#' @rdname wpplot
WPplot <- function(ftable, pathwayid='WP4172', high="red", low="blue", legend = TRUE, legend_x = 0.001, legend_y = 0.94, strip_out=0.01, fixed_range=F, lim=5) {
    wpplot(pathwayid) %>% WPbgfill(ftable %>% select(Gene, beta) %>% column_to_rownames('Gene') %>% {setNames(t(.), rownames(.))} %>% sort)
}

#' @param ftable Either a tibble in long format with columns 'multi', Gene and beta containing the conditions, the gene symbols and the value to plot respectively; or a matrix with the conditions as column names and the genes as row names.
#' @param multi The name of the colum with the conditions
#' @export
#' @rdname wpplot
WPmultiplot <- function(ftable, pathwayid='WP4172', multi='Cell_line', high="red", low="blue", legend = TRUE, legend_x = 0.001, legend_y = 0.94, strip_out=0.01, fixed_range=F, lim=5) {
    if (multi %in% colnames(ftable)) {
        wpplot(pathwayid) %>% WPbgfill(ftable %>% rename(multi=!!multi) %>% group_by(Gene, multi) %>% summarise(beta=mean(beta)) %>% ungroup %>% select(Gene, beta, multi) %>% pivot_wider(names_from=multi, values_from=beta) %>% column_to_rownames('Gene'))
    } else {
        message(paste0('No multi ("', multi, '") column in the data, assuming it is already in wide matrix format'))
        wpplot(pathwayid) %>% WPbgfill(ftable)
    }
}

    
#' Helper function to change the background for a gene
# https://stackoverflow.com/questions/61099149/svg-fill-color-not-working-with-hex-colors
# https://stackoverflow.com/questions/29335354/filling-an-svg-path-with-multiple-colors
# https://stackoverflow.com/questions/14051351/svg-gradient-using-css
replace_bg <- function(svg, gene, positions, colors) {
    if (is.null(positions) || all(is.na(positions)) || length(positions) == 0) {
        return(svg)
    }
    # Define the 'gradient'
    ii = 1
    fill = paste0('<defs><linearGradient id="color', gene,'">')
    while (ii <= length(colors)){
        fill = paste0(fill,
                     '<stop offset="', signif((ii-1)/length(colors), 2)*100,'%" style="stop-color:', substr(colors[ii],1,7), ';stop-opacity:1"/>',
                     '<stop offset="', signif(ii/length(colors), 2)*100,'%" style="stop-color:', substr(colors[ii],1,7), ';stop-opacity:1"/>'
               )
        ii=ii+1
    }
    fill = paste0(fill, '</linearGradient></defs>')
    svg[grep('<svg ', svg)+1] = sub('><', paste0('>', fill, "<"), svg[grep('<svg ', svg)+1])
    # Set all occurences of the gene to the gradient
    for (position in positions) {
        if (is.na(position)) next
        jj <- rev(grep("<g", svg[1:position]))[1]
        svg[jj] <- sub(
          "fill:.+;.+", # Look for the last fill before the gene name and change the background color
          paste("fill:url(#color", gene, "); text-rendering:geometricPrecision; stroke:white;\"", sep = ""),
        svg[jj] )
    }
    
    return(svg)
}


#' @title Fill the background of gene with color according to amount of gene expression.
#' @description Generate a color array.Fill the gene then generate the legend.
#' @param pp wpplot object as returned by wikiprofiler::wpplot or WPplot
#' @param value value is the amount of expression. Can be a matrix or a vector.
#' @param low The color of lowest gene.
#' @param high The color of highest gene.
#' @param legend Whether you need legend.
#' @param legend_x horizontal position of the legend
#' @param legend_y vertical position of the legend
#' @return A 'wpplot' object
# @import org.Hs.eg.db
# @import BiocGenerics
#' @export
WPbgfill <- function(pp, value, high="red", low="blue", legend = TRUE, legend_x = 0.001, legend_y = 0.94, strip_out=0.01, fixed_range=F, lim=5) {
#high="red"; low="blue"; legend = TRUE; legend_x = 0.001; legend_y = 0.94; strip_out=0.01; fixed_range=F; lim=5
#depmap_dependencies %>% filter(DepmapModelType=='BLL') %>% select(Gene, beta, Cell_line) %>% pivot_wider(names_from=Cell_line, values_from=beta) %>% column_to_rownames('Gene') -> value
    if (is.null(dim(value))) { return(wp_bgfill(pp, value, high, low, legend, legend_x, legend_y)) } # Use the basic version for vectors.
    # TODO Add legend for the samples plotted
    # TODO Add an option to rank the values, usefull to compare distribution
    
    if(legend_x < 0 || legend_x > 1 || legend_y < 0 || legend_y > 1){
      message('Parameters legend_x and legend_y must be numbers between 0 to 1!')
    }
    # Filter for genes in the pathway
    SYMBOLS <- sub('\\s+', '', sub('>', '', sub('</text', '', pp$svg[grep('</text', pp$svg)])))
    SYMBOLS <- SYMBOLS[is.na(suppressWarnings(as.numeric(SYMBOLS)))]
    if(!any(rownames(value) %in% SYMBOLS)){
        if (all(grepl('^ENSG', rownames(value)))) { # Auto detection Ensembl Ids
            message("Ensembl IDs detected, automatic conversion to SYMBOL and averaging over identical names.")
            value = value %>% as_tibble(rownames='Gene') %>% mutate(Gene=mapIds(org.Hs.eg.db, Gene, 'SYMBOL', 'ENSEMBL')) %>% group_by(Gene) %>% summarise_all(mean) %>% filter(!is.na(Gene)) %>% column_to_rownames('Gene')
        } else if(all(grep('^[A-Z0-9]{6,10}$', rownames(value)))){ # https://www.uniprot.org/help/accession_numbers
            message("Uniprot IDs detected, automatic conversion to SYMBOL and averaging over identical names.")
            value = value %>% as_tibble(rownames='Gene') %>% mutate(Gene=mapIds(org.Hs.eg.db, Gene, 'SYMBOL', 'UNIPROT')) %>% group_by(Gene) %>% summarise_all(mean) %>% filter(!is.na(Gene)) %>% column_to_rownames('Gene')
        } else {
            message("Please make sure the input gene ID type is 'SYMBOL'.")
            return(pp)
        }
    }
    value <- value[rownames(value) %in% SYMBOLS,]
    #value <- sort(value[names(value) %in% SYMBOLS])

    # Set the color scale
    if (strip_out >= 0.5) { stop("Cannot strip more than 50% of the data to generate the color scale") }
    if (!is.numeric(lim)) {stop("lim should be numeric")}
    lim = abs(lim)
    if (!fixed_range) {
        lowLim = max(-lim,quantile(value,strip_out, na.rm=T))
        upLim = min(lim,quantile(value,1-strip_out, na.rm=T))
    } else {
        lowLim = -lim
        upLim = lim
    }
    value[value < lowLim] = lowLim
    value[value > upLim] = upLim

    color_palette = circlize::colorRamp2(c(-1, 0, 1), c(low,'white',high))
    gene_colors = value %>% apply(2, color_palette)
    
    # Add the colors to the genes
    genes <- rownames(value)
    for (ii in seq_along(genes)) {
      pos <- grep(paste0('>',genes[ii],'<'), pp$svg)
      pp$svg <- replace_bg(pp$svg, genes[ii], pos, gene_colors[ii,])
    }
    
    svg_width <- as.numeric(strsplit(strsplit(pp$svg[4], 'width=\"')[[1]][2], '\" height=\"')[[1]][1])
    svg_height <- as.numeric(strsplit(strsplit(strsplit(pp$svg[4], 'width=\"')[[1]][2], '\" height=\"')[[1]][2], '\"')[[1]][1])
    
    incrementX <- svg_width * legend_x
    incrementY <- svg_height * (1 - legend_y)
    if(incrementX > svg_width - 48)
      incrementX <- svg_width - 48
    
    if(incrementY > svg_height - 122){
      incrementY <- svg_height - 122
    }else if(incrementY < 3)
      incrementY <- 3
    
    #### Breaks here
    prettyvalue = pretty(as.matrix(value), 4)
    textele <- rev(prettyvalue) 
    legendX <- 0 + incrementX
    legendY <- 0 + incrementY
    textX <- 40 + incrementX
    textY <- seq(from = 5,to = 120,length.out = length(textele)) + incrementY
    scalelineX <- 27 + incrementX
    scalelineY <- seq(from = 2,to = 118,length.out = length(textele)) +incrementY
    
    if(legend){
      zero_scale_line <- wikiprofiler:::find_zero_scale(as.matrix(value))
      proportion <- seq(from = 2,to = 118,length.out = length(textele)) / 120
      proportion <- proportion[length(which(prettyvalue >= zero_scale_line))]
      if(max(prettyvalue) == 0){
        proportion <- '0%'
      }
      if(min(prettyvalue) == 0){
        proportion <- '100%'
      }
      temp<-grep("</svg",pp$svg)
      pp$svg[temp]<-sub("</svg",paste("<defs><linearGradient id=\"grad1\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\"><stop offset=\"0%\" style=\"stop-color:",high,";stop-opacity:1\"></stop><stop offset=\"",proportion,"\" style=\"stop-color:","white",";stop-opacity:1\"></stop><stop offset=\"100%\" style=\"stop-color:",low,";stop-opacity:1\"></stop></linearGradient></defs><rect x=\"",legendX,"\" y=\"",legendY,"\" width =\"30\" height=\"120\" style=\"fill:url(#grad1 );stroke-width:0;stroke:black\"></rect></svg",sep = ""),pp$svg[temp])
      
      for (i in 1:length(prettyvalue)){
        temp<-grep("</svg",pp$svg)
        pp$svg[temp]<-sub("</svg",paste("<text x=\"",textX,"\" y=\"",textY[i],"\" style=\"font-size:10; fill:black; stroke:none\">",textele[i],"</text></svg",sep = ""),pp$svg[temp])
      }
      for (i in 1:length(prettyvalue)){
        temp<-grep("</svg",pp$svg)
        pp$svg[temp]<-sub("</svg",paste("<rect width=\"3\" height=\"1\" x=\"",scalelineX,"\" y=\"",scalelineY[i],"\" style=\"fill:white; stroke:none\"></rect></svg",sep = ""),pp$svg[temp])
      }
    }
    pp$geneExpr <- value
    return(pp)
}
