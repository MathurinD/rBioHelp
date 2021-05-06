############################################## venn.R ##################################################3
# Make euler plot rather than venn diagrams when possible

#' @import eulerr

#'
#' @title Generates combinations separated by '&'
#'
#' Generates all combinations or the elements in a vector.
#' @param vec A vector with names
#' @param low_order Whether to return the low order combinations or only the full
combiSet <- function(vec, low_order=TRUE) {
    if (length(vec)==1) { return(vec) }
    output = c()
    for (vv in seq_along(vec)) {
        others = vec[-seq(vv)] # The sets before have already been considered
        if (length(others) > 0) {
            output = c(output, vec[vv], paste0(vec[vv],"&",combiSet(others)))
        } else {
            output = c(output, vec[vv])
        }
    }
    return(output %>% strsplit('&') %>% sapply(length) %>% order(decreasing=TRUE) %>% {output[.]}) # Reorder so the low orders are first, usefull for customVenn
}

# Helper function to get the entries exclusive to the set
exclude <- function(set1, set2) { return(set1[!set1 %in% set2]) }

#' Euler diagram
#'
#' Make an euler diagram, always works exactly up to 3 sets, then has to approximate so use equisurface venn diagram
#' @param vennlist A named list where each entry is a vector of the elements of the set
#' @export
customVenn <- function(vennlist) {
    if (is.null(names(vennlist))) {
        conditions = paste0('set', 1:length(vennlist))
    } else {
        conditions = combiSet(names(vennlist))
    }
    clengths = conditions %>% strsplit('&') %>% sapply(length)
    cat = numeric(length(conditions))
    sets = list(length(conditions))
    sets = lapply(conditions, function(xx){ character(0) })
    names(sets) = conditions
    names(cat) = conditions
    for (cc in conditions) {
        singles = cc %>% str_split("&") %>% unlist
        setsize = length(singles)
        # Find the intersecting entries
        inter = vennlist[[singles[1]]]
        for (ss in singles[-1]) {
            inter = intersect(inter, vennlist[[ss]])
        }
        # Remove the entries already counted in bigger intersection sets
        # Relies on the sets being build in decreasing size order
        for (ss in conditions[clengths > setsize]) {
            inter = exclude(inter, sets[[ss]])
        }
        sets[[cc]] = inter
        cat[cc] = length(inter)
    }
    if (length(vennlist) >= 4) {
        efit = venn(cat)
    } else {
        efit = euler(cat, shape="ellipse")
        # TODO check that the fit is appropriate or try again to find a better one, importantly no category should not be represented
    }
    values = efit$original.values
    plot(efit, counts=TRUE, font=1, cex=2, quantities=paste0(values, "\n(", signif(values*100/sum(values), 2), "%)"), type="percent")
}
