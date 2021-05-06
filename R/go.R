#' @import tidyverse
#' @import GO.db

#slimgobp = read_csv("generic_go_slim.csv", col_names="GO")
# as.list(GOBPCHILDREN)["GO:0008150"] %>% unlist %>% as_tibble %>% rename(value="GO") -> GOBP1 # "biological process children"
# GOBP1$GO %>% as.list(GOBPCHILDREN)[.] %>% unlist %>% unique %>% as_tibble %>% rename(value="GO") -> GOBP2
# GOBP2$GO %>% as.list(GOBPCHILDREN)[.] %>% unlist %>% unique %>% as_tibble %>% rename(value="GO") -> GOBP3

#' @title Reduce a set of GO terms to higher order terms
#'
#' Reduce a set of GO terms to higher order terms defined in "slims" from http://geneontology.org/docs/download-ontology/#subsets. See http://geneontology.org/docs/go-subset-guide/ for more details
#' @param goset A semicolon-separated (";") set of GO terms
#' @export
findGenericGO <- function(goset, generic_set = GOBP1) {
    if (length(goset) > 1) {
        return( goset %>% mclapply(findGenericGO) %>% unlist)
    }
    if (is.na(goset) | length(goset)==0 | goset == "" | goset == "all") return("")
    goset %>% str_split(";") %>% unlist %>% unique %>% {.[. %in% generic_set$GO]} -> output
    if (length(output) == 0) {
        output = goset %>% str_split(";") %>% unlist %>% unique %>% {as.list(GOBPPARENTS)[.]} %>% unlist %>% {.[grepl("isa", names(.))]} %>% unique %>% paste0(collapse=";") %>% findGenericGO
    }
    return(paste(output, collapse=";"))
}

#' Convert GO ID to GO Term
#' @param Vector of GO IDs (in the form "GO:0123456")
#' @export
goToTerm <- function(goid) {
    return(Term(GOTERM)[goid])
}
