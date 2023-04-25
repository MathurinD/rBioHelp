# Set of functions to plot PCA results in tidyverse, 
# Work with the prcomp objects from the stats package

#' @import BiocGenerics

#' @export
# Should use the generic function defined by BiocGenerics
setMethod(plotPCA, signature(object='prcomp'), function(object, components=1:2) {
    sdev = signif(object$sdev^2/sum(object$sdev^2), 2)*100
    xpc = paste0('PC', components[1])
    ypc = paste0('PC', components[2])
    object$x %>% as_tibble(rownames='Samples') %>% ggplot(aes(get(xpc),get(ypc), label=Samples)) + geom_text() + xlab(paste0(xpc, ' (', sdev[components[1]], '%)')) + ylab(paste0(ypc, ' (', sdev[components[2]], '%)'))
    #%>% separate(Samples, into=c(''))
})

