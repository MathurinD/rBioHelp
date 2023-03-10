#' Reads gmt file
#'
#' Read GMT file in various format
readGMT <- function(gmt_file, format='hugo.') {
    outgmt = list()
    for (ll in readLines(gmt_file)) {
        line = unlist(str_split(ll, '\t'))
        hugo = unlist(str_split(line[2], '\\|'))[-1] %>% sapply(function(xx){ unlist(strsplit(xx, ';'))[1] })
        if (format == 'all') {
            outgmt[[line[1]]] = line[-1]
        } else if (format == 'hugo.') {
            outgmt[[line[1]]] = hugo %>% str_replace('_([TYS])', '.1\\1.') %>% str_split(':') %>% lapply(function(dd){dd[1]}) %>% unlist
        } else if (format == 'hugo_') {
            outgmt[[line[1]]] = hugo
        } else if (format == 'peptide') {
            outgmt[[line[1]]] = line[-c(1,2)]
        }
    }
    return(outgmt)
}

# PTMsigDB
#ptmsigdb = readGMT('PTMsigDB/ptm.sig.db.all.sitegrpid.human.v1.9.0.gmt')
# Only PhosphoSitePlus kinases
#kinase_sig = ptmsigdb[grepl('KINASE', names(ptmsigdb))]

#' KSEA
#'
#' Kinase Substrate Enrichment Analysis
#' @param sig One signature
ksea <- function(sig, signame='default') { # Implement from the methods of Wiredja 2017
    lfc_phospho %>% toGenes %>% filter(ID %in% sig) %>% summarise_at(vars(-ID), mean, na.rm=TRUE) -> mean_substrate
    lfc_phospho %>% toGenes %>% summarise_at(vars(-ID), mean, na.rm=TRUE) -> mean_all
    lfc_phospho %>% toGenes %>% summarise_at(vars(-ID), sd, na.rm=TRUE) -> sd_all
    lfc_phospho %>% toGenes %>% filter(ID %in% sig) %>% nrow -> matching_sites
    kzscore = sqrt(matching_sites)*(mean_substrate - mean_all)/sd_all
    samples_conditions %>% left_join(kzscore %>% t %>% as_tibble(rownames='value')) %>% group_by(Condition=gsub(".[ABC]$", "", value)) %>% summarise(KScore=mean(V1)) %>% column_to_rownames('Condition') %>% mutate(pvalue = pnorm(KScore)) %>% mutate(pvalue=sapply(pvalue, function(pp){min(pp, 1-pp)})) %>% mutate(Set=signame) %>% as_tibble(rownames='Condition')
}
