# Recalculate the estimates and SEs of PVEs and enrichments from MPH VC estimates
# vcfile: the .mq.vc.csv file produced by `mph --minque`.
# crossprod: the crossproduct of SNP incidence matrix and SNP weighting matrix. Its row names and columns should match the annotation categories of interest and the rows of vcfile, respectively.
# nsnps: the total number of SNPs. If `NA`, it is set to the first GRM's number of SNPs in vcfile.
# annot.size: a list of the number of SNPs in each annotation category listed in the row names of crossprod. It can be computed from the column-wise sum of the corresponding SNP incidence matrix. If `NA`, it is set to the `m` column of vcfile.
recalculate_enrichments <- function(vcfile, crossprod, nsnps=NA, annot.size=NA) {
    if(is.na(vcfile)) {
        stop("vcfile must be specifiled.")
    }
    if(length(crossprod) == 1 && is.na(crossprod)) {
        stop("crossprod must be specified.")
    }

    mq = read.csv(vcfile)
    mq = mq[-nrow(mq),]
    if(ncol(crossprod) > nrow(mq)) {
        stop(paste("ncol(crossprod) should equal the number of genetic VCs in", vcfile))
    }
    mq = mq[1:ncol(crossprod),]

    # PVE estimate
    var = crossprod %*% as.matrix(mq[,9:(8+nrow(mq))]) %*% t(crossprod)
    pve.est = crossprod %*% mq$enrichment / sum(mq$m)
    pve.se = sqrt(diag(var)) / sum(mq$m) 

    # enrichment estimate
    if(is.na(nsnps)) {
        nsnps = mq$m[1]
    }
    if(length(annot.size) == 1 && is.na(annot.size)) {
        annot.size = mq$m
    }
    if(length(annot.size) != nrow(crossprod)) {
        stop("annot.size must correspond to the rows of crossprod.")
    }
    e.est = pve.est / (annot.size/nsnps)
    e.se = pve.se / (annot.size/nsnps)

    out = cbind(annot.size/nsnps, pve.est, pve.se, e.est, e.se)
    rownames(out) = rownames(crossprod)
    colnames(out) = c("prop", "pve", "pve.se", "enrichment", "enrichment.se")
    return(out)
}

