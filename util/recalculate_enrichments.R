# Arguments
# vcfile: the .mq.vc.csv file produced by `mph --minque`.
# snpinfo: SNP info file. Igonored if crossprod is provided.
# crossprod: the crossproduct matrix of the weighting matrix in the SNP info file.
# nsnps: the total number of SNPs. If not provided, it is set to the first GRM's number of SNPs in vcfile.
# annot.size: a list of the number of SNPs in each output annotation category. If not provided, it is set to the "m" column of vcfile.
# index: a list of column indices in snpinfo or crossprod matching vcfile. The first category should be indexed as 1, the second as 2, and so on.
recalculate_enrichments <- function(vcfile, snpinfo=NA, crossprod=NA, nsnps=NA, annot.size=NA, index=NULL) {
    if(is.na(vcfile)) {
        stop("vcfile must be specifiled.")
    }

    if(is.na(snpinfo) && length(crossprod) == 1 && is.na(crossprod)) {
        stop("Either snpinfo or crossprod must be specified.")
    }

    if(! is.matrix(crossprod)) {
        library(data.table)
        sw = fread(snpinfo, head=T, sep=",")
        annot = colnames(sw)[-1]
        snp = sw$SNP
        sw = sw[,-1]
        sw[is.na(sw)] = 0
        sw = as.matrix(sw)
        crossprod = crossprod(sw)
    }

    mq = read.csv(vcfile)
    mq = mq[-nrow(mq),]
    if(length(index) == 0) {
        index = c(1:ncol(crossprod))
    }
    if(max(index) > ncol(crossprod)) {
        stop("Dim of snpinfo or crossprod NOT consistent with index specified.")
    }
    if(length(index) > nrow(mq)) {
        stop(paste("Length of index must be <= the number of genetic VCs in", vcfile))
    }
    mq = mq[1:length(index),]

    lhs = crossprod[, index]
    # PVE estimate
    var = lhs %*% as.matrix(mq[,9:(8+nrow(mq))]) %*% t(lhs)
    pve.est = lhs %*% mq$enrichment / sum(mq$m)
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

