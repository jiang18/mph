# Oct 27, 2023: initial release along with MPH v0.49.2
# Jan 27, 2024: revised recalculate_enrichments() along with MPH v0.52.0
# Oct 09, 2024: added the covariance matrix of PVE estimates in the output of recalculate_enrichments()
# Dec 10, 2025: added calculate_h2() function to calculate heritability and standard error from partitioned variance components using delta method

# Recalculate the estimates and SEs of PVEs and enrichments from MPH VC estimates
# vcfile: the .mq.vc.csv file produced by `mph --reml` or `mph --minque`.
# crossprod: the crossproduct of SNP incidence matrix and SNP weighting matrix. Its row names and columns should match the annotation categories of interest and the rows of vcfile, respectively.
# nsnps: the total number of SNPs. If `NA`, it is set to the first GRM's number of SNPs in vcfile.
# annot.size: a list of the number of SNPs in each annotation category listed in the row names of crossprod. It can be calculated from the column-wise sum of the corresponding SNP incidence matrix. If `NA`, it is set to the `m` column of vcfile.
# trait.x and trait.y: trait names matching the first two columns of vcfile. If both are `NA`, they are set to the trait names in the first line of vcfile. Otherwise, if only one is `NA`, it is set to the value of the other.
recalculate_enrichments <- function(vcfile, crossprod, nsnps=NA, annot.size=NA, trait.x=NA, trait.y=NA) {
    if(is.na(vcfile)) {
        stop("vcfile must be specifiled.")
    }
    if(length(crossprod) == 1 && is.na(crossprod)) {
        stop("crossprod must be specified.")
    }

    mq = read.csv(vcfile)

    if(is.na(trait.x) && is.na(trait.y)) {
      trait.x = mq$trait_x[1]
      trait.y = mq$trait_y[1]
    } else if (is.na(trait.x)) {
       trait.x = trait.y
    } else if (is.na(trait.y)) {
       trait.y = trait.x
    }

    mq = mq[mq$trait_x == trait.x & mq$trait_y == trait.y, -c(1,2)]
    mq = mq[-nrow(mq),]
    if(ncol(crossprod) > nrow(mq)) {
        stop(paste("ncol(crossprod) should equal the number of genomic VCs in", vcfile))
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

    out = cbind(annot.size/nsnps, pve.est, pve.se, e.est, e.se, var/(sum(mq$m)^2))
    rownames(out) = rownames(crossprod)
    colnames(out) = c("prop", "pve", "pve.se", "enrichment", "enrichment.se", rownames(crossprod))
    return(out)
}

# Write GRM into MPH format files (.iid and .bin)
write_grm=function(prefix, iid, grm){
  idfile = paste0(prefix, ".grm.iid")
  binfile = paste0(prefix, ".grm.bin")
  if(length(iid) != ncol(grm)) {
    print("Error: the size of iid is not equal to ncol of grm.")
    return(1)
  }
  write.table(iid, file=idfile, quote=F, row.names=F, col.names=F)
  print(paste0("Completed writing the iid file: ", prefix, ".grm.iid"))
  
  # sample size
  np = length(iid)
  # number of markers set to 1 for pseudo GRMs
  nm = 1
  con = file(binfile, "wb")
  writeBin(as.integer(np), con, size=4)
  writeBin(as.numeric( c(nm,grm[lower.tri(grm, diag=TRUE)]) ), con, size=4)
  close(con)
  print(paste0("Completed writing the bin file: ", prefix, ".grm.bin"))
}

# Read MPH GRM files (.iid and .bin)
read_grm=function(prefix){
  idfile = paste0(prefix, ".grm.iid")
  binfile = paste0(prefix, ".grm.bin")
  iid = read.table(idfile)
  iid = iid[,1]
  mph_file = file(binfile, "rb")
  np = readBin(mph_file, n=1, what=integer(), size=4)
  if(length(iid) != np) {
    print("Error: the size of iid is not equal to ncol of grm.")
    return(1)
  }
  nm = readBin(mph_file, n=1, what=numeric(), size=4)
  grm = matrix(0, nrow=np, ncol=np)
  grm[lower.tri(grm, diag=TRUE)] = readBin(mph_file, n=np*(np+1)/2, what=numeric(), size=4)
  close(mph_file)
  
  grm = grm / nm
  colnames(grm) = iid
  rownames(grm) = iid

  return(grm)
}

# Calculate heritability from partitioned variance components
#
# Inputs:
#   vc_estimates: Numeric vector [σ²_1, ..., σ²_N]
#   cov_matrix: Covariance matrix for vc_estimates (N × N)
#   n_genetic: Number of genetic VCs (default: N-1, assumes last component is error)
#
# Returns: List with
#   - Vg: Sum of genetic VCs
#   - Vg_se: SE of Vg
#   - Vp: Total variance (Vg + non-genetic components)
#   - Vp_se: SE of Vp
#   - h2: Heritability (Vg / Vp)
#   - h2_se: SE of heritability
calculate_h2 <- function(vc_estimates, cov_matrix, n_genetic = NULL) {
  n_total <- length(vc_estimates)
  
  # Default: first n-1 components are genetic
  if (is.null(n_genetic)) {
    n_genetic <- n_total - 1
  }
  
  cov_matrix = as.matrix(cov_matrix)
  
  # Step 1: Calculate sums
  Vg <- sum(vc_estimates[1:n_genetic])
  Vp <- sum(vc_estimates)
  
  # Gradient vectors for sums
  grad_Vg <- c(rep(1, n_genetic), rep(0, n_total - n_genetic))
  grad_Vp <- rep(1, n_total)
  
  # Variances of sums
  var_Vg <- t(grad_Vg) %*% cov_matrix %*% grad_Vg
  var_Vp <- t(grad_Vp) %*% cov_matrix %*% grad_Vp
  
  # Covariance between Vg and Vp
  cov_Vg_Vp <- t(grad_Vg) %*% cov_matrix %*% grad_Vp
  
  # Step 2: Calculate heritability and its SE using delta method
  h2 <- Vg / Vp
  
  # Gradient for ratio: ∂(Vg/Vp)/∂Vg = 1/Vp, ∂(Vg/Vp)/∂Vp = -Vg/Vp²
  grad_h2 <- c(1/Vp, -Vg/Vp^2)
  cov_sums <- matrix(c(var_Vg, cov_Vg_Vp, 
                       cov_Vg_Vp, var_Vp), nrow = 2)
  var_h2 <- t(grad_h2) %*% cov_sums %*% grad_h2
  
  return(list(
    Vg = Vg,
    Vg_se = sqrt(var_Vg),
    Vp = Vp,
    Vp_se = sqrt(var_Vp),
    h2 = h2,
    h2_se = sqrt(var_h2)
  ))
}

# Calculate ratio estimates and their standard errors
#
# Inputs:
#   estimates: Numeric vector [y, x1, x2, ..., xn]
#   cov_matrix: Square covariance matrix for estimates
#
# Returns: List with
#   - estimates: Ratio estimates [x1/y, x2/y, ..., xn/y]
#   - se: Standard errors of ratio estimates
#   - cov: Covariance matrix of ratio estimates
#
# Note: Uses multivariate Delta method. Assumes y not near zero.
ratio_estimates <- function(estimates, cov_matrix) {
  y <- estimates[1]
  x <- estimates[-1]
  n <- length(x)
  
  # Point estimates
  ratios <- x / y
  
  # Jacobian matrix
  J <- matrix(0, nrow = n, ncol = n + 1)
  J[, 1] <- -x / y^2
  diag(J[, -1]) <- 1 / y
  
  # Covariance matrix of ratios
  cov_ratios <- J %*% cov_matrix %*% t(J)
  
  # Standard errors
  se_ratios <- sqrt(diag(cov_ratios))
  
  return(list(estimates = ratios, se = se_ratios, cov = cov_ratios))
}
