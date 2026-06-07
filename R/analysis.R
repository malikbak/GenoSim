# =============================================================================
#  GenoSim R package ??? analysis.R
#  HWE, LD, ROH, PCA, FST, diversity metrics
# =============================================================================

#' Hardy-Weinberg Equilibrium Test
#'
#' @description
#' Performs a chi-squared goodness-of-fit test for Hardy-Weinberg equilibrium
#' at each SNP using observed genotype counts versus expected counts under HWE.
#' Loci that are fixed (monomorphic) are flagged as \code{NA} rather than
#' tested.
#'
#' @param geno_matrix Integer matrix of dosage values (0/1/2), dimensions
#'   \eqn{n_{individuals} \times n_{SNPs}}. Rownames = sample IDs,
#'   colnames = SNP IDs.
#' @param alpha Numeric. Significance level for the HWE test. Default
#'   \code{0.05}.
#'
#' @return A \code{data.frame} with one row per SNP and columns:
#' \describe{
#'   \item{\code{snp_id}}{SNP identifier.}
#'   \item{\code{chi2}}{Chi-squared test statistic (1 degree of freedom).}
#'   \item{\code{p_value}}{P-value from \eqn{\chi^2(1)} distribution.}
#'   \item{\code{hwe_pass}}{Logical. \code{TRUE} if \code{p_value >= alpha}.}
#'   \item{\code{sig_label}}{Character: \code{"HWE_pass"},
#'     \code{"HWE_deviation"}, or \code{"fixed"}.}
#' }
#'
#' @examples
#' sim <- simulate_population(n_founders = 100, n_snps = 200,
#'                            n_generations = 1, seed = 1)
#' hwe <- hwe_test(sim$genotypes[[1]])
#' table(hwe$sig_label)
#'
#' @seealso \code{\link{simulate_population}}, \code{\link{compute_ld}}
#' @export
hwe_test <- function(geno_matrix, alpha = 0.05) {
  n   <- nrow(geno_matrix)
  res <- apply(geno_matrix, 2, function(g) {
    n0 <- sum(g == 0L, na.rm = TRUE)
    n1 <- sum(g == 1L, na.rm = TRUE)
    n2 <- sum(g == 2L, na.rm = TRUE)
    p  <- (2*n2 + n1) / (2*n)
    q  <- 1 - p
    e0 <- n * q^2; e1 <- n * 2*p*q; e2 <- n * p^2
    if (any(c(e0, e1, e2) < 1e-10))
      return(c(chi2 = NA_real_, p_value = NA_real_))
    chi2 <- (n0-e0)^2/e0 + (n1-e1)^2/e1 + (n2-e2)^2/e2
    c(chi2 = chi2, p_value = stats::pchisq(chi2, df = 1, lower.tail = FALSE))
  })
  df <- as.data.frame(t(res))
  df$snp_id   <- rownames(df)
  df$hwe_pass <- !is.na(df$p_value) & df$p_value >= alpha
  df$sig_label <- ifelse(is.na(df$p_value), "fixed",
                  ifelse(df$p_value < alpha, "HWE_deviation", "HWE_pass"))
  df[, c("snp_id","chi2","p_value","hwe_pass","sig_label")]
}

#' Compute Pairwise Linkage Disequilibrium
#'
#' @description
#' Computes pairwise \eqn{r^2} and \eqn{D'} (normalised linkage disequilibrium)
#' between all SNP pairs within a specified distance on the same chromosome.
#' Fixed loci are skipped.
#'
#' @param geno_matrix Integer dosage matrix \eqn{n_{individuals} \times n_{SNPs}}.
#' @param snp_map A \code{data.frame} with columns \code{snp_id}, \code{chrom},
#'   and \code{pos_bp}.
#' @param max_snps Integer. If \code{ncol(geno_matrix) > max_snps}, a random
#'   subsample is used. Default \code{500}.
#' @param max_dist_bp Numeric. Maximum inter-SNP distance in base pairs for
#'   which LD is computed. Default \code{1e6} (1 Mb).
#' @param min_pair_n Integer. Minimum number of individuals genotyped at both
#'   SNPs of a pair for that pair to be evaluated. Default \code{2}.
#'
#' @return A \code{data.frame} with columns \code{snp_a}, \code{snp_b},
#'   \code{chrom}, \code{dist_bp}, \code{r2}, and \code{D_prime}.
#'   Returns an empty data frame if no pairs qualify.
#'
#' @details
#' Missing genotype calls (\code{NA}) are handled pairwise: for each SNP pair
#' only individuals genotyped at \emph{both} loci contribute to the estimate.
#' Pairs with fewer than \code{min_pair_n} shared non-missing individuals are
#' skipped. A message reports the overall missingness when any is present.
#'
#' @examples
#' sim <- simulate_population(n_founders = 100, n_snps = 300,
#'                            n_generations = 2, chromosomes = 1:3, seed = 1)
#' ld  <- compute_ld(sim$genotypes[[1]], sim$snp_map, max_snps = 200)
#' hist(ld$r2, main = "r-squared distribution", xlab = "r-sq")
#'
#' @seealso \code{\link{plot_ld_decay}}, \code{\link{hwe_test}}
#' @export
compute_ld <- function(geno_matrix, snp_map, max_snps = 500L,
                       max_dist_bp = 1e6, min_pair_n = 2L) {
  n_na <- sum(is.na(geno_matrix))
  if (n_na > 0)
    message(sprintf(paste0("compute_ld: %d missing genotype(s) (%.2f%%) ",
                           "handled pairwise (complete observations only)."),
                    n_na, 100 * n_na / length(geno_matrix)))
  if (ncol(geno_matrix) > max_snps) {
    idx        <- sample(ncol(geno_matrix), max_snps)
    geno_matrix <- geno_matrix[, idx, drop = FALSE]
    snp_map    <- snp_map[idx, ]
    message(sprintf("LD: subsampled to %d SNPs", max_snps))
  }
  ld_list <- list()
  for (chr in unique(snp_map$chrom)) {
    cidx  <- which(snp_map$chrom == chr)
    if (length(cidx) < 2) next
    g_chr <- geno_matrix[, cidx, drop = FALSE]
    pos   <- snp_map$pos_bp[cidx]
    ids   <- snp_map$snp_id[cidx]
    m     <- length(cidx)
    for (i in seq_len(m - 1)) {
      for (j in (i+1):m) {
        dist <- abs(pos[j] - pos[i])
        if (dist > max_dist_bp) next
        x  <- g_chr[, i]; y <- g_chr[, j]
        ok <- !is.na(x) & !is.na(y)        # complete observations only
        if (sum(ok) < min_pair_n) next
        x  <- x[ok]; y <- y[ok]
        px <- mean(x)/2;  py <- mean(y)/2
        if (px %in% c(0,1) || py %in% c(0,1)) next
        D    <- mean((x/2)*(y/2)) - px*py
        r2   <- D^2 / (px*(1-px)*py*(1-py))
        Dmax <- if (D > 0) min(px*(1-py), py*(1-px)) else max(-px*py, -(1-px)*(1-py))
        Dp   <- if (abs(Dmax) > 1e-12) D/Dmax else NA_real_
        ld_list[[length(ld_list)+1]] <- data.frame(
          snp_a   = ids[i], snp_b = ids[j], chrom = chr,
          dist_bp = dist, r2 = round(r2,5), D_prime = round(Dp,5),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(ld_list) == 0)
    return(data.frame(snp_a=character(),snp_b=character(),chrom=integer(),
                      dist_bp=numeric(),r2=numeric(),D_prime=numeric()))
  do.call(rbind, ld_list)
}

#' Detect Runs of Homozygosity
#'
#' @description
#' Scans each individual's genotype profile chromosome by chromosome to
#' identify contiguous runs of homozygous calls (dosage 0 or 2). A limited
#' number of heterozygous calls within a window are tolerated to account for
#' genotyping error.
#'
#' @param geno_matrix Integer dosage matrix \eqn{n_{individuals} \times n_{SNPs}}.
#' @param snp_map A \code{data.frame} with \code{snp_id}, \code{chrom},
#'   \code{pos_bp}.
#' @param min_snps Integer. Minimum number of SNPs in a qualifying ROH.
#'   Default \code{25}.
#' @param min_length_bp Numeric. Minimum ROH length in base pairs. Default
#'   \code{1e6}.
#' @param het_allowance Integer. Number of heterozygous calls tolerated within
#'   a run before it is broken. Default \code{1}.
#'
#' @details
#' Missing genotype calls (\code{NA}) are treated as non-homozygous, i.e. they
#' count against \code{het_allowance} and therefore break runs in the same way
#' as heterozygous calls. This makes the scan robust to incomplete matrices
#' (e.g. observed-pedigree generations carrying unimputed VCF missingness)
#' instead of erroring. A message reports the missingness when any is present.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{roh_segments}}{A \code{data.frame} of detected ROH segments
#'     with columns \code{ind_id}, \code{chrom}, \code{start_snp},
#'     \code{end_snp}, \code{start_bp}, \code{end_bp}, \code{length_bp},
#'     \code{n_snps_roh}.}
#'   \item{\code{roh_per_individual}}{Aggregated total ROH length and count
#'     per individual.}
#' }
#'
#' @examples
#' sim <- simulate_population(n_founders = 50, n_snps = 500,
#'                            inbreeding_F = 0.25,
#'                            n_generations = 2, chromosomes = 1:5, seed = 1)
#' roh <- detect_roh(sim$genotypes[[3]], sim$snp_map, min_snps = 10)
#' nrow(roh$roh_segments)
#'
#' @seealso \code{\link{plot_roh_per_individual}}
#' @export
detect_roh <- function(geno_matrix, snp_map, min_snps = 25L,
                       min_length_bp = 1e6, het_allowance = 1L) {
  n_na <- sum(is.na(geno_matrix))
  if (n_na > 0)
    message(sprintf(paste0("detect_roh: %d missing genotype(s) (%.2f%%) ",
                           "treated as run-breaking (non-homozygous) calls."),
                    n_na, 100 * n_na / length(geno_matrix)))
  roh_list <- list()
  for (i in seq_len(nrow(geno_matrix))) {
    ind <- rownames(geno_matrix)[i]
    for (chr in unique(snp_map$chrom)) {
      cidx   <- which(snp_map$chrom == chr)
      if (length(cidx) < min_snps) next
      g      <- geno_matrix[i, cidx]
      pos    <- snp_map$pos_bp[cidx]
      ids    <- snp_map$snp_id[cidx]
      # NA (missing) and heterozygous calls are both "not homozygous"; the
      # !is.na() guard prevents `if (!is_hom[k])` from erroring on NA.
      is_hom <- !is.na(g) & (g != 1L)
      in_run <- FALSE; rs <- 1L; hc <- 0L
      for (k in seq_along(g)) {
        if (!is_hom[k]) {
          hc <- hc + 1L
          if (hc > het_allowance) {
            if (in_run && (k-1L-rs+1L) >= min_snps &&
                (pos[k-1L]-pos[rs]) >= min_length_bp)
              roh_list[[length(roh_list)+1]] <- data.frame(
                ind_id=ind, chrom=chr, start_snp=ids[rs], end_snp=ids[k-1L],
                start_bp=pos[rs], end_bp=pos[k-1L],
                length_bp=pos[k-1L]-pos[rs], n_snps_roh=k-1L-rs+1L,
                stringsAsFactors=FALSE)
            in_run <- FALSE; hc <- 0L
          }
        } else {
          if (!in_run) { in_run <- TRUE; rs <- k; hc <- 0L }
        }
      }
      m <- length(g)
      if (in_run && (m-rs+1L) >= min_snps && (pos[m]-pos[rs]) >= min_length_bp)
        roh_list[[length(roh_list)+1]] <- data.frame(
          ind_id=ind, chrom=chr, start_snp=ids[rs], end_snp=ids[m],
          start_bp=pos[rs], end_bp=pos[m], length_bp=pos[m]-pos[rs],
          n_snps_roh=m-rs+1L, stringsAsFactors=FALSE)
    }
  }
  if (length(roh_list) == 0)
    return(list(roh_segments = data.frame(ind_id=character(),chrom=integer(),
                  start_snp=character(),end_snp=character(),start_bp=numeric(),
                  end_bp=numeric(),length_bp=numeric(),n_snps_roh=integer()),
                roh_per_individual = NULL))
  roh_df  <- do.call(rbind, roh_list)
  roh_sum <- stats::aggregate(cbind(total_roh_bp=length_bp, n_roh=length_bp) ~ ind_id,
                               data=roh_df,
                               FUN=function(x) c(sum=sum(x), count=length(x)))
  list(roh_segments=roh_df, roh_per_individual=roh_sum)
}

#' Principal Component Analysis on Genotype Matrix
#'
#' @description
#' Performs genomic PCA on a dosage genotype matrix after removing fixed
#' (zero-variance) loci. Returns PC scores, loadings, and variance explained.
#'
#' @param geno_matrix Integer dosage matrix \eqn{n_{individuals} \times n_{SNPs}}.
#' @param n_pc Integer. Number of principal components to compute. Default
#'   \code{10}. Capped at \eqn{\min(n_{individuals}-1, n_{SNPs})}.
#' @param center Logical. Mean-centre columns before PCA. Default \code{TRUE}.
#' @param scale Logical. Scale columns to unit variance. Default \code{FALSE}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{scores}}{Data frame of PC scores, \eqn{n_{individuals} \times n_{pc}}.}
#'   \item{\code{loadings}}{Data frame of SNP loadings.}
#'   \item{\code{variance_pct}}{Numeric vector of variance explained (\%) per PC.}
#'   \item{\code{cumvar_pct}}{Cumulative variance explained (\%).}
#'   \item{\code{n_snps_used}}{Number of variable SNPs used.}
#' }
#'
#' @examples
#' sim <- simulate_population(n_founders=100, n_snps=500, n_generations=3, seed=1)
#' all_g <- do.call(rbind, sim$genotypes)
#' pca   <- run_pca(all_g, n_pc = 5)
#' plot(pca$scores[,1], pca$scores[,2], xlab="PC1", ylab="PC2")
#'
#' @seealso \code{\link{plot_pca}}
#' @export
run_pca <- function(geno_matrix, n_pc = 10L, center = TRUE, scale = FALSE) {
  n_pc    <- min(n_pc, nrow(geno_matrix) - 1L, ncol(geno_matrix))
  keep    <- which(apply(geno_matrix, 2, stats::var) > 0)
  if (length(keep) < 2) stop("Fewer than 2 variable SNPs for PCA.")
  g_filt  <- geno_matrix[, keep, drop = FALSE]
  pca_res <- stats::prcomp(g_filt, center = center, scale. = scale, rank. = n_pc)
  var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  list(scores       = as.data.frame(pca_res$x[, seq_len(n_pc), drop=FALSE]),
       loadings     = as.data.frame(pca_res$rotation[, seq_len(n_pc), drop=FALSE]),
       variance_pct = round(var_exp[seq_len(n_pc)], 3),
       cumvar_pct   = round(cumsum(var_exp)[seq_len(n_pc)], 3),
       n_snps_used  = length(keep))
}

#' Compute Weir-Cockerham FST Between Generations
#'
#' @description
#' Estimates the fixation index \eqn{F_{ST}} between pairs of generations using
#' a simplified Weir-Cockerham (1984) formula. Reports mean and median across
#' loci.
#'
#' @param geno_list Named list of dosage matrices, as returned in the
#'   \code{genotypes} element of \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#' @param gen_a Integer. First generation index (0-based). Default \code{0}.
#' @param gen_b Integer or \code{NULL}. Second generation index. If \code{NULL}
#'   (default), all consecutive pairs are computed.
#'
#' @return A \code{data.frame} with columns \code{gen_from}, \code{gen_to},
#'   \code{mean_fst}, and \code{median_fst}.
#'
#' @examples
#' sim <- simulate_population(n_founders=100, n_snps=500, n_generations=4, seed=1)
#' fst <- compute_fst(sim$genotypes)
#' print(fst)
#'
#' @references
#' Weir BS, Cockerham CC (1984) Estimating F-statistics for the analysis of
#' population structure. \emph{Evolution} 38:1358???1370.
#'
#' @seealso \code{\link{diversity_metrics}}
#' @export
compute_fst <- function(geno_list, gen_a = 0L, gen_b = NULL) {
  non_empty <- Filter(function(g) !is.null(g) && nrow(g) > 0, geno_list)
  n         <- length(non_empty)
  pairs     <- if (is.null(gen_b)) cbind(0:(n-2), 1:(n-1)) else
                                    matrix(c(gen_a, gen_b), nrow=1)
  do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    g1  <- non_empty[[pairs[i,1]+1]]; g2 <- non_empty[[pairs[i,2]+1]]
    n1  <- nrow(g1); n2 <- nrow(g2)
    p1  <- colMeans(g1)/2; p2 <- colMeans(g2)/2
    pb  <- (n1*p1 + n2*p2)/(n1+n2)
    msp <- n1*(p1-pb)^2 + n2*(p2-pb)^2
    wcp <- (n1*p1*(1-p1) + n2*p2*(1-p2))/(n1+n2-2)
    fst <- ifelse(msp+wcp > 0, msp/(msp+wcp), 0)
    data.frame(gen_from=pairs[i,1], gen_to=pairs[i,2],
               mean_fst=round(mean(fst,na.rm=TRUE),5),
               median_fst=round(stats::median(fst,na.rm=TRUE),5),
               stringsAsFactors=FALSE)
  }))
}

#' Compute Population Diversity Metrics Per Generation
#'
#' @description
#' Computes Nei's gene diversity, proportion of segregating sites, and
#' Watterson's theta estimator for each generation in a simulation result.
#'
#' @param sim_result A list returned by \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#'
#' @return A \code{data.frame} with columns \code{generation},
#'   \code{nei_gene_div}, \code{seg_sites_frac}, and
#'   \code{watterson_theta}.
#'
#' @examples
#' sim <- simulate_population(n_founders=200, n_snps=1000, n_generations=5, seed=1)
#' div <- diversity_metrics(sim)
#' print(div)
#'
#' @references
#' Nei M (1973) Analysis of gene diversity in subdivided populations.
#' \emph{PNAS} 70:3321???3323.
#'
#' @seealso \code{\link{compute_fst}}
#' @export
diversity_metrics <- function(sim_result) {
  do.call(rbind, lapply(seq_along(sim_result$genotypes), function(i) {
    g  <- sim_result$genotypes[[i]]
    p  <- sim_result$allele_freqs[i, ]
    n  <- nrow(g)
    a1 <- sum(1/seq_len(max(n-1,1)))
    S  <- sum(apply(g, 2, stats::var) > 0)
    data.frame(
      generation      = i - 1L,
      nei_gene_div    = round(mean(2*p*(1-p), na.rm=TRUE), 5),
      seg_sites_frac  = round(mean(p > 0 & p < 1, na.rm=TRUE), 5),
      watterson_theta = round(if (a1 > 0) S/a1/ncol(g) else NA_real_, 5),
      stringsAsFactors = FALSE
    )
  }))
}
