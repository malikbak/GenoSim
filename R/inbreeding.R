# =============================================================================
#  GenoSim R package - inbreeding.R
#  Founder-referenced (temporal) inbreeding statistics.
#
#  Wright's within-generation F_IS is referenced to each generation's *own*
#  allele frequencies, so a single round of random mating restores
#  Hardy-Weinberg proportions and F_IS returns to ~0 - it cannot represent the
#  cumulative ("evolutionary") inbreeding that builds up over generations.
#  These helpers add founder-referenced statistics (baseline = the first
#  non-empty / founder generation) that accumulate as diversity is lost:
#
#    H_T          founder gene diversity (unbiased), the fixed baseline
#    F_IS(t)      1 - H_o(t)/H_e(t)      within-generation HWE departure (unbiased)
#    F_ST(t)      1 - H_e(t)/H_T         cumulative drift / diversity loss
#    F_IT(t)      1 - H_o(t)/H_T         total individual inbreeding vs founders
#    1 - F_IT = (1 - F_IS)(1 - F_ST)     (exact, over the same locus set)
#
#  Expected drift trajectory:  E[F_ST(t)] = 1 - (1 - 1/(2 Ne))^t
#  Realised Ne from decay:     Ne_hat   = 1 / (2 (1 - H_e(t)/H_e(t-1)))
#
#  All quantities use the unbiased Nei gene-diversity estimator
#  He = (2n / (2n - 1)) * 2 p (1 - p), are computed over loci polymorphic in the
#  founders, and use ratio-of-sums (Nei's G_ST form) so low-frequency loci do
#  not destabilise the estimates. Every step is NA-aware.
# =============================================================================

#' Per-locus observed het, unbiased gene diversity and called sample size
#' @param g integer dosage matrix (individuals x SNPs), possibly with NA
#' @return list with He (unbiased gene diversity), Ho (observed het) and n
#'   (called individuals), each a per-locus numeric vector
#' @noRd
.locus_diversity <- function(g) {
  if (is.null(g)) return(list(He = numeric(0), Ho = numeric(0), n = integer(0)))
  L <- ncol(g)
  if (is.null(L) || L == 0L || nrow(g) == 0L)
    return(list(He = rep(NA_real_, L), Ho = rep(NA_real_, L), n = rep(0L, L)))
  n  <- colSums(!is.na(g))
  s  <- colSums(g, na.rm = TRUE)
  p  <- ifelse(n > 0, s / (2 * n), NA_real_)
  nz <- ifelse(n > 0, n, NA_real_)
  Ho <- colSums(g == 1L, na.rm = TRUE) / nz
  # Unbiased Nei gene diversity (needs >= 2 diploid individuals for the
  # (2n)/(2n-1) correction; fixed loci correctly give 0).
  He <- ifelse(n >= 2, (2 * n / (2 * n - 1)) * 2 * p * (1 - p), NA_real_)
  list(He = as.numeric(He), Ho = as.numeric(Ho), n = as.integer(n))
}

#' Founder-referenced inbreeding trajectory for a list of genotype matrices
#'
#' @param geno_list list of per-generation dosage matrices (gen0 first)
#' @param n_eff effective population size used for the expected drift curve
#'   (or NULL to omit it)
#' @param ped_F_by_gen optional numeric vector (length = length(geno_list)) of
#'   the mean pedigree inbreeding coefficient per generation; NA where unknown
#' @return data.frame with one row per generation and columns
#'   \code{fis_unbiased}, \code{fst_vs_founder}, \code{fit_vs_founder},
#'   \code{expected_fst_drift}, \code{ne_estimate}, \code{mean_pedigree_F}
#' @noRd
.inbreeding_trajectory <- function(geno_list, n_eff = NULL, ped_F_by_gen = NULL) {
  G <- length(geno_list)
  empty_row <- function() data.frame(
    fis_unbiased = NA_real_, fst_vs_founder = NA_real_, fit_vs_founder = NA_real_,
    expected_fst_drift = NA_real_, ne_estimate = NA_real_, mean_pedigree_F = NA_real_,
    stringsAsFactors = FALSE)
  if (G == 0L) return(empty_row()[0, ])

  # Founder baseline = first generation that actually contains individuals.
  has_indiv <- vapply(geno_list, function(g) !is.null(g) && nrow(g) > 0L, logical(1))
  if (!any(has_indiv))
    return(do.call(rbind, replicate(G, empty_row(), simplify = FALSE)))
  founder_index <- which(has_indiv)[1]

  fl   <- .locus_diversity(geno_list[[founder_index]])
  poly <- which(!is.na(fl$He) & fl$He > 0)        # loci segregating in founders
  HT   <- if (length(poly)) sum(fl$He[poly]) else NA_real_

  prevHS <- NA_real_
  rows <- vector("list", G)
  for (i in seq_len(G)) {
    g <- geno_list[[i]]
    t <- i - founder_index                        # generations since founder
    ok <- has_indiv[i] && length(poly) > 0 && is.finite(HT) && HT > 0

    SumHS <- SumHo <- fis <- fst <- fit <- ne_est <- NA_real_
    if (ok) {
      dl    <- .locus_diversity(g)
      SumHS <- sum(dl$He[poly], na.rm = TRUE)
      SumHo <- sum(dl$Ho[poly], na.rm = TRUE)
      if (is.finite(SumHS) && SumHS > 0) fis <- 1 - SumHo / SumHS
      fst <- 1 - SumHS / HT
      fit <- 1 - SumHo / HT
      if (is.finite(prevHS) && prevHS > 0 && is.finite(SumHS)) {
        R <- SumHS / prevHS
        if (R < 1) ne_est <- 1 / (2 * (1 - R))    # realised Ne from He decay
      }
      prevHS <- SumHS
    }

    exp_fst <- if (!is.null(n_eff) && length(n_eff) == 1L && is.finite(n_eff) &&
                   n_eff > 0 && t >= 0)
      1 - (1 - 1 / (2 * n_eff))^t else NA_real_

    mpf <- if (!is.null(ped_F_by_gen) && i <= length(ped_F_by_gen))
      ped_F_by_gen[i] else NA_real_

    rows[[i]] <- data.frame(
      fis_unbiased       = round(fis, 5),
      fst_vs_founder     = round(fst, 5),
      fit_vs_founder     = round(fit, 5),
      expected_fst_drift = round(exp_fst, 5),
      ne_estimate        = round(ne_est, 3),
      mean_pedigree_F    = round(mpf, 5),
      stringsAsFactors   = FALSE)
  }
  do.call(rbind, rows)
}
