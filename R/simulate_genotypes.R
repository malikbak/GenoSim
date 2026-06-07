# =============================================================================
#  GenoSim R package ??? simulate_genotypes.R
#  Forward-time population genetics simulator (population mode)
# =============================================================================

# ---- Internal helpers -------------------------------------------------------

#' Validate simulation parameters
#' @noRd
.validate_params <- function(p) {
  stopifnot(
    "n_founders must be >= 10"            = p$n_founders >= 10,
    "n_snps must be between 1 and 100000" = p$n_snps >= 1 && p$n_snps <= 100000,
    "n_generations must be between 1-10"  = p$n_generations >= 1 && p$n_generations <= 10,
    "inbreeding_F must be in [0, 1)"      = p$inbreeding_F >= 0 && p$inbreeding_F < 1,
    "mut_rate must be in [0, 0.01]"       = p$mut_rate >= 0 && p$mut_rate <= 0.01,
    "selection_s must be in [-1, 1]"      = p$selection_s >= -1 && p$selection_s <= 1,
    "maf_min must be in (0, 0.5]"         = p$maf_min > 0 && p$maf_min <= 0.5,
    "n_offspring_per_couple must be >= 2" = p$n_offspring_per_couple >= 2,
    "chromosomes must be integers 1-22"   = all(p$chromosomes %in% 1:22)
  )
}

#' Assign SNP positions across selected chromosomes
#' @noRd
.assign_snp_positions <- function(n_snps, chromosomes) {
  chrom_lengths <- c(
    `1`=248,`2`=242,`3`=198,`4`=190,`5`=181,`6`=171,`7`=159,`8`=145,
    `9`=138,`10`=133,`11`=135,`12`=133,`13`=114,`14`=107,`15`=102,
    `16`=90,`17`=81,`18`=78,`19`=59,`20`=63,`21`=48,`22`=51
  )
  sel_lengths <- chrom_lengths[as.character(chromosomes)]
  total_len   <- sum(sel_lengths)
  n_per_chrom <- round(n_snps * sel_lengths / total_len)
  diff        <- n_snps - sum(n_per_chrom)
  n_per_chrom[which.max(n_per_chrom)] <-
    n_per_chrom[which.max(n_per_chrom)] + diff

  snp_list <- mapply(function(chr, n, len) {
    if (n == 0) return(NULL)
    data.frame(
      snp_id = paste0("chr", chr, "_", sort(sample.int(len * 1e6L, n))),
      chrom  = chr,
      pos_bp = sort(sample.int(len * 1e6L, n)),
      stringsAsFactors = FALSE
    )
  }, chromosomes, n_per_chrom, sel_lengths[as.character(chromosomes)],
  SIMPLIFY = FALSE)

  do.call(rbind, snp_list[!sapply(snp_list, is.null)])
}

#' Sample genotypes under Hardy-Weinberg with inbreeding
#' @noRd
.sample_genotypes <- function(p_vec, F_coef, n_indiv) {
  q_vec  <- 1 - p_vec
  n_snps <- length(p_vec)
  prob_0 <- q_vec^2 * (1 - F_coef) + q_vec * F_coef
  prob_1 <- 2 * p_vec * q_vec * (1 - F_coef)
  u      <- matrix(stats::runif(n_indiv * n_snps), nrow = n_indiv, ncol = n_snps)
  geno   <- matrix(0L, nrow = n_indiv, ncol = n_snps)
  for (j in seq_len(n_snps))
    geno[, j] <- ifelse(u[, j] < prob_0[j], 0L,
                 ifelse(u[, j] < prob_0[j] + prob_1[j], 1L, 2L))
  geno
}

#' Sample single individual genotype vector (internal)
#' @noRd
.sample_genotypes_vec <- function(p_vec, F_coef) {
  # Unknown frequencies (NA) would yield NA dosages; treat them as fixed
  # ancestral (0) and clamp to [0, 1] so the output is always a valid 0/1/2.
  p_vec[is.na(p_vec)] <- 0
  p_vec  <- pmin(pmax(p_vec, 0), 1)
  q_vec  <- 1 - p_vec
  prob_0 <- q_vec^2 * (1 - F_coef) + q_vec * F_coef
  prob_1 <- 2 * p_vec * q_vec * (1 - F_coef)
  u      <- stats::runif(length(p_vec))
  as.integer(ifelse(u < prob_0, 0L, ifelse(u < prob_0 + prob_1, 1L, 2L)))
}

#' Wright-Fisher drift
#'
#' NA-safe and boundary-safe: missing frequencies are treated as fixed (0), and
#' if \code{rbinom()} returns NA for a probability sitting exactly on (or a
#' floating-point hair past) the [0, 1] boundary, the expected count is used.
#' @noRd
.wf_drift <- function(p_vec, N_eff) {
  size <- 2L * N_eff
  prob <- pmin(pmax(p_vec, 0), 1)
  prob[is.na(prob)] <- 0
  k <- stats::rbinom(length(prob), size = size, prob = prob)
  na_k <- is.na(k)
  if (any(na_k)) k[na_k] <- round(prob[na_k] * size)
  k / size
}

#' Additive selection
#' @noRd
.apply_selection <- function(p_vec, s) {
  q_vec  <- 1 - p_vec
  mean_w <- p_vec^2 * (1+s) + 2*p_vec*q_vec*(1+s/2) + q_vec^2
  p_new  <- (p_vec^2 * (1+s) + p_vec*q_vec*(1+s/2)) / mean_w
  pmin(pmax(p_new, 0), 1)
}

#' Mutation
#' @noRd
.apply_mutation <- function(p_vec, mu)
  p_vec * (1 - mu) + (1 - p_vec) * mu

#' Recombination-aware gamete transmission
#' @noRd
.transmit_gamete <- function(hap1, hap2, snp_map, r_per_bp = 1e-8) {
  gamete <- integer(nrow(snp_map))
  for (chr in unique(snp_map$chrom)) {
    idx <- which(snp_map$chrom == chr)
    if (length(idx) == 1L) {
      gamete[idx] <- if (stats::runif(1) < 0.5) hap1[idx] else hap2[idx]
      next
    }
    r_probs <- 0.5 * (1 - exp(-2 * r_per_bp * diff(snp_map$pos_bp[idx])))
    current <- if (stats::runif(1) < 0.5) 1L else 2L
    for (k in seq_along(idx)) {
      gamete[idx[k]] <- if (current == 1L) hap1[idx[k]] else hap2[idx[k]]
      if (k < length(idx) && stats::runif(1) < r_probs[k])
        current <- 3L - current
    }
  }
  gamete
}

#' Simulate one generation via random mating
#' @noRd
.simulate_generation <- function(parent_genos, snp_map, F_coef, n_offspring) {
  n_parents <- nrow(parent_genos)
  n_snps    <- ncol(parent_genos)
  n_couples <- floor(n_parents / 2)
  pair_idx  <- matrix(sample(n_parents, n_couples * 2L), ncol = 2)

  do.call(rbind, lapply(seq_len(n_couples), function(i) {
    p1_dos <- parent_genos[pair_idx[i, 1], ]
    p2_dos <- parent_genos[pair_idx[i, 2], ]
    # Missing parental dosages (NA, e.g. from unimputed VCF calls carried in an
    # observed-pedigree generation) are treated as homozygous ancestral so that
    # transmitted gametes - and the resulting synthetic genotypes - never become
    # NA. Without the !is.na() guard, NA propagated straight into the offspring.
    h1_p1  <- ifelse(!is.na(p1_dos) & p1_dos == 2L, 1L, 0L)
    h2_p1  <- ifelse(!is.na(p1_dos) & p1_dos >= 1L, 1L, 0L)
    h1_p2  <- ifelse(!is.na(p2_dos) & p2_dos == 2L, 1L, 0L)
    h2_p2  <- ifelse(!is.na(p2_dos) & p2_dos >= 1L, 1L, 0L)
    t(sapply(seq_len(n_offspring), function(k) {
      g1 <- .transmit_gamete(h1_p1, h2_p1, snp_map)
      g2 <- .transmit_gamete(h1_p2, h2_p2, snp_map)
      if (stats::runif(1) < F_coef) g2 <- g1
      g1 + g2
    }))
  }))
}

# ---- Public API -------------------------------------------------------------

#' Simulate Diploid SNP Genotype Data Across Generations (Population Mode)
#'
#' @description
#' Performs forward-time simulation of diploid biallelic SNP genotypes for a
#' panmictic population over up to 10 generations. Biological processes
#' modelled include Hardy-Weinberg genotype sampling, Wright-Fisher genetic
#' drift, additive directional selection, per-locus mutation, and
#' recombination-aware gamete transmission. Inbreeding is implemented via the
#' standard \eqn{F}-coefficient deflation of heterozygosity.
#'
#' @param n_founders Integer. Number of founder individuals. Must be \eqn{\ge 10}.
#' @param n_snps Integer. Number of biallelic SNPs to simulate. Range 1???100,000.
#' @param n_generations Integer. Number of generations to simulate. Range 1???10.
#' @param inbreeding_F Numeric in \eqn{[0, 1)}. Wright's inbreeding coefficient
#'   \eqn{F}. A value of 0 gives Hardy-Weinberg proportions; 0.125 approximates
#'   a first-cousin mating system; 0.25 approximates half-sibling matings.
#' @param mut_rate Numeric in \eqn{[0, 0.01]}. Per-locus per-generation
#'   mutation rate (bidirectional: \eqn{A \to a} and \eqn{a \to A}).
#' @param selection_s Numeric in \eqn{[-1, 1]}. Additive selection coefficient.
#'   Positive values favour the alternative allele; 0 gives neutral evolution.
#' @param maf_min Numeric in \eqn{(0, 0.5]}. Minimum minor allele frequency
#'   assigned to founder SNPs. Founder allele frequencies are drawn uniformly
#'   from \code{[maf_min, 0.5]}.
#' @param n_eff Integer or \code{NULL}. Effective population size \eqn{N_e}
#'   governing genetic drift intensity. Defaults to \code{n_founders}.
#' @param n_offspring_per_couple Integer \eqn{\ge 2}. Number of offspring
#'   produced per mating pair per generation.
#' @param chromosomes Integer vector. Autosome numbers (subset of 1???22) across
#'   which SNPs are distributed proportionally to chromosome length.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#' @param verbose Logical. Whether to print progress messages. Default \code{TRUE}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{genotypes}}{Named list of dosage matrices (0/1/2) per
#'     generation. Each matrix is \eqn{n_{indiv} \times n_{SNPs}}, with
#'     rownames \code{gen0_ind1}, \code{gen1_ind1}, etc.}
#'   \item{\code{snp_map}}{A \code{data.frame} with columns \code{snp_id},
#'     \code{chrom}, \code{pos_bp}, and \code{founder_maf}.}
#'   \item{\code{allele_freqs}}{Matrix of alternative allele frequencies,
#'     dimensions \eqn{(n_{generations}+1) \times n_{SNPs}}.}
#'   \item{\code{summary_stats}}{A \code{data.frame} of per-generation
#'     statistics: observed/expected heterozygosity, \eqn{F_{IS}}, mean MAF,
#'     fraction of fixed loci, and mean dosage.}
#'   \item{\code{params}}{List of all input parameters used.}
#' }
#'
#' @details
#' \strong{Genotype model:}
#' Under inbreeding coefficient \eqn{F}, genotype probabilities are:
#' \deqn{P(AA) = p^2(1-F) + pF}
#' \deqn{P(Aa) = 2pq(1-F)}
#' \deqn{P(aa) = q^2(1-F) + qF}
#'
#' \strong{Drift:} Each generation, allele frequencies are updated via a
#' Binomial draw \eqn{Bin(2N_e, p)} (Wright-Fisher model).
#'
#' \strong{Selection:} Additive model with fitness \eqn{w_{AA}=1+s},
#' \eqn{w_{Aa}=1+s/2}, \eqn{w_{aa}=1}.
#'
#' \strong{Recombination:} Haldane mapping function applied between adjacent
#' SNPs on each chromosome, with default rate \eqn{10^{-8}} per bp per meiosis.
#'
#' @examples
#' # Basic neutral simulation
#' sim <- simulate_population(
#'   n_founders    = 100,
#'   n_snps        = 500,
#'   n_generations = 5,
#'   seed          = 42
#' )
#' # Inspect summary statistics
#' sim$summary_stats
#'
#' # High-inbreeding scenario (first-cousin system)
#' sim_inbred <- simulate_population(
#'   n_founders    = 80,
#'   n_snps        = 1000,
#'   n_generations = 5,
#'   inbreeding_F  = 0.125,
#'   n_eff         = 60,
#'   chromosomes   = 1:5,
#'   seed          = 2024
#' )
#'
#' # Positive selection scenario
#' sim_sel <- simulate_population(
#'   n_founders    = 200,
#'   n_snps        = 2000,
#'   n_generations = 8,
#'   selection_s   = 0.05,
#'   seed          = 999
#' )
#'
#' @seealso
#' \code{\link{simulate_from_pedigree}} for pedigree-constrained simulation,
#' \code{\link{hwe_test}}, \code{\link{compute_ld}}, \code{\link{run_pca}},
#' \code{\link{plot_dashboard}}, \code{\link{export_vcf}}
#'
#' @references
#' Wright S (1931) Evolution in Mendelian populations. \emph{Genetics} 16:97???159.
#'
#' @export
simulate_population <- function(
    n_founders             = 100,
    n_snps                 = 1000,
    n_generations          = 5,
    inbreeding_F           = 0.0,
    mut_rate               = 1e-4,
    selection_s            = 0.0,
    maf_min                = 0.05,
    n_eff                  = NULL,
    n_offspring_per_couple = 3,
    chromosomes            = 1:22,
    seed                   = 42,
    verbose                = TRUE
) {
  params <- as.list(environment())
  .validate_params(params)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(n_eff)) n_eff <- n_founders

  if (verbose) {
    message("=== GenoSim: Population Forward Simulator ===")
    message(sprintf("  Founders: %d | SNPs: %d | Generations: %d",
                    n_founders, n_snps, n_generations))
    message(sprintf("  F=%.3f | mu=%.2e | s=%.3f | Ne=%d",
                    inbreeding_F, mut_rate, selection_s, n_eff))
  }

  snp_map     <- .assign_snp_positions(n_snps, chromosomes)
  n_snps_act  <- nrow(snp_map)
  snp_map$founder_maf <- stats::runif(n_snps_act, maf_min, 0.5)
  flip        <- stats::runif(n_snps_act) < 0.5
  p_founder   <- ifelse(flip, snp_map$founder_maf, 1 - snp_map$founder_maf)

  if (verbose) message("[1/4] Simulating founder genotypes...")
  geno_founders <- .sample_genotypes(p_founder, inbreeding_F, n_founders)
  rownames(geno_founders) <- paste0("gen0_ind", seq_len(n_founders))
  colnames(geno_founders) <- snp_map$snp_id

  af_track <- matrix(NA_real_, nrow = n_generations + 1L, ncol = n_snps_act,
                     dimnames = list(paste0("gen", 0:n_generations), snp_map$snp_id))
  af_track[1L, ] <- p_founder

  geno_store      <- vector("list", n_generations + 1L)
  geno_store[[1]] <- geno_founders
  names(geno_store) <- paste0("gen", 0:n_generations)
  p_current <- p_founder

  if (verbose) message("[2/4] Running forward simulation...")
  for (gen in seq_len(n_generations)) {
    if (verbose) message(sprintf("  -> Generation %d/%d", gen, n_generations))
    new_genos  <- .simulate_generation(geno_store[[gen]], snp_map,
                                        inbreeding_F, n_offspring_per_couple)
    p_obs      <- colMeans(new_genos) / 2
    p_current  <- .apply_mutation(.apply_selection(.wf_drift(p_obs, n_eff),
                                                    selection_s), mut_rate)
    p_current  <- pmin(pmax(p_current, 0), 1)
    af_track[gen + 1L, ] <- p_current
    n_new <- nrow(new_genos)
    rownames(new_genos) <- paste0("gen", gen, "_ind", seq_len(n_new))
    colnames(new_genos) <- snp_map$snp_id
    geno_store[[gen + 1L]] <- new_genos
  }

  if (verbose) message("[3/4] Computing summary statistics...")
  summary_stats <- do.call(rbind, lapply(seq_along(geno_store), function(i) {
    g   <- geno_store[[i]]
    p   <- af_track[i, ]
    obs <- mean(g == 1L)
    exp <- mean(2 * p * (1 - p))
    data.frame(
      generation         = i - 1L,
      n_individuals      = nrow(g),
      n_snps             = ncol(g),
      obs_heterozygosity = round(obs, 5),
      exp_heterozygosity = round(exp, 5),
      inbreeding_fis     = round(if (exp > 0) 1 - obs/exp else NA_real_, 5),
      mean_maf           = round(mean(pmin(p, 1-p)), 5),
      frac_fixed         = round(mean(p <= 0 | p >= 1), 5),
      mean_dosage        = round(mean(g), 5),
      stringsAsFactors   = FALSE
    )
  }))

  if (verbose) {
    message("\n=== Simulation Complete ===")
    message(sprintf("  Total individuals: %d | SNPs: %d",
                    sum(sapply(geno_store, nrow)), n_snps_act))
  }

  list(genotypes = geno_store, snp_map = snp_map,
       allele_freqs = af_track, summary_stats = summary_stats, params = params)
}
