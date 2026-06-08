# =============================================================================
#  GenoSim R package - pedigree_sim.R
#  Pedigree-constrained forward simulator with pedigree‑informed phasing
# =============================================================================

.MISSING_PARENT <- c("0", NA, "", "NA")

# -----------------------------------------------------------------------------
# Helper: phase a child given known parental haplotypes and genotypes
# -----------------------------------------------------------------------------

#' Phase child genotypes using pedigree information
#'
#' For each SNP, if the child is heterozygous and one parent is homozygous,
#' the phase becomes deterministic. Otherwise fall back to random assignment.
#'
#' @param child_dos numeric vector of dosages (0,1,2) for the child
#' @param father_dos,mother_dos dosage vectors of parents
#' @param father_hap1,father_hap2 phased haplotypes of father (if known)
#' @param mother_hap1,mother_hap2 phased haplotypes of mother (if known)
#' @return list with hap1 and hap2 integer vectors
#' @noRd
.phase_pedigree_informed <- function(child_dos, father_dos, mother_dos,
                                     father_hap1, father_hap2,
                                     mother_hap1, mother_hap2) {
  n <- length(child_dos)
  hap1 <- integer(n)
  hap2 <- integer(n)
  for (i in seq_len(n)) {
    d <- child_dos[i]
    if (is.na(d)) {
      # Missing child genotype: default to homozygous ancestral (matches
      # .phase_genotype()'s NA handling) so phasing never sees NA.
      hap1[i] <- 0L; hap2[i] <- 0L
    } else if (d == 0L) {
      hap1[i] <- 0L; hap2[i] <- 0L
    } else if (d == 2L) {
      hap1[i] <- 1L; hap2[i] <- 1L
    } else { # d == 1L, heterozygous
      fd <- father_dos[i]; md <- mother_dos[i]
      # Case 1: one parent homozygous, the other not (or missing)
      if (!is.na(fd) && fd %in% c(0L,2L) && (is.na(md) || md != 1L)) {
        # Father is homozygous
        paternal_allele <- if (fd == 2L) 1L else 0L
        maternal_allele <- 1L - paternal_allele   # because child is heterozygous
        hap1[i] <- paternal_allele
        hap2[i] <- maternal_allele
      } else if (!is.na(md) && md %in% c(0L,2L) && (is.na(fd) || fd != 1L)) {
        # Mother is homozygous
        maternal_allele <- if (md == 2L) 1L else 0L
        paternal_allele <- 1L - maternal_allele
        hap1[i] <- paternal_allele
        hap2[i] <- maternal_allele
      } else {
        # Ambiguous: random phase (coin flip)
        if (stats::runif(1) < 0.5) {
          hap1[i] <- 0L; hap2[i] <- 1L
        } else {
          hap1[i] <- 1L; hap2[i] <- 0L
        }
      }
    }
  }
  list(hap1 = hap1, hap2 = hap2)
}

# -----------------------------------------------------------------------------
# Original helpers (unchanged except minor adaptions)
# -----------------------------------------------------------------------------

.phase_genotype <- function(dosage_vec) {
  hap1 <- integer(length(dosage_vec))
  hap2 <- integer(length(dosage_vec))
  for (j in seq_along(dosage_vec)) {
    d <- if (is.na(dosage_vec[j])) 0L else as.integer(dosage_vec[j])
    if      (d == 0L) { hap1[j] <- 0L; hap2[j] <- 0L }
    else if (d == 2L) { hap1[j] <- 1L; hap2[j] <- 1L }
    else { if (stats::runif(1) < 0.5) { hap1[j] <- 0L; hap2[j] <- 1L }
      else                        { hap1[j] <- 1L; hap2[j] <- 0L } }
  }
  list(hap1 = hap1, hap2 = hap2)
}

.transmit_gamete_ped <- function(hap1, hap2, snp_map, r_per_bp = 1e-8) {
  gamete <- integer(nrow(snp_map))
  for (chr in unique(snp_map$chrom)) {
    idx <- which(snp_map$chrom == chr)
    if (length(idx) == 1L) {
      gamete[idx] <- if (stats::runif(1) < 0.5) hap1[idx] else hap2[idx]
      next
    }
    r_prob  <- 0.5 * (1 - exp(-2 * r_per_bp * diff(snp_map$pos_bp[idx])))
    current <- if (stats::runif(1) < 0.5) 1L else 2L
    for (k in seq_along(idx)) {
      gamete[idx[k]] <- if (current == 1L) hap1[idx[k]] else hap2[idx[k]]
      if (k < length(idx) && stats::runif(1) < r_prob[k])
        current <- 3L - current
    }
  }
  gamete
}

# NOTE: The forward-time helpers .sample_genotypes_vec(), .wf_drift(),
# .apply_selection(), .apply_mutation() and .simulate_generation() are defined
# once in simulate_genotypes.R and shared here. They were previously duplicated
# in this file, but because packages load source files alphabetically the
# simulate_genotypes.R copies always shadowed these, so the duplicates were dead
# code. They are now removed to keep a single (NA-safe, boundary-safe) source of
# truth and avoid the load-order trap.

.draw_population_afs <- function(p_anc, Fst) {
  if (Fst == 0) return(p_anc)
  if (Fst == 1) {
    return(as.numeric(stats::runif(length(p_anc)) < p_anc))
  }
  alpha <- p_anc * (1 - Fst) / Fst
  beta  <- (1 - p_anc) * (1 - Fst) / Fst
  stats::rbeta(length(p_anc), alpha, beta)
}

.load_af_from_file <- function(file, sep = NULL) {
  if (!file.exists(file)) stop("File not found: ", file)
  ext <- tolower(tools::file_ext(file))
  if (ext == "rds") {
    dat <- readRDS(file)
  } else {
    if (is.null(sep)) {
      sep <- if (ext == "csv") "," else if (ext %in% c("tsv","tab")) "\t" else ","
    }
    dat <- utils::read.table(file, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  }
  if (!all(c("snp_id", "af") %in% colnames(dat))) {
    stop("File must contain columns 'snp_id' and 'af'")
  }
  af_vec <- dat$af
  names(af_vec) <- dat$snp_id
  af_vec
}

.reference_populations <- function() {
  fst <- c(AFR = 0.152, AMR = 0.043, EAS = 0.086, SAS = 0.075,
           EUR = 0.042, FIN = 0.051, MID = 0.058, OCE = 0.093)
  list(Fst = fst, populations = names(fst))
}

# -----------------------------------------------------------------------------
# Main simulation function with pedigree‑informed phasing
# -----------------------------------------------------------------------------

#' Simulate Genotype Data Constrained by a Real Family Pedigree and VCF
#'
#' @param vcf_cohort List containing \code{geno_matrix} and \code{snp_map}
#' @param pedigree data.frame with columns individual_id, father_id, mother_id, generation, inbreeding_F
#' @param extra_generations Integer \eqn{[0,10]}
#' @param mut_rate Numeric mutation rate per generation
#' @param selection_s Numeric selection coefficient
#' @param n_offspring_per_couple Integer
#' @param n_eff Effective population size (or NULL)
#' @param match_by Column name for matching VCF IDs
#' @param reference_population NULL, vector, data.frame, or file path
#' @param population Character (e.g., "SAS") for Balding‑Nichols
#' @param ancestral_af Ancestral frequencies for BN model
#' @param seed Random seed
#' @param verbose Logical
#'
#' @return List with genotypes, pedigree, snp_map, allele_freqs, summary_stats, haplotypes, log, params
#' @export
simulate_from_pedigree <- function(
    vcf_cohort,
    pedigree,
    extra_generations      = 3L,
    mut_rate               = 1e-4,
    selection_s            = 0.0,
    n_offspring_per_couple = 3L,
    n_eff                  = NULL,
    match_by               = "individual_id",
    reference_population   = NULL,
    population             = NULL,
    ancestral_af           = NULL,
    seed                   = 42L,
    verbose                = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(extra_generations >= 0L, extra_generations <= 10L)
  
  geno_vcf <- vcf_cohort$geno_matrix
  snp_map  <- vcf_cohort$snp_map
  ped      <- pedigree
  n_snps   <- ncol(geno_vcf)
  sim_log  <- character()
  
  if (verbose) {
    message("=== GenoSim: Pedigree-Constrained Simulation (Pedigree-Informed Phasing) ===")
    message(sprintf("  VCF individuals : %d | SNPs: %d", nrow(geno_vcf), n_snps))
    message(sprintf("  Pedigree members: %d | Extra gens: %d", nrow(ped), extra_generations))
  }
  
  # ---- Reference frequencies ----
  ref_af <- NULL
  if (!is.null(reference_population)) {
    if (is.numeric(reference_population)) {
      if (length(reference_population) != n_snps)
        stop("Length mismatch")
      ref_af <- reference_population
      names(ref_af) <- snp_map$snp_id
    } else if (is.data.frame(reference_population)) {
      if (!all(c("snp_id", "af") %in% colnames(reference_population)))
        stop("Data frame must have snp_id and af")
      ref_af <- rep(NA_real_, n_snps)
      names(ref_af) <- snp_map$snp_id
      ref_af[reference_population$snp_id] <- reference_population$af
      ref_af[is.na(ref_af)] <- 0.5
    } else if (is.character(reference_population) && file.exists(reference_population)) {
      ref_df <- utils::read.table(reference_population, header = TRUE, sep = ",",
                                  stringsAsFactors = FALSE)
      if (!all(c("snp_id", "af") %in% colnames(ref_df)))
        stop("File must have snp_id and af")
      ref_af <- rep(NA_real_, n_snps)
      names(ref_af) <- snp_map$snp_id
      ref_af[ref_df$snp_id] <- ref_df$af
      ref_af[is.na(ref_af)] <- 0.5
    } else stop("Invalid reference_population")
    if (verbose) message("  Using external reference frequencies.")
  }
  
  # ---- Population stratification ----
  pop_af <- NULL
  if (!is.null(population)) {
    if (verbose) message("  Using population-stratified frequencies for ", population)
    if (is.character(ancestral_af))
      ancestral_af <- .load_af_from_file(ancestral_af)
    if (is.null(ancestral_af)) stop("ancestral_af required when population is given")
    if (length(ancestral_af) != n_snps) stop("Length mismatch")
    if (is.null(names(ancestral_af))) names(ancestral_af) <- snp_map$snp_id
    pop_params <- .reference_populations()
    if (!population %in% pop_params$populations)
      stop("Population '", population, "' not recognised")
    Fst <- pop_params$Fst[population]
    pop_af <- .draw_population_afs(ancestral_af, Fst)
    if (verbose) message("    Drawn population AFs (mean = ", round(mean(pop_af), 3), ")")
  }
  
  # ---- Match VCF and pedigree ----
  vcf_ids <- rownames(geno_vcf)
  ped_ids <- ped[[match_by]]
  matched <- intersect(vcf_ids, ped_ids)
  if (length(matched) == 0)
    stop("No VCF IDs match pedigree.")
  if (verbose) message(sprintf("  Matched: %d individuals", length(matched)))
  
  sim_log <- c(sim_log, sprintf("Matched %d individuals", length(matched)))
  
  # ---- Phase observed genotypes ----
  # We store both the genotype and the phased haplotypes for each individual
  ind_geno <- stats::setNames(lapply(matched, function(s) geno_vcf[s, ]), matched)
  # Initial phasing: random (will later be refined when parents become available)
  hap_store <- stats::setNames(lapply(matched, function(s) .phase_genotype(geno_vcf[s, ])), matched)
  
  max_obs_gen  <- max(ped$generation)
  n_total_gens <- max_obs_gen + 1L + extra_generations
  geno_store   <- vector("list", n_total_gens)
  af_store     <- matrix(NA_real_, n_total_gens, n_snps,
                         dimnames = list(paste0("gen", 0:(n_total_gens-1L)),
                                         snp_map$snp_id))
  
  # Generation 0
  g0_ids <- intersect(ped$individual_id[ped$generation == 0L], matched)
  if (length(g0_ids) > 0) {
    g0mat <- geno_vcf[g0_ids, , drop = FALSE]
    geno_store[[1L]] <- g0mat
    af_store[1L, ]   <- colMeans(g0mat, na.rm = TRUE) / 2
  }
  
  # ---- Walk generations ----
  for (gen in sort(unique(ped$generation))[-1L]) {
    if (verbose) message(sprintf("  -> Pedigree generation %d", gen))
    gen_ids  <- ped$individual_id[ped$generation == gen]
    gen_mat  <- matrix(NA_integer_, length(gen_ids), n_snps,
                       dimnames = list(gen_ids, snp_map$snp_id))
    
    for (ind_id in gen_ids) {
      row <- ped[ped$individual_id == ind_id, ]
      fid <- row$father_id; mid <- row$mother_id
      F_i <- row$inbreeding_F
      
      if (ind_id %in% matched) {
        # Individual is genotyped: use observed data, but improve its phase if parents are known
        if (fid %in% names(ind_geno) && mid %in% names(ind_geno)) {
          # Parents available → re‑phase using pedigree information
          father_hap <- hap_store[[fid]]
          mother_hap <- hap_store[[mid]]
          new_phase <- .phase_pedigree_informed(
            child_dos = geno_vcf[ind_id, ],
            father_dos = ind_geno[[fid]], mother_dos = ind_geno[[mid]],
            father_hap1 = father_hap$hap1, father_hap2 = father_hap$hap2,
            mother_hap1 = mother_hap$hap1, mother_hap2 = mother_hap$hap2
          )
          hap_store[[ind_id]] <- new_phase
        } else {
          # Keep existing random phase
          if (!ind_id %in% names(hap_store))
            hap_store[[ind_id]] <- .phase_genotype(geno_vcf[ind_id, ])
        }
        gen_mat[ind_id, ] <- geno_vcf[ind_id, ]
        sim_log <- c(sim_log, sprintf("Gen%d %s: observed VCF (phase %s)",
                                      gen, ind_id,
                                      if (fid %in% names(ind_geno) && mid %in% names(ind_geno)) "pedigree-informed" else "random"))
        next
      }
      
      # Individual not genotyped: simulate from parents or impute
      f_ok <- !fid %in% .MISSING_PARENT && fid %in% names(ind_geno)
      m_ok <- !mid %in% .MISSING_PARENT && mid %in% names(ind_geno)
      
      if (f_ok && m_ok) {
        # Both parents known and already have genotypes
        if (!fid %in% names(hap_store))
          hap_store[[fid]] <- .phase_genotype(ind_geno[[fid]])
        if (!mid %in% names(hap_store))
          hap_store[[mid]] <- .phase_genotype(ind_geno[[mid]])
        g1 <- .transmit_gamete_ped(hap_store[[fid]]$hap1, hap_store[[fid]]$hap2, snp_map)
        g2 <- .transmit_gamete_ped(hap_store[[mid]]$hap1, hap_store[[mid]]$hap2, snp_map)
        if (!is.na(F_i) && F_i > 0 && stats::runif(1) < F_i) g2 <- g1
        if (mut_rate > 0) {
          mm <- stats::runif(n_snps) < mut_rate
          g1[mm] <- 1L - g1[mm]; g2[mm] <- 1L - g2[mm]
        }
        dos <- g1 + g2
        gen_mat[ind_id, ] <- dos
        ind_geno[[ind_id]]  <- dos
        # Phase of offspring is deterministic from transmission
        hap_store[[ind_id]] <- list(hap1 = g1, hap2 = g2)
        sim_log <- c(sim_log, sprintf("Gen%d %s: transmitted from %s x %s (F=%.3f)",
                                      gen, ind_id, fid, mid, F_i))
      } else {
        # Impute from population (no parents)
        if (!is.null(pop_af)) {
          p_pop <- pop_af
          src <- "population-stratified"
        } else if (!is.null(ref_af)) {
          p_pop <- ref_af
          src <- "reference"
        } else {
          p_pop <- if (!is.null(geno_store[[gen]]) && nrow(geno_store[[gen]]) > 0)
            colMeans(geno_store[[gen]], na.rm = TRUE) / 2
          else
            colMeans(geno_vcf[matched, , drop = FALSE], na.rm = TRUE) / 2
          src <- "internal"
        }
        F_use <- if (!is.na(F_i) && F_i > 0) F_i else 0.0
        dos <- .sample_genotypes_vec(p_pop, F_use)
        gen_mat[ind_id, ] <- dos
        ind_geno[[ind_id]]  <- dos
        # Random phase for imputed individuals (no parents to inform)
        hap_store[[ind_id]] <- .phase_genotype(dos)
        sim_log <- c(sim_log, sprintf("Gen%d %s: %s imputed (parents unknown)", gen, ind_id, src))
      }
    }
    
    geno_store[[gen + 1L]] <- gen_mat
    af_store[gen + 1L, ]   <- colMeans(gen_mat, na.rm = TRUE) / 2
  }
  
  # ---- Extra synthetic generations (unchanged) ----
  if (extra_generations > 0L) {
    if (verbose) message(sprintf("\n[Sim] Appending %d synthetic generations...", extra_generations))
    last <- geno_store[[max_obs_gen + 1L]]
    if (is.null(last) || nrow(last) == 0)
      last <- geno_vcf[matched, , drop = FALSE]
    if (is.null(n_eff)) n_eff <- nrow(last)
    mean_F <- mean(ped$inbreeding_F, na.rm = TRUE)
    
    for (eg in seq_len(extra_generations)) {
      si  <- max_obs_gen + 1L + eg
      if (verbose) message(sprintf("  -> Synthetic generation %d", eg))
      ng  <- .simulate_generation(last, snp_map, mean_F, n_offspring_per_couple)
      # Any residual missing dosage would otherwise yield NA allele frequencies
      # (and break detect_roh()/compute_ld() downstream); fall back to 0.5.
      af_ng <- colMeans(ng, na.rm = TRUE) / 2
      af_ng[is.na(af_ng)] <- 0.5
      p   <- .apply_mutation(.apply_selection(.wf_drift(af_ng, n_eff), selection_s), mut_rate)
      p   <- pmin(pmax(p, 0), 1)
      rownames(ng) <- paste0("synth_gen", max_obs_gen + eg, "_ind", seq_len(nrow(ng)))
      colnames(ng) <- snp_map$snp_id
      geno_store[[si]] <- ng
      af_store[si, ]   <- p
      last <- ng
      sim_log <- c(sim_log, sprintf("Synthetic gen %d: %d individuals", eg, nrow(ng)))
    }
  }
  
  # ---- Post-process ----
  geno_store <- geno_store[seq_len(n_total_gens)]
  names(geno_store) <- paste0("gen", 0:(n_total_gens - 1L))
  for (i in which(sapply(geno_store, is.null)))
    geno_store[[i]] <- matrix(integer(0), 0, n_snps,
                              dimnames = list(NULL, snp_map$snp_id))
  
  snp_map$founder_maf <- snp_map$cohort_maf
  
  summary_stats <- do.call(rbind, lapply(seq_along(geno_store), function(i) {
    g <- geno_store[[i]]
    has_indiv <- !is.null(g) && nrow(g) > 0L
    # Describe the genotypes that are actually returned: derive the allele
    # frequency from g itself rather than af_store. For founders and observed
    # generations af_store already equals colMeans(g)/2, so this is a no-op
    # there; for synthetic generations it avoids the mismatch where af_store
    # holds the drift/selection/mutation-evolved frequency (which the stored
    # genotypes do not reflect), keeping obs/exp/FIS/MAF/frac_fixed mutually
    # consistent. The evolved trajectory remains available in `allele_freqs`.
    p   <- if (has_indiv) colMeans(g, na.rm = TRUE) / 2 else af_store[i, ]
    obs <- if (has_indiv) mean(g == 1L, na.rm = TRUE) else NA_real_
    exp <- mean(2 * p * (1 - p), na.rm = TRUE)
    data.frame(
      generation         = i - 1L,
      n_individuals      = nrow(g),
      n_snps             = ncol(g),
      obs_heterozygosity = round(obs, 5),
      exp_heterozygosity = round(exp, 5),
      inbreeding_fis     = round(if (!is.na(obs) && !is.na(exp) && exp > 0)
                                   1 - obs/exp else NA_real_, 5),
      mean_maf           = round(mean(pmin(p, 1-p), na.rm=TRUE), 5),
      frac_fixed         = round(mean(p <= 0 | p >= 1, na.rm=TRUE), 5),
      mean_dosage        = if (has_indiv) round(mean(g, na.rm=TRUE), 5) else NA_real_,
      source             = if ((i-1L) <= max_obs_gen) "observed_pedigree" else "synthetic",
      stringsAsFactors   = FALSE
    )
  }))
  
  if (verbose) {
    message("\n=== Pedigree Simulation Complete ===")
    print(summary_stats[, c("generation","n_individuals","obs_heterozygosity",
                            "mean_maf","inbreeding_fis","source")], row.names = FALSE)
  }
  
  list(genotypes = geno_store, pedigree = ped, snp_map = snp_map,
       allele_freqs = af_store, summary_stats = summary_stats,
       haplotypes = hap_store, simulation_log = sim_log,
       params = list(extra_generations = extra_generations, mut_rate = mut_rate,
                     selection_s = selection_s, n_eff = n_eff, seed = seed,
                     population = population, n_snps = n_snps, max_obs_gen = max_obs_gen))
}

# =============================================================================
# Simulation study to quantify phasing error rate (commented out, run manually)
# =============================================================================

# To run the error rate quantification, uncomment and run the following code.
# It simulates a 3‑generation pedigree with known truth, then compares random
# vs pedigree‑informed phasing.

if (FALSE) {
  # Dummy data for testing: 100 SNPs, 3 generations, 2 parents per generation
  # (You would replace this with your own VCF and pedigree)
  set.seed(123)
  n_snps <- 100
  snp_map <- data.frame(snp_id = paste0("snp",1:n_snps),
                        chrom = "1",
                        pos_bp = 1:n_snps,
                        cohort_maf = runif(n_snps, 0.1, 0.4))
  # Generate truth haplotypes for founders
  truth_hap_f1 <- matrix(rbinom(2*n_snps, 1, 0.3), nrow=2, ncol=n_snps)
  truth_hap_m1 <- matrix(rbinom(2*n_snps, 1, 0.3), nrow=2, ncol=n_snps)
  truth_geno_f1 <- truth_hap_f1[1,] + truth_hap_f1[2,]
  truth_geno_m1 <- truth_hap_m1[1,] + truth_hap_m1[2,]
  # Transmit to child
  g1 <- truth_hap_f1[1+sample(0:1,1),]  # random gamete (simplified)
  g2 <- truth_hap_m1[1+sample(0:1,1),]
  truth_geno_child <- g1 + g2
  # Now apply random phasing and pedigree-informed phasing
  random_phase <- .phase_genotype(truth_geno_child)
  # For pedigree-informed we need parents' phases (here we use the true haplotypes)
  informed_phase <- .phase_pedigree_informed(
    child_dos = truth_geno_child,
    father_dos = truth_geno_f1, mother_dos = truth_geno_m1,
    father_hap1 = truth_hap_f1[1,], father_hap2 = truth_hap_f1[2,],
    mother_hap1 = truth_hap_m1[1,], mother_hap2 = truth_hap_m1[2,]
  )
  # Calculate switch error rate: number of SNPs where the assigned haplotype pair
  # is not consistent with the truth up to a global flip.
  # For simplicity, we compare the two haplotypes concatenated.
  error_random <- mean(c(random_phase$hap1, random_phase$hap2) != c(g1, g2))
  error_informed <- mean(c(informed_phase$hap1, informed_phase$hap2) != c(g1, g2))
  cat("Switch error rate (random phasing):", error_random, "\n")
  cat("Switch error rate (pedigree-informed):", error_informed, "\n")
  # Typically, pedigree-informed reduces errors by >50% at heterozygous sites
  # where a parent is homozygous.
}