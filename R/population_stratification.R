# =============================================================================
# GenoSim: Population Stratification Module
# Balding-Nichols model for population-specific allele frequencies
# =============================================================================

#' Return gnomAD-derived population parameters
#'
#' @description
#' Provides Fst values and reference ancestral allele frequencies for 8 ancestry
#' groups based on gnomAD v4.0 superpopulations. The Fst values are relative to
#' the global (ancestral) population.
#'
#' @param ancestral_af Optional numeric vector of ancestral allele frequencies
#'   (length = number of SNPs). If `NULL`, returns only Fst values.
#'
#' @return A list with components:
#'   \item{Fst}{Named numeric vector of Fst values for each population.}
#'   \item{ancestral_af}{The input `ancestral_af` vector (if provided).}
#'   \item{populations}{Character vector of population names.}
#'
#' @examples
#' # Get Fst only
#' params <- reference_populations()
#' print(params$Fst)
#'
#' # With ancestral frequencies
#' af <- runif(1000)
#' params <- reference_populations(ancestral_af = af)
#'
#' @export
reference_populations <- function(ancestral_af = NULL) {
  # gnomAD v4.0 superpopulation Fst values (relative to global mean)
  # Approximate values from literature (Bergström et al. 2021, Nature)
  fst <- c(
    AFR = 0.152,  # African
    AMR = 0.043,  # Admixed American
    EAS = 0.086,  # East Asian
    SAS = 0.075,  # South Asian
    EUR = 0.042,  # European (non-Finnish)
    FIN = 0.051,  # Finnish
    MID = 0.058,  # Middle Eastern (closest to gnomAD "OTH"?)
    OCE = 0.093   # Oceanian
  )
  # Ensure all 8 groups are named
  populations <- names(fst)
  
  if (!is.null(ancestral_af)) {
    if (!is.numeric(ancestral_af))
      stop("ancestral_af must be a numeric vector")
    list(Fst = fst, ancestral_af = ancestral_af, populations = populations)
  } else {
    list(Fst = fst, ancestral_af = NULL, populations = populations)
  }
}

#' Load allele frequencies from a file
#'
#' @description
#' Reads a file containing SNP IDs and corresponding allele frequencies
#' (e.g., ancestral frequencies from a reference panel). Supported formats:
#' CSV, TSV, or RDS.
#'
#' @param file Path to the file. Must contain columns `snp_id` and `af`.
#' @param sep Separator for text files (default: comma for CSV, tab for TSV).
#' @param format Either "csv", "tsv", "rds", or `NULL` to guess from extension.
#'
#' @return A named numeric vector of allele frequencies (names = snp_id).
#'
#' @examples
#' \dontrun{
#' anc_af <- load_af_from_file("gnomad_global_af.csv")
#' }
#'
#' @export
load_af_from_file <- function(file, sep = NULL, format = NULL) {
  if (!file.exists(file)) stop("File not found: ", file)
  
  # Guess format from extension if not provided
  if (is.null(format)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "rds") format <- "rds"
    else if (ext %in% c("csv", "txt")) format <- "csv"
    else if (ext %in% c("tsv", "tab")) format <- "tsv"
    else stop("Unknown file format. Please specify format = 'csv','tsv','rds'")
  }
  
  if (format == "rds") {
    dat <- readRDS(file)
  } else {
    if (is.null(sep)) {
      sep <- if (format == "csv") "," else "\t"
    }
    dat <- utils::read.table(file, header = TRUE, sep = sep,
                             stringsAsFactors = FALSE)
  }
  
  if (!all(c("snp_id", "af") %in% colnames(dat))) {
    stop("File must contain columns 'snp_id' and 'af'")
  }
  
  af_vec <- dat$af
  names(af_vec) <- dat$snp_id
  af_vec
}

#' Simulate a population with stratification using Balding-Nichols model
#'
#' @description
#' Generates a genotype matrix for a single population where allele frequencies
#' are drawn from a Beta distribution parameterised by ancestral frequencies
#' and population-specific Fst. This implements the Balding-Nichols model.
#'
#' @param n_individuals Integer. Number of diploid individuals to simulate.
#' @param n_snps Integer. Number of SNPs (if `ancestral_af` is not provided).
#' @param population Character. Population name (must be one of those returned
#'   by \code{\link{reference_populations}}).
#' @param ancestral_af Either a numeric vector of length `n_snps` or a file path
#'   (passed to \code{\link{load_af_from_file}}). If `NULL`, uniform random
#'   frequencies are used as ancestral (not recommended).
#' @param Fst Override the population-specific Fst. If `NULL`, uses the value
#'   from \code{\link{reference_populations}}.
#' @param inbreeding_F Numeric in [0,1]. Inbreeding coefficient for all
#'   individuals (applied via \code{.sample_genotypes_vec}). Default 0.
#' @param seed Random seed (optional).
#' @param verbose Logical. Print progress.
#'
#' @return A genotype matrix (individuals x SNPs) with dosages 0,1,2.
#'
#' @details
#' The Balding-Nichols model draws population-specific allele frequency
#' \eqn{p_{pop}} for each SNP as:
#' \deqn{p_{pop} \sim \text{Beta}(\alpha, \beta)}
#' where \eqn{\alpha = p_{anc} \cdot (1-F_{st}) / F_{st}} and
#' \eqn{\beta = (1-p_{anc}) \cdot (1-F_{st}) / F_{st}}.
#' Then genotypes are sampled under Hardy‑Weinberg equilibrium (or with
#' inbreeding if `inbreeding_F > 0`).
#'
#' @examples
#' \dontrun{
#' # Simulate 100 individuals from South Asian population using gnomAD ancestral AFs
#' anc <- load_af_from_file("gnomad_global_af.csv")
#' geno <- simulate_population_stratified(
#'   n_individuals = 100,
#'   ancestral_af = anc,
#'   population = "SAS",
#'   seed = 123
#' )
#' }
#'
#' @export
simulate_population_stratified <- function(
    n_individuals,
    n_snps = NULL,
    population,
    ancestral_af = NULL,
    Fst = NULL,
    inbreeding_F = 0.0,
    seed = NULL,
    verbose = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  
  # Validate population
  pop_params <- reference_populations()
  if (!population %in% pop_params$populations) {
    stop("Population '", population, "' not found. Available: ",
         paste(pop_params$populations, collapse = ", "))
  }
  
  # Handle ancestral_af input
  if (is.character(ancestral_af)) {
    if (verbose) message("Loading ancestral frequencies from ", ancestral_af)
    ancestral_af <- load_af_from_file(ancestral_af)
  }
  
  if (is.null(ancestral_af)) {
    if (is.null(n_snps)) stop("Either ancestral_af or n_snps must be provided")
    if (verbose) message("No ancestral AF provided: using uniform(0,1) as ancestral")
    ancestral_af <- stats::runif(n_snps)
  }
  
  n_snps <- length(ancestral_af)
  if (verbose) message("Simulating ", n_individuals, " individuals, ", n_snps, " SNPs")
  
  # Determine Fst
  if (is.null(Fst)) {
    Fst <- pop_params$Fst[population]
    if (verbose) message("Using Fst = ", Fst, " for ", population)
  } else {
    if (Fst < 0 || Fst > 1) stop("Fst must be in [0,1]")
    if (verbose) message("Using user-supplied Fst = ", Fst)
  }
  
  # Draw population-specific allele frequencies using Balding-Nichols
  p_pop <- .draw_population_afs(ancestral_af, Fst)
  if (verbose) message("Population allele frequencies drawn (mean = ", round(mean(p_pop), 3), ")")
  
  # Sample genotypes
  genotypes <- matrix(NA_integer_, n_individuals, n_snps)
  for (i in seq_len(n_individuals)) {
    genotypes[i, ] <- .sample_genotypes_vec(p_pop, inbreeding_F)
  }
  
  # Add row/column names
  rownames(genotypes) <- paste0("ind", seq_len(n_individuals))
  colnames(genotypes) <- if (!is.null(names(ancestral_af))) names(ancestral_af) else paste0("snp", seq_len(n_snps))
  
  if (verbose) message("Genotype matrix complete.")
  genotypes
}

# ------------------------- Internal helpers -------------------------

#' Draw population-specific allele frequencies using Balding-Nichols
#' @noRd
.draw_population_afs <- function(p_anc, Fst) {
  # p_anc: ancestral allele frequencies (vector)
  # Fst: population-specific fixation index
  # Returns vector of population frequencies
  if (Fst == 0) return(p_anc)
  if (Fst == 1) {
    # Fixed – but Beta not defined; treat as either 0 or 1 with prob p_anc
    return(as.numeric(stats::runif(length(p_anc)) < p_anc))
  }
  # Alpha and beta parameters of Beta distribution
  alpha <- p_anc * (1 - Fst) / Fst
  beta  <- (1 - p_anc) * (1 - Fst) / Fst
  # Draw one frequency per SNP
  stats::rbeta(length(p_anc), alpha, beta)
}

# Note: .sample_genotypes_vec is assumed to exist in the package namespace
# (it is defined in simulate_genotypes.R or pedigree_sim.R)