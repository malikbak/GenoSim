# =============================================================================
#  GenoSim R package ??? pedigree.R
#  Pedigree reader, validator, inbreeding computation
# =============================================================================

#' Read and Validate a Pedigree CSV File
#'
#' @description
#' Reads a pedigree file in CSV format, normalises column names, validates the
#' parent-offspring structure, infers generation depths via topological sort,
#' and computes per-individual Wright inbreeding coefficients using a recursive
#' kinship table.
#'
#' @param ped_path Character. Path to the pedigree CSV file. The file must
#'   contain at minimum three columns identifying each individual and its
#'   parents (see Details for accepted column names). Missing parents should
#'   be coded as \code{"0"}.
#' @param verbose Logical. Print a summary after loading. Default \code{TRUE}.
#'
#' @return A \code{data.frame} with one row per individual and the following
#'   columns (additional columns from the input file are preserved):
#' \describe{
#'   \item{\code{individual_id}}{Character. Unique individual identifier.}
#'   \item{\code{father_id}}{Character. Father's ID or \code{"0"}.}
#'   \item{\code{mother_id}}{Character. Mother's ID or \code{"0"}.}
#'   \item{\code{sex}}{Integer. 1 = male, 2 = female, 0 = unknown.}
#'   \item{\code{phenotype}}{Integer. 1 = unaffected, 2 = affected,
#'     \code{-9} = missing.}
#'   \item{\code{family_id}}{Character. Family identifier.}
#'   \item{\code{generation}}{Integer. Generation depth (founders = 0),
#'     inferred topologically if not present in the input file.}
#'   \item{\code{is_founder}}{Logical. \code{TRUE} if both parents are absent.}
#'   \item{\code{inbreeding_F}}{Numeric. Wright's inbreeding coefficient
#'     \eqn{F_i = \phi(father_i, mother_i)}, computed by recursive kinship.}
#' }
#'
#' @details
#' \strong{Accepted column name aliases} (case-insensitive):
#' \tabular{ll}{
#'   \strong{Canonical} \tab \strong{Accepted alternatives} \cr
#'   \code{individual_id} \tab \code{iid}, \code{id}, \code{sample_id}, \code{ind_id} \cr
#'   \code{father_id}     \tab \code{fid}, \code{father}, \code{pat}, \code{paternal_id} \cr
#'   \code{mother_id}     \tab \code{mid}, \code{mother}, \code{mat}, \code{maternal_id} \cr
#'   \code{sex}           \tab \code{gender} \cr
#'   \code{phenotype}     \tab \code{pheno}, \code{affection}, \code{status} \cr
#'   \code{family_id}     \tab \code{famid}, \code{fam_id}, \code{family} \cr
#'   \code{generation}    \tab \code{gen}, \code{generation_num}
#' }
#'
#' \strong{Missing parent coding:} \code{"0"}, \code{"NA"}, \code{""}, or
#' \code{"."} are all treated as absent parents.
#'
#' \strong{Missing sex / phenotype:} Blank cells in the \code{sex} and
#' \code{phenotype} columns are coerced to their documented missing codes
#' (\code{sex = 0} unknown, \code{phenotype = -9} missing) at load time, so
#' downstream functions never encounter \code{NA} in these fields.
#'
#' \strong{Unique IDs:} \code{individual_id} must be non-empty and unique.
#' Duplicate or blank IDs raise an informative error rather than failing later
#' with an opaque coercion message.
#'
#' \strong{Inbreeding coefficient:} Computed via the tabular kinship method.
#' \eqn{F_i = \phi(father_i, mother_i)} where \eqn{\phi} is the kinship
#' (coefficient of coancestry) computed recursively from founders upward.
#' Founders have \eqn{F = 0} and \eqn{\phi(i,i) = 0.5}.
#'
#' @examples
#' \dontrun{
#' ped <- read_pedigree("family_pedigree.csv")
#' # Show structure
#' summarise_pedigree(ped)
#' # Affected individuals with high inbreeding
#' ped[ped$phenotype == 2 & ped$inbreeding_F > 0.1, ]
#' }
#'
#' @seealso \code{\link{summarise_pedigree}}, \code{\link{read_vcf_cohort}},
#'   \code{\link{simulate_from_pedigree}}, \code{\link{plot_pedigree_tree}}
#' @export
read_pedigree <- function(ped_path, verbose = TRUE) {
  stopifnot(file.exists(ped_path))
  ped       <- utils::read.csv(ped_path, stringsAsFactors = FALSE,
                                na.strings = c("NA", ""))
  names(ped) <- tolower(trimws(names(ped)))

  colmap <- list(
    individual_id = c("individual_id","iid","id","sample_id","ind_id"),
    father_id     = c("father_id","fid","father","pat","paternal_id"),
    mother_id     = c("mother_id","mid","mother","mat","maternal_id"),
    sex           = c("sex","gender"),
    phenotype     = c("phenotype","pheno","affection","status"),
    family_id     = c("family_id","famid","fam_id","family"),
    generation    = c("generation","gen","generation_num")
  )

  for (canon in names(colmap)) {
    hit <- intersect(colmap[[canon]], names(ped))
    if (length(hit) > 0 && !canon %in% names(ped))
      ped[[canon]] <- ped[[hit[1]]]
  }

  missing_req <- setdiff(c("individual_id","father_id","mother_id"), names(ped))
  if (length(missing_req) > 0)
    stop("Pedigree missing required columns: ", paste(missing_req, collapse = ", "))

  ped$individual_id <- as.character(ped$individual_id)
  ped$father_id     <- as.character(ped$father_id)
  ped$mother_id     <- as.character(ped$mother_id)
  ped$father_id[is.na(ped$father_id)] <- "0"
  ped$mother_id[is.na(ped$mother_id)] <- "0"

  # ---- Validate individual IDs (fail fast with an informative message) ----
  blank_id <- is.na(ped$individual_id) | ped$individual_id %in% .MISSING_PARENT
  if (any(blank_id))
    stop("Pedigree contains ", sum(blank_id), " row(s) with a missing/blank ",
         "individual_id. Every individual must have a non-empty, unique ID.")
  dup_ids <- unique(ped$individual_id[duplicated(ped$individual_id)])
  if (length(dup_ids) > 0)
    stop("Duplicate individual_id value(s) in pedigree: ",
         paste(dup_ids, collapse = ", "),
         ". Each individual must appear on exactly one row with a unique ID.")

  if (!"sex"       %in% names(ped)) ped$sex       <- 0L
  if (!"phenotype" %in% names(ped)) ped$phenotype <- -9L
  if (!"family_id" %in% names(ped)) ped$family_id <- "FAM001"
  # Coerce to the documented integer codes. Blank cells (read as NA) become the
  # documented missing codes (sex = 0 unknown, phenotype = -9 missing) so that
  # downstream code (e.g. plot_pedigree_tree) never sees NA in these fields.
  ped$sex       <- suppressWarnings(as.integer(ped$sex))
  ped$phenotype <- suppressWarnings(as.integer(ped$phenotype))
  ped$sex[is.na(ped$sex)]             <- 0L
  ped$phenotype[is.na(ped$phenotype)] <- -9L

  ped$is_founder <- ped$father_id %in% .MISSING_PARENT &
                    ped$mother_id %in% .MISSING_PARENT

  self_ref <- ped$individual_id == ped$father_id |
              ped$individual_id == ped$mother_id
  if (any(self_ref, na.rm = TRUE))
    stop("Self-reference in pedigree: ",
         paste(ped$individual_id[self_ref], collapse = ", "))

  all_ids     <- ped$individual_id
  bad_fathers <- ped$father_id[!ped$father_id %in% .MISSING_PARENT &
                                !ped$father_id %in% all_ids]
  bad_mothers <- ped$mother_id[!ped$mother_id %in% .MISSING_PARENT &
                                !ped$mother_id %in% all_ids]
  if (length(bad_fathers) > 0)
    warning("Father IDs not in individual list: ",
            paste(unique(bad_fathers), collapse = ", "))
  if (length(bad_mothers) > 0)
    warning("Mother IDs not in individual list: ",
            paste(unique(bad_mothers), collapse = ", "))

  if (!"generation" %in% names(ped))
    ped$generation <- .topo_generation(ped)

  ped$inbreeding_F <- .compute_inbreeding(ped)

  if (verbose) {
    message(sprintf("[Pedigree] %d individuals | %d founders | %d generations | mean F=%.4f",
                    nrow(ped), sum(ped$is_founder),
                    max(ped$generation) + 1L,
                    mean(ped$inbreeding_F[!ped$is_founder], na.rm = TRUE)))
  }
  ped
}

#' Print a Structured Summary of a Pedigree
#'
#' @description
#' Prints a formatted summary of a pedigree data frame, including per-generation
#' counts, sex distribution, phenotype counts, inbreeding statistics, and
#' mating pair counts per generation.
#'
#' @param ped A \code{data.frame} returned by \code{\link{read_pedigree}}.
#'
#' @return The input \code{ped} data frame, invisibly.
#'
#' @examples
#' \dontrun{
#' ped <- read_pedigree("family_pedigree.csv")
#' summarise_pedigree(ped)
#' }
#'
#' @seealso \code{\link{read_pedigree}}
#' @export
summarise_pedigree <- function(ped) {
  cat("=== Pedigree Summary ===\n")
  cat(sprintf("  Individuals : %d\n", nrow(ped)))
  cat(sprintf("  Families    : %s\n", paste(unique(ped$family_id), collapse = ", ")))
  gen_tab <- table(ped$generation)
  cat("  Per generation:\n")
  for (g in names(gen_tab))
    cat(sprintf("    Gen %s: %d (%d founders)\n", g, gen_tab[g],
                sum(ped$is_founder & ped$generation == as.integer(g))))
  stab <- table(ped$sex)
  cat(sprintf("  Sex: %d male | %d female | %d unknown\n",
              stab["1"] %||% 0L, stab["2"] %||% 0L, stab["0"] %||% 0L))
  ptab <- table(ped$phenotype)
  cat(sprintf("  Phenotype: %d unaffected | %d affected | %d missing\n",
              ptab["1"] %||% 0L, ptab["2"] %||% 0L, ptab["-9"] %||% 0L))
  cat(sprintf("  Inbreeding F: mean=%.4f | max=%.4f\n",
              mean(ped$inbreeding_F, na.rm = TRUE),
              max(ped$inbreeding_F, na.rm = TRUE)))
  invisible(ped)
}

#' Extract Mating Pairs from a Pedigree
#'
#' @description
#' Returns a data frame of unique (father, mother) pairs optionally filtered
#' to a specific generation of offspring.
#'
#' @param ped A \code{data.frame} from \code{\link{read_pedigree}}.
#' @param generation Integer or \code{NULL}. If supplied, only pairs whose
#'   offspring belong to this generation are returned.
#'
#' @return A \code{data.frame} with columns \code{father_id}, \code{mother_id},
#'   \code{family_id}, and \code{offspring} (semicolon-separated IDs).
#'
#' @examples
#' \dontrun{
#' ped   <- read_pedigree("family_pedigree.csv")
#' pairs <- extract_mating_pairs(ped, generation = 1)
#' }
#'
#' @seealso \code{\link{read_pedigree}}
#' @export
extract_mating_pairs <- function(ped, generation = NULL) {
  if (!is.null(generation)) ped <- ped[ped$generation == generation, ]
  nf  <- ped[!ped$is_founder, ]
  if (nrow(nf) == 0) return(NULL)
  pairs <- unique(nf[, c("father_id","mother_id","family_id")])
  pairs <- pairs[!pairs$father_id %in% .MISSING_PARENT &
                 !pairs$mother_id %in% .MISSING_PARENT, ]
  pairs$offspring <- mapply(function(f, m)
    paste(nf$individual_id[nf$father_id == f & nf$mother_id == m],
          collapse = ";"),
    pairs$father_id, pairs$mother_id)
  pairs
}

# ---- Internal: topological generation assignment ----------------------------
#' @noRd
.topo_generation <- function(ped) {
  gen  <- stats::setNames(rep(NA_integer_, nrow(ped)), ped$individual_id)
  gen[ped$is_founder] <- 0L
  for (iter in seq_len(nrow(ped) + 1L)) {
    if (!any(is.na(gen))) break
    for (i in which(is.na(gen))) {
      fg <- if (ped$father_id[i] %in% .MISSING_PARENT) 0L else gen[ped$father_id[i]]
      mg <- if (ped$mother_id[i] %in% .MISSING_PARENT) 0L else gen[ped$mother_id[i]]
      if (!is.na(fg) && !is.na(mg)) gen[i] <- max(fg, mg) + 1L
    }
  }
  if (any(is.na(gen))) {
    warning("Could not resolve generation for: ",
            paste(names(gen)[is.na(gen)], collapse = ", "), ". Assigning max+1.")
    gen[is.na(gen)] <- max(gen, na.rm = TRUE) + 1L
  }
  unname(gen)
}

# ---- Internal: kinship-based inbreeding computation ------------------------
#' @noRd
.compute_inbreeding <- function(ped) {
  ids <- ped$individual_id
  n   <- length(ids)
  phi <- matrix(0.0, n, n, dimnames = list(ids, ids))
  diag(phi) <- 0.5

  for (k in order(ped$generation)) {
    i  <- ids[k]
    fi <- ped$father_id[k]; mi <- ped$mother_id[k]
    hf <- !fi %in% .MISSING_PARENT && fi %in% ids
    hm <- !mi %in% .MISSING_PARENT && mi %in% ids
    if (hf && hm) phi[i, i] <- 0.5 + 0.5 * phi[fi, mi]
    for (l in order(ped$generation)) {
      j <- ids[l]
      if (i == j) next
      val <- 0.0
      if (hf) val <- val + 0.5 * phi[fi, j]
      if (hm) val <- val + 0.5 * phi[mi, j]
      phi[i, j] <- phi[j, i] <- val
    }
  }

  vapply(seq_len(n), function(k) {
    fi <- ped$father_id[k]; mi <- ped$mother_id[k]
    if (fi %in% .MISSING_PARENT || mi %in% .MISSING_PARENT ||
        !fi %in% ids || !mi %in% ids) return(0.0)
    round(phi[fi, mi], 5)
  }, numeric(1))
}
