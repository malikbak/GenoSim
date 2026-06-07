# =============================================================================
#  GenoSim R package ??? vcf_reader.R
#  Read single multi-sample VCF or folder of per-individual VCFs
# =============================================================================

#' Parse a Single VCF File
#'
#' @description Internal helper that reads one VCF file and returns a dosage
#'   matrix plus SNP map. Handles phased and unphased GT fields, biallelic
#'   SNP filtering, and PASS/QUAL filtering.
#'
#' @param vcf_path Character. Full path to the VCF file.
#' @param min_qual Numeric or \code{NULL}. Minimum QUAL score; \code{NULL}
#'   skips quality filtering.
#' @param pass_filter_only Logical. If \code{TRUE} (default), only variants
#'   with \code{FILTER == "PASS"} or \code{FILTER == "."} are retained.
#' @param biallelic_only Logical. If \code{TRUE} (default), multi-allelic
#'   sites and indels are removed.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{geno_matrix}}{Integer matrix \eqn{n_{samples} \times n_{SNPs}}
#'     of dosage values (0/1/2). \code{NA} indicates missing genotype.}
#'   \item{\code{snp_map}}{A \code{data.frame} with columns \code{snp_id},
#'     \code{chrom}, \code{pos_bp}, \code{ref}, \code{alt}.}
#' }
#'
#' @keywords internal
#' @noRd
.parse_vcf_file <- function(vcf_path, min_qual = NULL,
                             pass_filter_only = TRUE,
                             biallelic_only   = TRUE) {
  stopifnot(file.exists(vcf_path))
  message(sprintf("  Reading VCF: %s", basename(vcf_path)))

  lines       <- readLines(vcf_path, warn = FALSE)
  header_line <- grep("^#CHROM", lines, value = TRUE)
  data_lines  <- lines[!grepl("^#", lines)]

  if (length(header_line) == 0) stop("No #CHROM header found in: ", vcf_path)
  if (length(data_lines)  == 0) stop("No variant records in: ", vcf_path)

  col_names    <- strsplit(header_line, "\t")[[1]]
  col_names[1] <- "CHROM"
  fixed_cols   <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  sample_ids   <- setdiff(col_names, fixed_cols)

  mat <- do.call(rbind, strsplit(data_lines, "\t"))
  if (ncol(mat) != length(col_names))
    stop("Column count mismatch in VCF body vs header: ", vcf_path)
  colnames(mat) <- col_names
  df <- as.data.frame(mat, stringsAsFactors = FALSE)

  if (pass_filter_only && "FILTER" %in% names(df))
    df <- df[df$FILTER %in% c("PASS", ".", ""), ]

  if (!is.null(min_qual) && "QUAL" %in% names(df)) {
    qn <- suppressWarnings(as.numeric(df$QUAL))
    df <- df[is.na(qn) | qn >= min_qual, ]
  }

  df$CHROM <- sub("^chr", "", df$CHROM)
  df       <- df[df$CHROM %in% as.character(1:22), ]

  if (biallelic_only) {
    df <- df[!grepl(",", df$ALT), ]
    df <- df[nchar(df$REF) == 1 & nchar(df$ALT) == 1, ]
  }

  if (nrow(df) == 0) stop("No variants passed filters in: ", vcf_path)

  snp_map <- data.frame(
    snp_id = ifelse(df$ID %in% c(".", ""),
                    paste0("chr", df$CHROM, "_", df$POS), df$ID),
    chrom  = as.integer(df$CHROM),
    pos_bp = as.integer(df$POS),
    ref    = df$REF,
    alt    = df$ALT,
    stringsAsFactors = FALSE
  )

  dup_pos  <- duplicated(paste0(snp_map$chrom, "_", snp_map$pos_bp))
  df       <- df[!dup_pos, ]
  snp_map  <- snp_map[!dup_pos, ]

  fmt_fields <- strsplit(df$FORMAT[1], ":")[[1]]
  gt_idx     <- which(fmt_fields == "GT")
  if (length(gt_idx) == 0) stop("No GT field in FORMAT: ", vcf_path)

  .gt_to_dosage <- function(gt_str) {
    gt_str <- gsub("\\|", "/", gt_str)
    gt_str <- vapply(strsplit(gt_str, ":"), `[[`, character(1), gt_idx)
    gt_map <- c("0/0" = 0L, "0/1" = 1L, "1/0" = 1L, "1/1" = 2L,
                "./." = NA_integer_, "."  = NA_integer_)
    out <- gt_map[gt_str]
    out[is.na(names(out))] <- NA_integer_
    as.integer(out)
  }

  if (length(sample_ids) > 0) {
    geno_mat <- vapply(sample_ids, function(s) .gt_to_dosage(df[[s]]),
                       integer(nrow(df)))
    if (is.vector(geno_mat)) geno_mat <- matrix(geno_mat, ncol = 1)
    rownames(geno_mat) <- snp_map$snp_id
    geno_mat <- t(geno_mat)
    rownames(geno_mat) <- sample_ids
  } else {
    sid      <- tools::file_path_sans_ext(basename(vcf_path))
    geno_mat <- matrix(NA_integer_, nrow = 1, ncol = nrow(snp_map),
                       dimnames = list(sid, snp_map$snp_id))
  }

  list(geno_matrix = geno_mat, snp_map = snp_map)
}

#' Merge Multiple VCF Parse Results on Common SNPs
#' @noRd
.merge_vcf_results <- function(vcf_results) {
  common_snps <- Reduce(intersect, lapply(vcf_results, function(r) r$snp_map$snp_id))
  if (length(common_snps) == 0)
    stop("No SNPs in common across all VCF files. Check genome build consistency.")
  message(sprintf("  Common SNPs across all VCFs: %d", length(common_snps)))
  merged_geno <- do.call(rbind, lapply(vcf_results, function(r)
    r$geno_matrix[, common_snps, drop = FALSE]))
  ref_map <- vcf_results[[1]]$snp_map
  list(geno_matrix = merged_geno,
       snp_map     = ref_map[match(common_snps, ref_map$snp_id), ])
}

#' Impute Missing Genotypes from Allele Frequency
#'
#' Fills every missing genotype by an HWE/inbreeding draw at the per-site
#' observed allele frequency. Sites with no observed calls fall back to the
#' global mean allele frequency (then to 0.5), so the returned matrix is
#' guaranteed to contain no \code{NA}.
#' @noRd
.impute_missing <- function(geno_matrix, F_coef = 0) {
  n_miss <- sum(is.na(geno_matrix))
  if (n_miss == 0) return(geno_matrix)
  message(sprintf("  Imputing %d missing genotypes (%.1f%%)...",
                  n_miss, 100 * n_miss / length(geno_matrix)))
  global_p <- mean(geno_matrix, na.rm = TRUE) / 2
  if (is.nan(global_p) || is.na(global_p)) global_p <- 0.5
  draw <- function(n, p) {
    q  <- 1 - p
    p0 <- q^2 * (1 - F_coef) + q * F_coef
    p1 <- 2 * p * q * (1 - F_coef)
    u  <- stats::runif(n)
    ifelse(u < p0, 0L, ifelse(u < p0 + p1, 1L, 2L))
  }
  for (j in seq_len(ncol(geno_matrix))) {
    miss <- which(is.na(geno_matrix[, j]))
    if (length(miss) == 0) next
    obs  <- geno_matrix[-miss, j]
    p    <- if (length(obs) > 0) mean(obs, na.rm = TRUE) / 2 else NA_real_
    if (is.nan(p) || is.na(p)) p <- global_p   # all-missing site -> global AF
    geno_matrix[miss, j] <- draw(length(miss), p)
  }
  # Guarantee completeness: if anything slipped through (e.g. a column that was
  # entirely NA and global_p was also undefined), fill from the global AF.
  resid <- which(is.na(geno_matrix))
  if (length(resid) > 0) {
    warning(sprintf("Imputation: %d genotype(s) had no usable frequency; ",
                    length(resid)),
            "filled from the global allele frequency to complete the matrix.")
    geno_matrix[resid] <- draw(length(resid), global_p)
  }
  storage.mode(geno_matrix) <- "integer"
  geno_matrix
}

# ---- Public API -------------------------------------------------------------

#' Read a Family VCF Cohort into a Dosage Matrix
#'
#' @description
#' Reads genotype data from either a single multi-sample VCF file or a
#' directory containing one VCF file per individual. Variant records are
#' filtered to biallelic autosomal SNPs, merged across files on their common
#' set of positions, and optionally imputed for missing genotypes.
#'
#' @param vcf_input Character. Either:
#'   \itemize{
#'     \item Path to a single multi-sample VCF file (all samples as columns), or
#'     \item Path to a directory containing one \code{.vcf} or \code{.vcf.gz}
#'       file per individual. The filename stem (without extension) is used as
#'       the sample identifier and must match \code{individual_id} in the
#'       pedigree.
#'   }
#' @param min_qual Numeric or \code{NULL}. Minimum QUAL score threshold.
#'   Variants below this value are excluded. Default \code{NULL} (no filter).
#' @param pass_filter_only Logical. Retain only \code{FILTER == "PASS"} or
#'   \code{"."} variants. Default \code{TRUE}.
#' @param biallelic_only Logical. Keep only biallelic single-nucleotide
#'   variants. Multi-allelic sites and indels are dropped. Default \code{TRUE}.
#' @param impute_missing Logical. If \code{TRUE} (default), missing genotypes
#'   (\code{./.}) are imputed by sampling from the observed allele frequency
#'   at that locus under the \code{F_coef} inbreeding model. When \code{TRUE}
#'   the returned matrix is guaranteed to contain no \code{NA}.
#' @param F_coef Numeric in \eqn{[0,1)}. Inbreeding coefficient used when
#'   imputing missing genotypes. Default \code{0}.
#' @param max_missing_rate Numeric in \eqn{(0,1]}. Sites whose per-SNP
#'   missingness exceeds this fraction are dropped \emph{before} imputation.
#'   This is the recommended way to handle raw WES/WGS callsets, where many
#'   sites are missing in most samples and imputing them from a handful of
#'   observed calls is unreliable. Default \code{1} (keep all sites; preserves
#'   prior behaviour). A typical quality threshold is \code{0.1}.
#' @param max_snps Integer. Maximum number of SNPs to retain. If the merged
#'   matrix exceeds this limit, SNPs are randomly subsampled. Default
#'   \code{100000}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{geno_matrix}}{Integer matrix \eqn{n_{individuals} \times
#'     n_{SNPs}} of dosage values (0 = hom-ref, 1 = het, 2 = hom-alt).
#'     Rownames are sample identifiers; colnames are SNP IDs.}
#'   \item{\code{snp_map}}{A \code{data.frame} with columns \code{snp_id},
#'     \code{chrom}, \code{pos_bp}, \code{ref}, \code{alt}, and
#'     \code{cohort_maf} (observed minor allele frequency in this cohort).}
#' }
#'
#' @section Expected input and size limits:
#' This reader is a pure-R parser intended for \strong{targeted / array-style}
#' callsets (family panels, gene panels, genotyping arrays) with relatively
#' complete genotypes and up to roughly \eqn{10^5} variants. It reads the whole
#' file into memory, so raw whole-exome (WES) or whole-genome (WGS) VCFs with
#' millions of records (hundreds of MB) are impractical to load directly and
#' will be slow or exhaust memory. For such inputs, pre-filter first with
#' \code{\link{prefilter_vcf}} (which streams the file without loading it whole)
#' or an external tool such as \code{bcftools}, and/or set
#' \code{max_missing_rate} to drop the very sparse sites typical of raw
#' WES/WGS before imputation. A warning is emitted when the input is large.
#'
#' @examples
#' \dontrun{
#' # From a directory of per-individual VCFs
#' vcf <- read_vcf_cohort("path/to/family_vcfs/")
#'
#' # From a single multi-sample VCF
#' vcf <- read_vcf_cohort("path/to/family.vcf", min_qual = 30)
#'
#' # With imputation using pedigree-derived F
#' vcf <- read_vcf_cohort("path/to/family_vcfs/",
#'                         impute_missing = TRUE, F_coef = 0.125)
#' }
#'
#' @seealso \code{\link{read_pedigree}}, \code{\link{simulate_from_pedigree}}
#' @export
read_vcf_cohort <- function(vcf_input,
                             min_qual         = NULL,
                             pass_filter_only = TRUE,
                             biallelic_only   = TRUE,
                             impute_missing   = TRUE,
                             F_coef           = 0,
                             max_missing_rate = 1,
                             max_snps         = 100000L,
                             verbose          = TRUE) {
  if (verbose) message("[VCF] Loading cohort genotypes...")

  # Warn about inputs that are too large for an in-memory pure-R parse.
  big_files <- character(0)
  if (file.exists(vcf_input) && !dir.exists(vcf_input))
    big_files <- vcf_input
  else if (dir.exists(vcf_input))
    big_files <- list.files(vcf_input, pattern = "\\.vcf(\\.gz)?$",
                            full.names = TRUE, ignore.case = TRUE)
  total_mb <- sum(file.info(big_files)$size, na.rm = TRUE) / 1024^2
  if (isTRUE(total_mb > 200))
    warning(sprintf(paste0("read_vcf_cohort(): input is ~%.0f MB. This pure-R ",
            "reader loads the whole file into memory and is impractical for ",
            "large WES/WGS VCFs. Consider prefilter_vcf() or bcftools first."),
            total_mb))

  if (dir.exists(vcf_input)) {
    vcf_files <- list.files(vcf_input, pattern = "\\.vcf(\\.gz)?$",
                            full.names = TRUE, ignore.case = TRUE)
    if (length(vcf_files) == 0)
      stop("No .vcf or .vcf.gz files found in: ", vcf_input)
    if (verbose) message(sprintf("  Found %d VCF files", length(vcf_files)))
    results <- lapply(vcf_files, .parse_vcf_file,
                      min_qual         = min_qual,
                      pass_filter_only = pass_filter_only,
                      biallelic_only   = biallelic_only)
    merged  <- .merge_vcf_results(results)
  } else if (file.exists(vcf_input)) {
    merged <- .parse_vcf_file(vcf_input,
                              min_qual         = min_qual,
                              pass_filter_only = pass_filter_only,
                              biallelic_only   = biallelic_only)
  } else {
    stop("vcf_input path does not exist: ", vcf_input)
  }

  geno_mat <- merged$geno_matrix
  snp_map  <- merged$snp_map

  # Drop sites that are missing in more than `max_missing_rate` of samples
  # (recommended for raw WES/WGS) BEFORE imputation, so imputed values come
  # from reasonably observed loci rather than a handful of calls.
  if (max_missing_rate < 1) {
    miss_rate <- colMeans(is.na(geno_mat))
    keep_site <- miss_rate <= max_missing_rate
    n_drop    <- sum(!keep_site)
    if (n_drop > 0) {
      geno_mat <- geno_mat[, keep_site, drop = FALSE]
      snp_map  <- snp_map[keep_site, , drop = FALSE]
      if (verbose)
        message(sprintf("  Dropped %d SNP(s) with >%.0f%% missingness; %d remain",
                        n_drop, 100 * max_missing_rate, ncol(geno_mat)))
    }
    if (ncol(geno_mat) == 0)
      stop("All sites exceeded max_missing_rate = ", max_missing_rate,
           ". Relax the threshold or check input genotype completeness.")
  }

  if (ncol(geno_mat) > max_snps) {
    idx      <- sort(sample(ncol(geno_mat), max_snps))
    geno_mat <- geno_mat[, idx, drop = FALSE]
    snp_map  <- snp_map[idx, ]
    if (verbose) message(sprintf("  Subsampled to %d SNPs", max_snps))
  }

  if (impute_missing)
    geno_mat <- .impute_missing(geno_mat, F_coef = F_coef)

  p_obs           <- colMeans(geno_mat, na.rm = TRUE) / 2
  snp_map$cohort_maf <- pmin(p_obs, 1 - p_obs)

  if (verbose)
    message(sprintf("  Cohort loaded: %d individuals, %d SNPs | mean MAF=%.4f",
                    nrow(geno_mat), ncol(geno_mat), mean(snp_map$cohort_maf)))

  list(geno_matrix = geno_mat, snp_map = snp_map)
}
