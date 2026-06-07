# =============================================================================
#  GenoSim R package - prefilter.R
#  Streaming pre-filter for large VCF files (does not load the whole file)
# =============================================================================

#' Pre-filter a Large VCF Without Loading It Into Memory
#'
#' @description
#' Streams a (optionally gzipped) VCF file in fixed-size chunks and writes a
#' smaller VCF containing only the records that pass the requested filters.
#' Because the file is processed line by line, it can handle inputs far larger
#' than available memory (e.g. multi-million-variant WES/WGS callsets) that
#' \code{\link{read_vcf_cohort}} cannot parse directly. The reduced output is
#' then suitable input for \code{read_vcf_cohort()}.
#'
#' Optionally, if \code{use_bcftools = TRUE} and a \code{bcftools} executable is
#' on the \code{PATH}, the filtering is delegated to \code{bcftools view} for
#' speed; otherwise a portable pure-R streaming implementation is used.
#'
#' @param vcf_in Character. Path to the input \code{.vcf} or \code{.vcf.gz}.
#' @param vcf_out Character. Path for the filtered output VCF. If it ends in
#'   \code{.gz} the output is gzip-compressed.
#' @param pass_only Logical. Keep only records with \code{FILTER} equal to
#'   \code{PASS}, \code{.}, or empty. Default \code{TRUE}.
#' @param biallelic_snp_only Logical. Keep only biallelic single-nucleotide
#'   variants (single-character \code{REF} and \code{ALT}, no comma in
#'   \code{ALT}). Default \code{TRUE}.
#' @param autosomes_only Logical. Keep only chromosomes 1-22 (with or without a
#'   \code{chr} prefix). Default \code{TRUE}.
#' @param min_qual Numeric or \code{NULL}. Minimum \code{QUAL}. Default
#'   \code{NULL} (no QUAL filter).
#' @param max_variants Integer or \code{NULL}. If set, at most this many passing
#'   records are written, chosen by uniform reservoir sampling across the whole
#'   file (representative, not just the first N). Default \code{NULL}.
#' @param chunk_lines Integer. Number of lines read per streaming chunk. Larger
#'   values are faster but use more memory. Default \code{1e5}.
#' @param use_bcftools Logical. Use an external \code{bcftools} if available.
#'   Default \code{FALSE}.
#' @param seed Integer or \code{NULL}. Seed for reservoir sampling
#'   reproducibility. Default \code{NULL}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return Invisibly, a list with \code{out} (the output path), \code{n_in}
#'   (data records scanned), and \code{n_out} (records written).
#'
#' @examples
#' \dontrun{
#' # Reduce a huge raw WGS VCF to <=50k PASS biallelic autosomal SNPs:
#' prefilter_vcf("big.vcf.gz", "small.vcf.gz",
#'               min_qual = 30, max_variants = 50000, seed = 1)
#' vcf <- read_vcf_cohort("small.vcf.gz")
#' }
#'
#' @seealso \code{\link{read_vcf_cohort}}
#' @export
prefilter_vcf <- function(vcf_in, vcf_out,
                          pass_only          = TRUE,
                          biallelic_snp_only = TRUE,
                          autosomes_only     = TRUE,
                          min_qual           = NULL,
                          max_variants       = NULL,
                          chunk_lines        = 1e5L,
                          use_bcftools       = FALSE,
                          seed               = NULL,
                          verbose            = TRUE) {
  if (!file.exists(vcf_in)) stop("Input VCF not found: ", vcf_in)
  if (!is.null(seed)) set.seed(seed)

  # ---- Optional fast path: delegate to bcftools --------------------------
  if (isTRUE(use_bcftools) && nzchar(Sys.which("bcftools"))) {
    if (verbose) message("[prefilter] Using bcftools view ...")
    args <- c("view")
    if (pass_only)          args <- c(args, "-f", "PASS,.")
    if (biallelic_snp_only) args <- c(args, "-m2", "-M2", "-v", "snps")
    if (autosomes_only)     args <- c(args, "-r", paste(1:22, collapse = ","))
    if (!is.null(min_qual)) args <- c(args, "-e", sprintf("QUAL<%g", min_qual))
    args <- c(args, shQuote(vcf_in))
    out_con <- if (grepl("\\.gz$", vcf_out)) "-Oz" else "-Ov"
    args <- c(args, out_con, "-o", shQuote(vcf_out))
    status <- system2("bcftools", args)
    if (status != 0)
      warning("bcftools returned status ", status,
              "; falling back to pure-R streaming.")
    else
      return(invisible(list(out = vcf_out, n_in = NA_integer_,
                            n_out = NA_integer_)))
  }

  # ---- Pure-R streaming implementation -----------------------------------
  con_in <- if (grepl("\\.gz$", vcf_in)) gzfile(vcf_in, "rt") else file(vcf_in, "rt")
  on.exit(close(con_in), add = TRUE)
  con_out <- if (grepl("\\.gz$", vcf_out)) gzfile(vcf_out, "wt") else file(vcf_out, "wt")
  on.exit(close(con_out), add = TRUE)

  autos <- as.character(1:22)
  n_in  <- 0L
  n_out <- 0L
  header_done <- FALSE

  # Reservoir state (only used when max_variants is set)
  reservoir <- if (!is.null(max_variants))
    character(max_variants) else NULL
  res_count <- 0L

  .keep_line <- function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < 8L) return(FALSE)
    if (autosomes_only) {
      chrom <- sub("^chr", "", f[1])
      if (!chrom %in% autos) return(FALSE)
    }
    if (biallelic_snp_only) {
      ref <- f[4]; alt <- f[5]
      if (grepl(",", alt, fixed = TRUE)) return(FALSE)
      if (nchar(ref) != 1L || nchar(alt) != 1L) return(FALSE)
    }
    if (pass_only && !f[7] %in% c("PASS", ".", "")) return(FALSE)
    if (!is.null(min_qual)) {
      q <- suppressWarnings(as.numeric(f[6]))
      if (!is.na(q) && q < min_qual) return(FALSE)
    }
    TRUE
  }

  repeat {
    lines <- readLines(con_in, n = chunk_lines, warn = FALSE)
    if (length(lines) == 0L) break

    # Header lines (always at the top of a VCF) are copied verbatim.
    if (!header_done) {
      is_hdr <- startsWith(lines, "#")
      hdr    <- lines[is_hdr]
      if (length(hdr) > 0L) writeLines(hdr, con_out)
      lines  <- lines[!is_hdr]
      if (length(lines) > 0L) header_done <- TRUE
      if (length(lines) == 0L) next
    }

    n_in <- n_in + length(lines)
    keep <- vapply(lines, .keep_line, logical(1), USE.NAMES = FALSE)
    kept <- lines[keep]
    if (length(kept) == 0L) next

    if (is.null(max_variants)) {
      writeLines(kept, con_out)
      n_out <- n_out + length(kept)
    } else {
      # Uniform reservoir sampling across the whole stream.
      for (ln in kept) {
        res_count <- res_count + 1L
        if (res_count <= max_variants) {
          reservoir[res_count] <- ln
        } else {
          r <- sample.int(res_count, 1L)
          if (r <= max_variants) reservoir[r] <- ln
        }
      }
    }
  }

  if (!is.null(max_variants)) {
    n_keep    <- min(res_count, max_variants)
    if (n_keep > 0L) writeLines(reservoir[seq_len(n_keep)], con_out)
    n_out <- n_keep
  }

  if (verbose)
    message(sprintf("[prefilter] %d records scanned -> %d written to %s",
                    n_in, n_out, vcf_out))
  invisible(list(out = vcf_out, n_in = n_in, n_out = n_out))
}
