# =============================================================================
#  GenoSim R package --- export.R
#  VCF, PLINK PED/MAP/RAW, tidy CSV writers
# =============================================================================

#' Export Simulated Genotypes to VCF Format
#'
#' @description
#' Writes one VCFv4.2 file per requested generation. Each file contains all
#' simulated individuals from that generation as sample columns, with
#' standard fixed fields and \code{GT} format.
#'
#' @param sim_result A list returned by \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#' @param generation Integer vector or \code{"all"}. Generation indices
#'   (0-based) to export. Default \code{"all"}.
#' @param out_dir Character. Output directory path (created if absent).
#'   Default \code{"."}.
#' @param compress Logical. Attempt gzip compression via \pkg{R.utils} if
#'   available. Default \code{FALSE}.
#'
#' @return Invisibly, a character vector of written file paths.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_population(n_founders=50, n_snps=100, n_generations=2, seed=1)
#' export_vcf(sim, generation = c(0,2), out_dir = "output/")
#' }
#'
#' @seealso \code{\link{export_plink}}, \code{\link{export_csv}}
#' @export
export_vcf <- function(sim_result, generation = "all",
                        out_dir = ".", compress = FALSE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  gens <- if (identical(generation, "all"))
    seq_along(sim_result$genotypes) - 1L
  else as.integer(generation)

  snp_map <- sim_result$snp_map
  written <- character()

  for (gen in gens) {
    geno    <- sim_result$genotypes[[gen + 1L]]
    if (is.null(geno) || nrow(geno) == 0) next
    af_row  <- sim_result$allele_freqs[gen + 1L, ]
    sids    <- rownames(geno)
    fname   <- file.path(out_dir, sprintf("genosim_gen%02d.vcf", gen))
    con     <- file(fname, "w")

    writeLines(c(
      "##fileformat=VCFv4.2",
      paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
      "##source=GenoSim_v1.2.0",
      "##reference=GRCh38",
      '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
      '##INFO=<ID=GEN,Number=1,Type=Integer,Description="Simulated Generation">',
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sids),
            collapse = "\t")
    ), con = con)

    gt_map <- c("0"="0/0","1"="0/1","2"="1/1")
    for (j in seq_len(nrow(snp_map))) {
      gt_vec <- gt_map[as.character(geno[, j])]
      writeLines(paste(c(snp_map$chrom[j], snp_map$pos_bp[j], snp_map$snp_id[j],
                         snp_map$ref[j] %||% "A", snp_map$alt[j] %||% "T",
                         ".", "PASS",
                         sprintf("AF=%.4f;GEN=%d", af_row[j], gen),
                         "GT", gt_vec), collapse = "\t"), con = con)
    }
    close(con)
    if (compress) {
      gz_fname <- paste0(fname, '.gz')
      con_in  <- file(fname, 'rb')
      con_out <- gzcon(file(gz_fname, 'wb'))
      writeBin(readBin(con_in, 'raw', file.info(fname)$size), con_out)
      close(con_in); close(con_out)
      file.remove(fname)
      fname <- gz_fname
    }
    written <- c(written, fname)
    message(sprintf("  VCF written: %s (%d samples)", fname, nrow(geno)))
  }
  invisible(written)
}

#' Export Simulated Genotypes to PLINK Format
#'
#' @description
#' Writes PLINK-format files (\code{.ped}, \code{.map}, and \code{.raw}) for a
#' specified generation. The \code{.raw} file uses additive dosage coding (0/1/2)
#' and is compatible with \code{--linear} and \code{--logistic} in PLINK.
#'
#' @param sim_result A list returned by \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#' @param generation Integer. Generation index (0-based) to export. Default
#'   \code{0}.
#' @param out_prefix Character. Filename prefix for the output files.
#'   Default \code{"genosim"}.
#' @param out_dir Character. Output directory. Default \code{"."}.
#'
#' @return Invisibly, a list with paths \code{ped}, \code{map}, \code{raw}.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_population(n_founders=50, n_snps=100, n_generations=2, seed=1)
#' export_plink(sim, generation = 1, out_prefix = "myfamily")
#' }
#'
#' @seealso \code{\link{export_vcf}}, \code{\link{export_csv}}
#' @export
export_plink <- function(sim_result, generation = 0L,
                          out_prefix = "genosim", out_dir = ".") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  geno    <- sim_result$genotypes[[generation + 1L]]
  if (is.null(geno) || nrow(geno) == 0)
    stop("No genotypes for generation ", generation)
  snp_map <- sim_result$snp_map
  n_ind   <- nrow(geno)
  stub    <- file.path(out_dir, paste0(out_prefix, "_gen", sprintf("%02d", generation)))

  # Allele letters per SNP, consistent with export_vcf() (REF default "A",
  # ALT default "T"): dosage 0 -> REF/REF, 1 -> REF/ALT, 2 -> ALT/ALT.
  ref <- if (!is.null(snp_map$ref)) snp_map$ref else rep("A", nrow(snp_map))
  alt <- if (!is.null(snp_map$alt)) snp_map$alt else rep("T", nrow(snp_map))
  ref[is.na(ref) | ref == ""] <- "A"
  alt[is.na(alt) | alt == ""] <- "T"

  # Assign SEX once and reuse for BOTH .ped and .raw so an individual's sex is
  # identical across files. Use the pedigree's recorded sex when available.
  sex_vec <- rep(0L, n_ind)
  ped_tab <- sim_result$pedigree
  if (!is.null(ped_tab) && all(c("individual_id","sex") %in% names(ped_tab))) {
    m <- match(rownames(geno), ped_tab$individual_id)
    sex_vec <- ifelse(is.na(m), 0L, as.integer(ped_tab$sex[m]))
    sex_vec[is.na(sex_vec)] <- 0L
  }

  # .map
  map_df <- data.frame(CHR=snp_map$chrom, SNP=snp_map$snp_id,
                        CM=0, BP=snp_map$pos_bp)
  utils::write.table(map_df, paste0(stub, ".map"),
                     quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

  # .ped  (missing dosages -> PLINK missing genotype "0 0")
  al_mat <- matrix("0 0", nrow = n_ind, ncol = nrow(snp_map))
  for (j in seq_len(nrow(snp_map))) {
    gj <- geno[, j]
    al_mat[, j] <- ifelse(is.na(gj), "0 0",
                   ifelse(gj == 0L, paste(ref[j], ref[j]),
                   ifelse(gj == 1L, paste(ref[j], alt[j]),
                                    paste(alt[j], alt[j]))))
  }
  ped_df <- data.frame(
    FID=paste0("GEN",generation), IID=rownames(geno),
    PAT=0, MAT=0, SEX=sex_vec, PHENO=-9L, al_mat, stringsAsFactors=FALSE)
  utils::write.table(ped_df, paste0(stub, ".ped"),
                     quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

  # .raw  (additive dosage; missing dosages written as NA)
  raw_df <- cbind(
    data.frame(FID=paste0("GEN",generation), IID=rownames(geno),
               PAT=0, MAT=0, SEX=sex_vec, PHENO=-9L),
    as.data.frame(geno))
  utils::write.table(raw_df, paste0(stub, ".raw"),
                     quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t",
                     na = "NA")

  message(sprintf("  PLINK: %s (.ped/.map/.raw) | %d individuals", stub, n_ind))
  invisible(list(ped=paste0(stub,".ped"), map=paste0(stub,".map"),
                 raw=paste0(stub,".raw")))
}

#' Export All Simulation Results to Tidy CSV Files
#'
#' @description
#' Writes four CSV files: a wide genotype table (all generations stacked),
#' the per-generation summary statistics, the allele frequency trajectory,
#' and the SNP map.
#'
#' @param sim_result A list returned by \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#' @param out_dir Character. Output directory. Default \code{"."}.
#' @param include_snp_map Logical. Whether to write the SNP map CSV. Default
#'   \code{TRUE}.
#'
#' @return Invisibly, a named list of written file paths.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_population(n_founders=50, n_snps=100, n_generations=2, seed=1)
#' export_csv(sim, out_dir = "output/csv/")
#' }
#'
#' @seealso \code{\link{export_vcf}}, \code{\link{export_plink}}
#' @export
export_csv <- function(sim_result, out_dir = ".", include_snp_map = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  geno_long <- do.call(rbind, lapply(seq_along(sim_result$genotypes), function(i) {
    g  <- sim_result$genotypes[[i]]
    df <- as.data.frame(g)
    df$individual_id <- rownames(g)
    df$generation    <- i - 1L
    df[, c("individual_id","generation", setdiff(names(df), c("individual_id","generation")))]
  }))

  geno_path  <- file.path(out_dir, "genosim_genotypes_all_generations.csv")
  stats_path <- file.path(out_dir, "genosim_summary_stats.csv")
  af_path    <- file.path(out_dir, "genosim_allele_freqs.csv")
  map_path   <- file.path(out_dir, "genosim_snp_map.csv")

  utils::write.csv(geno_long, geno_path, row.names = FALSE)
  utils::write.csv(sim_result$summary_stats, stats_path, row.names = FALSE)
  af_df <- as.data.frame(sim_result$allele_freqs)
  af_df$generation <- 0:(nrow(af_df)-1L)
  utils::write.csv(af_df, af_path, row.names = FALSE)
  if (include_snp_map) utils::write.csv(sim_result$snp_map, map_path, row.names=FALSE)

  message(sprintf("  CSV export complete -?? %s", out_dir))
  invisible(list(genotypes=geno_path, stats=stats_path,
                 allele_freq=af_path, snp_map=map_path))
}
