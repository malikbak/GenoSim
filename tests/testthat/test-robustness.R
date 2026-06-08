# Regression tests for the v1.1.2 robustness fixes.

# ---- Issue 1: NA-aware ROH / LD ---------------------------------------------

test_that("detect_roh tolerates NA genotypes without erroring", {
  sim <- simulate_population(n_founders = 30, n_snps = 200, n_generations = 2,
                             inbreeding_F = 0.25, chromosomes = 1:2,
                             seed = 11, verbose = FALSE)
  g <- sim$genotypes[[3]]
  g[sample(length(g), 50)] <- NA          # inject missingness
  expect_message(roh <- detect_roh(g, sim$snp_map, min_snps = 5,
                                   min_length_bp = 1e5),
                 "missing genotype")
  expect_type(roh, "list")
  expect_named(roh, c("roh_segments", "roh_per_individual"))
})

test_that("compute_ld tolerates NA genotypes (pairwise complete)", {
  sim <- simulate_population(n_founders = 60, n_snps = 100, n_generations = 2,
                             chromosomes = 1:2, seed = 12, verbose = FALSE)
  g <- sim$genotypes[[1]]
  g[sample(length(g), 100)] <- NA
  expect_message(ld <- compute_ld(g, sim$snp_map, max_snps = 80,
                                  max_dist_bp = 50e6),
                 "missing genotype")
  if (nrow(ld) > 0) {
    expect_false(any(is.na(ld$r2)))
    expect_true(all(ld$r2 >= 0 & ld$r2 <= 1))
  }
})

test_that("simulate_from_pedigree synthetic generations contain no NA", {
  vcf <- read_vcf_cohort(example_vcf_dir(), verbose = FALSE)
  ped <- read_pedigree(example_ped_path(), verbose = FALSE)
  sim <- simulate_from_pedigree(vcf, ped, extra_generations = 3,
                                seed = 7, verbose = FALSE)
  synth <- which(sim$summary_stats$source == "synthetic")
  for (i in synth) {
    g <- sim$genotypes[[i]]
    if (!is.null(g) && nrow(g) > 0) expect_false(anyNA(g))
  }
})

# ---- summary_stats consistency with stored genotypes ------------------------

test_that("simulate_from_pedigree summary_stats match the stored genotypes", {
  vcf <- read_vcf_cohort(example_vcf_dir(), verbose = FALSE)
  ped <- read_pedigree(example_ped_path(), verbose = FALSE)
  # Non-trivial selection/mutation make the evolved AF diverge from the
  # genotype AF, which previously corrupted exp_het/FIS/MAF/frac_fixed for
  # synthetic generations.
  sim <- suppressMessages(simulate_from_pedigree(
    vcf, ped, extra_generations = 3, selection_s = 0.1, mut_rate = 1e-3,
    seed = 42, verbose = FALSE))
  ss <- sim$summary_stats
  for (i in seq_along(sim$genotypes)) {
    g <- sim$genotypes[[i]]
    if (is.null(g) || nrow(g) == 0) next
    af <- colMeans(g, na.rm = TRUE) / 2
    ho <- mean(g == 1L, na.rm = TRUE)
    he <- mean(2 * af * (1 - af), na.rm = TRUE)
    expect_equal(ss$obs_heterozygosity[i], round(ho, 5), tolerance = 1e-5)
    expect_equal(ss$exp_heterozygosity[i], round(he, 5), tolerance = 1e-5)
    expect_equal(ss$inbreeding_fis[i], round(1 - ho/he, 5), tolerance = 1e-5)
    expect_equal(ss$mean_maf[i], round(mean(pmin(af, 1 - af)), 5),
                 tolerance = 1e-5)
    expect_equal(ss$mean_dosage[i], round(mean(g, na.rm = TRUE), 5),
                 tolerance = 1e-5)
  }
})

# ---- Issue 2: missing phenotype/sex codes -----------------------------------

test_that("read_pedigree coerces blank sex/phenotype to documented codes", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(c("individual_id,father_id,mother_id,sex,phenotype",
               "A,0,0,,",
               "B,0,0,1,2",
               "C,A,B,,"), tmp)
  ped <- read_pedigree(tmp, verbose = FALSE)
  expect_false(anyNA(ped$sex))
  expect_false(anyNA(ped$phenotype))
  expect_equal(ped$sex[ped$individual_id == "A"], 0L)
  expect_equal(ped$phenotype[ped$individual_id == "A"], -9L)
  # plot must not error on (now non-NA) phenotype/sex
  expect_invisible(plot_pedigree_tree(ped))
  unlink(tmp)
})

# ---- Issue 3: duplicate individual IDs --------------------------------------

test_that("read_pedigree rejects duplicate individual IDs", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(c("individual_id,father_id,mother_id",
               "A,0,0",
               "A,0,0",
               "B,A,0"), tmp)
  expect_error(read_pedigree(tmp, verbose = FALSE), "Duplicate individual_id")
  unlink(tmp)
})

# ---- Issue 5: imputation completes the matrix -------------------------------

test_that("read_vcf_cohort with imputation returns a complete matrix", {
  vcf <- read_vcf_cohort(example_vcf_dir(), impute_missing = TRUE,
                         verbose = FALSE)
  expect_false(anyNA(vcf$geno_matrix))
})

test_that("max_missing_rate drops sparse sites before imputation", {
  vcf_full <- read_vcf_cohort(example_vcf_dir(), max_missing_rate = 1,
                              verbose = FALSE)
  vcf_strict <- read_vcf_cohort(example_vcf_dir(), max_missing_rate = 0.05,
                                verbose = FALSE)
  expect_lte(ncol(vcf_strict$geno_matrix), ncol(vcf_full$geno_matrix))
  expect_false(anyNA(vcf_strict$geno_matrix))
})

# ---- Issue 4: streaming prefilter -------------------------------------------

test_that("prefilter_vcf streams and filters a VCF without bcftools", {
  vcf_files <- list.files(example_vcf_dir(), pattern = "\\.vcf$",
                          full.names = TRUE)
  skip_if(length(vcf_files) == 0)
  out <- tempfile(fileext = ".vcf")
  res <- prefilter_vcf(vcf_files[1], out, verbose = FALSE)
  expect_true(file.exists(out))
  expect_equal(res$out, out)
  # Output must be re-readable as a (single-sample) cohort
  vcf <- read_vcf_cohort(out, verbose = FALSE)
  expect_true(ncol(vcf$geno_matrix) >= 1)
  unlink(out)
})

test_that("prefilter_vcf reservoir subsampling caps the record count", {
  vcf_files <- list.files(example_vcf_dir(), pattern = "\\.vcf$",
                          full.names = TRUE)
  skip_if(length(vcf_files) == 0)
  out <- tempfile(fileext = ".vcf")
  res <- prefilter_vcf(vcf_files[1], out, max_variants = 3, seed = 1,
                       verbose = FALSE)
  expect_lte(res$n_out, 3)
  unlink(out)
})
