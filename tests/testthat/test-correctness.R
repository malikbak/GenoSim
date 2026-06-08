# Regression tests for the mathematical/logical correctness fixes (v1.1.3).

# ---- compute_ld: r^2 and D' scaling -----------------------------------------

test_that("compute_ld r2 equals dosage correlation^2 (perfect LD -> 1)", {
  set.seed(1)
  x  <- rbinom(500, 2, 0.4)
  g  <- cbind(s1 = x, s2 = x)                 # identical loci => r2 = 1
  sm <- data.frame(snp_id = c("s1","s2"), chrom = c(1,1), pos_bp = c(100,200))
  ld <- compute_ld(g, sm, max_dist_bp = 1e6)
  expect_equal(ld$r2, 1, tolerance = 1e-6)
  expect_equal(ld$D_prime, 1, tolerance = 1e-6)
})

test_that("compute_ld r2 matches cor(dosage)^2 for partial LD", {
  set.seed(2)
  a <- rbinom(800, 2, 0.3)
  b <- ifelse(rbinom(800, 1, 0.85) == 1, a, rbinom(800, 2, 0.3))
  g <- cbind(a = a, b = b)
  sm <- data.frame(snp_id = c("a","b"), chrom = c(1,1), pos_bp = c(1,2))
  ld <- compute_ld(g, sm, max_dist_bp = 1e6)
  expect_equal(ld$r2, cor(a, b)^2, tolerance = 1e-4)  # r2 is rounded to 5 dp
  expect_true(ld$r2 >= 0 && ld$r2 <= 1)
  expect_true(abs(ld$D_prime) <= 1 + 1e-9)
})

# ---- hwe_test: missing-data handling ----------------------------------------

test_that("hwe_test is not biased by missing genotypes", {
  set.seed(5)
  g <- matrix(sample(c(0L,1L,2L), 200, replace = TRUE, prob = c(.25,.5,.25)),
              ncol = 1)
  g_miss <- g; g_miss[1:100] <- NA           # half missing, same true freqs
  p_full <- hwe_test(g)$p_value
  p_miss <- hwe_test(g_miss)$p_value
  # In HWE -> both should be large; the old code returned ~0 for p_miss.
  expect_gt(p_full, 0.05)
  expect_gt(p_miss, 0.05)
})

# ---- selection & mutation actually change the genotypes ----------------------

test_that("selection_s shifts the genotype allele frequency (not just a table)", {
  neu <- simulate_population(n_founders = 150, n_snps = 400, n_generations = 6,
          selection_s = 0,   n_eff = 150, chromosomes = 1:4, seed = 7, verbose = FALSE)
  sel <- simulate_population(n_founders = 150, n_snps = 400, n_generations = 6,
          selection_s = 0.5, n_eff = 150, chromosomes = 1:4, seed = 7, verbose = FALSE)
  af_neu <- mean(colMeans(neu$genotypes[[7]]) / 2)
  af_sel <- mean(colMeans(sel$genotypes[[7]]) / 2)
  expect_gt(af_sel, af_neu)                  # positive selection raises alt AF
  expect_false(identical(neu$genotypes[[7]], sel$genotypes[[7]]))
})

test_that("mut_rate changes the genotypes", {
  m0 <- simulate_population(n_founders = 120, n_snps = 200, n_generations = 5,
          mut_rate = 0,    n_eff = 120, chromosomes = 1:2, seed = 3, verbose = FALSE)
  m1 <- simulate_population(n_founders = 120, n_snps = 200, n_generations = 5,
          mut_rate = 0.01, n_eff = 120, chromosomes = 1:2, seed = 3, verbose = FALSE)
  expect_false(identical(m0$genotypes[[6]], m1$genotypes[[6]]))
})

test_that("allele_freqs and summary_stats are consistent with genotypes", {
  sim <- simulate_population(n_founders = 100, n_snps = 300, n_generations = 4,
          selection_s = 0.2, mut_rate = 1e-3, n_eff = 100, chromosomes = 1:3,
          seed = 11, verbose = FALSE)
  for (i in seq_along(sim$genotypes)) {
    g  <- sim$genotypes[[i]]
    af <- colMeans(g) / 2
    expect_equal(unname(sim$allele_freqs[i, ]), unname(af), tolerance = 1e-9)
    ho <- mean(g == 1L); he <- mean(2 * af * (1 - af))
    expect_equal(sim$summary_stats$obs_heterozygosity[i], round(ho, 5), tolerance = 1e-5)
    expect_equal(sim$summary_stats$exp_heterozygosity[i], round(he, 5), tolerance = 1e-5)
    expect_equal(sim$summary_stats$inbreeding_fis[i], round(1 - ho/he, 5), tolerance = 1e-5)
  }
})

# ---- snp_id label matches pos_bp --------------------------------------------

test_that("snp_id-embedded position equals pos_bp", {
  sim <- simulate_population(n_founders = 20, n_snps = 80, n_generations = 1,
                             chromosomes = 1:4, seed = 9, verbose = FALSE)
  emb <- as.integer(sub(".*_", "", sim$snp_map$snp_id))
  expect_equal(emb, sim$snp_map$pos_bp)
})

# ---- export_plink: consistent SEX, NA-safe, VCF-consistent alleles ----------

test_that("export_plink writes identical SEX in .ped and .raw and tolerates NA", {
  sim <- simulate_population(n_founders = 20, n_snps = 20, n_generations = 1,
                             chromosomes = 1, seed = 1, verbose = FALSE)
  sim$genotypes[[1]][1, 1] <- NA              # inject a missing call
  od <- tempfile()
  p  <- suppressMessages(export_plink(sim, generation = 0, out_dir = od))
  ped_sex <- as.integer(read.table(p$ped)[, 5])
  raw_sex <- as.integer(read.table(p$raw, header = TRUE)[, 5])
  expect_identical(ped_sex, raw_sex)
  # First individual's first SNP was NA -> PLINK missing "0 0" (tokens 7 & 8).
  ped_line1 <- scan(p$ped, what = "character", nlines = 1, quiet = TRUE)
  expect_identical(ped_line1[7:8], c("0", "0"))
  unlink(od, recursive = TRUE)
})
