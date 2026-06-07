test_that("simulate_population returns correct structure", {
  sim <- simulate_population(
    n_founders    = 20,
    n_snps        = 50,
    n_generations = 3,
    chromosomes   = 1:3,
    seed          = 1,
    verbose       = FALSE
  )

  expect_type(sim, "list")
  expect_named(sim, c("genotypes","snp_map","allele_freqs","summary_stats","params"))

  # Correct number of generations stored
  expect_length(sim$genotypes, 4)  # gen0 .. gen3

  # Genotype matrices are integer with dosage 0/1/2
  g0 <- sim$genotypes[[1]]
  expect_true(is.matrix(g0))
  expect_true(all(g0 %in% 0:2))

  # SNP map has required columns
  expect_true(all(c("snp_id","chrom","pos_bp","founder_maf") %in% names(sim$snp_map)))

  # Allele freq matrix dimensions
  expect_equal(nrow(sim$allele_freqs), 4)
  expect_equal(ncol(sim$allele_freqs), nrow(sim$snp_map))

  # Summary stats has one row per generation
  expect_equal(nrow(sim$summary_stats), 4)
})

test_that("simulate_population respects seed reproducibility", {
  s1 <- simulate_population(n_founders=20, n_snps=30, n_generations=2,
                             chromosomes=1, seed=99, verbose=FALSE)
  s2 <- simulate_population(n_founders=20, n_snps=30, n_generations=2,
                             chromosomes=1, seed=99, verbose=FALSE)
  expect_identical(s1$genotypes[[1]], s2$genotypes[[1]])
})

test_that("inbreeding_F reduces observed heterozygosity", {
  sim_neutral <- simulate_population(n_founders=50, n_snps=200, n_generations=1,
                                      inbreeding_F=0,     chromosomes=1:5, seed=7, verbose=FALSE)
  sim_inbred  <- simulate_population(n_founders=50, n_snps=200, n_generations=1,
                                      inbreeding_F=0.25,  chromosomes=1:5, seed=7, verbose=FALSE)
  het_neutral <- sim_neutral$summary_stats$obs_heterozygosity[1]
  het_inbred  <- sim_inbred$summary_stats$obs_heterozygosity[1]
  expect_lt(het_inbred, het_neutral)
})

test_that("validate_params catches bad inputs", {
  expect_error(simulate_population(n_founders=5,  n_snps=100, n_generations=2, verbose=FALSE))
  expect_error(simulate_population(n_founders=20, n_snps=100, n_generations=11, verbose=FALSE))
  expect_error(simulate_population(n_founders=20, n_snps=100, n_generations=2,
                                    inbreeding_F=1.0, verbose=FALSE))
  expect_error(simulate_population(n_founders=20, n_snps=100, n_generations=2,
                                    chromosomes=c(23), verbose=FALSE))
})

test_that("selection_s shifts allele frequencies in expected direction", {
  sim_pos <- simulate_population(n_founders=100, n_snps=500, n_generations=5,
                                  selection_s=0.1, n_eff=500, chromosomes=1:5,
                                  seed=42, verbose=FALSE)
  sim_neu <- simulate_population(n_founders=100, n_snps=500, n_generations=5,
                                  selection_s=0.0, n_eff=500, chromosomes=1:5,
                                  seed=42, verbose=FALSE)
  # Positive selection should increase mean allele frequency relative to neutral
  p_pos <- mean(sim_pos$allele_freqs[nrow(sim_pos$allele_freqs), ])
  p_neu <- mean(sim_neu$allele_freqs[nrow(sim_neu$allele_freqs), ])
  expect_gt(p_pos, p_neu)
})
