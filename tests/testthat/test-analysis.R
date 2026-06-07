test_that("hwe_test returns correct structure", {
  sim <- simulate_population(n_founders=60, n_snps=100, n_generations=1,
                              chromosomes=1:3, seed=2, verbose=FALSE)
  res <- hwe_test(sim$genotypes[[1]], alpha=0.05)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("snp_id","chi2","p_value","hwe_pass","sig_label"))
  expect_equal(nrow(res), nrow(sim$snp_map))
  expect_true(all(res$sig_label %in% c("HWE_pass","HWE_deviation","fixed")))
})

test_that("hwe_test: inbreeding causes HWE deviations", {
  sim_n <- simulate_population(n_founders=80, n_snps=200, n_generations=1,
                                inbreeding_F=0,    chromosomes=1:5, seed=3, verbose=FALSE)
  sim_i <- simulate_population(n_founders=80, n_snps=200, n_generations=1,
                                inbreeding_F=0.35, chromosomes=1:5, seed=3, verbose=FALSE)
  pass_n <- mean(hwe_test(sim_n$genotypes[[1]])$hwe_pass, na.rm=TRUE)
  pass_i <- mean(hwe_test(sim_i$genotypes[[1]])$hwe_pass, na.rm=TRUE)
  expect_gt(pass_n, pass_i)
})

test_that("compute_ld returns valid r2 in [0,1]", {
  sim <- simulate_population(n_founders=60, n_snps=100, n_generations=2,
                              chromosomes=1:2, seed=5, verbose=FALSE)
  ld  <- compute_ld(sim$genotypes[[1]], sim$snp_map, max_snps=80, max_dist_bp=50e6)
  if (nrow(ld) > 0) {
    expect_true(all(ld$r2 >= 0 & ld$r2 <= 1, na.rm=TRUE))
    expect_true(all(abs(ld$D_prime) <= 1 + 1e-6, na.rm=TRUE))
  }
})

test_that("detect_roh returns correct structure", {
  sim <- simulate_population(n_founders=30, n_snps=300, n_generations=2,
                              inbreeding_F=0.25, chromosomes=1:3, seed=6, verbose=FALSE)
  roh <- detect_roh(sim$genotypes[[3]], sim$snp_map, min_snps=5, min_length_bp=1e5)
  expect_type(roh, "list")
  expect_named(roh, c("roh_segments","roh_per_individual"))
  segs <- roh$roh_segments
  expect_true(all(c("ind_id","chrom","length_bp","n_snps_roh") %in% names(segs)))
  if (nrow(segs) > 0) expect_true(all(segs$length_bp > 0))
})

test_that("run_pca returns explained variance summing to <= 100", {
  sim   <- simulate_population(n_founders=40, n_snps=200, n_generations=3,
                                chromosomes=1:3, seed=8, verbose=FALSE)
  all_g <- do.call(rbind, sim$genotypes)
  pca   <- run_pca(all_g, n_pc=5)
  expect_named(pca, c("scores","loadings","variance_pct","cumvar_pct","n_snps_used"))
  expect_true(all(pca$variance_pct >= 0))
  expect_lte(sum(pca$variance_pct), 100 + 1e-6)
  expect_equal(nrow(pca$scores), nrow(all_g))
})

test_that("compute_fst returns non-negative values", {
  sim <- simulate_population(n_founders=50, n_snps=200, n_generations=4,
                              chromosomes=1:3, seed=9, verbose=FALSE)
  fst <- compute_fst(sim$genotypes)
  expect_s3_class(fst, "data.frame")
  expect_true(all(fst$mean_fst >= 0))
  expect_equal(nrow(fst), length(sim$genotypes) - 1)
})

test_that("diversity_metrics returns one row per generation", {
  sim <- simulate_population(n_founders=50, n_snps=200, n_generations=3,
                              chromosomes=1:3, seed=10, verbose=FALSE)
  div <- diversity_metrics(sim)
  expect_equal(nrow(div), length(sim$genotypes))
  expect_true(all(div$nei_gene_div >= 0 & div$nei_gene_div <= 0.5))
})
