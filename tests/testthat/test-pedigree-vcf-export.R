test_that("read_pedigree loads and validates FAM_KHAN example", {
  ped_path <- example_ped_path()
  expect_true(file.exists(ped_path))
  ped <- read_pedigree(ped_path, verbose=FALSE)
  expect_s3_class(ped, "data.frame")
  expect_equal(nrow(ped), 25)
  expect_true(all(c("individual_id","father_id","mother_id","sex","phenotype",
                    "family_id","generation","is_founder","inbreeding_F") %in% names(ped)))
  expect_equal(sum(ped$is_founder), 4)
  expect_equal(max(ped$generation), 3)
  expect_true(all(ped$inbreeding_F >= 0))
  expect_true(any(ped$phenotype == 2))
})

test_that("read_pedigree detects self-reference", {
  tmp <- tempfile(fileext=".csv")
  df  <- data.frame(individual_id="A", father_id="A", mother_id="0",
                    stringsAsFactors=FALSE)
  write.csv(df, tmp, row.names=FALSE)
  expect_error(read_pedigree(tmp, verbose=FALSE), "Self-reference")
  unlink(tmp)
})

test_that("summarise_pedigree runs without error", {
  ped <- read_pedigree(example_ped_path(), verbose=FALSE)
  expect_invisible(summarise_pedigree(ped))
})

test_that("extract_mating_pairs returns expected pairs for gen 1", {
  ped   <- read_pedigree(example_ped_path(), verbose=FALSE)
  pairs <- extract_mating_pairs(ped, generation=1)
  expect_s3_class(pairs, "data.frame")
  expect_true(nrow(pairs) >= 1)
  expect_true(all(c("father_id","mother_id") %in% names(pairs)))
})

test_that("read_vcf_cohort loads the example VCF directory", {
  vcf_dir <- example_vcf_dir()
  expect_true(dir.exists(vcf_dir))
  vcf <- read_vcf_cohort(vcf_dir, verbose=FALSE)
  expect_named(vcf, c("geno_matrix","snp_map"))
  expect_equal(nrow(vcf$geno_matrix), 25)
  expect_true(ncol(vcf$geno_matrix) > 0)
  expect_true(all(vcf$geno_matrix %in% 0:2, na.rm=TRUE))
  expect_true("cohort_maf" %in% names(vcf$snp_map))
})

test_that("simulate_from_pedigree runs end-to-end", {
  vcf <- read_vcf_cohort(example_vcf_dir(), verbose=FALSE)
  ped <- read_pedigree(example_ped_path(), verbose=FALSE)
  sim <- simulate_from_pedigree(vcf, ped,
                                 extra_generations=2,
                                 seed=42, verbose=FALSE)
  expect_named(sim, c("genotypes","pedigree","snp_map","allele_freqs",
                       "summary_stats","haplotypes","simulation_log","params"))
  # Should have gen0..gen5 (3 observed + 2 extra)
  expect_gte(length(sim$genotypes), 3)
  # Source column distinguishes observed from synthetic
  expect_true("source" %in% names(sim$summary_stats))
  expect_true(any(sim$summary_stats$source == "observed_pedigree"))
  expect_true(any(sim$summary_stats$source == "synthetic"))
  # Haplotype store should have at least the matched founders
  expect_gte(length(sim$haplotypes), 4)
})

test_that("export_csv creates expected files", {
  sim     <- simulate_population(n_founders=20, n_snps=30, n_generations=2,
                                  chromosomes=1, seed=1, verbose=FALSE)
  out_dir <- tempfile()
  paths   <- export_csv(sim, out_dir=out_dir)
  expect_true(file.exists(paths$genotypes))
  expect_true(file.exists(paths$stats))
  expect_true(file.exists(paths$allele_freq))
  geno_csv <- read.csv(paths$genotypes)
  expect_true("generation" %in% names(geno_csv))
  unlink(out_dir, recursive=TRUE)
})

test_that("export_vcf creates VCF files with correct header", {
  sim     <- simulate_population(n_founders=20, n_snps=30, n_generations=2,
                                  chromosomes=1, seed=2, verbose=FALSE)
  out_dir <- tempfile()
  files   <- export_vcf(sim, generation=0, out_dir=out_dir)
  expect_length(files, 1)
  expect_true(file.exists(files[1]))
  header_lines <- readLines(files[1], n=10)
  expect_true(any(grepl("fileformat=VCFv4.2", header_lines)))
  expect_true(any(grepl("^#CHROM", header_lines)))
  unlink(out_dir, recursive=TRUE)
})

test_that("export_plink creates .ped .map .raw files", {
  sim     <- simulate_population(n_founders=20, n_snps=30, n_generations=1,
                                  chromosomes=1, seed=3, verbose=FALSE)
  out_dir <- tempfile()
  paths   <- export_plink(sim, generation=0, out_dir=out_dir)
  expect_true(file.exists(paths$ped))
  expect_true(file.exists(paths$map))
  expect_true(file.exists(paths$raw))
  map_df  <- utils::read.table(paths$map)
  expect_equal(ncol(map_df), 4)
  unlink(out_dir, recursive=TRUE)
})
