# GenoSim 1.2.0

## New feature: founder-referenced (cumulative) inbreeding statistics

The within-generation `inbreeding_fis` (Wright's *F*~IS~) is referenced to each
generation's own allele frequencies, so a single round of random mating restores
Hardy–Weinberg proportions and it returns to ~0 — it cannot represent the
cumulative inbreeding that accrues over generations, and small founder samples
also drive it spuriously negative. `simulate_population()` and
`simulate_from_pedigree()` now additionally report, per generation, founder-
referenced statistics that accumulate as drift erodes diversity:

* `fis_unbiased` — *F*~IS~ using the unbiased Nei gene-diversity estimator
  (`2n/(2n−1)` correction), removing the small-sample negativity.
* `fst_vs_founder` — *F*~ST~ = 1 − H~e~(t)/H~e~(0), the cumulative loss of gene
  diversity since the founder generation (the "evolutionary" inbreeding).
* `fit_vs_founder` — *F*~IT~ = 1 − H~o~(t)/H~e~(0), total individual inbreeding
  relative to the founders, satisfying 1 − *F*~IT~ = (1 − *F*~IS~)(1 − *F*~ST~).
* `expected_fst_drift` — the theoretical Wright–Fisher expectation
  1 − (1 − 1/2*N*~e~)^t for validation.
* `ne_estimate` — realised effective size from the heterozygosity decay.
* `mean_pedigree_F` — mean kinship-based pedigree inbreeding of the genotyped
  individuals in each observed generation (the input *F* in population mode).

All quantities are computed over loci polymorphic in the founders using
ratio-of-sums (Nei's *G*~ST~ form) and are missing-data aware. The derivation is
provided as a notebook in the project repository.

# GenoSim 1.1.3

## Correctness fixes (mathematical / logical)

* `compute_ld()` now reports the correct linkage-disequilibrium statistics.
  Previously `r^2` was under-estimated by a factor of ~4 and `D'` by a factor
  of ~2 (two perfectly correlated SNPs returned `r^2 = 0.25`). `r^2` is now the
  squared correlation of allele dosages and the haplotypic `D` is recovered
  correctly from the dosage covariance (`cov(x,y) = 2D`). New `min_pair_n`
  argument controls the minimum shared non-missing individuals per pair.

* `hwe_test()` no longer biases the chi-squared statistic when genotypes are
  missing: expected counts now use the number of *called* genotypes rather than
  the total number of individuals (previously 50% missingness could turn a true
  HWE locus into a spurious deviation).

* **Selection and mutation now actually affect the simulated genotypes.** In
  `simulate_population()` (and the synthetic generations of
  `simulate_from_pedigree()`) `selection_s` and `mut_rate` previously modified
  only a recorded allele-frequency vector, leaving the genotype output
  unchanged. Mutation is now applied to each transmitted gamete and selection
  as fitness-proportional resampling of offspring (preserving linkage), so both
  forces are reflected in the genotypes and all downstream analyses.

* `allele_freqs` and `summary_stats` from `simulate_population()` are now
  computed from the genotypes themselves, so `exp_heterozygosity`,
  `inbreeding_fis`, `mean_maf` and `frac_fixed` are mutually consistent with the
  returned matrices (same class of fix previously applied to
  `simulate_from_pedigree()`). `n_eff` now caps the per-generation breeding pool.

* `simulate_population()` SNP identifiers now match their coordinates: the
  position embedded in `snp_id` (`chrX_<pos>`) is the same value stored in
  `pos_bp` (previously they were two independent random draws).

* `export_plink()` now assigns each individual's sex once and writes it
  identically to the `.ped` and `.raw` files (previously the two files were
  randomised independently), tolerates missing genotypes (written as PLINK
  `0 0` / `NA`), and uses allele letters consistent with `export_vcf()`.

* Removed the now-unused allele-frequency helpers and a remaining same-name
  duplicate (`.draw_population_afs`) to keep a single source of truth and avoid
  the file-load-order shadowing that caused earlier bugs.

# GenoSim 1.1.2

## Bug fixes and robustness

* `simulate_from_pedigree()` no longer propagates `NA` allele frequencies (and
  hence `NA` genotypes) into synthetic generations. The Wright–Fisher drift
  helper now clips allele frequencies strictly inside `(0, 1)` before calling
  `rbinom()` (floating-point values landing on the `[0, 1]` boundary made
  `rbinom()` return `NA`), with a defensive fallback if `rbinom()` still
  returns `NA`. Synthetic-generation allele frequencies are also computed
  NA-safely.

* `detect_roh()` and `compute_ld()` are now NA-aware and no longer error or
  silently emit `NA` results on incomplete matrices:
  - `detect_roh()` treats missing calls as run-breaking (non-homozygous) calls
    instead of erroring in `if (!is_hom[k])`.
  - `compute_ld()` estimates each SNP pair from individuals genotyped at both
    loci (pairwise-complete), with a new `min_pair_n` argument.
  Both report any missingness via a message.

* `read_pedigree()` now coerces blank `sex` / `phenotype` cells to their
  documented missing codes (`sex = 0`, `phenotype = -9`) at load time, so
  `plot_pedigree_tree()` and others never see `NA` in these fields. The tree
  and ROH plots were additionally made NA-safe defensively.

* `read_pedigree()` now validates `individual_id`: duplicate or blank IDs raise
  an informative error at load time instead of the previous opaque
  `'length = 2' in coercion to 'logical(1)'` failure deep inside the simulator.

* Genotype imputation (`impute_missing = TRUE` in `read_vcf_cohort()`) now
  guarantees a complete matrix: sites with no observed calls fall back to the
  global allele frequency, any residual missing values are filled, and the
  result is returned as an integer matrix.

## New features

* `read_vcf_cohort()` gains `max_missing_rate`: sites whose per-SNP missingness
  exceeds this fraction are dropped before imputation. This is the recommended
  way to handle raw WES/WGS callsets (default `1` preserves prior behaviour).
  The function also warns when the input is too large for an in-memory pure-R
  parse and documents its expected input characteristics and size limits.

* New `prefilter_vcf()` streams a large (optionally gzipped) VCF in chunks and
  writes a filtered subset (PASS / biallelic SNP / autosomal / QUAL, with
  optional representative reservoir subsampling) without loading the whole file
  into memory. It can optionally delegate to `bcftools` when available. Use it
  to reduce multi-million-variant WES/WGS VCFs to an input `read_vcf_cohort()`
  can handle.
