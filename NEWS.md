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
