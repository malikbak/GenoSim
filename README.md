# GenoSim <img src="man/figures/logo.png" align="right" height="100" alt="" />

<!-- badges -->
![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Tests](https://img.shields.io/badge/tests-73%20passing-brightgreen)
![Status](https://img.shields.io/badge/status-stable-green)

**GenoSim** is an R package for forward-time simulation of diploid SNP
genotype data in clinical and population genetics settings where real NGS
human data are scarce due to ethical, budgetary, or regulatory constraints.

It provides two simulation modes:

| Mode | Function | Use case |
|---|---|---|
| **Population** | `simulate_population()` | Synthetic cohorts, power studies, pipeline testing |
| **Pedigree** | `simulate_from_pedigree()` | Family-based studies using real VCF data |

---

## Features

- Hardy-Weinberg equilibrium baseline with inbreeding coefficient F
- Wright-Fisher genetic drift (finite effective population size N_e)
- Additive directional selection coefficient *s*
- Per-locus per-generation mutation (bidirectional)
- Recombination-aware gamete transmission (Haldane mapping)
- Pedigree-constrained simulation from real family VCF files
- Per-individual inbreeding F computed by recursive kinship
- Up to **10 generations**, up to **100,000 SNPs**, all 22 autosomes
- Analyses: HWE test, LD (r^2 and D'), ROH detection, PCA, FST, Nei diversity
- Export: VCFv4.2, PLINK PED/MAP/RAW, tidy CSV
- 10 base-R plots -- no ggplot2 or Bioconductor dependencies

---

## Installation

### From a local zip (recommended for offline / clinical environments)

```r
# Download GenoSim_1.0.0.tar.gz then:
install.packages("path/to/GenoSim_1.0.0.tar.gz", repos = NULL, type = "source")
```

### From GitHub (requires remotes or devtools)

```r
# install.packages("remotes")
remotes::install_github("example/GenoSim")
```

### Dependencies

GenoSim uses **base R only** -- `stats`, `utils`, `graphics`, `grDevices`, `tools`.
No Bioconductor, no tidyverse, no external C libraries.

Optional suggested packages (only needed if you use them explicitly):
- `testthat >= 3.0.0` -- running the test suite
- `knitr`, `rmarkdown` -- building the vignette

---

## Quick Start

### Population mode (no real data needed)

```r
library(GenoSim)

# Simulate 100 founders, 1000 SNPs, 5 generations
# with first-cousin inbreeding (F = 0.125)
sim <- simulate_population(
  n_founders    = 100,
  n_snps        = 1000,
  n_generations = 5,
  inbreeding_F  = 0.125,
  n_eff         = 80,
  chromosomes   = 1:22,
  seed          = 42
)

# Summary statistics per generation
sim$summary_stats
#>   generation n_individuals n_snps obs_heterozygosity exp_heterozygosity inbreeding_fis mean_maf frac_fixed
#>           0           100   1000             0.26543            0.30451         0.12840  0.24781          0
#>           1           150   1000             0.25891            0.29820         0.13177  0.24263          0
#>  ...

# Six-panel dashboard (writes to PDF or plots inline)
all_g <- do.call(rbind, sim$genotypes)
pca   <- run_pca(all_g)
plot_dashboard(sim, pca_result = pca, out_file = "dashboard.pdf")

# Export to VCF
export_vcf(sim, generation = c(0, 5), out_dir = "output/")

# Export to PLINK
export_plink(sim, generation = 5, out_prefix = "cohort", out_dir = "output/")

# Export all generations to CSV
export_csv(sim, out_dir = "output/csv/")
```

### Pedigree mode (with real family VCF data)

```r
library(GenoSim)

# 1. Read your family VCFs (one per individual, or one multi-sample VCF)
vcf <- read_vcf_cohort("path/to/family_vcfs/")

# 2. Read and validate the pedigree
ped <- read_pedigree("path/to/family_pedigree.csv")
summarise_pedigree(ped)

# 3. Simulate: propagate through real family structure
#    then append 4 extra synthetic generations
sim <- simulate_from_pedigree(
  vcf_cohort        = vcf,
  pedigree          = ped,
  extra_generations = 4,
  mut_rate          = 1e-4,
  seed              = 42
)

# 4. Analyse
hwe    <- hwe_test(sim$genotypes[["gen0"]])
roh    <- detect_roh(do.call(rbind, sim$genotypes[1:4]), sim$snp_map)
fst    <- compute_fst(sim$genotypes)

# 5. Visualise
plot_pedigree_tree(ped, highlight_ids = ped$individual_id[ped$phenotype == 2])
plot_af_comparison(vcf, sim)
plot_kinship_heatmap(ped)

# 6. Export
export_vcf(sim, out_dir = "output/")
```

### Using bundled example data

```r
library(GenoSim)

# FAM_KHAN: 25 individuals, 4 generations, 200 SNPs, 4 affected
ped <- read_pedigree(example_ped_path())
vcf <- read_vcf_cohort(example_vcf_dir())

sim <- simulate_from_pedigree(vcf, ped, extra_generations = 3, seed = 1)
sim$summary_stats
```

---

## Pedigree CSV Format

Your pedigree CSV must contain at minimum three columns.
Column names are matched case-insensitively with common aliases.

| Column | Required | Values |
|---|---|---|
| `individual_id` | **Yes** | Unique ID -- must match VCF sample names |
| `father_id` | **Yes** | Father's `individual_id`, or `"0"` if unknown |
| `mother_id` | **Yes** | Mother's `individual_id`, or `"0"` if unknown |
| `sex` | No | `1` = male, `2` = female, `0` = unknown |
| `phenotype` | No | `1` = unaffected, `2` = affected, `-9` = missing |
| `family_id` | No | Family label, e.g. `"FAM001"` |
| `generation` | No | Integer depth (auto-inferred if absent) |

**Accepted aliases:** `iid`/`id`/`sample_id` for individual; `fid`/`father`/`pat` for father;
`mid`/`mother`/`mat` for mother; `gender` for sex; `pheno`/`affection` for phenotype.

Example (`family_pedigree.csv`):

```
individual_id,father_id,mother_id,sex,phenotype,family_id,generation
F001,0,0,1,1,FAM001,0
F002,0,0,2,1,FAM001,0
I101,F001,F002,1,1,FAM001,1
I102,F001,F002,2,2,FAM001,1
I201,I101,I102,1,2,FAM001,2
```

---

## VCF Format

Standard VCFv4.2. One sample column per file (per-individual VCFs) or multiple
sample columns (multi-sample VCF). The VCF sample name must match `individual_id`
in the pedigree.

```
##fileformat=VCFv4.2
##reference=GRCh38
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS       ID             REF  ALT  QUAL  FILTER  INFO  FORMAT  I101
1       2539584   chr1_2539584   A    T    .     PASS    .     GT      0/1
1       19261201  chr1_19261201  A    T    .     PASS    .     GT      0/0
```

Supported GT codes: `0/0`, `0/1`, `1/1`, `0|0`, `0|1`, `1|1`, `./.`
Indels and multi-allelic sites are automatically excluded.
Sex chromosomes (X, Y) and non-human contigs are silently skipped.

---

## All Functions

### Simulation
| Function | Description |
|---|---|
| `simulate_population()` | Forward-time population simulator |
| `simulate_from_pedigree()` | Pedigree-constrained simulator using real VCF data |

### Data input
| Function | Description |
|---|---|
| `read_vcf_cohort()` | Read VCF directory or multi-sample VCF |
| `read_pedigree()` | Read and validate pedigree CSV |
| `summarise_pedigree()` | Print pedigree structure summary |
| `extract_mating_pairs()` | Extract parent pairs per generation |

### Analysis
| Function | Description |
|---|---|
| `hwe_test()` | Chi-squared HWE test per SNP |
| `compute_ld()` | Pairwise r^2 and D' linkage disequilibrium |
| `detect_roh()` | Runs of homozygosity detection |
| `run_pca()` | PCA on dosage matrix |
| `compute_fst()` | Weir-Cockerham FST between generations |
| `diversity_metrics()` | Nei diversity, Watterson's theta |

### Export
| Function | Description |
|---|---|
| `export_vcf()` | Write VCFv4.2 per generation |
| `export_plink()` | Write PLINK .ped/.map/.raw |
| `export_csv()` | Write tidy CSV files |

### Visualisation
| Function | Description |
|---|---|
| `plot_dashboard()` | Six-panel summary dashboard |
| `plot_af_trajectory()` | Allele frequency drift across generations |
| `plot_heterozygosity()` | Observed vs expected heterozygosity |
| `plot_maf_distribution()` | MAF histograms per generation |
| `plot_fis_trajectory()` | Inbreeding Fis across generations |
| `plot_pca()` | PCA scatter coloured by generation |
| `plot_ld_decay()` | LD r^2 decay with distance |
| `plot_pedigree_tree()` | Family pedigree diagram |
| `plot_kinship_heatmap()` | Pairwise kinship heatmap |
| `plot_roh_per_individual()` | Total ROH per individual |
| `plot_af_comparison()` | Observed vs simulated MAF density |

### Example data
| Function | Description |
|---|---|
| `example_ped_path()` | Path to bundled FAM_KHAN pedigree CSV |
| `example_vcf_dir()` | Path to bundled FAM_KHAN VCF directory |
| `example_pedigree` | Pedigree data frame (via `data(example_pedigree)`) |

---

## Genetic Models

### Hardy-Weinberg + Inbreeding (F-model)

Genotype probabilities under inbreeding coefficient F:

```
P(AA) = p^2(1-F) + p?F
P(Aa) = 2pq(1-F)
P(aa) = q^2(1-F) + q?F
```

A value of F = 0 gives standard HWE. F = 0.125 approximates first-cousin
mating; F = 0.25 approximates half-sibling matings.

### Wright-Fisher Drift

Each generation, allele counts are drawn from:
`Binomial(2?N_e, p)`

Smaller N_e = stronger drift = faster fixation.

### Additive Selection

Fitness: w_AA = 1+s, w_Aa = 1+s/2, w_aa = 1.
- s > 0: positive selection on the alt allele
- s < 0: negative (purifying) selection
- s = 0: neutral

### Recombination

Haldane mapping function applied between adjacent SNPs on each chromosome,
with default rate 10?? per bp per meiosis:

`r = 0.5 * (1 - exp(-2 * 10?? * d_bp))`

---

## Test Suite

```r
library(testthat)
library(GenoSim)
test_dir(system.file("tests/testthat", package = "GenoSim"),
         package = "GenoSim")
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 73 ]
```

73 tests covering:
- Population simulation correctness and reproducibility
- Parameter validation
- HWE, LD, ROH, PCA, FST, diversity metrics
- Pedigree parsing, validation, inbreeding computation
- VCF reading, merging, imputation
- End-to-end pedigree simulation
- VCF, PLINK, and CSV export formats

---

---

## License

MIT -- see [LICENSE.md](LICENSE.md)
