## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 4.5,
  warning   = FALSE,
  message   = FALSE
)
library(GenoSim)

## ----install, eval=FALSE------------------------------------------------------
# # From a local zip file:
# install.packages("GenoSim_1.1.3.tar.gz", repos = NULL, type = "source")
# 
# # Or via remotes from GitHub (once published):
# # remotes::install_github("example/GenoSim")

## ----pop-basic----------------------------------------------------------------
sim <- simulate_population(
  n_founders    = 100,
  n_snps        = 1000,
  n_generations = 5,
  chromosomes   = 1:5,
  seed          = 42,
  verbose       = FALSE
)

# Each list element is a dosage matrix [n_indiv x n_SNPs]
names(sim$genotypes)
dim(sim$genotypes[["gen0"]])

## ----pop-stats----------------------------------------------------------------
sim$summary_stats[, c("generation", "n_individuals",
                       "obs_heterozygosity", "mean_maf", "inbreeding_fis")]

## ----pop-inbred---------------------------------------------------------------
sim_inbred <- simulate_population(
  n_founders    = 80,
  n_snps        = 1000,
  n_generations = 5,
  inbreeding_F  = 0.125,
  n_eff         = 60,
  chromosomes   = 1:5,
  seed          = 2024,
  verbose       = FALSE
)

sim_inbred$summary_stats[, c("generation","obs_heterozygosity",
                              "exp_heterozygosity","inbreeding_fis")]

## ----pop-sel------------------------------------------------------------------
sim_sel <- simulate_population(
  n_founders    = 200,
  n_snps        = 2000,
  n_generations = 8,
  selection_s   = 0.05,    # positive selection on alt allele
  n_eff         = 150,
  chromosomes   = 1:10,
  seed          = 999,
  verbose       = FALSE
)

## ----pop-viz, fig.height=5----------------------------------------------------
plot_dashboard(sim_inbred, pca_result = run_pca(
  do.call(rbind, sim_inbred$genotypes), n_pc = 5
))

## ----af-traj------------------------------------------------------------------
plot_af_trajectory(sim_inbred, n_snps_show = 30)

## ----hwe----------------------------------------------------------------------
hwe <- hwe_test(sim_inbred$genotypes[[1]], alpha = 0.05)
table(hwe$sig_label)

## ----ld-----------------------------------------------------------------------
ld <- compute_ld(
  sim$genotypes[[1]], sim$snp_map,
  max_snps = 200, max_dist_bp = 5e6
)
cat(sprintf("Pairs: %d | Mean r²: %.4f\n", nrow(ld), mean(ld$r2, na.rm=TRUE)))

## ----ld-plot------------------------------------------------------------------
plot_ld_decay(ld)

## ----roh----------------------------------------------------------------------
roh <- detect_roh(
  sim_inbred$genotypes[[6]],
  sim_inbred$snp_map,
  min_snps      = 10,
  min_length_bp = 500000
)
nrow(roh$roh_segments)

## ----roh-plot-----------------------------------------------------------------
plot_roh_per_individual(roh)

## ----pca----------------------------------------------------------------------
all_genos <- do.call(rbind, sim$genotypes)
pca <- run_pca(all_genos, n_pc = 10)
cat(sprintf("PC1: %.1f%%  PC2: %.1f%%  PC3: %.1f%%\n",
            pca$variance_pct[1], pca$variance_pct[2], pca$variance_pct[3]))

## ----pca-plot-----------------------------------------------------------------
gen_labels <- sub("_ind.*", "", rownames(pca$scores))
plot_pca(pca, color_by = gen_labels)

## ----fst----------------------------------------------------------------------
fst <- compute_fst(sim$genotypes)
print(fst)

## ----diversity----------------------------------------------------------------
div <- diversity_metrics(sim)
print(div)

## ----ped-load-----------------------------------------------------------------
ped <- read_pedigree(example_ped_path(), verbose = FALSE)
summarise_pedigree(ped)

## ----vcf-load-----------------------------------------------------------------
vcf <- read_vcf_cohort(example_vcf_dir(), verbose = FALSE)
cat(sprintf("Loaded: %d individuals × %d SNPs\n",
            nrow(vcf$geno_matrix), ncol(vcf$geno_matrix)))

## ----ped-sim------------------------------------------------------------------
ped_sim <- simulate_from_pedigree(
  vcf_cohort        = vcf,
  pedigree          = ped,
  extra_generations = 3,
  mut_rate          = 1e-4,
  seed              = 42,
  verbose           = FALSE
)

ped_sim$summary_stats[, c("generation","n_individuals",
                           "obs_heterozygosity","mean_maf","source")]

## ----ped-tree, fig.width=10, fig.height=7-------------------------------------
affected <- ped$individual_id[ped$phenotype == 2]
plot_pedigree_tree(ped, highlight_ids = affected,
                   title = "FAM_KHAN — Consanguineous Pedigree")

## ----kinship, fig.width=7, fig.height=6---------------------------------------
plot_kinship_heatmap(ped)

## ----af-compare---------------------------------------------------------------
plot_af_comparison(vcf, ped_sim)

## ----export-vcf, eval=FALSE---------------------------------------------------
# export_vcf(sim, generation = c(0, 5), out_dir = "output/vcf/")

## ----export-plink, eval=FALSE-------------------------------------------------
# export_plink(sim, generation = 5, out_prefix = "my_sim", out_dir = "output/plink/")

## ----export-csv, eval=FALSE---------------------------------------------------
# export_csv(sim, out_dir = "output/csv/")
# # Produces:
# #   genosim_genotypes_all_generations.csv
# #   genosim_summary_stats.csv
# #   genosim_allele_freqs.csv
# #   genosim_snp_map.csv

## ----session------------------------------------------------------------------
sessionInfo()

