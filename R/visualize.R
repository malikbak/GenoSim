# =============================================================================
#  GenoSim R package --- visualize.R
#  All plotting functions --- base R, no external dependencies
# =============================================================================

#' Plot Allele Frequency Trajectories Across Generations
#'
#' @description
#' Draws a line plot of alternative allele frequency over generations for a
#' random subsample of SNPs, overlaid with the population mean trajectory.
#'
#' @param sim_result List from \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#' @param n_snps_show Integer. Number of individual SNP trajectories to draw.
#'   Default \code{30}.
#' @param out_file Character or \code{NULL}. Path for PDF output. If
#'   \code{NULL} (default), the plot is drawn to the current device.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' sim <- simulate_population(n_founders=100, n_snps=200, n_generations=5, seed=1)
#' plot_af_trajectory(sim, n_snps_show=20)
#' @seealso \code{\link{plot_dashboard}}
#' @export
plot_af_trajectory <- function(sim_result, n_snps_show = 30L, out_file = NULL) {
  af   <- sim_result$allele_freqs
  gens <- 0:(nrow(af)-1L)
  idx  <- if (ncol(af) > n_snps_show) sample(ncol(af), n_snps_show) else seq_len(ncol(af))
  if (!is.null(out_file)) grDevices::pdf(out_file, width=9, height=5)
  graphics::par(mar=c(4,4.5,3,1), mgp=c(2.5,0.7,0))
  graphics::plot(NA, xlim=c(0, max(gens)), ylim=c(0,1),
                 xlab="Generation", ylab="Alt allele frequency",
                 main="Allele Frequency Trajectories", xaxt="n", las=1)
  graphics::axis(1, at=gens, labels=paste0("Gen ",gens), cex.axis=0.8)
  graphics::abline(h=c(0,0.5,1), lty=c(2,3,2), col="grey70", lwd=0.8)
  for (j in idx)
    graphics::lines(gens, af[,j], col=grDevices::adjustcolor(.SIM_COLORS[(j%%10)+1],0.4), lwd=0.9)
  graphics::lines(gens, rowMeans(af), col="#111827", lwd=2.5)
  graphics::legend("topright", legend=c("Individual SNPs","Mean AF"),
                   col=c(grDevices::adjustcolor(.SIM_COLORS[1],0.6),"#111827"),
                   lwd=c(1,2.5), bty="n", cex=0.85)
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}

#' Plot Observed vs Expected Heterozygosity Per Generation
#'
#' @description
#' Grouped bar chart comparing observed heterozygosity (deflated by inbreeding)
#' against expected HWE heterozygosity per generation, annotated with
#' \eqn{F_{IS}}.
#'
#' @param sim_result List from \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' sim <- simulate_population(n_founders=100, n_snps=500, n_generations=4,
#'                            inbreeding_F=0.125, seed=1)
#' plot_heterozygosity(sim)
#' @export
plot_heterozygosity <- function(sim_result, out_file = NULL) {
  ss <- sim_result$summary_stats
  if (!is.null(out_file)) grDevices::pdf(out_file, width=8, height=5)
  graphics::par(mar=c(4,4.5,3,1), mgp=c(2.5,0.7,0))
  bp <- graphics::barplot(rbind(ss$exp_heterozygosity, ss$obs_heterozygosity),
                          beside=TRUE, names.arg=paste0("G",ss$generation),
                          col=c("#BFDBFE","#2563EB"), border=NA,
                          ylim=c(0, max(ss$exp_heterozygosity,ss$obs_heterozygosity)*1.25),
                          ylab="Heterozygosity", main="Observed vs Expected Heterozygosity",
                          las=1, cex.names=0.85)
  graphics::legend("topright", legend=c("Expected (HWE)","Observed"),
                   fill=c("#BFDBFE","#2563EB"), border=NA, bty="n", cex=0.9)
  fis <- round(ss$inbreeding_fis, 3)
  graphics::text(colMeans(bp), ss$obs_heterozygosity + max(ss$obs_heterozygosity)*0.04,
                 labels=ifelse(is.na(fis),"",paste0("F=",fis)), cex=0.68, col="#374151")
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}

#' Plot MAF Distribution Per Generation
#'
#' @description
#' A grid of histograms showing the minor allele frequency distribution for
#' each generation, with mean MAF annotated.
#'
#' @param sim_result List from \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' sim <- simulate_population(n_founders=100, n_snps=500, n_generations=4, seed=1)
#' plot_maf_distribution(sim)
#' @export
plot_maf_distribution <- function(sim_result, out_file = NULL) {
  af    <- sim_result$allele_freqs
  n_gen <- nrow(af)
  ncols <- min(n_gen, 4L); nrows <- ceiling(n_gen/ncols)
  if (!is.null(out_file)) grDevices::pdf(out_file, width=ncols*3.5, height=nrows*3)
  graphics::par(mfrow=c(nrows,ncols), mar=c(3.5,3.5,2.5,1), mgp=c(2,0.6,0),
                oma=c(0,0,2,0))
  for (i in seq_len(n_gen)) {
    maf <- pmin(af[i,], 1-af[i,])
    maf <- maf[maf > 0 & maf < 1]
    graphics::hist(maf, breaks=30, col=.SIM_COLORS[((i-1)%%10)+1],
                   border="white", xlim=c(0,0.5),
                   main=paste0("Generation ", i-1L),
                   xlab="Minor allele frequency", ylab="Count", las=1, cex.main=0.95)
    graphics::abline(v=mean(maf), lty=2, col="#374151", lwd=1.5)
    graphics::legend("topright", legend=paste0("Mean=",round(mean(maf),3)),
                     bty="n", cex=0.75)
  }
  graphics::mtext("MAF Distribution by Generation", outer=TRUE, cex=1.1, font=2)
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}

#' Plot PCA of Simulated Genotypes
#'
#' @description
#' Scatter plot of PC1 vs PC2 from a \code{\link{run_pca}} result, with
#' optional colour grouping.
#'
#' @param pca_result List returned by \code{\link{run_pca}}.
#' @param color_by Character vector of length \eqn{n_{individuals}} for
#'   colour grouping (e.g. generation labels). Default \code{NULL}.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' sim <- simulate_population(n_founders=100, n_snps=300, n_generations=3, seed=1)
#' all_g <- do.call(rbind, sim$genotypes)
#' pca   <- run_pca(all_g)
#' gen_labels <- sub("_ind.*","", rownames(pca$scores))
#' plot_pca(pca, color_by = gen_labels)
#' @seealso \code{\link{run_pca}}
#' @export
plot_pca <- function(pca_result, color_by = NULL, out_file = NULL) {
  sc <- pca_result$scores; vp <- pca_result$variance_pct
  if (!is.null(out_file)) grDevices::pdf(out_file, width=7, height=6)
  graphics::par(mar=c(4.5,4.5,3,1), mgp=c(2.8,0.7,0))
  cols <- if (!is.null(color_by))
    .SIM_COLORS[(as.integer(as.factor(color_by))%%10)+1]
  else "#2563EB"
  graphics::plot(sc[,1], sc[,2],
                 xlab=sprintf("PC1 (%.1f%%)",vp[1]),
                 ylab=sprintf("PC2 (%.1f%%)",vp[2]),
                 main="PCA of Simulated Genotypes",
                 col=cols, pch=19, cex=0.7, las=1)
  if (!is.null(color_by)) {
    lvls <- levels(as.factor(color_by))
    graphics::legend("topright", legend=lvls,
                     col=.SIM_COLORS[seq_along(lvls)%%10+1],
                     pch=19, bty="n", cex=0.8)
  }
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}

#' Plot LD Decay
#'
#' @description
#' Scatter plot of per-pair \eqn{r^2} versus inter-SNP distance in kilobases,
#' overlaid with the binned mean \eqn{r^2} decay curve.
#'
#' @param ld_result Data frame returned by \code{\link{compute_ld}}.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' sim <- simulate_population(n_founders=100, n_snps=300, n_generations=3,
#'                            chromosomes=1:3, seed=1)
#' ld  <- compute_ld(sim$genotypes[[3]], sim$snp_map)
#' plot_ld_decay(ld)
#' @seealso \code{\link{compute_ld}}
#' @export
plot_ld_decay <- function(ld_result, out_file = NULL) {
  if (nrow(ld_result) == 0) { message("No LD data."); return(invisible(NULL)) }
  breaks  <- seq(0, max(ld_result$dist_bp), length.out=40)
  mids    <- (breaks[-1]+breaks[-length(breaks)])/2
  mean_r2 <- sapply(seq_len(length(breaks)-1), function(k) {
    sub <- ld_result$r2[ld_result$dist_bp>=breaks[k] & ld_result$dist_bp<breaks[k+1]]
    if (length(sub)==0) NA_real_ else mean(sub,na.rm=TRUE)
  })
  if (!is.null(out_file)) grDevices::pdf(out_file, width=8, height=5)
  graphics::par(mar=c(4.5,4.5,3,1), mgp=c(2.8,0.7,0))
  sub_idx <- if (nrow(ld_result)>5000) sample(nrow(ld_result),5000) else seq_len(nrow(ld_result))
  graphics::plot(ld_result$dist_bp[sub_idx]/1e3, ld_result$r2[sub_idx],
                 pch=".", col=grDevices::adjustcolor("#2563EB",0.3),
                 xlab="Distance (kb)", ylab=expression(r^2),
                 main="LD Decay", las=1, ylim=c(0,1))
  graphics::lines(mids/1e3, mean_r2, col="#DC2626", lwd=2.5)
  graphics::legend("topright", legend=c(expression(paste("Per-pair ",r^2)),"Mean"),
                   col=c(grDevices::adjustcolor("#2563EB",0.5),"#DC2626"),
                   pch=c(16,NA), lty=c(NA,1), lwd=c(NA,2.5), bty="n", cex=0.85)
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}

#' Plot Inbreeding Coefficient F_IS Across Generations
#'
#' @description
#' Line chart of the per-generation \eqn{F_{IS}} statistic with a reference
#' line for the input inbreeding coefficient \eqn{F}.
#'
#' @param sim_result List from \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' sim <- simulate_population(n_founders=100, n_snps=500, n_generations=5,
#'                            inbreeding_F=0.125, seed=1)
#' plot_fis_trajectory(sim)
#' @export
plot_fis_trajectory <- function(sim_result, out_file = NULL) {
  ss  <- sim_result$summary_stats
  fis <- ss$inbreeding_fis
  if (!is.null(out_file)) grDevices::pdf(out_file, width=7, height=4.5)
  graphics::par(mar=c(4,4.5,3,1), mgp=c(2.5,0.7,0))
  graphics::plot(ss$generation, fis, type="b", pch=19, col="#7C3AED", lwd=2,
                 xlab="Generation", ylab=expression(F[IS]),
                 main=expression(paste("Inbreeding Coefficient ",F[IS])),
                 ylim=c(min(c(fis,0),na.rm=TRUE)-0.01, max(c(fis,0.1),na.rm=TRUE)+0.01),
                 xaxt="n", las=1)
  graphics::axis(1, at=ss$generation, labels=paste0("Gen ",ss$generation))
  F_ref <- if (!is.null(sim_result$params$inbreeding_F)) sim_result$params$inbreeding_F else 0
  graphics::abline(h=F_ref, lty=2, col="#6B7280", lwd=1.5)
  graphics::legend("topright", legend=c(expression(F[IS]), paste0("Input F=",F_ref)),
                   col=c("#7C3AED","#6B7280"), lty=c(1,2), pch=c(19,NA), lwd=c(2,1.5),
                   bty="n", cex=0.85)
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}

#' Six-Panel Summary Dashboard
#'
#' @description
#' Produces a combined 2?-3 panel figure showing allele frequency drift,
#' heterozygosity change, \eqn{F_{IS}} trajectory, mean MAF, fraction of
#' fixed loci, and either a PCA scatter or individual counts per generation.
#'
#' @param sim_result List from \code{\link{simulate_population}} or
#'   \code{\link{simulate_from_pedigree}}.
#' @param pca_result Optional list from \code{\link{run_pca}} for Panel 6.
#'   If \code{NULL}, a bar chart of individual counts is shown instead.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' sim <- simulate_population(n_founders=100, n_snps=500, n_generations=5, seed=1)
#' all_g <- do.call(rbind, sim$genotypes)
#' pca <- run_pca(all_g)
#' plot_dashboard(sim, pca_result = pca)
#' @seealso \code{\link{plot_af_trajectory}}, \code{\link{plot_pca}}
#' @export
plot_dashboard <- function(sim_result, pca_result = NULL, out_file = NULL) {
  if (!is.null(out_file)) grDevices::pdf(out_file, width=14, height=9)
  graphics::par(mfrow=c(2,3), oma=c(0,0,3,0))
  ss   <- sim_result$summary_stats
  af   <- sim_result$allele_freqs
  gens <- 0:(nrow(af)-1)

  # P1: AF drift
  graphics::par(mar=c(3.5,4,2.5,1), mgp=c(2.2,0.6,0))
  idx <- if (ncol(af)>40) sample(ncol(af),40) else seq_len(ncol(af))
  graphics::plot(NA, xlim=c(0,max(gens)), ylim=c(0,1), xlab="Generation",
                 ylab="Alt allele freq", main="Allele Frequency Drift",
                 xaxt="n", las=1, cex.main=0.95)
  graphics::axis(1, at=gens, labels=gens, cex.axis=0.8)
  for (j in idx) graphics::lines(gens, af[,j],
                                   col=grDevices::adjustcolor(.SIM_COLORS[(j%%10)+1],0.3), lwd=0.8)
  graphics::lines(gens, rowMeans(af), col="#111827", lwd=2)

  # P2: Heterozygosity
  graphics::par(mar=c(3.5,4,2.5,1), mgp=c(2.2,0.6,0))
  graphics::plot(ss$generation, ss$exp_heterozygosity, type="b", pch=1, lty=2,
                 col="#6B7280", lwd=1.5, ylim=c(0,max(ss$exp_heterozygosity)*1.2),
                 xlab="Generation", ylab="Heterozygosity", main="Heterozygosity",
                 cex.main=0.95, xaxt="n", las=1)
  graphics::axis(1, at=ss$generation, labels=ss$generation, cex.axis=0.8)
  graphics::lines(ss$generation, ss$obs_heterozygosity, type="b", pch=19,
                  col="#2563EB", lwd=2)
  graphics::legend("topright", c("Expected","Observed"), col=c("#6B7280","#2563EB"),
                   lwd=c(1.5,2), pch=c(1,19), lty=c(2,1), bty="n", cex=0.75)

  # P3: Fis
  graphics::par(mar=c(3.5,4,2.5,1), mgp=c(2.2,0.6,0))
  graphics::plot(ss$generation, ss$inbreeding_fis, type="b", pch=19, col="#7C3AED",
                 lwd=2, xlab="Generation", ylab=expression(F[IS]), main="Inbreeding Fis",
                 las=1, cex.main=0.95, xaxt="n")
  graphics::axis(1, at=ss$generation, labels=ss$generation, cex.axis=0.8)
  graphics::abline(h=0, lty=2, col="#9CA3AF")

  # P4: Mean MAF
  graphics::par(mar=c(3.5,4,2.5,1), mgp=c(2.2,0.6,0))
  graphics::plot(ss$generation, ss$mean_maf, type="b", pch=19, col="#16A34A", lwd=2,
                 xlab="Generation", ylab="Mean MAF", main="Mean Minor Allele Frequency",
                 las=1, cex.main=0.95, xaxt="n")
  graphics::axis(1, at=ss$generation, labels=ss$generation, cex.axis=0.8)

  # P5: Fraction fixed
  graphics::par(mar=c(3.5,4,2.5,1), mgp=c(2.2,0.6,0))
  graphics::plot(ss$generation, ss$frac_fixed*100, type="b", pch=19, col="#DC2626",
                 lwd=2, xlab="Generation", ylab="% Fixed loci", main="Loci Fixed by Drift",
                 las=1, cex.main=0.95, xaxt="n")
  graphics::axis(1, at=ss$generation, labels=ss$generation, cex.axis=0.8)

  # P6: PCA or counts
  graphics::par(mar=c(3.5,4,2.5,1), mgp=c(2.2,0.6,0))
  if (!is.null(pca_result)) {
    sc <- pca_result$scores; vp <- pca_result$variance_pct
    gl <- sub("_ind.*","", rownames(sc))
    cols <- .SIM_COLORS[(as.integer(as.factor(gl))%%10)+1]
    graphics::plot(sc[,1], sc[,2],
                   xlab=sprintf("PC1 %.1f%%",vp[1]), ylab=sprintf("PC2 %.1f%%",vp[2]),
                   main="PCA (all generations)", las=1, col=cols, pch=19, cex=0.5, cex.main=0.95)
  } else {
    graphics::barplot(ss$n_individuals, names.arg=paste0("G",ss$generation),
                      col=.SIM_COLORS[seq_len(nrow(ss))%%10+1], border=NA,
                      main="Individuals per Generation", ylab="Count", cex.main=0.95, las=1)
  }
  graphics::mtext("GenoSim --- Population Genetics Dashboard", outer=TRUE, cex=1.2, font=2)
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}

#' Draw a Family Pedigree Tree
#'
#' @description
#' Renders a pedigree as a diagram using base R graphics. Squares represent
#' males, circles females, and diamonds individuals of unknown sex. Affected
#' individuals are filled red; unaffected are blue. Founders are marked with
#' an asterisk. Inbreeding coefficients are annotated on each symbol.
#'
#' @param ped A \code{data.frame} from \code{\link{read_pedigree}}.
#' @param highlight_ids Character vector of individual IDs to highlight with
#'   a purple border. Useful for marking probands or affected individuals.
#' @param title Character. Plot title. Default \code{"Family Pedigree"}.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' \dontrun{
#' ped <- read_pedigree("family_pedigree.csv")
#' affected <- ped$individual_id[ped$phenotype == 2]
#' plot_pedigree_tree(ped, highlight_ids = affected)
#' }
#' @seealso \code{\link{read_pedigree}}, \code{\link{plot_kinship_heatmap}}
#' @export
plot_pedigree_tree <- function(ped, highlight_ids = NULL,
                                title = "Family Pedigree", out_file = NULL) {
  gen_levels <- sort(unique(ped$generation))
  n_gen      <- length(gen_levels)
  max_in_gen <- max(table(ped$generation))
  if (!is.null(out_file)) grDevices::pdf(out_file, width=max(12, max_in_gen*1.4), height=max(8, n_gen*2))
  graphics::par(mar=c(1,2,3,1), bg="white")
  graphics::plot(NA, xlim=c(0,max_in_gen+1), ylim=c(n_gen-0.5,-0.5),
                 xaxt="n", yaxt="n", xlab="", ylab="", main=title, bty="n", cex.main=1.1)
  pos_map <- list()
  for (g in gen_levels) {
    gi  <- ped$individual_id[ped$generation==g]
    xp  <- seq((max_in_gen-length(gi))/2+1, by=1, length.out=length(gi))
    for (k in seq_along(gi)) pos_map[[gi[k]]] <- c(xp[k], g)
  }
  for (i in seq_len(nrow(ped))) {
    fid <- ped$father_id[i]; mid <- ped$mother_id[i]
    ind <- ped$individual_id[i]
    if (!fid %in% c("0","NA","") && fid %in% names(pos_map) &&
        !mid %in% c("0","NA","") && mid %in% names(pos_map)) {
      fp <- pos_map[[fid]]; mp <- pos_map[[mid]]; ip <- pos_map[[ind]]
      mx <- (fp[1]+mp[1])/2
      graphics::segments(fp[1],fp[2],mp[1],mp[2],col="#9CA3AF",lwd=1.2)
      graphics::segments(mx,fp[2],mx,(fp[2]+ip[2])/2,col="#6B7280",lwd=1)
      graphics::segments(mx,(fp[2]+ip[2])/2,ip[1],ip[2]+0.12,col="#6B7280",lwd=1)
    }
  }
  for (i in seq_len(nrow(ped))) {
    ind <- ped$individual_id[i]
    if (!ind %in% names(pos_map)) next
    p      <- pos_map[[ind]]
    ph     <- ped$phenotype[i]; if (is.na(ph)) ph <- -9L
    sx     <- ped$sex[i];       if (is.na(sx)) sx <- 0L
    bg_col <- if (ph==2) "#DC2626" else if (ph==1) "#DBEAFE" else "#F3F4F6"
    bd_col <- if (!is.null(highlight_ids) && ind %in% highlight_ids) "#7C3AED" else "#374151"
    lw     <- if (!is.null(highlight_ids) && ind %in% highlight_ids) 2.5 else 1.2
    if (sx==1)
      graphics::rect(p[1]-.3,p[2]-.25,p[1]+.3,p[2]+.25, col=bg_col, border=bd_col, lwd=lw)
    else if (sx==2)
      graphics::symbols(p[1],p[2], circles=.28, inches=FALSE, add=TRUE, fg=bd_col, bg=bg_col, lwd=lw)
    else
      graphics::polygon(c(p[1],p[1]+.3,p[1],p[1]-.3), c(p[2]-.25,p[2],p[2]+.25,p[2]),
                        col=bg_col, border=bd_col, lwd=lw)
    F_v <- ped$inbreeding_F[i]
    lbl <- if (!is.na(F_v) && F_v>.01) sprintf("%s\nF=%.3f",ind,F_v) else ind
    graphics::text(p[1],p[2],lbl,cex=0.5,col="#111827")
    if (ped$is_founder[i]) graphics::text(p[1]+.33,p[2]-.28,"*",cex=0.7,col="#16A34A")
  }
  for (g in gen_levels)
    graphics::text(.3,g,paste0("Gen ",g),cex=0.75,col="#6B7280",adj=0)
  graphics::legend("topright", bty="n", cex=0.72,
                   legend=c("Affected","Unaffected","Missing","Founder (*)","Highlighted"),
                   fill=c("#DC2626","#DBEAFE","#F3F4F6","white","white"),
                   border=c(rep("#374151",4),"#7C3AED"), lwd=c(1,1,1,1,2.5))
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}

#' Plot Kinship Coefficient Heatmap
#'
#' @description
#' Colour heatmap of a pairwise kinship proxy matrix derived from the pedigree
#' structure (parent-offspring = 0.25, full siblings = 0.25,
#' self = \eqn{(1+F)/2}).
#'
#' @param ped A \code{data.frame} from \code{\link{read_pedigree}}.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly, the kinship matrix.
#' @examples
#' \dontrun{
#' ped <- read_pedigree("family_pedigree.csv")
#' K   <- plot_kinship_heatmap(ped)
#' }
#' @seealso \code{\link{plot_pedigree_tree}}
#' @export
plot_kinship_heatmap <- function(ped, out_file = NULL) {
  ids <- ped$individual_id; n <- length(ids)
  K   <- matrix(0.0, n, n, dimnames=list(ids,ids))
  diag(K) <- 0.5 + ped$inbreeding_F/2
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i==j) next
    if (ids[i]==ped$father_id[j]||ids[i]==ped$mother_id[j]) K[i,j]<-K[j,i]<-0.25
    fi<-ped$father_id[i];mi<-ped$mother_id[i]
    fj<-ped$father_id[j];mj<-ped$mother_id[j]
    if (!fi%in%c("0","")&&!mi%in%c("0","")&&fi==fj&&mi==mj) K[i,j]<-K[j,i]<-0.25
  }
  if (!is.null(out_file))
    grDevices::pdf(out_file, width=max(7,n*.4+2), height=max(6,n*.4+2))
  pal    <- grDevices::colorRampPalette(c("#F0F9FF","#BAE6FD","#0284C7","#1E3A5F"))(100)
  graphics::par(mar=c(6,6,3,4))
  graphics::image(1:n,1:n,t(K[n:1,]), col=pal, breaks=seq(0,.5,length.out=101),
                  xaxt="n",yaxt="n",xlab="",ylab="",main="Pairwise Kinship Matrix")
  graphics::axis(1,at=1:n,labels=ids,las=2,cex.axis=0.6)
  graphics::axis(2,at=1:n,labels=rev(ids),las=1,cex.axis=0.6)
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(K)
}

#' Plot ROH Length Per Individual
#'
#' @description
#' Bar chart of total ROH length (in Mb) per individual, coloured by affection
#' status when a pedigree is supplied. Affected individuals are shown in red.
#'
#' @param roh_result List returned by \code{\link{detect_roh}}.
#' @param ped Optional \code{data.frame} from \code{\link{read_pedigree}} for
#'   colouring affected individuals.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' sim <- simulate_population(n_founders=50, n_snps=500, n_generations=2,
#'                            inbreeding_F=0.25, chromosomes=1:5, seed=1)
#' roh <- detect_roh(sim$genotypes[[3]], sim$snp_map, min_snps=5)
#' plot_roh_per_individual(roh)
#' @seealso \code{\link{detect_roh}}
#' @export
plot_roh_per_individual <- function(roh_result, ped = NULL, out_file = NULL) {
  segs <- roh_result$roh_segments
  if (nrow(segs)==0) { message("No ROH segments."); return(invisible(NULL)) }
  roh_mb <- sort(tapply(segs$length_bp, segs$ind_id, sum)/1e6, decreasing=TRUE)
  cols   <- rep("#2563EB", length(roh_mb))
  if (!is.null(ped)) {
    aff <- ped$individual_id[which(ped$phenotype==2)]   # which() drops NA safely
    cols[names(roh_mb) %in% aff] <- "#DC2626"
  }
  if (!is.null(out_file))
    grDevices::pdf(out_file, width=max(8,length(roh_mb)*.4+2), height=5)
  graphics::par(mar=c(7,5,3,1), mgp=c(3,0.6,0))
  graphics::barplot(roh_mb, col=cols, border=NA, las=2, cex.names=0.65,
                    ylab="Total ROH (Mb)", main="Total ROH Length per Individual")
  graphics::abline(h=mean(roh_mb), lty=2, col="#6B7280", lwd=1.5)
  graphics::legend("topright", bty="n", cex=0.8,
                   legend=c("Individual","Affected","Mean"),
                   fill=c("#2563EB","#DC2626",NA), border=NA, lty=c(NA,NA,2),
                   col=c(NA,NA,"#6B7280"), lwd=c(NA,NA,1.5))
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}

#' Plot Observed vs Simulated Allele Frequency Comparison
#'
#' @description
#' Two-panel plot: (1) scatter of observed (VCF founder) MAF vs simulated MAF
#' in the last generation with \eqn{R^2}; (2) overlaid density curves for the
#' two MAF distributions.
#'
#' @param vcf_cohort List from \code{\link{read_vcf_cohort}}.
#' @param ped_sim_result List from \code{\link{simulate_from_pedigree}}.
#' @param out_file Character or \code{NULL}. PDF output path.
#'
#' @return Invisibly \code{NULL}.
#' @examples
#' \dontrun{
#' vcf <- read_vcf_cohort("family_vcfs/")
#' ped <- read_pedigree("family_pedigree.csv")
#' sim <- simulate_from_pedigree(vcf, ped)
#' plot_af_comparison(vcf, sim)
#' }
#' @seealso \code{\link{simulate_from_pedigree}}
#' @export
plot_af_comparison <- function(vcf_cohort, ped_sim_result, out_file = NULL) {
  p_obs   <- colMeans(vcf_cohort$geno_matrix, na.rm=TRUE)/2
  maf_obs <- pmin(p_obs,1-p_obs)
  n_gen   <- length(ped_sim_result$genotypes)
  last    <- ped_sim_result$genotypes[[n_gen]]
  if (is.null(last)||nrow(last)==0) last <- ped_sim_result$genotypes[[max(1,n_gen-1)]]
  common  <- intersect(colnames(last), names(p_obs))
  maf_sim <- pmin(colMeans(last[,common,drop=FALSE],na.rm=TRUE)/2, 0.5)
  maf_obs_sub <- maf_obs[common]
  if (!is.null(out_file)) grDevices::pdf(out_file, width=10, height=5)
  graphics::par(mfrow=c(1,2), mar=c(4.5,4.5,3,1), mgp=c(2.8,0.7,0))
  graphics::plot(maf_obs_sub, maf_sim, pch=".", col=grDevices::adjustcolor("#2563EB",.4),
                 xlab="Observed MAF (VCF)", ylab="Simulated MAF (last gen)",
                 main="Observed vs Simulated MAF", las=1)
  graphics::abline(a=0,b=1,col="#DC2626",lwd=1.5,lty=2)
  r2 <- stats::cor(maf_obs_sub, maf_sim, use="complete.obs")^2
  graphics::legend("topleft",legend=sprintf("R?? = %.3f",r2),bty="n",cex=0.9)
  d1 <- stats::density(maf_obs_sub[maf_obs_sub>0&maf_obs_sub<.5],na.rm=TRUE)
  d2 <- stats::density(maf_sim[maf_sim>0&maf_sim<.5],na.rm=TRUE)
  graphics::plot(d1,col="#2563EB",lwd=2,ylim=c(0,max(d1$y,d2$y)*1.15),
                 main="MAF Density Comparison",xlab="MAF",las=1)
  graphics::lines(d2,col="#DC2626",lwd=2,lty=2)
  graphics::legend("topright",legend=c("Observed","Simulated"),
                   col=c("#2563EB","#DC2626"),lwd=2,lty=c(1,2),bty="n",cex=0.85)
  if (!is.null(out_file)) { grDevices::dev.off(); message("Saved: ", out_file) }
  invisible(NULL)
}
