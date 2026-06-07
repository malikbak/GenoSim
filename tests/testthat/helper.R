# tests/testthat/helper.R — loaded automatically before all test files
# When running via devtools::test() the package is loaded; when running
# plain test_dir() we source everything manually.
if (!isNamespaceLoaded("GenoSim")) {
  pkg_r <- list.files(
    file.path(dirname(dirname(dirname(getwd()))), "R"),
    full.names = TRUE, pattern = "\\.R$"
  )
  # fallback: look relative to testthat dir
  if (length(pkg_r) == 0) {
    pkg_r <- list.files(
      file.path(getwd(), "..", "..", "R"),
      full.names = TRUE, pattern = "\\.R$"
    )
  }
  invisible(lapply(pkg_r, source, local = globalenv()))
  # also expose example data helpers if not yet available
  if (!exists("example_ped_path")) {
    ext <- file.path(getwd(), "..", "..", "inst", "extdata")
    example_ped_path <- function() file.path(ext, "family_pedigree.csv")
    example_vcf_dir  <- function() file.path(ext, "family_vcfs")
    assign("example_ped_path", example_ped_path, envir = globalenv())
    assign("example_vcf_dir",  example_vcf_dir,  envir = globalenv())
  }
}
