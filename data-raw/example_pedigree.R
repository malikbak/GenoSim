## data-raw/example_pedigree.R
## Run this once to regenerate the bundled R data object.
## usethis::use_data() is not required — we write manually for portability.

example_pedigree <- read.csv(
  system.file("extdata", "family_pedigree.csv", package = "GenoSim"),
  stringsAsFactors = FALSE
)

# For building without package installed, read from inst/extdata directly
if (!nzchar(system.file(package = "GenoSim"))) {
  example_pedigree <- read.csv(
    file.path("inst", "extdata", "family_pedigree.csv"),
    stringsAsFactors = FALSE
  )
}

save(example_pedigree, file = "data/example_pedigree.rda", compress = "xz")
message("Saved data/example_pedigree.rda (", nrow(example_pedigree), " rows)")
