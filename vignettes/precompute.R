## vignettes/precompute.R
##
## Renders dbmm.Rmd.orig to dbmm.Rmd, evaluating all chunks and baking in the
## output. Run manually whenever the vignette or the package API changes; the
## MODGIRT fit takes several minutes per chain, so this cannot happen during
## R CMD check.
##
## Requires: cmdstanr with a working CmdStan, and estsubpop
##   (pak::pak("devincaughey/estsubpop"))
##
## Usage, from the package root:
##   Rscript vignettes/precompute.R

withr::with_dir("vignettes", {
    knitr::knit("dbmm.Rmd.orig", output = "dbmm.Rmd")
})

## Figures are written to vignettes/ and must be committed alongside dbmm.Rmd,
## since they are referenced by the built vignette.
cat("\nDone. Commit vignettes/dbmm.Rmd and vignettes/*.png.\n")
