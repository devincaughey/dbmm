## vignettes/precompute.R
##
## Renders modgirt.Rmd.orig to modgirt.Rmd, evaluating all chunks and baking in the
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
    knitr::knit("modgirt.Rmd.orig", output = "modgirt.Rmd")
})

## Figures are written to vignettes/ and must be committed alongside modgirt.Rmd,
## since they are referenced by the built vignette.
cat("\nDone. Commit vignettes/modgirt.Rmd and vignettes/*.png.\n")
