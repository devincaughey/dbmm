skip_if_no_cmdstan <- function() {
    testthat::skip_if_not_installed("cmdstanr")
    ok <- tryCatch(!is.null(cmdstanr::cmdstan_version()),
                   error = function(e) FALSE)
    testthat::skip_if_not(ok, "CmdStan not available")
}
