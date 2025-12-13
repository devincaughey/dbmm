.onAttach <- function(...) {
    pkgname <- "dbmm"
    packageStartupMessage("Package version: ",
                          as.character(utils::packageVersion(pkgname)))
}

