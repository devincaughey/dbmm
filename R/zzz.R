.onAttach <- function(...) {
  pkgname = "dynIRTtest"
  packageStartupMessage("Package version: ",
                        as.character(utils::packageVersion(pkgname)))
}
