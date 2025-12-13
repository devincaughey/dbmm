# a more than usually informative error message for handing in the
# wrong type to a function
# from keyATM package
check_arg_type <- function(arg, typename, message = NULL) {
  argname <- deparse(match.call()[['arg']])
  if (!inherits(arg, typename)) {
    if (is.null(message)) {
      cli::cli_abort(paste0('`', argname, '` is not a ', typename, ' object.'))
    } else {
      cli::cli_abort(message)
    }
  }
}
