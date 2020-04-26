#' run_quiet function
#'
#' Function to make other functions with cat() run quiet.
#' Credits to Hadley Wickham.
#' @param x function call that needs to run quiet.
#' @keywords run_quiet
#' @export

run_quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
