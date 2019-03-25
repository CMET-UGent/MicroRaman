#' Extract model fit statistics from a linear model
#'
#' This function extracts the R squared, adjusted R squared, predicted R squared,
#' and predictive residuals from a linear model.
#' @param linear.model an lm-type object
#' @importFrom stats summary.lm lm
#' @return A data.frame containing the R squared, adj. R squared, pred. R squared,
#'     and predictive residuals
#' @examples
#' ## Short example
#'
#' @export


model_fit_stats <- function(linear.model) {
  if(!class(linear.model)=="lm"){stop("No valid lm object was given")}
  r.sqr <- summary(linear.model)$r.squared
  adj.r.sqr <- summary(linear.model)$adj.r.squared
  pre.r.sqr <- pred_r_squared(linear.model)
  PRESS <- PRESS(linear.model)
  return.df <- data.frame(r.squared = r.sqr, adj.r.squared = adj.r.sqr,
                          pred.r.squared = pre.r.sqr, press = PRESS)
  return(return.df)
}
