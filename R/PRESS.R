#' Calculated predicted residuals from a linear model
#'
#' This function calculates the predicted resiiduals from a linear model
#' @param linear.model an lm-type object
#' @importFrom stats summary.lm lm anova residuals lm.influence
#' @return A numeric with the predicted residuals


PRESS <- function(linear.model) {
  # calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  # calculate the PRESS
  PRESS <- sum(pr^2)

  return(PRESS)
}
