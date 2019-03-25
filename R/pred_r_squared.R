#' Calculated predicted R squared from a linear model
#'
#' This function calculates the predicted R squared from a linear model
#' @param linear.model an lm-type object
#' @importFrom stats summary.lm lm anova
#' @return A numeric with the predicted R squared


pred_r_squared <- function(linear.model) {
  if(!class(linear.model)=="lm"){stop("No valid lm object was given")}
  # Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  # Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)

  return(pred.r.squared)
}
