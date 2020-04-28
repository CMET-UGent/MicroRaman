#' Train random forest model on hyperspec object
#'
#' @param hs.x Hyperspec object on which to perform predictions
#' @param model Model to use for predictions. Output from hs_RF.
#' @param ... additional parameters passed on to stats::predict
#' @importFrom tidyr spread
#' @importFrom randomForest randomForest
#' @examples
#' # Short example
#' data("hs_example")
#'
#' # Preprocess
#' hs_example <- hs_preprocess(hs_example)
#'
#' # Mock-up metadata
#' mock_meta <- data.frame(Spectrum_ID = rownames(hs_example@data$spc),
#' group = factor(c(rep(1,30),rep(2,34))))
#'
#' # Calculate metrics
#' hs.RF <- hs_RF(hs.x = hs_example, metadata = mock_meta, spectrumID_col= "Spectrum_ID",
#'  target_var = "group")
#'
#'  # Trained model
#'  print(hs.RF[[1]])
#'
#'  # Confusion matrix
#'  print(hs.RF[[2]])
#'
#'  # Variable importance metric
#'  caret::varImp(hs.RF[[1]])
#'
#'  # Perform predictions
#'  hs_RF_pred(hs.x = hs_example, model = hs.RF[[1]])
#' @export
#'
hs_RF_pred <- function(hs.x,
  model,
  ...) {

  # Check if hyperspec object
  if (is.null(hs.x) | class(hs.x) != "hyperSpec") {
    stop(
      "Error: you did not supply a valid hyperSpec object,
      and there is no default, please correct"
    )
  }

  # Predict on test data set for confusion matrix:
  RF_pred <- stats::predict(model, newdata = data.frame(hs.x@data$spc), ...)

  return(RF_pred)
  }
