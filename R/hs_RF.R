#' Train random forest model on hyperspec object
#'
#' @param hs.x Hyperspec object
#' @param metadata Dataframe containing the categorical variable/group to predict (target_var)
#' and also a column with the matching spectrum identifiers (spectrumID_col).
#' @param target_var Categorical variable/group to predict in metadata
#' @param spectrumID_col Column with the matching spectrum identifiers (spectrumID_col)
#' in metadata
#' @param ptrain Percentage of data to use in training model. Defaults to 0.75.
#' @param ntree Number of trees to build. Defaults to 500.
#'
#' @importFrom tidyr spread
#' @importFrom caret train
#' @importFrom randomForest randomForest
#' @examples
#' # Short example
#' data("hs_example")
#'
#' # Preprocess
#' hs_example <- hs_preprocess(hs_example)
#'
#' # Mock-up metadata
#' mock_meta <- metadata
#'
#' # Calculate metrics
#' hs_RF(hs.x = hs_example, metadata = metadata, spectrumID_col= "Spectrum_ID",
#'  target_var = "test_group")
#' @export
#'
hs_RF <- function(hs.x,
  metadata,
  target_var,
  spectrumID_col,
  ntree = 500,
  ptrain = 0.75,
  ...) {

  # Check if hyperspec object
  if (is.null(hs.x) | class(hs.x) != "hyperSpec") {
    stop(
      "Error: you did not supply a valid hyperSpec object,
      and there is no default, please correct"
    )
  }





  return(hs.RF)
}
