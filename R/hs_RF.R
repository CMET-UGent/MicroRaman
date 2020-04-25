#' Train random forest model on hyperspec object
#'
#' @param hs.x Hyperspec object
#' @param metadata Dataframe containing the categorical variable/group to predict (target_var)
#' and also a column with the matching spectrum identifiers (spectrumID_col).
#' @param target_var Categorical variable/group to predict in metadata
#' @param spectrumID_col Column with the matching spectrum identifiers (spectrumID_col)
#' in metadata
#' @param p_train Percentage of data to use in training model. Defaults to 0.75.
#' @param ntree Number of trees to build. Defaults to 500.
#' @param metric Metric to use to report/maximize performance of model (only for method_ML = "rf")
#' @importFrom tidyr spread
#' @importFrom caret trainControl createDataPartition train confusionMatrix
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
#' hs.RF <- hs_RF(hs.x = hs_example, metadata = metadata, spectrumID_col= "Spectrum_ID",
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
#' @export
#'
hs_RF <- function(hs.x,
  metadata,
  target_var,
  spectrumID_col,
  ntree = 500,
  p_train = 0.75,
  metric = "Accuracy",
  ...) {

  # Check if hyperspec object
  if (is.null(hs.x) | class(hs.x) != "hyperSpec") {
    stop(
      "Error: you did not supply a valid hyperSpec object,
      and there is no default, please correct"
    )
  }

  # Sort metadata ccording to spectral data
  metadata <- metadata[match(rownames(hs.x@data$spc), metadata[, spectrumID_col]), ]

  # Combine the data
  full_data <- data.frame(hs.x@data$spc, metadata)

  # Create data partitions
  trainIndex <- createDataPartition(full_data[, target_var], p = p_train)
  train_data <- full_data[trainIndex$Resample1, ]
  test_data <- full_data[-trainIndex$Resample1, ]

  # Set model parameters
  fitControl <- trainControl( ## 10-fold CV
    method = "repeatedcv",
    number = 10, ## repeated ten times
    repeats = 3)

  # Train Random Forest classifier on training set
  metric <- metric
  mtry <- round(base::sqrt(ncol(train_data)-1), 0)
  tunegrid <- base::expand.grid(.mtry = mtry)

  RF_train <- train(y = train_data[, target_var],
    x = train_data[, !colnames(train_data) %in% c(target_var, spectrumID_col)],
    method = "rf",
    metric = metric, tuneGrid = tunegrid,
    trControl = fitControl, ntree = ntree)

  # Predict on test data set for confusion matrix:
  RF_pred <- stats::predict(RF_train, newdata = test_data)

  # Make list containing model, confusion matrix, summary statistics and descision boundary
  results_list <- list()
  results_list[[1]] <- RF_train
  results_list[[2]] <- confusionMatrix(data = RF_pred, test_data[, target_var])

  return(results_list)
}
