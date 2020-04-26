#' Resampling function for hyperspec objects
#'
#' This function resamples the hyperspec object to an equal number of cells/spectra. Options are to resample
#' with or withour replacement.
#'
#' @param hs.x hyperSpec object.
#' @param sample Desired sample size. Defaults to minimum sample size.
#' @param replace Do you want to resample with or without replacement? Defaults to FALSE, which is without replacement.
#' @param seed Set seed for reproducibile sampling if that is necessary. Defaults to NULL.
#' @return Hyperspec object
#' @keywords resampling, hyperSpec
#' @examples
#' # Short example
#' data("hs_example")
#'
#' # Preprocess
#' hs_example <- hs_preprocess(hs_example)
#'
#' # Resample
#' hs_example_r <- hs_resample(hs_example, sample = 10, replace = TRUE)
#' @export

hs_resample <- function(hs.x,
  sample = 0,
  replace = FALSE,
  seed = NULL) {
  # Set seed if required
  if (!is.null(seed))
    set.seed(seed)

  # Check if hyperspec object
  if (is.null(hs.x) | class(hs.x) != "hyperSpec") {
    stop(
      "Error: you did not supply a valid hyperSpec object,
      and there is no default, please correct"
    )
  }

  if (sample > length(hs.x)) {
    stop(
      "Error: requested sample size is larger than number of spectra"
    )
  }

  # Resample
  s_r <- base::sample(1:length(hs.x), size = sample, replace = replace)
  mat.x.r <- hs.x@data$spc[s_r,]

  # Get file names
  if (!is.null(hs.x$filename))
    spec_names <-
    hs.x$filename
  else
    spec_names <- rownames(hs.x@data$spc)

  # Create new hyperSpec object
  hs.x.r <- new(
    "hyperSpec",
    spc = mat.x.r,
    wavelength = hs.x@wavelength,
    labels = spec_names[s_r]
  )

  rownames(hs.x.r@data$spc) <- spec_names[s_r]

  return(hs.x.r)
  }
