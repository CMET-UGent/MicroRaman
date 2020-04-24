#' Convert a MassSpectrum object to a hyperSpec object
#'
#' This function converts a list of MALDIquant::MassSpectrum  objects to a
#' hyperSpec::hyperSpec object for downstream processing. It assumes all
#' MassSpectrum objects have the same wavelengths as the first entry in the list.
#' @param mq.x list of MALDIquant objects
#' @importFrom MALDIquant isMassSpectrum mass
#' @importFrom methods new
#' @examples
#' # Load hyperSpec object
#' data("hs.example")
#'
#' # Convert to MassSpectrum object
#' mq.x <- hs_conv_mq(hs.example)
#'
#' # Convert to hyperspec object
#' hs.x <- mq_conv_hs(mq.x)
#' @export

mq_conv_hs <- function(mq.x) {
  if (is.null(mq.x) |
      sum(sapply(mq.x, isMassSpectrum)) != length(mq.x)) {
    stop(
      "Error: you did not supply a valid list of MassSpectrum objects,
      and there is no default, please correct"
    )
  }

  # Create matrix for holding spectral data from mq list
  matrix.spectra <-
    matrix(nrow = length(mq.x), ncol = length(mass(mq.x[[1]])))

  # Loop through mq list
  matrix.spectra <-
    matrix(unlist(lapply(mq.x, intensity)),
      ncol = length(mass(mq.x[[1]])), byrow = TRUE)
  cnames <- sapply(mq.x, function(x)
    x@metaData$name)
  hs.x <- new (
    "hyperSpec",
    spc = matrix.spectra,
    wavelength = mass(mq.x[[1]]),
    labels = cnames
  )
  rownames(hs.x@data$spc) <- cnames
  colnames(hs.x@data$spc) <- mass(mq.x[[1]])

  return(hs.x)
  }
