#' Run PCA on a hyperspec object
#'
#' @param hs.x hyperSpec object
#' @param ... Parameters to pass on to prcomp().
#' @importFrom stats prcomp
#' @examples
#' ## Short example
#'
#' # Load hyperSpec object
#' data("hs.example")
#'
#' # Convert to MassSpectrum object
#' hs.x.proc <- hs_preprocess(hs.x)
#'
#' df.PCA <- hs_PCA(hs.x.proc)
#' @export


hs_PCA <- function(hs.x, ...) {
  if (is.null(rownames(hs.x@data$spc)))
    spectra_ids <-
      hs.x@data$filename
  else
    spectra_ids <- rownames(hs.x@data$spc)
  hs.PCA <- prcomp(x = hs.x@data$spc, ...)
  hs.PCA.results <- data.frame(Spectrum_ID = spectra_ids,
    hs.PCA$x)
  return(hs.PCA.results)
}
