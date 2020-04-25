#' Run tSNE on a hyperspec object
#'
#' @param hs.x hyperSpec object or alternatively distance matrix
#' @param ... Parameters to pass on to tsne().
#' @importFrom tsne tsne
#' @examples
#' ## Short example
#'
#' # Load hyperSpec object
#' data("hs.example")
#'
#' # Convert to MassSpectrum object
#' hs.x.proc <- hs_preprocess(hs.x)
#'
#' df.tsne <- hs_tsne(hs.x.proc)
#' @export

hs_tsne <- function(hs.x, ...){
  # Extract spectra ID
  if (is.null(rownames(hs.x@data$spc)))
    spectra_ids <-
      hs.x@data$filename
  else
    spectra_ids <- rownames(hs.x@data$spc)
  # Run tSNE
  hs.tsne <- tsne(X = hs.x@data$spc, ...)
  hs.tsne.results <- data.frame(Spectrum_ID = spectra_ids,
    hs.tsne)
  return(hs.tsne.results)
}
