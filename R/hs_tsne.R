#' Run tSNE on a hyperspec object
#'
#' @param hs.x hyperSpec object or alternatively distance matrix
#' @param PCA Should PCA be used for initial dimension reduction? Defaults to FALSE.
#' @param ... Parameters to pass on to tsne().
#' @importFrom tsne tsne
#' @examples
#' ## Short example
#'
#' # Load hyperSpec object
#' data("hs_example")
#'
#' # Convert to MassSpectrum object
#' hs.x.proc <- hs_preprocess(hs_example)
#'
#' df.tsne <- hs_tsne(hs.x.proc)
#' @export

hs_tsne <- function(hs.x, PCA = FALSE, ...){

  # Extract spectra ID
  if (is.null(rownames(hs.x@data$spc)))
    spectra_ids <-
      hs.x@data$filename
  else
    spectra_ids <- rownames(hs.x@data$spc)

  # Run PCA if necessary

  # Run tSNE
  hs.tsne <- tsne(X = hs.x@data$spc, ...)
  hs.tsne.results <- data.frame(Spectrum_ID = spectra_ids,
    hs.tsne)
  return(hs.tsne.results)
}
