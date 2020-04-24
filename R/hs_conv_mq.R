#' Convert a hyperSpec object toa  MassSpectrum object
#'
#' This function converts a hyperSpec::hyperSpec objects (default of Raman spectra) to a
#' MALDIquant::MassSpectrum object for downstream processing.
#' @param hs.x hyperSpec object
#' @importFrom MALDIquant createMassSpectrum
#' @examples
#' ## Short example
#'
#' # Load hyperSpec object
#' data("hs.example")
#'
#' # Convert to MassSpectrum object
#' mq.x <- hs_conv_mq(hs.example)
#' @export

hs_conv_mq <- function(hs.x){
  # Sanity check
  if(is.null(hs.x)|class(hs.x)!="hyperSpec"){
    stop("Error: you did not supply a valid hyperSpec object,
         and there is no default, please correct")
  }

  mq <- list()
  for (i in 1:length(hs.x)) {
    # Extra necessary features
    wl <- hs.x@wavelength
    sc <- as.vector(hs.x@data$spc[i,])
    fn <- hs.x@data$filename[i]
    if(is.null(fn)){
      fn <- rownames(hs.x@data$spc)[i]
    }
    # Create mass spectrum
    mass.spectrum <-
      createMassSpectrum(mass = wl,
        intensity = sc,
        metaData = list(name = fn)
        )

    # Merge with other spectra in single object
    mq <- c(mq, mass.spectrum)
  }

  # Return list of maldiquant objects
  return(mq)
}
