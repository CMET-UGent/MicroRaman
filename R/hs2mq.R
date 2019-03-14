#' Convert a hyperSpec object toa  MassSpectrum object
#'
#' This function converts a hyperSpec::hyperSpec objects (default of Raman spectra) to a
#' MALDIquant::MassSpectrum object for downstream processing.
#' @param x hyperSpec object
#' @importFrom MALDIquant createMassSpectrum
#' @examples
#' ## Short example
#'
#' # Load hyperSpec object
#' data("spx.all")
#'
#' # Convert to MassSpectrum object
#' mdqs <- lapply(spx.all,hs2mq)
#' @export

hs2mq <- function(x){
  if(is.null(x)|class(x)!="hyperSpec"){
    stop("Error: you did not supply a valid hyperSpec object,
         and there is no default, please correct")
  }
  wl <- x@wavelength
  sc <- as.vector(x@data$spc)
  fn <- x@data$filename
  mass.spectrum <- createMassSpectrum(mass=wl,
                                      intensity = sc,
                                      metaData = list(name=fn))
  return(mass.spectrum)
}
