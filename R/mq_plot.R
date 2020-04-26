#' Convenience function to plot Raman spectra (from MALDIquant::MassSpectrum)
#'
#' This convenience function is a wrapper around MALDIquant plot-methods.
#' In essence it does nothing else than rename the axes.
#'
#' @param msspc a MALDIquant::MassSpectrum (trimmed to relevant wavelengths)
#' @param main a string with the plot header (defaults to "Raman Spectrum")
#' @param xlab a string with the horizontal axis label (defaults to inverse wavenumber)
#' @param ylab a string with the vertical axis label (devaults to Intensity (AU))
#' @param ... additional parameters passed on to plot
#' @importFrom MALDIquant isMassSpectrum
#' @importFrom graphics plot
#' @export

ramplot <- function(msspc,main="Raman spectrum",
                    xlab=expression("Wavenumber (cm"^-1*")"),
                    ylab="Intensity (AU)",...){
  if(!isMassSpectrum(msspc)){stop("You did not supply a valid MassSpecrum object")}
  plot(msspc,main=main,xlab=xlab,ylab=ylab,...)
}
