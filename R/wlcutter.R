#' Set intensities of a given bandpass to a preset value (anomalous peak removal)
#'
#' This function converts takes a list of MALDIquant::massSpectrum objects and
#' sets the intensitiy values between certain wavelengths to a given value.
#' By default, this value is zero, as this function can be used to remove errant
#' peaks from the spectra.
#'
#' @param x a list of MALDIquant::massSpectrum objects
#' @param lower_wl an integer of the lower boundary of wavelenghts to replace
#'     intensities for (defaults to 900)
#' @param upper_wl an integer of the upper boundary of wavelengths to replace
#'     intensities for (defaults to 1100)
#' @param new_intens new intensity value to assign to all peaks between lower_wl
#'     and upper_wl. Defaults to 0. Needs to be a single integer value (for now).
#' @return a list of  MALDIquant::massSpectrum objects with replaced intensities
#' @importFrom MALDIquant createMassSpectrum mass intensity metaData<-
#' @export

wlcutter <- function(x,lower_wl=900L,upper_wl=1100L,new_intens=0L){
  # TODO: check if x is a list of MassSpectrum objects
  wlvect <- mass(x[[1]])
  #operates under the assumption that the wavenumbers are identical in each
  #sample
  wl.lower <- wlvect[wlvect<lower_wl] #retain wavelengths below lower threshold
  wl.upper <- wlvect[wlvect>upper_wl] #retain wavelengths above upper threshold
  cut.length <- length(wlvect)-length(wl.lower)-length(wl.upper)
  # integer with the amount of observations to be contained in the cut-out area
  new.spectra <- list()
  for(i in 1:length(x)){
    init.intens <- intensity(x[[i]])
    intens.lower <- init.intens[wlvect<lower_wl]
    intens.upper <- init.intens[wlvect>upper_wl]
    cut.intens <- rep(new_intens,cut.length)
    new.intens <- append(append(intens.lower,cut.intens),intens.upper)
    massspeci <- createMassSpectrum(mass=wlvect,intensity = new.intens)
    metaData(massspeci) <- list(name=x[[i]]@metaData$name)
    new.spectra[[i]] <- massspeci
  }
  return(new.spectra)
}
