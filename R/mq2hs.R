#' Convert a MassSpectrum object to a hyperSpec object
#'
#' This function converts a list of MALDIquant::MassSpectrum  objects to a
#' hyperSpec::hyperSpec object for downstream processing. It assumes all
#' MassSpectrum objects have the same wavelengths as the first entry in the list.
#' @param x hyperSpec object
#' @param sampleNames an (optional) character vector of length(x) to designate sample names
#' @importFrom MALDIquant isMassSpectrum mass
#' @importFrom methods new
#' @examples
#' ## Short example
#'
#' # Load and pretransform MassSpectrum object
#' data("mass.spectra.baseline.corr")
#' data("spx.all")
#' mdqs <- lapply(spx.all,hs2mq)
#' labels <-  sub(pattern = ".*/(.*.spc)","\\1",
#' x =sapply(mdqs,function(x)x@metaData$name))
#' Medium <- unlist(lapply(strsplit(labels,split="_"),function(x) (x[2])))
#' Replicate <- unlist(lapply(strsplit(labels,split="_"),function(x) (x[3])))
#' smpnms <- paste(Medium,Replicate)
#' mass.spectra.baseline.corr <- wlcutter(mass.spectra.baseline.corr)
#' mq.norm <- calibrateIntensity(mass.spectra.baseline.corr, method="TIC",range=c(600, 1800))
#' # Convert to hyperSpec object
#' hs.norm <- mq2hs(mq.norm,sampleNames=smpnms)
#' @export

mq2hs <- function(x,sampleNames=NULL){
  if(is.null(x)|sum(sapply(x,isMassSpectrum))!=length(x)){
    stop("Error: you did not supply a valid list of MassSpectrum objects,
         and there is no default, please correct")
  }
  # in the current version of this script, we assume all spectra have the same
  # wavelengths as the first spectrum
  # TODO: remove this assumption in an intelligent way
  matrix.spectra <- matrix(nrow=length(x), ncol = length(mass(x[[1]])))
  for (i in 1:length(x)){
    matrix.spectra[i,] <- intensity(x[[i]])
  }
  cnames <- sapply(x,function(x)x@metaData$name)
  hsobj <- new ("hyperSpec", spc = matrix.spectra,
                wavelength = mass(x[[1]]), labels= cnames)
  rownames(hsobj@data$spc) <- cnames
  colnames(hsobj@data$spc) <- mass(x[[1]])
  if(!is.null(sampleNames)){
    if(!nrow(hsobj@data)==length(sampleNames)){
      stop(paste0("The supplied sample name vector length (",
                  length(sampleNames),
                  ") is not the same as the amount of samples (",
                  length(x),
                  ") in the list of MALDIquant objects"))
    }else{
    hsobj@label <- as.list(sampleNames)
    }
  }
  return(hsobj)
}
