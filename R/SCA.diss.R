#' Calculate the spectral contrast angle (SCA) between all cells in a hyperSpec object
#'
#'  The SCA measures the angle between two vectors corresponding to closely
#'  related spectra to measure whether they are the same or not
#'  (Wan *et al.* 2002). This function returns a dist object based on the SCA
#'  for downstream processing.
#'
#' @param x a hyperSpec object containing the cells in the rows, and wavenumbers (features) in columns
#' @importFrom stats as.dist
#' @return a distance object containing the SCA between each cell
#' @references Wan, K. X., Vidavsky, I., & Gross, M. L. (2002). Comparing similar spectra: from similarity index to spectral contrast angle. Journal of the American Society for Mass Spectrometry, 13(1), 85-88.
#' @examples
#' ## Short example
#'
#' # Load and pretransform MassSpectrum object
#' data("mass.spectra.baseline.corr")
#' mass.spectra.baseline.corr <- wlcutter(mass.spectra.baseline.corr)
#' mq.norm <- calibrateIntensity(mass.spectra.baseline.corr, method="TIC",range=c(600, 1800))
#' # Convert to hyperSpec object
#' hs.norm <- mq2hs(mq.norm)
#' disst <- SCA.diss(hs.norm)
#' @export

SCA.diss <- function(x){
  if(!class(x)=="hyperSpec"){
    stop("You did not supply a valid hyperSpec object, and there is no default.")
  }
  matr <- x@data$spc
  diss <- matrix( nrow = nrow(matr), ncol = nrow(matr))
  # there has to be a way to get rid of this double for-loop
  # TODO: have a look at cov and see how it's done there (apparently with C calls)
  for (i in 1:nrow(matr)){
    for (j in 1:nrow(matr)) {
      diss[i,j] <- SCA(matr[i,], matr[j,])
    }
  }
  row.names(diss) <- rownames(x@data$spc)
  diss <- as.dist(diss)
  attr(diss, "method") <- "SCA.diss"
  return(diss)
}
