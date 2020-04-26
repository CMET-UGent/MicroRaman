#' Calculate the spectral contrast angle (SCA) between two vectors
#'
#'  The SCA measures the angle between two vectors corresponding to closely
#'  related spectra to measure whether they are the same or not
#'  (Wan *et al.* 2002).
#' @param a vector 1
#' @param b vector 2
#' @importFrom NISTunits NISTradianTOdeg
#' @importFrom magrittr "%>%"
#' @references Wan, K. X., Vidavsky, I., & Gross, M. L. (2002). Comparing similar spectra: from similarity index to spectral contrast angle. Journal of the American Society for Mass Spectrometry, 13(1), 85-88.
#' @examples
#' # Load hyperSpec object
#' data("hs_example")
#'
#' # Convert to MassSpectrum object
#' hs.x.proc <- hs_preprocess(hs_example)
#'
#' # Calc SCA
#' disst <- SCA(hs.x.proc@data$spc[1,], hs.x.proc@data$spc[2,])
#' @export

SCA <- function(a, b) {
  numerator <- sum(a * b)
  quada <- sapply(
    a,
    FUN = function(x)
      x ^ 2
  ) %>% sum()
  quadb <- sapply(
    b,
    FUN = function(x)
      x ^ 2
  ) %>% sum()
  denominator <- sqrt(quada * quadb)
  costh <- numerator / denominator
  theta <- NISTradianTOdeg(acos(costh))
  theta <- theta / 90
  return(theta)
}
