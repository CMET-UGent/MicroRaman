#' Calculate the phenotypic heterogeneity indices for each cell
#'
#' @param hs.x Available hyperspec object
#' @param smooth Smoothing TRUE/FALSE
#' @param align Algn spectra TRUE/FALSE
#' @param niter Number of iterations.
#' @param peak_detection d
#' @param peak_window d
#' @param peak_method d
#' @param d d
#' @importFrom MALDIquant isMassSpectrum intensity createMassSpectrum mass calibrateIntensity trim estimateBaseline removeBaseline
#' @importFrom  MALDIquant alignSpectra
#' @importFrom tidyr spread
#' @importFrom graphics lines legend plot
#' @importFrom hyperSpec hy.setOptions spc.loess orderwl
#' @examples
#'
#' ## Short example
#'
#' # Load hyperSpec object
#'
#' @export

hs_preprocess <- function(hs.x,
  smooth = FALSE,
  align = FALSE,
  path = NULL,
  pattern = ".spc",
  trim.range = c(400, 1800),
  niter = 10) {

  hy.setOptions(file.keep.name = TRUE)

  # Load data if path is given
  if (!is.null(path)) {
    # DATA (spc files)
    filenames <- list.files(path = path, pattern = pattern)

    # read the first spectrum and order the wavelenghts with orderwl()
    hs.x <- orderwl(read.spc(paste0(path, filenames[1])))

    # rbind to add spectra to data
    for (i in 2:length(filenames)) {
      hs.x <- rbind(hs.x, orderwl(read.spc(paste0(
        path, filenames[i]
      ))))
    }

    labels <- filenames
    cell.name <-
      unlist(lapply(strsplit(labels, split = "_"), function(x)
        (x[1])))
  }

  # Preprocess data according to Garcia-Timermans et al. (2020)
  if (smooth) {
    # Smooth spectra
    hs.x <-
      spc.loess(hs.x, c(seq (0, hs.x@wavelength[length(hs.x@wavelength)])))
    plot(hs.x)
  }

  # hyperspec --> maldiquant
  mq <- list()
  for (i in 1:length(filenames)) {
    mass.spectrum <-
      createMassSpectrum(mass = hs.x@wavelength,
        intensity = hs.x@data$spc[i,])
    metaData(mass.spectrum) <- list(name = cell.name[i])
    mq <- c(mq, mass.spectrum)
  }

  ## TRIMMING
  # select the relevant part of the fingerprint
  mq.trim <- trim(mq, range = trim.range)
  wavelengths.trim <-
    mass(mq.trim[[1]]) #333 (intervals are unequally distributed)

  # BASELINE CORRECTION
  number.of.iterations <- niter
  plot(
    mq.trim[[1]],
    main = "SNIP baseline correction",
    xlab = expression("Wavenumber (cm" ^ -1 * ")"),
    ylab = "Intensity (AU)",
    ylim = c(min(intensity(mq.trim[[i]])), max(intensity(mq.trim[[i]])))
  )
  baseline <-
    estimateBaseline(mq.trim[[1]], method = "SNIP", iterations = number.of.iterations)
  lines(baseline, col = "blue", lwd = 2)
  mass.spectra.baseline.corr <-
    removeBaseline(mq.trim, method = "SNIP", iterations = number.of.iterations)

  # plot a spectrum to see the effect
  plot(
    mass.spectra.baseline.corr[[1]],
    main = "SNIP baseline correction",
    xlab = expression("Wavenumber (cm" ^ -1 * ")") ,
    ylab = "Intensity (AU)",
    ylim = c(min(intensity(
      mass.spectra.baseline.corr[[i]]
    )), max(intensity(
      mass.spectra.baseline.corr[[i]]
    )))
  )

  # NORMALIZATION:
  # surface = 1
  mq.norm <-
    calibrateIntensity(mass.spectra.baseline.corr,
      method = "TIC",
      range = trim.range)
  plot(
    mq.norm[[1]],
    main = "Normalization",
    xlab = expression("Wavenumber (cm" ^ -1 * ")") ,
    ylab = "Intensity (AU)",
    ylim = c(min(intensity(mq.norm[[i]])), max(intensity(mq.norm[[i]])))
  )
  par(mfrow = c(2, 1))

  # Aligned spectra
  #The first spectra is used as reference
  spectra_aligned <-
    alignSpectra(
      mq.norm ,
      tolerance = 5,
      halfWindowSize = 1,
      reference = mq.norm[[1]]
    )

  #### change MALDIquant (mq.norm) object in Hyperspec (hs.norm) object ####
  # get the intensity matrix from the mq.norm object
  matrix.spectra <-
    matrix(nrow = length(spectra_aligned),
      ncol = length(wavelengths.trim))
  for (i in 1:length(spectra_aligned)) {
    matrix.spectra[i,] <- intensity(spectra_aligned[[i]])
  }

  hs.x <- new(
    "hyperSpec",
    spc = matrix.spectra / max(matrix.spectra),
    wavelength = wavelengths.trim,
    labels = labels
  )

  return(hs.x)
}
