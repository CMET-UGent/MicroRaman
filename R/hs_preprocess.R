#' Calculate the phenotypic heterogeneity indices for each cell
#'
#' @param hs.x Available hyperspec object
#' @param smooth Smoothing TRUE/FALSE
#' @param align Align spectra TRUE/FALSE
#' @param align_ref Spectrum to use as reference for the alignment.
#' Defaults to 1, the first spectrum.
#' @param niter Number of iterations.
#' @param path Path contain spc files if no hyperSpec object is available
#' @param pattern Pattern of spc files to import from path
#' @param trim.range Wavenumber range to trim from raw spectra.
#' @importFrom MALDIquant isMassSpectrum intensity createMassSpectrum mass calibrateIntensity trim estimateBaseline removeBaseline
#' @importFrom  MALDIquant alignSpectra
#' @importFrom tidyr spread
#' @importFrom graphics lines legend plot
#' @importFrom hyperSpec hy.setOptions spc.loess orderwl
#' @importFrom graphics par
#' @examples
#' # Load hyperSpec object
#' data("hs_example")
#'
#' # Load hyperSpec object
#' hs.x.proc <- hs_preprocess(hs.x = hs_example, smooth = FALSE)
#' @export

hs_preprocess <- function(hs.x,
  smooth = FALSE,
  align = FALSE,
  align_ref = 1,
  path = NULL,
  pattern = ".spc",
  trim.range = c(400, 1800),
  niter = 10) {

  # Options
  hy.setOptions(file.keep.name = TRUE)
  # dev.off()
  par(mfrow = c(1, 1))

  # Load data if path is given
  if (!is.null(path)) {
    hs.x <- hs_import(path = path)
  }

  # Preprocess data according to Garcia-Timermans et al. (2020)
  if (smooth) {
    par(mfrow = c(2, 2))
    # Smooth spectra
    hs.x <- spc.loess(hs.x, c(seq (0, hs.x@wavelength[length(hs.x@wavelength)])))
    plot(hs.x)
  }

  # hyperspec --> maldiquant
  mq <- hs_conv_mq(hs.x)

  # TRIMMING
  mq.trim <- trim(mq, range = trim.range)
  wavelengths.trim <-
    mass(mq.trim[[1]]) #333 (intervals are unequally distributed)

  # BASELINE CORRECTION
  plot(
    mq.trim[[1]],
    main = "SNIP pre-baseline correction",
    xlab = expression("Wavenumber (cm" ^ -1 * ")"),
    ylab = "Intensity (AU)",
    ylim = c(min(intensity(mq.trim[[1]])), max(intensity(mq.trim[[1]])))
  )
  baseline <-
    estimateBaseline(mq.trim[[1]], method = "SNIP", iterations = niter)
  lines(baseline, col = "blue", lwd = 2)
  mass.spectra.baseline.corr <-
    removeBaseline(mq.trim, method = "SNIP", iterations = niter)

  # plot a spectrum to see the effect
  plot(
    mass.spectra.baseline.corr[[1]],
    main = "SNIP post-baseline correction",
    xlab = expression("Wavenumber (cm" ^ -1 * ")") ,
    ylab = "Intensity (AU)",
    ylim = c(min(intensity(
      mass.spectra.baseline.corr[[1]]
    )), max(intensity(
      mass.spectra.baseline.corr[[1]]
    )))
  )

  # NORMALIZATION:
  # surface = 1
  mq.norm <-
    calibrateIntensity(mass.spectra.baseline.corr,
      method = "TIC",
      range = trim.range)

  # Diagnostic plot
  plot(
    mq.norm[[1]],
    main = "Normalization",
    xlab = expression("Wavenumber (cm" ^ -1 * ")") ,
    ylab = "Intensity (AU)",
    ylim = c(min(intensity(mq.norm[[1]])), max(intensity(mq.norm[[1]])))
  )

  # Align spectra
  #The first spectra is used as reference
  if(align){
    spectra_aligned <-
      alignSpectra(
        mq.norm,
        tolerance = 5,
        halfWindowSize = 1,
        reference = mq.norm[[align_ref]]
      )
  }

  #### change MALDIquant (mq.norm) object in Hyperspec (hs.norm) object ####
  # get the intensity matrix from the mq.norm object
  hs.x <- mq_conv_hs(spectra_aligned)

  return(hs.x)
}
