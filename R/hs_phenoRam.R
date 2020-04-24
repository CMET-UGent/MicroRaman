#' Calculate the phenotypic heterogeneity indices for each cell
#'
#' @param hs.x d
#' @param preprocess d
#' @param smooth d
#' @param align d
#' @param niter d
#' @param return_spectra d
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
#'
PhenoRam <- function(hs.x,
  preprocess = FALSE,
  smooth = FALSE,
  align = FALSE,
  path = NULL,
  trim.range = c(400, 1800),
  pattern = ".spc",
  niter = 10,
  return_spectra = TRUE,
  peak_detection = TRUE,
  peak_window = 20,
  peak_method = c("MAD")) {

  # Make sure to keep file names on import
  hy.setOptions(file.keep.name = TRUE)

  # Load data if path is given
  if (!is.null(path) && preprocess) {
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
    par(mfrow = c(2,2))
    hs.x <-
      spc.loess(hs.x, c(seq (0, hs.x@wavelength[length(hs.x@wavelength)])))
    plot(hs.x)
    title(main = "Smoothed spectra")
  }

  if (preprocess) {
    plot(hs.x)
    title(main = "Before pre-processing")
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
      ylim = c(min(
        intensity(mass.spectra.baseline.corr[[i]])
      ), max(
        intensity(mass.spectra.baseline.corr[[i]])
      ))
    )

    # NORMALIZATION:
    # surface = 1
    mq.norm <-
      calibrateIntensity(mass.spectra.baseline.corr,
        method = "TIC",
        range = trim.range)
    plot(
      mq.norm[[1]],
      main = "Normalisation",
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
    hs.x <-
      new(
        "hyperSpec",
        spc = matrix.spectra / max(matrix.spectra),
        wavelength = wavelengths.trim,
        labels = labels
      )
  } else
    labels <- labels(hs.x@data)[[1]]

  # Peak detection algorithm
  if (peak_detection) {
    # hyperspec --> maldiquant
    mq <- list()
    for (i in 1:length(hs.x)) {
      mass.spectrum <-
        createMassSpectrum(mass = hs.x@wavelength,
          intensity = hs.x@data$spc[i,])
      metaData(mass.spectrum) <- list(name = cell.name[i])
      mq <- c(mq, mass.spectrum)
    }

    mq.peaks <-
      detectPeaks(mq, halfWindowSize = peak_window, method = peak_method)

    # maldiquant --> hyperspec
    wavelengths.trim <-
      mass(mq.peaks[[1]])
    for (i in 1:length(mq.peaks)) {
      tmp <- data.frame(Spectrum_ID = i,
        wavenumber = mass(mq.peaks[[i]]),
        intensity = intensity(mq.peaks[[i]])
      )
      if(i == 1) {
        df_peaks <- tmp
      } else {
        df_peaks <- rbind(df_peaks, tmp)
      }
    }
    df_peaks <- tidyr::spread(df_peaks, wavenumber, intensity)
    wavelengths.trim <- as.numeric(colnames(df_peaks)[-1])
    hs.x <-
      new(
        "hyperSpec",
        spc = as.matrix(df_peaks[-1]),
        wavelength = wavelengths.trim,
        labels = labels
      )
  }

  # Extract spectra
  x <- hs.x@data$spc

  # Normalize spectra by sum of peak intensities
  # This assumes that each peak width is identical
  x <- x / rowSums(x, na.rm = TRUE)

  # Calculate Hill numbers
  D2.fun <- function(x) {
    x <- x[!is.na(x)]
    1 / sum((x) ^ 2)
  }
  D1.fun <- function(x) {
    x <- x[x != 0]
    x <- x[!is.na(x)]
    exp(-sum((x) * log(x)))
  }
  D0.fun <- function(x) {
    x <- x[!is.na(x)]
    sum(x != 0)
  }

  ### D0
  D0 <- apply(x,
    1,
    FUN = D0.fun)

  ### D1
  D1 <- apply(x,
    1,
    FUN = D1.fun)

  ### D2
  D2 <- apply(x,
    1,
    FUN = D2.fun)

  results <-
    data.frame(
      labels = labels,
      D0 = D0,
      D1 = D1,
      D2 = D2
    )

  colnames(results) = c("Spectrum_ID", "D0", "D1", "D2")

  if (return_spectra) {
    results_comb <- list(results, hs.x)
    return(results_comb)
  } else {
    return(results)
  }

}
