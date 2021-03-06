#' Calculate the phenotypic heterogeneity indices for each cell
#'
#' @param hs.x HyperSpec object
#' @param preprocess Preprocess? TRUE/FALSE. Defaults to FALSE.
#' @param smooth Option for preprocessing. Defaults to FALSE.
#' @param align Align option for preprocessing. Defaults to FALSE.
#' @param niter Number of iteration Option for preprocessing. Defaults to FALSE.
#' @param return_spectra Should hyperSpec object be returned as final object? Defaults to TRUE.
#' @param peak_detection Should peak detection be used instead of raw spectra for Hill calculations? Defaults to FALSE
#' @param peak_window Peak windows size for a peak to be considered a signal. Defaults to 20.
#' @param peak_method Which peak detection method should be used (requires peak_detection == TRUE).
#' @param path path where Raman data is located, only useful if preprocess==TRUE
#' @param trim.range Wavenumber range to trim from raw spectra. Only used if preprocess==TRUE
#' @param pattern (POSIX) regular expression used to select files (defaults to ".spc"), only used if preprocess==TRUE
#' Options are "MAD" and "SuperSmoother"
#' @importFrom MALDIquant isMassSpectrum intensity createMassSpectrum mass calibrateIntensity trim estimateBaseline removeBaseline
#' @importFrom MALDIquant alignSpectra detectPeaks
#' @importFrom tidyr spread
#' @importFrom graphics lines legend plot
#' @importFrom hyperSpec hy.setOptions spc.loess orderwl
#' @examples
#' # Short example
#' data("hs_example")
#'
#' # Preprocess
#' hs_example <- hs_preprocess(hs_example)
#'
#' # Calculate metrics
#' hs_phenoRam(hs.x = hs_example, preprocess = TRUE, peak_detection = TRUE)
#' @export
#'
hs_phenoRam <- function(hs.x,
  preprocess = FALSE,
  smooth = FALSE,
  align = FALSE,
  path = NULL,
  trim.range = c(400, 1800),
  pattern = ".spc",
  niter = 10,
  return_spectra = TRUE,
  peak_detection = FALSE,
  peak_window = 20,
  peak_method = c("MAD")) {

  # Make sure to keep file names on import
  hy.setOptions(file.keep.name = TRUE)

  # Check if hyperspec object
  if (is.null(hs.x) | class(hs.x) != "hyperSpec") {
    stop(
      "Error: you did not supply a valid hyperSpec object,
      and there is no default, please correct"
    )
  }

  # Preprocess data according to Garcia-Timermans et al. (2020)
  if(preprocess){
    hs.x <- hs_preprocess(hs.x,
      smooth = smooth,
      align = align,
      path = path,
      pattern = pattern,
      trim.range = c(400, 1800),
      niter = 10)
  }

  # Peak detection algorithm
  if (peak_detection) {
    # hyperspec --> maldiquant
    mq <- hs_conv_mq(hs.x)
    mq.peaks <- detectPeaks(mq, halfWindowSize = peak_window, method = peak_method)

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
    df_peaks <- tidyr::spread(df_peaks, .data$wavenumber, .data$intensity)
    wavelengths.trim <- as.numeric(colnames(df_peaks)[-1])
    hs.x <-
      new(
        "hyperSpec",
        spc = as.matrix(df_peaks[-1]),
        wavelength = wavelengths.trim,
        labels =  sapply(mq, function(x)
          x@metaData$name)
      )
    rownames(hs.x@data$spc) <- sapply(mq, function(x)
      x@metaData$name)
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
      Spectrum_ID = rownames(hs.x@data$spc),
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
