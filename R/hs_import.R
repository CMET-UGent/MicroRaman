#' Import Raman spectral data
#'
#' This function imports raman data (Thermo Galactic's spc file format) into the R environment
#' @param hs.x hyperSpec object
#' @importFrom hyperSpec orderwl read.spc
#' @examples
#' # Import Raman spectral data
#' hs.x <- hs_import(path = "")
#' @export

hs_import <- function(path = NULL) {
  # Sanity check
  if (is.null(path)) {
    stop("Please specify the correct path containing the spc files")
  }
  # Load data if path is given
  if (!is.null(path)) {
    # filenames (spc files)
    filenames <- list.files(path = path, pattern = pattern)

    # read the first spectrum and order the wavelenghts with orderwl()
    hs.x <- orderwl(read.spc(paste0(path, filenames[1])))

    # rbind to add spectra to data
    for (i in 2:length(filenames)) {
      hs.x <- rbind(hs.x, orderwl(read.spc(paste0(
        path, filenames[i]
      ))))
    }
  }

  # Return hyperspec object
  return(hs.x)
}
