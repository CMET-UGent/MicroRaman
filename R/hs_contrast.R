#' Contrast plot for Raman spectra
#'
#' A way to visualize contrasts between Raman spectra from different samples in an
#'  existing hyperSpec object
#'
#' @param hs.x a hyperSpec object
#' @param comp1 numeric (any length) or boolean (TRUE/FALSE - length equal to number of spectra in hs.x) vector
#'  with sample names from one side of the comparisons
#' @param comp2 numeric (any length) or boolean (TRUE/FALSE - length equal to number of spectra in hs.x) vector
#'  with sample names from the other side of the comparisons
#' @param plot a logical indicating whether or not the constrasts should be plotted (defaults to TRUE)
#' @return a dataframe containing the contrast values between the averaged spectra
#' @importFrom ggplot2 geom_line ggplot geom_point scale_fill_distiller theme_bw
#' @examples
#' ## Short example
#' data("hs_example")
#'
#' # Preprocess
#' hs_example <- hs_preprocess(hs_example)
#'
#' # Run contrast comparison
#' hs.cont <- hs_contrast(hs_example, comp1 = c(1:4), comp2 = c(8:20))
#' @export

hs_contrast <- function(hs.x, comp1, comp2, plot = TRUE) {
  if (is.null(hs.x) | class(hs.x) != "hyperSpec") {
    stop(
      "Error: you did not supply a valid hyperSpec object,
      and there is no default, please correct"
    )
  }

  # Extract spectra
  x <- hs.x@data$spc

  # Remove clusters column
  # Samples <- factor(unlist(hs.x@label))
  # Samples <-
  #   Samples[Samples != "clusters"]
  # Samples <- droplevels(Samples)
  # lev1 <- comp1
  # lev2 <- comp2

  # Make contrasts
  max.total <- max(x)
  if (length(comp1) == 1 & length(comp2) == 1)
    tmp <- (x[comp1, ] - x[comp2, ])
  if (length(comp1) == 1 & length(comp2) != 1)
    tmp <- x[comp1, ] - colMeans(x[comp2, ])
  if (length(comp1) != 1 & length(comp2) == 1)
    tmp <- colMeans(x[comp1, ]) - x[comp2, ]
  if (length(comp1) != 1 & length(comp2) != 1)
    tmp <- colMeans(x[comp1, ]) - colMeans(x[comp2,])

  # normalize for max intensity in comparison
  tmp <- tmp / max.total

  # Results to dataframe
  hs.contrast <-
    data.frame(Density = tmp, Wavenumber = hs.x@wavelength)

  # Sanity check and plot
  # cat(paste0(
  #   "\t Returning contrasts between mean spectra for ",
  #   length(comp1),
  #   " cells of\n ",
  #   rownames(hs.x@data$spc)[comp1]
  # ))
  # cat(paste0("\n\t and ", length(comp2) , " cells of\n ",
  #   rownames(hs.x@data$spc)[comp1], "\n"))
  # cat(
  #   paste0(
  #     "-----------------------------------------------------------------------------------------------------",
  #     "\n \n"
  #   )
  # )
  if (plot) {
    # Plot
    v <-
      ggplot2::ggplot(
        hs.contrast,
        ggplot2::aes(
          x = .data$Wavenumber,
          y = .data$Density,
          fill = .data$Density
        )
      ) +
      ggplot2::geom_point(
        shape = 21,
        colour = "black",
        alpha = 0.7,
        size = 3
      ) +
      geom_line() +
      ggplot2::scale_fill_distiller(palette = "RdBu", na.value = "white") +
      ggplot2::theme_bw()
    print(v)
  }
  return(hs.contrast)
  }
