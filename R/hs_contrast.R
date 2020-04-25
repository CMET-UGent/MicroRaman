#' Contrast plot for Raman spectra
#'
#' A way to visualize contrasts between Raman spectra from different samples an existing hyperSpec object
#'
#' @param hs.x a hyperspec object
#' @param comp1 character vector with sample names from one side of the comparisons
#' @param comp2 character vector with sample names from the other side of the comparisons
#' @param plot a logical indicating whether or not the constrasts should be plotted (defaults to TRUE)
#' @return a ggplot2 object with the contrast plot
#' @importFrom dendextend is.hclust cutree
#' @importFrom stats rect.hclust
#' @importFrom rlang .data
#' @importFrom ggplot2 geom_line ggplot geom_point scale_fill_distiller theme_bw
#' @examples
#' ## Short example
#' ### Load and pretransform MassSpectrum object
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
#' ### Convert to hyperSpec object
#' hs.norm <- mq2hs(mq.norm,sampleNames=smpnms)
#' hs_contrast(hs.x = hs.norm, comp1 = c("LB rep1", "LB rep2", "LB rep3"),
#' comp2 = c("NB rep1","NB rep2","NB rep3"))
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
  Samples <- factor(unlist(hs.x@label))
  Samples <-
    Samples[Samples != "clusters"]
  Samples <- droplevels(Samples)
  lev1 <- comp1
  lev2 <- comp2
  comp1 <- Samples %in% comp1
  comp2 <- Samples %in% comp2

  # Print table of samples
  cat(
    paste0(
      "-----------------------------------------------------------------------------------------------------",
      "\n \n"
    )
  )
  cat(paste0("\t Your cells are distributed over these samples:\n\n "))
  print(table(Samples))
  cat(
    paste0(
      "-----------------------------------------------------------------------------------------------------",
      "\n \n"
    )
  )
  max.total <- max(x)
  if (length(comp1) == 1 & length(comp2) == 1)
    tmp <- (x[comp1,] - x[comp2,])
  if (length(comp1) == 1 & length(comp2) != 1)
    tmp <- x[comp1,] - colMeans(x[comp2,])
  if (length(comp1) != 1 & length(comp2) == 1)
    tmp <- colMeans(x[comp1,]) - x[comp2,]
  if (length(comp1) != 1 & length(comp2) != 1)
    tmp <- colMeans(x[comp1,]) - colMeans(x[comp2, ])
  tmp <- tmp / max.total
  results <- data.frame(Density = tmp, Wavenumber = hs.x@wavelength)
  cat(paste0(
    "\t Returning contrasts between mean spectra for ",
    sum(comp1),
    " cells of\n ",
    list(lev1)
  ))
  cat(paste0("\n\t and ", sum(comp2) , " cells of\n ",
    list(lev2), "\n"))
  cat(
    paste0(
      "-----------------------------------------------------------------------------------------------------",
      "\n \n"
    )
  )
  if (plot) {
    # Plot
    v <-
      ggplot2::ggplot(
        results,
        ggplot2::aes(
          x = .data$Wavenumber,
          y = .data$Density,
          fill = .data$Density
        )
      ) +
      ggplot2::geom_point(
        shape = 21,
        colour = "gray",
        alpha = 0.4,
        size = 4
      ) +
      geom_line() +
      ggplot2::scale_fill_distiller(palette = "RdBu", na.value = "white") +
      ggplot2::theme_bw()
    print(v)
  }
  return(results)
  }
