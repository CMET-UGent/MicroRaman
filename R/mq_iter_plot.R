#' Plot the effect of iterations on baseline correction
#'
#' This function plots baseline correction for several iterations on the
#' spectral data to determine the optimal number of iterations
#'
#' @param mdq a MALDIquant::MassSpectrum (trimmed to relevant wavelengths)
#' @param itervect a numeric vector of different iteration options (max. 11 values)
#' @param plotlegend a boolean indicating wether or not to include a legend in lower left corner
#'     (defaults to TRUE)
#' @param method one of "SNIP", "TopHat", "ConvexHull" or "median" for baseline
#'     estimation by MALDIquant::estimateBaseline (see ?MALDIquant::estimateBaseline).
#'     defaults to "SNIP".
#' @param ... additional parameters passed on to plot or lines
#' @importFrom RColorBrewer brewer.pal
#' @importFrom hyperSpec plot
#' @importFrom MALDIquant intensity estimateBaseline
#' @importFrom graphics lines legend
#' @export


iterationsplot <- function(mdq,itervect,method="SNIP",plotlegend=TRUE,...){
  method <- match.arg(method,c("SNIP", "TopHat", "ConvexHull","median"))
  ncol <- length(itervect)
  if(ncol>11){
    stop("Currently the maximum number of different iteration options is 11.
         Please revise.")
  }
  if(class(mdq)!="MassSpectrum"){
    stop("You did not supply a valid MALDIquant::MassSpectrum object,
         please correct.")
  }
  if(ncol < 3){
    #Rcolorbrewer palettes accept a minimum of 3 colors only
    spectrvect <- brewer.pal(11, "Spectral")
    colors <- spectrvect[1:ncol]
  } else {
    colors <- brewer.pal(ncol, "Spectral")
  }
  smpnm <- mdq@metaData$name
  plot(mdq,
       main = paste0("SNIP baseline correction: influence of iterations for ",
                     smpnm),
       xlab=expression("Wavenumber (cm"^-1*")"), ylab="Intensity (AU)",
       ylim=c(min(intensity(mdq)), max(intensity(mdq))),...)
  for (i in 1:length(itervect)){
    baseline <- estimateBaseline(mdq, method=method,
                                 iterations=itervect[i])
    lines(baseline,col=colors[i],lwd=2,...)
  }
  if(plotlegend){
    legend("bottomleft",
           title="Number of iterations",
           legend = itervect, lwd=2,
           lty=c(1,1,1,1,1,1,1),
           col=colors, ncol=3, bty ="o")
  }
  }
