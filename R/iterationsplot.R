#' Plot the effect of iterations on baseline correction
#'
#' This function plots baseline correction for several iterations on the
#' spectral data to determine the optimal number of iterations
#'
#' @param mdq a MALDIquant::MassSpectrum (trimmed to relevant wavelengths)
#' @param itervect a numeric vector of different iteration options (max. 11 values)
#' @param plotlegend a boolean indicating wether or not to include a legend in lower left corner
#'     (defaults to TRUE)
#' @param ... additional parameters passed on to plot or lines
#' @importFrom RColorBrewer brewer.pal
#' @importFrom hyperSpec plot
#' @importFrom MALDIquant intensity estimateBaseline
#' @importFrom graphics lines legend
#' @examples
#' ## Short example
#'
#' # Load hyperSpec object
#' data("mdqs")
#'
#' # Extract wavelengths and plot
#' mdqs.trim <- trim(mdqs, range=c(600, 1800))
#' selms <- mdqs.trim[[1]]
#' iteration.options <- c(5,10,20,30,40,50,100)
#' iterationsplot(mdq=selms,itervect=iteration.options)
#' # if you want to check for all
#' # pdf("iterationplots.pdf")
#' # lapply(mdqs,iterationsplot,iteration.options)
#' # dev.off()
#' @export


iterationsplot <- function(mdq,itervect,plotlegend=TRUE,...){
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
       ylim=c(min(intensity(mdq)), max(intensity(mdq))))
  for (i in 1:length(itervect)){
    baseline <- estimateBaseline(mdq, method="SNIP",
                                 iterations=itervect[i])
    lines(baseline,col=colors[i],lwd=2)
  }
  if(plotlegend){
    legend("bottomleft",
           title="Number of iterations",
           legend = itervect, lwd=2,
           lty=c(1,1,1,1,1,1,1),
           col=colors, ncol=3, bty ="o")
  }
  }
