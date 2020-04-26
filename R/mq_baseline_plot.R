#' Plot the result of optimal baseline correction
#'
#' This function plots a selected baseline-corrected Raman spectrum.
#' Optionally, the original spectrum can be added (for reference). At the moment
#' the function assumes that the cells are in the same order in the baseline-
#' corrected and original data.
#'
#' @param mscorr a list of MALDIquant::massSpectrum objects after baseline
#'     correction
#' @param msorig an optional list of MALDIquant::massSpectrum objects before
#'     baseline correction (defaults to NULL)
#' @param i an integer (defaults to 1) selecting the cell number. Cannot be
#'     larger than the number of cells in the mscorr list
#' @param addorig a boolean indicating wether or not to add the spectrum before
#'     baseline correction to the plot. If the list containing this data is not
#'     supplied, an error will be thrown (in the future this should be a warning).
#'     Defaults to FALSE
#' @param methodused a string indicating which method was used for baseline
#'     correction. Only used in the plot title. Defaults to "SNIP".
#' @param ... additional parameters passed on to plot
#' @importFrom MALDIquant isMassSpectrum intensity
#' @importFrom graphics lines legend plot
#' @export

mq_baseline_plot <- function(mscorr,msorig=NULL,i=1L,addorig=FALSE,methodused="SNIP",...){
  # TODO: allow for a vector of i to be passed
  if(i>length(mscorr)){
    stop(paste0("You selected cell number ",i," but there are only ",
                length(mscorr)," cells in the dataset"))
  }
  if(!MALDIquant::isMassSpectrum(mscorr[[i]])){
    stop("The baseline corrected data is not a valid MassSpectrum object")
  }
  if(!is.null(msorig)){
    if(!MALDIquant::isMassSpectrum(msorig[[i]])){
      stop("The original data you supplied is not a valid MassSpectrum object")
    }
  }
  if(addorig){
    if(is.null(msorig)){
      stop("You specified to overlay with the uncorrected spectrum, but did not supply it")
    }else{
      if(length(msorig)!=length(mscorr)){
        stop("The specified original and baseline corrected dataset do not have the same length")
      }else{
        plot(mscorr[[i]], main = paste0(methodused," baseline correction"),
             xlab=expression("Wavenumber (cm"^-1*")") ,
             ylab="Intensity (AU)",
             ylim=c(min(intensity(mscorr[[i]])), max(intensity(msorig[[i]]))),...)
        lines(x=msorig[[1]]@mass,y=msorig[[1]]@intensity,lty=3,col="red")
        legend("right",
               legend=c("Original spectrum","Baseline corrected spectrum"),
               lwd=2,lty=c(3,1),col=c("red","black"),bty="n")
      }
    }
  }else{
    plot(mscorr[[i]], main = paste0(methodused," baseline correction"),
         xlab=expression("Wavenumber (cm"^-1*")") ,
         ylab="Intensity (AU)",
         ylim=c(min(intensity(mscorr[[i]])), max(intensity(mscorr[[i]]))),...)
  }
}
