#' Plot a range of wavelengths and their intervals
#'
#' This function plots the intervals between wavelengths to assess whether
#' or not they were constant.
#' @param trimmedwl a numeric vector of wavelengths (trimmed to the biological relevant part)
#' @importFrom graphics plot axis
#' @examples
#' data("mdqs")
#'
#' # Extract wavelengths and plot
#' mdqs.trim <- trim(mdqs, range=c(600, 1800))
#' wavelengths.trim <-  mass(mdqs.trim[[1]])
#' intervalplot(trimmedwl=wavelengths.trim)
#' @export
intervalplot <- function(trimmedwl){
  interval <- vector()
  name.interval <- vector()
  for (i in 1:length(trimmedwl)-1){
    interval[i] <- trimmedwl[i+1]-trimmedwl[i]
    name.interval[i] <- paste(round(trimmedwl[i], digits = 4),
                              round(trimmedwl[i+1], digits=4), sep=" - ")
  }
  plot(interval, pch = 20, cex = 0.8, bty="l",
       main = "Interval between wavenumbers",
       ylab=expression("Size interval [cm"^-1*"]"),
       xlab=expression("Range wavenumbers [cm"^-1*"]"),
       xaxt="n")
  axis(at = c(1,length(trimmedwl)),
       side = 1, labels = c(name.interval[1],name.interval[length(trimmedwl)]),
       las = 0)
}
