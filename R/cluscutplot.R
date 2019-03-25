#' Plot clusters on hclust object at given height.
#'
#' An easy exploratory way to draw clusters at a given height on a dendrogram.
#'
#' @param dend a hclust object
#' @param h cutoff distance for the clusters to be plotted (defaults to 0.75)
#' @param bordercol the color of the border
#' @param ... other parameters passed on to plot.hclust
#' @return a matrix with the assigned clusters for each leaf (in the case of Raman data each cell)
#' @importFrom dendextend is.hclust cutree
#' @importFrom stats rect.hclust
#' @examples
#' ## Short example
#' data("SCA.similarity")
#' dendrogram <- hclust(SCA.similarity, method="ward.D2")
#' cluscutplot(dendrogram)
#' @export

cluscutplot <- function(dend,h=0.75,bordercol="red",...){
  if(!is.hclust(dend)){stop("dend is not a valid hclust object")}
  clusters <- as.matrix(cutree(dend,h=h))
  plot(dend,...)
  rect.hclust(dend,k=max(clusters),border=bordercol)
  return(clusters)
}
