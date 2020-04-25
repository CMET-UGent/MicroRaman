#' Perform hierarchical clustering on Raman spectra
#'
#' @param hs.x hyperSpec object
#' @param dist_method Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark",
#' "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn",
#' "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis".
#' @param ... Parameters to pass on to hclust().
#' @importFrom stats hclust dist
#' @importFrom vegan vegdist
#' @examples
#' ## Short example
#'
#' # Load hyperSpec object
#' data("hs.example")
#'
#' # Convert to MassSpectrum object
#' hs.x.proc <- hs_preprocess(hs.x)
#'
#' hs_hclust <- hs_hclust(hs.x.proc)
#' @export

hs_hclust <- function(hs.x, dist_method = "bray", ...) {
  hs.dist <- vegan::vegdist(hs.x@data$spc, method = dist_method)
  hs.hclust <- hclust(d = hs.dist, ...)
  return(hs.hclust)
}
