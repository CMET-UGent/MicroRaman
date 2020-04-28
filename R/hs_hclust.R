#' Perform hierarchical clustering on Raman spectra
#'
#' @param hs.x hyperSpec object
#' @param dist_method Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark",
#' "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn",
#' "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis". Can be specified to "SCA" for Raman data.
#' @param clust_method Choose regular hierarchical clustering or bootstrap supported hierarchical clustering with pvclust
#' @param nboot Numnber of bootstraps for pvclust. Defaults to 1000.
#' @param method.hclust Clustering method to use for pvclust. Either hclust or pvclust. Defaults to "pvclust"..
#' @param aggl_method the agglomerative method used in hierarchical clustering. This should be (an abbreviation of) one of "average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid". The default is "average".
#' @param ... Parameters to pass on to hclust().
#' @importFrom stats hclust dist
#' @importFrom pvclust pvrect pvclust
#' @importFrom vegan vegdist
#' @examples
#' ## Short example
#'
#' # Load hyperSpec object
#' data("hs_example")
#'
#' # Preprocess spectra
#' hs.x.proc <- hs_preprocess(hs_example)
#' hs.x.proc <- hs_resample(hs.x.proc, sample = 10)
#'
#' # Cluster
#' hclust_obj <- hs_hclust(hs.x.proc, dist_method = "manhattan",
#' clust_method = "pvclust")
#' @export

hs_hclust <- function(hs.x,
  dist_method = "bray",
  clust_method = "pvclust",
  aggl_method = "ward.D2",
  nboot = 1000,
  ...) {
  if (clust_method == "pvclust") {
    if (dist_method == "SCA") {
      dist_method <- function(x) {
        x <- t(x)
        diss <- matrix(nrow = nrow(x), ncol = nrow(x))
        for (i in 1:nrow(x)) {
          for (j in 1:nrow(x)) {
            diss[i, j] <- SCA(x[i,], x[j,])
          }
        }
        diss <- as.dist(diss)
        attr(diss, "method") <- "SCA.raman"
        return(diss)
      }
    }
    hs.hclust <-
      pvclust(
        data = t(hs.x@data$spc),
        method.dist = dist_method,
        method.hclust = aggl_method,
        nboot = nboot
      )
  } else {
    hs.dist <- vegan::vegdist(hs.x@data$spc, method = dist_method)
    hs.hclust <- hclust(d = hs.dist, aggl_method, ...)
  }
  return(hs.hclust)
}
