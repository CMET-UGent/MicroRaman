#' Perform clustering of spectra using partitioning around medoids (PAM)
#'
#' Number of clusters is optimized based on Silhouette index.
#' @param hs.x hyperSpec object
#' @param PCA.var Variance threshold for PCA preprocessing prior to clustering.
#' Defaults to NULL which means that no PCA will be performed.
#' @param ... Parameters to pass on to pam().
#' @importFrom cluster pam
#' @importFrom stats prcomp
#' @examples
#' ## Short example
#'
#' # Load hyperSpec object
#' data("hs.example")
#'
#' # Convert to MassSpectrum object
#' hs.x.proc <- hs_preprocess(hs.x)
#'
#' hclust_obj <- hs_type(hs.x.proc, PCA.var = 0.8)
#' @export

hs_type <- function(hs.x, PCA.var = NULL, ...) {

  # Perform PCA to reduce number of features in fingerprint
  if(!is.null(PCA.var)){
    hs_pca <- prcomp(x = hs.x@data$spc)
    # Only retain PC which explain PCA.var% of the variance
    nr_pc <-
      min(which((
        cumsum(vegan::eigenvals(hs_pca) / sum(vegan::eigenvals(hs_pca))) >
          PCA.var
      ) == TRUE))
    pc_cluster <- hs_pca$x[, 1:nr_pc]
    if (nr_pc < 2) {
      stop(
        "Error: only one PC remained - consider increasing PCA.var"
      )
    }
  } else pc_cluster <- hs.x@data$spc


  # Evaluate number of robust clusters by means of silhouette index
  tmp.si <- c()
  for (i in 2:(nrow(pc_cluster) - 1)) {
    tmp.si[i] <- pam(pc_cluster, k = i)$silinfo$avg.width
  }
  nr_clusters <- which(tmp.si == max(tmp.si, na.rm = TRUE))

  # Cluster samples and export cluster labels
  hs.type <- pam(pc_cluster, k = nr_clusters)

  # Extract cluster labels
  if(!is.null(rownames(hs.x@data$spc))) sample_id <- rownames(hs.x@data$spc) else {
    sample_id <- hs.x@data$filename
  }

  # Merge results
  hs.type <- data.frame(Spectrum_ID = sample_id,
    cluster_label = hs.type$clustering)

  return(hs.type)
}
