#' Calculate the spectral conrast angle (SCA) and generate an iTol-compatible tree
#'
#'  The SCA measures the angle between two vectors corresponding to closely
#'  related spectra to measure whether they are the same or not
#'  (Wan *et al.*, 2002). This function returns an iTol (Letunic & Bork, 2016)
#'
#' @param x a hyperSpec object containing the cells in the rows, and wavenumbers (features) in columns
#' @param Factor1 a character vector containing treatment levels (defaults to different media)
#' @param Factor1Name a string naming Factor 1 (defaults to "Medium")
#' @param Factor2 a character vector containing "treatment" levels (defaults to different replicates)
#' @param Factor2Name a string naming Factor 2 (defautls to "Replicate")
#' @param clusMethod a string showing which stats::hclust linkage method should be chosen
#' @importFrom stats as.dist hclust
#' @importFrom ape as.phylo write.tree
#' @return an ape::phylo object based upon hclust of SCA distances of the hyperSpec object's cells
#' @references Wan, K. X., Vidavsky, I., & Gross, M. L. (2002). Comparing similar spectra: from similarity index to spectral contrast angle. Journal of the American Society for Mass Spectrometry, 13(1), 85-88.
#' @references Letunic, I., & Bork, P. (2016). Interactive tree of life (iTOL) v3: an online tool for the display and annotation of phylogenetic and other trees. Nucleic acids research, 44(W1), W242-W245.
#' @export

hs_SCA_conv_itol <- function(x,Factor1=c("LB","NB","LB"),Factor1Name="Medium",
                      Factor2=c("rep1","rep2","rep3"),Factor2Name="Replicate",
                      clusMethod="ward.D2"){
  if(class(x)!="hyperSpec"){
    stop("You did not supply a valid hyperSpec object, and there is no default.")
  }
  if(length(Factor1)!=length(x)){
    stop(paste0("The length of Factor 1 (",
                Factor1Name,
                ") is not equal to that of the number of observations in the ",
                "hyperSpec object : ",length(Factor1)," in stead of ",
                length(x)))}
  if(length(Factor2)!=length(x)){
    stop(paste0("The length of Factor 2 (",
                Factor2Name,
                ") is not equal to that of the number of observations in the ",
                "hyperSpec object :",length(Factor2)," in stead of ",
                length(x)))}
  cnvecttmp <- paste(Factor1,Factor2,collapse=NULL)
  uniqcnvecctmp <- make.unique(cnvecttmp,sep="_")
  rownames(x@data$spc) <- uniqcnvecctmp
  hyperSpec::rownames(x) <- uniqcnvecctmp
  simitol <- SCA.diss(x)
  denditol <- hclust(simitol,method="ward.D2")
  denditol$tip.label <- uniqcnvecctmp
  dendrogram.itol <- as.phylo(denditol)
  return(dendrogram.itol)
}
