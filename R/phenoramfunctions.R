# Function for creating a common legend for 2 ggplot2 figures.
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)

}

### Functions for extracting the predicted R squared from lm models.
model_fit_stats <- function(linear.model) {
  r.sqr <- summary(linear.model)$r.squared
  adj.r.sqr <- summary(linear.model)$adj.r.squared
  pre.r.sqr <- pred_r_squared(linear.model)
  PRESS <- PRESS(linear.model)
  return.df <- data.frame(r.squared = r.sqr, adj.r.squared = adj.r.sqr, pred.r.squared = pre.r.sqr, press = PRESS)
  return(return.df)
}

pred_r_squared <- function(linear.model) {
  #' Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #' Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)

  return(pred.r.squared)
}

PRESS <- function(linear.model) {
  #' calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #' calculate the PRESS
  PRESS <- sum(pr^2)

  return(PRESS)
}

##### Normalization #######

# Better rounding function than R's base round
myround <- function(x) { trunc(x + 0.5) }


# Scales reads by
# 1) taking proportions
# 2) multiplying by a given library size of n
# 3) rounding
# Default for n is the minimum sample size in your library
# Default for round is floor
scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {

  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq,
                                          function(x) {(n * x/sum(x))}
  )

  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- myround(otu_table(physeq.scale))
  }

  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

# Extract Legend from ggplot2 object
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Function to pool FCS files based on sample name patterns

FCS_pool <- function(x, stub){
  if(length(stub) == length(x)) cat("-- No samples to merge --")
  for(i in 1:length(stub)){
    index <- grep(stub[i], flowCore::sampleNames(x))
    temp <- flowCore::flowSet(as(x[index], "flowFrame"))
    flowCore::sampleNames(temp) <- as.character(stub[i])
    if(stub[i] == stub[1]){
      concat_x <- temp
    } else {
      concat_x <- flowCore::rbind2(concat_x, temp)
    }
  }
  return(concat_x)
}

# Calculate correlations between bins and OTUs

FPcorSeq <- function(fp, phy, cor_thresh = 0.5, fp_thresh = 1e-10,
                     d = 10, param = c("FL1-H", "FL3-H"),
                     cor_m = "pearson", npoint = 10){

  # Make sure the phyloseq and flowbasis object are in same order
  phy@otu_table <- phy@otu_table[match(rownames(fp@basis),
                                       phyloseq::sample_names(phy)), ]
  cat(paste0("\t", "Check if sample names match between fbasis and phyloseq objects:",
             "\n"))
  cat(paste(rownames(fp@basis),"-", rownames(phy@otu_table), "\n"))
  # Normalize fingerprint
  fp@basis <- fp@basis/apply(fp@basis, 1, max)

  # Give bins coordinates
  nbin <- fp@nbin
  Y <- c()
  for (i in 1:nbin) Y <- c(Y, rep(i, nbin))
  X <- rep(1:nbin, nbin)

  # Position of data in @basis
  npos <- seq(1:nrow(fp@param))[fp@param[, 1] == param[1] & fp@param[, 2] ==
                                  param[2]]
  region <- ((npos - 1) * nbin * nbin + 1):(npos * nbin * nbin)

  fp@basis <- fp@basis[, region]

  # Additional filtering to reduce number of bins
  X <- X[colSums(fp@basis) > fp_thresh]
  Y <- Y[colSums(fp@basis) > fp_thresh]
  fp@basis <- fp@basis[, colSums(fp@basis) > fp_thresh]

  cat(paste0("\t",sum(colSums(fp@basis) > fp_thresh)," bins used to calculate correlation with OTUs: ", min(region), " - ", max(region),
             "\n"))

  # Calculate correlations
  cat(date(), paste0("---- Returning correlations >",
                     cor_thresh, "\n"))
  for(i in 1:ncol(fp@basis)){
    if(i%%100 == 0) cat(date(), paste0("---- at bin ", i, "/",  ncol(fp@basis), "\n"))
    cor_temp <- c(); p_temp <- c(); npoint_temp <- c()
    for(j in 1:dim(phy@otu_table)[2]){
      cor_temp[j] <- cor(fp@basis[,i], phy@otu_table[,j], method = cor_m)
      p_temp[j] <- cor.test(fp@basis[,i], phy@otu_table[,j], method = cor_m)$p.value
      npoint_temp[j] <- sum((phy@otu_table[, j] > 0)&(round(fp@basis[,i], d) > fp_thresh))
    }
    results_m <- data.frame(cor = cor_temp, cor.pval = p_temp, taxa = taxa_names(phy),
                            X = X[i], Y = Y[i], npoint = npoint_temp)

    if(i == 1) results_cat <- results_m
    results_cat <- rbind(results_cat, results_m)
  }

  ### Filter out low correlation values
  ### and correlation values calculated with low number of points
  results_cat <- results_cat[abs(results_cat$cor) > cor_thresh, ]
  colnames(results_cat)[4:5] <- param
  cat(date(), paste0("---- Done!\n"))
  return(results_cat)
}

# Example usage:
# ram_contrast(hyprs = hs.norm, comp1 = c("LB rep1", "LB rep2", "LB rep3"),
# comp2 = c("NB rep1","NB rep2","NB rep3"))

ram_contrast <- function(hyprs, comp1, comp2, plot = TRUE){
  # Extract spectra
  x <- hyprs@data$spc

  # Remove clusters column
  Samples <- factor(unlist(hyprs@label))
  Samples <- Samples[Samples!="clusters"]; Samples <- droplevels(Samples)
  lev1 <- comp1
  lev2 <- comp2
  comp1 <- Samples %in% comp1
  comp2 <- Samples %in% comp2

  # Print table of samples
  cat(paste0("-----------------------------------------------------------------------------------------------------",
             "\n \n"))
  cat(paste0("\t Your cells are distributed over these samples:\n\n "))
  print(table(Samples))
  cat(paste0("-----------------------------------------------------------------------------------------------------",
             "\n \n"))
  max.total <- max(x)
  if (length(comp1) == 1 & length(comp2) == 1)
    tmp <- (x[comp1, ] - x[comp2, ])
  if (length(comp1) == 1 & length(comp2) != 1)
    tmp <- x[comp1, ] - colMeans(x[comp2, ])
  if (length(comp1) != 1 & length(comp2) == 1)
    tmp <- colMeans(x[comp1, ]) - x[comp2, ]
  if (length(comp1) != 1 & length(comp2) != 1)
    tmp <- colMeans(x[comp1, ]) - colMeans(x[comp2,])
  tmp <- tmp/max.total
  results <- data.frame(Density = tmp, Wavenumber = hyprs@wavelength)
  cat(paste0("\t Returning contrasts between mean spectra for ", sum(comp1),
             " cells of\n ", list(lev1)))
  cat(paste0("\n\t and ", sum(comp2) ," cells of\n ",
             list(lev2), "\n"))
  cat(paste0("-----------------------------------------------------------------------------------------------------",
             "\n \n"))
  if(plot){
    # Plot
    v <- ggplot2::ggplot(results, ggplot2::aes(x = Wavenumber, y = Density, fill = Density))+
      ggplot2::geom_point(shape = 21, colour="gray", alpha = 0.4,
                          size = 4)+
      geom_line()+
      ggplot2::scale_fill_distiller(palette="RdBu", na.value="white") +
      ggplot2::theme_bw()
    print(v)
  }
  return(results)
}

# This function takes as input a hyperspec or maldiquant object and
# clusters the cells using partition around medioids (PAM) optimized through
# the silhouette index and model-based clustering using MClust optimized
# through the BIC criterion.
ram_clust <- function(sp.x, PCA = FALSE, ram_object = c("hs", "mq"),
                      PCA.var = 0.9, nclust = 50){

  # Get cell labels
  cell.labels <- rownames(sp.x)

  # Special formatting for maldiquant data
  if(ram_object == "mq"){
    wv_mq <- mass(sp.x[[1]])
    cell.labels <- paste(do.call(rbind, lapply(mq.norm, FUN = function(x) metaData(x)$name)),
                         seq(1:length(mq.norm)), sep = ".")
    matrix.spectra <- matrix(nrow=length(mq.norm), ncol = length(wv_mq))
    for (i in 1:length(mq.norm)){
      matrix.spectra[i,] <- intensity(sp.x[[i]])
    }
    sp.x <- new("hyperSpec", spc = matrix.spectra, wavelength = wv_mq, labels = cell.labels)
    rownames(sp.x) <- cell.labels
  }

  # Do PCA if requested
  if(PCA == TRUE){
    # Perform PCA to reduce number of features in fingerprint
    hs.x <- prcomp(sp.x)
    # Only retain PC which explain 90% of the variance
    nr_pc <- min(which((cumsum(vegan::eigenvals(hs.x)/sum(vegan::eigenvals(hs.x)))>PCA.var) == TRUE))
    pc_cluster <- hs.x$x[, 1:nr_pc]
  } else {
    pc_cluster <- sp.x
  }
  # Start PAM clustering
  # Evaluate number of robust clusters by means of silhouette index
  # We limit the search to "nclust" clusters
  tmp.si <- c()
  for(i in 2:nclust){
    if(i%%10 == 0) cat(date(), paste0("---- at k =  ", i, "/",  nrow(pc_cluster), "\n"))
    tmp.si[i] <- cluster::pam(pc_cluster, k = i)$silinfo$avg.width
  }
  nr_clusters_bacteria <- which(tmp.si == max(tmp.si, na.rm = TRUE))

  # Plot Silhouette index distribution
  plot(tmp.si, type = "l", ylab = "Silhouette index",
       xlab = "Number of clusters")
  lines(x = c(nr_clusters_bacteria, nr_clusters_bacteria),
        y = c(0,100), lty = 2, col = "red")

  # Cluster samples and export cluster labels
  clusters_x <- cluster::pam(pc_cluster, k = nr_clusters_bacteria)

  # Extract cluster labels
  cluster_labels_pam <- data.frame(Sample_label = rownames(sp.x),
                                   OPU = clusters_x$clustering)

  # Start MClust
  mc_fit <- Mclust(as.matrix(pc_cluster))
  summary(mc_fit) # display the best model

  cluster_labels_mc <- data.frame(Sample_label =  rownames(sp.x),
                                  OPU = mc_fit$classification)

  # Merge cluster outputs in long format df
  OPU_mq_merged <- rbind(cluster_labels_pam, cluster_labels_mc)
  OPU_mq_merged <- data.frame(OPU_mq_merged, method =
                                c(rep("PAM", nrow(cluster_labels_pam)),
                                  rep("Mclust", nrow(cluster_labels_mc)))
  )
  return(OPU_mq_merged)
}
