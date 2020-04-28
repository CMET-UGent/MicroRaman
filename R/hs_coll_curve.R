#' Collection curves for phenotypic heterogeneity measurements
#'
#' Calculates collection curves that describe the sensitivity of the diversity metrics to sample size.
#' This can be used to check whether the Raman spectra in their totality cover the phenotypic heterogeneity of the sample.
#'
#' @param hs.x HyperSpec object
#' @param intervals Decide in which sample size intervals you want to inspect the curves.
#' Defaults to 10 cells.
#' @param nboot Number of bootstraps to condut for each sample size. Defaults to 10.
#' @param peak_detection Should peak detection be used instead of raw spectra for Hill calculations? Defaults to FALSE
#' @param peak_window Peak windows size for a peak to be considered a signal. Defaults to 20.
#' @param peak_method Which peak detection method should be used (requires peak_detection == TRUE).
#' Options are "MAD" and "SuperSmoother"
#' @param plot_fig Should figures of collectors curves be shown? Defaults to TRUE.
#' @importFrom dplyr group_by summarize
#' @importFrom magrittr %>%
#' @importFrom ggplot2 geom_point theme_bw geom_smooth geom_errorbar labs theme ylim ggplot aes element_text
#' @importFrom cowplot plot_grid
#' @importFrom rlang .data
#' @importFrom stats sd
#' @examples
#' # Short example
#' data("hs_example")
#'
#' # Preprocess
#' hs_example <- hs_preprocess(hs_example)
#'
#' # Calculate metrics
#' hs_coll_curve(hs_example, intervals = 3, nboot = 10, plot_fig = FALSE)
#' @export
#'
hs_coll_curve <- function(hs.x,
  intervals = 10,
  nboot = 10,
  peak_detection = FALSE,
  peak_window = 20,
  peak_method = c("MAD"),
  plot_fig = TRUE){

  # Create sequence of cells to check
  seq_int <- c(1, seq(from = intervals, to = length(hs.x), by = intervals))

  # Loop through bootstraps and sample sizes
  for(i in seq_int){
    for(j in 1:nboot){
      hs.x.boot <- hs_resample(hs.x, sample = i, replace = TRUE)
      hs.x.boot.ram <- hs_phenoRam(hs.x.boot,
        peak_detection = peak_detection,
        peak_window = peak_window,
        peak_method = peak_method,
        return_spectra = FALSE)
      # add colmeans
      if (j == 1 & i == seq_int[1])
        tmp <- data.frame(nboot = j, nspec = i, hs.x.boot.ram)
      else
        tmp <- rbind(tmp, data.frame(nboot = j, nspec = i, hs.x.boot.ram))
      #
    }
  }

  # Calc means and st dev for each cell size
  results_coll <- tmp %>%
    group_by(.data$nspec) %>%
    summarize(D0.mean = mean(.data$D0), D1.mean = mean(.data$D1),
              D2.mean = mean(.data$D2),
              D0.sd = sd(.data$D0), D1.sd = sd(.data$D1),
              D2.sd = sd(.data$D2))

  # Plot the results
  if(plot_fig){
    p1 <- ggplot(results_coll, aes(x = .data$nspec, y = .data$D0.mean))+
      geom_point(size = 3, shape = 21)+
      theme_bw()+
      geom_smooth(se = FALSE, col = "black")+
      geom_errorbar(aes(ymin = .data$D0.mean - .data$D0.sd,
                        ymax = .data$D0.mean + .data$D0.sd), width = 0.02)+
      labs(y = "D0", x = "number of spectra")+
      theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
      ylim(0, 1.1*max(results_coll$D0.mean))

    p2 <- ggplot(results_coll, aes(x = .data$nspec, y = .data$D1.mean))+
      geom_point(size = 3, shape = 21)+
      theme_bw()+
      geom_smooth(se = FALSE, col = "black")+
      geom_errorbar(aes(ymin = .data$D1.mean - .data$D1.sd,
                        ymax = .data$D1.mean + .data$D1.sd), width = 0.02)+
      labs(y = "D1", x = "number of spectra")+
      theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
      ylim(0, 1.1*max(results_coll$D0.mean))

    p3 <- ggplot(results_coll, aes(x = .data$nspec, y = .data$D2.mean))+
      geom_point(size = 3, shape = 21)+
      theme_bw()+
      geom_smooth(se = FALSE, col = "black")+
      geom_errorbar(aes(ymin = .data$D2.mean - .data$D2.sd,
                        ymax = .data$D2.mean + .data$D2.sd), width = 0.02)+
      labs(y = "D2", x = "number of spectra")+
      theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
      ylim(0, 1.1*max(results_coll$D0.mean))

    print(plot_grid(p1, p2, p3, nrow = 1))
  }

  return(results_coll)
}
