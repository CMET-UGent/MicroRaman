#' Function for creating a common legend for 2 ggplot2 figures.
#'
#' @param ... ggplot2 objects which could have a common legend
#' @param ncol number of columns in the grid (defaults to the number of plots)
#' @param nrow number of rows in the grid (defaults to 1)
#' @param position character vector (of length 2) defining the position of the legend. Defaults to bottom and right.
#' @importFrom ggplot2 ggplotGrob theme
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid unit.c unit grid.newpage grid.draw
#' @export

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
