#' Similarity matrix based upon the Spectral contrast Angle (SCA) between all cells
#'
#' E. coli 2092 cultured in LB and NB, in biological triplicates at 28°C and
#' 120 rpm on an orbital shaker.
#' After 24h, cells are measured using FCM and fixed with PFA.
#' They are measured that day (day1, all samples).
#' Here we check this biological variability and try to identify how many
#' phenotypes were induced. As a means of calculating the distance between two
#' spectra, a distance matrix of all pairwise comparisons was generated using 
#' the SCA.diss function.
#' @format a dist object containing the SCA between each cell in the dataset
"SCA.similarity"