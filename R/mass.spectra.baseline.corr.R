#' Baseline-corrected Raman spectrometry test single-cell data of Escherichia coli LMG2092 biological replicates
#'
#' E. coli 2092 cultured in LB and NB, in biological triplicates at 28°C and 
#' 120 rpm on an orbital shaker. 
#' After 24h, cells are measured using FCM and fixed with PFA. 
#' They are measured that day (day1, all samples). 
#' In this case, the spectra were trimmed to a biologically relevant region and
#' subsequently baseline correction using MALDIquant::estimateBaseline was 
#' executed with method="SNIP".
#' @format A list of MALDIquant::massSpectrum objects after baseline correction
"mass.spectra.baseline.corr"