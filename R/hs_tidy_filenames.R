#' hs_tidy_filenames function
#'
#' Function to tidy/clean filenames in hyperspec objects (i.e. remove absolute path in name, and/or replace filenames).
#' @param hs.x HyperSpec object
#' @param remove_pattern What should be adjusted in filenames.
#' "path" means remove the absolute path of the file and only keep the filename
#' "rename" means that you supply a dataframe with matching old and new filenames.
#' In this case you must supply a dataframe in rename_df and designate the old filename
#' column (old_name_col) and new filename column (new_name_col).
#' @param rename_df dataframe containing old and new filenames for the hyperspec objects
#' @param old_name_col column in rename_df with old filenames
#' @param new_name_col column in rename df with new filenames
#' @importFrom dplyr left_join
#' @keywords hyperspec rename
#' @examples
#' data("hs_example")
#'
#' # Remove path from spectral ID
#' hs_tidy_filenames(hs_example, remove_pattern = "path")
#'
#' @export

hs_tidy_filenames <- function(hs.x,
  remove_pattern = "path",
  rename_df = NULL,
  old_name_col,
  new_name_col) {
  if (remove_pattern == "path") {
    # Edit spectra IDs
    if (is.null(rownames(hs.x@data$spc))) {
      hs.x@data$filename <-
        gsub(pattern = ".*/",
          replacement = "",
          x = hs.x@data$filename)
    }  else {
      rownames(hs.x@data$spc) <-
        gsub(pattern = ".*/",
          replacement = "",
          x = rownames(hs.x@data$spc))
    }
  } else if (remove_pattern == "rename" && !is.null(rename_df)) {
    # Edit spectra IDs
    if (is.null(rownames(hs.x@data$spc))) {
      hs_order <- match(hs.x@data$filename, rename_df[, old_name_col])
      hs.x@data$filename <- rename_df[hs_order, new_name_col]
    }  else {
      hs_order <-
        match(rownames(hs.x@data$spc), rename_df[, old_name_col])
      rownames(hs.x@data$spc) <- rename_df[hs_order, new_name_col]
    }
  }
  return(hs.x)
}
