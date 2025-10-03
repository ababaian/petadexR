#' Initialize Colony Data from File
#'
#' Reads colony data from a tab-separated values (TSV) file, automatically
#' skipping comment lines that start with "#" and properly formatting the
#' first column as "row". The function handles files with variable numbers
#' of header comment lines.
#'
#' The coordinates in the input file are assumed to already be in the top-left
#' origin coordinate system, so no coordinate transformation is performed.
#'
#' @param path_name `character(1)`. Path to the TSV file containing colony data.
#'   The file may contain comment lines starting with "#" at the beginning.
#' @param image_obj `cimg`. Image object from the imager package (currently
#'   unused but included for potential future functionality).
#'
#' @return A `data.frame` containing the colony data from the file with the
#'   first column renamed to "row" and comment lines properly skipped.
#' @export
initialize_colony_data <- function(path_name, image_obj) {
  con <- file(path_name, "r")
  lines_to_skip <- 0
  while(TRUE) {
    line <- readLines(con, n = 1)
    if(length(line) == 0) break
    if(startsWith(line, "#")) {
      lines_to_skip <- lines_to_skip + 1
    } else {
      break
    }
  }
  close(con)
  data <- read_tsv(path_name, skip = lines_to_skip - 1, show_col_types = FALSE)
  colnames(data)[1] <- "row"

  # No coordinate flipping needed - coordinates are already in top-left origin system
  return(data)
}
