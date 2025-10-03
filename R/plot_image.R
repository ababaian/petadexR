#' Plot Image with Colony Overlays
#'
#' Creates a visualization of an image with optional overlays showing colony
#' grids, center points, and row-column labels. Uses base R graphics for
#' memory-efficient plotting of large images.
#'
#' @param img `cimg`. Image object loaded by `imager::load.image()`.
#' @param colony_data `data.frame` or `NULL`. Colony data with coordinates and
#'   bounding box information. If `NULL`, only the base image is plotted.
#' @param show_colony_grid `logical(1)`. If `TRUE`, displays colony bounding
#'   boxes as red rectangles. Default is `FALSE`.
#' @param use_new_boxes `logical(1)`. If `TRUE`, uses normalized bounding boxes
#'   (new_xl, new_xr, new_yt, new_yb) from `create_grid_boxes()`. If `FALSE`,
#'   uses original boxes (xl, xr, yt, yb). Default is `FALSE`.
#' @param show_center_points `logical(1)`. If `TRUE`, displays colony center
#'   points as yellow dots. Default is `FALSE`.
#' @param show_rowcol_labels `logical(1)`. If `TRUE`, displays row-column
#'   labels for each colony. Default is `FALSE`.
#'
#' @return Invisibly returns `NULL`. Plots are created as side effects.
#' @export
plot_image <- function(img,
                       colony_data = NULL,
                       show_colony_grid = FALSE,
                       use_new_boxes = FALSE,
                       show_center_points = FALSE,
                       show_rowcol_labels = FALSE) {

  # Plot the base image using imager's plot method
  plot(img, axes = FALSE, main = NULL)

  # Add overlays if colony_data is provided
  if (!is.null(colony_data)) {

    # Add colony grid overlay
    if (show_colony_grid) {
      if (use_new_boxes) {
        rect(xleft = colony_data$new_xl,
             ybottom = colony_data$new_yb,
             xright = colony_data$new_xr,
             ytop = colony_data$new_yt,
             border = "red",
             lwd = 2)
      } else {
        rect(xleft = colony_data$xl,
             ybottom = colony_data$yb,
             xright = colony_data$xr,
             ytop = colony_data$yt,
             border = "red",
             lwd = 2)
      }
    }

    # Add center points
    if (show_center_points) {
      points(x = colony_data$x,
             y = colony_data$y,
             col = "yellow",
             pch = 16,
             cex = 1.2)
    }

    # Add row-column labels
    if (show_rowcol_labels) {
      text(x = colony_data$x,
           y = colony_data$y,
           labels = paste0(colony_data$row, "-", colony_data$col),
           col = "white",
           cex = 0.8,
           pos = 3,  # Position above the point
           font = 2)  # Bold font
    }
  }

  invisible(NULL)
}
