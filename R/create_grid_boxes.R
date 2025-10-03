#' Create Normalized Grid Boxes for Colonies
#'
#' Creates square bounding boxes for each colony based on their size (area) or
#' based on inter-colony distances. The function supports two modes:
#' - "inside": Creates the largest square that fits within each colony (assuming circular)
#' - "outside": Creates boxes based on median inter-colony distance with optional gaps
#'
#' If there is no colony, the function will label the new coordinates as NA.
#'
#' @param colony_data `data.frame`. Colony data containing at minimum columns
#'   `x`, `y`, `size`, `row`, and `col`. The `size` column should represent
#'   the area of each colony.
#' @param mode `character`. Either "inside" (default) or "outside".
#'   - "inside": Boxes fit within colony boundaries
#'   - "outside": Boxes based on inter-colony spacing
#' @param gap_factor `numeric`. Only used when mode = "outside".
#'   Fraction of median distance to use as gap between boxes (default: 0.1)
#'
#' @return A `data.frame` containing the original colony data plus new columns:
#'   \itemize{
#'     \item `new_xl`, `new_xr`: Left and right x-coordinates of bounding box
#'     \item `new_yt`, `new_yb`: Top and bottom y-coordinates of bounding box
#'     \item `grid_box_size`: Area of the created square bounding box
#'     \item `box_mode`: The mode used to create the boxes
#'   }
#' @export
create_grid_boxes <- function(colony_data, mode = "inside", gap_factor = 0.05) {

  # Validate inputs
  if (!mode %in% c("inside", "outside")) {
    stop("mode must be either 'inside' or 'outside'")
  }

  if (gap_factor < 0 || gap_factor >= 1) {
    stop("gap_factor must be between 0 and 1 (exclusive)")
  }

  if (mode == "inside") {

    colony_data$new_xl <- NA
    colony_data$new_xr <- NA
    colony_data$new_yt <- NA
    colony_data$new_yb <- NA
    colony_data$grid_box_size <- NA

    # Original functionality: boxes fit within colony boundaries
    for(i in 1:nrow(colony_data)) {
      # Get the size (area) for this specific colony
      colony_size <- colony_data$size[i]

      # Skip if size is NA or invalid
      if(is.na(colony_size) || colony_size <= 0) {
        next
      }

      # Back-calculate radius assuming colony size is area of a circle
      colony_radius <- sqrt(colony_size / pi)

      # Find largest square that fits within this circle
      # For a square inscribed in a circle: diagonal of square = 2 * radius
      # side_length = (2 * radius) / √2 = radius * √2
      square_side_length <- colony_radius * sqrt(2)
      half_side <- square_side_length / 2

      # Create bounding box centered on colony's x, y coordinates
      colony_data$new_xl[i] <- colony_data$x[i] - half_side
      colony_data$new_xr[i] <- colony_data$x[i] + half_side
      colony_data$new_yt[i] <- colony_data$y[i] - half_side
      colony_data$new_yb[i] <- colony_data$y[i] + half_side
      colony_data$grid_box_size[i] <- square_side_length^2
    }

  } else if (mode == "outside") {
      # Calculate grid spacing based on colony center positions
      rows <- sort(unique(colony_data$row))
      cols <- sort(unique(colony_data$col))

      # Calculate average spacing between colonies
      x_spacing <- mean(diff(sort(unique(colony_data$x))))
      y_spacing <- mean(diff(sort(unique(colony_data$y))))

      # Use the smaller spacing to ensure boxes don't overlap
      box_size <- min(x_spacing, y_spacing) * 0.95  # 95% to leave small gap
      half_box <- box_size / 2

      # Create new bounding boxes centered on colony positions
      colony_data$new_xl <- colony_data$x - half_box
      colony_data$new_xr <- colony_data$x + half_box
      colony_data$new_yt <- colony_data$y - half_box
      colony_data$new_yb <- colony_data$y + half_box
  }

  return(colony_data)
}
