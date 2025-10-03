#' Calculate Edge Intensity Difference Between Two Images
#'
#' Compares the mean pixel intensities in the edge regions of two images.
#' Extracts pixels within a specified distance from the image borders and
#' calculates the difference in their mean grayscale intensities. Useful
#' for background normalization and image comparison tasks.
#'
#' @param img1 `cimg`. First image object from the imager package.
#' @param img2 `cimg`. Second image object from the imager package.
#' @param edge_width `numeric(1)`. Width of edge region in pixels from the
#'   outer border. Default is 150.
#'
#' @return A `list` containing:
#'   \itemize{
#'     \item `mean_intensity1`: Mean intensity of edge pixels in first image
#'     \item `mean_intensity2`: Mean intensity of edge pixels in second image
#'     \item `difference`: Difference between the two mean intensities
#'   }
#' @export
calculate_edge_intensity <- function(img1, img2, edge_width = 150) {

  # Helper function to extract edge pixels and calculate median intensity
  calculate_edge_median <- function(img, width) {
    # Convert to grayscale if needed
    if (dim(img)[4] > 1) {
      gray_img <- grayscale(img)
    } else {
      gray_img <- img
    }

    # Get image dimensions
    img_width <- width(gray_img)
    img_height <- height(gray_img)

    # Create edge mask - TRUE for edge pixels, FALSE for interior
    edge_mask <- imfill(img_width, img_height, val = FALSE)

    # Left edge
    edge_mask[1:width, , , ] <- TRUE
    # Right edge
    edge_mask[(img_width - width + 1):img_width, , , ] <- TRUE
    # Top edge
    edge_mask[, 1:width, , ] <- TRUE
    # Bottom edge
    edge_mask[, (img_height - width + 1):img_height, , ] <- TRUE

    # Extract edge pixels
    edge_pixels <- gray_img[edge_mask]

    return(median(edge_pixels, na.rm = TRUE))
  }

  # Validate inputs
  if (missing(img1) || missing(img2)) {
    stop("Both img1 and img2 must be provided")
  }

  if (!inherits(img1, "cimg") || !inherits(img2, "cimg")) {
    stop("Both img1 and img2 must be cimg objects from the imager package")
  }

  # Get dimensions for validation
  dims1 <- c(width(img1), height(img1))
  dims2 <- c(width(img2), height(img2))

  if (!identical(dims1, dims2)) {
    warning("Images have different dimensions")
  }

  # Check if edge width is reasonable
  min_dim <- min(dims1)
  max_edge_width <- floor(min_dim / 2) - 1

  if (edge_width > max_edge_width) {
    warning(paste("Edge width", edge_width, "is large for image size.",
                  "Consider using", max_edge_width, "or smaller."))
  }

  median_intensity1 <- calculate_edge_median(img1, edge_width)
  median_intensity2 <- calculate_edge_median(img2, edge_width)

  intensity_difference <- median_intensity1 - median_intensity2

  return(list(
    median_intensity1 = median_intensity1,
    median_intensity2 = median_intensity2,
    difference = intensity_difference
  ))
}
