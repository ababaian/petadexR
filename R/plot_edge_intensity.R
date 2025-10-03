#' Create Edge Intensity Histograms with Optional Normalization (for cimg objects)
#'
#' Generates histograms comparing edge pixel intensities between two images
#' with optional background normalization. Works with cimg objects from the
#' imager package loaded with load_image(). Can either compare original images
#' directly or show the effect of normalization on the second image. Useful
#' for visualizing background correction and intensity distributions.
#'
#' @param img_cimg1 `cimg`. First image object (typically background)
#'   loaded with `load_image()` from the imager package.
#' @param img_cimg2 `cimg`. Second image object (typically washed)
#'   loaded with `load_image()` from the imager package.
#' @param edge_width `numeric(1)`. Width of edge region in pixels from the
#'   outer border. Default is 150.
#' @param img1_name `character(1)`. Display name for first image.
#'   Default is "Background".
#' @param img2_name `character(1)`. Display name for second image.
#'   Default is "Washed".
#' @param normalize_background `logical(1)`. If `FALSE`, compares original images.
#'   If `TRUE`, shows normalization effect on second image. Default is `FALSE`.
#'
#' @return A `list` containing:
#'   \itemize{
#'     \item `plot`: ggplot2 comparison histogram
#'     \item `statistics`: Statistical summaries of intensities
#'     \item `normalized_image`: Normalized image cimg object (only if normalize_background = TRUE)
#'   }
#' @export
plot_edge_intensity <- function(img_cimg1, img_cimg2, edge_width = 150,
                                img1_name = "Background", img2_name = "Washed",
                                normalize_background = FALSE) {

  # Check if objects are cimg
  if (!inherits(img_cimg1, "cimg") || !inherits(img_cimg2, "cimg")) {
    stop("Both inputs must be cimg objects from the imager package")
  }

  # Helper function to extract edge pixels and return intensity data
  extract_edge_intensities <- function(img_cimg, width) {
    # Convert to grayscale using imager's built-in function
    gray_img <- imager::grayscale(img_cimg)

    # Get image dimensions
    img_width <- dim(gray_img)[1]
    img_height <- dim(gray_img)[2]

    # Create logical mask for edge pixels
    edge_mask <- array(FALSE, dim = c(img_width, img_height))

    # Mark edge regions
    edge_mask[1:width, ] <- TRUE                           # Left edge
    edge_mask[(img_width - width + 1):img_width, ] <- TRUE # Right edge
    edge_mask[, 1:width] <- TRUE                           # Bottom edge
    edge_mask[, (img_height - width + 1):img_height] <- TRUE # Top edge

    # Extract grayscale intensities for edge pixels
    intensities <- as.vector(gray_img[,,1,1])[edge_mask]

    return(intensities)
  }

  # Extract edge pixel intensities
  intensities1 <- extract_edge_intensities(img_cimg1, edge_width)
  intensities2 <- extract_edge_intensities(img_cimg2, edge_width)

  # Calculate statistics for original images
  stats1 <- list(
    mean = mean(intensities1, na.rm = TRUE),
    median = median(intensities1, na.rm = TRUE),
    sd = sd(intensities1, na.rm = TRUE),
    n = length(intensities1)
  )

  stats2 <- list(
    mean = mean(intensities2, na.rm = TRUE),
    median = median(intensities2, na.rm = TRUE),
    sd = sd(intensities2, na.rm = TRUE),
    n = length(intensities2)
  )

  if (!normalize_background) {
    # Compare background vs washed (before normalization)
    cat("=== BACKGROUND vs WASHED COMPARISON ===\n")
    cat(sprintf("%s Image: Mean=%.2f, Median=%.2f, SD=%.2f, n=%d\n",
                img1_name, stats1$mean, stats1$median, stats1$sd, stats1$n))
    cat(sprintf("%s Image: Mean=%.2f, Median=%.2f, SD=%.2f, n=%d\n",
                img2_name, stats2$mean, stats2$median, stats2$sd, stats2$n))
    cat(sprintf("Median difference: %.2f\n\n", stats1$median - stats2$median))

    # Create individual plots
    plot1 <- ggplot(data.frame(intensity = intensities1), aes(x = intensity)) +
      geom_histogram(bins = 30, fill = "lightblue", alpha = 0.7, color = "black") +
      geom_vline(xintercept = stats1$mean, color = "darkblue", linetype = "dashed", size = 1) +
      geom_vline(xintercept = stats1$median, color = "blue", linetype = "dotted", size = 1) +
      labs(title = paste(img1_name, "Image"),
           subtitle = paste("Mean:", round(stats1$mean, 2), "| Median:", round(stats1$median, 2),
                            "| SD:", round(stats1$sd, 2)),
           x = "Intensity (0-1)", y = "Pixel Count") +
      theme_minimal() +
      theme(plot.title = element_text(size = 12))

    plot2 <- ggplot(data.frame(intensity = intensities2), aes(x = intensity)) +
      geom_histogram(bins = 30, fill = "lightcoral", alpha = 0.7, color = "black") +
      geom_vline(xintercept = stats2$mean, color = "darkred", linetype = "dashed", size = 1) +
      geom_vline(xintercept = stats2$median, color = "red", linetype = "dotted", size = 1) +
      labs(title = paste(img2_name, "Image"),
           subtitle = paste("Mean:", round(stats2$mean, 2), "| Median:", round(stats2$median, 2),
                            "| SD:", round(stats2$sd, 2)),
           x = "Intensity (0-1)", y = "Pixel Count") +
      theme_minimal() +
      theme(plot.title = element_text(size = 12))

    comparison_plot <- grid.arrange(plot1, plot2, ncol = 2,
                                    top = paste("Background vs Washed Comparison | Median difference:",
                                                round(stats1$median - stats2$median, 2)))

    return(list(
      plot = comparison_plot,
      statistics = list(
        background_stats = stats1,
        washed_stats = stats2,
        difference = stats1$median - stats2$median
      )
    ))

  } else {
    # Perform normalization and compare original vs normalized washed

    # Calculate normalization factor
    normalization_factor <- stats1$median - stats2$median

    # Create normalized washed image
    normalized_img2 <- img_cimg2

    # Apply normalization factor to all channels
    if (dim(img_cimg2)[4] >= 3) {
      # Color image - normalize all RGB channels
      for (channel in 1:3) {
        channel_data <- normalized_img2[,,1,channel]

        # Convert to 0-255 if needed
        if (max(channel_data, na.rm = TRUE) <= 1) {
          channel_data <- channel_data * 255
          channel_data <- pmax(0, pmin(255, channel_data + normalization_factor))
          channel_data <- channel_data / 255  # Convert back to 0-1
        } else {
          channel_data <- pmax(0, pmin(255, channel_data + normalization_factor))
        }

        normalized_img2[,,1,channel] <- channel_data
      }
    } else {
      # Grayscale image
      channel_data <- normalized_img2[,,1,1]

      # Convert to 0-255 if needed
      if (max(channel_data, na.rm = TRUE) <= 1) {
        channel_data <- channel_data * 255
        channel_data <- pmax(0, pmin(255, channel_data + normalization_factor))
        channel_data <- channel_data / 255  # Convert back to 0-1
      } else {
        channel_data <- pmax(0, pmin(255, channel_data + normalization_factor))
      }

      normalized_img2[,,1,1] <- channel_data
    }

    # Extract normalized edge intensities
    normalized_intensities2 <- extract_edge_intensities(normalized_img2, edge_width)

    # Calculate normalized statistics
    normalized_stats2 <- list(
      mean = mean(normalized_intensities2, na.rm = TRUE),
      median = median(normalized_intensities2, na.rm = TRUE),
      sd = sd(normalized_intensities2, na.rm = TRUE),
      n = length(normalized_intensities2)
    )

    # Create individual plots
    plot1 <- ggplot(data.frame(intensity = intensities2), aes(x = intensity)) +
      geom_histogram(bins = 30, fill = "lightcoral", alpha = 0.7, color = "black") +
      geom_vline(xintercept = stats2$mean, color = "darkred", linetype = "dashed", size = 1) +
      geom_vline(xintercept = stats2$median, color = "red", linetype = "dotted", size = 1) +
      labs(title = paste("Original", img2_name),
           subtitle = paste("Mean:", round(stats2$mean, 2), "| Median:", round(stats2$median, 2),
                            "| SD:", round(stats2$sd, 2)),
           x = "Intensity (0-1)", y = "Pixel Count") +
      theme_minimal() +
      theme(plot.title = element_text(size = 12))

    plot2 <- ggplot(data.frame(intensity = normalized_intensities2), aes(x = intensity)) +
      geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7, color = "black") +
      geom_vline(xintercept = normalized_stats2$mean, color = "darkgreen", linetype = "dashed", size = 1) +
      geom_vline(xintercept = normalized_stats2$median, color = "green", linetype = "dotted", size = 1) +
      labs(title = paste("Normalized", img2_name),
           subtitle = paste("Mean:", round(normalized_stats2$mean, 2), "| Median:", round(normalized_stats2$median, 2),
                            "| SD:", round(normalized_stats2$sd, 2)),
           x = "Intensity (0-1)", y = "Pixel Count") +
      theme_minimal() +
      theme(plot.title = element_text(size = 12))

    comparison_plot <- grid.arrange(plot1, plot2, ncol = 2,
                                    top = paste("Washed Image: Original vs Normalized | Shift:",
                                                round(normalization_factor, 2)))

    return(list(
      plot = comparison_plot,
      normalized_image = normalized_img2,
      statistics = list(
        original_washed_stats = stats2,
        normalized_washed_stats = normalized_stats2,
        normalization_factor = normalization_factor
      )
    ))
  }
}
