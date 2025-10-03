#' Classify Pixels by Colony Bounding Box
#'
#' Classifies pixels within each colony's bounding box as "clearing", "background",
#' or "halo" based on intensity thresholds derived from each individual box's background statistics.
#' Uses median ± 4 standard deviations as classification thresholds. Requires
#' colony data with background intensity statistics.
#'
#' @param img `cimg`. Image object loaded by `imager::load.image()`.
#' @param colony_data_with_bg `data.frame`. Colony data that includes background
#'   intensity statistics (bg_median_intensity, bg_sd_intensity) from
#'   `calculate_average_intensity()`.
#' @param background_normalization `numeric(1)`. Value to add to all pixel
#'   intensities for normalization purposes. Calculated by
#'   `edge_intensity_difference()$difference`.
#'
#' @return A `list` containing three elements:
#'   \itemize{
#'     \item `classified_data`: Image data frame with added `label` and `colony_id` columns
#'     \item `colony_summary`: Summary statistics for each colony including pixel counts and percentages
#'     \item `total_stats`: Overall statistics across all colonies
#'   }
#' @export
classify_pixels <- function(img, colony_data_with_bg, background_normalization) {

  # Convert imager object to data frame more efficiently
  img_df <- as.data.frame(img)

  # Extract RGB values efficiently using matrix operations
  rgb_matrix <- matrix(img_df$value, ncol = 3, byrow = FALSE)
  intensity_vals <- rgb_to_gray(rgb_matrix[,1], rgb_matrix[,2], rgb_matrix[,3]) + background_normalization

  # Create processed image data frame with pre-allocated vectors
  n_pixels <- nrow(img_df) / 3  # 3 color channels
  processed_img <- data.frame(
    x = img_df$x[img_df$cc == 1],
    y = img_df$y[img_df$cc == 1],
    intensity = intensity_vals,
    label = rep("outside", n_pixels),
    colony_id = rep(NA_character_, n_pixels),
    stringsAsFactors = FALSE
  )

  # Pre-allocate summary data frame with exact size
  n_colonies <- nrow(colony_data_with_bg)
  classification_summary <- data.frame(
    colony_id = character(n_colonies),
    row = character(n_colonies),
    col = character(n_colonies),
    clearing_count = integer(n_colonies),
    background_count = integer(n_colonies),
    halo_count = integer(n_colonies),
    total_pixels = integer(n_colonies),
    clearing_pct = numeric(n_colonies),
    background_pct = numeric(n_colonies),
    halo_pct = numeric(n_colonies),
    median_intensity = numeric(n_colonies),
    sd_intensity = numeric(n_colonies),
    clearing_median_intensity = numeric(n_colonies),
    background_median_intensity = numeric(n_colonies),
    halo_median_intensity = numeric(n_colonies),
    stringsAsFactors = FALSE
  )

  # Pre-calculate all thresholds to avoid repeated calculations
  valid_rows <- !is.na(colony_data_with_bg$new_xl) &
    !is.na(colony_data_with_bg$new_xr) &
    !is.na(colony_data_with_bg$new_yt) &
    !is.na(colony_data_with_bg$new_yb)

  bg_medians <- colony_data_with_bg$bg_median_intensity
  bg_sds <- colony_data_with_bg$bg_sd_intensity
  lower_thresholds <- bg_medians - 4 * bg_sds
  upper_thresholds <- bg_medians + 4 * bg_sds

  total_classified <- 0
  summary_idx <- 0


  # Vectorized processing of all colonies
  for(i in 1:n_colonies) {
    summary_idx <- summary_idx + 1

    # Set default values for this row
    colony_id <- paste0(colony_data_with_bg$row[i], "-", colony_data_with_bg$col[i])
    classification_summary$colony_id[summary_idx] <- colony_id
    classification_summary$row[summary_idx] <- colony_data_with_bg$row[i]
    classification_summary$col[summary_idx] <- colony_data_with_bg$col[i]
    classification_summary$median_intensity[summary_idx] <- bg_medians[i]
    classification_summary$sd_intensity[summary_idx] <- bg_sds[i]

    # Skip invalid bounding boxes
    if(!valid_rows[i]) {
      next
    }

    # Vectorized bounding box mask
    box_mask <- processed_img$x >= colony_data_with_bg$new_xl[i] &
      processed_img$x <= colony_data_with_bg$new_xr[i] &
      processed_img$y >= colony_data_with_bg$new_yt[i] &
      processed_img$y <= colony_data_with_bg$new_yb[i]

    box_pixels_count <- sum(box_mask)

    if(box_pixels_count == 0) {
      next
    }

    # Vectorized assignment and classification
    processed_img$colony_id[box_mask] <- colony_id

    box_intensities <- processed_img$intensity[box_mask]

    # Vectorized classification
    labels <- ifelse(box_intensities < lower_thresholds[i], 1L,  # clearing = 1
                     ifelse(box_intensities > upper_thresholds[i], 3L, 2L))  # halo = 3, background = 2

    processed_img$label[box_mask] <- c("clearing", "background", "halo")[labels]

    # Fast counting using tabulate
    label_counts <- tabulate(labels, nbins = 3)
    clearing_count <- label_counts[1]
    background_count <- label_counts[2]
    halo_count <- label_counts[3]

    # Update summary efficiently
    classification_summary$clearing_count[summary_idx] <- clearing_count
    classification_summary$background_count[summary_idx] <- background_count
    classification_summary$halo_count[summary_idx] <- halo_count
    classification_summary$total_pixels[summary_idx] <- box_pixels_count
    classification_summary$clearing_pct[summary_idx] <- round(clearing_count / box_pixels_count * 100, 1)
    classification_summary$background_pct[summary_idx] <- round(background_count / box_pixels_count * 100, 1)
    classification_summary$halo_pct[summary_idx] <- round(halo_count / box_pixels_count * 100, 1)

    # Vectorized median calculations
    if(clearing_count > 0) {
      classification_summary$clearing_median_intensity[summary_idx] <- median(box_intensities[labels == 1])
    } else {
      classification_summary$clearing_median_intensity[summary_idx] <- NA
    }

    if(background_count > 0) {
      classification_summary$background_median_intensity[summary_idx] <- median(box_intensities[labels == 2])
    } else {
      classification_summary$background_median_intensity[summary_idx] <- NA
    }

    if(halo_count > 0) {
      classification_summary$halo_median_intensity[summary_idx] <- median(box_intensities[labels == 3])
    } else {
      classification_summary$halo_median_intensity[summary_idx] <- NA
    }

    total_classified <- total_classified + box_pixels_count
  }

  # Convert labels to factor efficiently
  processed_img$label <- factor(processed_img$label,
                                levels = c("clearing", "background", "halo", "outside"))

  # Calculate overall statistics using vectorized operations
  overall_totals <- colSums(classification_summary[c("clearing_count", "background_count", "halo_count")], na.rm = TRUE)

  return(list(
    classified_data = processed_img,
    colony_summary = classification_summary,
    total_stats = list(
      total_pixels = total_classified,
      clearing_total = overall_totals["clearing_count"],
      background_total = overall_totals["background_count"],
      halo_total = overall_totals["halo_count"],
      overall_clearing_median_intensity = median(classification_summary$clearing_median_intensity, na.rm = TRUE),
      overall_background_median_intensity = median(classification_summary$background_median_intensity, na.rm = TRUE),
      overall_halo_median_intensity = median(classification_summary$halo_median_intensity, na.rm = TRUE)
    )
  ))
}


rgb_to_gray <- function(r, g, b) {
  # Standard luminance formula
  return(0.299 * r + 0.587 * g + 0.114 * b)
}
