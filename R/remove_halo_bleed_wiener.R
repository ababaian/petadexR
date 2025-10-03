# Estimate Point Spread Function from Peak Intensity Centering
#
# Identifies highly active colonies and extracts PSF by centering on peak
# intensities within each colony region. Accumulates and averages PSF patterns
# from multiple colonies to create a robust estimate of the optical blur.
#
# @param img `cimg`. Washed image object from the imager package.
# @param colony_data `data.frame`. Colony data with bounding box coordinates.
# @param background_img `cimg`. Background image for activity calculation.
# @param activity_percentile `numeric(1)`. Percentile threshold for selecting
#   highly active colonies. Default is 0.9.
# @param psf_radius `numeric(1)`. Radius of the PSF in pixels. Default is 15.
# @param invert_for_clearing `logical(1)`. Whether to invert intensities for
#   clearing assays. Default is TRUE.
#
# @return A normalized `matrix` representing the estimated point spread function.
# @export

estimate_psf_from_peaks <- function(img, colony_data, background_img,
                                    activity_percentile = 0.9,
                                    psf_radius = 15,
                                    invert_for_clearing = TRUE) {

  coords <- if (all(c("new_xl", "new_yt", "new_xr", "new_yb") %in% names(colony_data))) {
    c("new_xl", "new_yt", "new_xr", "new_yb")
  } else {
    c("xl", "yt", "xr", "yb")
  }

  # Calculate activities for each colony
  colony_data$activity <- sapply(1:nrow(colony_data), function(i) {
    colony <- colony_data[i, ]
    colony_coords <- c(colony[[coords[1]]], colony[[coords[2]]],
                       colony[[coords[3]]], colony[[coords[4]]])

    if (any(is.na(colony_coords))) return(NA_real_)
    if (colony_coords[1] >= colony_coords[3] || colony_coords[2] >= colony_coords[4]) {
      return(NA_real_)
    }

    tryCatch({
      region <- img[colony_coords[1]:colony_coords[3],
                    colony_coords[2]:colony_coords[4], 1, 1]
      bg_region <- background_img[colony_coords[1]:colony_coords[3],
                                  colony_coords[2]:colony_coords[4], 1, 1]

      if (invert_for_clearing) {
        activity <- mean(bg_region - region, na.rm = TRUE)
      } else {
        activity <- mean(region - bg_region, na.rm = TRUE)
      }

      return(activity)
    }, error = function(e) NA_real_)
  })

  # Find highly active colonies
  valid_activities <- colony_data$activity[!is.na(colony_data$activity)]

  if (length(valid_activities) == 0) {
    stop("Cannot estimate PSF without valid colonies")
  }

  threshold <- quantile(valid_activities, activity_percentile, na.rm = TRUE)
  highly_active <- colony_data[!is.na(colony_data$activity) &
                                 colony_data$activity > threshold, ]

  if (nrow(highly_active) == 0) {
    return(create_gaussian_psf(psf_radius, sigma = psf_radius / 3))
  }

  # Build PSF by finding peaks and extracting centered regions
  psf_size <- 2 * psf_radius + 1
  psf_accumulator <- matrix(0, psf_size, psf_size)
  psf_count <- 0

  for (i in 1:nrow(highly_active)) {
    colony <- highly_active[i, ]
    colony_coords <- c(colony[[coords[1]]], colony[[coords[2]]],
                       colony[[coords[3]]], colony[[coords[4]]])

    if (any(is.na(colony_coords))) next

    tryCatch({
      # Get colony region
      region <- img[colony_coords[1]:colony_coords[3],
                    colony_coords[2]:colony_coords[4], 1, 1]
      bg_region <- background_img[colony_coords[1]:colony_coords[3],
                                  colony_coords[2]:colony_coords[4], 1, 1]

      # Process based on assay type
      if (invert_for_clearing) {
        processed <- bg_region - region
      } else {
        processed <- region - bg_region
      }

      # Find the peak location within this colony
      peak_idx <- which(processed == max(processed), arr.ind = TRUE)[1, ]

      # Convert to image coordinates
      peak_x <- colony_coords[1] + peak_idx[1] - 1
      peak_y <- colony_coords[2] + peak_idx[2] - 1

      # Extract PSF region centered on this peak
      x_start <- peak_x - psf_radius
      x_end <- peak_x + psf_radius
      y_start <- peak_y - psf_radius
      y_end <- peak_y + psf_radius

      # Check bounds
      if (x_start < 1 || x_end > width(img) ||
          y_start < 1 || y_end > height(img)) {
        next
      }

      # Extract PSF from raw image
      psf_region <- img[x_start:x_end, y_start:y_end, 1, 1]

      if (invert_for_clearing) {
        psf_processed <- max(psf_region) - psf_region
      } else {
        psf_processed <- psf_region
      }

      # Normalize
      psf_processed <- psf_processed - min(psf_processed)

      if (max(psf_processed) > 0) {
        psf_processed <- psf_processed / max(psf_processed)

        # Verify the peak is centered
        center_idx <- ceiling(psf_size / 2)
        peak_check <- which(psf_processed == max(psf_processed), arr.ind = TRUE)[1, ]
        peak_offset <- sqrt((peak_check[1] - center_idx)^2 + (peak_check[2] - center_idx)^2)

        if (peak_offset < 3) {
          psf_accumulator <- psf_accumulator + psf_processed
          psf_count <- psf_count + 1
        }
      }
    }, error = function(e) {})
  }

  if (psf_count == 0) {
    return(create_gaussian_psf(psf_radius, sigma = psf_radius / 3))
  }

  # Average the accumulated PSFs
  psf <- psf_accumulator / psf_count

  # Final normalization
  psf <- psf / sum(psf)

  return(psf)
}

# Create Gaussian Point Spread Function
#
# Generates a 2D Gaussian PSF as a fallback when empirical PSF estimation fails
# or when a simple analytical PSF is preferred.
#
# @param radius `numeric(1)`. Radius of the PSF in pixels.
# @param sigma `numeric(1)`. Standard deviation of the Gaussian distribution.
#   Default is radius/3.
#
# @return A normalized `matrix` representing a Gaussian point spread function.
# @export

create_gaussian_psf <- function(radius, sigma = NULL) {
  if (is.null(sigma)) sigma <- radius / 3

  size <- 2 * radius + 1
  center <- radius + 1

  psf <- matrix(0, size, size)
  for (i in 1:size) {
    for (j in 1:size) {
      dist_sq <- (i - center)^2 + (j - center)^2
      psf[i, j] <- exp(-dist_sq / (2 * sigma^2))
    }
  }

  return(psf / sum(psf))
}

# Wiener Deconvolution in Fourier Domain
#
# Applies Wiener deconvolution to remove optical blur from an image using
# the estimated PSF. Operates in the Fourier domain for computational efficiency
# and applies optimal filtering based on noise characteristics.
#
# @param img `cimg`. Image to be deconvolved.
# @param psf `matrix`. Point spread function estimated from the image.
# @param noise_power `numeric(1)`. Estimated noise power for Wiener filter
#   regularization. Default is 0.01.
#
# @return A `cimg` object containing the deconvolved image.
# @export

wiener_deconvolve <- function(img, psf, noise_power = 0.01) {

  img_array <- as.array(img)[,,1,1]

  # Store original dimensions
  orig_rows <- nrow(img_array)
  orig_cols <- ncol(img_array)

  # Pad PSF to image size
  psf_padded <- matrix(0, orig_rows, orig_cols)
  psf_center_x <- floor(nrow(psf) / 2) + 1
  psf_center_y <- floor(ncol(psf) / 2) + 1

  for (i in 1:nrow(psf)) {
    for (j in 1:ncol(psf)) {
      x_idx <- (i - psf_center_x + orig_rows) %% orig_rows + 1
      y_idx <- (j - psf_center_y + orig_cols) %% orig_cols + 1
      psf_padded[x_idx, y_idx] <- psf[i, j]
    }
  }

  # Apply Wiener filter in Fourier domain
  img_fft <- fft(img_array)
  psf_fft <- fft(psf_padded)

  psf_conj <- Conj(psf_fft)
  psf_power <- Mod(psf_fft)^2

  wiener_filter <- psf_conj / (psf_power + noise_power)
  result_fft <- img_fft * wiener_filter
  result <- Re(fft(result_fft, inverse = TRUE)) / length(result_fft)
  result <- pmax(result, 0)

  return(as.cimg(result))
}

#' Remove Halo Bleed Using Wiener Deconvolution
#'
#' Complete pipeline for removing optical blur (halo bleed) from colony images
#' using Wiener deconvolution. Estimates PSF from highly active colonies,
#' performs background subtraction, applies deconvolution, and generates
#' comparison statistics between original and deconvolved intensities.
#'
#' @param img `cimg`. Washed image object from the imager package.
#' @param background_img `cimg`. Background image for subtraction.
#' @param colony_data `data.frame`. Colony data with bounding box coordinates.
#' @param psf_radius `numeric(1)`. Radius of the PSF in pixels. Default is 15.
#' @param noise_power `numeric(1)`. Noise power for Wiener filter. Default is 0.01.
#' @param activity_percentile `numeric(1)`. Percentile for selecting active colonies.
#'   Default is 0.9.
#' @param test_mode `logical(1)`. Whether to process only first 4x4 colonies for
#'   testing. Default is FALSE.
#' @param invert_for_clearing `logical(1)`. Whether to invert intensities for
#'   clearing assays. Default is TRUE.
#'
#' @return A `list` containing deconvolved images, PSF, comparison statistics,
#'   and processing settings.
#' @export

remove_halo_bleed_wiener <- function(img, background_img, colony_data,
                                     psf_radius = 15,
                                     noise_power = 0.01,
                                     activity_percentile = 0.9,
                                     test_mode = FALSE,
                                     invert_for_clearing = TRUE) {

  coords <- if (all(c("new_xl", "new_yt", "new_xr", "new_yb") %in% names(colony_data))) {
    c("new_xl", "new_yt", "new_xr", "new_yb")
  } else {
    c("xl", "yt", "xr", "yb")
  }

  # Handle test mode (first 4x4 colonies)
  if (test_mode) {
    test_colonies <- colony_data[colony_data$row <= 4 & colony_data$col <= 4, ]

    x_min <- max(1, min(test_colonies[[coords[1]]], na.rm = TRUE) - 50)
    x_max <- min(width(img), max(test_colonies[[coords[3]]], na.rm = TRUE) + 50)
    y_min <- max(1, min(test_colonies[[coords[2]]], na.rm = TRUE) - 50)
    y_max <- min(height(img), max(test_colonies[[coords[4]]], na.rm = TRUE) + 50)

    img_use <- as.cimg(img[x_min:x_max, y_min:y_max, 1, 1])
    bg_use <- as.cimg(background_img[x_min:x_max, y_min:y_max, 1, 1])
    colony_data_use <- test_colonies
  } else {
    img_use <- img
    bg_use <- background_img
    colony_data_use <- colony_data
    x_min <- y_min <- NULL
  }

  # Estimate PSF using peak-based centering
  psf <- estimate_psf_from_peaks(img_use, colony_data_use, bg_use,
                                 activity_percentile, psf_radius,
                                 invert_for_clearing)

  # Convert to arrays for processing
  img_array <- as.array(img_use)[,,1,1]
  bg_array <- as.array(bg_use)[,,1,1]

  if (invert_for_clearing) {
    img_bg_sub <- pmax(bg_array - img_array, 0)
  } else {
    img_bg_sub <- pmax(img_array - bg_array, 0)
  }

  img_bg_subtracted <- as.cimg(img_bg_sub)

  # Deconvolve
  deconvolved <- wiener_deconvolve(img_bg_subtracted, psf, noise_power)

  # Convert back to original intensity space if needed
  if (invert_for_clearing) {
    deconv_array <- as.array(deconvolved)[,,1,1]

    if (nrow(deconv_array) != nrow(bg_array) || ncol(deconv_array) != ncol(bg_array)) {
      min_rows <- min(nrow(deconv_array), nrow(bg_array))
      min_cols <- min(ncol(deconv_array), ncol(bg_array))
      deconv_array <- deconv_array[1:min_rows, 1:min_cols]
      bg_array <- bg_array[1:min_rows, 1:min_cols]
    }

    result_array <- bg_array - deconv_array
    deconvolved_original_space <- as.cimg(result_array)
  } else {
    deconvolved_original_space <- deconvolved
  }

  comparison_data <- data.frame(
    row = colony_data$row,
    col = colony_data$col,
    original_median = NA_real_,
    original_mean = NA_real_,
    deconvolved_median = NA_real_,
    deconvolved_mean = NA_real_,
    median_change = NA_real_,
    mean_change = NA_real_,
    median_change_percent = NA_real_,
    mean_change_percent = NA_real_
  )

  for (i in 1:nrow(colony_data)) {
    colony <- colony_data[i, ]

    if (test_mode && (colony$row > 4 || colony$col > 4)) next

    colony_coords <- c(colony[[coords[1]]], colony[[coords[2]]],
                       colony[[coords[3]]], colony[[coords[4]]])

    if (any(is.na(colony_coords))) next

    tryCatch({
      if (test_mode) {
        adj_coords <- pmax(colony_coords - c(x_min - 1, y_min - 1, x_min - 1, y_min - 1), 1)
        adj_coords[3] <- min(adj_coords[3], dim(img_use)[1])
        adj_coords[4] <- min(adj_coords[4], dim(img_use)[2])

        region_orig <- img_use[adj_coords[1]:adj_coords[3],
                               adj_coords[2]:adj_coords[4], 1, 1]
        region_deconv <- deconvolved_original_space[adj_coords[1]:adj_coords[3],
                                                    adj_coords[2]:adj_coords[4], 1, 1]
      } else {
        region_orig <- img_use[colony_coords[1]:colony_coords[3],
                               colony_coords[2]:colony_coords[4], 1, 1]
        region_deconv <- deconvolved_original_space[colony_coords[1]:colony_coords[3],
                                                    colony_coords[2]:colony_coords[4], 1, 1]
      }

      orig_mean <- mean(region_orig, na.rm = TRUE)
      orig_median <- median(region_orig, na.rm = TRUE)
      deconv_mean <- mean(region_deconv, na.rm = TRUE)
      deconv_median <- median(region_deconv, na.rm = TRUE)

      comparison_data$original_mean[i] <- orig_mean
      comparison_data$original_median[i] <- orig_median
      comparison_data$deconvolved_mean[i] <- deconv_mean
      comparison_data$deconvolved_median[i] <- deconv_median
      comparison_data$mean_change[i] <- deconv_mean - orig_mean
      comparison_data$median_change[i] <- deconv_median - orig_median

      if (orig_mean > 0) {
        comparison_data$mean_change_percent[i] <-
          ((deconv_mean - orig_mean) / orig_mean) * 100
      }
      if (orig_median > 0) {
        comparison_data$median_change_percent[i] <-
          ((deconv_median - orig_median) / orig_median) * 100
      }
    }, error = function(e) {})
  }

  return(list(
    deconvolved_img = deconvolved_original_space,
    deconvolved_inverted_space = deconvolved,
    original_img = img_use,
    original_bg_subtracted = img_bg_subtracted,
    psf = psf,
    comparison_data = comparison_data,
    test_mode = test_mode,
    test_region = if(test_mode) list(x_min=x_min, x_max=x_max,
                                     y_min=y_min, y_max=y_max) else NULL,
    settings = list(
      noise_power = noise_power,
      psf_radius = psf_radius,
      invert_for_clearing = invert_for_clearing,
      activity_percentile = activity_percentile
    )
  ))
}
