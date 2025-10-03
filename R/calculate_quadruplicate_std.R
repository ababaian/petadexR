#' Calculate Standard Deviation for Quadruplicates After Background Subtraction
#'
#' Compares original versus deconvolved variability across quadruplicate groups.
#' Groups colonies into 2x2 quadruplicates and calculates statistics for both
#' original and deconvolved intensities after background subtraction.
#'
#' @param result `list`. Output from `remove_halo_bleed_wiener()`.
#' @param data_with_bg `data.frame`. Colony data with background intensities.
#'   Must have bg_median_intensity column.
#'
#' @return A `list` containing quadruplicate summary statistics, individual
#'   colony details, and overall statistics comparing original vs deconvolved.
#' @export
calculate_quadruplicate_std <- function(result, data_with_bg) {

  # Extract comparison data
  comparison <- result$comparison_data

  # Standardize column names to lowercase
  names(comparison) <- tolower(names(comparison))
  names(data_with_bg) <- tolower(names(data_with_bg))

  # Check if bg_median_intensity is already in comparison data
  if ("bg_median_intensity" %in% names(comparison)) {
    merged_data <- comparison
  } else {
    # Merge with background data to get bg_median_intensity
    merged_data <- merge(comparison,
                         data_with_bg[, c("row", "col", "bg_median_intensity")],
                         by = c("row", "col"),
                         all.x = TRUE)
  }

  # Remove rows with missing data
  merged_data <- merged_data[!is.na(merged_data$deconvolved_median) &
                               !is.na(merged_data$bg_median_intensity), ]

  # Calculate background-subtracted values (both original and deconvolved)
  merged_data$bg_subtracted_original <- merged_data$original_median - merged_data$bg_median_intensity
  merged_data$bg_subtracted_deconv <- merged_data$deconvolved_median - merged_data$bg_median_intensity

  # Define quadruplicate groups
  # Rows: groups of 2 (1-2, 3-4, 5-6, ..., 15-16)
  # Cols: groups of 2 (1-2, 3-4, 5-6, ..., 23-24)
  merged_data$row_group <- ceiling(merged_data$row / 2)
  merged_data$col_group <- ceiling(merged_data$col / 2)

  # Create unique quadruplicate ID
  merged_data$quad_id <- paste0("R", merged_data$row_group, "_C", merged_data$col_group)

  # Calculate statistics for ORIGINAL (pre-deconvolution) quadruplicates
  quad_summary_orig <- aggregate(
    bg_subtracted_original ~ quad_id + row_group + col_group,
    data = merged_data,
    FUN = function(x) {
      c(
        mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        n = sum(!is.na(x)),
        cv = if(mean(x, na.rm = TRUE) != 0) {
          sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE)) * 100
        } else NA
      )
    }
  )

  # Calculate statistics for DECONVOLVED quadruplicates
  quad_summary_deconv <- aggregate(
    bg_subtracted_deconv ~ quad_id + row_group + col_group,
    data = merged_data,
    FUN = function(x) {
      c(
        mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        n = sum(!is.na(x)),
        cv = if(mean(x, na.rm = TRUE) != 0) {
          sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE)) * 100
        } else NA
      )
    }
  )

  # Unpack the results and combine
  quad_summary <- data.frame(
    quad_id = quad_summary_deconv$quad_id,
    row_group = quad_summary_deconv$row_group,
    col_group = quad_summary_deconv$col_group,
    original_mean = quad_summary_orig$bg_subtracted_original[, "mean"],
    original_median = quad_summary_orig$bg_subtracted_original[, "median"],
    original_std_dev = quad_summary_orig$bg_subtracted_original[, "sd"],
    original_cv_percent = quad_summary_orig$bg_subtracted_original[, "cv"],
    deconv_mean = quad_summary_deconv$bg_subtracted_deconv[, "mean"],
    deconv_median = quad_summary_deconv$bg_subtracted_deconv[, "median"],
    deconv_std_dev = quad_summary_deconv$bg_subtracted_deconv[, "sd"],
    deconv_cv_percent = quad_summary_deconv$bg_subtracted_deconv[, "cv"],
    n_values = quad_summary_deconv$bg_subtracted_deconv[, "n"]
  )

  # Calculate improvement metrics
  quad_summary$std_dev_reduction <- quad_summary$original_std_dev - quad_summary$deconv_std_dev
  quad_summary$std_dev_reduction_pct <- (quad_summary$std_dev_reduction / quad_summary$original_std_dev) * 100
  quad_summary$cv_reduction <- quad_summary$original_cv_percent - quad_summary$deconv_cv_percent

  # Sort by row and column groups
  quad_summary <- quad_summary[order(quad_summary$row_group, quad_summary$col_group), ]

  # Add individual colony values for reference
  colony_details <- merged_data[, c("row", "col", "row_group", "col_group",
                                    "quad_id", "original_median", "deconvolved_median",
                                    "bg_median_intensity", "bg_subtracted_original",
                                    "bg_subtracted_deconv")]
  colony_details <- colony_details[order(colony_details$row, colony_details$col), ]

  # Calculate overall statistics
  overall_stats <- data.frame(
    metric = c("Original: Mean CV", "Original: Median CV",
               "Original: Mean Std Dev", "Original: Median Std Dev",
               "Deconvolved: Mean CV", "Deconvolved: Median CV",
               "Deconvolved: Mean Std Dev", "Deconvolved: Median Std Dev",
               "Mean CV Reduction", "Median CV Reduction",
               "Mean Std Dev Reduction", "Median Std Dev Reduction",
               "Mean Std Dev Reduction %", "Median Std Dev Reduction %"),
    value = c(
      mean(quad_summary$original_cv_percent, na.rm = TRUE),
      median(quad_summary$original_cv_percent, na.rm = TRUE),
      mean(quad_summary$original_std_dev, na.rm = TRUE),
      median(quad_summary$original_std_dev, na.rm = TRUE),
      mean(quad_summary$deconv_cv_percent, na.rm = TRUE),
      median(quad_summary$deconv_cv_percent, na.rm = TRUE),
      mean(quad_summary$deconv_std_dev, na.rm = TRUE),
      median(quad_summary$deconv_std_dev, na.rm = TRUE),
      mean(quad_summary$cv_reduction, na.rm = TRUE),
      median(quad_summary$cv_reduction, na.rm = TRUE),
      mean(quad_summary$std_dev_reduction, na.rm = TRUE),
      median(quad_summary$std_dev_reduction, na.rm = TRUE),
      mean(quad_summary$std_dev_reduction_pct, na.rm = TRUE),
      median(quad_summary$std_dev_reduction_pct, na.rm = TRUE)
    )
  )

  return(list(
    quadruplicate_summary = quad_summary,
    colony_details = colony_details,
    overall_stats = overall_stats
  ))
}
