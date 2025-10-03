#' Plot Point Spread Function
#'
#' Creates visualizations of the PSF including a 2D heatmap and radial profile,
#' along with diagnostic metrics for quality assessment.
#'
#' @param psf `matrix`. Point spread function to visualize.
#'
#' @return A `list` containing ggplot objects (heatmap, radial profile) and a
#'   data frame of diagnostic metrics.
#' @export
#'
plot_psf <- function(psf) {
  library(ggplot2)
  library(reshape2)

  psf_df <- melt(psf)
  colnames(psf_df) <- c("x", "y", "value")

  p1 <- ggplot(psf_df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    theme_minimal() +
    labs(title = "Point Spread Function (2D)", fill = "Intensity")

  # Radial profile
  center_x <- ceiling(nrow(psf) / 2)
  center_y <- ceiling(ncol(psf) / 2)

  radial_data <- do.call(rbind, lapply(1:nrow(psf), function(i) {
    lapply(1:ncol(psf), function(j) {
      data.frame(
        distance = sqrt((i - center_x)^2 + (j - center_y)^2),
        intensity = psf[i, j]
      )
    })
  }))
  radial_data <- do.call(rbind, radial_data)
  radial_data$dist_bin <- round(radial_data$distance)
  radial_summary <- aggregate(intensity ~ dist_bin, radial_data, mean)

  p2 <- ggplot(radial_summary, aes(x = dist_bin, y = intensity)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 2) +
    theme_minimal() +
    labs(title = "PSF Radial Profile",
         x = "Distance from Center (pixels)",
         y = "Mean Intensity",
         subtitle = "Should decay smoothly from center")

  diagnostics <- data.frame(
    Metric = c("Peak Value", "Center Value", "Sum (should be ~1)",
               "FWHM (pixels)", "Symmetry Score"),
    Value = c(
      sprintf("%.4f", max(psf)),
      sprintf("%.4f", psf[center_x, center_y]),
      sprintf("%.4f", sum(psf)),
      sprintf("%.1f", calculate_fwhm(radial_summary)),
      sprintf("%.3f", assess_psf_symmetry(psf))
    )
  )

  return(list(heatmap = p1, radial = p2, diagnostics = diagnostics))
}

calculate_fwhm <- function(radial_summary) {
  max_int <- max(radial_summary$intensity)
  half_max <- max_int / 2
  above_half <- radial_summary[radial_summary$intensity >= half_max, ]

  if (nrow(above_half) > 0) {
    return(max(above_half$dist_bin) * 2)
  }
  return(NA)
}

assess_psf_symmetry <- function(psf) {
  center_x <- ceiling(nrow(psf) / 2)
  center_y <- ceiling(ncol(psf) / 2)

  left_sum <- sum(psf[, 1:center_y])
  right_sum <- sum(psf[, center_y:ncol(psf)])
  lr_ratio <- min(left_sum, right_sum) / max(left_sum, right_sum)

  top_sum <- sum(psf[1:center_x, ])
  bottom_sum <- sum(psf[center_x:nrow(psf), ])
  tb_ratio <- min(top_sum, bottom_sum) / max(top_sum, bottom_sum)

  return((lr_ratio + tb_ratio) / 2)
}

assess_psf_quality <- function(psf) {
  center_x <- ceiling(nrow(psf) / 2)
  center_y <- ceiling(ncol(psf) / 2)

  issues <- c()
  score <- 100

  # Check peak location
  max_idx <- which(psf == max(psf), arr.ind = TRUE)[1, ]
  dist_from_center <- sqrt((max_idx[1] - center_x)^2 + (max_idx[2] - center_y)^2)

  if (dist_from_center > 3) {
    issues <- c(issues, sprintf("Peak is %.1f pixels from center (should be <3)",
                                dist_from_center))
    score <- score - 20
  } else {
    issues <- c(issues, "Peak is centered")
  }

  # Check normalization
  psf_sum <- sum(psf)
  if (abs(psf_sum - 1) > 0.01) {
    issues <- c(issues, sprintf("PSF sum is %.4f (should be ~1.0)", psf_sum))
    score <- score - 10
  } else {
    issues <- c(issues, "PSF is properly normalized")
  }

  # Check symmetry
  symmetry <- assess_psf_symmetry(psf)
  if (symmetry < 0.8) {
    issues <- c(issues, sprintf("PSF asymmetry detected (score: %.2f, should be >0.8)",
                                symmetry))
    score <- score - 20
  } else {
    issues <- c(issues, sprintf("PSF is symmetric (score: %.2f)", symmetry))
  }

  # Check smooth decay
  radial_data <- do.call(rbind, lapply(1:nrow(psf), function(i) {
    lapply(1:ncol(psf), function(j) {
      data.frame(
        distance = sqrt((i - center_x)^2 + (j - center_y)^2),
        intensity = psf[i, j]
      )
    })
  }))
  radial_data <- do.call(rbind, radial_data)
  radial_data$dist_bin <- round(radial_data$distance)
  radial_summary <- aggregate(intensity ~ dist_bin, radial_data, mean)
  radial_summary <- radial_summary[order(radial_summary$dist_bin), ]

  diffs <- diff(radial_summary$intensity)
  increasing_count <- sum(diffs > 0)
  if (increasing_count > length(diffs) * 0.2) {
    issues <- c(issues, sprintf("PSF doesn't decay smoothly (%d increases)",
                                increasing_count))
    score <- score - 15
  } else {
    issues <- c(issues, "PSF decays smoothly from center")
  }

  # Check FWHM
  fwhm <- calculate_fwhm(radial_summary)
  if (!is.na(fwhm)) {
    if (fwhm < 2 || fwhm > nrow(psf) * 0.8) {
      issues <- c(issues, sprintf("FWHM is %.1f pixels (unusual)", fwhm))
      score <- score - 15
    } else {
      issues <- c(issues, sprintf("FWHM is %.1f pixels (reasonable)", fwhm))
    }
  }

  quality <- if (score >= 90) {
    "EXCELLENT"
  } else if (score >= 70) {
    "GOOD"
  } else if (score >= 50) {
    "FAIR"
  } else {
    "POOR"
  }

  return(list(score = score, quality = quality, issues = issues))
}
