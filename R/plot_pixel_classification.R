#' Visualize Pixel Classification Results
#'
#' Creates a visualization showing the spatial distribution of pixel
#' classifications (clearing, background, halo) overlaid as colored tiles.
#' Filters out "outside" pixels for cleaner visualization and provides
#' summary information about the number of colonies analyzed.
#'
#' @param classification_result `list`. Result object from
#'   `classify_pixels_by_box()` containing `classified_data`, `colony_summary`,
#'   and `total_stats` elements.
#'
#' @return A `ggplot2` object showing pixel classifications as colored tiles:
#'   \itemize{
#'     \item Dark blue: clearing pixels
#'     \item Light gray: background pixels
#'     \item Orange: halo pixels
#'   }
#' @export
plot_pixel_classification <- function(classification_result) {

  classified_data <- classification_result$classified_data

  # Filter out "outside" pixels for cleaner visualization
  colony_pixels <- classified_data[classified_data$label != "outside", ]

  # Count total colonies analyzed
  total_colonies <- length(unique(colony_pixels$colony_id[!is.na(colony_pixels$colony_id)]))

  # Create the plot
  class_plot <- ggplot(colony_pixels, aes(x = x, y = y, fill = label)) +
    geom_tile() +
    scale_fill_manual(values = c("clearing" = "darkblue",
                                 "background" = "lightgray",
                                 "halo" = "orange"),
                      name = "Pixel Type") +
    scale_y_reverse() +  # Flip y-axis to match image orientation
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(title = "Pixel Classification by Individual Colony Background",
         subtitle = paste("Total colonies analyzed:", total_colonies))

  return(class_plot)
}
