#' Plot enzymatic activity
#' 
#' Generates a plot for enzymatic activity (Raw, normalized activity; Total 
#' Activity; Net Activity; Net/Total Activity) based on the user input.
#' 
#' The data needs to be aggregated first before graphing.
#' 
#' @param graph.type, restricted to four different strings ('Raw' for 
#' raw, normalized colony.value data, 'Total' for total activity, 'Net' for net 
#' activity, 'Net/Total' for net/total enzymatic activity)
#' @param input.plateID, a list of plateIDs. If multiply plate IDs are given,
#' the graphs will show a trendline corresponding to the mean of the
#' replicates across plates while showing individual replicates as data points.
#' @param input.media, a string which corresponds to a media condition within 
#' the data
#' @param input.gene, a list of strings, can also be empty. Regardless, even if 
#' empty, the script will automatically include "EV_ctrl", "IsPETase_ctrl", and 
#' "FAST-PETase_ctrl".
#' @param input.plasmid a string which corresponds to a specific plasmids within
#' the data

#' @return a 'ggplot2' object.
#' @export

plotActivity <- function(plate.tbl,
                         graph.type = "Raw",
                         input.plateID = c(),
                         input.media = "BHET12.5",
                         input.gene = c(),
                         input.plasmid = "152") {
  
  ##---- sanity check --------------------------------------------------------##
  stopifnot(is.data.frame(plate.tbl))

  if (!graph.type %in% c("Raw", "Total", "Net", "Net/Total")) {
    usethis::ui_stop("Graph type must be one of 'Raw', 'Total', 'Net', or 'Net/Total'")
  }
  
  if(length(input.plateID) == 0) {
    usethis::ui_stop("input.plateID must contain at least 1 plateID")
  }
  
  if(!all(input.plateID %in% plate.tbl$plateID)) {
    missing <- input.plateID[!input.plateID %in% plate.tbl$plateID]
    usethis::ui_stop("The following plateID(s) are not present in the table: {missing}")
  }
  
  if(!input.media %in% plate.tbl$media) {
    usethis::ui_stop("The media provided is not present in the table")
  }
  
  if(!all(input.gene %in% plate.tbl$gene)) {
    missing <- input.gene[!input.gene %in% plate.tbl$gene]
    usethis::ui_stop("The following gene(s) are not present in the table: {missing}")
  }
  
  if(!input.plasmid %in% plate.tbl$plasmid) {
    usethis::ui_stop("The plasmid provided is not present in the table")
  }
  
  ## ---- prepare data ------------------------------------------------------ ##
  df <- plate.tbl %>%
        group_by(plateID, media, gene, plasmid) %>%
        arrange(timepoint, .by_group = TRUE)
  
  first_timepoint = df_changes$timepoint[1]
  
  # Calculating total activity (Summation of absolute value of the difference 
  # between timepoints)
  df <- df %>%
    mutate(change = ifelse(timepoint == first_timepoint, 0, colony.value_mean-lag(colony.value_mean))) %>%
    mutate(abs_change = abs(change)) %>%
    mutate(total_activity = {
      tmp <- numeric(n())
      
      for(i in seq_along(tmp)) {
        if(timepoint[i] == first_timepoint) {
          tmp[i] <- abs(colony.value_mean[i])     # reset at timepoint 6
        } else if(i == 1) {
          tmp[i] <- abs(colony.value_mean[i])     # first row uses colony.value_mean
        } else {
          tmp[i] <- tmp[i-1] + abs_change[i]      # cumulative addition
        }
      }
      
      tmp 
    })
  
  # Calculating net activity (Summation of the difference between timepoints)
  df = df %>%              
    mutate(net_activity = {
      tmp <- numeric(n())
      
      for(i in seq_along(tmp)) {
        if(timepoint[i] == first_timepoint) {
          tmp[i] <- colony.value_mean[i]       # reset at timepoint 6
        } else if(i == 1) {
          tmp[i] <- colony.value_mean[i]       # first row uses colony.value_mean
        } else {
          tmp[i] <- tmp[i-1] + change[i]       # cumulative addition
        }
      }
      
      tmp
    })
  
  # Calculating the activity ratio (Net activity / Total Activity)
  df = df %>%
    mutate(NetTotal_Ratio = net_activity/total_activity)
  
  input.gene <- c(input.gene, list("EV_ctrl", 
                                   "IsPETase_ctrl", 
                                   "FAST-PETase_ctrl"))
  small_subset <- filter(df,
                         plateID %in% input.plateID,
                         media   == input.media,
                         gene %in% input.gene,
                         plasmid == input.plasmid)
  
  ## ---- graphing ---------------------------------------------------------- ##
  
  small_subset_summary <- small_subset %>%
    group_by(gene, timepoint)
  
  
  if (graph.type == "Raw") {
    small_subset_summary <- small_subset_summary %>%
                            summarise (mean_raw_activity = mean(colony.value_mean, 
                                                                na.rm = TRUE), 
                                                           .groups = 'drop')
    
    ggplot(small_subset, aes(x = timepoint, 
                             y = colony.value_mean, 
                             color = gene, group = gene)) +
    geom_line(data = small_subset_summary, aes(x = timepoint, 
                                               y = mean_raw_activity, 
                                               color = gene),
                                               size = 0.6) +
    geom_point(position = position_jitter(width = 0, height = 0)) +
    labs(title = "Raw, Normalized Activity",
         x = "Time (hr)",
         y = "Background corrected pixel intensities (a.u.)") +
    theme_classic()
    
  } else if (graph.type == "Total") {
    small_subset_summary <- small_subset_summary %>%
                            summarise (mean_total_activity = mean(total_activity, 
                                                                  na.rm = TRUE), 
                                                             .groups = 'drop')
  
    ggplot(small_subset, aes(x = timepoint, 
                             y = total_activity, 
                             color = gene, group = gene)) +
    geom_line(data = small_subset_summary, aes(x = timepoint, 
                                               y = mean_total_activity, 
                                               color = gene),
                                               size = 0.6) +
    geom_point(position = position_jitter(width = 0, height = 0)) +
    labs(title = "Total Enzymatic Activity",
         x = "Time (hr)",
         y = "Σ |Δ background corrected pixel intensities| (a.u.)") +
    theme_classic()
  
  } else if (graph.type == "Net") {
      small_subset_summary <- small_subset_summary %>%
                              summarise (mean_net_activity = mean(net_activity, 
                                                                  na.rm = TRUE), 
                                                             .groups = 'drop')
    
      ggplot(small_subset, aes(x = timepoint, 
                               y = net_activity, 
                               color = gene, group = gene)) +
      geom_line(data = small_subset_summary, aes(x = timepoint, 
                                                 y = mean_net_activity, 
                                                 color = gene),
                                                 size = 0.6) +
      geom_point(position = position_jitter(width = 0, height = 0)) +
      labs(title = "Net Enzymatic Activity",
           x = "Time (hr)",
           y = "Σ Δ background corrected pixel intensities (a.u.)") +
      theme_classic()
    
  } else if (graph.type == "Net/Total") {
      small_subset_summary <- small_subset_summary %>%
                              summarise (mean_NetTotal_Ratio = mean(NetTotal_Ratio, 
                                                                    na.rm = TRUE), 
                                                               .groups = 'drop')
    
      ggplot(small_subset, aes(x = timepoint, 
                               y = NetTotal_Ratio, 
                               color = gene, group = gene)) +
      geom_line(data = small_subset_summary, aes(x = timepoint, 
                                                 y = mean_NetTotal_Ratio, 
                                                 color = gene),
                                                 size = 0.6) +
      geom_point(position = position_jitter(width = 0, height = 0)) +
      labs(title = "Net/Total Enzymatic Activity",
           x = "Time (hr)",
           y = "Σ |Δ| / Σ Δ (a.u.)") +
      theme_classic()
  }
  
}