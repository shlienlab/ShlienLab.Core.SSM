create_final_filtered_dataframe <- function(data = NULL, pon_threshold = 1, flagged_threshold = 1, soft_threshold = 0, hard_threshold = 0) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  data <- data %>%
    filter(pon_count <= pon_threshold) %>%
    filter(flagged_count <= flagged_threshold) %>%
    filter(in_hard_clips <= hard_threshold) %>%
    filter(in_soft_clips <= soft_threshold)

  return(data)
  }
