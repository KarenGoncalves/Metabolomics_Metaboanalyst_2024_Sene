devtools::source_gist("https://gist.github.com/KarenGoncalves/0db105bceff4ff69547ee25460dda978")

install_from_dif_sources(
  cran_packages = c("tidyverse")
)

get_positive_ids <- function(data, group_names, IDcolName, threshold ) {
  stopifnot(exprs = is.numeric(threshold))
  
  positive_ids <- data.frame(sapply(group_names, function(group_name) {
    group_columns <- grep(group_name, colnames(data), value = TRUE)
    # Get the averages
    if (length(group_columns) < 2) {
      data[, group_columns] > threshold
    } else {
      average_value <- apply(data[, group_columns], 2, as.numeric) %>%
        as.data.frame %>% rowMeans()
      # Get the number of replicates in which they appear
      numberOfReplicatesPresent <- apply(data[, group_columns], 2, function(x) {
        x > threshold
      }) %>% rowSums
      # Return the positive IDs
      average_value > threshold & numberOfReplicatesPresent >= 2
    }
  })) %>% cbind(ID = data[, IDcolName]) %>% as.data.frame
  return(positive_ids)
}
