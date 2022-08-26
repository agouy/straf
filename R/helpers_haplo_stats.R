#' Get haplotype specific statistics for all populations
#' @inheritParams getHaploStatsFromGenind
#' @export
#' @keywords internal
getHaploStatsAllPop <- function(data) {
  ind <- list()
  ind$all <- getHaploStatsFromGenind(data)
  for(popu in unique(data$pop)) {
    ind <- c(ind, x = NA)
    mat <- getHaploStatsFromGenind(data[data@pop == popu, ])
    ind$x <- mat
    names(ind)[length(ind)] <- popu
  }
  return(ind)
}

#' Get haplotype specific statistics
#' @param data A genind object.
#' @export
#' @keywords internal
getHaploStatsFromGenind <- function(data) {
  # check ploidy
  df <- adegenet::genind2df(data)
  
  hap_df <- tidyr::unite_(df, "haplotype", colnames(df)[-1])
  
  n_loc <- ncol(df) - 1
  
  h_count <- table(hap_df$haplotype)
  h_freq <- prop.table(h_count)
  h_n <- sum(h_count)
  h_div <- (h_n / (h_n - 1)) * sum(h_freq ^ 2)
  
  ## haplotype level table
  hap_sq <- seq_along(h_count)
  hap_labels <- paste0(
    "H",
    formatC(
      hap_sq,
      width = max(nchar(hap_sq)),
      flag = "0"
    )
  )
  hap_data <- data.frame(
    haplotype_id = hap_labels,
    h_count = h_count,
    h_freq = h_freq,
    h_value = gsub("_", "|", names(h_count))
  )
  
  ## population level table
  hap_pop_data <- data.frame(
    h_div = h_div
  )
  
  ## number of pairwise differences
  compdist <- function(a, b) {
    x <- strsplit(a, "|")[[1]]
    y <- strsplit(b, "|")[[1]]
    return(sum(x != y))
  }
  d_mat <- outer(X = hap_data$h_value, Y = hap_data$h_value, Vectorize(compdist))
  colnames(d_mat) <- hap_data$haplotype_id
  rownames(d_mat) <- hap_data$haplotype_id
  
  ## proportion of pairwise differences
  d_mat_prop <- d_mat / n_loc
  
  list_out <- list(
    hap_data = hap_data,
    hap_pop_data = hap_pop_data,
    d_mat = d_mat,
    d_mat_prop = d_mat_prop
  )
  return(list_out)
}
