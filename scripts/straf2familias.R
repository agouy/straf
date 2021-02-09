### straf to familias

straf2familias <- function(f.name) {
  
  df <- readLines(f.name)
  
  spt <- do.call("rbind", strsplit(df, "\t"))
  colnames(spt) <- spt[1, ]
  df <- as.data.frame(spt[-1, ])
  
  df_tmp <- df[, -1:-2]
  
  # add leading zeros
  # concatenate
  idx <- seq_len(length(df_tmp))
  ids <- as.logical(idx %% 2)
  
  df_out <- list()
  for(i in idx[ids]) {
    nm <- names(df_tmp[i])
    nm <- gsub(" ", "", nm)
    df_out[[nm]] <- c(df_tmp[[i]], df_tmp[[i + 1]])
  }
  
  
  tbs <- lapply(df_out, table)
  prop.tb <- lapply(tbs, prop.table)
  
  str.list <- lapply(prop.tb, function(x) {
    vec <- x
    str_loc <- paste0(names(vec), "\t", unname(vec), collapse = "\n")
    return(str_loc)
  })
  
  output <- paste(names(prop.tb), str.list, sep = "\n")
  out <- paste(output, collapse = "\n\n")
  out <- paste0(out, "\n")
  return(out)
}
str <- straf2familias("C:/Users/alexa/Desktop/Ausgangsdatei_STRAF_32markers_6pop.txt")
cat(str)


, file = "./scripts/testfamilias.txt")

    