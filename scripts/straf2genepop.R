### straf to genepop

straf2genepop <- function(f.name) {
  df <- readLines(f.name)
  
  spt <- do.call("rbind", strsplit(df, "\t"))
  colnames(spt) <- spt[1, ]
  df <- as.data.frame(spt[-1, ])
  
  df_tmp <- lapply(df[, -1:-2], function(x) gsub("[.]", "", x))
  
  # add leading zeros
  df_tmp2 <- lapply(df_tmp, function(x) {
    x[nchar(x) == 1] <- paste0(x[nchar(x) == 1], "00")
    x[nchar(x) == 2] <- paste0(x[nchar(x) == 2], "0")
    if(any(nchar(x) != 3)) stop("Error while converting allele labels.")
    return(x)
  })
  
  # concatenate
  idx <- seq_len(length(df_tmp2))
  ids <- idx %% 2
  
  df_out <- list()
  for(i in idx[as.logical(ids)]) {
    nm <- names(df_tmp2[i])
    nm <- gsub(" ", "", nm)
    df_out[[nm]] <- paste0(df_tmp2[[i]], df_tmp2[[i + 1]])
  }
  
  df_out <- as.data.frame(df_out)
  
  
  first.line <- "STRAF-generated GENEPOP input file."
  loci <- colnames(df_out)
  
  str_out <- apply(df_out, 1, paste0, collapse = "\t")
  
  ## get pops
  populations <- unique(df$pop)
  vec_out <- c()
  for(i in populations) {
    idx <- df$pop %in% i
    vec_out <- c(
      vec_out,
      paste0("Pop\t", i),
      paste(df[idx, ]$ind, str_out[idx], sep = "\t,\t")
    )
  }
  
  ## write file
  output <- c(first.line, loci, vec_out)
  output <- paste(output, "\n", collapse = "")
  
  return(output)
}

str <- straf2genepop("./www/exampleSTRAFdiplo.txt")
cat(str, file = "scripts/test.txt")
