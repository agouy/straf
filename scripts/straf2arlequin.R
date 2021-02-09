f.name <- "./www/exampleSTRAFdiplo.txt"
straf2arlequin <- function(f.name) {
  df <- readLines(f.name)
  
  spt <- do.call("rbind", strsplit(df, "\t"))
  colnames(spt) <- spt[1, ]
  df <- as.data.frame(spt[-1, ])
  
  df_tmp2 <- lapply(df[, -1:-2], function(x) {
    dot_idx <- grep("[.]", x)
    x <- gsub("[.]", "", x)
    if(length(dot_idx) > 0) x[-dot_idx] <- paste0(x[-dot_idx], "0")
    return(x)
  })
  df_tmp2 <- as.data.frame(df_tmp2)
  
  # concatenate
  idx <- seq_len(length(df_tmp2))
  ids <- as.logical(idx %% 2)
  
  out_str <- c()
  for(pp in unique(df$pop)) {
    df_pop <- df[df$pop == pp, ]
    out_str <- c(out_str, paste0('SampleName="', pp, '"\nSampleSize=',nrow(df_pop),'\nSampleData={\n'))
    
    for(i in seq_len(nrow(df_pop))) {
      samp_nm <- df_pop[i, 1]
      l1 <- paste0(c(samp_nm, "1", unname(unlist(df_tmp2[i, c(idx[ids])]))), collapse = "\t")
      l2 <- paste0(c("", "", unname(unlist(df_tmp2[i, c(idx[!ids])]))), collapse = "\t")
      out_str <- c(out_str, l1, l2)
    }
    out_str <- c(out_str, "}\n\n")
  }
  npop <- length(unique(df$pop))
  header <- paste0('[Profile]\nTitle="STRAF-generated Arlequin file."\nNbSamples=',npop,'\nDataType=MICROSAT\n
GenotypicData=1\nGameticPhase=0\nMissingData="?"\nLocusSeparator=WHITESPACE\n\n[Data]\n\n[[Samples]]\n\n')
  
  output <- c(header, out_str)
  ## write file
  output <- paste(output, "\n", collapse = "")
  
  return(output)
}

str <- straf2arlequin("./www/exampleSTRAFdiplo.txt")
cat(str, file = "scripts/testoutput.arp")



