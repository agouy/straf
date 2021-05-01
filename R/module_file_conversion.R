#' Generate file conversion tab UI.
#' @export
#' @noRd
file_conv_UI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "File conversion",
    h3("Genepop"),
    downloadButton(ns('dlGenepop'), 'Download file in the Genepop format'),
    tags$br(),
    h3("Familias"),
    downloadButton(ns('dlFamilias'), 'Download file in the Familias format'),
    tags$br(),
    h3("Arlequin (diploid data only)"),
    downloadButton(ns('dlArlequin'), 'Download file in the Arlequin format'),
    tags$br()
  )
}

#' Generate file conversion tab Server
#' @export
#' @noRd
file_conv_Server <- function(id, fpath, ploidy) {
  moduleServer(
    id,
    function(input, output, session) {
      
      #### FILE CONVERSION
      output$dlGenepop <- downloadHandler(
        filename = function() { 
          paste('straf2genepop.txt', sep='') 
        },
        content = function(file) {
          straf2genepop(input_file = fpath(), output_file = file, ploidy = ploidy())
        }
      )
      
      output$dlArlequin <- downloadHandler(
        filename = function() { 
          paste('straf2arlequin.arp', sep='') 
        },
        content = function(file) {
          straf2arlequin(fpath(), file)
        }
      )
      
      output$dlFamilias <- downloadHandler(
        filename = function() { 
          paste('straf2familias.txt', sep='')
        },
        content = function(file) {
          straf2familias(fpath(), output_file = file)
        }
      )
    }
  )}

#' Convert STRAF file to the Genepop format.
#' @export
#' @noRd
straf2genepop <- function(input_file, output_file, ploidy = 2) {
  df <- readLines(input_file)
  
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
  
  if(ploidy == 2) {
    ids <- idx %% 2
    
    df_out <- list()
    for(i in idx[as.logical(ids)]) {
      nm <- names(df_tmp2[i])
      nm <- gsub(" ", "", nm)
      df_out[[nm]] <- paste0(df_tmp2[[i]], df_tmp2[[i + 1]])
    }
  } else if (ploidy == 1) {
    df_out <- list()
    for(i in idx) {
      nm <- names(df_tmp2[i])
      nm <- gsub(" ", "", nm)
      df_out[[nm]] <- df_tmp2[[i]]
    }
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
      "Pop",
      paste(df[idx, ]$ind, str_out[idx], sep = "\t,\t")
    )
  }
  
  ## write file
  output <- c(first.line, loci, vec_out)
  output <- paste(output, "\n", collapse = "")
  cat(output, file = output_file)  
  return(NULL)
}

#' Convert STRAF file to the Familias format.
#' @export
#' @noRd
straf2familias <- function(input_file, output_file) {
  
  df <- readLines(input_file)
  
  spt <- do.call("rbind", strsplit(df, "\t"))
  colnames(spt) <- spt[1, ]
  df <- as.data.frame(spt[-1, ])
  
  df_tmp <- df[, -1:-2]
  
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
  cat(out, file = output_file)
  return(NULL)
}

#' Convert STRAF file to the Arlequin format.
#' This function
#' @param fname Input file name in the STRAF format.
#' @return NULL
#' @export
#' @noRd
straf2arlequin <- function(input_file, output_file) {
  df <- readLines(input_file)
  
  spt <- do.call("rbind", strsplit(df, "\t"))
  colnames(spt) <- spt[1, ]
  df <- as.data.frame(spt[-1, ])
  
  df_tmp2 <- list()
  
  for(i in seq_len(ncol(df[, -1:-2]))) {
    x <- df[, -1:-2][, i]
    if(i %% 2 != 0) {
      x2 <- df[, -1:-2][, i + 1]
      dot_idx <- grep("[.]", x)
      dot_idx2 <- grep("[.]", x2)
      if(length(dot_idx) > 0 | length(dot_idx2) > 0) {
        x <- gsub("[.]", "", x)
        x[-dot_idx] <- paste0(x[-dot_idx], "0")
        x2 <- gsub("[.]", "", x2)
        x2[-dot_idx2] <- paste0(x2[-dot_idx2], "0")
        if(length(dot_idx) == 0) x <- paste0(x, "0")
        if(length(dot_idx2) == 0) x2 <- paste0(x2, "0")
      }
      df_tmp2[[i]] <- x
      df_tmp2[[i+1]] <- x2
      
    } 
    
  }
  df_tmp2 <- as.data.frame(df_tmp2)
  
  # concatenate
  idx <- seq_len(length(df_tmp2))
  ids <- as.logical(idx %% 2)
  
  out_str <- c()
  for(pp in unique(df$pop)) {
    df_pop <- df[df$pop == pp, ]
    out_str <- c(out_str, paste0('SampleName="', pp, '"\nSampleSize=',nrow(df_pop),'\nSampleData={\n'))
    
    for(i in which(df$pop == pp)) {
      samp_nm <- df[i, 1]
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
  cat(output, file = output_file)
  return(NULL)
}

