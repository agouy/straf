#' Generate file conversion tab UI.
#' @export
#' @keywords internal
file_conv_UI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "File conversion",
    h3("Population genetics software"),
    downloadButton(ns('dlGenepop'), 'Download file in the Genepop format'),
    tags$br(),
    tags$br(),
    downloadButton(ns('dlArlequin'), 'Download file in the Arlequin format (diploid data only)'),
    tags$br(),
    tags$br(),
    h3("Allele frequencies - reference datasets"),
    uiOutput(ns("selectPopAF")),
    tags$br(),
    downloadButton(ns('dlFamilias'), 'Download file in the Familias format'),
    tags$br(),
    tags$br(),
    downloadButton(ns('dlEuroformix'), 'Download file in the Euroformix format'),
    tags$br(),
    tags$br(),
    downloadButton(ns('dlSTRmix'), 'Download file in the STRmix format'),
    tags$br(),
    tags$br(),
    downloadButton(ns('dlLRmix'), 'Download file in the LRmix format'),
    tags$br(),
    tags$br()
    
  )
}

#' Generate file conversion tab Server
#' @export
#' @keywords internal
file_conv_Server <- function(id, fpath, ploidy, getgenind, popnames) {
  moduleServer(
    id,
    function(input, output, session) {
      
      ns <- session$ns
      
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
          paste(input$selectPopAF, '_', 'Familias.csv', sep = '') 
        },
        content = function(file) {
          straf2familias(fpath(), output_file = file, pop = input$selectPopAF)
        }
      )
      
      alleleFreqTabsAF <- reactive({
        req(getgenind())
        getFreqAllPop(getgenind())
      })
      
      output$selectPopAF <- renderUI({
        selectInput(ns("selectPopAF"), "Select a population:", popnames())
      })
      
      output$dlEuroformix <- downloadHandler(
        filename = function() { 
          paste(input$selectPopAF, '_', 'Euroformix.csv', sep = '') 
        },
        content = function(file) {
          if(input$selectPopAF == "") matr <- alleleFreqTabsAF()[[1]]
          else matr <- alleleFreqTabsAF()[[input$selectPopAF]]
          write.table(matr[-nrow(matr), ], file, sep = ",", na = "", row.names = FALSE, quote = FALSE)
        }
      )
      # 
      output$dlSTRmix <- downloadHandler(
        filename = function() { 
          paste(input$selectPopAF, '_', 'STRmix.csv', sep = '') 
        },
        content = function(file) {
          if(input$selectPopAF == "") matr <- alleleFreqTabsAF()[[1]]
          else matr <- alleleFreqTabsAF()[[input$selectPopAF]]
          write.table(matr, file, sep = ",", na = "0", row.names = FALSE, quote = FALSE)
        }
      )
      
      output$dlLRmix <- downloadHandler(
        filename = function() { 
          paste(input$selectPopAF, '_', 'LRmix.csv', sep = '') 
        },
        content = function(file) {
          if(input$selectPopAF == "") matr <- alleleFreqTabsAF()[[1]]
          else matr <- alleleFreqTabsAF()[[input$selectPopAF]]
          write.table(matr[-nrow(matr), ], file, sep = ",", na = "", row.names = FALSE, quote = FALSE)
        }
      )
      
    }
  )}

#' Convert STRAF file to the Genepop format.
#' @export
#' @keywords internal
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
#' @keywords internal
straf2familias <- function(input_file, output_file, pop = "all") {
  
  df <- readLines(input_file)
  df <- df[df != ""]
  
  spt <- do.call("rbind", strsplit(df, "\t"))
  colnames(spt) <- spt[1, ]
  df <- as.data.frame(spt[-1, ])
  
  if(pop != "all") {
    df <- df[df$pop == pop, ]
  }
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
#' @keywords internal
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



#' Convert POPTREE file to the STRAF allele frequencies format.
#' This function
#' @param fname Input file name in the STRAF format.
#' @return NULL
#' @export
#' @keywords internal
poptree2straf <- function(input_file, output_file) {
  
  fname <- input_file
  if(tools::file_ext(fname) %in% c("xls", "xlsx")) {
    df <- readxl::read_excel(fname, sheet = 1, skip = 0, col_names = FALSE)
  } else {
    no_col <- max(utils::count.fields(fname, sep = ""))
    df <- read.table(fname, 
                     comment.char = "",
                     header = T, sep = "\t", col.names = 1:no_col, fill = TRUE)
    df[df == ""] <- NA
    df <- df[, !colSums(is.na(df)) == nrow(df)]
  }
  
  popnames <- df[[1]]
  popnames <- popnames[2:(which(is.na(popnames))[1] - 1)]
  n_pops <- length(popnames)
  
  a <- grep("locus", df[[1]])
  b <- grep("#", df[[1]])
  if(length(a) != length(b)) stop("There has been an issue with file conversion.")
  str_out <- c()
  for(i in seq_along(a)) {
    
    df_tmp <- df[a[i]:(b[i]-1), ]
    pop_names <- df_tmp[1, 4:(4+n_pops-1)]
    Locus <- df_tmp[1, 2][[1]]
    freq <- df_tmp[, 4:(4+n_pops-1)]
    all_names <- c("Allele", unname(unlist(df_tmp[-1, 2])))
    df_tmp2 <- cbind(all_names, freq)
    df_tmp3 <- apply(df_tmp2, 1, paste0, collapse = ",")
    str_out <- c(str_out, Locus, df_tmp3)
  }
  str_out <- gsub(",[.]", ",0." , str_out)
  
  output <- paste(str_out, "\n", collapse = "")
  cat(output, file = output_file)
  return(NULL)
}

