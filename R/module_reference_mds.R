#' Generate the reference population tab UI.
#' @export
#' @keywords internal
ref_mds_UI <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "Reference population",
    h3("Custom allele frequency database"),
    fileInput(
      ns("refdata"),
      "Import allele frequencies (if empty, the STRidER database will be used.)",
      width = "300px"
    ),
    awesomeCheckbox(ns("add_current_ref"), "Include uploaded data to the MDS", FALSE),
    uiOutput(ns('plotMDS_ref')),
    uiOutput(ns('plotMDStree_ref')),
    uiOutput(ns('select_ref_pops')),
    uiOutput(ns('common_all')),
    tags$hr()
  )
}

#' Generate the reference population tab Server.
#' @export
#' @keywords internal
ref_mds_Server <- function(id, getgenind) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- NS(id)
      
      getRefData <- reactive({
        if(is.null(input$refdata)) { 
          fr <- suppressWarnings({freq_to_mds("./www/STRidER_frequencies_2019-08-02.csv")})
        } else { 
          fr <- freq_to_mds(input$refdata$datapath)
        }
        idx_in <- colSums(is.na(fr)) == 0
        fr <- fr[, idx_in]
        return(fr)
      })

      output$common_all <- renderUI({
        textOutput(ns("common_all"))
      })
      
      output$plotMDS_strider <- renderUI({
        plotOutput(ns('runMDS_strider'))
      })
      common_alleles <- reactive({
        req(input$refpops, getRefData(), getgenind())
        X <- getRefData()[input$refpops, ]
        X <- X[, colSums(is.na(X)) == 0]
        obj <- genind2genpop(getgenind(), quiet = TRUE)
        cn <- colnames(obj@tab)
        cn <- gsub("[.]", "_", cn)
        cn <- gsub("[-]", ".", cn)
        colnames(obj@tab) <- cn
        common_cols <- intersect(colnames(obj@tab), colnames(X))
        if(length(common_cols) == 0) stop("No alleles in common.") 
        X <- rbind(
          subset(obj@tab, select = common_cols), 
          subset(X, select = common_cols)
        )
        return(X)
      })
      output$common_all <- renderText({
        if(!input$add_current_ref) {          
          return(NULL)
        }
        req(getRefData())
        X <- common_alleles()
        loci <- unique(unlist(lapply(strsplit(colnames(X), "_"), "[", 1)))
        str_out <- paste0(
          "Loci in common used to perform MDS: ",
          paste0(loci, collapse = "; "),
          collapse = ""
        )
        return(str_out)
      })
    
      output$plotMDS_ref <- renderUI({
        plotOutput(ns('runMDS_ref'))
      })
      ref_pops <- reactive({
        refdat <- getRefData()
        req(refdat)
        rownames(refdat)
      })
      output$select_ref_pops <- renderUI({
        req(ref_pops())
        all <- ref_pops()
        req(all)
        suppressMessages(pickerInput(
          ns('refpops'), 'Select populations',
          choices = all,
          selected = all[!all %in% "THAILAND"],
          multiple = TRUE,
          width = "100%",
          options = list(
            `actions-box` = TRUE
          )
        ))
      })
      output$runMDS_ref <- renderPlot({
        req(do.dist_ref())
        d <- do.dist_ref()
        if(is.null(d)) return(NULL)
        mds <- cmdscale(d)
        MDS <- data.frame(ax1 = mds[, 1], ax2 = mds[, 2], pop = rownames(mds))
        .data <- NA
        p <- ggplot(MDS, aes(x=.data$ax1, y=.data$ax2, color = pop, label = pop)) +
          geom_point() +
          geom_text_repel(max.overlaps = 50) + 
          labs( x = "MDS Axis 1", y = "MDS Axis 2", title = "MDS based on Nei's distance")  +
          theme_minimal()
        plot(p)
        
      })
      
      output$plotMDStree_ref <- renderUI({
        plotOutput(ns('runMDStree_ref'))
      })
      output$runMDStree_ref <- renderPlot({
        req(do.dist_ref())
        dst <- do.dist_ref()
        if(is.null(dst)) return(NULL)
        hc <- hclust(dst)
        plot(ape::as.phylo(hc), cex = 0.9)        
      })
      do.dist_ref <- reactive({
        req(getRefData())
        X <- getRefData()
        if(is.null(X)) return(NULL)
        if(!any(rownames(X) %in% input$refpops)) return(NULL)
        X <- X[input$refpops, ]
        X <- X[, colSums(is.na(X)) == 0]
        if(input$add_current_ref) {          
          X <- common_alleles()
        }
        d <- X %*% t(X)
        vec <- sqrt(diag(d))
        d <- d / vec[col(d)]
        d <- d / vec[row(d)]
        d <- -log(d)
        d <- as.dist(d)

        return(d)
      })
      

    }
  )
}


#' Convert allele frequencies to a format suitable for MDS analysis.
#' @export
#' @keywords internal
freq_to_mds <- function(fname) {
  ln <- readLines(fname)
  ln2 <- lapply(ln, function(x) strsplit(x, ",")[[1]])
  ln3 <- lapply(ln2, function(x) {
    if(sum(nchar(x[-1]) == 0) == length(x[-1])) return(x[1])
    else return(x)
  })
  hd <- lengths(ln3)
  names_idx <- which(hd == 1)
  st_idx <- names_idx + 1
  en_idx <- names_idx - 1
  en_idx <- c(en_idx[-1], length(ln2))
  df <- lapply(seq_along(names_idx), function(i) {
    loc_id <- names_idx[i]
    loc_name <- ln3[[loc_id]]
    if(en_idx[i] - st_idx[i] > 1) {
      mat <- do.call(rbind, ln2[st_idx[i]:en_idx[i]])
      colnames(mat) <- mat[1, ]
      mat[mat == ""] <- "0"
      df <- as.data.frame(mat[-1:-2, ])
      colnames(df) <- gsub(pattern = " ", replacement = "_", colnames(df))
      colnames(df) <- gsub(pattern = "\"", replacement = "", colnames(df))
      Allele <- NA
      location <- NA
      df_long <- tidyr::gather(df, location, frequency, -Allele, factor_key=TRUE)
      df_long$locus <- loc_name
      return(df_long)
      
    } else {
      return(NULL)
    }
  })
  df_l <- do.call(rbind, df)
  df_l$frequency <- as.numeric(df_l$frequency)
  df_l$location  <- as.character(df_l$location)
  tt <- reshape2::acast(df_l, location ~ locus + Allele, value.var = 'frequency', fun.aggregate = mean, fill = -1)
  ct <- rownames(tt)
  tt <- tt %>% dplyr::as_tibble()
  df_f <- tt %>%
    dplyr::as_tibble() %>%
    dplyr::mutate_all(~ifelse(.x == -1, NA, .x)) #mean(.x[.x != -1], na.rm = TRUE)
  matt <- (as.matrix(df_f))
  rownames(matt) <- ct
  return(matt)
}
