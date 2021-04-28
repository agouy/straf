#' Generate the referecne population tab UI.
#' @export
#' @noRd
#' @importFrom shinyWidgets awesomeCheckbox 
ref_mds_UI <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "Reference population",
    h3("Custom allele frequency database"),
    fileInput(ns("refdata"), "Import allele frequencies (if empty, STRidER database will be used.)"),
    awesomeCheckbox(ns("add_current_ref"), "Include uploaded data to the MDS", FALSE),
    uiOutput(ns('plotMDS_ref')),
    uiOutput(ns('select_ref_pops')),
    div(
      "If no reference frequencies are uploaded, the MDS is performed on STRidER allele frequencies. Missing frequencies are imputed per allele as the mean allele frequency of other populations.", 
      tags$a(href = "https://strider.online/frequencies", "Link to source data.")
    ),
    tags$hr()
  )
}

#' Generate the referecne population tab Server.
#' @export
#' @noRd
#' @importFrom adegenet genind2genpop
ref_mds_Server <- function(id, getgenind) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- NS(id)
      
      getRefData <- reactive({
        if(is.null(input$refdata)) { 
          fr <- freq_to_mds("./www/STRidER_frequencies_2019-08-02.csv")
        } else { 
          fr <- freq_to_mds(input$refdata$datapath)
        }
        idx_in <- colSums(is.na(fr)) == 0
        fr <- fr[, idx_in]
        return(fr)
      })
      output$plotMDS_strider <- renderUI({
        plotOutput(ns('runMDS_strider'))
      })
      common_alleles <- reactive({
        X <- getRefData()[input$refpops, ]
        X <- X[, colSums(is.na(X)) == 0]
        obj <- genind2genpop(getgenind(), quiet = FALSE)
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
        X <- common_alleles()
        str_out <- paste0(
          "Alleles in common used to perform MDS: ",
          paste0(colnames(X), collapse = "; "),
          collapse = ""
        )
        return(str_out)
      })
    
      output$plotMDS_ref <- renderUI({
        plotOutput(ns('runMDS_ref'))
      })
      ref_pops <- reactive({rownames(getRefData())})
      output$select_ref_pops <- renderUI({
        all <- ref_pops()
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
        
        X <- getRefData()[input$refpops, ]
        if(is.null(X)) return(NULL)
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
        
        mds <- cmdscale(d)
        MDS <- data.frame(ax1 = mds[, 1], ax2 = mds[, 2], pop = rownames(mds))
        
        p <- ggplot(MDS, aes(x=.data$ax1, y=.data$ax2, color = pop, label = pop)) +
          geom_point() +
          geom_text_repel(max.overlaps = 50) + 
          labs( x = "MDS Axis 1", y = "MDS Axis 2", title = "MDS based on Nei's distance")  +
          theme_minimal()
        plot(p)
        
      })  
      

    }
  )
}