#' Generate input data tab UI.
#' @export
#' @noRd
data_UI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Data",
    h3("Dataset"),
    div(DT::dataTableOutput(ns('contents')), style = "font-size:70%"),
    
    tags$hr(),
    h3("Allele frequencies per locus"),
    shinyWidgets::awesomeCheckbox(
      ns('displayAlleleFreq'), 'Plot the distribution of allele frequencies',
      FALSE
    ),
    conditionalPanel(
      condition = 'input.displayAlleleFreq == true',
      ns = ns,
      uiOutput(ns('plotAF'))
    ),
    
    tags$hr(),
    shinyWidgets::awesomeCheckbox(
      ns('displayAlleleTable'),
      'Display a table of allele frequencies',
      FALSE
    ),
    conditionalPanel(
      condition = "input.displayAlleleTable == true",
      ns = ns,
      uiOutput(ns("selectPop")),
      div(DT::dataTableOutput(ns('tableFreq')), style = "font-size:70%"),
      downloadButton(ns('dlTabfreq'), 'Download as text (.tsv)'),
      downloadButton(ns('dlTabfreqXL'), 'Download as Excel (.xlsx)')
    ),
    
    tags$hr()
  )
}

#' Generate input data tab Server.
#' @export
#' @noRd
data_Server <- function(id, getgenind, getData, barplotcolor, transparency, width, height, popnames) {
  moduleServer(
    id,
    function(input, output, session) {
      
      ns <- session$ns
      
      output$contents <- DT::renderDataTable({
        if (is.null(getData())) return(NULL)
        X <- read.table(
          getData()$datapath,
          header = TRUE, sep = "\t",
          colClasses = "character"
        )
        DT::datatable(X)
      })
      
      output$alleleFreq <- renderPlot({ # barplots of allele freq distributions
        
        if (!input$displayAlleleFreq)  return(NULL)
        dat2 <- getgenind()
        freq <- apply(dat2@tab, 2, sum, na.rm = TRUE) #counts number of alleles
        nam <- strsplit(names(freq), split = "[.]") #split locus and allele name
        loc <- as.factor(unlist(
          lapply(nam, function(x) x[1])
        ))
        alle <- as.numeric(unlist(
          lapply(nam, function(x) sub("-", ".", x[2]))
        ))
        DAT <- data.frame(freq, loc, alle)
        DAT <- DAT[order(DAT$alle), ]
        
        ###depending on the number of loci, different number of columns:
        nL <- length(unique(DAT$loc))
        if(nL <= 5) n_col <- 2
        if(nL >= 6) n_col <- 3
        if(nL >= 10) n_col <- 4
        
        par(mfrow = c(ceiling(nL / n_col), n_col), mar = rep(2, 4))
        for(i in unique(DAT$loc)) {
          barplot(
            DAT$freq[DAT$loc == i],
            names.arg = DAT$alle[DAT$loc == i],
            main = i,
            col = adegenet::transp(barplotcolor(), transparency()),
            border = 0
          )
        }
      })
      
      output$plotAF <- renderUI({ #display UI only if allele freq is checked
        plotOutput(
          ns('alleleFreq'),
          width = paste(width(), "%", sep = ""),
          height = height()
        )
      })

      alleleFreqTabs <- reactive({
        if (!input$displayAlleleTable)  return(NULL)
        dat2 <- getgenind()
        matr <- getFreqAllPop(dat2)
        return(matr)
      })
      
      output$selectPop <- renderUI({
        selectInput(ns("selectPop"), "Select a population:", popnames())
      })
      
      output$tableFreq <- DT::renderDataTable({
        if (!input$displayAlleleTable | is.null(input$selectPop)) return(NULL)
        if(input$selectPop == "") matr <- alleleFreqTabs()[[1]]
        else matr <- alleleFreqTabs()[[input$selectPop]]
        DT::datatable(matr) %>% DT::formatRound(columns = colnames(matr), digits = 3)
      })
      
      output$dlTabfreq <- downloadHandler(
        filename = function() { paste('allele_frequencies.tsv', sep='') },
        content = function(file) {
          if (!input$displayAlleleTable) return(NULL)
          if(input$selectPop == "") matr <- alleleFreqTabs()[[1]]
          else matr <- alleleFreqTabs()[[input$selectPop]]
          write.table(matr, file, sep = "\t", na = "", row.names = TRUE)
        }
      )
      
      output$dlTabfreqXL <- downloadHandler(
        filename = function() { paste('allele_frequencies.xlsx', sep='') },
        content = function(file) {
          if (!input$displayAlleleTable) return(NULL)
          if(input$selectPop == "") matr <- alleleFreqTabs()[[1]]
          else matr <- alleleFreqTabs()[[input$selectPop]]
          openxlsx::write.xlsx(list(allele_frequencies = matr), file = file, rowNames = TRUE)
        }
      )
      
    }
  )
}