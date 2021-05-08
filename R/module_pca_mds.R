#' Generate the PCA UI.
#' @export
#' @noRd
pca_mds_UI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "PCA - MDS",
    h3("Principal Component Analysis (PCA)"),
    awesomeCheckbox(
      ns('displayPCA'),
      'Run and plot a PCA (Principal Component Analysis)',
      FALSE
    ),
    conditionalPanel(
      condition = "input.displayPCA == true", ns = ns,
      checkboxGroupInput(ns('PCAaxis'), 'PCA axis', c(1,2,3), c(1,2), inline = TRUE)
    ),
    conditionalPanel(
      condition = "input.displayPCA == true", ns = ns,
      uiOutput(ns('plotPCA')),
      verbatimTextOutput(ns('info')),
      downloadButton(ns('dlPCAeigen'), 'Download PCA eigenvectors'),
      downloadButton(ns('dlPCAcoord'), 'Download PCA coordinates'),
      awesomeCheckbox(ns('displayloadings'), 'Plot loadings (alleles contributions)', FALSE)
    ),
    conditionalPanel(
      condition = "input.displayPCA == true & input.displayloadings == true",
      ns = ns,
      uiOutput(ns('plotLoadings'))
    ),
    
    tags$hr(),
    h3("Multidimensional Scaling (MDS) based on Nei's distance"),
    awesomeCheckbox(
      ns('displayMDS'),
      "Compute Nei's genetic distance between populations and run MDS",
      FALSE
    ),
    conditionalPanel(
      condition = "input.displayMDS == true",
      ns = ns,
      uiOutput(ns('plotMDS'))
    ),
    
    tags$hr()
  )
}

#' Generate the PCA server.
#' @export
#' @noRd
pca_mds_Server <- function(id, getgenind) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      output$runPCA <- renderPlot({
        if (!input$displayPCA)  return(NULL)
        dat2 <- getgenind()
        pca.obj <- do.pca()
        
        if(length(input$PCAaxis)==2){
          par(mfrow=c(1,1))
          coul <- transp(funky(length(unique(pop(dat2)))),.6)
          plotPCA(
            pca = pca.obj,
            popus = pop(dat2),
            coul = coul,
            axis = c(as.numeric(input$PCAaxis[1]), as.numeric(input$PCAaxis[2]))
          )
        }
        if(length(input$PCAaxis)==3){
          
          par(mfrow = c(1,2))
          coul <- transp(funky(length(unique(pop(dat2)))),.6)
          plotPCA(
            pca = pca.obj,
            popus = pop(dat2),
            coul = coul,
            axis = c(as.numeric(input$PCAaxis[1]), as.numeric(input$PCAaxis[2]))
          )
          plotPCA(
            pca = pca.obj,
            popus = pop(dat2),
            coul = coul,
            axis = c(as.numeric(input$PCAaxis[1]), as.numeric(input$PCAaxis[3]))
          )
          
        }
      })
      
      output$plotPCA <- renderUI({
        plotOutput(ns('runPCA'), click = ns("plot_click"))
      })
      output$plotMDS <- renderUI({
        plotOutput(ns('runMDS'))
      })
      output$runMDS <- renderPlot({
        if (!input$displayMDS)  return(NULL)
        dat2 <- getgenind()
        if(length(levels(pop(dat2))) < 2) return(NULL)
        obj <- genind2genpop(dat2, quiet = TRUE)
        dst <- dist.genpop(obj, method = 1)
        MDS <- cmdscale(dst)
        MDS <- data.frame(ax1 = MDS[, 1], ax2 = MDS[, 2], pop = rownames(MDS))
        .data <- NA
        p <- ggplot(MDS, aes(x = .data$ax1, y = .data$ax2, color = pop, label = pop)) +
          geom_point() +
          geom_text_repel() + 
          labs( x = "MDS Axis 1", y = "MDS Axis 2", title = "MDS based on Nei's distance")  +
          theme_minimal()
        plot(p)
      })
      
      do.pca <- reactive({
        freq.tab <- makefreq(getgenind(), missing = "mean", quiet = TRUE)
        pca.obj <- dudi.pca(freq.tab, nf = 3, scannf = FALSE)
        
        return(pca.obj)
      })
      # DL principal components
      output$dlPCAeigen <- downloadHandler(
        filename = function() {
          paste('PCA_eigenvectors.tsv', sep = '')
        },
        content = function(file) {
          if (!input$displayPCA)  return(NULL)
          pca.obj <- do.pca()
          write.table(pca.obj$c1, file, sep = "\t", na = "",row.names = TRUE)
        }
      )
      output$dlPCAcoord <- downloadHandler(
        filename = function() {
          paste('PCA_coordinates.tsv', sep = '')
        },
        content = function(file) {
          if (!input$displayPCA)  return(NULL)
          pca.obj <- do.pca()
          write.table(pca.obj$li, file, sep = "\t", na = "",row.names = TRUE)
        }
      )
      output$info <- renderPrint({
        if (!input$displayPCA)  return(NULL)
        pca.obj <- do.pca()
        if(length(input$PCAaxis) == 2){
          cat("Click on a point to get its ID and coordinates\n\n")
          ta <- c("Axis1","Axis2","Axis3")
          if(!is.null(input$plot_click)){
            nearPoints(pca.obj$li, input$plot_click, 
                       xvar = ta[as.numeric(input$PCAaxis[1])], 
                       yvar = ta[as.numeric(input$PCAaxis[2])])
          }
        }
      })
      output$loadings <- renderPlot({
        if (!input$displayPCA) return(NULL)
        pca.obj <- do.pca()
        loadingplot(pca.obj$c1 ^ 2)
      })
      output$plotLoadings <- renderUI({
        plotOutput(ns('loadings'))
      })
      
    }
  )
}