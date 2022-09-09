#' Generate the PCA UI.
#' @export
#' @keywords internal
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
      selectInput(ns('PCAaxis1'), 'PC on x axis', 1:5, 1),
      selectInput(ns('PCAaxis2'), 'PC on y axis', 1:5, 2)
    ),
    conditionalPanel(
      condition = "input.displayPCA == true", ns = ns,
      uiOutput(ns('plotPCA')),
      downloadButton(ns('dlPCAeigen'), 'Download PCA eigenvectors'),
      downloadButton(ns('dlPCAcoord'), 'Download PCA coordinates')
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
      uiOutput(ns('plotMDS')),
      uiOutput(ns('plotMDStree'))
    ),
    
    tags$hr()
  )
}

#' Generate the PCA server.
#' @export
#' @keywords internal
pca_mds_Server <- function(id, getgenind) {
  
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      output$runPCA <- renderPlotly({
        if (!input$displayPCA)  return(NULL)
        dat2 <- getgenind()
        pc <- do.pca()
        
        df_tmp <- pc$x
        
        pc_var <- summary(pc)$importance["Proportion of Variance", ]
        
        ax_1 <- as.numeric(input$PCAaxis1)
        ax_2 <- as.numeric(input$PCAaxis2)
        
        axes_labels <- paste0(colnames(df_tmp), " (", round(pc_var * 100, 2), " %)")
        df <- as.data.frame(df_tmp) %>%
          dplyr::mutate(Population = pop(dat2), ID = rownames(dat2$tab))
        
        plt <- suppressWarnings(
          ggplot2::ggplot(
            df,
            ggplot2::aes_string(
              x = colnames(df)[ax_1],
              y = colnames(df)[ax_2],
              color = "Population"
            )
          ) +
            ggplot2::geom_point(ggplot2::aes(text = ID)) +
            ggplot2::stat_ellipse() + 
            ggplot2::theme_minimal() +
            ggplot2::labs(x = axes_labels[ax_1], y = axes_labels[ax_2]) 
        )
        
        plotly::ggplotly(plt)
      })
      
      output$plotPCA <- renderUI({
        plotlyOutput(ns('runPCA'))
      })
      output$plotMDS <- renderUI({
        plotOutput(ns('runMDS'))
      })
      output$runMDS <- renderPlot({
        if (!input$displayMDS)  return(NULL)
        req(do.dist())
        dst <- do.dist()
        MDS <- cmdscale(dst)
        MDS <- data.frame(ax1 = MDS[, 1], ax2 = MDS[, 2], pop = rownames(MDS))
        .data <- NA
        p <- ggplot(MDS, aes(x = .data$ax1, y = .data$ax2, color = pop, label = pop)) +
          geom_point() +
          geom_text_repel() + 
          labs(x = "MDS Axis 1", y = "MDS Axis 2", title = "MDS based on Nei's distance")  +
          theme_minimal()
        plot(p)
      })
      
      output$plotMDStree <- renderUI({
        plotOutput(ns('runMDStree'))
      })
      
      output$runMDStree <- renderPlot({
        req(do.dist())
        if (!input$displayMDS)  return(NULL)
        dst <- do.dist()
        hc <- stats::hclust(dst)
        plot(ape::as.phylo(hc), cex = 0.9)        
      })
      
      do.dist <- reactive({
        if (!input$displayMDS)  return(NULL)
        dat2 <- getgenind()
        if(length(levels(pop(dat2))) < 2) stop("Multiple populations are required for the MDS.")
        obj <- genind2genpop(dat2, quiet = TRUE)
        dst <- dist.genpop(obj, method = 1)

        return(dst)
      })
      
      
      do.pca <- reactive({
        freq.tab <- makefreq(getgenind(), missing = "mean", quiet = TRUE)
        pca.obj <- prcomp(freq.tab, scale = TRUE, center = TRUE)
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
          write.table(pca.obj$rotation, file, sep = "\t", na = "",row.names = TRUE)
        }
      )
      output$dlPCAcoord <- downloadHandler(
        filename = function() {
          paste('PCA_coordinates.tsv', sep = '')
        },
        content = function(file) {
          if (!input$displayPCA)  return(NULL)
          pca.obj <- do.pca()
          write.table(pca.obj$x, file, sep = "\t", na = "",row.names = TRUE)
        }
      )
      
    }
  )
}
