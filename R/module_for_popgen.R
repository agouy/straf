### Module
for_UI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Forensic parameters",
    h3("Forensic parameters"),
    awesomeCheckbox(
      ns('displayForensics'),
      'Compute forensics statistics (H, GD, PIC, PD, PE & TPI)',
      FALSE
    ),
    
    conditionalPanel(
      condition = "input.displayForensics == true",
      ns = ns,
      uiOutput(ns("selectPop2")),
      div(tableOutput(ns('forensics')), style = "font-size:75%"),
      downloadButton(ns('dlForensics'), 'Download as text (.tsv)'),
      downloadButton(ns('dlForensicsXL'), 'Download as Excel (.xlsx)'),
      tags$hr(),
      uiOutput(ns("uiFOR")),
      uiOutput(ns("plotFOR"))
    ),
    tags$hr()
  )
}
popgen_UI <- function(id) {
  ns <- NS(id)    
  tabPanel(
    "Population genetics",
    h3("Summary statistics"),
    awesomeCheckbox(
      ns('displayDiv'),
      'Compute heterozygosities and F-statistics',
      FALSE
    ),
    
    conditionalPanel(
      condition = "input.displayDiv == true",
      ns = ns,
      uiOutput(ns("selectPop3"))
    ),
    
    conditionalPanel(
      condition = "input.ploidy == 'Diploid'",
      ns = ns,
      awesomeCheckbox(
        ns('computeHW'), 'Test for Hardy-Weinberg equilibrium',
        FALSE
      ),
      numericInput(
        ns('hw_nperm'), 'Number of permutations for HW test',
        1000, min = 100, max = 10000, step = 100
      )
    ),
    conditionalPanel(
      condition = "input.displayDiv == true",
      ns = ns,
      div(tableOutput(ns('diversity')), style = "font-size:75%"),
      downloadButton(ns('dlPopgen'), 'Download as text (.txt)'),
      downloadButton(ns('dlPopgenCL'), 'Download as Excel (.xlsx)'),
      tags$hr(),
      uiOutput(ns("uiPG")),
      uiOutput(ns("plotPG"))
    ),
    
    tags$hr(),
    h3("Linkage disequilibrium"),
    awesomeCheckbox(
      ns('displayLDtable'), 'Display pairwise LD p-values matrix',
      FALSE
    ),
    conditionalPanel(
      condition = "input.displayLDtable == true",
      ns = ns,
      div(tableOutput(ns('LDtable')), style = "font-size:75%"),
      downloadButton(ns('dlLDtable'), 'Download as text'),
      downloadButton(ns('dlLDtableXL'), 'Download as Excel')
    ),
    conditionalPanel(
      condition = "output.LD30",
      ns = ns,
      awesomeCheckbox(
        ns('displayLDplot'),
        'Plot pairwise LD p-values matrix',
        FALSE),
      conditionalPanel(
        condition = "input.displayLDplot == true",
        ns = ns,
        uiOutput(ns("plotLD2"))
      ),
      conditionalPanel(
        condition = "input.displayLDplot == true | input.displayLDtable == true",
        ns = ns,
        awesomeCheckbox(
          ns('displayLDpvalplot'),
          'Plot LD p-values distribution',
          FALSE
        ),
        conditionalPanel(
          condition = "input.displayLDpvalplot == true",
          ns = ns,
          uiOutput(ns("plotLDpval2"))
        )
      )
    ),
    
    tags$hr(),
    h3("Pairwise Fst"),
    awesomeCheckbox(
      ns('displayFstMat'),
      'Compute pairwise Fst matrix',
      FALSE
    ),
    conditionalPanel(
      condition = "input.displayFstMat == true",
      ns = ns,
      div(tableOutput(ns('FstMat')), style = "font-size:75%"),
      downloadButton(ns('dlFstMat'), 'Download as text (.txt)'),
      downloadButton(ns('dlFstMatXL'), 'Download as Excel (.xlsx)')
    ),
    tags$hr()
  )
}
for_popgen_Server <- function(id, getgenind, popnames, ploidy) {
  moduleServer(
    id,
    function(input, output, session) { 
      ns <- session$ns
      output$selectPop2 <- renderUI({
        selectInput(ns("selectPop2"), "Select a population:", popnames())
      })
      output$selectPop3 <- renderUI({
        selectInput(ns("selectPop3"), "Select a population:", popnames())
      })
      reacIndices <- reactive({
        if(is.null(getgenind())) return(NULL)
        
        DF <- getIndicesAllPop(
          getgenind(),
          hw = input$computeHW,
          hwperm = input$hw_nperm,
          ploidy = ploidy()
        )
        return(DF)
      })
      
      ### forensics table display
      output$forensics <- renderTable({
        
        if (is.null(getgenind())) return(NULL)
        if(is.null(input$selectPop2)) taB <- reacIndices()[[1]]
        else taB <- reacIndices()[[input$selectPop2]]
        
        if(!is.null(input$selectPop2)) {
          dat2 <- getgenind()
          if(length(unique(dat2@pop)) > 1 & input$selectPop2 == "all") {
            taB[, ! colnames(taB) %in% c("Ht", "Fis", "Fst")]
          } else {
            taB
          }
        } else {
          taB
        }
        
      }, digits = 4)
      
      ### forensics table DL
      output$dlForensics <- downloadHandler(
        filename = function() { paste('forensics_parameters.tsv', sep = '') },
        content = function(file) {
          if(is.null(input$selectPop2)) taB <- reacIndices()[[1]]
          else taB <- reacIndices()[[input$selectPop2]]
          write.table(
            taB[, ! colnames(taB) %in% c("Ht", "Fis", "Fst")],
            file, sep = "\t", row.names = FALSE
          )
        }
      )
      output$dlForensicsXL <- downloadHandler(
        filename = function() { paste('forensics_parameters.xlsx', sep = '') },
        content = function(file) {
          if(is.null(input$selectPop2)) taB <- reacIndices()[[1]]
          else taB <- reacIndices()[[input$selectPop2]]
          out <- taB[, ! colnames(taB) %in% c("Ht", "Fis", "Fst")]
          openxlsx::write.xlsx(list(forensics_parameters = out), file, row.names = FALSE)
        }
      )
      
      ### popgen table display + download
      output$diversity <- renderTable({
        if(is.null(input$selectPop3)) taB <- reacIndices()[[1]]
        else taB <- reacIndices()[[input$selectPop3]]
        taB2 <- taB[, !colnames(taB) %in% c("PIC", "PD", "PE", "TPI", "PM")]
        taB2
      }, digits = 4)
      
      output$dlPopgen <- downloadHandler(
        filename = function() {
          paste('popgen_statistics.tsv', sep='') 
        },
        content = function(file) {
          if(is.null(input$selectPop3)) taB <- reacIndices()[[1]]
          else taB <- reacIndices()[[input$selectPop3]]
          taB2 <- taB[, !colnames(taB) %in% c("PIC","PD","PE","TPI", "PM")]
          write.table(taB2, file, sep = "\t", row.names = FALSE)
        }
      )
      output$dlPopgenXL <- downloadHandler(
        filename = function() {
          paste('popgen_statistics.xlsx', sep='') 
        },
        content = function(file) {
          if(is.null(input$selectPop3)) taB <- reacIndices()[[1]]
          else taB <- reacIndices()[[input$selectPop3]]
          taB2 <- taB[, !colnames(taB) %in% c("PIC","PD","PE","TPI", "PM")]
          openxlsx::write.xlsx(
            list(popgen_parameters = taB2),
            file,
            row.names = FALSE
          )
        }
      )  
      ### plot forensics parameters + popgen indices
      
      # reactive UI: displayed if indices are computed
      output$uiFOR <- renderUI({
        if(is.null(input$selectPop2)) return(NULL)
        else taB <- reacIndices()[[input$selectPop2]]
        selectInput(ns("plotIndicesFOR"), "Select the statistic to plot:",
                    choices = colnames(taB)[-1:-2],
                    selected = "Nall"
        )
      })
      
      # reactive UI: displayed if indices are computed
      output$uiPG <- renderUI({
        if(is.null(input$selectPop3)) return(NULL)
        else taB <- reacIndices()[[input$selectPop3]]
        selectInput(ns("plotIndicesPG"), "Select the statistic to plot:",
                    choices = colnames(taB)[-1:-2],
                    selected = "Nall"
        )
      })
      
      # plot in the reactive UI
      output$plotIndices <-  renderPlotly({
        
        if(is.null(input$selectPop2)) return(NULL)
        else taB <- reacIndices()[[input$selectPop2]]
        
        if(is.null(taB) || dim(taB)[1] == 0) return(NULL)
        dat <- taB
        
        if(is.null(input$plotIndicesFOR)) return(NULL)
        
        fig <- plot_ly(
          x = dat[, input$plotIndicesFOR],
          y = dat[, "locus"],
          name = "bp_for",
          type = "bar",
          text = paste0(
            "Locus: ", dat[, "locus"],
            "\nValue: ", round(dat[, input$plotIndicesFOR], 4)
          ),
          hoverinfo = 'text',
          marker = list(color = transp(input$barplotcolor, input$transparency))
        ) %>% layout(
          xaxis = list(title = input$plotIndicesFOR, zeroline = FALSE),
          yaxis = list(title = "Locus")
        )
        
        (fig)
      })
      
      output$plotFOR <- renderUI({
        plotlyOutput(ns('plotIndices'))
      })
      
      output$plotIndicesPopgen <-  renderPlot({
        
        if(is.null(input$selectPop3)) return(NULL)
        else taB <- reacIndices()[[input$selectPop3]]
        
        if(is.null(taB) | dim(taB)[1] == 0) return(NULL)
        datpl <- taB
        
        if(is.null(input$plotIndicesPG)) return(NULL)
        
        par(mar = rep(input$margin, 4))
        barplot(
          datpl[, input$plotIndicesPG],
          names.arg = datpl[, "locus"],
          horiz = TRUE,
          las = 1,
          col = transp(input$barplotcolor, input$transparency),
          border = as.numeric(input$borderbarplot),
          cex.axis = input$cexaxis, 
          cex.names = input$cexaxis,
          xlab = input$plotIndicesPG
        )
      })
      
      output$plotPG <- renderUI({
        plotOutput(ns('plotIndicesPopgen'))
      })

      # PAIRWISE FST -----------------------------------------------------------------
      
      #reactive function to compute FST matrix
      fstMatInput <- reactive({
        if (!input$displayFstMat)  return(NULL)
        dat2 <- getgenind()
        
        if (length(unique(dat2@pop)) == 1)
          stop("Multiple populations are required to perform this analysis")
        matFST <- pairwise.WCfst(genind2hierfstat(dat2))
        matFST[lower.tri(matFST)] <- NA
        matFST
      })
      
      #display FST matrix
      output$FstMat <- renderTable({
        if(!is.null(fstMatInput()))  matFST <- apply(fstMatInput(),1,rev)
      }, digits = 4,include.rownames =TRUE,na="")
      
      #DL fst matrix
      output$dlFstMat <- downloadHandler(
        filename = function() { 
          paste('pairwise_fst.tsv', sep='') 
        },
        content = function(file) {
          matFST <- apply(fstMatInput(),1,rev)
          write.table(matFST, file, sep="\t", na = "",row.names = TRUE)
        }
      )
      output$dlFstMatXL <- downloadHandler(
        filename = function() { 
          paste('pairwise_fst.xlsx', sep='') 
        },
        content = function(file) {
          matFST <- apply(fstMatInput(),1,rev)
          openxlsx::write.xlsx(list(fst=matFST), file, keepNA = FALSE,row.names = TRUE)
        }
      )
      # LD ---------------------------------------------------------------------------
      
      #compute LD table
      reacLDtable <- reactive({
        if (is.null(getgenind()))
          return(NULL)
        
        D <- getgenind()
        datLD <- genind2loci(D)
        loci <- (unique(D@loc.fac))
        nloc <- length(unique(D@loc.fac))
        LDmat <- matrix(NA, nrow = nloc, ncol = nloc)
        
        npairs <- nloc * (nloc - 1) / 2
        withProgress(message = 'LD computation', value = 0, {
          for(i in 1:nloc){
            for(ii in 1:nloc){
              if(i<ii){
                if(ploidy()=="Diploid") lx <- LD2(datLD,
                                                      locus=c(loci[i],loci[ii]))$T2
                if(ploidy()=="Haploid") lx <- LD(datLD,
                                                     locus=c(loci[i],loci[ii]))$T2
                LDmat[i,ii] <- lx[3]
                LDmat[ii,i] <- lx[1]
                incProgress(1/npairs,message="Computing LD...")
              } 
            }
          }
        })
        rownames(LDmat) <- colnames(LDmat) <- loci
        LDmat
      })
      
      output$LDtable <- renderTable({
        if (is.null(getgenind()) | !input$displayLDtable)
          return(NULL)
        
        M <- reacLDtable()
        M[lower.tri(M)] <- NA
        M <- apply(M,1,rev)
      }, digits = 4,include.rownames =TRUE,na="")
      
      output$LD30 <- reactive({
        if (is.null(getgenind()) | !input$displayLDtable) return(FALSE)
        M <- reacLDtable()
        return(length(M)>30)
      })
      outputOptions(output, 'LD30', suspendWhenHidden=FALSE)
      
      #plot heatmap LP p-values
      output$plotLD <- renderPlot({
        if (is.null(getgenind()) | !input$displayLDplot) return(NULL)
        
        M <-  -log10(reacLDtable())
        M[lower.tri(M)] <- NA
        col <- redpal(100)
        
        image(
          M,
          col = col,
          frame = F,
          xaxt = "n",
          yaxt = "n"
        )
        axis(
          2,
          at = seq(0, 1, length.out = ncol(M)),
          labels = colnames(M),
          las = 2,
          cex.axis = 0.8
        )
        axis(
          3,
          at = seq(0, 1, length.out = nrow(M)),
          labels = rownames(M),
          las = 2,
          cex.axis = 0.8
        )
        
      })
      
      output$plotLD2 <- renderUI({
        plotOutput(ns('plotLD'))
      })
      
      #plot p-values distribution
      output$plotLDpval <- renderPlot({
        if (is.null(getgenind()))
          return(NULL)
        
        M <-  (reacLDtable())
        
        par(mfrow=c(1,2))
        pv <- M[upper.tri(M)]
        
        hist(pv, xlab = "P-value", main = "LD p-values distribution")
        
        qqplot(
          pv,
          qunif(seq(0, 1, 0.01)),
          pch = 16,
          xlab = "Observed quantiles",
          ylab = "Expected quantiles",
          main = paste(
            "KS test p-value:",
            round(ks.test(pv, qunif(seq(0, 1, 0.01)))$p.value, digits = 4)
          )
        )
        abline(0, 1, lty = 2, lwd = 2, col = "grey")
        
      })
      
      output$plotLDpval2 <- renderUI({
        plotOutput(ns('plotLDpval'))
      })
      
      #DL LD matrix
      output$dlLDtable <- downloadHandler(
        
        filename = function() { 
          paste('LD_pvalues.tsv', sep='') 
        },
        
        content = function(file) {
          pairLD <- reacLDtable()
          pairLD[lower.tri(pairLD)] <- NA
          pairLD <- apply(pairLD, 1, rev)
          write.table(pairLD, file, sep = "\t", na = "", row.names = TRUE)
        }
      )
      output$dlLDtableXL <- downloadHandler(
        
        filename = function() { 
          paste('LD_pvalues.xlsx', sep='') 
        },
        
        content = function(file) {
          pairLD <- reacLDtable()
          pairLD[lower.tri(pairLD)] <- NA
          pairLD <- apply(pairLD, 1, rev)
          openxlsx::write.xlsx(list(LD=pairLD), file, keepNA = FALSE, row.names = TRUE)
        }
      )
    }
  )
}