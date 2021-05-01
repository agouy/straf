#' Generate the forensic parameters UI.
#' @export
#' @noRd
#' @importFrom shinyWidgets awesomeCheckbox
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

#' Generate the population genetics indices UI.
#' @export
#' @noRd
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
      condition = "input.ploidy == 2",
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

#' Generate the forensic parameters and population genetics tabs Server.
#' @export
#' @noRd
#' @importFrom  hierfstat pairwise.WCfst genind2hierfstat
for_popgen_Server <- function(id, getgenind, popnames, ploidy, barplotcolor, transparency, cexaxis) {
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
            taB[, !colnames(taB) %in% c("Ht", "Fis", "Fst")],
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
      output$plotIndices <-  plotly::renderPlotly({
        
        if(is.null(input$selectPop2)) return(NULL)
        else taB <- reacIndices()[[input$selectPop2]]
        
        if(is.null(taB) || dim(taB)[1] == 0) return(NULL)
        dat <- taB
        
        if(is.null(input$plotIndicesFOR)) return(NULL)
        
        fig <- plotly::plot_ly(
          x = dat[, input$plotIndicesFOR],
          y = dat[, "locus"],
          name = "bp_for",
          type = "bar",
          text = paste0(
            "Locus: ", dat[, "locus"],
            "\nValue: ", round(dat[, input$plotIndicesFOR], 4)
          ),
          hoverinfo = 'text',
          marker = list(color = adegenet::transp(barplotcolor(), transparency()))
        ) %>% plotly::layout(
          xaxis = list(title = input$plotIndicesFOR, zeroline = FALSE),
          yaxis = list(title = "Locus")
        )
        
        (fig)
      })
      
      output$plotFOR <- renderUI({
        plotly::plotlyOutput(ns('plotIndices'))
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
          col = transp(barplotcolor(), transparency()),
          border = 0,
          cex.axis = cexaxis(), 
          cex.names = cexaxis(),
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
        matFST <- hierfstat::pairwise.WCfst(hierfstat::genind2hierfstat(dat2))
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
        datLD <- pegas::genind2loci(D)
        loci <- (unique(D@loc.fac))
        nloc <- length(unique(D@loc.fac))
        LDmat <- matrix(NA, nrow = nloc, ncol = nloc)
        
        npairs <- nloc * (nloc - 1) / 2
        withProgress(message = 'LD computation', value = 0, {
          for(i in 1:nloc){
            for(ii in 1:nloc){
              if(i<ii){
                if(ploidy() == 2) lx <- pegas::LD2(datLD,
                                                      locus=c(loci[i],loci[ii]))$T2
                if(ploidy() == 1) lx <- pegas::LD(datLD,
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
      }, digits = 4, include.rownames =TRUE,na="")
      
      output$LD30 <- reactive({
        if (is.null(getgenind()) | !input$displayLDtable) return(FALSE)
        M <- reacLDtable()
        return(length(M) > 30)
      })
      outputOptions(output, 'LD30', suspendWhenHidden=FALSE)
      
      #plot heatmap LP p-values
      output$plotLD <- renderPlot({
        if (is.null(getgenind()) | !input$displayLDplot) return(NULL)
        
        M <-  -log10(reacLDtable())
        M[lower.tri(M)] <- NA
        col <- adegenet::redpal(100)
        
        image(M, col = col, frame = F, xaxt = "n", yaxt = "n")
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
      
      output$plotLD2 <- renderUI({plotOutput(ns('plotLD'))})
      
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

#' Get forensics and popgen indices for all populations
#' @inheritParams getIndicesFromGenind
#' @export
#' @noRd
getIndicesAllPop <- function(data, hw = FALSE, hwperm = 1000, ploidy = 2) {
  ind <- list()
  ind$all <- getIndicesFromGenind(data, hw, hwperm, ploidy)
  for(popu in unique(data$pop)) {
    ind <- c(ind, x = NA)
    mat <- getIndicesFromGenind(data[data@pop == popu, ], hw, hwperm, ploidy)
    ind$x <- mat
    names(ind)[length(ind)] <- popu
  }
  return(ind)
}

#' Get forensics and popgen indices for a given population
#' @param data A genind object.
#' @param hw boolean, whether to compute Hardy-Weinberg test p-values (default = FALSE).
#' @param hwperm numeric, number of permutations for Hardy-Weinberg test.
#' @param ploidy character, ploidy, either "Diploid" or "Haploid".
#' @export
#' @noRd
getIndicesFromGenind <- function(data,
                                 hw = FALSE,
                                 hwperm = 1000,
                                 ploidy = "Diploid") {
  
  freq <- apply(data@tab, 2, sum, na.rm = TRUE)
  
  nam <- strsplit(
    names(freq),
    split="[.]"
  )
  
  loc <- as.factor(unlist(
    lapply(nam, function(x) x[1])
  ))
  
  alle <- as.numeric(unlist(
    lapply(nam, function(x) sub("-", ".", x[2]))
  ))
  DAT <- data.frame(freq, loc, alle)
  N <- tapply(DAT$freq, DAT$loc, sum)
  DAT$frequency <- DAT$freq / N[DAT$loc]
  
  PIC <- rep(NA, length(unique(loc)))
  for(i in unique(loc)) {
    FR <- c(DAT$frequency[DAT$loc == i])
    p2q2 <- outer(FR, FR, get_2p2q2)
    sum_2p2q2 <- sum(p2q2[lower.tri(p2q2)])
    PIC[i] <-  1 - sum(FR ^ 2) - sum_2p2q2
  } 
  Nall <- tapply(
    DAT[DAT$freq > 0, ]$freq,
    DAT[DAT$freq > 0, ]$loc,
    length
  )
  
  GD <- tapply(DAT$frequency, DAT$loc, function(x) 1 - sum(x ^ 2))
  
  GD <- GD * N / (N - 1) 
  PIC <- PIC[names(GD)]
  D2 <- pegas::genind2loci(data)
  sumloc <- summary(D2)[names(GD)]
  PM1 <- lapply(sumloc, function(x) {
    sum((x$genotype / sum(x$genotype)) ^ 2)
  })
  PM <- unlist(PM1)
  
  DF <- data.frame(
    locus = names(GD),
    N = N,
    Nall = Nall,
    GD = GD,
    PIC = PIC,
    PM = PM,
    PD = 1 - PM
  )
  
  if(ploidy == 2) {
    DF$Hobs <- adegenet::summary(data)$Hobs[names(GD)]
    DF$PE <- (DF$Hobs ^ 2) * (1 - 2 * (DF$Hobs) * ((1 - DF$Hobs) ^ 2))
    DF$TPI <- 1 / (2 * (1 - DF$Hobs))
  }
  
  if(length(unique(data@pop)) > 1 & length(locNames(data)) > 1) {
    basicstat <- hierfstat::basic.stats(
      data,
      diploid = (ploidy == 2),
      digits = 4
    )$perloc
    rownames(basicstat) <- as.character(unique(data@loc.fac))
    Fst <- hierfstat::wc(
      data,
      diploid = (ploidy == 2)
    )$per.loc$FST
    names(Fst) <- as.character(unique(data@loc.fac))
    
    DF$Fst <- Fst[names(GD)]
    DF$Ht <- basicstat[names(GD), "Ht"]
    DF$Fis <- basicstat[names(GD), "Fis"]
  }
  
  if(ploidy == 2 & hw) {
    withProgress(message = 'Performing HW test...', value = 0, {
      DF$pHW <- pegas::hw.test(data, B = hwperm)[names(GD), 4]
    })
  } 
  return(DF)
}

#' 2p2q2
#' @param a first value
#' @param b second value
#' @return Computes 2 * a^2 * b^2
#' @export
#' @noRd
get_2p2q2 <- function(a, b) {2 * (a ^ 2) * (b ^ 2)}

