#' Generate the forensic parameters UI.
#' @export
#' @keywords internal
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
      div(tableOutput(ns('forensics')) %>% shinycssloaders::withSpinner(color="#dd4814"), style = "font-size:75%"),
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
#' @keywords internal
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
    
    uiOutput(ns("uiPG_HW")),
    
    conditionalPanel(
      condition = "input.ploidy == 2",
      ns = ns,
      uiOutput(ns("HW_comp"))
    ),
    conditionalPanel(
      condition = "input.displayDiv == true",
      ns = ns,
      div(tableOutput(ns('diversity')) %>% shinycssloaders::withSpinner(color="#dd4814"), style = "font-size:75%"),
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
      uiOutput(ns("selectPopLD")),
      div(tableOutput(ns('LDtable')) %>% shinycssloaders::withSpinner(color="#dd4814"), style = "font-size:75%"),
      downloadButton(ns('dlLDtable'), 'Download as text'),
      downloadButton(ns('dlLDtableXL'), 'Download as Excel'),
      uiOutput(ns("plotLD2"))
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
    tags$hr(),
    uiOutput(ns("haplo_UI"))
  )
}

#' Generate the forensic parameters and population genetics tabs Server.
#' @export
#' @keywords internal
for_popgen_Server <- function(id, input_file, getgenind, popnames, ploidy, hw_perm, barplotcolor, transparency, cexaxis) {
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
      output$selectPopLD <- renderUI({
        selectInput(ns("selectPopLD"), "Select a population:", popnames()[!popnames() %in% "all"])
      })
      output$uiPG_HW <- renderUI({
        pl <- ploidy()
        if(pl == 2) {
          return(list(
            numericInput(
              ns('hw_nperm'), 'Number of permutations for HW test',
              1000, min = 100, max = 10000, step = 100
            )
          ))
        } 
        else {
          return(NULL)
        } 
      })
      
      reacIndices <- reactive({
        if(is.null(getgenind())) return(NULL)
        
        DF <- getIndicesAllPop(getgenind(), ploidy())
        
        if(ploidy() == 2 & !is.null(input$hw_nperm)) {
          if(input$hw_nperm > 10000) stop("Number of permutations should not be greater than 10000.")
          hw_results <- getGenepopHW(input_file(), ploidy(), input$hw_nperm)
          for(pop in names(DF)) { # 
            df_tmp <- DF[[pop]]
            rownames(df_tmp) <- df_tmp$locus
            
            
            if(pop == "all") hw_tmp <- hw_results[hw_results$Population == pop, ]
            else {
              popu <- unique(hw_results$Population)[grep(pop, names(DF)[names(DF) != "all"])]
              hw_tmp <- hw_results[hw_results$Population == popu, ]
              hw_tmp$Population <- popu
            }
            
            hw_tmp$locus <- hw_tmp$Locus
            df_tmp <- dplyr::left_join(df_tmp, hw_tmp[, c("locus", "P_value")], "locus")
            
            DF[[pop]] <- df_tmp
          }
        }

        return(DF)
      })
      
      
      ### forensics table display
      output$forensics <- renderTable({
        req(reacIndices_pop_for())
        taB <- reacIndices_pop_for()
        
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
          req(reacIndices_pop_for())
          taB <- reacIndices_pop_for()
          write.table(
            taB[, !colnames(taB) %in% c("Ht", "Fis", "Fst")],
            file, sep = "\t", row.names = FALSE
          )
        }
      )
      output$dlForensicsXL <- downloadHandler(
        filename = function() { paste('forensics_parameters.xlsx', sep = '') },
        content = function(file) {
          req(reacIndices_pop_for())
          taB <- reacIndices_pop_for()
          out <- taB[, ! colnames(taB) %in% c("Ht", "Fis", "Fst")]
          openxlsx::write.xlsx(list(forensics_parameters = out), file, row.names = FALSE)
        }
      )
      
      ### popgen table display + download
      output$diversity <- renderTable({
        req(reacIndices_pop_pg())
        taB <- reacIndices_pop_pg()
        taB2 <- taB[, !colnames(taB) %in% c("PIC", "PD", "PE", "TPI", "PM")]
        taB2
      }, digits = 4)
      
      output$dlPopgen <- downloadHandler(
        filename = function() {
          paste('popgen_statistics.tsv', sep='') 
        },
        content = function(file) {
          req(reacIndices_pop_pg())
          taB <- reacIndices_pop_pg()
          taB2 <- taB[, !colnames(taB) %in% c("PIC","PD","PE","TPI", "PM")]
          write.table(taB2, file, sep = "\t", row.names = FALSE)
        }
      )
      output$dlPopgenXL <- downloadHandler(
        filename = function() {
          paste('popgen_statistics.xlsx', sep='') 
        },
        content = function(file) {
          req(reacIndices_pop_pg())
          taB <- reacIndices_pop_pg()
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
        taB <- reacIndices_pop_for()
        selectInput(
          ns("plotIndicesFOR"), "Select the statistic to plot:",
          choices = colnames(taB)[-1:-2],
          selected = "Nall"
        )
      })
      reacIndices_pop_for <- reactive({
        req(reacIndices())
        req(input$selectPop2)
        reacIndices()[[input$selectPop2]]
      })
      reacIndices_pop_pg <- reactive({
        req(reacIndices())
        req(input$selectPop3)
        reacIndices()[[input$selectPop3]]
      })
      
      # reactive UI: displayed if indices are computed
      output$uiPG <- renderUI({
        taB <- reacIndices_pop_for()
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
        
        fig
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
        datpl[, input$plotIndicesPG] <- as.numeric(datpl[, input$plotIndicesPG])
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
        if(!input$displayFstMat)  return(NULL)
        data <- getgenind()
        
        if(length(unique(data@pop)) == 1) {
          stop("Multiple populations are required to perform this analysis.")
        }
        fst_mat <- hierfstat::pairwise.WCfst(hierfstat::genind2hierfstat(data))
        fst_mat[lower.tri(fst_mat)] <- NA
        fst_mat <- apply(fst_mat, 1, rev)
        fst_mat
      })
      
      #display FST matrix
      output$FstMat <- renderTable({
        req(fstMatInput())
        fstMatInput()
      }, digits = 4, include.rownames = TRUE, na = "")
      
      #DL fst matrix
      output$dlFstMat <- downloadHandler(
        filename = function() { 
          paste('pairwise_fst.tsv', sep='')
        },
        content = function(file) {
          matFST <- fstMatInput()
          write.table(matFST, file, sep="\t", na = "", row.names = TRUE)
        }
      )
      output$dlFstMatXL <- downloadHandler(
        filename = function() { 
          paste('pairwise_fst.xlsx', sep='') 
        },
        content = function(file) {
          matFST <- fstMatInput()
          openxlsx::write.xlsx(list(fst = matFST), file, keepNA = FALSE, row.names = TRUE)
        }
      )
      # LD ---------------------------------------------------------------------------
      
      #compute LD table
      LD_genepop_table <- reactive({
        ld_results <- getGenepopLD(input_file(), ploidy(), 5000)
        return(ld_results)
      })

      LD_genepop_pairwise_pop <- reactive({
        req(LD_genepop_table_pop())
        ld_results_sub <- LD_genepop_table_pop()

        ld_results_sub <- ld_results_sub[, c("Locus_1", "Locus_2", "P_value")]
        ld_results_sub <- tidyr::pivot_wider(
          ld_results_sub, 
          names_from = "Locus_2", 
          values_from  = "P_value"
        )
        ld_results_sub <- as.data.frame(ld_results_sub)
        rownames(ld_results_sub) <- ld_results_sub$Locus_1
        ld_results_sub$Locus_1 <- NULL
        
        ld_results_sub <- ld_results_sub[order(rownames(ld_results_sub)), ]
        ld_results_sub <- ld_results_sub[, order(colnames(ld_results_sub), decreasing = TRUE)]
        
        return(ld_results_sub)
        
      })        
      LD_genepop_table_pop <- reactive({
        req(LD_genepop_table())
        req(input$selectPopLD)

        pop <- input$selectPopLD
        ld_results <- LD_genepop_table()
        
        popu <- unique(ld_results$Population)[grep(pop, popnames()[popnames() != "all"])]
        
        ld_results_sub <- ld_results[ld_results$Population == popu, ]
        ld_results_sub$Population <- pop
        
        ld_results_sub <- ld_results_sub[, c("Locus_1", "Locus_2", "P_value", "plab")]
        return(ld_results_sub)
      })

            
      reacLDtable <- reactive({
        if (is.null(getgenind()))
          return(NULL)
        LD_genepop_pairwise_pop()
      })
      
      output$LDtable <- renderTable({
        if (is.null(getgenind()) | !input$displayLDtable)
          return(NULL)
        LD_genepop_pairwise_pop()
      }, digits = 4, include.rownames = TRUE,na="")

      #plot heatmap LP p-values
      output$plotLD <- plotly::renderPlotly({
        df <- LD_genepop_table_pop()
        plt <- plot_ld(df)
        plt
      })
      # 
      output$plotLD2 <- renderUI({plotly::plotlyOutput(ns('plotLD'))})

      #DL LD matrix
      output$dlLDtable <- downloadHandler(
        
        filename = function() { 
          paste('LD_pvalues.tsv', sep='') 
        },
        
        content = function(file) {
          pairLD <- LD_genepop_pairwise_pop()
          write.table(pairLD, file, sep = "\t", na = "", row.names = TRUE)
        }
      )
      output$dlLDtableXL <- downloadHandler(
        
        filename = function() { 
          paste('LD_pvalues.xlsx', sep='') 
        },
        
        content = function(file) {
          pairLD <- LD_genepop_pairwise_pop()
          openxlsx::write.xlsx(list(LD=pairLD), file, keepNA = FALSE, row.names = TRUE)
        }
      )
      
      
      output$haplo_UI <- renderUI({
        if(ploidy() == 2) {
            list()
        } else {
          list(
            h3("Haploid data statistics"),
            uiOutput(ns("selectPopHaplo")),
            tableOutput(ns("haplo.stats")),
            tags$hr()
          )
        }
      })

      output$selectPopHaplo <- renderUI({
        selectInput(ns("selectPopHaplo"), "Select a population:", popnames())
      })
      
      haploStats_pop <- reactive({
        req(haploStats())
        req(input$selectPopHaplo)
        haploStats()[[input$selectPopHaplo]]$hap_data
      })
      
      haploStats <- reactive({
        data <- getgenind()
        req(data)
        df <- getHaploStatsAllPop(data)
        return(df)
      })
      
      output$haplo.stats <- renderTable({
        req(haploStats_pop())
        haplo_stats <- haploStats_pop()
        haplo_stats
      }, digits = 4)
      
    }
  )
}

#' Get forensics and popgen indices for all populations
#' @inheritParams getIndicesFromGenind
#' @export
#' @keywords internal
getIndicesAllPop <- function(data, ploidy = 2) {
  ind <- list()
  ind$all <- getIndicesFromGenind(data,ploidy)
  for(popu in unique(data$pop)) {
    ind <- c(ind, x = NA)
    mat <- getIndicesFromGenind(data[data@pop == popu, ],ploidy)
    ind$x <- mat
    names(ind)[length(ind)] <- popu
  }
  return(ind)
}

#' Get forensics and popgen indices for a given population
#' @param data A genind object.
#' @param ploidy character, ploidy, either "Diploid" or "Haploid".
#' @export
#' @keywords internal
getIndicesFromGenind <- function(data, ploidy = 2) {
  
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
  
  ploidy <- as.numeric(ploidy)
  return(DF)
}

#' 2p2q2
#' @param a first value
#' @param b second value
#' @return Computes 2 * a^2 * b^2
#' @export
#' @keywords internal
get_2p2q2 <- function(a, b) {2 * (a ^ 2) * (b ^ 2)}

