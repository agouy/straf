options(stringsAsFactors = FALSE)
source("./scripts/helpers.R")
load("./www/strider_frequencies.rda")

shinyServer(function(input, output) {

  getData <- reactive({
    ### check the input: returns NULL if problem
    ### if NULL, the function checkinputfile will print a message
    if(is.null(input$file1)) { return(NULL) } else {
      
      X <- read.table(
        input$file1$datapath,
        header = TRUE,
        sep = "\t",
        colClasses = "character"
      )
      
      if(class(X) != "data.frame" |
         colnames(X)[1] != "ind" |
         colnames(X)[2] != "pop" |
         dim(X)[1] < 2 |
         dim(X)[2] < 3) {
        
        warning("Input file error.")
        return(NULL)
        
      }
      
      testGeno <- try(getgenind(), silent = TRUE)

      if(class(testGeno) == "try-error") {
        warning("Input file error.")
        return(NULL)
      }
      
      if(length(unique(testGeno@pop)) > 1 & length(locNames(testGeno)) > 1) {
        
        testGeno3 <- try(
          wc(
            testGeno,
            diploid = switch(
              input$ploidy,
              Diploid = TRUE,
              Haploid = FALSE
            )
          ),
          silent = TRUE
        )
        if(class(testGeno3) == "try-error") {
          warning("Input file error.")
          return(NULL)
        }
      }
      return(input$file1)
    }
  })
  
  output$fileUploaded <- reactive({ return(!is.null(getData())) })
  
  output$checkInputFile <- renderText({ 
    
    if(is.null(input$file1)) {
      
      return("Please import a data set.")
      
    } else {
      
      X <- read.table(input$file1$datapath,
                      header = TRUE,
                      sep = "\t",
                      colClasses="character")
      
      if(class(X) != "data.frame") {
        return("Input file is not a data frame. Please read the documentation for more information.")
      }
      if(colnames(X)[1] != "ind" | colnames(X)[2] != "pop") {
        return("Input file error. Please read the documentation for more information.")
      }
      if(dim(X)[1] < 2) {
        return("Input file error: incorrect number of rows. Please read the documentation for more information.")
      }
      if(dim(X)[2] < 3) {
        return("Input file error: incorrect number of columns. Please read the documentation for more information.")
      } 
      
      testGeno <- try(getgenind(), silent = TRUE)
      
      if(class(testGeno) == "try-error") {
        return("Input file error. Wrong number of columns per locus (please check ploidy and number of digits).")
      }

      testGeno3 <- try(
        wc(
          testGeno,
          diploid = switch(
            input$ploidy,
            Diploid = TRUE,
            Haploid = FALSE
          )
        ),
        silent = TRUE
      )
      
      if(class(testGeno3) == "try-error") {
        return("Input file error. Wrong ploidy (please check the left tab).")
      }
      
    }
    
  })
  
  outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)
  #makes the fileUploaded variable "more reactive"
  
  popnames <- reactive({
    
    if(is.null(input$file1)) { 
      return("") 
    } else {
      X <- read.table(
        input$file1$datapath,
        header = TRUE,
        sep = "\t",
        colClasses = "character"
      )
      pop.names <- c("all", unique(X$pop))
      return(pop.names)
    }
  })
  
  # DISPLAY DATASET --------------------------------------------------------------
  output$contents <- renderDataTable({
    if (is.null(getData())) return(NULL)
    X <- read.table(
      getData()$datapath,
      header = TRUE, sep = "\t",
      colClasses = "character"
    )
    DT::datatable(X)
  })
  
  # ALLELE FREQUENCIES -----------------------------------------------------------
  
  output$alleleFreq <- renderPlot({ #barplots of allele freq distributions
    
    if (!input$displayAlleleFreq)  return(NULL)
    dat2 <- getgenind()##ts
    
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
        col = transp(input$barplotcolor, input$transparency),
        border = as.numeric(input$borderbarplot)
      )
    }
  })
  
  output$plotAF <- renderUI({ #display UI only if allele freq is checked
    
    plotOutput(
      'alleleFreq',
      width = paste(input$width, "%", sep = ""),
      height = input$height
    )
    
  })
  
  # ALLELE FREQUENCIES TABLE -----------------------------------------------------
  
  alleleFreqTabs <- reactive({
    if (!input$displayAlleleTable)  return(NULL)
    dat2 <- getgenind()
    matr <- getFreqAllPop(dat2)
    return(matr)
  })
  
  output$selectPop <- renderUI({
    selectInput("selectPop", "Select a population:", popnames())
  })
  
  output$tableFreq <- renderDataTable({
    
    if (!input$displayAlleleTable) return(NULL)
    
    if(input$selectPop == "") matr <- alleleFreqTabs()[[1]]
    else matr <- alleleFreqTabs()[[input$selectPop]]
    
    DT::datatable(matr) %>% formatRound(columns = colnames(matr), digits = 3)
    
  })
  
  output$dlTabfreq <- downloadHandler(
    
    filename = function() { 
      paste('allele_frequencies.tsv', sep='') 
    },
    
    content = function(file) {
      
      if (!input$displayAlleleTable) return(NULL)
      
      if(input$selectPop == "") matr <- alleleFreqTabs()[[1]]
      else matr <- alleleFreqTabs()[[input$selectPop]]
      
      write.table(matr, file, sep = "\t", na = "", row.names = TRUE)
    }
  )
  output$dlTabfreqXL <- downloadHandler(
    
    filename = function() { 
      paste('allele_frequencies.xlsx', sep='') 
    },
    
    content = function(file) {
      
      if (!input$displayAlleleTable) return(NULL)
      
      if(input$selectPop == "") matr <- alleleFreqTabs()[[1]]
      else matr <- alleleFreqTabs()[[input$selectPop]]
      
      openxlsx::write.xlsx(list(allele_frequencies = matr), file = file, rowNames = TRUE)
    }
  )
  
  
  # STATISTICS TABLE -------------------------------------------------------------
  
  ### compute forensics indices, reactive function
  output$selectPop2 <- renderUI({
    
    selectInput("selectPop2", "Select a population:", popnames())
    
  })
  
  output$selectPop3 <- renderUI({
    
    selectInput("selectPop3", "Select a population:", popnames())
    
  })
  
  reacIndices <- reactive({
    inFile <- input$file1
    if(is.null(inFile)) return(NULL)

    DF <- getIndicesAllPop(
      getgenind(),
      hw = input$computeHW,
      hwperm = input$hw_nperm,
      ploidy = input$ploidy
    )
    return(DF)
  })
  
  ### forensics table display
  output$forensics <- renderTable({
    
    if (is.null(input$file1)) return(NULL)
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
      write.table(
        taB2,
        file,
        sep = "\t",
        row.names = FALSE
      )
      
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
    
    selectInput("plotIndicesFOR", "Select the statistic to plot:",
                choices = colnames(taB)[-1:-2],
                selected = "Nall"
    )
  })
  
  # reactive UI: displayed if indices are computed
  output$uiPG <- renderUI({
    if(is.null(input$selectPop3)) return(NULL)
    else taB <- reacIndices()[[input$selectPop3]]
    
    selectInput("plotIndicesPG", "Select the statistic to plot:",
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
    plotlyOutput('plotIndices', width=paste(input$width,"%",sep=""),
                 height=input$height)
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
    plotOutput('plotIndicesPopgen', width=paste(input$width,"%",sep=""),
               height=input$height)
  })
  
  # create genind function
  
  getgenind <- reactive({
    inFile <- input$file1
    
    if(is.null(inFile)) return(NULL)
    
    df_out <- createGenind(
      Ifile = inFile,
      Imicrovariants = input$microvariants,
      Incode = input$ncode,
      Iploidy = input$ploidy
    )
    return(df_out)
  })
  # create genind function
  
  # getgenind4ld <- reactive({
  #   inFile <- input$file1
  #   
  #   if(is.null(inFile)) return(NULL)
  #   
  #   df_out <- genind4LD(
  #     Ifile = inFile,
  #     Imicrovariants = input$microvariants,
  #     Incode = input$ncode,
  #     Iploidy = input$ploidy
  #   )
  #   return(df_out)
  # })
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
    inFile <- input$file1
    
    if (is.null(inFile))
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
            if(input$ploidy=="Diploid") lx <- LD2(datLD,
                                                  locus=c(loci[i],loci[ii]))$T2
            if(input$ploidy=="Haploid") lx <- LD(datLD,
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
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayLDtable)
      return(NULL)
    
    M <- reacLDtable()
    M[lower.tri(M)] <- NA
    M <- apply(M,1,rev)
  }, digits = 4,include.rownames =TRUE,na="")
  
  output$LD30 <- reactive({
    if (is.null(input$file1) | !input$displayLDtable) return(FALSE)
    M <- reacLDtable()
    return(length(M)>30)
  })
  outputOptions(output, 'LD30', suspendWhenHidden=FALSE)
  
  #plot heatmap LP p-values
  output$plotLD <- renderPlot({
    
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayLDplot) return(NULL)
    
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
    plotOutput(
      'plotLD',
      width = paste(input$width, "%", sep = ""),
      height = input$height
    )
  })
  
  #plot p-values distribution
  output$plotLDpval <- renderPlot({
    inFile <- input$file1
    
    if (is.null(inFile))
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
    plotOutput(
      'plotLDpval',
      width = paste(input$width, "%", sep = ""),
      height = input$height)
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
  # PCA --------------------------------------------------------------------------
  
  output$runPCA <- renderPlot({
    if (!input$displayPCA)  return(NULL)
    dat2 <- getgenind()
    
    FREQ <- makefreq(dat2,missing="mean")
    pcaY <- dudi.pca(FREQ, nf = 3, scannf = FALSE)
    
    
    if(length(input$PCAaxis)==2){
      par(mfrow=c(1,1))
      coul <- transp(funky(length(unique(pop(dat2)))),.6)
      plotPCA(
        pca = pcaY,
        popus = pop(dat2),
        coul = coul,
        axis = c(as.numeric(input$PCAaxis[1]), as.numeric(input$PCAaxis[2]))
      )
    }
    if(length(input$PCAaxis)==3){
      
      par(mfrow = c(1,2))
      coul <- transp(funky(length(unique(pop(dat2)))),.6)
      plotPCA(
        pca = pcaY,
        popus = pop(dat2),
        coul = coul,
        axis = c(as.numeric(input$PCAaxis[1]), as.numeric(input$PCAaxis[2]))
      )
      plotPCA(
        pca = pcaY,
        popus = pop(dat2),
        coul = coul,
        axis = c(as.numeric(input$PCAaxis[1]), as.numeric(input$PCAaxis[3]))
      )
      
    }
  })
  
  output$plotPCA <- renderUI({
    plotOutput('runPCA', width = paste(input$width,"%",sep = ""), 
               height = input$height, click = "plot_click")
  })
  output$plotMDS <- renderUI({
    plotOutput('runMDS', width = paste(input$width,"%",sep = ""), 
               height = input$height)
  })
  output$runMDS <- renderPlot({
    if (!input$displayMDS)  return(NULL)
    
    dat2 <- getgenind()
    
    if(length(levels(pop(dat2))) < 2) return(NULL)
    
    obj <- genind2genpop(dat2, quiet = FALSE)
    dst <- dist.genpop(obj, method = 1)
    MDS <- cmdscale(dst)
    MDS <- data.frame(ax1 = MDS[, 1], ax2 = MDS[, 2], pop = rownames(MDS))
    
    p <- ggplot(MDS, aes(x=ax1, y=ax2, color = pop, label = pop)) +
      geom_point() +
      geom_text_repel() + 
      labs( x = "MDS Axis 1", y = "MDS Axis 2", title = "MDS based on Nei's distance")  +
      theme_minimal()
    plot(p)
  })
  
  output$plotMDS_strider <- renderUI({
    plotOutput('runMDS_strider', width = paste(input$width,"%",sep = ""), 
               height = input$height)
  })
  common_alleles <- reactive({
    X <- strider_frequencies[input$location, ]
    X <- X[, colSums(is.na(X)) == 0]
    dat2 <- getgenind()
    obj <- genind2genpop(dat2, quiet = FALSE)
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
  output$runMDS_strider <- renderPlot({
    #if (!input$displayMDS)  return(NULL)
    
    X <- strider_frequencies[input$location, ]
    X <- X[, colSums(is.na(X)) == 0]
    
    if(input$add_current) {
      X <- common_alleles()
    }
    
    d <- X %*% t(X)
    vec <- sqrt(diag(d))
    d <- d/vec[col(d)]
    d <- d/vec[row(d)]
    d <- -log(d)
    d <- as.dist(d)
    
    mds <- cmdscale(d)
    MDS <- data.frame(ax1 = mds[, 1], ax2 = mds[, 2], pop = rownames(mds))
    
    p <- ggplot(MDS, aes(x=ax1, y=ax2, color = pop, label = pop)) +
      geom_point() +
      geom_text_repel(max.overlaps = 50) + 
      labs( x = "MDS Axis 1", y = "MDS Axis 2", title = "MDS based on Nei's distance")  +
      theme_minimal()
    plot(p)
    
  })  
  
  # DL principal components
  output$dlPCAeigen <- downloadHandler(
    filename = function() {
      paste('PCA_eigenvectors.tsv', sep = '')
    },
    content = function(file) {
      
      if (!input$displayPCA)  return(NULL)
      dat2 <- getgenind()
      
      FREQ <- makefreq(dat2,missing = "mean")
      pcaY <- dudi.pca(FREQ, nf = 3, scannf = FALSE)
      write.table(pcaY$c1, file, sep = "\t", na = "",row.names = TRUE)
    }
  )
  
  output$dlPCAcoord <- downloadHandler(
    filename = function() {
      paste('PCA_coordinates.tsv', sep = '')
    },
    content = function(file) {
      if (!input$displayPCA)  return(NULL)
      dat2 <- getgenind()
      FREQ <- makefreq(dat2,missing = "mean")
      pcaY <- dudi.pca(FREQ, nf = 3, scannf = FALSE)
      write.table(pcaY$li, file, sep = "\t", na = "",row.names = TRUE)
    }
  )
  
  # plot eigen values
  
  output$info <- renderPrint({
    if (!input$displayPCA)  return(NULL)
    dat2 <- getgenind()
    
    FREQ <- makefreq(dat2,missing = "mean",quiet = TRUE)
    pcaY <- dudi.pca(FREQ, nf = 3, scannf = FALSE)
    
    if(length(input$PCAaxis) == 2){
      cat("Click on a point to get its ID and coordinates\n\n")
      
      ta <- c("Axis1","Axis2","Axis3")
      
      if(!is.null(input$plot_click)){
        nearPoints(pcaY$li, input$plot_click, 
                   xvar = ta[as.numeric(input$PCAaxis[1])], 
                   yvar = ta[as.numeric(input$PCAaxis[2])])
      }
    }
  })
  
  output$loadings <- renderPlot({
    if (!input$displayPCA)  return(NULL)
    dat2 <- getgenind()
    FREQ <- makefreq(dat2,missing = "mean")
    pcaY <- dudi.pca(FREQ, nf = 3, scannf = FALSE)
    loadingplot(pcaY$c1^2)
  })
  output$plotLoadings <- renderUI({
    plotOutput('loadings', width = paste(input$width,"%",sep = ""), height = input$height)
  })
  
  
  #### FILE CONVERSION
  output$dlGenepop <- downloadHandler(
    filename = function() { 
      paste('straf2genepop.txt', sep='') 
    },
    content = function(file) {
      gp <- straf2genepop(f.name = input$file1$datapath, ploidy = switch(input$ploidy, Diploid = 2, Haploid = 1))
      cat(gp, file = file)
    }
  )
  
  output$dlArlequin <- downloadHandler(
    filename = function() { 
      paste('straf2arlequin.arp', sep='') 
    },
    content = function(file) {
      gp <- straf2arlequin(input$file1$datapath)
      cat(gp, file = file)
    }
  )
  
  getInput <- reactive(input$file1$datapath)
  
  output$dlFamilias <- downloadHandler(
    filename = function() { 
      paste('straf2familias.txt', sep='')
    },
    content = function(file) {
      fmi <- straf2familias(getInput())
      cat(fmi, file = file)
    }
  )
  
})
