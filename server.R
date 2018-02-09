# SERVER -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)
  library(ade4)
  library(adegenet)
  library(pegas)
  library(hierfstat)
  library(DT)
  library(car)
})

shinyServer(function(input, output) {

# GET INPUT FILE -------------------------------------------------------------
    
  getData <- reactive({ #get file from the input
    
    ### check the input: returns NULL if problem
    ### if NULL, the function checkinputfile will print a message
    if(is.null(input$file1)) {
      
      return(NULL)
      
    } else {
      
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
      
      testGeno <- try(
        createGenind(
          Ifile = input$file1,
          Imicrovariants = input$microvariants,
          Incode = input$ncode,
          Iploidy = input$ploidy
        ),
        silent = TRUE
      )
      
      testGeno2 <- try(
        genind4LD(
          Ifile = input$file1,
          Imicrovariants = input$microvariants,
          Incode = input$ncode,
          Iploidy = input$ploidy
        ),
        silent = TRUE
      )
      
      if(class(testGeno) == "try-error" | class(testGeno2) == "try-error") {
        
        warning("Input file error.")
        return(NULL)
        
      }

      if(length(unique(testGeno2@pop)) > 1 & length(locNames(testGeno2)) > 1) {

        testGeno3 <- try(
          wc(
            testGeno2,
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
  
  output$fileUploaded <- reactive({
    
    return(!is.null(getData()))
    
  })
  
  output$checkInputFile <- renderText({ 
    
    if(is.null(input$file1)) {
      
      return("Please import a data set.")
      
    } else {
      
      X <- read.table(input$file1$datapath,
                      header = TRUE,
                      sep = "\t",
                      colClasses="character")
      
      if(class(X) != "data.frame") {
        return("Input file is not a data frame. Please read the documentation 
               for more information.")
      }
      if(colnames(X)[1] != "ind" | colnames(X)[2] != "pop") {
        return("Input file error. Please read the documentation for more 
               information.")
      }
      if(dim(X)[1] < 2) {
        return("Input file error: incorrect number of rows. Please read the 
               documentation for more information.")
      }
      if(dim(X)[2] < 3) {
        return("Input file error: incorrect number of columns. Please read 
               the documentation for more information.")
      } 

      testGeno <- try(
        createGenind(
          Ifile=input$file1,
          Imicrovariants=input$microvariants,
          Incode=input$ncode,
          Iploidy=input$ploidy
        ),
        silent=TRUE)
      
      if(class(testGeno) == "try-error") {
        return("Input file error. Wrong number of columns per locus 
               (please check the left tab).")
      }
      
      testGeno2 <- try(
        genind4LD(
          Ifile = input$file1,
          Imicrovariants = input$microvariants,
          Incode = input$ncode,
          Iploidy = input$ploidy
        ),
        silent = TRUE
      )

      if(class(testGeno2) == "try-error") {
        return("Input file error. Wrong ploidy or number of digits (please check
               the left tab).")
      }

      testGeno3 <- try(
        wc(
          testGeno2,
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
    
    if (is.null(getData()) | !input$displayTable) return(NULL)
    
    X <- read.table(
      getData()$datapath,
      header = TRUE,
      sep = "\t",
      colClasses = "character"
    )
    
    DT::datatable(X) #datatable returns an interactive table
    
  })
  
# ALLELE FREQUENCIES -----------------------------------------------------------
  
  output$alleleFreq <- renderPlot({ #barplots of allele freq distributions
    
    if (is.null(getData()) | !input$displayAlleleFreq)  return(NULL)
    
    dat2 <- createGenind(
      Ifile = getData(),
      Imicrovariants = input$microvariants,
      Incode = input$ncode,
      Iploidy = input$ploidy
    )
    
    freq <- apply(dat2@tab, 2, sum, na.rm = TRUE) #counts number of alleles
    
    nam <- strsplit(
      names(freq),
      split = "[.]"
    ) #split locus and allele name
    
    loc <- as.factor(unlist(
      lapply(nam, function(x) x[1])
    ))
    alle <- as.numeric(unlist(
      lapply(nam, function(x) sub("-", ".", x[2]))
    ))
    
    DAT <- data.frame(freq,loc,alle)
    DAT <- DAT[order(DAT$alle),]
    
    ###depending on the number of loci, different number of columns:
    nL <- length(unique(DAT$loc))
    if(nL <= 5) ncolo <- 2
    if(nL >= 6) ncolo <- 3
    if(nL >= 10) ncolo <- 4
    par(mfrow = c(ceiling(nL / ncolo), ncolo))
    
    par(mar=c(2,2,2,2))
    for(i in unique(DAT$loc)) {
      
      barplot(
        DAT$freq[DAT$loc == i],
        names.arg = DAT$alle[DAT$loc == i],
        main = i,
        col = transp(
          input$barplotcolor,
          input$transparency
        ),
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
    
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayAlleleTable) return(NULL)
    
    dat2 <- createGenind(
      Ifile = inFile,
      Imicrovariants = input$microvariants,
      Incode = input$ncode,
      Iploidy = input$ploidy
    )

    matr <- getFreqAllPop(dat2)
    
    return(matr)
    
  })

  output$selectPop <- renderUI({
    
    selectInput("selectPop", "Select a population:", popnames())
    
  })

  output$tableFreq <- renderDataTable({
    
    inFile <- input$file1

    if (is.null(inFile) | !input$displayAlleleTable)
      return(NULL)
    
    if(input$selectPop == "") matr <- alleleFreqTabs()[[1]]
    else matr <- alleleFreqTabs()[[input$selectPop]]
    
    DT::datatable(matr) %>% formatRound(columns = colnames(matr), digits = 3)
    
  })

  output$dlTabfreq <- downloadHandler(
    
    filename = function() { 
      paste('allele_frequencies.tsv', sep='') 
    },
    
    content = function(file) {
      inFile <- input$file1
      
      if (is.null(inFile) | !input$displayAlleleTable)
        return(NULL)
      
      if(input$selectPop == "") matr <- alleleFreqTabs()[[1]]
      else matr <- alleleFreqTabs()[[input$selectPop]]
      
      write.table(matr, file, sep = "\t", na = "", row.names = TRUE)
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
    
    if (is.null(inFile))
      return(NULL)
    
    dat2 <- genind4LD(
      Ifile = inFile,
      Imicrovariants = input$microvariants,
      Incode = input$ncode,
      Iploidy = input$ploidy
    )

    hw <- input$computeHW
    ploidy <- input$ploidy

    DF <- getIndicesAllPop(
      dat2,
      hw = hw,
      hwperm = input$hw_nperm,
      ploidy = ploidy
    )
    return(DF)
    
  })
  
  ### forensics table display
  output$forensics <- renderTable({
    
    inFile <- input$file1
    
    if (is.null(inFile)) return(NULL)
    
    dat2 <- genind4LD(
      Ifile = inFile,
      Imicrovariants = input$microvariants,
      Incode = input$ncode,
      Iploidy = input$ploidy
    )

    if(is.null(input$selectPop2)) taB <- reacIndices()[[1]]
    else taB <- reacIndices()[[input$selectPop2]]

    if(! is.null(input$selectPop2)) {
      
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
    
    filename = function() { 
      
      paste('forensics_parameters.tsv', sep = '') 
      
    },
    
    content = function(file) {
      
      if(is.null(input$selectPop2)) taB <- reacIndices()[[1]]
      else taB <- reacIndices()[[input$selectPop2]]

      write.table(
        taB[, ! colnames(taB) %in% c("Ht", "Fis", "Fst")],
        file,
        sep = "\t",
        row.names = FALSE
      )
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
  output$plotIndices <-  renderPlot({
    
    if(is.null(input$selectPop2)) return(NULL)
    else taB <- reacIndices()[[input$selectPop2]]
    
    if(is.null(taB) | dim(taB)[1]==0) return(NULL)
    dat <- taB

    if(is.null(input$plotIndicesFOR)) return(NULL)
    
    par(mar=rep(input$margin,4))
    barplot(
      dat[,input$plotIndicesFOR],
      names.arg = dat[,"locus"],
      horiz = TRUE,
      las = 1,
      col = transp(input$barplotcolor, input$transparency),
      border = as.numeric(input$borderbarplot),
      cex.axis = input$cexaxis,
      cex.names = input$cexaxis,
      xlab = input$plotIndicesFOR
    )
  })
  
  output$plotFOR <- renderUI({
    plotOutput('plotIndices', width=paste(input$width,"%",sep=""),
               height=input$height)
  })
  
  output$plotIndicesPopgen <-  renderPlot({
    
    if(is.null(input$selectPop3)) return(NULL)
    else taB <- reacIndices()[[input$selectPop3]]
    
    if(is.null(taB) | dim(taB)[1]==0) return(NULL)
    datpl <- taB
    
    if(is.null(input$plotIndicesPG)) return(NULL)
    
    par(mar=rep(input$margin,4))
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
  
  
# PAIRWISE FST -----------------------------------------------------------------
  
  #reactive function to compute FST matrix
  fstMatInput <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayFstMat)
      return(NULL)
    
    dat2 <- createGenind(Ifile=inFile,Imicrovariants=input$microvariants,
                       Incode=input$ncode,Iploidy=input$ploidy)
    
    if (length(unique(dat2@pop)) == 1)
      stop("Multiple populations are required to perform this analysis")
    matFST <- pairwise.fst(dat2, res.type = "matrix")
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
  
# LD ---------------------------------------------------------------------------
  
  #compute LD table
  reacLDtable <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    D <- genind4LD(Ifile=inFile,Imicrovariants=input$microvariants,
                 Incode=input$ncode,Iploidy=input$ploidy)
    
    datLD <- genind2loci(D)
    loci <- (unique(D@loc.fac))
    nloc <- length(unique(D@loc.fac))
    LDmat <- matrix(NA,nrow=nloc,ncol=nloc)
    
    npairs <- nloc*(nloc-1)/2
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
  
# PCA --------------------------------------------------------------------------
  
  output$runPCA <- renderPlot({
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayPCA)
      return(NULL)
    
    dat2 <- createGenind(Ifile=inFile,Imicrovariants=input$microvariants,
                       Incode=input$ncode,Iploidy=input$ploidy)
    
    FREQ <- makefreq(dat2,missing="mean")
    pcaY <- dudi.pca(FREQ, nf = 3, scannf = FALSE)
    
    
    if(length(input$PCAaxis)==2){
      par(mfrow=c(1,1))
      # s.class(pcaY$li, fac = pop(dat2),
      #         col = transp(funky(length(unique(pop(dat2)))),.6),
      #         axesel = FALSE, cstar = 0, cpoint = 2,xax = as.numeric(input$PCAaxis[1]),
      #         yax = as.numeric(input$PCAaxis[2]))
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
      # s.class(pcaY$li, fac = pop(dat2),
      #         col = transp(funky(length(unique(pop(dat2)))),.6),
      #         axesel = FALSE, cstar = 0, cpoint = 2,xax = as.numeric(input$PCAaxis[1]),
      #         yax = as.numeric(input$PCAaxis[2]))
      # s.class(pcaY$li, fac = pop(dat2),
      #         col = transp(funky(length(unique(pop(dat2)))),.6),
      #         axesel = FALSE, cstar = 0, cpoint = 2,xax = as.numeric(input$PCAaxis[1]),
      #         yax = as.numeric(input$PCAaxis[3]))
    }
  })
  
  output$plotPCA <- renderUI({
    plotOutput('runPCA', width = paste(input$width,"%",sep = ""), 
               height = input$height, click = "plot_click")
  })
  
  # DL principal components
  output$dlPCAeigen <- downloadHandler(
      filename = function() {
        paste('PCA_eigenvectors.tsv', sep = '')
      },
      content = function(file) {
    
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayPCA)
      return(NULL)
    
    dat2 <- createGenind(Ifile = inFile,Imicrovariants = input$microvariants,
                         Incode = input$ncode,Iploidy = input$ploidy)
  
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
      
      inFile <- input$file1
      
      if (is.null(inFile) | !input$displayPCA)
        return(NULL)
      
      dat2 <- createGenind(Ifile = inFile,Imicrovariants = input$microvariants,
                           Incode = input$ncode,Iploidy = input$ploidy)
      
      FREQ <- makefreq(dat2,missing = "mean")
      pcaY <- dudi.pca(FREQ, nf = 3, scannf = FALSE)
      write.table(pcaY$li, file, sep = "\t", na = "",row.names = TRUE)
    }
  )
  
  # plot eigen values
  
  output$info <- renderPrint({
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayPCA & length(input$PCAaxis)==2)
      return(NULL)
    
    dat2 <- createGenind(Ifile = inFile,Imicrovariants = input$microvariants,
                       Incode = input$ncode,Iploidy = input$ploidy)
    
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
    inFile <- input$file1
    if(is.null(inFile) | !input$displayPCA) return(NULL)
    
    dat2 <- createGenind(Ifile = inFile,Imicrovariants = input$microvariants,
                       Incode = input$ncode,Iploidy = input$ploidy)
    FREQ <- makefreq(dat2,missing = "mean")
    pcaY <- dudi.pca(FREQ, nf = 3, scannf = FALSE)
    loadingplot(pcaY$c1^2)
  })
  output$plotLoadings <- renderUI({
    plotOutput('loadings', width = paste(input$width,"%",sep = ""), height = input$height)
  })
})

# EXTERNAL FUNCTIONS -----------------------------------------------------------

### convert input to genind object
createGenind <- function(Ifile, Imicrovariants, Incode, Iploidy) {
  
  if(Imicrovariants == 2) {
    
    readLines(Ifile$datapath) -> matrix
    matrix <- strsplit(matrix,"[\t]")
    
    mat <- matrix(unlist(matrix), nrow = length(matrix), ncol = length(matrix[[1]]),
                byrow = TRUE)
    mat[mat == "0"] <- NA ###
    colnames(mat) <- mat[1,]
    rownames(mat) <- mat[,1]
    
    mat <- mat[-1,]
    loci <- unique(colnames(mat[,-1:-2]))
    freqTAB <- NULL
    mat2 <- mat
    mat2 <- sub("[.]","-",mat2)
    for(i in 1:length(loci)){
      ids <- which(colnames(mat)==loci[i])
      alleles <- unique(c(mat[,ids]))
      alleles <- sub("[.]","-",alleles)
      alleles <- alleles[!is.na(alleles)] ###
      nameCol <- paste(loci[i],".",alleles,sep = "")
      
      newmat <- matrix(NA,ncol = length(nameCol),nrow = dim(mat)[1])
      for(ii in 1:length(alleles)){
        newmat[,ii] <- apply(mat2[,ids]==alleles[ii],1,sum)
        colnames(newmat) <- nameCol
      }
      freqTAB <- cbind(freqTAB,newmat)
    }
    rownames(freqTAB) <- mat[,1]
    colnames(freqTAB) <- sub(" ","",colnames(freqTAB))
    dat2 <- as.genind(tab = freqTAB)
    pop(dat2) <- mat[,"pop"]
    
  } else {
    dat <- read.table(Ifile$datapath, header = TRUE,
                    sep = "\t",colClasses = "character")
    rownames(dat) <- dat$ind
    dat2 <- df2genind(dat[,-1:-2],ncode = switch(
      Incode, "2" = 2, "3" = 3),ploidy = switch(
        Iploidy, Diploid = 2, Haploid = 1));pop(dat2) <- dat$pop
  }
  return(dat2)
}

### convert input to genind object adapted to LD computation
genind4LD <- function(Ifile, Imicrovariants, Incode, Iploidy) {
  
  if(Imicrovariants == 2) {
    readLines(Ifile$datapath) -> matrix
    matrix <- strsplit(matrix, "[\t]")

  mat <- matrix(
    unlist(matrix),
    nrow = length(matrix),
    ncol = length(matrix[[1]]),
    byrow = TRUE
  )
  mat[mat=="0"] <- NA
  
  colnames(mat) <- mat[1, ]
  rownames(mat) <- mat[ ,1]
  
  mat <- mat[-1, ]
  mat <- sub("[.]", "", mat)
  mat[nchar(mat) == 2 & !is.na(mat)] <- paste(
    "0",mat[nchar(mat) == 2 & !is.na(mat)],
    
    sep=""
  )
  mat[nchar(mat) == 1 & !is.na(mat)] <- paste(
    "00",mat[nchar(mat) == 1 & !is.na(mat)],
    
    sep = ""
  )
  loci <- unique(colnames(mat[, -1:-2]))
  freqTAB <- NULL
  
  for(i in 1:length(loci)) {
    ids <- which(colnames(mat)==loci[i])
    alleles <- unique(c(mat[,ids]))
    alleles <- alleles[!is.na(alleles)] ###
    
    nameCol <- paste(loci[i],".",alleles,sep="")
    
    newmat <- matrix(NA,ncol=length(nameCol),nrow=dim(mat)[1])
    for(ii in 1:length(alleles)){
      newmat[,ii] <- apply(mat[,ids]==alleles[ii],1,sum)
      colnames(newmat) <- nameCol
    }
    freqTAB <- cbind(freqTAB,newmat)
  }
  rownames(freqTAB) <- mat[, 1]
  colnames(freqTAB) <- sub(" ", "", colnames(freqTAB))
  
  D <- genind(tab = freqTAB, pop = mat[, "pop"])
  
  } else {
    
    dat <- read.table(
      Ifile$datapath,
      header = TRUE,
      sep = "\t",
      colClasses = "character"
    )
    rownames(dat) <- dat$ind
    
    D <- df2genind(
      dat[, -1:-2],
      ncode = switch(
        Incode,
        "2" = 2,
        "3" = 3
      ),
      ploidy = switch(
        Iploidy,
        Diploid = 2,
        Haploid = 1
      )
    )
    pop(D) <- dat$pop
  }
  
  return(D)
}


## getFreqAllPop
## returns allele frequencies for each population
# input: genind object
# output: a list of allele frequencies tables
getFreqAllPop <- function(data) {
  
  freq <- list()
  freq$all <- getFreqFromGenind(data)
  
  for(popu in unique(data$pop)) {
    
    freq <- c(freq, x = NA)
    mat <- getFreqFromGenind(data[data@pop == popu, ])
    freq$x <- mat
    
    names(freq)[length(freq)] <- popu
    
  }
  
  return(freq)
}

## getFreqFromGenind
## returns allele frequencies for a given population
# input: genind object
# output: an allele frequencies tables
getFreqFromGenind <- function(data) {
  
  freq <- apply(data@tab, 2, sum, na.rm = TRUE)
  nam <- strsplit(
    names(freq),
    split = "[.]"
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
  
  DAT$alle[is.na(DAT$alle)] <- 0
  matr <- matrix(
    NA,
    ncol = length(unique(DAT$loc)),
    nrow = length(unique(DAT$alle))
  )
  rownames(matr) <- sort(unique(DAT$alle))
  colnames(matr) <- sort(unique(DAT$loc))
  
  for(i in sort(unique(DAT$alle))) {
    matr[
      as.character(i),
      DAT[DAT$alle == i, ]$loc
    ] <- DAT[DAT$alle==i, ]$frequency
  }
  
  return(matr)
}


## getIndicesAllPop
## returns indices for each population
# input: genind object
# output: a list of indices tables
getIndicesAllPop <- function(data,
                             hw = FALSE,
                             hwperm = 1000,
                             ploidy = "Diploid") {
  
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

## getIndicesFromGenind
## returns indices for a given population
# input: genind object
# output: an indices table
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

  PIC <- NULL
  for(i in unique(loc)) {
    
    FR <- c(DAT$frequency[names(DAT$frequency) == i])
    xu <- outer(FR, FR, "fu")
    som <- sum(xu[lower.tri(xu)])
    PIC[i] <-  1 - sum(FR ^ 2) - som
    
  } 
  
  Nall <- tapply(
    DAT[DAT$freq>0, ]$freq,
    DAT[DAT$freq>0, ]$loc,
    length
  )
  
  GD <- tapply(
    DAT$frequency,
    DAT$loc,
    function(x) 1 - sum(x ^ 2)
  )
  
  GD <- GD * N / (N - 1) 
  
  PIC <- PIC[names(GD)]
  
  D2 <- genind2loci(data)
  
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
  
  if(ploidy == "Diploid") {
    
    DF$Hobs <- summary(data)$Hobs[names(GD)]
    DF$PE <- (DF$Hobs ^ 2) * (1 - 2 * (DF$Hobs) * ((1 - DF$Hobs) ^ 2))
    DF$TPI <- 1 / (2 * (1 - DF$Hobs))
    
  }
  
  
  if(length(unique(data@pop)) > 1 & length(locNames(data)) > 1) {
    
    basicstat <- basic.stats(
      data,
      diploid = switch(
        ploidy,
        Diploid = TRUE,
        Haploid = FALSE
      ),
      digits = 4
    )
    
    Fst <- wc(
      data,
      diploid = switch(
        ploidy,
        Diploid = TRUE,
        Haploid = FALSE
      )
    )$per.loc$FST
    
    names(Fst) <- as.character(unique(data@loc.fac))
    DF$Fst <- Fst[names(GD)]
    
    DF$Ht <- basicstat$perloc[names(GD), "Ht"]
    DF$Fis <- basicstat$perloc[names(GD), "Fis"]

  }
  
  if(ploidy == "Diploid" & hw) {
    withProgress(message = 'Performing HW test...', value = 0, {
      DF$pHW <- hw.test(data, B = hwperm)[names(GD), 4]
    })
  } 
    
  return(DF)
}


fu <- function(a, b){
  2 * (a ^ 2) * (b ^ 2)
}


plotPCA <- function(pca, popus, coul, axis) {
  
  var1 <- round(100*(pca$eig/sum(pca$eig))[axis[1]], 2)
  var2 <- round(100*(pca$eig/sum(pca$eig))[axis[2]], 2)
  
  plot(pca$li[, axis[1]],
       pca$li[, axis[2]],
       col = transp(coul[popus], 0.5),
       pch = 16, cex = 1.7,
       xlab = paste0("PC", axis[1], " (", var1, " %)"),
       ylab = paste0("PC", axis[2], " (", var2, " %)"),
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5,
       bty = "l",
       main = "PCA projection"
  )
  
  sapply(unique(popus), function(x) {
    ellipse(
      c(mean(pca$li[popus %in% x, axis[1]]), mean(pca$li[popus %in% x, axis[2]])),
      cov(pca$li[, axis[1:2]]),
      0.95,
      col = transp(coul[x], 1),
      lty = 1,
      lwd = 2,
      center.pch = 0,
      fill = TRUE,
      fill.alpha = 0.2,
      segments = 100
    )
  })
  
  invisible(sapply(unique(popus), function(x) {
    legend(
      mean(pca$li[popus %in% x, axis[1]]),
      mean(pca$li[popus %in% x, axis[2]]),
      x,
      xjust = 0.5,
      yjust = 0.5,
      text.col=coul[x],
      text.font = 2,
      box.col = NA,
      bg = transp("white", 0.5),
      adj = 0.2
    )
  }))
  
}

