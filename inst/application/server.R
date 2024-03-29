shinyServer(function(input, output) {
  
  getData <- reactive({
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
      
      return(input$file1)
    }
  })
  
  output$fileUploaded <- reactive({ return(!is.null(getData())) })
  
  output$checkInputFile <- renderText({
    
    if(is.null(input$file1)) { return("Please import a data set.") }
    else {
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
        return("Input file error. Wrong number of columns per locus (please check ploidy).")
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
  
  getgenind <- reactive({
    inFile <- input$file1
    if(is.null(inFile)) return(NULL)
    df_out <- createGenind(
      Ifile = inFile,
      ploidy = input$ploidy
    )
    return(df_out)
  })
  
  output$dl_example <- downloadHandler(
    filename = function() { paste('allele_frequencies.tsv', sep='') },
    content = function(file) {
      if (!input$displayAlleleTable) return(NULL)
      if(input$selectPop == "") matr <- alleleFreqTabs()[[1]]
      else matr <- alleleFreqTabs()[[input$selectPop]]
      write.table(matr, file, sep = "\t", na = "", row.names = TRUE)
    }
  )
  
  ploidy <- reactive({input$ploidy})
  hw_perm <- reactive({input$hw_nperm})
  barplotcolor <- reactive({input$barplotcolor})
  transparency <- reactive({input$transparency})
  cexaxis <- reactive({input$cexaxis})
  width <- reactive({input$width})
  height <- reactive({input$height})
  input_file <- reactive({input$file1$datapath})
  
  data_Server("data_ns", getgenind, getData, barplotcolor, transparency, width, height, popnames)
  for_popgen_Server(
    "for_popgen", input_file, getgenind,
    popnames, ploidy, hw_perm, barplotcolor, transparency, cexaxis
  )
  pca_mds_Server("pca_mds", getgenind)
  ref_mds_Server("ref_mds", getgenind)
  file_conv_Server("file_conv", reactive({input$file1$datapath}), ploidy, getgenind, popnames)
  
  autoInvalidate <- reactiveTimer(30000)
  timeout <- reactiveValues(counter = 0)
  observe({
    autoInvalidate()
    timeout$counter <- isolate(timeout$counter) + 1
    if(timeout$counter < 10) {
      cat(".")
    } 
  })
  
})