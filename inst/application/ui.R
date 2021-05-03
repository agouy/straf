library(straf)
shinyUI(
  navbarPage(
    paste0("STRAF ", packageVersion("straf"), ": STR Analysis for Forensics"),
    
    ##### ANALYSIS TAB ---------------------------------------------------------
    tabPanel(
      "Analysis",
      fluidPage(
        theme = "bootstrap.css",
        tags$head(includeHTML(("./www/googleanalytics.html"))),
        tags$head(tags$style(type="text/css", ".container-fluid {max-width: 1200px}")),
        tags$head(tags$style('body {font-family: Arial;}')),
        tags$head(tags$style('table {font-family: Arial;}')),
        tags$head(tags$style('h1 {font-family: Arial;}')),
        tags$head(tags$style('h2 {font-family: Arial;}')),
        tags$head(tags$style('h3 {font-family: Arial;}')),
        tags$head(tags$style('h4 {font-family: Arial;}')),
        
        sidebarLayout(
          sidebarUI(),
          mainPanel(
            width = 9,
            conditionalPanel(
              condition = "!output.fileUploaded",
              uiOutput("checkInputFile")
            ),
            
            conditionalPanel(
              condition = "output.fileUploaded",
              
              tabsetPanel(
                type = "tabs",
                data_UI("data_ns"),
                for_UI("for_popgen"),
                popgen_UI("for_popgen"),
                pca_mds_UI("pca_mds"),
                ref_mds_UI("ref_mds"),
                file_conv_UI("file_conv")
              )
            )
          )
        )
      )
    ),
    documentation_tab(),
    about_tab()
  )
)

