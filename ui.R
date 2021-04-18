# UI ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(colourpicker)
  library(DT)
  library(plotly)
})
strider_pop <- c('Asia', 'AUSTRIA', 'BELGIUM', 'BOSNIA_AND_HERZEGOWINA', 'CZECH_REPUBLIC', 'DENMARK', 'Entire_Database', 'Europe', 'FINLAND', 'FRANCE', 'GERMANY', 'GREECE', 'HUNGARY', 'IRELAND', 'MONTENEGRO', 'NORWAY', 'POLAND', 'SAUDI_ARABIA', 'SLOVAKIA', 'SLOVENIA', 'SPAIN', 'SWEDEN', 'SWITZERLAND', 'THAILAND')

shinyUI(
  navbarPage(
    "STRAF 1.4.3: STR Analysis for Forensics",
    
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
        tags$head(tags$script(type="text/javascript", src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML")),

        sidebarLayout(
          sidebarPanel(
            p('STRAF performs forensics and population genetics analysis
              of STR data. Please read the documentation for details about input
              files and analyses.'),
            
            h4('Input'),
            radioButtons('microvariants', "Number of columns per locus:", c('2', '1'), inline = TRUE),
            radioButtons('ploidy', "Ploidy:", c('Diploid', 'Haploid'), inline = TRUE),
            conditionalPanel(
              condition="input.microvariants == 1",
              radioButtons('ncode', 'Number of digits for allele sizes:', c('2', '3'), inline = TRUE)
            ),
            fileInput(
              'file1', 'Choose file to upload:',
              accept = c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values', 'text/plain', '.csv', '.tsv')
            ),
            
            h4('Graphical parameters'),
            checkboxInput("hidegraph", "Display graphical parameters", FALSE),
            conditionalPanel(
              condition = "input.hidegraph",
              p("Barplot color"),
              colourInput("barplotcolor", NULL, "#36648B", showColour = "background"),
              checkboxInput("borderbarplot", "Bar border", FALSE),
              sliderInput("transparency", "Tranparency", 0, 1, 0.8, ticks = FALSE),
              sliderInput("width", "Plot width", 40, 100, 100, ticks = FALSE, post = "%"),
              sliderInput("height", "Plot height", 300, 800, 500, ticks = FALSE, post = "px"),
              sliderInput("cexaxis", "Axis label size", 0.2, 1.5, 1, ticks = FALSE),
              sliderInput("margin", "Margin", 1, 10, 7, ticks = FALSE)
            ),
            
            tags$hr(),
            h4('Contact'),
            p('Please address your questions and bug reports to Alexandre Gouy
              (alexandre.gouy [at] protonmail.com). Any suggestions are welcome!'),
            
            tags$hr(),
            h4('Citation'),
            p("Gouy, A., & Zieger, M. (2017). STRAF - A convenient online tool for STR data evaluation in forensic genetics. Forensic Science International: Genetics, 30, 148-151.")
          ),
          
          mainPanel(
            
            conditionalPanel(
              condition = "!output.fileUploaded",
              uiOutput("checkInputFile")
            ),
            
            conditionalPanel(
              condition = "output.fileUploaded",
              
              tabsetPanel(
                type = "tabs",
                
                tabPanel(
                  "Data",
                  
                  h3("Dataset"),
                  div(dataTableOutput('contents'), style = "font-size:70%"),
                  
                  tags$hr(),
                  h3("Allele frequencies per locus"),
                  checkboxInput(
                    'displayAlleleFreq', 'Plot the distribution of allele frequencies',
                    FALSE
                  ),
                  conditionalPanel(
                    condition = 'input.displayAlleleFreq == true',
                    uiOutput('plotAF')
                  ),
                  
                  tags$hr(),
                  checkboxInput(
                    'displayAlleleTable',
                    'Display a table of allele frequencies',
                    FALSE
                  ),
                  conditionalPanel(
                    condition = "input.displayAlleleTable == true",
                    uiOutput("selectPop"),
                    div(dataTableOutput('tableFreq'), style = "font-size:70%"),
                    downloadButton('dlTabfreq', 'Download as text (.tsv)'),
                    downloadButton('dlTabfreqXL', 'Download as Excel (.xlsx)')
                  ),
                  
                  tags$hr()
                ),
                
                tabPanel(
                  "Forensic parameters",
                  h3("Forensic parameters"),
                  checkboxInput(
                    'displayForensics',
                    'Compute forensics statistics (H, GD, PIC, PD, PE & TPI)',
                    FALSE
                  ),
                  
                  conditionalPanel(
                    condition = "input.displayForensics == true",
                    
                    uiOutput("selectPop2"),
                    div(tableOutput('forensics'), style = "font-size:75%"),
                    downloadButton('dlForensics', 'Download as text (.tsv)'),
                    downloadButton('dlForensicsXL', 'Download as Excel (.xlsx)'),
                    tags$hr(),
                    uiOutput("uiFOR"),
                    uiOutput("plotFOR")
                  ),
                  tags$hr()
                ),
                
                tabPanel(
                  "Population genetics indices",
                  h3("Summary statistics"),
                  checkboxInput(
                    'displayDiv',
                    'Compute heterozygosities and F-statistics',
                    FALSE
                  ),
                  
                  conditionalPanel(
                    condition = "input.displayDiv == true",
                    uiOutput("selectPop3")
                  ),
                  
                  conditionalPanel(
                    condition = "input.ploidy == 'Diploid'",
                    checkboxInput(
                      'computeHW', 'Test for Hardy-Weinberg equilibrium',
                      FALSE
                    ),
                    numericInput(
                      'hw_nperm', 'Number of permutations for HW test',
                      1000, min = 100, max = 10000, step = 100
                    )
                  ),
                  conditionalPanel(
                    condition = "input.displayDiv == true",
                    div(tableOutput('diversity'), style = "font-size:75%"),
                    downloadButton('dlPopgen', 'Download as text (.txt)'),
                    downloadButton('dlPopgenCL', 'Download as Excel (.xlsx)'),
                    tags$hr(),
                    uiOutput("uiPG"),
                    uiOutput("plotPG")
                  ),
                  
                  tags$hr(),
                  h3("Linkage disequilibrium"),
                  checkboxInput(
                    'displayLDtable', 'Display pairwise LD p-values matrix',
                    FALSE
                  ),
                  conditionalPanel(
                    condition = "input.displayLDtable == true",
                    div(tableOutput('LDtable'), style = "font-size:75%"),
                    downloadButton('dlLDtable', 'Download as text'),
                    downloadButton('dlLDtableXL', 'Download as Excel')
                  ),
                  conditionalPanel(
                    condition = "output.LD30",
                    checkboxInput(
                      'displayLDplot',
                      'Plot pairwise LD p-values matrix',
                      FALSE),
                    conditionalPanel(
                      condition = "input.displayLDplot == true",
                      uiOutput("plotLD2")
                    ),
                    conditionalPanel(
                      condition = "input.displayLDplot == true | input.displayLDtable == true",
                      checkboxInput(
                        'displayLDpvalplot',
                        'Plot LD p-values distribution',
                        FALSE
                      ),
                      conditionalPanel(
                        condition = "input.displayLDpvalplot == true",
                        uiOutput("plotLDpval2")
                      )
                    )
                  ),
                  
                  tags$hr(),
                  h3("Population structure"),
                  h4("Pairwise Fst"),
                  checkboxInput(
                    'displayFstMat',
                    'Compute pairwise Fst matrix',
                    FALSE
                  ),
                  conditionalPanel(
                    condition = "input.displayFstMat == true",
                    div(tableOutput('FstMat'), style = "font-size:75%"),
                    downloadButton('dlFstMat', 'Download as text (.txt)'),
                    downloadButton('dlFstMatXL', 'Download as Excel (.xlsx)')
                  ),
                  tags$hr()
                ),
                
                tabPanel(
                  "PCA - MDS",
                  h4("Principal Component Analysis (PCA)"),
                  checkboxInput(
                    'displayPCA',
                    'Run and plot a PCA (Principal Component Analysis)',
                    FALSE
                  ),
                  conditionalPanel(
                    condition = "input.displayPCA == true",
                    checkboxGroupInput('PCAaxis','PCA axis',c(1,2,3),c(1,2),inline = TRUE)
                  ),
                  conditionalPanel(
                    condition = "input.displayPCA == true",
                    uiOutput('plotPCA'),
                    verbatimTextOutput("info"),
                    downloadButton('dlPCAeigen', 'Download PCA eigenvectors'),
                    downloadButton('dlPCAcoord', 'Download PCA coordinates'),
                    checkboxInput('displayloadings', 'Plot loadings (alleles contributions)', FALSE)
                  ),
                  conditionalPanel(
                    condition = "input.displayPCA == true & input.displayloadings == true",
                    uiOutput('plotLoadings')
                  ),
                  
                  tags$hr(),
                  h4("Multidimensional Scaling (MDS) based on Nei's distance"),
                  checkboxInput(
                    'displayMDS',
                    "Compute Nei's genetic distance between populations and run MDS",
                    FALSE
                  ),
                  conditionalPanel(
                    condition = "input.displayMDS == true",
                    uiOutput('plotMDS')
                  ),
                  
                  tags$hr(),
                  h4("MDS on STRidER allele frequency database"),
                  checkboxInput("add_current", "Include uploaded data to the MDS", FALSE),
                  uiOutput('plotMDS_strider'),
                  checkboxGroupInput(
                    'location', 'Select populations',
                    choices = strider_pop, select = strider_pop,
                    inline = TRUE
                  ),
                  tags$hr(),
                  div("This MDS is performed on STRidER allele frequencies. Missing frequencies are imputed per allele as the mean allele frequency of other populations."), 
                  tags$a(href = "https://strider.online/frequencies", "Link to source data.")
                ),
                
                tabPanel(
                  "File conversion (beta)",
                  h4("NB: these new features have not been extensively tested yet."),
                  h3("Genepop (diploid data only)"),
                  downloadButton('dlGenepop', 'Download file in the Genepop format'),
                  h3("Arlequin (diploid data only)"),
                  downloadButton('dlArlequin', 'Download file in the Arlequin format'),
                  h3("Familias"),
                  downloadButton('dlFamilias', 'Download file in the Familias format')
                )
              )
            )
          )
        )
      )
    ),
    
    ##### DOCUMENTATION TAB ----------------------------------------------------
    tabPanel(
      "Documentation",
      fluidRow(
        column(width = 3),
        column(width = 6, includeMarkdown("./ui_files/doc.md")),
        column(width = 3)
      )
    ),
    
    ##### ABOUT STRAF TAB ------------------------------------------------------
    tabPanel(
      "About STRAF",
      fluidRow(
        column(width = 3),
        column(
          width = 6,
          p('STRAF is a browser-based application that allows to perform forensics 
  and population genetics analysis of STR data.'),
          includeMarkdown("./ui_files/changelog.md"),
          includeMarkdown("./ui_files/license.md"),
          includeMarkdown("./ui_files/acknowledgments.md"),
        ),
        column(width = 3)
      )
    )
  )
)
