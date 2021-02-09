# UI ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(colourpicker)
  library(DT)
  library(plotly)
  
})

shinyUI(

  navbarPage(

    "STRAF 1.4.0: STR Analysis for Forensics",
    
    ##### ANALYSIS TAB ----------------------------------------------------------
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
          
          sidebarPanel(
            
            p('STRAF performs forensics and population genetics analysis
              of STR data. Please read the documentation for details about input
              files and analysis.'),
            
            h4('Input'),
            
            radioButtons('microvariants', "Number of columns per locus:", c('2', '1'), inline = TRUE),
            radioButtons('ploidy', "Ploidy:", c('Diploid', 'Haploid'), inline = TRUE),

            conditionalPanel(
              
              condition="input.microvariants == 1",
              
              radioButtons(
                'ncode',
                'Number of digits for allele sizes:',
                c('2', '3'),
                inline = TRUE
              )
            ),

            fileInput(
              'file1', 'Choose file to upload:',
              accept = c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values', 'text/plain', '.csv', '.tsv')
            ),

            h4('Graphical parameters'),
            
            checkboxInput("hidegraph", "Display graphical parameters", FALSE),

            conditionalPanel(
              
              condition="input.hidegraph",
              
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
            p("Gouy, A., & Zieger, M. (2017). STRAF - A convenient online tool 
              for STR data evaluation in forensic genetics. 
              Forensic Science International: Genetics, 30, 148-151.")
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
                  
                  checkboxInput(
                    'displayTable',
                    'Display the dataset',
                    FALSE
                  ),
                  conditionalPanel(
                    condition = "input.displayTable == true",
                    div(dataTableOutput('contents'), style = "font-size:70%"),
                  ),

                  tags$hr(),
                  h3("Allele frequencies per locus"),
                  
                  checkboxInput(
                    'displayAlleleFreq',
                    'Plot the distribution of allele frequencies',
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
                  "Forensics analysis",

                  h3("Forensics statistics"),
                  
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
                  "Population genetics analysis",

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
                      'computeHW',
                      'Test for Hardy-Weinberg equilibrium',
                      FALSE
                    ),
                    numericInput(
                      'hw_nperm',
                      'Number of permutations for HW test',
                      1000,
                      min = 100,
                      max = 10000,
                      step = 100
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
                    'displayLDtable',
                    'Display pairwise LD p-values matrix',
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
                  )
                  
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
    
    ##### DOCUMENTATION TAB --------------------------------------------------------
    
    tabPanel(
      "Documentation",
      
      fluidRow(
        
        column(width=3),
        column(width=6,
               
               p('STRAF is a browser-based application that allows to perform forensics 
        and population genetics analysis of STR data.'),

               h3("Input file format"),

               p('STRAF accepts tab delimited txt-tables of different content. 
        The first column needs to contain the sample ID and the second the population 
        ID (if several populations are studied). Most convenient for the analysis of 
        forensically relevant autosomal STR data, STRAF accepts point alleles. Allele 
        data for haploid samples is entered with one column per locus and for diploid 
        data with two columns per locus. Missing data (e.g. null alleles) must be 
        indicated with a “0”. An example of a diploid input file is provided in the 
        supplementary data for this article.

        This format is designed to facilitate the input file generation from
        a typical Excel file (Save as > Text (Tab-delimited) (*.txt)).
        
        Examples of diploid and haploid input files can be downloaded using the 
        following links:'),

               tags$a('Download STRAF diploid example file.', 
                      target="_blank",
                      download="exampleSTRAFdiplo.txt",
                      href = 'exampleSTRAFdiplo.txt'),

               p(),
               tags$a('Download STRAF haploid example file.', 
                      target="_blank",
                      download="exampleSTRAFhaplo.txt",
                      href = 'exampleSTRAFhaplo.txt'),
               
               p(),
               p('It is also possible to use a format similar to the Genepop 
        input format, 
        with both alleles (for diploid data) coded in one column, either with 2 or 3 
        digits. Note that in order to use this format, no point alleles should be 
        present in the data set. Here is an example input file:'),

               tags$a('Download Genepop-like example file.', 
                      target="_blank",
                      download="exampleGenepop.txt",
                      href = 'exampleGenepop.txt'),
               
               h3("Using STRAF"),
               
               p("STRAF computes standard forensics parameters. Some standard population
        genetics analysis can be achieved if the samples are assigned to different 
        populations. STRAF generates downloadable tables, and plots can be personnalized 
        (using the Graphical parameters section on the left panel)
        before saving them. Details about the methods can be found in Gouy & Zieger (2017)."),

               p("Use the left panel to choose and upload your file. If no error appear
        once the file is uploaded, three tabs appear on the right page: 
        Data, Forensics analysis and Population genetics analysis."),

               p("On the Data tab, three checkboxes allow to 1. display the dataset, 
        2. plot the distribution of alleles frequencies per locus, and 3. display a 
        table of allele frequencies. This table is formatted as in most forensics data
        reports (rows = alleles; columns = loci). You can download this table as a 
        TSV file readable in Excel by clicking the Download button."),
               
               p("On the Forensics analysis tab, you can compute all the standard
                   forensics parameters for your dataset. You can download the table
                   and plot the results."),
               
               p("On the Population genetics analysis tab, standard population
                   genetics statistics can be computed (F-statistics, HWE, LD, ...).
                   A PCA can also be performed to study population structure or
                   discover outliers."),
               
               h3("Reference"),
               
               p("Please cite STRAF if you use it for your project
        using the following reference:"),
               p("Gouy, A., & Zieger, M. (2017). STRAF - A convenient online tool 
                  for STR data evaluation in forensic genetics. 
                   Forensic Science International: Genetics, 30, 148-151.")

        ),
        column(width=3)
      )
    ),
    
    ##### DOCUMENTATION TAB --------------------------------------------------------
    
    tabPanel(
      "About STRAF",
      
      fluidRow(
        
        column(width=3),
        column(width=6,
               
               p('STRAF is a browser-based application that allows to perform forensics 
                 and population genetics analysis of STR data.'),

               h3("Updates"),
               
               tags$ul(
                 tags$li("1.4.0 (09/02/2021) – STRAF can now perform an MDS. All results can be downloaded as Excel files. Minor bug fixes in file conversion."), 
                 tags$li("1.3.3 (03/02/2021) – Minor bug fixes in file conversion and PIC computation. Improved graphics."), 
                 tags$li("1.3.0 (01/02/2021) – STRAF can convert files to the Genepop and Familias formats. A File conversion tab has been added."), 
                 tags$li("1.2.2 (09/01/2021) – STRAF has moved to an AWS server (without any changes in the License)."), 
                 tags$li("1.1.2 (03/10/2020) – Preparation for integration to Tercen; few minor updates"), 
                 tags$li("1.0.5 (18/09/2018) – a few bug fixes; new page with
                         updates, license, data usage"), 
                 tags$li("1.0.4 (30/01/2018) - percentages of explained variance 
                         appear now on the PCA plot + style improvements"), 
                 tags$li("1.0.3 (18/12/2017) - a few bugs fixed, and it is now 
                         possible to analyse a single locus"),
                 tags$li("1.0.2 (20/10/2017) - forensics parameters and 
                         population genetics indices can be computed for 
                         each population separately; PCA coordinates and 
                         eigenvectors can be downloaded"), 
                 tags$li("1.0.1 (22/09/2017) - allele frequencies can now be 
                         computed and downloaded for each population separately"), 
                 tags$li("1.0.0 (02/06/2017) - STRAF release")
               ),
               
               h3("Data policy"),
               
               p("Any file uploaded to STRAF is deleted when the session is
                 closed. No data is retained on the server."),

               h3("License"),
               
               p("The STRAF software as a whole is distributed under GPL-3 (GNU 
                GENERAL PUBLIC LICENSE version 3). This program is free 
                software: you can redistribute it and/or modify it under the 
                terms of the GNU General Public License version 3 (GPL-3) as 
                published by the Free Software Foundation."),
               
               p("This program is distributed in the hope that it will be 
                useful, but WITHOUT ANY WARRANTY; without even the implied 
                warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
                PURPOSE.  See the GNU General Public License for more details 
                (https://www.gnu.org/licenses/gpl.txt)."),
               
               h3("Acknowledgments"),
               
               p("STRAF benefited of the comments of Martin Zieger, 
                 Peter Vallone, Peter de Knijff, Guanglin He and Martin Bodner.")
               
               ),
        column(width=3)
               )
        )
    
  )
)
