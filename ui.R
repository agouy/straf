library(shiny)
library(DT)
library("shinythemes")
library("colourpicker")


shinyUI(
  navbarPage("STRAF: STR Analysis for Forensics",
    tabPanel("Analysis",
  
  fluidPage(theme = "bootstrap.css",
  tags$head(tags$style(type="text/css", ".container-fluid {max-width: 1200px}")),

  sidebarLayout(
    sidebarPanel(
      p('STRAF performs forensics and population genetics analysis
        of STR data. Please read the documentation for details about input files and analysis.
        A walkthrough example is also described.'),

      h4('Input'),
      # checkboxInput('microvariants', 'Presence of microvariants', TRUE),
      radioButtons('microvariants', "Number of columns per locus:" , c('2', '1'),inline = TRUE),
      radioButtons('ploidy', "Ploidy:", c('Diploid', 'Haploid'),inline = TRUE),
      
      conditionalPanel(
        condition="input.microvariants == 1",
        radioButtons('ncode', 'Number of digits for allele sizes:', c('2', '3'),inline = TRUE)
      ),
      fileInput('file1', 'Choose file to upload:',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )
      ),
      
      h4('Graphical parameters'),
      checkboxInput("hidegraph","Display graphical parameters",FALSE),
      conditionalPanel(
        condition="input.hidegraph",
        p("Barplot color"),
        colourInput(
          "barplotcolor", NULL, "#AEA79F",
          showColour = "background"),
        checkboxInput("borderbarplot","Bar border",TRUE),
        sliderInput("transparency","Tranparency",0,1,0.8,ticks=FALSE),
        sliderInput("width","Plot width",40,100,100,ticks=FALSE,post="%"),
        sliderInput("height","Plot height",300,800,500,ticks=FALSE,post="px")
        # sliderInput("cexaxis","Axis labels size",0.5,2,1,ticks=FALSE)
      ),

      tags$hr(),
      h4('Contact'),
      p('Please address your questions and bug reports to Alexandre Gouy
        (alexandre.gouy [at] iee.unibe.ch). Any suggestions are welcome!')
    ),
    
    mainPanel(
      conditionalPanel(
        condition = "output.fileUploaded",
        
      tabsetPanel(type = "tabs",
        tabPanel("Data",
                 
          h3("Dataset"),
          
          checkboxInput('displayTable', 'Display the dataset', FALSE),
          conditionalPanel(
            condition = "input.displayTable == true",
            div(dataTableOutput('contents'), style = "font-size:70%")
          ),
          
          tags$hr(),
          h3("Allele frequencies per locus"),
          
          checkboxInput('displayAlleleFreq', 'Plot the distribution of allele frequencies', FALSE),
          conditionalPanel(
            condition = "input.displayAlleleFreq == true",
            uiOutput("plotAF")
          ),
          
          tags$hr(),

          
          checkboxInput('displayAlleleTable', 'Display a table of allele frequencies', FALSE),
          conditionalPanel(
            condition = "input.displayAlleleTable == true",
            div(dataTableOutput('tableFreq'), style = "font-size:70%"),
            downloadButton('dlTabfreq', 'Download')
          ),
          tags$hr()
        ),
        
        tabPanel("Forensics analysis",
           h3("Forensics statistics"),
           checkboxInput('displayForensics', 'Compute forensics statistics (H, GD, PIC, PD, PE & TPI)', FALSE),


          conditionalPanel(
            condition = "input.displayForensics == true",
            
            div(tableOutput('forensics'), style = "font-size:75%"),
            downloadButton('dlForensics', 'Download'),
            tags$hr(),
            uiOutput("uiFOR"),
            uiOutput("plotFOR")
          ),
          tags$hr()

          
        ),
        
        tabPanel("Population genetics analysis",
                 
          h3("Summary statistics"),
          
          checkboxInput('displayDiv', 'Compute heterozygosities and F-statistics', FALSE),
          conditionalPanel(
            condition = "input.displayDiv == true & input.ploidy == 'Diploid'",
            checkboxInput('computeHW', 'Test for Hardy-Weinberg equilibrium', FALSE)
          ),
           
          conditionalPanel(
            condition = "input.displayDiv == true",
            
            div(tableOutput('diversity'), style = "font-size:75%"),
            downloadButton('dlPopgen', 'Download'),
            tags$hr(),
            uiOutput("uiPG"),
            uiOutput("plotPG")
          ),
          

          tags$hr(),
          h3("Linkage disequilibrium"),
          
            checkboxInput('displayLDtable', 'Display pairwise LD p-values matrix', FALSE),
          conditionalPanel(
            condition = "input.displayLDtable == true",
            div(tableOutput('LDtable'), style = "font-size:75%"),
            downloadButton('dlLDtable', 'Download')
            ),
          checkboxInput('displayLDplot', 'Plot pairwise LD p-values matrix', FALSE),
            conditionalPanel(
              condition = "input.displayLDplot == true",
              uiOutput("plotLD2")
          ),
          conditionalPanel(
            condition = "input.displayLDplot == true | input.displayLDtable == true",
            checkboxInput('displayLDpvalplot', 'Plot LD p-values distribution', FALSE),
            conditionalPanel(
              condition = "input.displayLDpvalplot == true",
              uiOutput("plotLDpval2")
            )
          ),
          
          tags$hr(),
          h3("Population structure"),
          
          h4("Pairwise Fst"),
          checkboxInput('displayFstMat', 'Compute pairwise Fst matrix', FALSE),
          conditionalPanel(
            condition = "input.displayFstMat == true",
            
            div(tableOutput('FstMat'), style = "font-size:75%"),
            downloadButton('dlFstMat', 'Download')
          ),
          
          tags$hr(),
          h4("Principal Component Analysis"),
           checkboxInput('displayPCA', 'Run and plot a PCA (Principal Component Analysis)', FALSE),
           conditionalPanel(
             condition = "input.displayPCA == true",
             checkboxGroupInput('PCAaxis','PCA axis',c(1,2,3),c(1,2),inline = TRUE)
           ),

          conditionalPanel(
            condition = "input.displayPCA == true",
            uiOutput('plotPCA'),
            checkboxInput('displayloadings', 'Plot loadings (alleles contributions)', FALSE)
            
          ),
          conditionalPanel(
            condition = "input.displayPCA == true & input.displayloadings == true",
            uiOutput('plotLoadings')
          ),
          tags$hr()
        )
      )
      )
    )
  )
)
),
  tabPanel("Documentation",
           fluidRow(
             column(width=3),
             column(width=6,
             p('STRAF is a browser-based application that allows to perform forensics and population genetics analysis
          of STR data.'),
             h3("Input file"),
             p('The input file must be a text file (.txt) with tab separated values. The first two colums
          include respectively the individual identifier and the population. The other columns contain the genotypes
          (one locus per column).'),
             p('Genotypes can be coded with 2 or 3 digits.'),
             p("It is usual in forensics to observe microvariants (i.e. partial repeats) for which the allele size
            is coded with a decimal (e.g. 12.1). In this case, Presence of microvariants must be checked
            and the input format changes (2 columns per loci)."),
             
             h3("Walkthrough example"),
             
             p("Example input file + tutorial coming soon."),
             p("+ Example file microvariants."),
             tags$a('Download example file.',target="_blank", download="example.txt", href = 'exampleFile.txt'),
             p(),
             tags$img(src = "tuto1.png", width = "200px"),
             p(),
             p("Once the file is uploaded, three tabs appear on the Analysis page: Data, Forensics analysis
               and Population genetics analysis."),
             
             p("1. On the Data tab, three checkboxes allow to 1. display the dataset, 2. plot the distribution
               of alleles frequencies per locus, and 3. display a table of allele frequencies. This table is
               formatted as in most forensics data reports (rows = alleles; columns = loci). 
               You can download this table as a TSV file readable in Excel by clicking the Download button."),
             
             p("2. Forensics analysis"),
             
             p("3. Population genetics analysis"),
             
             h3("Methods"),
             p("STRAF computes standard forensics parameters. Some standard population genetics analysis
             can be achieved if the samples are assigned to different populations."),
             
             h4("Forensics parameters"),
             withMathJax(),
             p('For each locus, with n the number of alleles at a given locus and pi the frequency of allele i, 
              the following statistics are computed:'),

            p('Genetic diversity, or expected heterozygosity:
            
                          $$ H_{\\textrm{exp}} = GD = \\sum\\limits_{i=1}^{n} (p_{i})^{2}$$
            
            Polymorphism information content:
            
                           $$ PIC = 1 - \\sum\\limits_{i=1}^{n} (p_{i})^{2} - (\\sum\\limits_{i=1}^{n} (p_{i})^{2})^2 + \\sum\\limits_{i=1}^{n} (p_{i})^{4}$$
            
            Power of discrimination:
            
            $$ PD = 1 - 2(\\sum\\limits_{i=1}^{n} (p_{i})^{2})^2 - \\sum\\limits_{i=1}^{n} (p_{i})^{4}$$
            
            Power of exclusion:
            
            $$ PE = H_{\\textrm{exp}}^2(1-(1-H_{\\textrm{exp}})H_{\\textrm{exp}}^2)$$
            
            Paternity index:
            
            $$ PI = \\frac{1}{2\\sum\\limits_{i=1}^{n}(p_{i})^{2}}$$
            
            '),
            
            h4("Population genetics analysis"),
            
            tags$b("Per locus heterozygozities and F-statistics"),
            p("H and F-stats
              $$ H_{\\textrm{exp}} = \\sum\\limits_{i=1}^{n} (p_{i})^{2}$$"),
            
            tags$b("Hardy-Weinberg equilibrium testing"),
            p("Hardy-Weinberg equilibrium is tested by computing the p-value of an exact test
            based on 1,000 Monte Carlo permutations of alleles."),
            
            tags$b("Linkage disequilibrium"),
            p("Linkage disequilibrium is computed
              using the T2 statistic (Zaykin et al., 2008). Statistical testing is achieved
              using the chi-square approximation of this statistic. No permutation test is yet
              implemented."),
            tags$b("Pairwise Fst"),
            p("Pairwise Fst matrix"),
            tags$b("PCA"),
            p("PCA")
             ),
             column(width=3))
  )
)
)