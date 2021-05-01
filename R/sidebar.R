#' Generate the sidebar UI.
#' @export
#' @noRd
#' @importFrom shinyWidgets awesomeCheckbox pickerInput
#' @importFrom colourpicker colourInput
sidebarUI <- function() {
  sidebarPanel(
    width = 3,
    fluidRow(
      column(width = 6,
             div(tags$img(src='STRAF_logo.png', height = "120"), style="text-align: center;"),
      ), 
      column( width = 6,
              p(strong('Welcome!'), br(), br(), 'STRAF is an STR data analysis application.'),
      )),
    h4('Input'),
    p('Please read the documentation for details about input files and analyses.'),
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
    # tags$hr(),
    # actionButton("load_example", "Download example dataset"),
    tags$hr(),
    
    h4('Graphical parameters'),
    awesomeCheckbox("hidegraph", "Display graphical parameters", FALSE),
    conditionalPanel(
      condition = "input.hidegraph",
      p("Barplot color"),
      colourInput("barplotcolor", NULL, "#36648B", showColour = "background"),
      awesomeCheckbox("borderbarplot", "Bar border", FALSE),
      sliderInput("transparency", "Tranparency", 0, 1, 0.8, ticks = FALSE),
      sliderInput("width", "Plot width", 40, 100, 100, ticks = FALSE, post = "%"),
      sliderInput("height", "Plot height", 300, 800, 500, ticks = FALSE, post = "px"),
      sliderInput("cexaxis", "Axis label size", 0.2, 1.5, 1, ticks = FALSE),
      sliderInput("margin", "Margin", 1, 10, 7, ticks = FALSE)
    ),
    
    tags$hr(),
    h4('The STRAF Book'),
    p('Click on the image below to open our online book with a lot more details about the software!'),
    HTML("<div align='center'><a href='https://agouy.github.io/straf/' target='_blank'><img src='cover.png' align='center' width='169' /></a></div>"),
    tags$hr(),
    
    h4('Contact'),
    p('Please address your questions and bug reports to Alexandre Gouy
              (alexandre.gouy [at] protonmail.com). Any suggestions are welcome!'),
    h4('Citation'),
    p("Gouy, A., & Zieger, M. (2017). STRAF - A convenient online tool for STR data evaluation in forensic genetics. Forensic Science International: Genetics, 30, 148-151.")
  )
}