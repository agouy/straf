#' straf: STR Analysis for Forensics
#'
#' straf is a Shiny application to perform Short Tandem Repeats (STRs, also 
#' known as microsatellites) data analysis. The application allows one to 
#' compute forensic parameters, population genetics indices, and investigate 
#' population structure through various methods and generate relevant data 
#' visualisations. It also implements file conversion to other popular formats.
#' 
#' @section Running the app:
#' One simply needs to call the runStraf() function to start the application.
#'
#' @docType package
#' @name straf
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom adegenet as.genind df2genind pop<- locNames genind2genpop dist.genpop makefreq pop transp
#' @importClassesFrom adegenet genind
#' @importFrom ape as.phylo
#' @importFrom colourpicker colourInput
#' @importFrom DT dataTableOutput renderDataTable datatable
#' @importFrom genepop clean_workdir test_HW test_LD
#' @importFrom ggplot2 ggplot geom_point labs theme_minimal aes geom_tile
#' @importFrom ggrepel geom_text_repel
#' @importFrom graphics abline axis barplot hist image legend par
#' @importFrom hierfstat pairwise.WCfst genind2hierfstat wc
#' @importFrom magrittr "%>%"
#' @importFrom openxlsx write.xlsx
#' @importFrom pegas genind2loci
#' @importFrom plotly renderPlotly plotlyOutput ggplotly
#' @importFrom readxl read_excel
#' @importFrom reshape2 acast
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyWidgets awesomeCheckbox pickerInput
#' @importFrom stats as.dist cmdscale cov frequency hclust ks.test median qqplot qunif p.adjust prcomp
#' @importFrom tidyr gather
#' @importFrom utils read.table write.table count.fields
NULL
#> NULL