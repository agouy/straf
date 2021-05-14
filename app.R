library(straf)
library(shiny)
library(markdown)  # required for shinyapps.io deployment
library(rmarkdown) #
dir <- system.file("application", package = "straf")
setwd(dir)
shiny::shinyAppDir(".")
