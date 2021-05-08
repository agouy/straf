library(straf)
dir <- system.file("application", package = "straf")
setwd(dir)
shiny::shinyAppDir(".")
