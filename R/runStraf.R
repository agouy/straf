#' @import shiny
#' @import colourpicker
#' @import plotly
#' @import shinyWidgets
#' @import ade4
#' @import adegenet
#' @import pegas
#' @import hierfstat
#' @import car
#' @import openxlsx
#' @import reshape2
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import ggrepel

runStraf <- function() {
  options(
    warn = -1,
    shiny.sanitize.errors = FALSE,
    stringsAsFactors = FALSE
  )
  shiny::runApp(appDir = system.file("application", package = "straf"), port = 80)
}
