#' Run the STRAF application
#' @description Main function to be called in order to start STRAF,
#' @return Runs the shiny application.
#' @examples
#' # runStraf()
#' 
#' @export
#' 
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @import dplyr
#' @import ggplot2
#' @importFrom graphics abline axis barplot hist image legend par
#' @importFrom stats as.dist cmdscale cov frequency ks.test qqplot qunif
#' @importFrom utils read.table write.table
runStraf <- function() {
  options(
    warn = -1,
    shiny.sanitize.errors = FALSE,
    stringsAsFactors = FALSE
  )
  shiny::runApp(appDir = system.file("application", package = "straf"), port = 80)
}
