#' Run the STRAF application
#' @description Main function to be called in order to start STRAF,
#' @param port TCP port.
#' @param host IP address.
#' @return Runs the shiny application.
#' @examples
#' library(straf)
#' # runStraf()
#' 
#' @export
#' 
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom graphics abline axis barplot hist image legend par
#' @importFrom stats as.dist cmdscale cov frequency ks.test qqplot qunif
#' @importFrom utils read.table write.table
runStraf <- function(port = 3838, host = "0.0.0.0") {
  options(
    warn = -1,
    shiny.sanitize.errors = FALSE,
    stringsAsFactors = FALSE
  )
  shiny::runApp(
    appDir = system.file("application", package = "straf"),
    port = port,
    host = host,
    display.mode = "normal"
  )
}
