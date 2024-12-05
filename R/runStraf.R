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
    quiet = TRUE,
    display.mode = "normal"
  )
}
