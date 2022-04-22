#' Generate "Documentation" panel.
#' @export
#' @keywords internal
documentation_tab <- function() {
  tabPanel(
    "Documentation",
    fluidRow(
      column(width = 3),
      column(width = 6, includeMarkdown("./ui_files/doc.md")),
      column(width = 3)
    )
  )
}

#' Generate "About" panel.
#' @export
#' @keywords internal
about_tab <- function() {
  tabPanel(
    "About STRAF",
    fluidRow(
      column(width = 3),
      column(
        width = 6,
        p('STRAF is a browser-based application that allows to perform forensics 
    and population genetics analysis of STR data.'),
        includeMarkdown("./ui_files/changelog.md"),
        includeMarkdown("./ui_files/license_ui.md"),
        includeMarkdown("./ui_files/acknowledgments.md"),
      ),
      column(width = 3)
    )
  )
}
