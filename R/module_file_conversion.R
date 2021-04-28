## File conversion module
file_conv_UI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "File conversion",
    h3("Genepop"),
    downloadButton(ns('dlGenepop'), 'Download file in the Genepop format'),
    tags$br(),
    h3("Familias"),
    downloadButton(ns('dlFamilias'), 'Download file in the Familias format'),
    tags$br(),
    h3("Arlequin (diploid data only)"),
    downloadButton(ns('dlArlequin'), 'Download file in the Arlequin format'),
    tags$br()
  )
}

file_conv_Server <- function(id, fpath, ploidy) {
  moduleServer(
    id,
    function(input, output, session) {
      
      #### FILE CONVERSION
      output$dlGenepop <- downloadHandler(
        filename = function() { 
          paste('straf2genepop.txt', sep='') 
        },
        content = function(file) {
          gp <- straf2genepop(f.name = fpath(), ploidy = switch(ploidy(), Diploid = 2, Haploid = 1))
          cat(gp, file = file)
        }
      )
      
      output$dlArlequin <- downloadHandler(
        filename = function() { 
          paste('straf2arlequin.arp', sep='') 
        },
        content = function(file) {
          gp <- straf2arlequin(fpath())
          cat(gp, file = file)
        }
      )
      
      output$dlFamilias <- downloadHandler(
        filename = function() { 
          paste('straf2familias.txt', sep='')
        },
        content = function(file) {
          fmi <- straf2familias(fpath())
          cat(fmi, file = file)
        }
      )
    }
  )}