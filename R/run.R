#' @title Run mclexplorer
#' @description function to lunch the mclexplorer shiny application.
#'
#' @param 
#'
#' @author Pernice Simone
#' @import shiny bslib ggplot2 dplyr tidyr patchwork readxl sas7bdat survival survminer shinydashboard shinybusy shinyalert connector latex2exp labeling viridis ggmosaic purrr car cowplot
#' @rawNamespace import(DT, except=c(dataTableOutput,renderDataTable))
#' @rawNamespace import(shinyWidgets,except=alert)
#' 
#' @examples
#'\dontrun{
#' mclexplorer.run()
#' }
#' @export

mclexplorer.run <-function()
{
  x = T
  
  Appui <- system.file("Shiny","ui.R", package = "mclexplorer")
  Appserver <- system.file("Shiny","server.R", package = "mclexplorer")
  
  source(Appui)
  source(Appserver)
  
  app <-shinyApp(ui, server,
                 options =  options(shiny.maxRequestSize=1000*1024^2,
                                    shiny.launch.browser = .rs.invokeShinyWindowExternal)
  )
  
  app$staticPaths <- list(
    `/` = httpuv::staticPath(system.file("Shiny","www", package = "mclexplorer"), indexhtml = FALSE, fallthrough = TRUE)
  )
  
  shiny::runApp(app)
  
}
