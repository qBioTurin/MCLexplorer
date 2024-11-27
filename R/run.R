#' @title Run ShinyMCL
#' @description function to lunch the ShinyMCL shiny application.
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
#' ShinyMCL.run()
#' }
#' @export

ShinyMCL.run <-function()
{
  x = T
  
  Appui <- system.file("Shiny","ui.R", package = "ShinyMCL")
  Appserver <- system.file("Shiny","server.R", package = "ShinyMCL")
  
  source(Appui)
  source(Appserver)
  
  app <-shinyApp(ui, server,
                 options =  options(shiny.maxRequestSize=1000*1024^2,
                                    shiny.launch.browser = .rs.invokeShinyWindowExternal)
  )
  
  app$staticPaths <- list(
    `/` = httpuv::staticPath(system.file("Shiny","www", package = "ShinyMCL"), indexhtml = FALSE, fallthrough = TRUE)
  )
  
  shiny::runApp(app)
  
}
