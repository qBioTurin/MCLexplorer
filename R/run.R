#' @title Run MCLexplorer
#' @description Function to launch the MCLexplorer Shiny application.
#'
#' @author Pernice Simone
#' @import shiny bslib sas7bdat ggplot2 dplyr tidyr patchwork readxl sas7bdat survival survminer shinydashboard shinybusy shinyalert connector latex2exp labeling viridis ggmosaic purrr car cowplot
#' @rawNamespace import(DT, except=c(dataTableOutput,renderDataTable))
#' @rawNamespace import(shinyWidgets,except=alert)
#'
#' @examples
#' \dontrun{
#' mclexplorer.run()
#' }
#' @export

mclexplorer.run <- function() {
  Appui <- system.file("Shiny", "ui.R", package = "MCLexplorer")
  Appserver <- system.file("Shiny", "server.R", package = "MCLexplorer")
  if (Appui == "" || Appserver == "") {
    stop("Could not find the Shiny application files. Make sure 'MCLexplorer' is installed correctly.", call. = FALSE)
  }

  source(Appui)
  source(Appserver)

  app <- shinyApp(ui, server,
    options = list(shiny.maxRequestSize = 1000 * 1024^2)
  )

  app$staticPaths <- list(
    `/` = httpuv::staticPath(system.file("Shiny", "www", package = "MCLexplorer"), indexhtml = FALSE, fallthrough = TRUE)
  )

  shiny::runApp(app)
}
