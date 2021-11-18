#' Launches the NxtIRF Graphics User Interface using Shiny Dashboard
#' 
#' This function launches the NxtIRF interactive app using Shiny Dashboard
#' @param mode (default `"dialog"`) `"dialog"` displays NxtIRF in a dialog box 
#'   with specified width and height. `"browser"` opens NxtIRF in a browser-
#'   like resizable window.
#' @param width,height If `mode` is set to `"dialog"`, the specified width
#'   and height of the NxtIRF app.
#' @return An interactive shinydashboard NxtIRF app runs.
#' @examples
#' \dontrun{
#' # Launches interactive ShinyDashboard NxtIRF app as fixed-size dialog box
#' nxtIRF(mode = "dialog", width = 1600, height = 900) 
#'
#' # Launches interactive ShinyDashboard NxtIRF app as browser window
#' nxtIRF(mode = "browser") 
#' }
#' @md
#' @export
nxtIRF <- function(mode = c("dialog", "browser"), 
        width = 1600, height = 900) {
    if(!interactive()) {
        .log(paste("In nxtIRF(),",
            "NxtIRF App can only be run in interactive mode (i.e. RStudio)."))
    }
    mode = match.arg(mode)
    ui_dash <- dashboardPage(
        dashboardHeader(title = "NxtIRF"),
        ui_sidebar(),
        dashboardBody(
            tabItems(
                ui_tab_title(),     # "navTitle"
                
                ui_tab_system(),    # "navSystem"
                
                ui_tab_ref_new(),   # "navRef_New"    

                ui_tab_expr(),      # "navExpr"
                
                ui_tab_expr_load(), # "navExprLoad"
                ui_tab_qc(),        # "navQC"
                ui_tab_filter(),    # "navFilter"
                ui_tab_analyse(),   # "navAnalyse"

                ui_tab_diag(),      # "navDiag"
                ui_tab_volcano(),   # "navVolcano"
                ui_tab_heatmap(),   # "navHeatmap"
                ui_tab_coverage()   # "navCoverage"
            )
        )
    )
    if(mode == "dialog") {
        runGadget(
            shinyApp(ui_dash, dash_server),
            viewer = dialogViewer('NxtIRF', width = width, height = height)
        )
    } else {
        runApp(shinyApp(ui_dash, dash_server))
    }
}