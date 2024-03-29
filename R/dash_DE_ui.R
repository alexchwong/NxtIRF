ui_DE <- function(id) {
    ns <- NS(id)
    DE_opts <- c()
    if(.check_package_installed("limma", "3.44.0", "silent")) {
        DE_opts <- "limma"
    }
    if(.check_package_installed("DESeq2", "1.30.0", "silent")) {
        DE_opts <- c(DE_opts, "DESeq2")
    }
    if(.check_package_installed("DoubleExpSeq", "1.1", "silent")) {
        DE_opts <- c(DE_opts, "DoubleExpSeq")
    }
    if(is.null(DE_opts)) {
        DE_opts <- ("(none) - please install limma, DESeq2 or DoubleExpSeq")
    }
    fluidRow(
        .ui_notice(),
        column(4,
            textOutput(ns("warning_DE")),
            selectInput(ns('method_DE'), 'Method', 
                DE_opts),
            selectInput(ns('variable_DE'), 'Variable', 
                c("(none)")),
            selectInput(ns('nom_DE'), 'Nominator', 
                c("(none)")),
            selectInput(ns('denom_DE'), 'Denominator', 
                c("(none)")),
            selectInput(ns('batch1_DE'), 'Batch Factor 1', 
                c("(none)")),
            selectInput(ns('batch2_DE'), 'Batch Factor 2', 
                c("(none)")),
            shinyWidgets::switchInput(ns("adjP_DE"), 
                label = "Sort by Adjusted P Values", 
                value = TRUE, labelWidth = "120px"),                    
            actionButton(ns("perform_DE"), "Perform DE"),
            shinyFilesButton(ns("load_DE"), label = "Load DE", 
                    title = "Load DE from RDS", multiple = FALSE),
            shinySaveButton(ns("save_DE"), "Save DE", "Save DE as...", 
                filetype = list(RDS = "Rds")),
        ),
        column(8,
            actionButton(ns("clear_selected_DE"), "Clear Selected Events"),
            div(style = 'overflow-x: scroll',  
                DT::dataTableOutput(ns('DT_DE'))
            )
        )
    )
}