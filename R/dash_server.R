dash_server = function(input, output, session) {
    # Volumes / storage drives
    default_volumes <- c("Working Directory" = getwd(), 
        "Home" = "~", getVolumes()())
    volumes = reactive(default_volumes)
    # Defines for Reactives
    settings_expr <- settings_expr_load <- settings_filtered_SE <- 
        settings_DE <- settings_Diag <- settings_Volc <- c()
    # Defines for Page Refreshers
    settings_refresh <- reactiveValues(
        new_ref = c(), 
        expr = c(), 
        expr_load = c(), qc = c(), filters = c(), DE = c(),
        diag = c(), volc = c(), heat = c(), cov = c()
    )
    
    refresh_newref      <- reactive(settings_refresh$new_ref)
    
    refresh_expr        <- reactive(settings_refresh$expr)
    
    refresh_exprload    <- reactive(settings_refresh$expr_load)
    refresh_QC          <- reactive(settings_refresh$qc)
    refresh_filters     <- reactive(settings_refresh$filters)
    refresh_DE          <- reactive(settings_refresh$DE)
    
    refresh_diag        <- reactive(settings_refresh$diag)
    refresh_volc        <- reactive(settings_refresh$volc)
    refresh_heat        <- reactive(settings_refresh$heat)
    refresh_cov         <- reactive(settings_refresh$cov)
    
        # Reactives for shared objects
    
    # NxtSE object
    get_se_reactive <- reactive(settings_expr_load$se)
    # Current differential ASE results
    get_de_reactive <- reactive(settings_DE$res)  
    # Filtered NxtSE object
    get_filtered_se_reactive <- reactive({
        req(settings_expr_load$se)
        if(is_valid(settings_filtered_SE$filterSummary)) {
            settings_expr_load$se[settings_filtered_SE$filterSummary,]
        } else {
            settings_expr_load$se
        }
    })  
    # Collate path
    get_se_path_reactive <- reactive(settings_expr_load$collate_path) 
    # Annotation object
    get_df_anno_reactive <- reactive(settings_expr_load$df.anno) 
    # Filters from Filter tab
    get_filters_reactive <- reactive(settings_filtered_SE$filters)
    # Filters from DE tab (loading saved DE)    
    get_filters_DE_reactive <- reactive(settings_DE$filters) 
    # Filtered DE rows
    get_rows_all <- reactive(settings_DE$DT_DE_rows_all)    
    
    # Two-way talk between selected rows 
    get_rows_selected <- reactive(settings_DE$DT_DE_rows_selected)
    get_rows_selected_diag <- reactive(settings_Diag$selected)
    get_rows_selected_volc <- reactive(settings_Volc$selected)
        
    # Reactive that returns the number of threads to use
    get_threads_reactive <- reactive(.dash_get_threads(
        input$thread_option, input$cores_numeric))
    
    # Tie module data to their server objects
    settings_system <- setreactive_system()
    settings_newref <- server_ref_new("new_ref", refresh_newref, volumes)
    settings_expr <- server_expr("build_expr", refresh_expr, volumes, 
        get_threads_reactive)
    settings_expr_load <- server_expr("load_expr", refresh_exprload, volumes, 
        get_threads_reactive, limited = TRUE)
    settings_QC <- server_qc("qc", refresh_QC, get_se_path_reactive, 
        get_df_anno_reactive)
    settings_filtered_SE <- server_filters("filters", refresh_filters, volumes, 
        get_se_reactive, get_filters_DE_reactive)
    settings_DE <- server_DE("DE", refresh_DE, volumes, get_threads_reactive,
        get_filtered_se_reactive, get_filters_reactive,
        get_rows_selected_diag, get_rows_selected_volc)
    settings_Diag = server_vis_diag("diag", refresh_diag, volumes, 
        get_filtered_se_reactive, get_de_reactive,
        get_rows_all, get_rows_selected)
    settings_Volc = server_vis_volcano("volcano", refresh_diag, volumes, 
        get_filtered_se_reactive, get_de_reactive,
        get_rows_all, get_rows_selected)
    settings_Heat = server_vis_heatmap("heatmap", refresh_diag, volumes, 
        get_filtered_se_reactive, get_de_reactive,
        get_rows_all, get_rows_selected)
    settings_Cov <- server_cov("cov", refresh_cov, volumes, 
        get_filtered_se_reactive, get_de_reactive,
        get_rows_all, get_rows_selected)

# tabEvent Observer
    observeEvent(input$navSelection, {
        if(input$navSelection == "navTitle") {

        } else if(input$navSelection == "navRef_New") {
            settings_refresh$new_ref = runif(1)
        } else if(input$navSelection == "navSystem") {

        } else if(input$navSelection == "navExpr") {
            settings_refresh$expr = runif(1)
        } else if(input$navSelection == "navExprLoad") {
            settings_refresh$expr_load = runif(1)
        } else if(input$navSelection == "navQC") {
            settings_refresh$qc = runif(1)
        } else if(input$navSelection == "navFilter") {
            settings_refresh$filters = runif(1)
        } else if(input$navSelection == "navAnalyse") {
            settings_refresh$DE = runif(1)
        } else if(input$navSelection == "navDiag") {
            settings_refresh$diag = runif(1)
        } else if(input$navSelection == "navDiag") {
            settings_refresh$volc = runif(1)
        } else if(input$navSelection == "navHeatmap") {
            settings_refresh$heat = runif(1)
        } else if(input$navSelection == "navCoverage") {
            settings_refresh$cov = runif(1)
        }
    })
# End of server function
}

.dash_get_threads <- function(thread_option, cores_numeric) {
    if(thread_option == "Single-Thread"){
        n_threads = 1
    } else if(thread_option == "Multi-Thread (Low)") {
        n_threads = min(parallel::detectCores(), 4)
        if(n_threads < 1) n_threads = 1
    } else if(thread_option == "Multi-Thread (High)") {
        n_threads = min(parallel::detectCores(), 16)
        if(n_threads < 1) n_threads = 1
    } else if(thread_option == "Custom") {
        n_threads = cores_numeric
    }
    n_threads
}