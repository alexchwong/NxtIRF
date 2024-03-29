server_vis_diag <- function(
        id, refresh_tab, volumes, get_se, get_de,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Diag <- setreactive_Diag()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
            output$warning_diag <- renderText({
                validate(need(get_se(), 
                "Please load Differential Expression via 'Analysis' tab"))
                
                "Differential Expression Loaded"
            })
            req(get_se())
            colData <- colData(get_se())
            if(
                    is_valid(input$variable_diag) && 
                    input$variable_diag %in% colnames(colData)
            ) {
                selected <- isolate(input$variable_diag)
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = selected
                )
            } else {
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = "(none)"
                )
            }

        })
        observeEvent(rows_selected(), {
            settings_Diag$selected = rows_selected()
        })
    
        output$plot_diag <- renderPlotly({
            # settings_Diag$plot_ini = FALSE
            validate(need(get_se(), "Load Experiment first"))
            validate(need(get_de(), "Load DE Analysis first"))
            validate(need(input$variable_diag, 
                "Select conditions and contrasts"))
            validate(need(input$nom_diag, 
                "Select conditions and contrasts"))
            validate(need(input$denom_diag, 
                "Select conditions and contrasts"))
            validate(need(input$variable_diag != "(none)", 
                "Select conditions and contrasts"))
            validate(need(input$nom_diag != "(none)", 
                "Select conditions and contrasts"))
            validate(need(input$denom_diag != "(none)", 
                "Select conditions and contrasts"))

            selected <- settings_Diag$selected      

            num_events <- input$number_events_diag
            res <- as.data.table(get_de()[rows_all(),])
            if(is_valid(input$EventType_diag)) {
                res <- res[get("EventType") %in% input$EventType_diag]
            }
            if(num_events < nrow(res)) {
                res <- res[seq_len(num_events)]
            }
            df.diag <- make_diagonal(
                get_se(), res$EventName, 
                input$variable_diag, input$nom_diag, input$denom_diag
            )

            if(is_valid(settings_Diag$selected)) {
                df.diag$selected <- 
                    (df.diag$EventName %in% get_de()$EventName[selected])
            } else {
                df.diag$selected <- FALSE
            }
            df.diag$NMD_direction <- get_de()$NMD_direction[
                match(df.diag$EventName, get_de()$EventName)]
            
            settings_Diag$plot_ini = TRUE
            if(input$NMD_diag == TRUE) {
                df.diag             <- df.diag[df.diag$NMD_direction != 0, ]
                df.diag$nom_NMD     <- ifelse(df.diag$NMD_direction == 1, 
                                        df.diag$nom, df.diag$denom)
                df.diag$denom_NMD   <- ifelse(df.diag$NMD_direction == -1, 
                                        df.diag$nom, df.diag$denom)
                                        
                p <- ggplot(df.diag, 
                        aes(
                            x = get("nom_NMD"), y = get("denom_NMD"), 
                            key = get("EventName"), text = get("EventName"), 
                            colour = get("selected")
                        )
                    ) + geom_point() + 
                    scale_color_manual(values = c("black", "red")) +
                    labs(
                        x = paste(input$nom_diag, "NMD substrate"),
                        y = paste(input$denom_diag, "NMD substrate")
                    )
            } else {
                p <- ggplot(df.diag, 
                        aes(
                            x = get("nom"), y = get("denom"), 
                            key = get("EventName"), text = get("EventName"), 
                            colour = get("selected")
                        )
                    ) + geom_point() + 
                    scale_color_manual(values = c("black", "red")) +
                    labs(
                        x = paste(input$nom_diag),
                        y = paste(input$denom_diag)
                    )         
            }
            p <- p + labs(color = "Selected")
            settings_Diag$final_plot <- ggplotly(
                p, tooltip = "text",
                source = "plotly_diagonal") %>% 
                layout(
                    dragmode = "lasso",
                    yaxis = list(scaleanchor="x", scaleratio=1)
                )      
            print(
                settings_Diag$final_plot
            )
        })
    
        observe({
            shinyFileSave(
                input, "saveplot_diag", 
                roots = volumes(), session = session,
                filetypes = c("pdf")
            )
        })
        observeEvent(input$saveplot_diag, {
            req(settings_Diag$final_plot)
            selectedfile <- parseSavePath(volumes(), input$saveplot_diag)
            req(selectedfile$datapath)

            obj <- isolate(settings_Diag$final_plot)
            plotly::orca(obj, .make_path_relative(getwd(), 
                selectedfile$datapath))
        })

        settings_Diag$plotly_click = reactive({
            plot_exist <- settings_Diag$plot_ini
            if(plot_exist) 
                event_data("plotly_click", source = "plotly_diagonal")
        })
    
        observeEvent(settings_Diag$plotly_click(), {
            req(settings_Diag$plotly_click())
            click <- settings_Diag$plotly_click()
            # print(click)
            click.id <- which(get_de()$EventName == click$key)
            req(click.id)

            selected <- settings_Diag$selected

            if(click.id %in% selected) {
                selected <- selected[-which(selected == click.id)]
            } else {
                selected <- c(selected, click.id)
            }
            settings_Diag$selected <- selected
            # DT::dataTableProxy("DT_DE") %>% DT::selectRows(selected)
        })

        settings_Diag$plotly_brush = reactive({
            plot_exist <- settings_Diag$plot_ini
            if(plot_exist)
                event_data("plotly_selected", source = "plotly_diagonal")
        })
    
        observeEvent(settings_Diag$plotly_brush(), {
            req(settings_Diag$plotly_brush())
            brush <- settings_Diag$plotly_brush()
            # print(brush)
            brush.id <- which(get_de()$EventName %in% brush$key)
            req(brush.id)

            selected <- settings_Diag$selected
            selected <- unique(c(selected, brush.id))
            settings_Diag$selected <- selected
            # DT::dataTableProxy("DT_DE") %>% DT::selectRows(selected)
        })
    
        observeEvent(input$variable_diag, {
            req(get_se())
            req(input$variable_diag != "(none)")
            colData <- colData(get_se())
            req(input$variable_diag %in% colnames(colData))

            if(!is(colData[,input$variable_diag], "factor")) {
                output$warning_diag <- renderText(
                    "Contrast must be performed on discrete categories")
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = "(none)"
                )
            } else {
                output$warning_diag <- renderText("")
                col_levels <- levels(colData[,input$variable_diag])
                updateSelectInput(
                    session = session, inputId = "nom_diag", 
                    choices = c("(none)", col_levels), 
                    selected = "(none)"
                )
                updateSelectInput(
                    session = session, inputId = "denom_diag", 
                    choices = c("(none)", col_levels), 
                    selected = "(none)"
                )
            }
        })

        observeEvent(input$clear_diag, {
            updateSelectInput(session = session, 
                "EventType_diag", selected = NULL)
            shinyWidgets::updateSliderTextInput(session = session, 
                "number_events_diag", selected = 10000)
            
            if(is_valid(get_se())) {
                colData <- colData(get_se())
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = "(none)"
                )
            } else {
                updateSelectInput(session = session, inputId = "variable_diag", 
                    choices = c("(none)"), selected = "(none)")
            }
            
            updateSelectInput(session = session, inputId = "nom_diag", 
                choices = c("(none)"), selected = "(none)")
            updateSelectInput(session = session, inputId = "denom_diag", 
                choices = c("(none)"), selected = "(none)")
        })
    
        return(settings_Diag)        
    })
}

server_vis_volcano <- function(
        id, refresh_tab, volumes, get_se, get_de,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Volc <- setreactive_Diag()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
        })
        observeEvent(rows_selected(), {
            settings_Volc$selected <- rows_selected()
        })

        settings_Volc$plotly_click <- reactive({
            plot_exist <- settings_Volc$plot_ini
            if(plot_exist) 
                event_data("plotly_click", source = "plotly_volcano")
        })
    
        observeEvent(settings_Volc$plotly_click(), {
            req(settings_Volc$plotly_click())
            click <- settings_Volc$plotly_click()
            click.id <- which(get_de()$EventName == click$key)
            req(click.id)

            selected <- settings_Volc$selected

            if(click.id %in% selected) {
                selected <- selected[-which(selected == click.id)]
            } else {
                selected <- c(selected, click.id)
            }
            settings_Volc$selected <- selected
        })

        settings_Volc$plotly_brush <- reactive({
            plot_exist <- settings_Volc$plot_ini
            if(plot_exist)
                event_data("plotly_selected", source = "plotly_volcano")
        })

        observeEvent(settings_Volc$plotly_brush(), {
            req(settings_Volc$plotly_brush())
            brush <- settings_Volc$plotly_brush()
            brush.id <- which(get_de()$EventName %in% brush$key)
            req(brush.id)

            selected <- settings_Volc$selected
            selected <- unique(c(selected, brush.id))
            settings_Volc$selected <- selected
        })


        output$plot_volc <- renderPlotly({
            validate(need(get_se(), "Load Experiment first"))
            validate(need(get_de(), "Load DE Analysis first"))

            selected <- settings_Volc$selected

            num_events <- input$number_events_volc
            res <- as.data.table(get_de()[rows_all(),])
            if(is_valid(input$EventType_volc)) {
                res <- res[get("EventType") %in% input$EventType_volc]
            }
            if(num_events < nrow(res)) {
                res <- res[seq_len(num_events)]
            }

            df.volc <- data.frame(
                EventName = res$EventName, 
                EventType = res$EventType, 
                NMD_direction = res$NMD_direction,
                log2FoldChange = res$log2FoldChange
            )
            if("pvalue" %in% colnames(res)) {
                df.volc$pvalue <- res$pvalue
                df.volc$padj <- res$padj
            } else {
                df.volc$pvalue <- res$P.Value
                df.volc$padj <- res$adj.P.Val
            }

            if(is_valid(selected)) {
                df.volc$selected <- 
                    (df.volc$EventName %in% get_de()$EventName[selected])
            } else {
                df.volc$selected <- FALSE
            }
            if(input$NMD_volc) {
                df.volc <- df.volc[df.volc$NMD_direction != 0, ]
                df.volc$log2FoldChange <- 
                    df.volc$log2FoldChange * df.volc$NMD_direction
            }

            settings_Volc$plot_ini <- TRUE
            if(input$adjP_volc) {
                p <- ggplot(df.volc, aes(
                        x = get("log2FoldChange"), y = -log10(get("padj")),
                        key = get("EventName"), text = get("EventName"), 
                        colour = get("selected")))           
            } else {
                p <- ggplot(df.volc, aes(
                        x = get("log2FoldChange"), y = -log10(get("pvalue")),
                        key = get("EventName"), text = get("EventName"), 
                        colour = get("selected")))               
            }

            p <- p + geom_point() + 
                scale_color_manual(values = c("black", "red"))

            if(input$facet_volc) {
                p <- p + facet_wrap(vars(get("EventType")))
            }
            if(input$NMD_volc) {
                p <- p + labs(x = "Log2 Fold Change NMD substrate")
            } else {
                p <- p + labs(x = "Log2 Fold Change")            
            }
            if(input$adjP_volc) {
                p <- p + labs(y = "Adjusted P Value (-log10)")
            } else {
                p <- p + labs(x = "Nominal P Value (-log10)")            
            }
            
            p <- p + labs(color = "Selected")
            settings_Volc$final_plot <- ggplotly(
                p, tooltip = "text",
                source = "plotly_volcano"
            ) %>% layout(dragmode = "lasso")
            
            print(
                settings_Volc$final_plot
            )
        })
        
        observe({
            shinyFileSave(input, "saveplot_volc", 
                roots = volumes(), session = session,
                filetypes = c("pdf"))
        })
        observeEvent(input$saveplot_volc, {
            req(settings_Volc$final_plot)
            selectedfile <- parseSavePath(volumes(), input$saveplot_volc)
            req(selectedfile$datapath)

            obj <- isolate(settings_Volc$final_plot)
            plotly::orca(obj, .make_path_relative(
                getwd(), selectedfile$datapath))
        })
        
        observeEvent(input$clear_volc, {
            updateSelectInput(session = session, "EventType_volc", 
                selected = NULL)
            shinyWidgets::updateSliderTextInput(session = session, 
                "number_events_volc", selected = 10000)
        })
        
        return(settings_Volc)
    })
}

server_vis_heatmap <- function(
        id, refresh_tab, volumes, get_se, get_de,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Heat <- setreactive_Diag()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
            req(get_se())
            updateSelectInput(session = session, inputId = "anno_col_heat", 
                choices = colnames(colData(get_se())), selected = NULL)
        })

        output$plot_heat <- renderPlotly({
            
            validate(need(get_se(), "Load Experiment first"))
            validate(need(get_de(), "Load DE Analysis first"))

            if(input$select_events_heat == "Highlighted") {
                selected <- rows_selected()
            } else if(input$select_events_heat == "Top N Filtered Results") {
                selected <- rows_all()
                if(length(selected) > input$slider_num_events_heat) {
                    selected <- selected[seq_len(input$slider_num_events_heat)]
                }
            } else {
                selected = seq_len(min(input$slider_num_events_heat, 
                    nrow(get_de())))
            }

            validate(need(length(selected) > 0, "Select some Events first"))

            colData <- as.data.frame(colData(get_se()))

            if(input$mode_heat == "PSI") {
                mat <- make_matrix(get_se(), get_de()$EventName[selected],
                rownames(colData), "PSI")
            } else if(input$mode_heat == "Logit") {
                mat <- make_matrix(get_se(), get_de()$EventName[selected],
                rownames(colData), "logit")
            } else {
                mat <- make_matrix(get_se(), get_de()$EventName[selected],
                rownames(colData), "Z-score")
            }

            validate(need(nrow(mat) > 0 & ncol(mat) > 0, 
                "No data after filtering results"))

            colors.df <- RColorBrewer::brewer.pal.info
            color.index <- which(rownames(colors.df) == input$color_heat)
            color <- grDevices::colorRampPalette(
                rev(RColorBrewer::brewer.pal(
                    n = colors.df$maxcolors[color.index],
                    name = rownames(colors.df)[color.index])
                )
            )

            # Hopefully the fixed filtering in limma pipeline will also fix the 
            #   NA issues here:
            na.exclude <- (rowSums(!is.na(mat)) == 0)
            if(any(na.exclude == TRUE)) {
                output$warning_heat <- renderText({
                    cat("The following have been excluded due to NA values:")
                    paste(rownames(mat)[which(na.exclude)])
                })
                mat <- mat[-which(na.exclude),]
            }

            if(
                    is_valid(input$anno_col_heat) && 
                    all(input$anno_col_heat %in% colnames(colData))
            ) {
                settings_Heat$final_plot <- heatmaply::heatmaply(
                    mat, color = color, 
                    col_side_colors = colData[, input$anno_col_heat, drop=FALSE]
                )
            } else {
                settings_Heat$final_plot <- heatmaply::heatmaply(
                    mat, color = color)
            }      
            print(
                settings_Heat$final_plot
            )
        })
        
        observe({
            shinyFileSave(input, "saveplot_heat", 
                roots = volumes(), session = session,
                filetypes = c("pdf"))
        })
        observeEvent(input$saveplot_heat, {
            req(settings_Heat$final_plot)
            selectedfile <- parseSavePath(volumes(), input$saveplot_heat)
            req(selectedfile$datapath)
            
            obj <- isolate(settings_Heat$final_plot)
            plotly::orca(obj, 
                .make_path_relative(getwd(), selectedfile$datapath))
        })
        
    })
}