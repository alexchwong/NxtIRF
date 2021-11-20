server_cov <- function(
        id, refresh_tab, volumes, get_se, get_de,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Cov <- setreactive_Cov()
        get_ref <- function(){
            se = get_se()
            ref(se)
        }
        observeEvent(refresh_tab(), {
            req(refresh_tab())
            output$warning_cov <- renderText({
                validate(need(get_se(), "Please build experiment first"))
            })
            req(get_se())
            .server_cov_refresh(
                session, get_ref()$gene_list,
                get_de(), rows_all(), rows_selected(),
                input$slider_num_events_cov, input$events_cov,
                input$select_events_cov
            )
            
            se <- get_se()
            settings_Cov$event.ranges <- as.data.table(
                CoordToGR(rowData(se)$EventRegion))
            settings_Cov$event.ranges$EventName <- rowData(se)$EventName
            
            .server_cov_refresh_tracks_cov(session, input$mode_cov, 
                input$condition_cov, se)
        })
        
        # Delayed (debounced) reactives
        chr_r <- reactive({
            req(is_valid(input$chr_cov))
            req(input$chr_cov %in% names(get_ref()$seqInfo))
            
            input$chr_cov
        })
        start_r <- reactive({
            req(input$start_cov)
            
            input$start_cov
        })
        end_r <- reactive({
            req(input$end_cov)
            
            input$end_cov
        })
        tracks_r <- reactive({
            server_cov_get_all_tracks(input)
        })
        trigger_r <- reactive({
            settings_Cov$trigger
        })
        
        chr_rd <- chr_r %>% debounce(1000)
        start_rd <- start_r %>% debounce(1000)
        end_rd <- end_r %>% debounce(1000)
        tracks_rd <- tracks_r %>% debounce(3000)
        trigger_rd <- trigger_r %>% debounce(1000)
    
        observeEvent(list(chr_rd(), start_rd(), end_rd()), {
            .server_cov_update_norm_event(
                input, session, settings_Cov$event.ranges)
        })
        observeEvent(list(
            input$refresh_coverage, input$stack_tracks_cov,
            input$graph_mode_cov, input$pairwise_t_cov, input$condense_cov
        ), {
            settings_Cov$trigger <- runif(1)
        })
        
        observeEvent(list(trigger_rd(), tracks_rd()), {
            tracks <- tracks_r()
            settings_Cov$plot_params <- .server_cov_refresh_plot_args(
                get_se(), get_ref(), 
                input$event_norm_cov, 
                input$chr_cov, input$start_cov, 
                input$end_cov, tracks, 
                settings_Cov$plot_params, input
            )
            if(.server_cov_check_plot_args(settings_Cov$plot_params)) {
                obj <- do.call(Plot_Coverage, settings_Cov$plot_params)            
            
                req(obj)
                settings_Cov$final_plot <- obj$final_plot
                settings_Cov$final_plot$x$source <- "plotly_ViewRef"
                output$plot_cov <- renderPlotly({
                    settings_Cov$plot_ini <- TRUE      
                    print(
                        settings_Cov$final_plot
                    )
                })
            }
        })
        
        observeEvent(input$graph_mode_cov, {
            req(settings_Cov$plot_ini == TRUE)
            .server_cov_plot_change_mode(session, input$graph_mode_cov)
        })
        observeEvent(input$mode_cov, {
            .server_cov_refresh_track_condition(
                session, input$mode_cov, get_se())
            .server_cov_refresh_tracks_cov(
                session, input$mode_cov, input$condition_cov, get_se())
        })
        observeEvent(input$condition_cov, {
            .server_cov_refresh_tracks_cov(session, input$mode_cov, 
                input$condition_cov, get_se())
        })
        settings_Cov$plotly_relayout = reactive({
            req(settings_Cov$plot_ini == TRUE)
            event_data("plotly_relayout", source = "plotly_ViewRef")
        })
        observeEvent(settings_Cov$plotly_relayout(), {
            print(settings_Cov$plotly_relayout())
            req(length(settings_Cov$plotly_relayout()) == 2)
            req(all(c("xaxis.range[0]", "xaxis.range[1]") %in% 
                names(settings_Cov$plotly_relayout())))
                
            updateTextInput(
                session = session, inputId = "start_cov", 
                value = max(1, 
                    round(settings_Cov$plotly_relayout()[["xaxis.range[0]"]]))
                )
            updateTextInput(session = session, inputId = "end_cov", 
                value = round(
                    settings_Cov$plotly_relayout()[["xaxis.range[1]"]]))
        })
        observeEvent(input$zoom_out_cov, {
            req(input$zoom_out_cov, input$chr_cov) 
            req(input$chr_cov %in% names(get_ref()$seqInfo))
            
            seqInfo <- get_ref()$seqInfo[input$chr_cov]
            .server_cov_zoom_out(input, session, seqInfo)
        })
        observeEvent(input$zoom_in_cov, {
            req(input$zoom_in_cov)
            
            .server_cov_zoom_in(input, session)
        })
        observeEvent(input$events_cov, {
            req(input$events_cov)
            req(input$events_cov != "(none)")
            
            events_id_view <- settings_Cov$event.ranges[
                get("EventName") == input$events_cov]
            .server_cov_locate_events(input, session, events_id_view)
        })
        observeEvent(input$genes_cov, {
            req(input$genes_cov)
            req(input$genes_cov != "(none)")
            
            gene_id_view <- get_ref()$gene_list[
                get("gene_display_name") == input$genes_cov]
            .server_cov_locate_genes(input, session, gene_id_view)
        })
        
        observeEvent(chr_rd(), {
            seqInfo <- get_ref()$seqInfo[chr_rd()]
            seqmax <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            req(seqmax > 50)
            .server_cov_change_start_end(input, session, output, seqmax)
            settings_Cov$trigger <- runif(1)
        })
        observeEvent(list(start_rd(), end_rd()), {
            req(input$chr_cov, input$chr_cov %in% names(get_ref()$seqInfo))
            seqInfo <- get_ref()$seqInfo[input$chr_cov]
            seqmax <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            req(seqmax > 50)
            output <- .server_cov_change_start_end(
                input, session, output, seqmax)
            settings_Cov$trigger <- runif(1)
        })
        
        # Populate events
        observeEvent(input$select_events_cov, {
            req(rows_all())
            req(get_de())
            
            .server_cov_change_event_list(
                session, input$select_events_cov, 
                input$slider_num_events_cov,
                get_de(), rows_all(), rows_selected()
            )
        })
        observe({
            shinyFileSave(input, "saveplot_cov", roots = volumes(), 
                session = session, filetypes = c("pdf"))    
        })
        observeEvent(input$saveplot_cov, {    
            req(settings_Cov$final_plot)
            selectedfile <- parseSavePath(volumes(), input$saveplot_cov)
            req(selectedfile$datapath)
            plotly::orca(settings_Cov$final_plot, 
                .make_path_relative(getwd(), selectedfile$datapath),
                width = 1920, height = 1080)
        })
    
    })
}

# Get the i'th track as character
.server_cov_get_track_selection <- function(input, i) {
    return(input[[paste0("track", as.character(i), "_cov")]])
}

# Get a list of all tracks
server_cov_get_all_tracks <- function(input) {
    tracks = list()
    for(i in seq_len(4)) {
        tracks[[i]] = .server_cov_get_track_selection(input, i)       
    }
    tracks <- Filter(is_valid, tracks)
    return(tracks)
}

# Updates drop-downs
.server_cov_refresh <- function(
        session, gene_list, DE,
        rows_all, rows_selected, num_events, selected_event, mode
) {
    if(!is.null(gene_list)) {
        message("Populating drop-down box with ", 
            length(unique(gene_list$gene_display_name)), " genes")
        updateSelectInput(session = session, inputId = "chr_cov", 
            choices = c("(none)", 
                as.character(sort(unique(gene_list$seqnames)))),
            selected = "(none)")
        updateSelectizeInput(session = session, inputId = "genes_cov",
            server = TRUE, choices = c("(none)", gene_list$gene_display_name), 
            selected = "(none)")
    } else {
        updateSelectInput(session = session, inputId = "chr_cov", 
            choices = c("(none)"), selected = "(none)")
        updateSelectizeInput(session = session, inputId = "genes_cov", 
            server = TRUE, choices = c("(none)"), selected = "(none)") 
    }
    if(is_valid(DE)) {
        if(mode == "Highlighted") {
            selected <- rows_selected
        } else if(mode == "Top N Filtered Results") {
            selected <- rows_all
        } else {
            selected <- seq_len(nrow(DE))
        }
        if(length(selected) > num_events) {
            selected <- selected[seq_len(num_events)]
        }
        if(length(selected) > 0 & is_valid(DE)) {
            if(is_valid(selected_event)) {
                if(!(selected_event %in% DE$EventName[selected])) {
                    selected_event = "(none)"
                }
            } else {
                selected_event = "(none)"
            }
            updateSelectizeInput(session = session, 
                inputId = "events_cov", server = TRUE,
                choices = c("(none)", 
                    DE$EventName[selected]), 
                    selected = selected_event)
        } else {
            updateSelectizeInput(session = session, 
                inputId = "events_cov", server = TRUE,
                choices = c("(none)"), selected = "(none)")
        }
    }
}

# Get a list of all EventName in the given genomic region
.server_cov_get_inrange_events <- function(
        view_chr, view_start, view_end, event.ranges
) {
    req(event.ranges)
    DT <- event.ranges[get("seqnames") == view_chr &
        get("end") > view_start & get("start") < view_end]
    DT$EventName
}

# Updates dropdown of Event Norm options
.server_cov_update_norm_event <- function(input, session, event.ranges) {
    view_chr <- isolate(input$chr_cov)
    view_start <- isolate(input$start_cov)
    view_end <- isolate(input$end_cov)
    selected_event <- isolate(input$events_cov)
    cur_event <- isolate(input$event_norm_cov)
    req(view_chr)
    req(view_start)
    req(view_end)
    
    event_choices <- c("(none)")
    if(is_valid(selected_event)) {
        event_choices <- c(event_choices, selected_event)
    } else if(is_valid(cur_event)) {
        event_choices <- c(event_choices, cur_event)        
    }
    event_choices = unique(c(event_choices, 
        .server_cov_get_inrange_events(view_chr, view_start, view_end,
            event.ranges)))
            
    if(is_valid(selected_event)) {
        updateSelectInput(session = session, inputId = "event_norm_cov", 
            choices = event_choices, selected = selected_event)
    } else if(is_valid(cur_event)) {
        updateSelectInput(session = session, inputId = "event_norm_cov", 
            choices = event_choices, selected = cur_event)
    } else {
        updateSelectInput(session = session, inputId = "event_norm_cov", 
            choices = event_choices, selected = "(none)")        
    }
}

# Compiles a list of arguments to pass into Plot_Coverage
.server_cov_refresh_plot_args <- function(
        se, ref, norm_event, 
        view_chr, view_start, view_end, 
        tracks, plot_params, input
) {
    req(view_chr, view_start, view_end, se)

    args <- list(
        se = se, 
        seqname = view_chr, start = view_start, end = view_end, 
        strand = input$strand_cov, 
        norm_event = norm_event,
        zoom_factor = 0, bases_flanking = 0,
        tracks = tracks, # track_names = "",
        condition = input$condition_cov, 
        
        condense_tracks = input$condense_cov,
        stack_tracks = input$stack_tracks_cov,
        t_test = input$pairwise_t_cov

        # graph_mode = input$graph_mode_cov,
    )
    args <- Filter(is_valid, args)
    
    return(args)
}

.server_cov_check_plot_args <- function(args) {
    if(length(args$tracks) == 0) return(FALSE)
    if(!is_valid(args$condition)) {
        check <- tryCatch({
            .plot_cov_validate_args(
                se = args$se, 
                tracks = args$tracks, 
                condition = args$condition,
                seqname = args$seqname, 
                start = args$start, 
                end = args$end
            )
            TRUE
        }, error = function(e) FALSE)
    } else {
        check <- tryCatch({
            .plot_cov_validate_args(
                se = args$se, 
                tracks = args$tracks, 
                # condition = args$condition,
                seqname = args$seqname, 
                start = args$start, 
                end = args$end
            )
            TRUE
        }, error = function(e) FALSE)

    }
    return(check)
}

# Change plotly mode (Pan / Zoom / Movable Labels)
.server_cov_plot_change_mode <- function(session, mode) {
    if(mode == "Pan") {
        plotlyProxy("plot_cov", session) %>% 
            plotlyProxyInvoke("relayout", list(dragmode = "pan")) %>%
            plotlyProxyInvoke("reconfig", editable = FALSE)
    } else if(mode == "Zoom") {
        plotlyProxy("plot_cov", session) %>% 
            plotlyProxyInvoke("relayout", list(dragmode = "zoom")) %>%
            plotlyProxyInvoke("reconfig", editable = FALSE)
    } else if(mode == "Movable Labels") {
        plotlyProxy("plot_cov", session) %>% 
            plotlyProxyInvoke("relayout", list(dragmode = FALSE)) %>%
            plotlyProxyInvoke("reconfig", editable = TRUE)
    }
}

# Refresh list of available conditions
.server_cov_refresh_track_condition <- function(session, mode, se) {
    req(se)
    if(mode == "By Condition") {
        colData <- colData(se)
        updateSelectInput(session = session, inputId = "condition_cov", 
            choices = c("(none)", colnames(colData)))
        updateSelectInput(session = session, inputId = "track1_cov", 
            choices = c("(none)"), selected = "(none)")     
        updateSelectInput(session = session, inputId = "track2_cov", 
            choices = c("(none)"), selected = "(none)")  
        updateSelectInput(session = session, inputId = "track3_cov", 
            choices = c("(none)"), selected = "(none)")    
        updateSelectInput(session = session, inputId = "track4_cov", 
            choices = c("(none)"), selected = "(none)")             
    } else {
        updateSelectInput(session = session, inputId = "condition_cov", 
            choices = c("(none)"))
    }
}

# Refresh list of available tracks
.server_cov_refresh_tracks_cov <- function(session, mode, condition, se) {
    req(se)
    if(mode == "By Condition") {
        if(is_valid(condition)) {
            colData <- colData(se)
            updateSelectInput(session = session, inputId = "track1_cov", 
                choices = c("(none)", 
                unique(as.character(unlist(colData[, condition])))))    
            updateSelectInput(session = session, inputId = "track2_cov", 
                choices = c("(none)", 
                unique(as.character(unlist(colData[, condition])))))    
            updateSelectInput(session = session, inputId = "track3_cov", 
                choices = c("(none)", 
                unique(as.character(unlist(colData[, condition])))))    
            updateSelectInput(session = session, inputId = "track4_cov", 
                choices = c("(none)", 
                unique(as.character(unlist(colData[, condition])))))    
        }  else {
            updateSelectInput(session = session, inputId = "track1_cov", 
                choices = "(none)")
            updateSelectInput(session = session, inputId = "track2_cov", 
                choices = "(none)")
            updateSelectInput(session = session, inputId = "track3_cov", 
                choices = "(none)")
            updateSelectInput(session = session, inputId = "track4_cov", 
                choices = "(none)")
        }
    } else {
        avail_samples <- names(covfile(se)[file.exists(covfile(se))])
        updateSelectInput(session = session, inputId = "track1_cov", 
            choices = c("(none)", avail_samples), selected = "(none)")
        updateSelectInput(session = session, inputId = "track2_cov", 
            choices = c("(none)", avail_samples), selected = "(none)")
        updateSelectInput(session = session, inputId = "track3_cov", 
            choices = c("(none)", avail_samples), selected = "(none)")
        updateSelectInput(session = session, inputId = "track4_cov", 
            choices = c("(none)", avail_samples), selected = "(none)")
    }
}

# Zoom out
.server_cov_zoom_out <- function(input, session, seqInfo) {
    view_start  <- input$start_cov
    view_end    <- input$end_cov
    req(view_start, view_end, view_end - view_start >= 50)
    
    seqmax      <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
    # get center of current range
    center      <- round((view_start + view_end) / 2)
    span        <- view_end - view_start
    # zoom range is 50 * 3^z
    cur_zoom    <- floor(log(span/50) / log(3))

    new_span <- round(span * 3)
    # if(new_span > seqmax - 1) new_span = seqmax - 1
    new_start <- max(1, center - round(new_span / 2))
    
    updateTextInput(session = session, inputId = "start_cov", 
        value = new_start)
    updateTextInput(session = session, inputId = "end_cov", 
        value = new_start + new_span)
}

# Zoom in
.server_cov_zoom_in <- function(input, session) {
    view_start  <- input$start_cov
    view_end    <- input$end_cov
    req(view_start, view_end, view_end - view_start >= 50)
    
    # get center of current range
    center      <- round((view_start + view_end) / 2)
    span        <- view_end - view_start
    # zoom range is 50 * 3^z
    cur_zoom    <- floor(log(span/50) / log(3))

    new_span <- round(span / 3)
    if(new_span < 50) new_span <- 50
    new_zoom <- floor(log(new_span / 50) / log(3))
    new_start <- max(1, center - round(new_span / 2))

    updateTextInput(session = session, inputId = "start_cov", 
        value = new_start)
    updateTextInput(session = session, inputId = "end_cov", 
        value = new_start + new_span)
}

# Sets start and end given an Event
.server_cov_locate_events <- function(input, session, events_id_view) {
        
    span        <- events_id_view$end[1] - events_id_view$start[1]
    view_start  <- max(1, events_id_view$start[1] - span)
    view_end    <- view_start + 3 * span
    
    updateSelectInput(session = session, inputId = "chr_cov", 
        selected = events_id_view$seqnames[1])
    updateTextInput(session = session, inputId = "start_cov", 
        value = view_start)
    updateTextInput(session = session, inputId = "end_cov", 
        value = view_end)
}

# Sets start and end given a Gene
.server_cov_locate_genes <- function(input, session, gene_id_view) {
    updateSelectInput(session = session, inputId = "chr_cov", 
        selected = gene_id_view$seqnames[1])
    updateTextInput(session = session, inputId = "start_cov", 
        value = gene_id_view$start[1])
    updateTextInput(session = session, inputId = "end_cov", 
        value = gene_id_view$end[1])
}

# Changes start and end coordinates, if over seqmax
# also if changing to shorter chromosome and prev end > chrom length
.server_cov_change_start_end <- function(input, session, output, seqmax) {
    target_start    <- input$start_cov
    target_end      <- input$end_cov
    req(target_end, target_start)
    
    if(target_end > seqmax) target_end <- seqmax
    
    # assumes chromosome length > 50
    if(target_end - target_start < 50) target_start <- target_end - 50
    
    span <- target_end - target_start
    cur_zoom = floor(log(span/50) / log(3))
    
    output$label_zoom_cov <- renderText({16 - cur_zoom})
    updateTextInput(session = session, inputId = "end_cov", 
        value = target_end)
    updateTextInput(session = session, inputId = "start_cov", 
        value = target_start)
    return(output)
}

# Updates Event List
.server_cov_change_event_list <- function(
        session, mode, num_events,
        DE, rows_all, rows_selected
) {
    if(mode == "Highlighted") {
        selected <- rows_selected
    } else if(mode == "Top N Filtered Results") {
        selected <- rows_all
        if(length(selected) > num_events) {
            selected <- selected[seq_len(num_events)]
        }
    } else {
        selected <- seq_len(min(num_events, nrow(DE)))
    }
    
    if(length(selected) > 0 & is_valid(DE)) {
        updateSelectizeInput(
            session = session, inputId = "events_view", 
            server = TRUE, 
            choices = c("(none)", DE$EventName[selected]), 
            selected = "(none)"
        )
        updateSelectizeInput(
            session = session, inputId = "events_cov", 
            server = TRUE, 
            choices = c("(none)", DE$EventName[selected]), 
            selected = "(none)"
        )
    } else {
        updateSelectizeInput(session = session, inputId = "events_view", 
            server = TRUE, choices = c("(none)"), selected = "(none)")
        updateSelectizeInput(session = session, inputId = "events_cov", 
            server = TRUE, choices = c("(none)"), selected = "(none)")
    }
}

################################################################################

# Validate given arguments in Plot_Coverage()
.plot_cov_validate_args <- function(se, tracks, condition,
        Event, Gene,
        seqname, start, end, bases_flanking
) {
    .plot_cov_validate_args_se(se, tracks, condition)
    cov_data <- ref(se)
    checked <- .plot_cov_validate_args_loci(
        cov_data, Event, Gene, seqname, start, end)
    if (!checked) .plot_cov_validate_args_event(se, Event, bases_flanking)
}

# Check se, tracks, conditions
.plot_cov_validate_args_se <- function(se, tracks, condition) {
    if (missing(se) || !is(se, "NxtSE"))
        .log("In Plot_Coverage, `se` must be a valid `NxtSE` object")

    if (!all(file.exists(covfile(se))))
        .log(paste("In Plot_Coverage,",
            "COV files are not defined in se.",
            "Please supply the correct paths of the COV files",
            "using covfile(se) <- vector_of_correct_COVfile_paths"))

    # Check condition and tracks
    if (length(tracks) < 1 | length(tracks) > 4)
        .log(paste("In Plot_Coverage,", "tracks must be of length 1-4"))

    if (!missing(condition)) {
        if (length(condition) != 1)
            .log(paste("In Plot_Coverage,", "condition must be of length 1"))

        if (!(condition %in% names(colData(se))))
            .log(paste("In Plot_Coverage,",
                "condition must be a valid column name in colData(se)"))

        condition_options <- unique(colData(se)[, condition])
        if (!all(tracks %in% condition_options))
            .log(paste("In Plot_Coverage,",
                "some tracks do not match valid condition names in",
                args[["condition"]]))

    } else {
        if (!all(tracks %in% colnames(se)))
            .log(paste("In Plot_Coverage,",
                "some tracks do not match valid sample names in se"))
    }
}

# Checks Gene and loci. Only this is run if Plot_Genome is run
.plot_cov_validate_args_loci <- function(cov_data,
    Event, Gene, seqname, start, end, bases_flanking = 0
) {
    if (!all(c("seqInfo", "gene_list", "elem.DT", "transcripts.DT") %in%
            names(cov_data)))
        .log(paste("In Plot_Coverage,",
            "cov_data must be a valid object",
            "created by prepare_covplot_data()"))

    # Check we know where to plot
    if (missing(Event) & missing(Gene) &
            (missing(seqname) | missing(start) | missing(end))
    ) {
        .log(paste("In Plot_Coverage,",
            "Event or Gene cannot be empty, unless coordinates are provided"))
    } else if ((is_valid(seqname) & is_valid(start) & is_valid(end))) {
        view_chr <- as.character(seqname)
        view_start <- start
        view_end <- end
    } else if (is_valid(Gene)) {
        if (!(Gene %in% cov_data$gene_list$gene_id) &
                !(Gene %in% cov_data$gene_list$gene_name)) {
            .log(paste("In Plot_Coverage,",
                Gene, "is not a valid gene symbol or Ensembl gene id"))
        }
        if (!(Gene %in% cov_data$gene_list$gene_id)) {
            gene.df <- as.data.frame(
                cov_data$gene_list[get("gene_name") == get("Gene")])
            if (nrow(gene.df) != 1) {
                .log(paste("In Plot_Coverage,", Gene,
                    "is an ambiguous name referring to 2 or more genes.",
                    "Please provide its gene_id instead"))
            }
        } else {
            gene.df <- as.data.frame(
                cov_data$gene_list[get("gene_id") == get("Gene")])
        }
        view_chr <- as.character(gene.df$seqnames)
        view_start <- gene.df$start
        view_end <- gene.df$end
    } else {
        return(FALSE)
    }
    view_center <- (view_start + view_end) / 2
    view_length <- view_end - view_start
    if (!(view_chr %in% names(cov_data$seqInfo)))
        .log(paste("In Plot_Coverage,", view_chr,
            "is not a valid chromosome reference name in the given genome"))

    if (is_valid(bases_flanking) &&
            (!is.numeric(bases_flanking) || bases_flanking < 0))
        .log(paste("In Plot_Coverage,",
            "bases_flanking must be a non-negative number"))

    if (!is.numeric(view_length) || view_length < 0)
        .log(paste("In Plot_Coverage,",
            "view_length must be a non-negative number"))

    return(TRUE)
}

# Checks whether Event given is valid.
.plot_cov_validate_args_event <- function(se, Event, bases_flanking) {
    cov_data <- ref(se)
    rowData <- as.data.frame(rowData(se))
    if (!(Event %in% rownames(rowData))) {
        .log(paste("In Plot_Coverage,", Event,
            "is not a valid IR or alternate splicing event in rowData(se)"))
    }
    rowData <- rowData[Event, ]
    view_chr <- tstrsplit(rowData$EventRegion, split = ":")[[1]]
    temp1 <- tstrsplit(rowData$EventRegion, split = "/")
    temp2 <- tstrsplit(temp1[[1]], split = ":")[[2]]
    view_start <- as.numeric(tstrsplit(temp2, split = "-")[[1]])
    view_end <- as.numeric(tstrsplit(temp2, split = "-")[[2]])

    view_center <- (view_start + view_end) / 2
    view_length <- view_end - view_start
    if (!(view_chr %in% names(cov_data$seqInfo)))
        .log(paste("In Plot_Coverage,", view_chr,
            "is not a valid chromosome reference name in the given genome"))

    if (is_valid(bases_flanking) &&
            (!is.numeric(bases_flanking) || bases_flanking < 0))
        .log(paste("In Plot_Coverage,",
            "bases_flanking must be a non-negative number"))

    if (!is.numeric(view_length) || view_length < 0)
        .log(paste("In Plot_Coverage,",
            "view_length must be a non-negative number"))

    return(TRUE)
}
