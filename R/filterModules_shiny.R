filterModule_UI <- function(id, label = "Counter") {
    ns <- NS(id)
    wellPanel(
        h5(label),  # e.g. "Filter #1"
        selectInput(ns("filterClass"), "Filter Class", 
            width = '100%', choices = c("(none)", "Annotation", "Data")),
        selectInput(ns("filterType"), "Filter Type", 
            width = '100%', choices = c("(none)")),
        conditionalPanel(ns = ns,
        condition = paste0("['TSL'].",
            "indexOf(input.filterType) >= 0"),
            shinyWidgets::sliderTextInput(ns("slider_TSL_min"), 
                "TSL Threshold", 
                choices = seq_len(5), selected = 1)
        ),
        conditionalPanel(ns = ns,
            condition = paste0("['Consistency'].",
                "indexOf(input.filterType) >= 0"),
            shinyWidgets::sliderTextInput(ns("slider_cons_max"), 
                "log-fold maximum", choices = seq(0.2, 5, by = 0.2), 
                selected = 1)
        ),
        conditionalPanel(ns = ns,
            condition = "['Coverage'].indexOf(input.filterType) >= 0",
            sliderInput(ns("slider_cov_min"), "Percent Coverage", 
                min = 0, max = 100, value = 80)
        ),
        conditionalPanel(ns = ns,
            condition = "['Depth'].indexOf(input.filterType) >= 0",
            shinyWidgets::sliderTextInput(ns("slider_depth_min"), 
                "Minimum", choices = c(1,2,3,5,10,20,30,50,100,200,300,500), 
                selected = 20),
        ),
        conditionalPanel(ns = ns,
            condition = "['Depth', 'Coverage'].indexOf(input.filterType) >= 0",
            tagList(
                shinyWidgets::sliderTextInput(ns("slider_mincond"), 
                    "Minimum Conditions Satisfy Criteria", 
                    choices = c(as.character(seq_len(8)), "All"), 
                    selected = "All"),
                selectInput(ns("select_conds"), "Condition", width = '100%',
                    choices = c("(none)")),
                sliderInput(ns("slider_pcTRUE"), 
                    "Percent samples per condition satisfying criteria", 
                    min = 0, max = 100, value = 80)
            )
        ),
        conditionalPanel(ns = ns,
            condition = paste0("['Coverage', 'Consistency'].",
                "indexOf(input.filterType) >= 0"),
            shinyWidgets::sliderTextInput(ns("slider_minDepth"), 
                "Signal Threshold to apply criteria", 
                choices = c(1,2,3,5,10,20,30,50,100,200,300,500), 
                selected = 20),
        ),
        conditionalPanel(ns = ns,
            condition = "['(none)'].indexOf(input.filterClass) < 0",
            selectInput(ns("EventType"), "Splice Type", width = '100%', 
                multiple = TRUE,
                choices = c("IR", "MXE", "SE", "A5SS", "A3SS",
                    "AFE", "ALE", "RI"))
        )
    )
}

filterModule_server <- function(id, filterdata, conditionList) {
    moduleServer(id, function(input, output, session) {
        final <- reactiveValues(
            filterObj = NxtFilter() # initialize to defaults
        )

        # Observe whether colData of NxtSE changes
        fCond <- final$filterObj@condition
        observeEvent(conditionList(), {
            choices_conds <- c("(none)", conditionList())
            if(
                    # Valid condition
                    length(choices_conds) > 1 && is_valid(fCond) && 
                    fCond %in% choices_conds[-1]
            ) {
                updateSelectInput(
                    session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = fCond
                )
            } else if(
                is_valid(fCond) && 
                !(fCond %in% choices_conds)                
            ){
                # If condition is valid but not in column, reset it and return
                updateSelectInput(
                    session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = "(none)"
                )
                return()
            } else {
                updateSelectInput(
                    session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = "(none)"
                )            
            }
        })

        # inputs from final -> UI
        observeEvent(filterdata(), {
            final <- filterdata()
            class_choices <- c("(none)", "Annotation", "Data")
            type_choices <- c("(none)")

            fClass <- final$filterObj@filterClass
            if(is_valid(fClass) && fClass %in% class_choices) {
                if(fClass == "Annotation") {
                    type_choices <- c("Protein_Coding", "NMD", "TSL", 
                        "Terminus", "ExclusiveMXE")
                } else if(fClass == "Data") {
                    type_choices <- c("Depth", "Coverage", "Consistency")
                }
                updateSelectInput(session = session, 
                    inputId = "filterClass", choices = class_choices, 
                    selected = fClass)
            } else {
                # fClass == "" | fClass == "(none)"
                updateSelectInput(session = session, 
                    inputId = "filterClass", choices = class_choices)
                updateSelectInput(session = session, 
                    inputId = "filterType", choices = type_choices)
                return()
            }
            
            fType <- final$filterObj@filterType
            if(is_valid(fType) && fType %in% type_choices) {
                updateSelectInput(session = session, inputId = "filterType", 
                    choices = type_choices, selected = fType)
            } else if(is_valid(fClass) && fClass %in% class_choices) {
                # fClass != "" & fClass != "(none)"
                updateSelectInput(session = session, inputId = "filterType", 
                    choices = type_choices) # Sets default fType if not set
                return()
            } else {
                # Invalid fClass
                updateSelectInput(session = session, 
                    inputId = "filterClass", choices = class_choices)
                updateSelectInput(session = session, 
                    inputId = "filterType", choices = type_choices)
                return()
            }
            
            fMin <- final$filterObj@minimum # always valid
            if(fType == "Depth") {
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_depth_min", 
                    selected = fMin)
            } else if(final$filterType == "Coverage"){
                updateSliderInput(session = session, 
                    inputId = "slider_cov_min", 
                    value = fMin)
            } else  if(final$filterType == "TSL"){
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_TSL_min", 
                    selected = fMin)
            }
            
            fMax <- final$filterObj@maximum # always valid
            shinyWidgets::updateSliderTextInput(
                session = session, inputId = "slider_cons_max", 
                selected = fMax)

            fmDepth <- final$filterObj@minDepth # always valid
            updateSelectInput(session = session, 
                inputId = "slider_minDepth", 
                selected = fmDepth)
            
            fmCond <- final$filterObj@minCond # always valid
            shinyWidgets::updateSliderTextInput(
                session = session, inputId = "slider_mincond", 
                selected = fmCond)
                
            choices_conds = c("(none)", conditionList())
            fCond <- final$filterObj@condition
            if(is_valid(fCond) && fCond %in% choices_conds) {
                updateSelectInput(session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = fCond)
            } else {
                updateSelectInput(session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = "(none)")            
            }
            
            fpcTRUE <- final$filterObj@pcTRUE
            updateSliderInput(session = session, 
                inputId = "slider_pcTRUE", 
                value = fpcTRUE)
            
            feType <- final$filterObj@EventTypes
            eOptions <- c("IR", "MXE", "SE", "A3SS", "A5SS", "ALE", "AFE", "RI")
            
            # make sure feType is always valid
            if(length(feType) > 0) feType <- feType[feType %in% eOptions]
            if(length(feType) == 0) feType <- eOptions
            updateSelectInput(session = session, 
                inputId = "EventType", 
                selected = feType)
        })

        # outputs from UI -> final
        observeEvent(input$filterClass, {
            final$filterObj@filterClass <- input$filterClass
            if(input$filterClass == "Annotation") {
                type_choices <- c("Protein_Coding", "NMD", "TSL", 
                    "Terminus", "ExclusiveMXE")
            } else if(input$filterClass == "Data") {
                type_choices <- c("Depth", "Coverage", "Consistency")
            } else {
                type_choices <- "(none)"
            }
            cur_choice <- isolate(final$filterType)
            if(is_valid(cur_choice) && cur_choice %in% type_choices) {
                updateSelectInput(session = session, 
                    inputId = "filterType", 
                    choices = type_choices, selected = cur_choice)
            } else {
                final$filterObj@filterType <- type_choices[1]
                updateSelectInput(session = session, 
                    inputId = "filterType", 
                    choices = type_choices)
            }
        })
        observeEvent(input$filterType, {
            # final$trigger = NULL
            req(input$filterType)
            fType <- input$filterType
            final$filterObj@filterType <- fType

            fMin <- final$filterObj@minimum
            if(fType == "Depth") {
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_depth_min", 
                    selected = fMin)
            } else if(final$filterType == "Coverage"){
                updateSliderInput(session = session, 
                    inputId = "slider_cov_min", 
                    value = fMin)
            } else  if(final$filterType == "TSL"){
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_TSL_min", 
                    selected = fMin)
            }
        })
        observeEvent(input$slider_depth_min, {
            if(final$filterObj@filterType == "Depth") {
                final$filterObj@minimum = input$slider_depth_min
            }
        })
        observeEvent(input$slider_cov_min, {
            if(final$filterObj@filterType == "Coverage"){
                final$filterObj@minimum = input$slider_cov_min
            }
        })
        observeEvent(input$slider_TSL_min,{        
            if(final$filterObj@filterType == "TSL"){
                final$filterObj@minimum = as.numeric(input$slider_TSL_min)
            }
        })
        observeEvent(input$slider_cons_max,{        
            final$filterObj@maximum = input$slider_cons_max
        })
        observeEvent(input$slider_minDepth,{        
            final$filterObj@minDepth = input$slider_minDepth
        })
        observeEvent(input$slider_mincond,{        
            final$filterObj@minCond = input$slider_mincond
        })
        observeEvent(input$select_conds,{        
            final$filterObj@condition = input$select_conds
        })
        observeEvent(input$slider_pcTRUE,{        
            final$filterObj@pcTRUE = input$slider_pcTRUE
        })
        observeEvent(input$EventType,{        
            final$filterObj@EventTypes = input$EventType
        })

        # Returns filter list from module
        return(final)
    })
}