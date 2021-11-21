ui_ref_new <- function(id) {
    ns <- NS(id)
    fluidRow(
        .ui_notice(),
        column(6,
            # h4("Select Reference Directory"),
            tags$div(title = paste("Specify (or create) a directory for NxtIRF",
                    "to create its IRFinder/NxtIRF reference"),
                wellPanel(
                    shinyDirButton(ns("dir_reference_path"), 
                        label = "Select Reference Directory", 
                        title = "Select Reference Directory",
                        buttonType = "primary"
                    ), textOutput(ns("txt_reference_path"))
                )
            ), br(),
            wellPanel(
                h4("Select Ensembl Reference"),
                selectInput(ns('release'), 'Select Release', width = '100%',
                    choices = c("")),
                selectInput(ns('species'), 'Select Species', width = '100%',
                    choices = c("")),
                selectInput(ns('fasta'), 'Select Genome', width = '100%',
                    choices = c("")),
                selectInput(ns('gtf'), 'Select Annotation', width = '100%',
                    choices = c(""))
            ),
            wellPanel(
                h4("or select Reference from File"),
                tags$div(title = paste("Choose a user-supplied genome fasta file"),
                    shinyFilesButton(ns("file_genome"), 
                    label = "Choose genome FASTA File", 
                    title = "Choose genome FASTA File", 
                    buttonType = "primary",
                    multiple = FALSE)
                ),
                textOutput(ns("txt_genome")), br(),
                tags$div(title = paste("Choose a user-supplied transcript",
                    "reference gtf file"),
                    shinyFilesButton(ns("file_gtf"), 
                    label = "Choose transcriptome GTF File", 
                    title = "Choose transcriptome GTF File", 
                    buttonType = "primary",
                    multiple = FALSE)
                ),
                textOutput(ns("txt_gtf")),
            ),
            actionButton(ns("load_ref_example"), "Load Example FASTA / GTF")
        ),
        column(6,
            wellPanel(
                tags$div(title = paste("NxtIRF will auto-populate default",
                        "mappability and non-polyA reference files for",
                        "hg38, hg19, mm10 and mm9 genomes"),
                    selectInput(ns('newref_genome_type'),
                        'Select Genome Type to set Mappability and non-PolyA files',
                        c("(custom)", "hg38", "mm10", "hg19", "mm9"))
                ),
                tags$div(title = paste("Select Mappability Exclusion file.",
                        "This is typically a 3 columns",
                        "of values containing seqnames,",
                        "start and end coordinates of low-mappability regions"),
                    shinyFilesButton(ns("file_mappa"), 
                        label = "Choose Mappability Exclusion file", 
                        title = "Choose Mappability Exclusion file", 
                        buttonType = "success",
                        multiple = FALSE),
                    actionButton(ns("clear_mappa"), "Clear",
                        class = "btn-outline-danger")
                ), 
                textOutput(ns("txt_mappa")), br(),
            
                tags$div(title = paste("Select Non-PolyA reference file.",
                        "This is used by IRFinder",
                        "to calculate reads from known non-polyadenylated",
                        "transcripts to assess",
                        "quality of poly-A enrichment in sample QC"),
                    shinyFilesButton(ns("file_NPA"), 
                        label = "Choose non-PolyA BED file", 
                        title = "Choose non-PolyA BED file", 
                        buttonType = "success",
                        multiple = FALSE),
                    actionButton(ns("clear_NPA"), "Clear",
                        class = "btn-outline-danger")
                ), 
                textOutput(ns("txt_NPA")), br(),
                
                tags$div(title = paste("Select Blacklist file.",
                        "This is typically a 3 columns",
                        "of values containing seqnames, start and end coordinates",
                        "of regions to exclude from IRFinder analysis"),
                    shinyFilesButton(ns("file_bl"), 
                        label = "Choose blacklist BED file", 
                        title = "Choose blacklist BED file", 
                        buttonType = "success",
                        multiple = FALSE),
                    actionButton(ns("clear_bl"), "Clear",
                        class = "btn-outline-danger")
                ), 
                textOutput(ns("txt_bl"))
            ),
            actionButton(ns("buildRef"), "Build Reference", 
                class = "btn-success"),
            actionButton(ns("clearNewRef"), "Clear settings",
                class = "btn-outline-danger"),
            br(),
            uiOutput(ns("refStatus"))
        )
    )
}
