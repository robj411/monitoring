shinyUI(
  bootstrapPage(
    theme = shinythemes::shinytheme("spacelab"),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      tags$link(rel = "shortcut icon", type = "image/x-icon",
                href = "favicon.ico")
    ),
    # include the js code
    shinyjs::useShinyjs(),
    includeScript("scripts.js"),
    column(12,
           HTML(
             "
             <meta name='keywords' content='infectious,disease,epidemiology,genome,coronavirus,SARS-CoV-2,phylodynamic,phylogenetics,lineage,mutation,Spike,UK'>
             <meta name='author' content='Rob Johnson'>
             <h1>COG-UK variants</h1>
             "
           ),
      #tabsetPanel(
        #tabPanel("About",
          HTML(
             "
             <p>Recent frequencies of variants with prevalence >2.5% (at positions with missingness <5%) in the last six weeks.</p>
             "
             #)
        #)
        ),
    tags$hr(style="background: #cccccc; height: 4px;")
    ),
    column(12, id = "plot",
           tabsetPanel(
             tabPanel("Variants summarised",
                      #downloadButton("save_variant_table", "Save table"),
                      HTML(
                        "
             <p>Log odds = log odds of the new variant. Hospital log odds ratio = log odds in hospital - log odds in lighthouse. Female log odds ratio = female log odds - male log odds.</p>
             <p>Mean recent growth = smoothed log odds 3 weeks ago - smoothed log odds 13 weeks ago. </p>
             "
                      ),
                      #tableOutput("estimated_r_output")
                      fluidRow(
                        #checkboxInput(inputId='show_all_sequences', label='Show all sequences', value=F),
                        DT::dataTableOutput(outputId = 'variants')
                        , style = "font-size:15px", align="left")
             ),
             tabPanel("Trends",
                                 fluidRow(
                                   column(3,
                                          selectInput(inputId='ti_genes', label='Gene', choices = NULL, selected = NULL)
                                   ),
                                   column(3,
                                          selectInput(inputId='ti_variants', label='Position', choices = NULL, selected = NULL)
                                   )),
                      tabsetPanel(
                        tabPanel("Log odds",
                                 column(3, id = "menu",
                                        div(id = "1.2",
                                            radioButtons(inputId='ti_lh', label='Labs', choices = c('Hospital','Lighthouse','Both'), selected ='Both')
                                            , selectInput(inputId='ti_regions', label='Region', choices = NULL, selected = NULL)
                                        )
                                        
                                 ),
                                 #         downloadButton("save_plot", "Save Image"),
                                 column(6, fluidRow( plotOutput( 'logodds', width = "100%", height = "300px")),
                                        fluidRow( plotOutput( 'trend', width = "100%", height = "300px")),
                                        fluidRow( plotOutput( 'start', width = "100%", height = "300px"))),
                                 column(3,plotOutput('legend'))
                        ),
                        tabPanel("Frequencies by region",
                                 column(5,fluidRow(
                                   checkboxGroupInput('ti_region2', label = 'Regions', choices=NULL,selected=NULL)
                                 )),
                                 column(6,
                                        fluidRow( plotOutput( 'trend2', width = "100%", height = "300px")))
                                 ))),
             
             tabPanel("Genes summarised",
                      #downloadButton("save_gene_table", "Save table"),
                      #tableOutput("estimated_r_output")
                      fluidRow(
                        #checkboxInput(inputId='show_all_sequences', label='Show all sequences', value=F),
                        DT::dataTableOutput(outputId = 'genes')
                        , style = "font-size:15px", align="left")
             )
           )
    )
  )
)
