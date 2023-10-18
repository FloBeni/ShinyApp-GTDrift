ui <- shinyUI(fluidPage(style = "background-color:#a7c2da;height: 100%;font-family: Economica; font-size: 21px",tags$style(HTML( #couleur des onglets
  "
    .tabbable > .nav > li > a {background-color: #529bde; font-size: 23px; font-family: 'Economica' ; color:white}
    .tabbable > .nav > li[class=active]  > a {background-color: #136dc0; font-size: 25px; font-family: 'Economica' ; color:white}
    .my-class {font-size: 21px;}
    ")),
  setSliderColor(c("#529bde","#529bde","#529bde","#529bde","#529bde","#529bde"), c(1, 2, 3, 4, 5,6)),shinyjs::useShinyjs(),
  tabsetPanel(id = "tabs",
              tabPanel("Inter-species graphics",
                       column(12,
                              column(1,offset = 1,h2("Y axis :")),
                              column(10,selectInput("y_inter", "",choices = axisInter_list_quantitative)),
                              column(10,class = "well",offset = 0,withSpinner(type = 4,color ="#136dc0",plotlyOutput("plot_inter", width = "100%", height = "100%"))),
                              column(2,class="well", 
                                     fluidRow(column(4,
                                                     materialSwitch(inputId = "boxplot_inter", label = h3("Boxplot"),status="primary")),
                                              column(5,
                                                     prettyCheckboxGroup( "scale_inter", shape = "round",animation="pulse",bigger=T,
                                                                          status = "primary",
                                                                          fill = TRUE,h3("Scales"),choices = list("X log10","Y log10"),
                                                                          selected = list(""),inline = T))),
                                     prettyRadioButtons("busco_inter",shape = "round",animation="pulse",
                                                        status = "primary",bigger=T,
                                                        fill = TRUE,h3("Busco dataset"),
                                                        choices = list("Eukaryota"="eukaryota","Metazoa"="metazoa","Embryophyta"="embryophyta",
                                                                       "None"="None"),selected = "None",inline = T,),
                                     pickerInput(inputId = "clades_inter",label = h3("Select/deselect clades"), 
                                                 choices = levels(data_by_species$clade.qual),selected = levels(data_by_species$clade.qual), 
                                                 choicesOpt = list(
                                                   style = rep_len("font-size: 21px", length(levels(data_by_species$clade.qual)))
                                                 ),
                                                 options = list(style = "my-class",`actions-box` = TRUE,`selected-text-format` = "count > 3"
                                                 ), 
                                                 multiple = TRUE ),
                                     materialSwitch(inputId = "pgls_inter", label = h3("PGLS"),status="primary",inline = T),
                                     selectInput("tree_inter", h3("Tree for PGLS"), choices = phylogenetic_trees, selected = "1"),
                                     dropdown(
                                       tags$h2("Parameters"),
                                       pickerInput("svr_class",h3("Introns AS classes"), choices = c( "All major introns" = "all",
                                                                                                       "AS rate >= 5%" = "high_SV",
                                                                                                       "AS rate < 5%" = "low_SV"
                                       )),
                                       
                                       sliderInput("coverage_inter",h3("Minimal coverage (reads/bp)"),min = 0, max = 1000, value = 200),
                                       fileInput("upload_data", h3("Choose data File"),
                                                 accept = c(
                                                   "text/csv/tab",
                                                   "text/comma-separated-values,text/plain",
                                                   ".tab")),
                                       
                                       style = "unite", icon = icon("gears"),
                                       status = "primary", 
                                     )
                              ),
                              column(1,offset = 3,h2("X axis :")),
                              column(5,selectInput("x_inter", "",choices = axisInter_list_quantitative))
                       )
                       
              ),
              tabPanel("Inter-species Axis",dataTableOutput('tableInter')),
              ## INTRA SPECIES GRAPHIC
              tabPanel("Intra-species graphics",
                       column(12,
                              column(2,offset = 1,selectInput("species_selected_intra", h2("Species studied"),choices =  listNomSpecies , selected = "Drosophila melanogaster")),
                              column(2, style='padding:10px;',imageOutput('species_image_intra',height=5 , width=5), div(style = "height:120px;")),
                              column(3,selectInput("y_intra", h2("Y axis"),choices = axisIntra_list , selected = 1)),
                              column(3,selectInput("x_intra", h2("X axis or Histogram"),choices = axisIntra_list , selected = 1)),
                              column(9,class = "well",offset = 1,withSpinner(type = 4,color ="#136dc0", plotlyOutput("plot_intra",height = "100%",width="100%"))),
                              column(width =2,class="well",  
                                     
                                     fluidRow(column(4,
                                                     materialSwitch(inputId = "histogram_intra", label = h3("Histogram"),status="primary")),
                                              column(5,
                                                     prettyCheckboxGroup( "scale_intra", shape = "round",animation="pulse",bigger=T,
                                                                          status = "primary",
                                                                          fill = TRUE,h3("Scales"),choices = list("X log10","Y log10"),
                                                                          selected = list(""),inline = T))),
                                     prettyRadioButtons("busco_intra",shape = "round",animation="pulse",
                                                        status = "primary",bigger=T,
                                                        fill = TRUE,h3("Busco dataset"),
                                                        choices = list("Eukaryota"="eukaryota","Metazoa"="metazoa","Embryophyta"="embryophyta",
                                                                       "None"="None"),selected = "None",inline = T,),
                                     dropdown( 
                                       sliderInput("svr_range_intra",h3("AS range of the introns studied"),min = 0, max = 0.5, value =  c(0,0.5)),
                                       sliderInput("bin_intra",h3("Proportion of N by points (%)"),min = 0, max = 100, value =  10),style = "unite", icon = icon("gears"),
                                       status = "primary", 
                                     ),
                                     downloadButton("download_fpkm", "Genes Expression", style = "color: #fff; background-color: #27ae60; border-color: #fff;font-size: 21px;"),
                                     downloadButton("download_busco_id", "Busco identification", style = "color: #fff; background-color: #27ae60; border-color: #fff;font-size: 21px;"),
                                     downloadButton("download_svr", "Introns Alternative Splicing", style = "color: #fff; background-color: #27ae60; border-color: #fff;font-size: 21px;")
                              )
                       )
              ),
              tabPanel("Intra-species Axis",dataTableOutput('tableIntra')),
              
              ### GENE STRUCTURE
              tabPanel("Gene structure",
                       column(10,offset=1,class = "well",fluidRow(
                         column(2,offset = 1,
                                selectInput("species_gene_struct", h4("Species studied"), 
                                            choices = listNomSpecies , selected = "Drosophila melanogaster")),
                         column(2, imageOutput('species_image',height=10,width=10),
                                div(style = "height:100px;")),
                         
                         column(2,
                                prettyRadioButtons("gene_list",shape = "round",animation="pulse",
                                                   status = "primary",bigger=T,
                                                   fill = TRUE,h4("ID choice"),
                                                   choices = c("gene id"="gene_id","metazoa busco id" = "busco_id_metazoa",
                                                               "eukaryota busco id" = "busco_id_eukaryota",
                                                               "embryophyta busco id" = "busco_id_embryophyta"),selected = "gene_id",inline = T)),
                         
                         column(2, selectInput("studied_gene", h4("Gene studied"),choices="")),
                         column(2,  dropdown(   
                           sliderInput("sliderscale", h4("Bp Bins on x axis"),
                                       min = 0, max = 100000, value = 1000),
                           style = "unite", icon = icon("gears"),
                           status = "primary", 
                         ))))
                       ,column(10,offset = 1,class = "well",withSpinner(type = 4,color ="#136dc0",
                                                                        plotlyOutput("structureGene",width = "100%", height = "100%")))),
              
              ### PHYLOGENETIC TREE
              tabPanel("Phylogenetic tree",
                       fluidRow(
                         column(2,offset = 1,selectInput("select_tree", h3("Tree selected"),choices = phylogenetic_trees,
                                                         selected = phylogenetic_trees[1])),
                         column(2,selectInput("layout_tree", h3("Layout"),
                                              choices = c("roundrect","ellipse","circular","equal_angle","daylight")
                         )),
                         column(2,sliderInput("spacing", h3("Spacing"),min = 0, max = 1, value = 0.05
                         )),
                         column(2,sliderInput("tip_size", h3("Tips size"),min = 0, max = 10, value = 2
                         ))),
                       fluidRow(
                         column(10,offset = 1, class = "well",
                                column(6,
                                       withSpinner(type = 4,color ="#136dc0",plotOutput("plot_principal", height = 900,
                                                                                        brush = brushOpts(
                                                                                          id = "plot2_brush",
                                                                                          resetOnNew = TRUE
                                                                                        )
                                       )))
                                ,
                                column(6,h4("Zoom"),
                                       # dropdown(
                                       
                                       withSpinner(type = 4,color ="#136dc0",plotOutput("plot_zoomed", height = 600)),
                                       column(3,offset = 3,sliderInput("tip_size_zoom", h4("Tips size zoomed"),min = 0, max = 10, value = 2
                                       )),
                                       column(3,sliderInput("spacing_zoom", h4("Spacing zoomed"),min = 0, max = 1, value = 0.05
                                       ))
                                       # ,                                  tags$h3("Parameters"),
                                       
                                       # style = "unite", icon = icon("zoom-in", lib = "glyphicon"),
                                       # status = "primary")
                                )
                         )
                       )
              )
  )
)
)

