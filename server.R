
server <- function(input, output,session) {
  
  output$tableIntra <- renderDataTable(ExplicationsIntra) #tableau des explications des axes du 2eme onglet
  
  output$tableInter <- renderDataTable(ExplicationsInter) #tableau des explications des axes du 5eme onglet
  
  
  # Block mount image
  addResourcePath(prefix = "imgResources", directoryPath = "www/species_images/")
  
  hover_event <- reactive({
    event_data(event = "plotly_hover", source = "hoverplotsource")
  })
  
  unhover_event <- reactive({
    event_data(event = "plotly_unhover", source = "hoverplotsource")
  })
  
  hoverplotlyProxy <- plotlyProxy("plot_inter", session)
  
  observeEvent(unhover_event(), {
    hoverplotlyProxy %>%
      plotlyProxyInvoke("relayout", list(images = list(NULL)))
  })
  
  observeEvent(input$upload_data,ignoreNULL = T, {
    inFile <- input$upload_data
    dt = read.delim(inFile$datapath,header = T)
    data_by_species <<- merge(dt,data_by_species, by.x = "species", by.y = "species", all.x = TRUE, all.y = TRUE)
    
    axisInter <<- rbind(axisInter,data.frame(group = "Uploaded" , "display_label" = colnames(dt),name_label=colnames(dt),description=NA,quantitative=T))
    axisInter_quantitative <<- axisInter[axisInter$quantitative,]
    axisInter_list_quantitative <<- tapply(axisInter_quantitative$display_label,axisInter_quantitative$group,list)
    updateSelectizeInput(session, "y_inter", choices = axisInter_list_quantitative,options = list(maxOptions = 3000))
    updateSelectizeInput(session, "x_inter", choices = axisInter_list_quantitative,options = list(maxOptions = 3000))
  })
  
  observeEvent(hover_event(), {
    dt = unlist(str_split(hover_event()$customdata,"_;_"))
    if ("Y log10" %in% input$scale_inter  & "X log10" %in% input$scale_inter){
      hoverplotlyProxy %>%
        plotlyProxyInvoke("relayout", list(images = list(
          list(
            source = dt[5],
            xref = "x",
            yref = "y",
            x = hover_event()$x,
            y = hover_event()$y,
            sizex = (log10(as.numeric(dt[2])) - log10(as.numeric(dt[1])))/6,
            sizey = (log10(as.numeric(dt[4])) - log10(as.numeric(dt[3])))/6,
            opacity = 1
          )
        )))
    } else {
      hoverplotlyProxy %>%
        plotlyProxyInvoke("relayout", list(images = list(
          list(
            source = dt[5],
            xref = "x",
            yref = "y",
            x = hover_event()$x,
            y = hover_event()$y,
            sizex = (as.numeric(dt[2]) - as.numeric(dt[1]))/6,
            sizey = (as.numeric(dt[4]) - as.numeric(dt[3]))/6,
            opacity = 1
          )
        )))
    }
  })
  
  ###### End block
  observeEvent(input$species_gene_struct,ignoreNULL = FALSE,{
    species = input$species_gene_struct
    species = str_replace_all(species," ","_")
    if (dt_species[species,]$clade_group == "Embryophyta" ){
      updatePrettyRadioButtons(session,
                               prettyOptions = list(shape = "round",animation="pulse",
                                                    status = "primary",bigger=T,
                                                    fill = TRUE), "gene_list",
                               choices = c("gene id"="gene_id",
                                           "eukaryota busco id" = "busco_id_eukaryota",
                                           "embryophyta busco id" = "busco_id_embryophyta"),selected = "gene_id",inline = T)
    } else {
      updatePrettyRadioButtons(session,
                               prettyOptions = list(shape = "round",animation="pulse",
                                                    status = "primary",bigger=T,
                                                    fill = TRUE), "gene_list",
                               choices = c("gene id"="gene_id","metazoa busco id" = "busco_id_metazoa",
                                           "eukaryota busco id" = "busco_id_eukaryota"
                               ),selected = "gene_id")
    }
  })
  
  observeEvent(input$boxplot_inter,ignoreNULL = FALSE,{
    if (input$boxplot_inter){
      updateSelectizeInput(session, "x_inter", choices = axisInter_list_qualitative,options=list(maxOptions=3000))
      updatePrettyCheckboxGroup(session, 
                                prettyOptions = list(shape = "round",animation="pulse",bigger=T,
                                                     status = "primary",
                                                     fill = TRUE), "scale_inter", choices = list("Y log10"))
      updateMaterialSwitch(session, "pgls_inter", value = F)
      shinyjs::disable("pgls_inter")
    } else {    
      updateSelectizeInput(session, "x_inter", choices = axisInter_list_quantitative,options=list(maxOptions=3000))
      updatePrettyCheckboxGroup(session, 
                                prettyOptions = list(shape = "round",animation="pulse",bigger=T,
                                                     status = "primary",
                                                     fill = TRUE),"scale_inter", choices = list("X log10","Y log10"))
      shinyjs::enable("pgls_inter")
    }
  })
  
  
  observeEvent(input$pgls_inter, {
    if (!input$pgls_inter){
      shinyjs::disable("tree_inter")
    } else {
      shinyjs::enable("tree_inter")}
  })
  
  
  output$plot_inter <- renderPlotly({
    ylabel = axisInter[axisInter$display_label == input$y_inter,]$name_label 
    xlabel = axisInter[axisInter$display_label == input$x_inter,]$name_label
    
    if ( grepl("buscodataset_",xlabel)){ xlabel = paste(xlabel,input$busco_inter,".quant",sep="")}
    if ( grepl("buscodataset_",ylabel)){ ylabel = paste(ylabel,input$busco_inter,".quant",sep="")}
    
    if ( grepl("svr_class_",xlabel)){ xlabel = str_replace(xlabel,"all",input$svr_class)}
    if ( grepl("svr_class_",ylabel)){ ylabel = str_replace(ylabel,"all",input$svr_class)}
    if ( input$busco_inter != "None"){ 
      minimum_coverage = paste("median_coverage_exon.buscodataset_",input$busco_inter,".quant",sep="")
      data_by_species = data_by_species[data_by_species[minimum_coverage] > input$coverage_inter , ]
    }
    
    if ( !all(input$clades_inter == levels(data_by_species$clade.qual))){
      data_by_species = data_by_species[data_by_species$clade.qual %in% input$clades_inter,]
    }
    data_by_species = data_by_species[!is.na(data_by_species[,xlabel]) & !is.na(data_by_species[,ylabel]),]
    if ( input$pgls_inter ){
      arbrePhylo = read.tree(input$tree_inter)
      data_by_species = data_by_species[data_by_species$species %in% arbrePhylo$tip.label,]
    }
    data_by_species$img_local <- paste("imgResources/",data_by_species$species,".png",sep="")
    
    if ( !input$boxplot_inter  ){
      data_by_species$custom_data = paste(min(data_by_species[,xlabel]),max(data_by_species[,xlabel]),
                                          min(data_by_species[,ylabel]),max(data_by_species[,ylabel]),data_by_species$img_local,sep="_;_")
    } else {
      data_by_species$custom_data = "custom_data"
    }
    data_by_species$speciesname = str_replace(data_by_species$species,"_"," ")
    p = ggplot(
      data_by_species,aes_string(x=xlabel,y=ylabel , text="speciesname",
                                 customdata="custom_data"))+ 
      scale_fill_manual(values=Clade_color) +
      xlab(input$x_inter) + 
      ylab(input$y_inter) + theme_bw() + theme(
        axis.title.x = element_text(color="black", size=25,family="economica"),
        axis.title.y = element_text(color="black", size=25, family="economica"),
        axis.text.y =  element_text(color="black", size=20, family="economica"),
        axis.text.x =  element_text(color="black", size=20, family="economica"),
        title =  element_text(color="black", size=15, family="economica"),
        legend.text =  element_text(color="black", size=20, family="economica")
      ) + geom_point(aes(fill=clade.qual),size=4,alpha=0.7)
    p
    
    if ( "Y log10" %in% input$scale_inter  ){
      lm_y = log10(data_by_species[,ylabel])
      p = p + scale_y_log10()
    } else {
      lm_y = data_by_species[,ylabel]
    }
    
    if ( "X log10" %in% input$scale_inter  ){ 
      lm_x = log10(data_by_species[,xlabel])
      p = p + scale_x_log10()
    } else {
      lm_x = data_by_species[,xlabel]
    }
    
    if ( input$boxplot_inter  ){
      p = p + geom_boxplot(alpha=.1,fill="grey") + 
        geom_point(aes(fill=clade.qual),size=3,alpha=0.7)+ theme_bw() + theme(
          axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
          axis.title.y = element_text(color="black", size=25, family="economica"),
          axis.text.y =  element_text(color="black", size=20, family="economica"),
          axis.text.x =  element_text(color="black", size=25,angle = 50, family="economica"),
          title =  element_text(color="black", size=15, family="economica"),
          legend.text =  element_text(color="black", size=20, family="economica")
        ) + ggtitle(paste("N=",nrow(data_by_species))) + theme(legend.position='none')
    } else if ( !input$boxplot_inter ){
      if (input$pgls_inter){
        shorebird <- comparative.data(arbrePhylo, 
                                      data.frame(species=data_by_species$species,
                                                 pgls_x=lm_x,
                                                 pgls_y=lm_y), species, vcv=TRUE)
        
        gls = GLS(shorebird)
        p = p + geom_abline(slope=1,intercept=0,alpha = .6) +
          ggtitle(paste("N=",nrow(data_by_species)," / LM:",lm_eqn(lm(lm_y ~ lm_x)),
                        " / PGLS:",lm_eqn(pgls(pgls_y~pgls_x,shorebird))
                        # "/ Best",gls[[3]],":",lm_eqn(gls[[2]])
          ))
      } else { 
        p = p + geom_abline(slope=1,intercept=0,alpha = .6)+ ggtitle(paste("N=",nrow(data_by_species)," / LM:",lm_eqn(lm(lm_y~lm_x))))
      }
    }
    
    ggplotly( p, tooltip = c("text") , 
              height = 1000*.8, 
              hoverinfo = 'none',
              source = "hoverplotsource"
    ) %>%
      event_register('plotly_hover') %>%
      event_register('plotly_unhover')
    
  })
  
  
  observe({
    if (input$tabs == "Intra-species graphics"){
      species = input$species_selected_intra
      
      species = str_replace_all(species," ","_")
      species_genes = read.delim(paste("www/database/Transcriptomic/",dt_species[species,]$path_db,"/by_gene_analysis.tab.gz",sep="") , header=T , sep="\t",comment.char = "#")
      rownames(species_genes) = species_genes$gene_id
      species_intron = read.delim(paste("www/database/Transcriptomic/",dt_species[species,]$path_db,"/by_intron_analysis.tab.gz",sep="") , header=T , sep="\t",comment.char = "#")
      species_intron$median_fpkm = species_genes[species_intron$gene_id,]$median_fpkm
      species_intron$svr = species_intron$splice_variant_rate
      species_intron$nsvr = species_intron$nonsplice_variant_rate
      species_intron$sum_ns = species_intron$ns
      species_intron$sum_na = species_intron$na_spl3 + species_intron$na_spl5
      species_intron$sum_nu = species_intron$nu_spl3 + species_intron$nu_spl5
      species_intron$geneprop_tabid = species_intron$gene_id
      species_intron$length = abs(species_intron$splice5 - species_intron$splice3)
      species_intron <<- species_intron
    }
  })
  
  
  output$species_image_intra <- renderImage({
    secies_name <- str_replace(input$species_selected_intra," ","_")
    list(src=paste("www/species_images/",
                   secies_name,".png",sep=""))
  },deleteFile=FALSE)
  
  
  output$plot_intra <- renderPlotly({
    species = str_replace(input$species_selected_intra," ","_")
    color_species = Clade_color[dt_species[species,]$clade_group]
    
    species_intron = species_intron[ species_intron$into_cds == "True",]
    
    
    if ("None" != input$busco_intra){
      busco_gene = read.delim(paste("www/database/BUSCO_annotations/",dt_species[species,]$path_db,"/busco_to_gene_id_",input$busco_intra,sep=""))
      
      species_intron = species_intron[species_intron$gene_id %in% busco_gene$gene_id, ]
    }
    
    
    species_intron = species_intron[!is.na(species_intron$svr) & 
                                      species_intron$svr >= input$svr_range_intra[1] & 
                                      species_intron$svr < input$svr_range_intra[2],]
    
    
    gene_n1 =  tapply(species_intron$sum_ns,species_intron$geneprop_tabid,sum)
    gene_n2 =   tapply(species_intron$sum_na,species_intron$geneprop_tabid,sum) 
    SVRgene = 1 - (1 - (gene_n2 / (gene_n1 + gene_n2)))^(tapply(species_intron$sum_ns,species_intron$geneprop_tabid,length))
    averageLength = tapply(species_intron$length,species_intron$geneprop_tabid,mean) # average length intron
    FPKMgene = tapply(species_intron$median_fpkm,species_intron$geneprop_tabid,mean) # FPKM per gene
    
    listeAxis = list("FPKMgene"=FPKMgene,"SVRgene"=SVRgene,"IntronPerGene"=table(species_intron$geneprop_tabid),"averageLength"=averageLength,
                     "Ns"=species_intron$sum_ns,"Ns+Na"=species_intron$sum_ns+species_intron$sum_na,
                     "SVRintron"=species_intron$svr,"FPKMintron"=species_intron$median_fpkm,"IntronLength"=species_intron$length)
    
    yintra = axisIntra[axisIntra$display_label == input$y_intra,]$name_label 
    xintra = axisIntra[axisIntra$display_label == input$x_intra,]$name_label
    if ( !input$histogram_intra ){
      xaxis = unlist(listeAxis[xintra])
      proportion = input$bin_intra / 100
      quantile = unique(quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T))
      intervalle = cut(xaxis, quantile,include.lowest = T,include.higher=T)
      
      X = tapply(xaxis, intervalle, mean)
      XerrorBar = tapply(xaxis, intervalle, std)
      yaxis = unlist(listeAxis[yintra])
      Y = tapply(yaxis, intervalle, mean)
      YerrorBar = tapply(yaxis, intervalle, std)
      
      data_sp = data.frame(X=X ,Y=Y ,XerrorBar=XerrorBar ,YerrorBar=YerrorBar)
      table(intervalle)
      
      p9 = ggplot(data_sp,aes(x=X,y=Y,text=paste("Nb of samples by group",table(intervalle))))  + theme_bw() +
        ylab(input$y_intra) + xlab(axisIntra==input$x_intra) +
        geom_errorbar(aes(ymin=Y-YerrorBar, ymax=Y+YerrorBar),size=0.1) +
        geom_errorbarh(aes(xmin=X-XerrorBar, xmax=X+XerrorBar),size=0.1)+
        ggtitle(paste("No genes or introns studied=",sum(table(intervalle)))) + 
        geom_point(pch=21,size=5,fill=color_species) + theme(
          axis.title.x = element_text(color="black", size=25,family="economica"),
          axis.title.y = element_text(color="black", size=25, family="economica"),
          axis.text.y =  element_text(color="black", size=20, family="economica"),
          axis.text.x =  element_text(color="black", size=20, family="economica"),
          title =  element_text(color="black", size=15, family="economica"),
          legend.text =  element_text(color="black", size=20, family="economica")
        )
    } else {
      
      data_sp = data.frame(X=unlist(listeAxis[xintra]))
      
      p9 = ggplot(data_sp,aes(x=X)) + geom_histogram( position="dodge",bins=200,alpha=0.8,fill=color_species,col="black")+ theme_bw() +
        xlab(input$x_intra) + ggtitle(paste("N=",nrow(data_sp))) + ylab("Density") + theme(
          axis.title.x = element_text(color="black", size=25,family="economica"),
          axis.title.y = element_text(color="black", size=25, family="economica"),
          axis.text.y =  element_text(color="black", size=20, family="economica"),
          axis.text.x =  element_text(color="black", size=20, family="economica"),
          title =  element_text(color="black", size=15, family="economica"),
          legend.text =  element_text(color="black", size=20, family="economica")
        )
    }
    
    if ( "Y log10" %in% input$scale_intra  ){ 
      p9 = p9 + scale_y_log10()
    }
    
    if ( "X log10" %in% input$scale_intra  ){
      p9 = p9 + scale_x_log10()
    } 
    
    ggplotly(p9,height = 700, width = 1100) 
  })
  
  ### ONGLET STRUCTURE GENE
  observe({
    if (input$tabs == "Gene structure"){
      species = input$species_gene_struct
      species = str_replace_all(species," ","_")
      species_genes = read.delim(paste("www/database/Transcriptomic/",dt_species[species,]$path_db,"/by_gene_analysis.tab.gz",sep=""), header=T , sep="\t",comment.char = "#")
      
      if ( grepl("busco_id_",input$gene_list) ){
        domain = str_replace(input$gene_list,"busco_id_","")
        
        busco_gene = read.delim(paste("www/database/BUSCO_annotations/",dt_species[species,]$path_db,"/busco_to_gene_id_",domain,sep=""))
        busco_gene = busco_gene[!(duplicated(busco_gene$busco_id,fromLast = FALSE) | duplicated(busco_gene$busco_id,fromLast = TRUE)) &
                                  !(duplicated(busco_gene$gene_id,fromLast = FALSE) | duplicated(busco_gene$gene_id,fromLast = TRUE)) ,]
        
        rownames(busco_gene) = busco_gene$gene_id
        
        species_genes = species_genes[species_genes$gene_id %in% busco_gene$gene_id, ]
        species_genes$busco_id = busco_gene[species_genes$gene_id,]$busco_id
        nameGene <<- species_genes$gene_id
        names(nameGene) <<- paste(species_genes$busco_id,species_genes$gene_name,sep=" | ")
        
      } else {
        nameGene <<- species_genes$gene_id
        names(nameGene) <<- paste(species_genes$gene_id,species_genes$gene_name,sep=" | ")
      }
      updateSelectizeInput(session, "studied_gene", choices=nameGene,selected = nameGene[1],
                           options=list(maxOptions=100), server = T)
    }
  })
  
  output$species_image <- renderImage({
    secies_name <- str_replace(input$species_gene_struct," ","_")
    list(src=paste("www/species_images/",
                   secies_name,".png",sep=""))
  } , deleteFile=FALSE)
  
  output$structureGene <- renderPlotly({
    species = input$species_gene_struct
    species = str_replace_all(species," ","_")
    id_selected = input$studied_gene
    if (id_selected != ""){
      if (id_selected != ""){
        start.time <- Sys.time()
        header = readLines(paste("www/database/Transcriptomic/",dt_species[species,]$path_db,"/by_intron_analysis.tab.gz",sep=""),n=16)
        header = read.table(text=header)
        intron = system(paste("zgrep '",id_selected,"\t' ","www/database/Transcriptomic/",dt_species[species,]$path_db,"/by_intron_analysis.tab.gz",sep=""),intern=T)
        intron = read.table(text=intron,sep="\t")
        colnames(intron) = header
        intron$length = abs(intron$splice5 - intron$splice3)
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        as.data.frame(intron,sep="\t")
        
        intron$category_intron = paste("Class:",intron$intron_class," in CDS:",intron$into_cds,sep="" )
        
        intron$position = paste("Sp3:",intron$splice3,"Sp5:",intron$splice5)
        intron$sum_ns = as.numeric(intron$ns)
        intron = intron[,colnames(intron)[colnames(intron)!= "id"]]
        
        p2 = ggplot(intron,aes(group=category_intron)) +
          geom_rect( aes(label=position , xmin=splice3,ymin=sum_ns*1.1,xmax=splice5,ymax=sum_ns,fill=category_intron),
                     size=0.5,alpha=1,col="black" ) +
          geom_segment(aes(x=splice3,y=sum_ns*1.1,xend=splice5+length/2*(splice3-splice5)/abs(splice5-splice3),
                           yend=sum_ns*1.4,fill=category_intron),
                       size=0.5,alpha=1,col="black")+
          geom_segment(aes(x=splice5,y=sum_ns*1.1,xend=splice3-length/2*(splice3-splice5)/abs(splice5-splice3),
                           yend=sum_ns*1.4,fill=category_intron),
                       size=0.5,alpha=1,col="black") +
          theme_bw()  +
          ylab("Sum of Ns") +
          scale_fill_manual(name="Intron\ngroup",values=set_color) +
          xlab("Position on chromosome (bp)") +
          ggtitle(paste("Chromosome:",intron[1,"seqname"],"and Strand:",intron[1,"strand"])) +  theme(
            axis.title.x = element_text(color="black", size=25,family="economica"),
            axis.title.y = element_text(color="black", size=25, family="economica"),
            axis.text.y =  element_text(color="black", size=20, family="economica"),
            axis.text.x =  element_text(color="black", size=20, family="economica"),
            title =  element_text(color="black", size=15, family="economica"),
            text =  element_text(color="black", size=20, family="economica"),
            legend.text =  element_text(color="black", size=20, family="economica")
          ) + scale_y_log10()  +
          geom_vline(aes(xintercept=splice3,col=category_intron),alpha=0.1,fill="black") +
          geom_vline(aes(xintercept=splice5,col=category_intron),alpha=0.1,fill="black")+
          guides(col = guide_legend(override.aes = list( size = 6,alpha=1),
                                    label.theme = element_text(color="black",
                                                               size=26,face="italic", family="economica",vjust = 1.5,margin = margin(t = 5))))
        
        if (length(seq(min(intron$splice3,intron$splice5),
                       max(intron$splice3,intron$splice5),
                       input$sliderscale)) < 10){
          p2 = p2 + scale_x_continuous(breaks = seq(min(intron$splice3,intron$splice5),
                                                    max(intron$splice3,intron$splice5),
                                                    input$sliderscale) )
        }
        
        ggplotly(p2,height = 1000 * .8, )
      }
    }
  })
  
  ## PHYLOGENETIC TREE
  ranges2 <- reactiveValues(x = NULL, y = NULL, text_size=NULL)
  
  observe({
    brush <- input$plot2_brush
    if (!is.null(brush)) {
      ranges2$x <- c(brush$xmin, brush$xmax)
      ranges2$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })
  
  observe({
    tree_name <- input$select_tree
    tree <- read.tree(tree_name)
    # tree$tip.label[grepl("dingo",tree$tip.label)] = "Canis_lupus"
    tree$tip.label <- str_replace_all(tree$tip.label,"_"," ")
    edge_group <- str_replace_all(tree$tip.label,"_"," ")
    edge_clade <- rep("branch",length(tree$edge[,2]))
    for (group in unique(edge_group)){
      if (group %in% unlist(all_listNomSpecies)){
        edge_clade[tree$edge[,2] %in% grep(group,edge_group)] =
          names(all_listNomSpecies[unlist(lapply(all_listNomSpecies,function(x) group %in% x))])
      }
    }
    
    for (clade in names(all_listNomSpecies)){
      edge_clade[ which.edge(tree,  tree$edge[,2][edge_clade == clade] ) ] = clade
    }
    node_metadata = data.frame(node=tree$edge[,2],color=edge_clade)
    
    node_metadata$color = factor(node_metadata$color, levels = c("Embryophyta","Lepido Diptera","Hymenoptera",
                                                                 "Other Insecta","Nematoda","Other Invertebrates",
                                                                 "Mammalia","Aves","Teleostei","Other Vertebrates"))
    
    p = ggtree(tree, layout=input$layout_tree,size=1)  
    p <- p %<+% node_metadata  + aes(color=color) + 
      scale_color_manual("Clade",values=Clade_color[unique(edge_clade)]) +    theme(
        panel.background = element_rect(fill = "#f5f5f5", linetype = "dashed")
      ) + theme(
        title =  element_text(color="black", size=25, family="economica"),
        legend.text =  element_text(color="black", size=20, family="economica")
      )+ guides(color = guide_legend(override.aes = list(lwd = 3)))
    tree_plot <<- p
  })
  
  output$plot_principal <- renderPlot({
    input$layout_tree
    input$select_tree
    tree_plot + geom_tiplab(size=input$tip_size,nudge_x = input$spacing,colour="black")
  })
  
  output$plot_zoomed <- renderPlot({
    input$layout_tree
    input$select_tree
    tree_plot + geom_tiplab(size=input$tip_size_zoom,nudge_x = input$spacing_zoom,colour="black") +
      coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE) +  theme(legend.position="none")
  })
  
  
  output$download_fpkm <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), sep="")
    },
    content = function(file) {
      species = str_replace(input$species_selected_intra," ","_")
      data = read.delim(paste("www/database/Transcriptomic/",dt_species[species,]$path_db,"/by_gene_analysis.tab.gz",sep=""))
      write.table(data, file=file,row.names=F, col.names=T, sep="\t", quote=F)
    }
  )
  output$download_busco_id <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), sep="")
    },
    content = function(file) {
      species = str_replace(input$species_selected_intra," ","_")
      domain = input$busco_intra
      data = read.delim(paste("www/database/BUSCO_annotations/",dt_species[species,]$path_db,"/busco_to_gene_id_",domain,sep=""))
      write.table(data, file=file,row.names=F, col.names=T, sep="\t", quote=F)
    }
  )
  
  output$download_svr <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), sep="")
    },
    content = function(file) {
      species = str_replace(input$species_selected_intra," ","_")
      data = read.delim(paste("www/database/Transcriptomic/",dt_species[species,]$path_db,"/by_intron_analysis.tab.gz",sep=""))
      write.table(data, file=file,row.names=F, col.names=T, sep="\t", quote=F)
    }
  )
}
