# Prepare for applications
<<<<<<< HEAD
# system("wget -O www/database_ShyniApp.tar.gz https://zenodo.org/records/10529523/files/database_ShyniApp.tar.gz?download=1")
# system("tar -xvzf www/database_ShyniApp.tar.gz")
# system("mv database www/database")
=======
system("wget -O www/database_ShyniApp.tar.gz https://zenodo.org/records/10527331/files/database_ShyniApp.tar.gz?download=1")
system("cd www")
system("tar -xvzf database_ShyniApp.tar.gz")
>>>>>>> a98af614b5a0019f6ef7e9e4d1bbd9e0d0996f4b


## Tree preparation
library(ape)
library(stringr)
list_files = list.files("database/dNdS/phylogeny/",recursive = T,pattern = "root.nwk" ,full.names = T)

phylo=data.frame()
for (file in list_files){
  print(file)
  tree = read.tree(file)
  set = str_split_1(str_replace(file,"database/dNdS/phylogeny//",""),"_")[1]
  nb_genes = ""
  if (file.exists(paste("data/dnds_phylo/",set,"/readme",sep=""))){
    readme = read.table(paste("data/dnds_phylo/",set,"/readme",sep=""))
    rownames(readme) = readme$V1
    nb_genes = readme["nb_genes_raxml","V2"]
  }
  nb_species = length(tree$tip.label)
  
  phylo = rbind(phylo,
                data.frame(
                  file,
                  name = paste("database/dNdS/phylogeny/",set,"_root.nwk",sep=""),
                  description = paste(set,nb_species,"species",nb_genes,"genes")
                ))
}

write.table(phylo[c(nrow(phylo):1),c("name" , "description")] , paste("phylogenetic_trees_description.tab",sep=""),quote=F,row.names = F,col.names = T,sep="\t")



## Informations table
options(stringsAsFactors = F, scipen = 999)
library(rgbif)
library(ggpubr)
library(ggrepel)
library(readxl)
library(seqinr)
library(stringr)
library(dplyr)
library(tidyr)

add_charac <- function(data_summary,label,description,value){
  data_summary = rbind(data_summary,data.frame(
    label = label,
    description = description,
    value=value
  ))
  return(data_summary)
}

get_rsquared_slope = function(prop.quantile = 0.1,Xaxis,Yaxis){
  quantile = unique(quantile(Xaxis, probs = seq(0, 1,prop.quantile),na.rm=T))
  intervalle = cut(Xaxis, quantile,include.lowest = T,include.higher=T)
  X = tapply(Xaxis, intervalle, median)
  if ( !any(is.na(X)) ){
    Y = tapply(Yaxis, intervalle, mean)
    X = log10(X)
    pearson_method = cor.test(X, Y,method="pearson")
    cor = pearson_method$estimate
    pval_cor = pearson_method$p.value
    return( c(cor,pval_cor) )
  } else {
    return(
      c(NA,NA)
    )
  }
}


list_species = read.delim(paste("www/database/list_species.tab",sep=""))
rownames(list_species) = list_species$species
list_species$path_db = paste(list_species$species,"_NCBI.taxid",list_species$NCBI.taxid,"/",list_species$assembly_accession,sep="")

data_lht = read.delim("www/database/life_history_traits.tab")
rownames(data_lht) = paste(data_lht$species,data_lht$life_history_traits)

list_species$max_lifespan_days = data_lht[paste(list_species$species,"lifespan_days"),]$value
list_species$max_length_cm = data_lht[paste(list_species$species,"length_cm"),]$value
list_species$max_weight_kg = data_lht[paste(list_species$species,"weight_kg"),]$value

all_dt = data.frame()
for (species in rev(list_species$species) ){print(species)
  # pathData = "/home/fbenitiere/data/Projet-SplicedVariants/"
  pathData = "/beegfs/data/fbenitiere/Projet-SplicedVariants/"
  
  gff_path = paste(pathData , "Annotations/",species,"/data_source/annotation.gff",sep="")
  gc_table_path = paste(pathData, "Annotations/",species,"/GC_content.tab",sep="")
  by_gene_analysis_path = paste("www/database/Transcriptomic/",list_species[species,]$path_db,"/by_gene_analysis.tab.gz",sep="")
  by_intron_analysis_path = paste("www/database/Transcriptomic/",list_species[species,]$path_db,"/by_intron_analysis.tab.gz",sep="")
  prot_path = paste(pathData , "Annotations/",species,"/data_source/protein.faa",sep="")
  
  data_summary = data.frame()
  data_summary = add_charac(data_summary,'species',"",species)
  
  con <- file(gff_path,"r")
  first_line <- readLines(con,n=10)
  close(con)
  print(first_line[5])
  genome_assembly = first_line[5]
  genome_assembly = str_replace(genome_assembly,"#!genome-build-accession NCBI_Assembly:","")
  data_summary = add_charac(data_summary,'genome_assembly',"",genome_assembly)
  
  data_summary = add_charac(data_summary,'clade;qual',"", list_species[species,]$clade_group)
  data_summary = add_charac(data_summary,'lifespan_days;quant',"", list_species[species,]$max_lifespan_days)
  data_summary = add_charac(data_summary,'length_cm;quant',"", list_species[species,]$max_length_cm)
  data_summary = add_charac(data_summary,'weight_kg;quant',"", list_species[species,]$max_weight_kg)
  
  key = name_backbone(name=str_replace(species,"_"," "),rank="species")$usageKey
  RGBIF_count = occ_search(key,limit=0)$meta$count
  
  data_summary = add_charac(data_summary,'RGBIF_observation;quant',"",RGBIF_count)
  
  
  genome_character = read.table(gc_table_path,header=T)
  genome_size_var =  sum(genome_character[genome_character$genome_character %in% c("A","T","G","C"), ]$Freq) / 1000000
  GC_content_proportion_var = round( sum(genome_character[genome_character$genome_character %in% c("G","C"), ]$Freq) /
                                       sum(genome_character[genome_character$genome_character %in% c("A","T","G","C"), ]$Freq) ,2)
  
  data_summary = add_charac(data_summary,'genome_size;quant',"Mb",genome_size_var)
  data_summary = add_charac(data_summary,'GC_content_proportion;quant',"",GC_content_proportion_var)
  data_summary = add_charac(data_summary,'No_prot_annot;quant',"",length(read.fasta(prot_path)))
  
  #### Data_set
  print(file.exists(by_gene_analysis_path))
  if ( file.exists(by_gene_analysis_path) ){
    by_gene =  read.delim(by_gene_analysis_path , header=T , sep="\t",comment.char = "#")
    by_gene = by_gene[by_gene$type == "gene" & grepl("gene_biotype=protein_coding" , by_gene$attributes),]
    rownames(by_gene) = by_gene$gene_id
    by_gene = by_gene[by_gene$type == "gene",] # FILTRE PSEUDOGENE
    
    by_intron =  read.delim(by_intron_analysis_path , header=T , sep="\t",comment.char = "#")
    by_intron = by_intron[by_intron$gene_id %in% by_gene$gene_id,] # FILTRE PSEUDOGENE
    by_intron$median_fpkm = by_gene[by_intron$gene_id,]$median_fpkm
    
    for (busco_group in c("metazoa","embryophyta","eukaryota","None")){ 
      can_analyse = T
      if ( busco_group != "None" ){
        if (file.exists(paste("www/database/BUSCO_annotations/",list_species[species,]$path_db,"/busco_to_gene_id_",busco_group,".gz",sep=""))){
          busco_to_gene = read.delim(paste("www/database/BUSCO_annotations/",list_species[species,]$path_db,"/busco_to_gene_id_",busco_group,".gz",sep=""))
          
          by_gene$busco_id = by_gene$gene_id %in% busco_to_gene$gene_id
          by_gene_selected = by_gene[by_gene$busco_id,]
          
          data_summary = add_charac(data_summary,paste("median_coverage_exon;buscodataset_",busco_group,";quant",sep=""),"",median(by_gene_selected$exon_coverage,na.rm = T))
          
          busco_to_gene = busco_to_gene[!(duplicated(busco_to_gene$busco_id,fromLast = FALSE) | duplicated(busco_to_gene$busco_id,fromLast = TRUE)) &
                                          !(duplicated(busco_to_gene$gene_id,fromLast = FALSE) | duplicated(busco_to_gene$gene_id,fromLast = TRUE)) ,]
          
          rownames(busco_to_gene) = busco_to_gene$gene_id
          by_gene$busco_id = by_gene$gene_id %in% busco_to_gene$gene_id
          by_intron$busco_id = by_intron$gene_id %in% busco_to_gene$gene_id
          
          by_gene_selected = by_gene[by_gene$busco_id,]
          by_intron_selected = by_intron[by_intron$busco_id & by_intron$intron_class == "major" & by_intron$into_cds == "True",]
          
        } else { can_analyse = F }
      } else {
        by_gene_selected = by_gene
        by_intron_selected = by_intron[ by_intron$intron_class == "major" & by_intron$into_cds == "True",]
        data_summary = add_charac(data_summary,paste("median_coverage_exon;buscodataset_",busco_group,";quant",sep=""),"",median(by_gene_selected$exon_coverage,na.rm = T))
      }
      
      if (can_analyse){
        data_summary = add_charac(data_summary,paste("No_genes;buscodataset_",busco_group,";quant",sep=""),"",nrow(by_gene_selected))
        data_summary = add_charac(data_summary,paste("No_introns;buscodataset_",busco_group,";quant",sep=""),"",nrow(by_intron_selected))
        data_summary = add_charac(data_summary,paste("median_length_introns;buscodataset_",busco_group,";quant",sep=""),"",median(abs(by_intron_selected$splice3-by_intron_selected$splice5)))
        data_summary = add_charac(data_summary,paste("No_multi_exonic_genes;buscodataset_",busco_group,";quant",sep=""),"",length(unique(by_intron_selected$gene_id)))
        data_summary = add_charac(data_summary,paste("No_intron_per_gene;buscodataset_",busco_group,";quant",sep=""),"",nrow(by_intron_selected) / length(unique(by_intron_selected$gene_id)))
        data_summary = add_charac(data_summary,paste("prop_intron_annotated;buscodataset_",busco_group,";quant",sep=""),"",sum(by_intron_selected$Annotation) / nrow(by_intron_selected))
        
        
        data_summary = add_charac(data_summary,paste("median_gene_fpkm;buscodataset_",busco_group,";quant",sep=""),"",median(by_gene_selected$median_fpkm,na.rm = T))
        
        for (svr_class in c("all" , "high_SV" , "low_SV")){
          by_intron_selected_svr = by_intron_selected
          if ( svr_class == "high_SV" ){
            by_intron_selected_svr = by_intron_selected[  by_intron_selected$splice_variant_rate >= 0.05 ,]
          } else if (svr_class == "low_SV") {
            by_intron_selected_svr = by_intron_selected[  by_intron_selected$splice_variant_rate < 0.05 ,]
          }
          
          data_summary = add_charac(data_summary,paste("average_svr;svr_class_",svr_class,";buscodataset_",busco_group,";quant",sep=""),"",mean(by_intron_selected_svr$splice_variant_rate))
        }
      }
    }
  }
  data_summary$species = species
  
  all_dt = rbind(all_dt,data_summary)
}

columns = unique(all_dt$label)


species_table = data.frame()
for (species in unique(all_dt$species)){print(species)
  dt = all_dt[all_dt$species==species,]
  rownames(dt) = dt$label
  species_table = rbind(species_table,dt[columns,]$value)
}

colnames(species_table) = columns


list_method = list.files("www/database/dNdS",full.names = F,recursive = F,pattern = ".tab")
list_method = str_replace(list_method,".tab","")

for (method in list_method){ print(method)
  data_dNdS = read.delim(paste("www/database/dNdS/",method,".tab",sep=""))
  rownames(data_dNdS) = data_dNdS$species
  
  species_table[,paste(method,"_dn;quant",sep="")] = data_dNdS[species_table$species,]$dN
  species_table[,paste(method,"_ds;quant",sep="")] = data_dNdS[species_table$species,]$dS
  species_table[,paste(method,"_dnds;quant",sep="")] = data_dNdS[species_table$species,]$dNdS
}

write.table(species_table,paste("www/species_informations_tables/data_by_species.tab",sep=""),sep="\t",quote=F,row.names = F, col.names = T)

# inter_axis = data.frame( name_label = colnames(species_table))
# inter_axis[,"group"] = "Default"
# group	display_label	name_label	description	quantitative
# write.table(all_dt,paste("www/species_informations_tables/data_by_species.tab",sep=""),sep="\t",quote=F,row.names = F, col.names = T)

