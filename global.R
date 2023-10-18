options(stringsAsFactors = F, scipen = 999) 
library(ape)
library(stringr)
library(plotly)
library(ggtree)
library(shinythemes)
library(caper)
library(shinyWidgets)
library(shinyjs)
library(phylolm)
library(shinycssloaders)
library(RColorBrewer)
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))


std <- function(x) sd(x)/sqrt(length(x))

lm_eqn <- function(m=lm(Y ~ X,data)){
  paste( "R2 = " , round( summary(m)$r.squared , 2) , "; p-value = " , formatC(summary(m)$coefficients[2,4], format = "e" , digits = 0 ))
}


GLS <- function(dataframe=shorebird){
  # aic = 1000000
  # dt = data.frame()
  # for (model in c("LM","lambda","OUfixedRoot","OUrandomRoot","BM")){
  #   for (measurement_error in c(T,F)){
  #     if (model == "LM"){
  #       fit = lm(pgls_y~pgls_x, data = dataframe$data)
  #       measurement_error = NA
  #     } else if (model != "lambda"){
  #       fit <- phylolm(pgls_y~pgls_x, phy = dataframe$phy, data = dataframe$data, model = model,measurement_error=measurement_error)
  #     } else{ fit <- phylolm(pgls_y~pgls_x, phy = dataframe$phy, data = dataframe$data, model = model)
  #     measurement_error = NA}
  #     a = summary(fit)
  #     if (length(a$optpar)==0){a$optpar=NA}
  #     if (length(a$aic)==0){a$aic=NA
  #     a$logLik=NA
  #     a$optpar=NA
  #     a$sigma2=NA}
  #     
  #     dt = rbind(dt,data.frame(
  #       model,
  #       measurement_error,
  #       p_val_slope = a$coefficients[2,4],
  #       r.squared = a$r.squared,
  #       adj.r.squared = a$adj.r.squared,
  #       aic = a$aic,
  #       logLik = a$logLik,
  #       optpar = a$optpar,
  #       sigma2 = a$sigma2
  #     ))
  #     if ( !is.na(a$aic < aic) & a$aic < aic ){ best_fit_model = fit
  #     best_model = model
  #     aic = a$aic}
  #   }
  # }
  # dt = dt[!duplicated(dt$aic),]
  # dt = dt[order(dt$aic),]
  # return(list(dt,best_fit_model,best_model))
  return(list(0,0,0))
}


Clade_color = c("Other Invertebrates"="#f5b48a","Lepido Diptera"="red","Other Vertebrates"="#A6CEE3","Other Insecta"="#FF7F00",
                Nematoda="#B2DF8A",Teleostei="#1F78B4",Hymenoptera="#ba8e18",Aves="#5b5b5b",Mammalia="#66281A",Embryophyta="#33A02C"
)

table_phylo = read.delim("www/phylogenetic_trees_description.tab")
phylogenetic_trees = paste("www/",table_phylo$name,sep="")
names(phylogenetic_trees) = table_phylo$description

data_by_species = data.frame(species="")
for (file in list.files("www/species_informations_tables",full.names = T,pattern = "data_by_species.tab")){
  dt = read.delim(file,header = T)
  data_by_species = merge(dt,data_by_species, by.x = "species", by.y = "species", all.x = TRUE, all.y = TRUE)
}
data_by_species$clade.qual = factor(data_by_species$clade.qual, levels = c("Embryophyta","Lepido Diptera","Hymenoptera",
                                                                           "Other Insecta","Nematoda","Other Invertebrates",
                                                                           "Mammalia","Aves","Teleostei","Other Vertebrates"))

dt_species = read.delim("www/database/list_species.tab",header=T)
rownames(dt_species) = dt_species$species
all_listNomSpecies = tapply(dt_species$species,dt_species$clade_group,function(x)  str_replace_all(x,"_"," "))
dt_species = dt_species[dt_species$expression_data,]
dt_species$path_db = paste(dt_species$species,"_NCBI.taxid",dt_species$NCBI.taxid,"/",dt_species$assembly_accession,sep="")
listNomSpecies = tapply(dt_species$species,dt_species$clade_group,function(x)  str_replace_all(x,"_"," "))

axisInter = read.delim("www/inter_axis.tab",sep="\t")
ExplicationsInter = axisInter[axisInter$group != "",c("display_label","description")]

axisInter_quantitative = axisInter[axisInter$quantitative,]
axisInter_qualitative = axisInter[!axisInter$quantitative,]


axisInter_list_qualitative = tapply(axisInter_qualitative$display_label,axisInter_qualitative$group,list)
axisInter_list_quantitative = tapply(axisInter_quantitative$display_label,axisInter_quantitative$group,list)



axisIntra = read.delim("www/intra_axis.tab",sep="\t")
axisIntra_list = tapply(axisIntra$display_label,axisIntra$group,list)
ExplicationsIntra = axisIntra[axisIntra$group != "",c("display_label","description")]

