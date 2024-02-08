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

## Source codesetSliderColor
setSliderColor <- function(color, sliderId) {
  
  # some tests to control inputs
  stopifnot(!is.null(color))
  stopifnot(is.character(color))
  stopifnot(is.numeric(sliderId))
  stopifnot(!is.null(sliderId))
  
  # the css class for ionrangeslider starts from 0
  # therefore need to remove 1 from sliderId
  sliderId <- sliderId - 1
  
  # create custom css background for each slider
  # selected by the user
  sliderCol <- lapply(sliderId, FUN = function(i) {
    paste0(
      ".js-irs-", i, " .irs-single,",
      " .js-irs-", i, " .irs-from,",
      " .js-irs-", i, " .irs-to,",
      " .js-irs-", i, " .irs-bar-edge,",
      " .js-irs-", i,
      " .irs-bar{  border-color: transparent;background: ", color[i+1],
      "; border-top: 1px solid ", color[i+1],
      "; border-bottom: 1px solid ", color[i+1],
      ";}"
    )
  })
  
  # insert this custom css code in the head
  # of the shiy app
  custom_head <- tags$head(tags$style(HTML(as.character(sliderCol))))
  return(custom_head)
}
##

fitted_model <- function(x_value=x_axis,y_value=y_axis,species_label=species_label,tree = NA){
  dt_fit = data.frame()
  if ( length(tree) != 1){ # If tree exist
    shorebird <- comparative.data(tree, data.frame(label=species_label, x=x_value, y=y_value), label, vcv=TRUE)
    fit = pgls(y~x,shorebird)
    summ_fit = summary(fit)
    dt_fit = rbind(dt_fit,data.frame(
      model="PGLS",
      p_val_slope = summ_fit$coefficients[2,4],
      r.squared = summ_fit$r.squared,
      adj.r.squared = summ_fit$adj.r.squared,
      aic = AIC(fit),
      slope = coef(fit)[2],
      intercept = coef(fit)[1]
    ))
    
    fit <- phylolm(y_value~x_value, phy = shorebird$phy, data = shorebird$data, model = "lambda")
    summ_fit = summary(fit)
    dt_fit = rbind(dt_fit,data.frame(
      model="Pagel's λ",
      p_val_slope = summ_fit$coefficients[2,4],
      r.squared = summ_fit$r.squared,
      adj.r.squared = summ_fit$adj.r.squared,
      aic = AIC(fit),
      slope = coef(fit)[2],
      intercept = coef(fit)[1]
    ))
  }
  
  fit = lm(y_value~x_value)
  summ_fit = summary(fit)
  dt_fit = rbind(dt_fit,data.frame(
    model="LM",
    p_val_slope = summ_fit$coefficients[2,4],
    r.squared = summ_fit$r.squared,
    adj.r.squared = summ_fit$adj.r.squared,
    aic = AIC(fit),
    slope = coef(fit)[2],
    intercept = coef(fit)[1]
  ))
  
  
  dt_fit = dt_fit[order(dt_fit$aic),]
  
  title_graphic = paste(paste("/ ",dt_fit$model,":", 
                              # AIC=",round(dt_fit$aic),
                              " R2=",round(dt_fit$r.squared, 2),
                              "; p-value=",formatC(dt_fit$p_val_slope, format = "e", digits = 0),sep=""),collapse = " ")
  if ( nrow(dt_fit) == 1 ){
    title_graphic = paste(paste("/ ",dt_fit$model,": R2=",round(dt_fit$r.squared, 2),
                                "; p-value=",formatC(dt_fit$p_val_slope, format = "e", digits = 0),sep=""),collapse = " ")
  }
  
  dt_fit = dt_fit[dt_fit$aic == min(dt_fit$aic),]
  
  return(list(no_species=length(x_value),title_graphic=title_graphic,slope=dt_fit$slope,intercept=dt_fit$intercept))
  
}


Clade_color = c(Embryophyta="#33A02C",Diptera="red",Lepidoptera="#FB9A99",Coleoptera="#e2cc1a",Hymenoptera="#ba8e18","Other Insects"="#FF7F00",
                Nematoda="#B2DF8A",Teleostei="#1F78B4",Mammalia="#66281A",Aves="#5b5b5b","Other Vertebrates"="#A6CEE3","Other Metazoans"="#f5b48a", "NA"="grey"
)

table_phylo = read.delim("www/phylogenetic_trees_description.tab")
phylogenetic_trees = paste("www/",table_phylo$name,sep="")
names(phylogenetic_trees) = table_phylo$description

data_by_species = data.frame(species="")
for (file in list.files("www/species_informations_tables",full.names = T,pattern = "data_by_species.tab")){
  dt = read.delim(file,header = T)
  data_by_species = merge(dt,data_by_species, by.x = "species", by.y = "species", all.x = TRUE, all.y = TRUE)
}
data_by_species$clade.qual = factor(data_by_species$clade.qual, levels = names(Clade_color))

data_by_species_original = data_by_species

dt_species = read.delim("www/database/list_species.tab",header=T)
rownames(dt_species) = dt_species$species
all_listNomSpecies = tapply(dt_species$species,dt_species$clade_group,function(x)  str_replace_all(x,"_"," "))
dt_species = dt_species[dt_species$expression_data,]
dt_species$path_db = paste(dt_species$species,"_NCBI.taxid",dt_species$NCBI.taxid,"/",dt_species$assembly_accession,sep="")
listNomSpecies = tapply(dt_species$species,dt_species$clade_group,function(x)  str_replace_all(x,"_"," "))

axisInter = read.delim("www/inter_axis.tab",sep="\t")
axisInter_original = axisInter
ExplicationsInter = axisInter[axisInter$group != "",c("display_label","description")]

axisInter_quantitative = axisInter[axisInter$quantitative,]
axisInter_qualitative = axisInter[!axisInter$quantitative,]


axisInter_list_qualitative = tapply(axisInter_qualitative$display_label,axisInter_qualitative$group,list)
axisInter_list_quantitative = tapply(axisInter_quantitative$display_label,axisInter_quantitative$group,list)



axisIntra = read.delim("www/intra_axis.tab",sep="\t")
axisIntra_list = tapply(axisIntra$display_label,axisIntra$group,list)
ExplicationsIntra = axisIntra[axisIntra$group != "",c("display_label","description")]

set_color_structure_gene = c("Class:major in CDS"="#33A02C","Class:major "="#B2DF8A",
                             "Class:minor in CDS"="#E31A1C","Class:minor "="#FB9A99",
                             "Class:Unclassified in CDS"="#1F78B4","Class:Unclassified "="#A6CEE3")
