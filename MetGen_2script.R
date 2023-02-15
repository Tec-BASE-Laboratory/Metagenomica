## Correr MetGen_Utils antes de correr este código 

# Instalación y carga de librerias
packages <- c("data.table", "tidyverse", "vegan", "RColorBrewer","colorspace")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
  }
invisible(lapply(packages, library, character.only = TRUE))

#Preparación de los datos 
setwd("C:/Users/benja/Desktop/Shared/TeBase_2023") 
raw_data <- read.table(file.choose(), header = T, sep = "\t",quote = "\"", stringsAsFactors = F, fill = F) ## Dos Warnings, aún hay que checar el porque de estos y porque no sale el output como se desearía
raw_data
raw_data[,-c(1,2,3)] <-lapply(raw_data[, -c(1,2,3)], as.integer)

str(raw_data)


taxa_data <- list()
taxa_long_data <- list()
Metatop <- list()
taxa <- c("phylum", "family","genus", "species") 
taxa_data 
taxa 

# En esta primera función se hace el filtrado de los datos por Rank, por Abundance y el top (10,20, etc). Primero se filtran los datos comparandolos contra la columna "Rank" y los resultados son depositados en la lista

for (i in 1:length(taxa)) {  
  taxa_data[[i]] <- raw_data %>% filter(Rank == taxa[i]) 
  taxa_data[[i]] <- taxa_data[[i]][,-1] 
  taxa_long_data[[i]] <- taxa_data[[i]] %>% 
    as.data.table() %>% 
    melt(id = c("Scientific.Name"), variable.name = "Sample", value.name = "Abundance") %>% 
    filter(Abundance != 0) 
  df <- taxa_long_data[[i]] 
  colnames(df) <- c("Scientific Name","Sample","Abundance") 
  #df$Sample <- df$Sample %>% str_remove_all('[r\\d]')  ## Opcional 
  taxa_long_data[[i]] <- df 
  treatments <- unique(df$Sample) 
  for (k in treatments) { 
    dt <- df[Sample == k] 
    dt <- aggregate(Abundance ~ ., dt, sum) 
    dt <- dt %>% arrange(desc(Abundance)) %>% slice(1:10) 
    RelAb <- sum(dt$Abundance) 
    dt <- as.data.table(dt)
    dt[,Abundance := 100*Abundance/RelAb]
    if (k == treatments[1]){
      top_df <- rbind(dt)
    }
    else {
      top_df <- rbind(top_df,dt)
    }
    
  }
  Metatop[[i]] <- top_df
}

names(taxa_data) <- taxa 
names(taxa_long_data) <- taxa
names(Metatop) <- taxa
dt
df
taxa_data
taxa_long_data
Metatop

for (i in 1:4) {
  name = paste("MetaGen",taxa[i])
  write_csv(taxa_data[[i]],file=paste("taxa_data",name,".csv"))
  write_csv(taxa_long_data[[i]],file=paste("taxa_long",name,".csv"))
  write_csv(Metatop[[i]],file=paste("MetaTop",name,".csv"))

}


#color palette and filepath
library("colorspace")
pal <- choose_palette()
mypal<- c(pal(7))

#bcol.pal <- 
  #mypal <- c("#f29080", "#f27933", "#ffba4a", "#fbff91", "#78a644", "#91ffd9", "#77def2", "#396799", "#91a0ff", "#c936ff", "#99268f", "#ff369e", "#ff363c", "#8c3b1d", "#a6754b", "#99844b", "#b9ff40", "#30e676", "#188c7d", "#2bc0ff", "#1a4099", "#3643ff", "#e89cff", "#cc74a6", "#bf284b")


#create barplots, function from utils_kks.R

for (i in 1:length(taxa)){
  x  <- Metatop[[i]]
  genbarplot(x,1,taxa[i],"MetaGen",mypal)
}

#diversity tables and boxplots
boxplotpath <- ""

diversity <- list()

for (i in 1:length(taxa)) {
  diversity[[i]] <- genDiversityIndexTable(raw_data, tax = taxa[i])
}

names(diversity) <- taxa
diversityT <- list()

for (i in 1:length(diversity)) {
  df <- diversity[[i]]
  df$Sample <- df$Sample %>% str_remove_all('[r\\d]') 
  diversityT[[i]] <- df
  rm(df)
}
names(diversityT) <- taxa

setwd(boxplotpath)
for (i in 1:length(diversity)){
  name = paste("MetaGen",taxa[i])
  genboxplots(diversityT[[i]],name)
  write.csv(diversity[[i]],file=paste(name,".csv"))
}

