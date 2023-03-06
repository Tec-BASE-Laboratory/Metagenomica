## Función para determinar abundancia (normalizar datos) 


otherize <- function(dt,limite,other) {
  totAb <- sum(dt$Abundance)
  dt[,Abundance := 100*Abundance/totAb]
  dt[Abundance < limite, `Scientific Name` := other]
  dt <- aggregate(Abundance ~ ., dt, sum)
}

#Función para generación de barplots

genbarplot <- function(x,limite,tax,name,pal){
  
  other <- paste("Other <", limite,"%",sep = "")
  
  df <- x 
  samples <- unique(df$Sample)
  for (k in samples) {
    dt <- df[Sample == k]
    otherize(dt=dt,limite=limite,other=other)
    if (k == samples[1]) {
      data_barplot <- rbind(dt)
    }
    else {
      data_barplot <- rbind(data_barplot,dt)
    }
  }
  data_barplot <- aggregate(Abundance ~ `Scientific Name` + Sample, data_barplot,sum)
  
  
  data_barplot$`Scientific Name` <- reorder(data_barplot$`Scientific Name`, data_barplot$Abundance)
  
  data_barplot$`Scientific Name` <- factor(data_barplot$`Scientific Name`, levels = rev(c(other,levels(data_barplot$`Scientific Name`)[!levels(data_barplot$`Scientific Name`) == other])))
  
  n <- length(unique(data_barplot$`Scientific Name`))
  
  if(n > length(pal)){
    colpal <- colorRampPalette(pal)(n)
    colpal[n] <- "#BBBBBB"
  } else {
    colpal <- pal
    colpal[n] <- "#BBBBBB"
  }
  
  barplot <- ggplot(data_barplot, 
                    aes(x=Sample, y=Abundance, fill=`Scientific Name`)) + 
    geom_bar(stat="identity") + 
    theme(axis.text.x=element_text(), axis.title.x=element_blank(), panel.background = element_blank()) + 
    ylab("Relative Abundance (%)") +
    scale_fill_manual(values = colpal)
  
  barplot$labels$fill <- str_to_title(tax)
  
  outname <- paste(name,"Barplot",tax,".png",sep = "_")
  
  ggsave(filename=outname,plot=barplot, device = "png", width = 10, height = 7, dpi = 300)
  
}

## Funciones para hacer boxplots de los índices de diversidad Shannon y Simpson 

genDiversityIndexTable <- function(x, tax) {
  shannonS <- c()
  simpsonS <- c()
  df <- x %>% filter(Rank == tax)
  df <- df[,c(-1,-2)]
  df <- na.omit(df)
  samples <- names(df[,2:length(df)])
  
  for (col in 2:length(df)){
    #shannon
    dfcol <- df[,col]
    shannon <- diversity(dfcol, index="shannon")
    shannonS <- c(shannonS,shannon)
    
    #simpson
    dfcol <-df[,col]
    simpson <- diversity(dfcol, index="simpson")
    simpsonS <- c(simpsonS,simpson)
  }
  #tabla samples
  indexes <- data.frame(samples,shannonS,simpsonS)
  names(indexes) <- c("Sample", "Shannon", "Simpson")
  
  
  indexes
}


genboxplots <- function(indexesT,name){
  #Shannon Index
  shannonboxname <- paste(name,"_ShannonIndexBoxplot", ".png", sep="")
  shannonbox <- ggplot(indexesT, aes(x= Sample, y= Shannon)) + 
    ylab("Shannon Index") + 
    geom_boxplot(fatten=1, outlier.shape = NA) +
    theme_light() + 
    theme(axis.title.x = element_blank())
  ggsave(filename = shannonboxname, plot = shannonbox, device="png", width = 10, height = 10, dpi = 300)
  
  #Simpson Index
  simpsonboxname <- paste(name,"_SimpsonIndexBoxplot", ".png", sep="")
  simpsonbox <- ggplot(indexesT, aes(x= Sample, y= Simpson)) + 
    ylab("Simpson Index") + 
    geom_boxplot(fatten=1, outlier.shape = NA) + 
    theme_light() + 
    theme(axis.title.x = element_blank())
  ggsave(filename = simpsonboxname, plot = simpsonbox, device="png", width = 10, height = 10, dpi = 300)
}

## Paleta de colores general

bcol.pal <- 
  mypal <- c("#f29080", "#f27933", "#ffba4a", "#fbff91", "#78a644", "#91ffd9", "#77def2", "#396799", "#91a0ff", "#c936ff", "#99268f", "#ff369e", "#ff363c", "#8c3b1d", "#a6754b", "#99844b", "#b9ff40", "#30e676", "#188c7d", "#2bc0ff", "#1a4099", "#3643ff", "#e89cff", "#cc74a6", "#bf284b")

## Preparar librerías

packages <- c("data.table", "tidyverse", "vegan", "RColorBrewer","colorspace")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

## Espacio destinado para dar colores específicos a las variables 

# RB1 = "#00FE08"
# RB2 = "#0B6B0E"
# MF1 = "#04E2F8"
# MF2 = "#1114DC"
# CMF1 = "#CB72F0"
# CRB1 = "#6007D2"

## Preparación de datos 

setwd("~/Tec_BASE/ProyectosTCB/MetGen_Totoaba//") 
raw_data <- read.table(file.choose(), header = T, sep = ",",quote = "\"", stringsAsFactors = F, fill = F) 
raw_data
raw_data[,-c(1,2,3)] <-lapply(raw_data[, -c(1,2,3)], as.integer)
#raw_data <- raw_data[,c(-22,-23,-24)]
str(raw_data)

taxa_data <- list()
taxa_long_data <- list()
Metatop <- list()
taxa <- c("phylum", "family","genus", "species") 

## Filtrado de datos 
for (i in 1:length(taxa)) {  
  taxa_data[[i]] <- raw_data %>% filter(Rank == taxa[i]) 
  taxa_data[[i]] <- taxa_data[[i]][,c(-1,-2)] 
  taxa_long_data[[i]] <- taxa_data[[i]] %>% 
    as.data.table() %>% 
    melt(id = c("Scientific.Name"), variable.name = "Sample", value.name = "Abundance") %>% 
    filter(Abundance != 0) 
  colnames(taxa_long_data[[i]]) <- c("Scientific Name","Sample","Abundance") 
  #df$Sample <- df$Sample %>% str_remove_all('[r\\d]')  ## Opcional ##
  treatments <- unique(taxa_long_data[[i]]$Sample) 
  for (k in treatments) { 
    dt <- taxa_long_data[[i]][Sample == k] 
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

## Generación de documentos a partir del filtrado de los datos 

for (i in 1:4) {
  name = paste("MetaGen",taxa[i],sep = "_")
  write_csv(taxa_data[[i]],file=paste("Taxa_data",name,".csv",sep="_"))
  write_csv(taxa_long_data[[i]],file=paste("Taxa_long",name,".csv",sep="_"))
  write_csv(Metatop[[i]],file=paste("MetaTop",name,".csv",sep="_"))
}


## Visualización de datos
##Tablas y Boxplots 

for (i in 1:length(taxa)){
  x  <- taxa_long_data[[i]]
  genbarplot(x,1,taxa[i],"MetaGen",mypal)
}

#diversity tables and boxplots

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

for (i in 1:length(diversity)){
  name = paste("MetaGen",taxa[i])
  genboxplots(diversityT[[i]],name)
  write.csv(diversity[[i]],file=paste(name,".csv"))
}

