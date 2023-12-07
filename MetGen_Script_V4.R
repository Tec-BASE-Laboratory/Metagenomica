## Función para agrupar los datos que restan (otro 1%)

other <- function(data) {
  for (group in taxa) {
    current_df <- data[[group]]
    current_df$`Scientific Name`[current_df$Abundance < 1] <- "Other" # El valor se puede cambiar a criterio
    data[[group]] <- current_df
    # Optionally, if you want to print the modified current_df
    print(current_df)
  }
  return(data)
}


#Función para caluclar los índices de diversidad

calculate_diversity_indices <- function(data) {
  # Shannon
  shannon_diversity <- function(x) {
    p <- x / sum(x)
    -sum(p * log2(p), na.rm = TRUE)
  }
  
  #Simpson
  simpson_diversity <- function(x) {
    p <- x / sum(x)
    1 - sum(p^2, na.rm = TRUE)
  }
  
  # Agrupar los datos por muestra y calcular los índices
  
  result <- data %>%
    group_by(Sample) %>%
    summarize(
      Shannon_Index = shannon_diversity(Abundance),
      Simpson_Index = simpson_diversity(Abundance)
    )
  
  return(result)
  
}

## Funciones para hacer boxplots de los índices de diversidad Shannon y Simpson

calculate_abundance <- function(abundance_data, index_data) {
  merged_data <- merge(abundance_data, index_data, by = "Sample", all.x = TRUE)
  merged_data <- merged_data %>%
    mutate(
      Simpson_Diversity = Simpson_Index * Abundance,
      Shannon_Diversity = Shannon_Index * Abundance
    )
  return(merged_data)
}



# --------------------------Librerías ------------------------------------------

packages <- c("data.table", "dplyr","tidyverse", "vegan", "RColorBrewer","colorspace","viridis","tidyr")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))


## -------------------  Preparación de datos  ------------------------------ 

dir <- "/some/directory/path" # <- Directorio, archivos de metagenómica
barplotpath <- "/some/barplot/path" # <- Directorio para guardar los barplots
boxplotpath <- "/some/boxplot/path" # <- Directorio para guardar boxplots

setwd(dir) 

raw_data <- read.table("", header = T, sep = ",",quote = "", stringsAsFactors = F, fill = F) # <- Archivo de Metagenómica

raw_data$TaxId = NULL
str(raw_data)

#treatments 

#En esta primera parte se creará una lista vacía llamada "taxa_data" en la cual se van a depositar los valores filtrados de "raw_data". Este filtro será "taxa", el cual está separando por grupo taxonómico. 

taxa_data <- list()
taxa_long_data <- list()# Se crean listas vacías 
Metatop <- list()
taxa <- c("phylum", "family","genus", "species") #Se crea una lista llamada "taxa" contra la cual se va a comparar los datos para poder hacer el filtrado 

taxa_data # visualización de taxa data,  está vacía 

taxa #visualización de la lista, contiene 4 elementos 

# -------------


# En esta primera función se hace el filtrado de los datos por Rank, por Abundance y el top (10,20, etc). Primero se filtran los datos comparandolos contra la columna "Rank" y los resultados son depositados en la lista 

for (i in 1:length(taxa)) { #por cada elemento en taxa se va a hacer un filtrado 
  taxa_data[[i]] <- raw_data %>% filter(Rank == taxa[i]) # Se agragan datos de raw_data a taxa_data a partir de un filtrado usando taxa
  taxa_data[[i]] <- taxa_data[[i]][,c(-1)] #por cada elemento de taxa data, se elimina la columna 1 (columna de "Rank") (en caso de tener otras columnas como TaxId es necesario concatenar las columnas que se desean eliminar: [,c(-1,-2,-...)])
  taxa_long_data[[i]] <- taxa_data[[i]] %>% #Se crean elementos dentro de taxa_treatments_data mediante un filtrado de taxa_data 
    as.data.table() %>% # Esta funcion agrupa los datos de la lista taxa_data y son convertidos en columnas y celdas agrupadas 
    melt(id = c("Scientific.Name"), variable.name = "Sample", value.name = "Abundance") %>% #Esta función ayuda a reorganizar los datos del filtro anterior y los agrupa en varias columnas, en este caso evalúa Scientific.name y los agrupa por sample y abundance. 
    filter(Abundance != 0) #Aquí se eliminan todos los datos de abundance que sean 0, de esta manera se evitan datos irrelevantes. 
  df <- taxa_long_data[[i]] # Se define una variable que ayude a usar colnames por cada elemento dentro de taxa_treatments_data
  colnames(df) <- c("Scientific Name","Sample","Abundance") #Se le agregan los nombres a las columnas de cada elemento de taxa_treatment_data
  #df$Sample <- df$Sample %>% str_remove_all('[r\\d]')  <---- esta parte del código es útil si es que se tienen varias réplicas y se requiere agrupar todos los elementos de una muestra en una sola barra, en este caso, al ser 7 muestras diferentes no es necesario usarlo
  taxa_long_data[[i]] <- df #Se agregan los nombres del data frame a la lista taxa_long_data
  #esta siguiente parte del código es para llenar la lista de Metatop, en donde se obtienen los topX de cada grupo taxonómico y se compara su abundancia para poder 
  treatments <- unique(df$Sample) #aquí se obtienen todas las variables únicas que existen en la columna sample de taxa_long_data
  for (k in treatments) { 
    dt <- df[Sample == k] #Se compara cada elemento de la variable treatments contra sample y se define la variable dt
    dt <- aggregate(Abundance ~ ., dt, sum) #se utiliza la función aggregate para poder obtener la abundancia relativa  
    dt <- dt %>% arrange(desc(Abundance)) %>% slice(1:10) #dt se filtra utilizando la función arrange y desc para poder acomodar los datos de mayor a menor y después slice hace el corte del topX (en este caso, un slice para el top10 de 1:10)
    dt <- dt %>% filter(!grepl("TaxID",Sample))
    RelAb <- sum(dt$Abundance) #se define una variable de la abundancia relativa que suma la abundancia de la variable dt
    dt <- as.data.table(dt) # se hace una tabla de datos 
    dt[,Abundance := 100*Abundance/RelAb] #se convierte la columna de abundancia a una columna con los porcentajes de abundancia relativa 
    if (k == treatments[1]){
      top_df <- rbind(dt)
    }
    else {
      top_df <- rbind(top_df,dt)
    }
    
  }
  Metatop[[i]] <- top_df
}

names(taxa_data) <- taxa # Se le dan nombres a los elementos de taxa_data, taxa_long_data y Metatop a partir de taxa 
names(taxa_long_data) <- taxa
names(Metatop) <- taxa
dt
df
taxa_data # Las listas ya no se encuentran vacías y ahora son listas de 4 elementos con data frames diferentes, acomodados dependiendo del tipo de filtro que se agregó.
taxa_long_data
Metatop

## ----------------  Generación de documentos a partir del filtrado de los datos --------------
setwd()
for (i in 1:4) {
  name = paste("MetaGen_5DIC",taxa[i],sep = "_")
  write_csv(taxa_data[[i]],file=paste("Taxa_dataS",name,".csv",sep="_"))
  write_csv(taxa_long_data[[i]],file=paste("Taxa_LD",name,".csv",sep="_"))
  write_csv(Metatop[[i]],file=paste("MetaTop",name,".csv",sep="_"))
}

#  -------------------- Preparación de los datos: TOP por cada muestra ----------------- 

long_taxa_data <- list()
for (i in 1:4) {
  long_taxa_data[[i]] <- taxa_data[[i]] %>% 
    as.data.table() %>% 
    melt.data.table(id = c("Scientific.Name"), 
                    variable.name = "Sample", 
                    value.name = "Abundance")
  df <- long_taxa_data[[i]]
  colnames(df) <- c("Scientific Name","Sample","Abundance")
  df$Sample <- gsub("_(\\d+)$", "", df$Sample)
  long_taxa_data[[i]] <- df
  rm(df)
}

#for (i in 1:4) {
long_taxa_data[[i]] <- taxa_data[[i]] %>% 
  as.data.table() %>% 
  melt.data.table(id = c("Scientific.Name"), 
                  variable.name = "Sample", 
                  value.name = "Abundance")
df <- long_taxa_data[[i]]
colnames(df) <- c("Scientific Name","Sample","Abundance")
df$Sample <-  df$Sample %>% str_remove_all('[r\\d]') 
long_taxa_data[[i]] <- df
rm(df)
}
names(long_taxa_data) <- taxa

top10 <- list() # <-  lISTA DE TOP 10 por grupo taxonómico
for (i in 1:4) {
  df <- long_taxa_data[[i]]
  treatments <- unique(df$Sample)
  for (k in treatments) {
    dt <- df[Sample == k] 
    dt <- aggregate(Abundance ~ ., dt, sum)
    dt <- dt %>% arrange(desc(Abundance)) %>% slice(1:10)
    dt <- dt %>% filter(!grepl("TaxID",Sample))
    
    totAb <- sum(dt$Abundance)
    dt <- as.data.table(dt)
    dt[,Abundance := 100*Abundance/totAb]
    if (k == treatments[1]){
      top_df <- rbind(dt)
    }
    else {
      top_df <- rbind(top_df,dt)
    }
    
  }
  top10[[i]] <- top_df
}
names(top10) <- taxa


## BARPLOTS.

setwd(barplotpath) # Directorio para guardar los barplots

# Barplots por cada grupo taxonómico.

for (i in 1:4) {
  #color palette
  dt <- top10[[i]]
  num <- length(unique(dt$`Scientific Name`))
  if(num > length(mypal)){
    colpal <- colorRampPalette(mypal)(num)
  } else {
    colpal <- mypal
  }
  barplot <- ggplot(dt, aes(x=Sample, y=Abundance, fill=fct_reorder(`Scientific Name`,Abundance, .desc=TRUE)))+
    geom_bar(stat="identity") +
    
    theme(axis.text.x=element_text(), 
          axis.title.x=element_blank(), 
          panel.background = element_blank()) + 
    ylab("Relative Abundance (%)") +
    scale_fill_manual(values = colpal)
  barplot$labels$fill <- str_to_title(taxa[i])
  print(barplot)
  outname <- paste("Datos_taxonómicos_2",taxa[i],".png",sep = "")
  ggsave(filename = outname, plot = barplot, device = "png", width = 30, height = 20, dpi = 300)
  rm(dt)
}

#  -------------------------------   BOXPLOTS  ------------------------------------

# Preparación de datso: Cálculo de los índices y lista de objetos
other(top10) # <- Función para clasificar especies de baja abundancia como "Other"
df_family <- calculate_abundance(top10$family, indices_fa)
df_genus <- calculate_abundance(top10$genus, indices_g)
df_especie <-  calculate_abundance(top10$species, indices_sp)
df_phylum <- calculate_abundance(top10$phylum,indices_ph)
boxplots_df <- list(df_family,df_genus,df_especie,df_phylum)
titulos <- c("Familia", "Género", "Especie", "Filo")

#-------------- Boxplots usando los índices de Simpson -----------
setwd(boxplotpath)
getwd()
for (i in (1:4)){
boxplot_simpson <- ggplot(boxplots_df[[i]], aes(x = Sample, y = Simpson_Diversity, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = mypal) +
  labs(
    x = "Muestra",
    y = "Simpson Abundance",
    fill = "Sample",
    title = paste("Abundancia Simpson,",titulos[i])
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "white", size = 0.2), # Add major grid lines
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
print(boxplot_simpson)
ggsave(filename = paste0(boxplotpath, "plot_", i, ".png"), plot = boxplot_simpson)

}

# -------- Boxplots usando los índices de Shannon -----------------

for (i in (1:4)){
  boxplot_shannon <- ggplot(boxplots_df[[i]], aes(x = Sample, y = Shannon_Diversity, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = mypal) +
  labs(
    x = "Muestra",
    y = "Shannon Abundance",
    fill = "Sample",
    title = paste("Abundancia Shannon,",titulos[i])
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "white", size = 0.2), # Add major grid lines
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
print(boxplot_shannon)
ggsave(filename = paste0(boxplotpath, "plot_", i, ".png"), plot = boxplot_shannon)
}

