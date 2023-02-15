#utils kks

##barplots

RB1 = "#00FE08"
RB2 = "#0B6B0E"
MF1 = "#04E2F8"
MF2 = "#1114DC"
CMF1 = "#CB72F0"
CRB1 = "#6007D2"

otherize <- function(dt,limite,other) {
  #dt = data table, limite = threshold for changing scientific name to other
  totAb <- sum(dt$Abundance)
  dt[,Abundance := 100*Abundance/totAb]
  dt[Abundance < limite, `Scientific Name` := other]
  dt <- aggregate(Abundance ~ ., dt, sum)
}


genbarplot <- function(x,limite,tax,name,pal){
  ## x = data table long form
  ## limite = limite en porcentaje para cambiar scientific name a other
  ## tax = str of tax to make barplot of
  ## name = distinct name for file

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
    theme(legend.text = element_text(size=8), legend.key.size = unit(8, 'mm'), axis.text.x=element_text(size = 7), axis.title.x=element_blank(), panel.background = element_blank()) + 
    ylab("Relative Abundance (%)") +
    scale_fill_manual(values = colpal) 
  
  barplot$labels$fill <- str_to_title(tax)
  
  outname <- paste(name,"Barplot",tax,".png",sep = "_")
  
  ggsave(filename=outname,plot=barplot, device = "png", width = 15, height = 7, dpi = 300)
  
}

##diversity

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
    scale_fill_brewer(palette = "Spectral") + 
    theme_light() + 
    theme(axis.title.x = element_blank())
  ggsave(filename = shannonboxname, plot = shannonbox, device="png", width = 10, height = 10, dpi = 300)
  
  #Simpson Index
  simpsonboxname <- paste(name,"_SimpsonIndexBoxplot", ".png", sep="")
  simpsonbox <- ggplot(indexesT, aes(x= Sample, y= Simpson)) + 
    ylab("Simpson Index") + 
    geom_boxplot(fatten=1, outlier.shape = NA) + 
    scale_fill_brewer(palette = "Spectral") +
    theme_light() + 
    theme(axis.title.x = element_blank())
  ggsave(filename = simpsonboxname, plot = simpsonbox, device="png", width = 10, height = 10, dpi = 300)
}

