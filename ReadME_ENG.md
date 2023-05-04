# Metagenomics  
This project was developed by for academic use in the Core Lab Genomics investigation group at the Instituto Tecnológico y de Estudios Superiores de Monterrey (ITESM).

This program was designed for making a metagenomic analysis from 16s and shotgun sequencing. It is made from TXT or CSV format data from a taxonomic analysis from [KRAKEN 2](https://github.com/DerrickWood/kraken2.).

This program aims to create statiscal analyses and allow the user the visualization of taxonomical data from KRAKEN 2. 

# Demo Data

If you want to see how the code works with some demo data, you can download it on the repository. The document is called taxonomic_clasiffication_ReadMe.csv

# Libraries and Utils 

The firts part of the code is designed to declare different functions that will be used for the creation of barplots and boxplots that will be part of the results and analysis. 

It's important to run these first lines because some functions are not part of any library and they are necessary for the program to work correctly.

Aside from that, some libraries are installed too, these will be used to make different statistical anlysis. 


On the line:
```Rscript
geom_boxplot(fill=sample(mypal, length(treatcol)),fatten=1, outlier.shape = NA)
```
Random colors will be generated for the boxplots from Shannon and Simpson indexes. It is possible to change the parameter `fill =` in order to assign specific colors depending on the number of variables. This must be done following this format `fill = c(C1,C2,C3,...`. If the user wants only one color for the boxplots, it must be done like this  `fill = *hexadexcimal color code*`. If the user wants to have blank boxplots just needs to eliminate the `fill =` portion of the code. This will leave the code looking like this.

```Rscript
geom_boxplot(fill=sample(mypal, length(treatcol)),fatten=1, outlier.shape = NA)
```

It is recommended to use the hexadecimal code for the colors in order to get more specific and personalized colors depending on the study. If specificity is not required, the names of the colors can be used instead. 

### Nota 
Is very important to mention that the modification of this parameters would depend solely on the requirements of the study. It depends on what needs to be shown or analyzed. 

# Data
In the second section of the code, the pathway for the directory in which the analysis is going to take place is delcared. It's important to mention that every output of the analysis in the form of tables and images is goint to be saved inside this directory.  That's why it's recommended to have different Directories for each individual analysis. 

This section of the code is destined to treat and format the data in order to avoid any type of errors or warnings. The input for this code must be a **.csv** or a **.txt** document from a taxonomic analysis using [KRAKEN 2](https://github.com/DerrickWood/kraken2.). 


The input must have the following layout: 


| Rank | TaxId | Scientific Name | Sample 1 | Sample 2 | Sample 3 | 
| --- | --- | --- | --- | --- | --- |
| unclassified  | 0 | Uknown | 68 | 1232 | 1696 |
| superkingdom | 2 | Bacteria (bacteria) | 66708 | 76666 | 64937 | 
| genus | 131079 | Limnobacter| 0 | 0 | 0 |
| species | 2060312 | Altererythrobacter sp. B11 | 0 | 0 | 2 |

It is recommended that only the names of the samples are modified to avoid any issues or inconvenients with the program. The sample names can be modified freely depending on the analysis and what needs to be observed. To access more information regarding on how the data is treated and filtered, go to [**Filtering**](#Filtering).

In the line:
```Rscript
raw_data <- read.table(file.choose(), header = T, sep = "\t",quote = "", stringsAsFactors = F, fill = F)
```

The `"sep ="` parameter must be modified depending on the type of document that is used: for **.txt** `"\t"` and for **.csv** `","`. 

The following code lines are directed to assure that the data is inputted correctly and there's no data loss or change in values (numbers to strings) and to avoid the existence of NA values insde the table.

The line: 
```Rscript
raw_data <- raw_data[,c(n1,n2,n3)]
```
it's used to eliminate columns that might be added in previous steps, wheter it was when loadind the data or converting the numbers. This line can be modified depending on how many columns were added. E.g. if the data frame had 15 columns and now has 17 columns, columns 16 and 17 were added. In this case, the line must be modified with the following `raw_data <- raw_data[,-c(16,17)]`. In this way, the two columns that were added by coercion are eliminated. This does not mean that there is data loss, sometimes KRAKEN creates documents with invisible rows that R reads as nameless empty rows filled with NA's. If this columns are not eliminated, the program won't work as intended. 

The lines:
  
```Rscript
raw_data
str(raw_data)
```

are meant to allow the user to visualize the data and how was read by R and see if there's any problem with the data. If everything was loaded sucessfully, when running `srt(raw_data)` the data must have the following format: 
            
```
 $ Rank           : chr  "unclassified" "superkingdom" "superkingdom" "superkingdom" ...
 $ TaxId          : int  0 2 10239 2157 2759 200783 200795 200918 28890 200930 ...
 $ Scientific.Name: chr  "Unknown" "Bacteria <bacteria>" "Viruses" "Archaea" ...
 $ Sample 1       : int  33369 25062 0 1 0 0 1 1 1 0 ...
 $ Sample 2       : int  184966 73818 0 1 0 1 39 20 1 3 ...
 $ Sample 3       : int  329303 75747 0 0 0 0 29 0 0 0 ...
```
The following lines are related to the creation of empty lists for making the next step in the analysis, the data filtering based on 4 taxonomic groups. 
            
# Filtering

In this part of the code the user can find a loop that filters the data depending on the taxonomic group they belong to **"phyllum"**, **"family"**, **"genus"**, or **"species"**. In this section, there are various personalizable parameters to satisfy the requirements of various studies. One of these parameters is the one in the next line:

```Rscript
#taxa_long_data[[i]]$Sample <- taxa_long_data[[i]]$Sample  %>% str_remove_all('[r\\d]') 
```
This line allows the user to group the samples with its own replicates. This line is commented by default to avoid grouping of samples with similar names. 
The following data frame is is used to give an example on how data is grouped. 

| Rank | TaxId | Scientific Name | SP1 | SP2 | SP3 | SP4
| --- | --- | --- | --- | --- | --- | --- |
| unclassified | 0 | Uknown | 68 | 1232 | 1696 | 22013 |
| superkingdom | 2 | Bacteria (bacteria) | 66708 | 76666 | 64937 | 25496 |
| genus | 131079 | Limnobacter| 0 | 0 | 0 | 90 |
| species | 2060312 | Altererythrobacter sp. B11 | 0 | 0 | 2 | 10 |

When filtering this dataframe by grouping, the results are not going to be shown for individual samples as "SP1, SP2, SP3, SP4", instead, they are going to be grouped into sample SP and replicates. In this case, the line must be kept as default in order to show each sample as an individual. For grouping samples, the user must define the names of the samples and its replicates using the follwoing format "A1, A2, B1, B2". In this case the user has two samples with two replicates. The code will show results from samples "A" and "B" taking each replicate with its sample. For further information and more examples on how data can be grouped go to [**Boxplot Creation**](#Boxplot-Creation)


Inside this loop the data is filtered to obtain a **Top n** of the results. Which means the program obtains a Top "n" from the different taxonomical data from species, genus, family and phyllum. The defalut is a **Top 10** but it can be modified depending on the study and what needs to be analyzed. 

This modification can be done on the line:
```Rscript
dt <- dt %>% arrange(desc(Abundance)) %>% slice(1:10)  
```
The user must change `slice(1:n)`, to obtain a different Top. 

# Data Visualization 

In this section different outputs are created in the form of **.csv** documents and tables, barplots and boxplots that allow the user the visualization of the data. This outputs are helpful to make further statistical analysis because they have been filtered into different groups.

## Table Creation 

In this part of the code, different tables are created after the filtering of the data. The code can be modified to make personalized outputs to give specific names to the individual analysis. 

```Rscript
name = paste("MetaGen_ReadME",taxa[i],sep = "_") 
```
This line can be modified by adding specific identifiers to make unique documents for each analysis. This can be done modifying `"MetaGen_ReadME"` to make something personalized for the research. E.g. `"MetaGen_Analysis1"`.

After modifying this part of the code and running the loop, twelve **.csv** documents are going to be created, 3 documents per taxonomic group (**"phyllum"**, **"family"**, **"genus"**, and **"species"**). The documents called `"Taxa_Data"` have the full table of counts for each taxonomic group on each sample. the documents called "Taxa_LD" have the abundance of each taxonomic group for each sample. The documents called `"MetaTop"` is similar to `"Taxa_LD"` but it only has the **Top n** of each taxonomic group based on their abundance.  

The outputs of these documents should look like the following tables:

**For MetaTop (Family)**

| Scientific Name | Sample | Abundace
| --- | --- | ---
| Pseudomonadaceae | SP1 | 61.0684020452676
| Oxalobacteraceae | SP1 | 14.191725962741800
| Micrococcaceae | SP1 | 8.821110594629390
| Planococcaceae | SP1 | 4.172775056316370

**For Taxa_Data (Family)**

| Scientific.Name | SP1 | SP2 | SP3 | SP4
| --- | --- | --- | --- | ---
| Chthonomonadaceae | 0 | 0 | 0 | 2
| Salinivirgaceae | 0 | 0 | 0 | 0
| Hyphomonadaceae | 0 | 0 | 0 | 0
| Archangiaceae | 0 | 1 | 0 | 0

**For Taxa_LD (Family)**

| Scientific Name | Sample | Abundace
| --- | --- | --- 
| Leptospiraceae | SP1 | 2
| Methylobacteriaceae | SP1 | 2
| Burkholderiaceae | SP1 | 893
| Selenomonadaceae | SP1 | 5

Further down in the code, there are four documents more that get created. In the last loop, the program performs an statistical analysis to determine Shannon and Simpson indexes. The **.csv** documents that are created have tables with the Shannon and Simpson indexes for each sample in each taxonomic group.

These documents should be presented with the following format:

| Sample | Shannon | Simpson
| --- | --- | --- 
| SP1 | 2.49878742275665 | 0.826907066675482
| SP2 | 2.25058905458163 | 0.774696960815804
| SP3 | 1.97042119751292 | 0.678785150456466
| SP4 | 3.24981277477337 | 0.929375693849251

These tables are meant to allow the user use the data for further statistical analysis. 

## Boxplot Creation 

After the first twelve **.csv** documents are created, the following parts of the code are dedicated to the creation of barplots and boxplots. 

The barplots are created using the fuction barplots from the utils section of the code. The plots are generated for each taxonomic group, at the end four **.png** files should be generated for the barplots.

Here's an example of a barplot generated using the program: 

![MetaGen_ReadME_species_Barplot_family_](https://user-images.githubusercontent.com/106264605/236337607-828736a6-ed3b-463b-8c50-4b7581c6d3ca.png)

This barplot shows the relative abundance of each family in each sample. 

The boxplot fuction can be found in the last section of the code. This loop creates the Shannon and Simpson indexes tables and two boxplots for each taxonomic group (one for Shannon and one for Simpson). 

Here are examples of the boxplots generated by the program for Shannon (first) and Simpson (second).

![S SDiversity_ReadME_family_ShannonIndexBoxplot](https://user-images.githubusercontent.com/106264605/236338339-5b07228b-9875-4b88-8e45-61c2fca29c0b.png)

![S SDiversity_ReadME_family_SimpsonIndexBoxplot](https://user-images.githubusercontent.com/106264605/236338358-1b2b73b8-2500-47fb-a235-bb59c6463e82.png)

The colors for the boxplots are generated randomly, this can be changed following the steps that can be found in [**Libraries and Utils**](#Libraries-and-utils)
