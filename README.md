# Metagenómica 
Proyecto desarrollado por alumnos y docentes del Laboratorio Nacional de Secuenciación Genomica del Instituto Tecnológico y de Estudios superiores de Monterrey.

Este programa fue diseñado para hacer un análisis metagenómico de 16S y shotgun, a partir de datos en formato TXT o CSV provenientes de un análisis taxonómico de [KRAKEN 2](https://github.com/DerrickWood/kraken2.).

Este programa busca hacer análisis estadísticos y permitir al usuario la visualización de los datos taxonómicos arrojados por KRAKEN 2. 

# Utils y librerías 
La primera parte del código está diseñada para declarar las diferentes funciones que serán utilizadas para la creación de barplots y boxplots que formarán parte de los resultados de este análisis.

Es importante correr estas lineas del código ya que algunas de las funciones no se encuentran en ninguna librería y son necesarias para el funcionamiento del programa.

A su vez, también se instalan (en caso de no tenerlas) y se llama a las librerías que serán utilizadas para hacer los diferentes análisis estadísticos.  

En la línea:
```Rscript
geom_boxplot(fill=sample(mypal, length(treatcol)),fatten=1, outlier.shape = NA)
```
se generan colores aleatorios para los boxplots de los indices shannon y simpson. Es posible cambiar el parámetro `fill =` para poder poner colores específicos a las variables. Esto se debe de hacer de la siguiente forma: `fill = c(C1,C2,C3,...`. Se recomienda usar los códigos de hexadecimales de color para poder obtener colores más específicos y personalizados. 

# Datos
En la segunda sección del código se declara el pathway de la carpeta donde se encuentran los datos con los cuales va a trabajar el programa. Es importante mencionar que dentro de esta carpeta se van a guardar los diferentes documentos e imágenes resultantes de este análisis. Por esto se recomienda destinar Carpetas específicas para cada análisis indívidual. 

Dentro de esta sección se hace un tratado de los datos para asegurar que el formato es el correcto y R o Rstudio no arroje códigos de error y códigos de advertencia. El input de este código debe de ser un archivo **.csv** o **.txt** proveniente del análisis taxonómico de [KRAKEN 2](https://github.com/DerrickWood/kraken2.). 

El input se debe de ver acomodado de la siguiente manera: 


| Rank | TaxId | Scientific Name | Sample 1 | Sample 2 | Sample 3 | 
| --- | --- | --- | --- | --- | --- |
| unclassified  | 0 | Uknown | 68 | 1232 | 1696 |
| superkingdom | 2 | `Bacteria <bacteria>` | 66708 | 76666 | 64937 | 
| genus | 131079 | Limnobacter| 0 | 0 | 0 |
| species | 2060312 | Altererythrobacter sp. B11 | 0 | 0 | 2 |

Es recomendable modificar solamente los nombres de las muestras (samples) para evitar cualquier inconveniente con el programa. Los nombres de las muestras pueden ser modificados libremente dependiendo de lo que se busque observar en el estudio. Más información de como se tratan, filtran y agrupan las muestras puede ser observado en la sección de [**Filtrado**](#Filtrado)
  
En la línea:
```Rscript
raw_data <- read.table(file.choose(), header = T, sep = "\t",quote = "", stringsAsFactors = F, fill = F)
```

Se debe de modificar el parámetro `"sep ="` dependiendo del tipo de documento que se esté usando: para **.txt** se usa `"\t"` y para **.csv** se usa `","`.

Las siguientes líneas de este código están dirigidas para asegurar que los datos se hayan introducido de manera correcta y no haya una pérdida de datos o se conviertan valores numéricos a strings y evitar que existan NA dentro de la tabla con la que se trabajará. Las líneas:
  
Rscript
raw_data
str(raw_data)

funcionan para poder visualizar como fueron cargados los datos a R y ver si es que existe algún problema con los mismos. Si todo fue cargado de manera correcta se deberían de ver los datos de la misma manera que se ve el documento de origen. En el caso de `srt(raw_data)` se deben de ver los datos de la siguiente manera: 
            
```
 $ Rank           : chr  "unclassified" "superkingdom" "superkingdom" "superkingdom" ...
 $ TaxId          : int  0 2 10239 2157 2759 200783 200795 200918 28890 200930 ...
 $ Scientific.Name: chr  "Unknown" "Bacteria <bacteria>" "Viruses" "Archaea" ...
 $ Sample 1       : int  33369 25062 0 1 0 0 1 1 1 0 ...
 $ Sample 2       : int  184966 73818 0 1 0 1 39 20 1 3 ...
 $ Sample 3       : int  329303 75747 0 0 0 0 29 0 0 0 ...
```
Las siguientes líneas están destinadas a la creación de listas vacías para poder hacer el siguiente paso que es el filtrado de los datos a partir de 4 grupos taxonómicos. 

            
# Filtrado

En esta parte del programa se encuentra un loop que se encarga de filtrar los datos dependiendo del grupo taxonómico al que pertenecen **"phyllum"**, **"family"**, **"genus"**, o **"species"**. 
            
Dentro de este mismo loop se filtran los datos para poder obtener un ** Top n** de los datos, es decir, obtener un Top n de especies, generos, familias y filos. Se puede modificar el tamaño del top dependiendo de las necesidades del estudio. El default de este programa es un Top 10. 
