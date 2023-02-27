# Metagenómica 
Proyecto desarrollado por alumnos y docentes del Laboratorio Nacional de Secuenciación Genomica del Instituto Tecnológico y de Estudios superiores de Monterrey.

Este programa fue diseñado para hacer un análisis metagenómico de 16S y shotgun, a partir de datos en formato TXT o CSV provenientes de un análisis taxonómico de [KRAKEN 2](https://github.com/DerrickWood/kraken2.) 

Este programa busca hacer análisis estadísticos y permitir al usuario la visualización de los datos taxonómicos arrojados por KRAKEN 2. 
# Utils y librerías 
La primera parte del código está diseñada para declarar las diferentes funciones que serán utilizadas para la creación de barplots y boxplots que formarán parte de los resultados de este análisis.

Es importante correr estas lineas del código ya que las funciones no se encuentran en ninguna librería y son necesarias para el funcionamiento del programa.

A su vez, también se instalan (en caso de no tenerlas) y se llama a las librerías que serán utilizadas para hacer los diferentes análisis estadísticos.  
# Datos
En la segunda sección del código se declara el pathway de la carpeta donde se encuentran los datos con los cuales va a trabajar el programa. También, dentro de esta carpeta se van a guardar los diferentes documentos e imágenes resultantes de este análisis. 

Igual dentro de esta sección se hace un tratado de los datos para asegurar que el formato es el correcto y R o Rstudio no arroje códigos de error y códigos de advertencia 

En la línea `raw_data <- read.table(file.choose(), header = T, sep = "\t",quote = "\"", stringsAsFactors = F, fill = F)` se debe de modificar el parámetro "sep" dependiendo del tipo de documento que se esté usando: para **.txt** se usa **"\t"** y para **.csv** se usa **","**.


# Filtrado

En esta parte del programa se encuentra un 
