#!/bin/bash
echo "Instalando paquetes de R necesarios..."

# Instalamos httr y jsonlite usando R desde la terminal
Rscript -e '
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr", repos="http://cran.us.r-project.org")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite", repos="http://cran.us.r-project.org")
'
echo "Instalación completada."
