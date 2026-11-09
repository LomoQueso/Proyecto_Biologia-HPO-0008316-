#!/bin/bash
echo "Instalando paquetes de R necesarios..."

Rscript -e '

# Ya que uso windows puedo tener problemas con el acceso a los servidores, por ello:
Sys.setenv(CURL_SSL_BACKEND = "openssl")
install.packages("curl", force = TRUE, dependencies = TRUE)

# Instalamos httr, jsonlite, BiocManager y STRINGdb usando R desde la terminal.
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("STRINGdb", quietly = TRUE)) {
  BiocManager::install("STRINGdb", update = FALSE, ask = FALSE)
}
'
echo "Instalación completada."
