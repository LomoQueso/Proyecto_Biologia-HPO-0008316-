#!/bin/bash

mkdir -p code/R_LIBS
export R_LIBS_USER=code/R_LIBS

echo "Instalando paquetes de R necesarios..."

Rscript -e '

# Configurar carpeta personal de librerías de usuario con permisos de escritura
user_lib <- Sys.getenv("R_LIBS_USER")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

# Ya que uso windows puedo tener problemas con el acceso a los servidores, por ello:
Sys.setenv(CURL_SSL_BACKEND = "openssl")
options(download.file.method = "libcurl", repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("curl", quietly = TRUE)) {
  install.packages("curl", dependencies = TRUE)
}

# Instalamos httr, jsonlite, BiocManager y STRINGdb usando R desde la terminal.
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib)
}
if (!requireNamespace("STRINGdb", quietly = TRUE)) {
  BiocManager::install("STRINGdb", update = FALSE, ask = FALSE)
}

# Instalamos igraph para cambiar el formato de la red de STRINGdb
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

# Instalamos clusterProfiler, org.Hs.eg.db y enrichplot para el enriquecimiento funcional
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot", update = FALSE, ask = FALSE)
}

'
echo "Instalación completada"
