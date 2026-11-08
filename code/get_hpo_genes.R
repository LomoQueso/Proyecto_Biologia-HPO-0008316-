# Cargar librerías necesarias
library(httr)
library(jsonlite)

# Definir el fenotipo HPO asignado al grupo
hpo_id <- "HP:0008316" # Abnormal Mitochondria in muscle tissue
url <- paste0("https://ontology.jax.org/api/network/annotation/", hpo_id)

# Consultar la API de la HPO
response <- GET(url)
data <- content(response, as = "text", encoding = "UTF-8")
genes <- fromJSON(data)$genes

# Número total de genes
cat("Total de genes obtenidos:", nrow(genes), "\n\n")
# Primeros genes
cat("Primeros genes:\n")
print(head(genes))
