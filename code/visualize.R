# Cargar las librerías necesarias
suppressPackageStartupMessages(library(igraph))

# Leer el archivo 'string_network.txt' para posterior visualización
interaction_network_file <- "results/data/string_network.txt"
if (!file.exists(interaction_network_file)) {
  stop("Error: No se encuentra 'results/data/string_network.txt'. Ejecuta primero get_string_network.R.")
}

# Pasamos el archivo a una tabla de R
interaction_network <- read.table(interaction_network_file, header = TRUE)

# Conversión a dataframe
g <- graph_from_data_frame(interaction_network, directed = FALSE)

# Guardamos el archivo en la carpeta results/data
pdf("results/Basic_interaction_network.pdf", width = 10, height = 10)

# Visualización básica
plot(g, vertex.size = 5, vertex.label.cex = 0.7, main = "Basic Interaction Network")


