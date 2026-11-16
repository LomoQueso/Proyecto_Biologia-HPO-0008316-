# Cargar las librerías necesarias
suppressPackageStartupMessages(library(igraph))

# Leer el archivo 'string_network.txt' para posteriormente construir la red
interaction_network_file <- "results/data/string_network.txt"
if (!file.exists(interaction_network_file)) {
  stop("Error: No se encuentra 'results/data/string_network.txt'. Ejecuta primero get_string_network.R.")
}

# Una vez detectado el archivo, se convierte en una tabla para que pueda entenderlo R
interaction_network <- read.table(interaction_network_file, header = TRUE)

# Leer el archivo 'string_ids.txt'
string_ids_file <- "results/data/string_ids.txt"
if (!file.exists(string_ids_file)) {
  stop("Error: No se encuentra 'results/data/string_ids.txt'. Ejecuta primero get_string_network.R.")
}

# Transformar a tabla 'string_ids.txt' para R
string_ids <- read.table(string_ids_file, header = TRUE)

# Usamos 'graph_from_data_frame()' para construir la red a partir de un data frame
g <- graph_from_data_frame(interaction_network, directed = FALSE)
cat("Red construida con éxito \n")
cat("\n")

# Genes mapeados a STRING pero sin interacciones (no están en la red)
no_interactions <- setdiff(string_ids$name, V(g)$name)
cat("Genes sin interacciones (no están en la red):", length(no_interactions), "\n")

cat("\n--- PROPIEDADES DE LA RED ---\n")
cat("\n")

# Número de nodos en la red
num_nodes <- gorder(g)
cat("Número de nodos:", num_nodes, "\n")

# Número de aristas en la red
num_edges <- gsize(g)
cat("Número de aristas:", num_edges, "\n")

# Grado medio de la red
degree_values <- degree(g)
average_degree <- mean(degree_values)
cat("Grado medio:", average_degree, "\n")

# Densidad de la red
network_density <- edge_density(g)
cat("Densidad:", network_density, "\n")

# Dispersión de la red
sparsity <- 1 - edge_density(g)
cat("Dispersión:", sparsity, "\n")

# Centralidad de cercanía en la red
centrality_closeness <- closeness(g)
cat("Centralidad de cercanía (5 primeros nodos):", head(centrality_closeness, 5), "\n")

# Modularidad de la red
modularidad <- modularity(cluster_louvain(g))
cat("Modularidad:", modularidad, "\n")

# Longitud promedio de la red
avg_distance <- mean_distance(g)
cat("Longitud Promedio:", avg_distance, "\n")

# Diámetro de la red
diametro <- diameter(g)
cat("Diámetro:", diametro, "\n")

# Grado de asortatividad de la red
assortativity <- assortativity_degree(g)
cat("Grado de asortatividad:", assortativity, "\n")

