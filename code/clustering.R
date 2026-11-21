# Cargar la librería necesaria
suppressPackageStartupMessages(library(igraph))

# Leer el archivo 'string_network.txt' generado previamente
interaction_network_file <- "results/data/string_network.txt"
if (!file.exists(interaction_network_file)) {
  stop("Error: No se encuentra 'results/data/string_network.txt'. Ejecuta primero get_string_network.R.")
}

interaction_network <- read.table(interaction_network_file, header = TRUE)

# Construir la red en formato igraph
g <- graph_from_data_frame(interaction_network, directed = FALSE)

# Clustering de la red con el método de Louvain (fijando semilla)
set.seed(123)
communities <- cluster_louvain(g)

# Mostrar el número de clusters y sus tamaños en la terminal
cat("\n--- CLUSTERING DE LA RED (LOUVAIN) ---\n")
cat("Number of clusters:", length(communities), "\n")
cat("Community sizes:\n")
print(sizes(communities))

# Guardar la asignación de cada gen a su cluster en results/data/
clusters_df <- data.frame(
  gene = V(g)$name,
  cluster = membership(communities)
)
write.table(clusters_df, file = "results/data/louvain_clusters.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Crear la visualización del clustering
pdf("results/Clustering_Louvain.pdf", width = 10, height = 10)

plot(communities, g, vertex.size = 5, vertex.label = NA, main = "Clustering with Louvain")

cat("Gráfico guardado en 'results/Clustering_Louvain.pdf' y tabla en 'results/data/louvain_clusters.txt'\n")
