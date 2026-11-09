# Cargar la librería necesaria
library(STRINGdb)

# Leer los genes descargados previamente
gene_file <- "results/data/hpo_gene.txt"
if (!file.exists(gene_file)) {
  stop("Error: No se encuentra 'results/data/hpo_gene.txt'. Ejecuta primero get_hpo_genes.R.")
}
genes <- read.table(gene_file, header = TRUE, stringsAsFactors = FALSE)

# Inicializar STRINGdb
dir.create("results/data", recursive = TRUE, showWarnings = FALSE)
string_db <- STRINGdb$new(
  version = "12.5",
  species = 9606,
  score_threshold = 800,
  input_directory = "results/data"
)

# Mapear los genes a identificadores de STRING
string_ids <- string_db$map(genes, "name", removeUnmappedRows = TRUE)

# Genes Perdidos
cat("\n--- CONTROL DE GENES PERDIDOS EN STRING ---\n")
cat("HPO genes originales:", nrow(genes), "\n")
cat("Genes mapeados a STRING:", nrow(string_ids), "\n")
genes_no_mapeados <- setdiff(genes$name, string_ids$name)
if(length(genes_no_mapeados) > 0) {
  cat("Genes no mapeados:", paste(genes_no_mapeados, collapse = ", "), "\n")
} else {
  cat("Genes no mapeados: Ninguno. Todos se han mapeado con éxito.\n")
}
cat("-------------------------------------------\n\n")

# Obtener la red de interacción entre proteínas
interaction_network <- string_db$get_interactions(string_ids$STRING_id)
interaction_network <- unique(interaction_network[, 1:2]) # Eliminar interacciones duplicadas

# Reemplazar identificadores de STRING por nombres oficiales de los genes
interaction_network$from <- string_ids$name[match(interaction_network$from, string_ids$STRING_id)]
interaction_network$to <- string_ids$name[match(interaction_network$to, string_ids$STRING_id)]

# Guardar la red en results/data/string_network.txt
write.table(interaction_network, file = "results/data/string_network.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Resumen de ejecución en la terminal
cat("Total de interacciones encontradas:", nrow(interaction_network), "\n\n")
cat("Primeras 6 interacciones:\n")
print(head(interaction_network))
