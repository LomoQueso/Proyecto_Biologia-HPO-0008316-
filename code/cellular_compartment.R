# Cargar las librerías necesarias (ya instaladas en setup.sh)
suppressPackageStartupMessages({
  library(igraph)
  library(org.Hs.eg.db)
})

# 1. Comprobar y leer la red de interacciones y los genes mapeados
network_file <- "results/data/string_network.txt"
ids_file     <- "results/data/string_ids.txt"

if (!file.exists(network_file) || !file.exists(ids_file)) {
  stop("Error: Faltan archivos en 'results/data/'. Ejecuta primero get_string_network.R.")
}

interaction_network <- read.table(network_file, header = TRUE, stringsAsFactors = FALSE)
string_ids          <- read.table(ids_file, header = TRUE, stringsAsFactors = FALSE)

# Construir el objeto igraph incluyendo TODOS los genes mapeados (también nodos aislados como MYH7)
vertices_df <- data.frame(name = unique(string_ids$name), stringsAsFactors = FALSE)
g <- graph_from_data_frame(interaction_network, directed = FALSE, vertices = vertices_df)

# 2. Consultar anotaciones de Componente Celular (GO:0005739 = Mitochondrion)
# GO:0005739 es el identificador estándar del Gene Ontology Consortium (extraído de org.Hs.eg.db)
# para el término 'mitochondrion'. Se emplea como referencia única y sin ambigüedades para
# clasificar todos los genes cuyas proteínas operan en este orgánulo (incluyendo tanto los
# codificados en el ADNmt 'MT-' como los de origen nuclear), mientras que la columna GOALL
# engloba automáticamente todos sus subcompartimentos (membrana interna, matriz, cadena respiratoria, etc.).
go_cc <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(c(string_ids$name, V(g)$name)),
  keytype = "SYMBOL",
  columns = c("GOALL", "ONTOLOGYALL")
))

mito_genes <- unique(go_cc$SYMBOL[!is.na(go_cc$GOALL) & go_cc$GOALL == "GO:0005739"])

# Clasificar los nodos de la red según su compartimento celular
V(g)$is_mitochondrial <- V(g)$name %in% mito_genes
V(g)$compartment      <- ifelse(V(g)$is_mitochondrial, "Mitocondria", "Otro compartimento")

# 3. Resumen por pantalla y guardado de tabla en results/data/
cat("\n--- ANÁLISIS POR COMPARTIMENTO CELULAR (GO:0005739) ---\n")
cat("Nodos totales en la red:", gorder(g), "\n")
cat("Genes mitocondriales en la red:", sum(V(g)$is_mitochondrial), "\n")
cat("Genes NO mitocondriales en la red:", sum(!V(g)$is_mitochondrial), "\n")

non_mito_net <- V(g)$name[!V(g)$is_mitochondrial]
if (length(non_mito_net) > 0) {
  cat("Genes fuera de mitocondria:", paste(non_mito_net, collapse = ", "), "\n")
}

dir.create("results/data", recursive = TRUE, showWarnings = FALSE)
compartment_df <- data.frame(
  gene        = V(g)$name,
  compartment = V(g)$compartment,
  stringsAsFactors = FALSE
)
write.table(
  compartment_df,
  file      = "results/data/cellular_compartment.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# 4. Visualización destacando los genes mitocondriales sobre el resto
V(g)$color       <- ifelse(V(g)$is_mitochondrial, "#E64B35", "#4DBBD5")
V(g)$size        <- ifelse(V(g)$is_mitochondrial, 8, 6)
V(g)$label.color <- ifelse(V(g)$is_mitochondrial, "black", "gray30")
V(g)$label.font  <- ifelse(V(g)$is_mitochondrial, 2, 1)

pdf("results/Network_Mitochondrial_Compartment.pdf", width = 11, height = 10)
set.seed(123)
layout_fr <- layout_with_fr(g)

plot(
  g,
  layout           = layout_fr,
  vertex.label.cex = 0.75,
  edge.color       = "gray80",
  main             = "Red de Interacción por Compartimento Celular (HP:0008316)"
)

legend(
  "topleft",
  legend = c(
    paste0("Mitocondria (n = ", sum(V(g)$is_mitochondrial), ")"),
    paste0("Otro compartimento (n = ", sum(!V(g)$is_mitochondrial), ")")
  ),
  col    = c("#E64B35", "#4DBBD5"),
  pch    = 19,
  pt.cex = 1.8,
  bty    = "n",
  cex    = 1.05
)

cat("Gráfico guardado en 'results/Network_Mitochondrial_Compartment.pdf' y tabla en 'results/data/cellular_compartment.txt'\n")
