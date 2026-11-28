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

# 3. Consultar localización cromosómica en org.Hs.eg.db para separar ADN mitocondrial (CHR == "MT") vs nuclear
chr_info <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = V(g)$name,
  keytype = "SYMBOL",
  columns = "CHR"
))
mt_dna_genes <- unique(chr_info$SYMBOL[!is.na(chr_info$CHR) & chr_info$CHR == "MT"])

# Clasificar los nodos de la red según su compartimento celular y origen genómico
V(g)$is_mitochondrial <- V(g)$name %in% mito_genes
V(g)$is_mt_dna        <- V(g)$name %in% mt_dna_genes
V(g)$chromosome       <- chr_info$CHR[match(V(g)$name, chr_info$SYMBOL)]

V(g)$compartment      <- ifelse(V(g)$is_mitochondrial, "Mitocondria", "Otro compartimento")
V(g)$subgroup    <- ifelse(
  !V(g)$is_mitochondrial,
  "Otro compartimento",
  ifelse(V(g)$is_mt_dna, "Mitocondria (ADN mitocondrial)", "Mitocondria (ADN nuclear)")
)

# Conteos por grupo
n_mt_dna  <- sum(V(g)$is_mitochondrial & V(g)$is_mt_dna)
n_nuc_dna <- sum(V(g)$is_mitochondrial & !V(g)$is_mt_dna)
n_other   <- sum(!V(g)$is_mitochondrial)

# 4. Resumen por pantalla y guardado de tabla en results/data/
cat("\n--- ANÁLISIS POR COMPARTIMENTO CELULAR Y ORIGEN GENÓMICO (GO:0005739) ---\n")
cat("Nodos totales en la red:", gorder(g), "\n")
cat("  - Mitocondria (ADN mitocondrial, CHR = MT):", n_mt_dna, "\n")
cat("  - Mitocondria (ADN nuclear):", n_nuc_dna, "\n")
cat("  - Otro compartimento:", n_other, "\n")

non_mito_net <- V(g)$name[!V(g)$is_mitochondrial]
if (length(non_mito_net) > 0) {
  cat("Genes fuera de mitocondria:", paste(non_mito_net, collapse = ", "), "\n")
}

dir.create("results/data", recursive = TRUE, showWarnings = FALSE)
compartment_df <- data.frame(
  gene        = V(g)$name,
  compartment = V(g)$compartment,
  subgroup    = V(g)$subgroup,
  chromosome  = V(g)$chromosome,
  grado       = degree(g),
  stringsAsFactors = FALSE
)
write.table(
  compartment_df,
  file      = "results/data/cellular_compartment.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# 5. Visualización destacando los genes mitocondriales(Rojo/Salmón) sobre el resto
V(g)$color       <- ifelse(!V(g)$is_mitochondrial, "#4DBBD5", ifelse(V(g)$is_mt_dna, "#D73027", "#F4A582"))
V(g)$size        <- 8
V(g)$label.color <- "black"
V(g)$label.font  <- 2

pdf("results/Network_Mitochondrial_Compartment.pdf", width = 11, height = 10)
set.seed(123)
layout_fr <- layout_with_fr(g)

plot(
  g,
  layout           = layout_fr,
  vertex.label.cex = 0.75,
  edge.color       = "gray80",
  main             = "Red de Interacción por Compartimento Celular y Origen Genómico (HP:0008316)"
)

legend(
  "topleft",
  legend = c(
    paste0("Mitocondria - ADN mitocondrial (n = ", n_mt_dna, ")"),
    paste0("Mitocondria - ADN nuclear (n = ", n_nuc_dna, ")"),
    paste0("Otro compartimento (n = ", n_other, ")")
  ),
  col    = c("#D73027", "#F4A582", "#4DBBD5"),
  pch    = 19,
  pt.cex = 1.8,
  bty    = "n",
  cex    = 1.05
)

cat("Gráfico guardado en 'results/Network_Mitochondrial_Compartment.pdf' y tabla en 'results/data/cellular_compartment.txt'\n")
