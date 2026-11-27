# Cargar las librerías necesarias
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

# Comprobar y leer los genes descargados de la HPO (Parte 1)
gene_file <- "results/data/hpo_gene.txt"
if (!file.exists(gene_file)) {
  stop("Error: No se encuentra 'results/data/hpo_gene.txt'. Ejecuta primero get_hpo_genes.R.")
}

genes <- read.table(gene_file, header = TRUE, stringsAsFactors = FALSE)

# 1. Enriquecimiento funcional (ORA) de la lista completa de genes HPO

cat("\n--- ENRIQUECIMIENTO FUNCIONAL (GO - BP) LISTA COMPLETA ---\n")

enrich_results <- enrichGO(
  gene          = genes$name,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

print(head(enrich_results))

# Guardar tabla de resultados en CSV
write.csv2(as.data.frame(enrich_results), file = "results/data/GO_enrichment.csv", row.names = FALSE)

# Guardar dotplot y barplot en PDF
if (nrow(as.data.frame(enrich_results)) > 0) {
  pdf("results/GO_Enrichment_Global.pdf", width = 10, height = 8)
  print(dotplot(enrich_results, showCategory = 10, title = "GO Enrichment (BP)"))
  print(barplot(enrich_results, showCategory = 10, title = "GO Enrichment (BP)"))
}


# 2. Enriquecimiento funcional usando como fondo (universe) los genes de la HPO

cat("\n--- ENRIQUECIMIENTO FUNCIONAL CON UNIVERSO HPO ---\n")

g2p <- tryCatch(
  read.delim("https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt"),
  error = function(e) {
    message("Aviso: No se pudo descargar genes_to_phenotype.txt: ", e$message)
    NULL
  }
)

if (!is.null(g2p)) {
  enrich_results_hpo <- enrichGO(
    gene     = genes$name,
    universe = unique(g2p$gene_symbol),
    OrgDb    = org.Hs.eg.db,
    keyType  = "SYMBOL",
    ont      = "BP"
  )

  print(head(enrich_results_hpo))
  write.csv2(as.data.frame(enrich_results_hpo), file = "results/data/GO_enrichment_hpo_universe.csv", row.names = FALSE)
}

# 3. Enriquecimiento ORA de cada Cluster de Louvain (tamaño >= 5)

cat("\n--- ENRIQUECIMIENTO FUNCIONAL POR CLUSTER (LOUVAIN) ---\n")

clusters_file <- "results/data/louvain_clusters.txt"
if (!file.exists(clusters_file)) {
  stop("Error: No se encuentra 'results/data/louvain_clusters.txt'. Ejecuta primero clustering.R.")
}

clusters_df <- read.table(clusters_file, header = TRUE, stringsAsFactors = FALSE)

# Agrupar genes por cluster y filtrar clusters con al menos 5 genes
genes_clusters <- split(clusters_df$gene, clusters_df$cluster)

min_size <- 5
genes_clusters <- genes_clusters[lengths(genes_clusters) >= min_size]

cat("Clusters with at least", min_size, "genes:", length(genes_clusters), "\n")
print(lengths(genes_clusters))

# Enriquecimiento ORA de cada cluster usando la red completa como fondo (universe)
enrichment_results_clusters <- list()

for (i in names(genes_clusters)) {
  enrichment_results_clusters[[i]] <- enrichGO(
    gene          = genes_clusters[[i]],
    universe      = clusters_df$gene,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
}

# Término más significativo para cada cluster
cluster_summary <- data.frame(
  cluster  = names(genes_clusters),
  genes    = lengths(genes_clusters),
  top_term = sapply(enrichment_results_clusters, function(x) {
    res <- as.data.frame(x)
    if (nrow(res) == 0) "(none)" else res$Description[1]
  })
)
print(cluster_summary)
write.csv2(cluster_summary, file = "results/data/enrichment_clusters_summary.csv", row.names = FALSE)

# Guardar los resultados de cada cluster en un CSV (exacto al manual)
for (i in names(enrichment_results_clusters)) {
  write.csv2(
    as.data.frame(enrichment_results_clusters[[i]]),
    file = paste0("results/data/enrichment_cluster_", i, ".csv"),
    row.names = FALSE
  )
}

# Guardar dotplot de todos los clusters que tengan términos significativos
clusters_with_terms <- names(enrichment_results_clusters)[
  sapply(enrichment_results_clusters, function(x) nrow(as.data.frame(x)) > 0)
]

if (length(clusters_with_terms) > 0) {
  pdf("results/GO_Enrichment_Clusters.pdf", width = 10, height = 8)
  for (i in clusters_with_terms) {
    print(dotplot(
      enrichment_results_clusters[[i]],
      showCategory = 10,
      title = paste("GO Enrichment, cluster", i)
    ))
  }

}


# 4. Refinamiento funcional (simplificación semántica y análisis multiontología BP, MF y CC)

cat("\n--- REFINAMIENTO FUNCIONAL: SIMPLIFICACIÓN GO-BP Y MULTIONTOLOGÍA (BP, MF, CC) ---\n")

if (!is.null(g2p)) {
  # 4.1. Eliminación de redundancia jerárquica en GO-BP (universo HPO)
  enrich_bp_simp <- simplify(enrich_results_hpo, cutoff = 0.7, by = "p.adjust", select_fun = min)
  write.csv2(as.data.frame(enrich_bp_simp), file = "results/data/GO_enrichment_hpo_simplified.csv", row.names = FALSE)

  # 4.2. Enriquecimiento en Función Molecular (MF) y Componente Celular (CC) con universo HPO
  enrich_mf <- enrichGO(
    gene          = genes$name,
    universe      = unique(g2p$gene_symbol),
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  enrich_mf_simp <- simplify(enrich_mf, cutoff = 0.7, by = "p.adjust", select_fun = min)
  write.csv2(as.data.frame(enrich_mf_simp), file = "results/data/GO_MF_enrichment_hpo.csv", row.names = FALSE)

  enrich_cc <- enrichGO(
    gene          = genes$name,
    universe      = unique(g2p$gene_symbol),
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  enrich_cc_simp <- simplify(enrich_cc, cutoff = 0.7, by = "p.adjust", select_fun = min)
  write.csv2(as.data.frame(enrich_cc_simp), file = "results/data/GO_CC_enrichment_hpo.csv", row.names = FALSE)

  # 4.3. Exportación de gráficos (dotplot y barplot) en PDF
  pdf("results/GO_Enrichment_Refined.pdf", width = 10, height = 8)

  if (nrow(as.data.frame(enrich_bp_simp)) > 0) {
    print(dotplot(enrich_bp_simp, showCategory = 10, title = "GO-BP Simplificado (Universo HPO)"))
    print(barplot(enrich_bp_simp, showCategory = 10, title = "GO-BP Simplificado (Universo HPO)"))
  }

  if (nrow(as.data.frame(enrich_mf_simp)) > 0) {
    print(dotplot(enrich_mf_simp, showCategory = 10, title = "GO-MF Simplificado (Universo HPO)"))
  }

  if (nrow(as.data.frame(enrich_cc_simp)) > 0) {
    print(dotplot(enrich_cc_simp, showCategory = 10, title = "GO-CC Simplificado (Universo HPO)"))
  }

 
}

cat("\nEnriquecimiento funcional completado con éxito.\n")
