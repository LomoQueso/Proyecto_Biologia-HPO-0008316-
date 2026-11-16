#!/bin/bash
echo "=========================================="
echo " INICIANDO ANÁLISIS - FENOTIPO HP:0008316 "
echo "=========================================="

# Ejecutar obtención de genes
echo ""
echo "[INFO] Ejecutando la obtención de genes desde HPO..."
echo ""
Rscript code/get_hpo_genes.R

# Ejecutar la construcción de red con STRINGdb
echo ""
echo "[INFO] Ejecutando la construcción de red con STRINGdb..."
echo ""
Rscript code/get_string_network.R

# Ejecutar la conversión de la red de STRINGdb a igraph
echo ""
echo "[INFO] Ejecutando la conversión de la red de STRINGdb a igraph..."
echo ""
Rscript code/stringdb_to_igraph.R

echo ""
echo "Pipeline finalizado con éxito."
