#!/bin/bash
echo "=========================================="
echo " INICIANDO ANÁLISIS - FENOTIPO HP:0008316 "
echo "=========================================="

# Ejecutar obtención de genes
echo ""
echo "[INFO] Ejecutando la obtención de genes desde HPO..."
Rscript code/get_hpo_genes.R

# Ejecutar la construcción de red con STRINGdb
echo ""
echo "[INFO] Ejecutando la construcción de red con STRINGdb..."
Rscript code/get_string_network.R

echo "Pipeline finalizado con éxito."
