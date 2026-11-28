
export R_LIBS_USER=code/R_LIBS
mkdir -p ../results/data

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

# Ejecutar la visualización de la red
echo ""
echo "[INFO] Ejecutando la visualización de la red..."
echo ""
Rscript code/visualize.R
echo "Archivo con la visualización de la red creado con éxito."

# Ejecutar análisis por compartimento celular
echo ""
echo "[INFO] Ejecutando el análisis por compartimento celular..."
echo ""
Rscript code/cellular_compartment.R
echo "Análisis por compartimento celular completado con éxito."

# Ejecutar el clustering de la red
echo ""
echo "[INFO] Ejecutando el clustering de la red con Louvain..."
echo ""
Rscript code/clustering.R
echo "Ejecución de la red completada con éxito"

# Ejecutar el enriquecimiento funcional
echo ""
echo "[INFO] Ejecutando el análisis de enriquecimiento funcional..."
echo ""
Rscript code/functional_enrichment.R
echo "Ejecución del análisis de enriquecimiento funcional completado con éxito" 

echo ""
echo "Pipeline finalizado con éxito."
