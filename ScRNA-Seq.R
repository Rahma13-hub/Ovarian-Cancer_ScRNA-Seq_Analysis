#needed packages
library(tidyverse)
library(Seurat)
library(hdf5r)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(enrichR)
#Load the scRNA-seq dataset from 10X Genomics HDF5 file
data = Read10X_h5(filename = "C:/Users/3nany/Downloads/17k_Ovarian_Cancer_scFFPE_count_filtered_feature_bc_matrix.h5")
str(data)
# Create Seurat object for downstream analysis, filtering genes detected in at least 3 cells
# and cells expressing at least 200 genes
data.seurat.obj = CreateSeuratObject(counts = data, project = "OvarianCancer",
                                    min.cells = 3, min.features = 200)
str(data.seurat.obj)
#Quality Control - Calculate percent of mitochondrial genes to assess cell quality
View(data.seurat.obj@meta.data)
data.seurat.obj[["percent.mt"]] = PercentageFeatureSet(data.seurat.obj, pattern =  "^MT-")
# Visualize QC metrics
VlnPlot(data.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Plot relationship between total counts and detected features per cell
FeatureScatter(data.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")
# Filter cells based on QC thresholds to exclude low quality or dying cells
data.seurat.obj = subset(data.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                           percent.mt < 5)
# Normalize the data to correct for sequencing depth differences across cells
data.seurat.obj = NormalizeData(data.seurat.obj)
#Identify highly variable genes for dimensionality reduction and clustering
data.seurat.obj = FindVariableFeatures(data.seurat.obj, selection.method = "vst", nfeatures = 2000)
# Visualize top 10 variable genes
top10 = head(VariableFeatures(data.seurat.obj), 10)
plot1 = VariableFeaturePlot(data.seurat.obj)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
#scale the data
data.seurat.obj = ScaleData(data.seurat.obj, features = VariableFeatures(data.seurat.obj))
#Perform PCA for dimensionality reduction
data.seurat.obj = RunPCA(data.seurat.obj, features = VariableFeatures(object = data.seurat.obj))
# # Visualize PCA results for first 5 components
print(data.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(data.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
# determine dimensionality of the data
ElbowPlot(data.seurat.obj)
#Construct nearest neighbor graph based on selected PCs and perform clustering
data.seurat.obj = FindNeighbors(data.seurat.obj, dims = 1:15)
data.seurat.obj = FindClusters(data.seurat.obj, resolution = 0.5)
view(data.seurat.obj@meta.data)
# Run UMAP for nonlinear dimensionality reduction and visualize clusters
data.seurat.obj = RunUMAP(data.seurat.obj, dims = 1:15)
DimPlot(data.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
# setting identity of clusters
Idents(data.seurat.obj)
Idents(data.seurat.obj) = "RNA_snn_res.0.5"
# Identify marker genes for each cluster with thresholds on minimum expression and log fold change
markers = FindAllMarkers(data.seurat.obj, only.pos = TRUE, min.pct = 0.25,
                         logfc.threshold = 0.25)
head(markers)
# Extract top 10 markers per cluster based on average log2 fold change
top10.markers = markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
# Print top 10 marker genes per cluster
print(top10.markers)
view(top10.markers)
# Extract scaled expression data from Seurat object for marker genes
expr_data = GetAssayData(data.seurat.obj, assay = "RNA", slot = "scale.data")
# Prepare gene list for heatmap plotting by filtering genes present in expression matrix
genes_to_plot = unique(top10.markers$gene)
genes_to_plot = genes_to_plot[genes_to_plot %in% rownames(expr_data)]
# Subset expression data to marker genes only
expr_subset = expr_data[genes_to_plot, ]
# Get cluster identities for ordering cells
clusters = Idents(data.seurat.obj)
# Order cells by cluster identity for heatmap visualization
ordered_cells = names(sort(clusters[colnames(expr_subset)]))
expr_subset = expr_subset[, ordered_cells]
clusters = clusters[ordered_cells]
# Create annotation for heatmap columns with cluster identities
cluster_annotation = HeatmapAnnotation(Cluster = clusters)
# Plot heatmap of scaled expression data for marker genes across clusters
Heatmap(expr_subset,
        name = "Expression",
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = FALSE,    
        top_annotation = cluster_annotation,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 4))
# Rename cluster identities to biologically meaningful cell types based on marker analysis
data.seurat.obj = RenameIdents(
  data.seurat.obj,
  "0" = "Cancer-associated fibroblasts",
  "1" = "Ciliated epithelial cells",
  "2" = "Proliferating tumor cells",
  "3" = "Antigen-presenting myeloid cells",
  "4" = " inflammatory stromal cells",
  "5" = "Endothelial cells",
  "6" = "Smooth muscle cells",
  "7" = "Epithelial (tumor) cells",
  "8" = "Tumor-associated macrophages",
  "9" = "Activated cytotoxic T cells",
  "10" = "Ovarian tumor epithelial cells",
  "11" = "Pericytes / vascular smooth muscle cells",
  "12" = "Basal epithelial cells",
  "13" = "Neuroendocrine-like epithelial cells",
  "14" = "Interferon-responsive immune cells",
  "15" = "Multi-ciliated epithelial cells"
)
# Plot UMAP embedding with clusters labeled by cell type
DimPlot(data.seurat.obj, reduction = "umap", label = T, pt.size = 0.5) +
  NoLegend()
# Save the renamed Seurat object for future use or downstream analysis
saveRDS(data.seurat.obj, file = "data_seurat_renamed.rds")
# Specify databases for enrichment analysis
dbs = c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")
# Split marker genes by cluster to perform enrichment per cluster
cluster_genes = split(top10.markers$gene, top10.markers$cluster)
# Initialize list to store enrichment results
results_list = list()
# Perform enrichment analysis for each cluster
for (cl in names(cluster_genes)) {
  genes = cluster_genes[[cl]]
  enriched = enrichr(genes, dbs)
  results_list[[cl]] = enriched
}
# Export enrichment results as CSV file
for (cl in names(results_list)) {
  for (db in dbs) {
    file_name = paste0("Cluster_", cl, "_", db, "_Enrichment.csv")
    write.csv(results_list[[cl]][[db]], file_name, row.names = FALSE)
  }
}
# Export top 10 markers results as CSV file
write.csv(top10.markers, file = "top10_markers_all_clusters.csv", row.names = FALSE)

