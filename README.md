## Functional Enrichment Analysis of Ovarian Cancer scRNA-seq Data
#Project Overview
This project analyzed single-cell RNA sequencing (scRNA-seq) data from ovarian cancer tissue samples. The aim was to identify cellular heterogeneity by clustering cells and to explore functional enrichment of each cluster through Gene Ontology, KEGG pathways, and Reactome databases. This analysis provides insights into the biological roles of diverse cell populations within the tumor microenvironment and highlights mechanisms relevant to ovarian cancer progression, immune response, and therapy resistance.
#Dataset and Analysis Pipeline
Dataset: scRNA-seq data from ovarian cancer samples (filtered and quality controlled).
Tools: Seurat (clustering, dimensionality reduction, visualization), enrichR (functional enrichment analysis).
Clustering: Identified 16 clusters (0 to 15) representing distinct cell types within ovarian cancer tissue.
Annotation: Clusters renamed based on marker genes.
##Visualization
Below are key visualizations representing the main steps and findings from the single-cell RNA-seq analysis of the ovarian cancer sample:
1. UMAP Plot of Cell Clusters: shows the clustering of single cells into biologically meaningful groups identified in the ovarian cancer tissue. Each cluster is labeled with its inferred cell type.
2. Quality Control (QC) Plot: Violin plot showing the distribution of number of detected genes per cell (nFeature_RNA), total counts per cell (nCount_RNA), and percentage of mitochondrial gene expression (%MT).
3. Highly Variable Genes Plot: This figure displays the top highly variable genes used for downstream dimensionality reduction and clustering.
4. PCA Elbow Plot: The elbow plot guides the selection of the number of principal components to include in the clustering analysis.
5. Complex Heatmap of Top Marker Genes per Cluster: Heatmap illustrating the scaled expression of the top 10 marker genes in each cluster. 


##Cluster-wise Functional Interpretation
#Cluster 0 – Cancer-associated Fibroblasts:
-Involved in extracellular matrix remodeling, immune regulation, and tumor stroma formation.
-Cancer Relevance: Supports tumor growth and invasion by modifying the tumor microenvironment and mediating interactions with immune cells.
#Cluster 1 – Ciliated Epithelial Cells
-Cells specialized in fluid movement and mucociliary clearance.
-Cancer Relevance: Reflects presence of normal or reactive epithelial cells within ovarian tissue, potentially affected by tumor-induced changes.
#Cluster 2 – Proliferating Tumor Cells
-Actively dividing cancer cells with high cell cycle and DNA replication activity
-Cancer Relevance: Represents aggressive tumor cell populations driving tumor growth.
#Cluster 3 – Antigen-presenting Myeloid Cells
-Immune cells presenting antigens to T cells, activating immune responses.
-Cancer Relevance: Reflects tumor-infiltrating immune cells involved in immune surveillance or immune escape.
#Cluster 4 – Inflammatory Stromal Cells
-Cells involved in inflammation and cytokine signaling within the tumor stroma.
-Cancer Relevance: Contributes to inflammatory tumor microenvironment, potentially promoting cancer progression or immune suppression.
#Cluster 5 – Endothelial Cells
-Cells forming blood vessels, supporting angiogenesis.
-Cancer Relevance: Supports tumor vascularization, essential for nutrient supply and metastasis.
#Cluster 6 – Smooth Muscle Cells
-Contractile cells involved in tissue structure and vessel support.
-Cancer Relevance: Present in tumor-associated vasculature and stromal remodeling.
#Cluster 7 – Epithelial (Tumor) Cells with Stress Response
-Tumor cells with active heat shock response, ER stress adaptation, and senescence regulation.
-Cancer Relevance: Enables tumor cell survival under proteotoxic stress and contributes to therapy resistance.
#Cluster 8 – Tumor-associated Macrophages and Nutrient Metabolism
-Macrophages modulating immune responses and remodeling extracellular matrix.
-Cancer Relevance: Supports immune evasion, tumor invasion, and metabolic adaptation.
#Cluster 9 – Activated Cytotoxic T Cells and Immune Regulation.
-Effector T cells and NK cells mediating anti-tumor cytotoxicity and immune regulation.
-Cancer Relevance: Central in tumor immune surveillance and immune escape mechanisms.
#Cluster 10 – Ovarian Tumor Epithelial Cells with Developmental and Hormonal Signatures
-Epithelial tumor cells expressing developmental transcription factors and hormonal regulation genes.
-Cancer Relevance: Reflects tumor cells maintaining epithelial identity and remodeling immune and metabolic pathways for survival.
#Cluster 11 – Pericytes / Vascular Smooth Muscle Cells with Apoptosis and ECM Remodeling
-Cells involved in vessel stabilization, apoptosis regulation, and extracellular matrix metabolism.
-Cancer Relevance: Supports angiogenesis, tumor-stroma interactions, and survival under therapy stress.
#Cluster 12 – Basal Epithelial Cells with Barrier and Signaling Functions
-Epithelial basal cells involved in skin and epithelial development, water transport, and cell signaling.
-Cancer Relevance: Linked to epithelial barrier integrity, signaling pathways that contribute to tumor growth and immune interactions.
#Cluster 13 – Neuroendocrine-like Epithelial Cells and GPCR Signaling
-Cells expressing prostanoid and retinoic acid receptors, involved in cell signaling and tissue homeostasis.
-Cancer Relevance: May influence tumor microenvironment communication, cell adhesion, and metabolic regulation.
#Cluster 14 – Interferon-responsive Immune Cells
-Immune cells with active antiviral and interferon-stimulated gene expression.
-Cancer Relevance: Reflects innate immune activation against viral components and tumor immune modulation.
#Cluster 15 – Multi-ciliated Epithelial Cells and Immune Regulation
-Epithelial cells with cilia involved in immune signaling and tissue integrity.
-Cancer Relevance: Supports immune cell proliferation and tissue barrier functions, possibly influencing tumor microenvironment stability.
##Conclusion
This comprehensive scRNA-seq analysis of ovarian cancer tissue reveals a complex cellular ecosystem involving tumor cells, immune populations, stromal cells, and vascular components. Functional enrichment highlights diverse biological processes including proliferation, immune regulation, metabolic adaptation, extracellular matrix remodeling, and therapy resistance mechanisms. These insights provide a valuable resource for identifying potential therapeutic targets.




