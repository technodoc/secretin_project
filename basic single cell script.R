library(scCATCH)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Pre-Processing
# 1. Expresison Matrix
# 2. Filter Cells /QC
# 3. Normalization

# Clustering
# 4. Identify Variable genes
# 5. Dimensionality Reduction
# 6. Explore Marker Genes

# Biology
# 7. DEG
# 8. Assign Cell Type
# 9. Functional Annotation

# Seurat STEPS
# Read10X
# CreateSeurateObject
# NormalizeData
# FindVariableFeatures
# ScaleData
# RunPCA
# FindNeighbors
# FindClusters
# RunTSNE
# DimPlot

#TODO Build CLI for user argument parsing
#TODO Multiple tools for each step using switch

# Set data directory
# Path to parent directory for tabula muris dataset
data_dir <-
  '/home/technodoc/Desktop/NAS/Tabula_Muris/01_droplet_raw_data/droplet'
# List of subdirectories of parent directory
subdirectories_name_list <-
  list.dirs(path = data_dir,
            full.names = FALSE,
            recursive = FALSE)
number_subdirectories <- length(subdirectories_name_list)

if (number_subdirectories == 1) {
  folder_path <-
    paste(data_dir, subdirectories_name_list[[1]], sep = "/")
  expression_matrix <- Read10X(data.dir = folder_path)
  all_seurat <-
    CreateSeuratObject(
      counts = expression_matrix,
      project = subdirectories_name_list[[1]],
      min.cells = 4,
      min.features = 250
    )
  
} else if (number_subdirectories > 1) {
  # List of Seurat objects, one for each subdirectory
  seurat_object_list <-
    vector(mode = "list", length = number_subdirectories)
  
  for (i in 1:number_subdirectories) {
    # Construct full path to subdirectories
    folder_path <-
      paste(data_dir, subdirectories_name_list[[i]], sep = "/")
    # Read subdirectories
    expression_matrix <- Read10X(data.dir = folder_path)
    # Create Seurat Object for read data and store in Seurat object list.
    seurat_object_list[[i]] <-
      CreateSeuratObject(
        counts = expression_matrix,
        project = subdirectories_name_list[[i]],
        min.cells = 3,
        min.features = 200
      )
  }
  
  # Rename cells using subdirectories names as prefix to enforce unique cell names
  for (i in 1:length(seurat_object_list)) {
    newID <- seurat_object_list[[i]]@project.name
    seurat_object_list[[i]] <-
      RenameCells(seurat_object_list[[i]], add.cell.id = newID)
  }
  # Merge all Seurat Objects in one.
  all_seurat <-
    merge(seurat_object_list[[1]], y = c(seurat_object_list[2:number_subdirectories]))
}


all_seurat[["fraction.mito"]] <-
  PercentageFeatureSet(all_seurat, pattern = "^mt-")

if (sum(all_seurat[["fraction.mito"]] != 0) > 0) {
  VlnPlot(
    all_seurat,
    features = c("nFeature_RNA", "nCount_RNA", "fraction.mito"),
    ncol = 3
  )
  plot_count_feature <-
    FeatureScatter(all_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot_count_mito <-
    FeatureScatter(all_seurat, feature1 = "nCount_RNA", feature2 = "fraction.mito")
  plot_count_feature + plot_count_mito
} else {
  VlnPlot(all_seurat,
          features = c("nFeature_RNA", "nCount_RNA"),
          ncol = 2)
  plot_count_feature <-
    FeatureScatter(all_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot_count_feature
}

# Filter cell using fraction.mito and nFeature_RNA
filter.list <- c("nFeature_RNA", "fraction.mito")
all_seurat <- subset(x = all_seurat,
                     nFeature_RNA > 200 & fraction.mito > -Inf &
                       nFeature_RNA < 2500 & fraction.mito < 0.1)

# Normalize Data
# TODO Multiple normalization methods.
all_seurat <-
  NormalizeData(all_seurat, normalization.method = "LogNormalize")

# Find variable genes/features
all_seurat <-
  FindVariableFeatures(all_seurat,
                       selection.method = "vst",
                       nfeatures = 2000)

# Extract top 10 most variable genes
top10 <- head(VariableFeatures(all_seurat), 10)
# Plot variable features without labels
plot_unlabeled_vf <- VariableFeaturePlot(all_seurat)
# Plot variable features with labels
plot_labeled_vf <- LabelPoints(
  plot = plot_unlabeled_vf,
  points = top10,
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)
plot_labeled_vf


all.genes <- rownames(all_seurat)
# Scale Data
all_seurat <- ScaleData(all_seurat, features = all.genes)

# Run PCA
all_seurat <-
  RunPCA(all_seurat, features = VariableFeatures(object = all_seurat))

print(all_seurat[["pca"]], dims = 1:2, nfeatures = 5)
VizDimLoadings(
  all_seurat,
  dims = 1:2,
  nfeatures = 15,
  reduction = "pca"
)
DimPlot(all_seurat,
        reduction = "pca",
        label = TRUE,
        repel = TRUE)
DimHeatmap(all_seurat,
           dims = 1,
           cells = 500,
           balanced = TRUE)


all_seurat <- JackStraw(all_seurat, num.replicate = 100)
all_seurat <- ScoreJackStraw(all_seurat, dims = 1:20)
JackStrawPlot(all_seurat, dims = 1:15)
ElbowPlot(all_seurat, ndims = 20)

all_seurat <- FindNeighbors(all_seurat, dims = 1:12)
all_seurat <- FindClusters(all_seurat, resolution = 0.5)

all_seurat <- RunUMAP(all_seurat, dims = 1:20)
all_seurat <- RunTSNE(all_seurat, dims = 1:20)

# Calculatge all marker genes and group them by clusters
all_seurat.markers <-
  FindAllMarkers(
    all_seurat,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
x <-
  all_seurat.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)

# TODO Multiple cell-type annotation techniques
# Marker gene database-based cell-type annotation
clu_ann <- scCATCH(
  object = x,
  species = 'Mouse',
  cancer = NULL,
  tissue = c(
    'Kidney',
    'Mesonephros',
    'Lung',
    'Trachea',
    'Bronchiole',
    'Liver',
    'Spleen',
    'Thymus'
  )
)

# TODO Loop for FeaturePlot
FeaturePlot(all_seurat, features = x$gene[1:4])
FeaturePlot(all_seurat, features = x$gene[5:8])
FeaturePlot(all_seurat, features = x$gene[9:12])
FeaturePlot(all_seurat, features = x$gene[13:16])
FeaturePlot(all_seurat, features = x$gene[17:20])

# Rename clusters to maker gene database-based cell-type annotation
new.cluster.ids <- c(clu_ann$cell_type)
names(new.cluster.ids) <- levels(all_seurat)
all_seurat <- RenameIdents(all_seurat, new.cluster.ids)

# Dimensionality reduction plot tSNE
DimPlot(
  all_seurat,
  reduction = "tsne",
  label = TRUE,
  pt.size = 0.5
)

# Dimensionality reduction plot umap
DimPlot(
  all_seurat,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5
)
