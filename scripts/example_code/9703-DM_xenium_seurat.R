library(Seurat)
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)

path <- "/nfs/mm-isilon/bioinfcore/labshare/bfxcore/damki_lab/20231116__185301__9703-DM/output-XETG00077__0006740__9301-DM-1_ROI_B__20231116__185322/"
# Load the Xenium data
# ?? FOV in morphology_fov_locations.json
# ?? Negative control codewords vs Unassigned codewords vs Deprecated codewords
# ?? decoding assigns a unique codeword for each transcript? (deduplicate transcripts?)
xenium.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

# normalize and cluster
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.1)
DimPlot(xenium.obj)
ImageDimPlot(xenium.obj, cols = "polychrome", axes=T)

# new fov by cropping
crop <- Crop(xenium.obj[["fov"]], x = c(2500, 3000), y = c(1000, 1600))
xenium.obj[["crop"]] <- crop
DefaultBoundary(xenium.obj[["crop"]]) <- "segmentation"
ImageDimPlot(xenium.obj, fov='crop', cols = "polychrome", axes=T, molecules=c("MYH11","ACTA2"))


ImageDimPlot(xenium.obj, fov = "fov", molecules = c("MYH11","ACTA2"), nmols = 20000)

# molecules vs features
ImageFeaturePlot(xenium.obj, features = c("MYH11","ACTA2"),molecules = c("MYH11","ACTA2"), axes = TRUE)
ImageFeaturePlot(xenium.obj, features = c("MYH11","ACTA2"), axes = TRUE, size = 0.75, cols = c("white", "red"))


cropped.coords <- Crop(xenium.obj[["fov"]], x = c(1200, 1300), y = c(3750, 4000), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords

# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",coord.fixed = FALSE, molecules = c("MYH11","ACTA2"), nmols = 10000,alpha=0.1)



## annotation

### singleR
aggexp <- AggregateExpression(xenium.obj, group.by = "seurat_clusters", return.seurat = T)@assays$SCT@data
write.csv(aggexp,'xenium.aggexp.csv')

BiocManager::install(c("SingleR", "celldex"))

library(SingleR)
library(celldex)

hpca.se <- HumanPrimaryCellAtlasData()

pred.hpca.se <- SingleR(test = aggexp, ref = hpca.se, labels = hpca.se$label.main)
View(as.data.frame(pred.hpca.se) %>% select(pruned.labels))

## using ref data
# count matrix and annotation data downloaded from: https://singlecell.broadinstitute.org/single_cell/study/SCP1909/aortic-cellular-diversity-and-quantitative-genome-wide-association-study-trait-prioritization-through-single-nuclear-rna-sequencing-of-the-aneurysmal-human-aorta#study-download
ref.expression_matrix <- ReadMtx(
  mtx = "AorticAneurysm_Expression_Matrix_raw_counts_V1.mtx", features = "AorticAneurysm_Expression_Matrix_genes_V1.tsv",
  cells = "AorticAneurysm_Expression_Matrix_barcodes_V1.tsv")
ref.obj <- CreateSeuratObject(counts = ref.expression_matrix)
ref.meta = read.delim('AorticAneurysm_MetaData_V1.txt')[-1,] %>% select(NAME,celltype)
ref.meta$cluster = sapply(ref.meta$celltype, function(x){strsplit(gsub('\\.','',x),' ')[[1]][1]})
ref.meta$cell_type = sapply(ref.meta$celltype, function(x){strsplit(gsub('\\.','',x),' ')[[1]][2]})


ref.obj <- SCTransform(ref.obj, assay = "RNA")
ref.obj <- RunPCA(ref.obj)
ref.obj <- RunUMAP(ref.obj, dims = 1:30)
ref.obj <- FindNeighbors(ref.obj, dims = 1:30)
ref.obj <- FindClusters(ref.obj, resolution = 0.1)

ref.meta = ref.obj@meta.data %>% rownames_to_column('NAME') %>% left_join(ref.meta, by='NAME')
ref.obj@meta.data$cell_type = ref.meta$cell_type
ref.obj@meta.data$ref_clusters = ref.meta$cluster
Idents(ref.obj) = ref.obj@meta.data$ref_clusters

counts <- GetAssayData(ref.obj, assay = "RNA", slot = "counts")
cluster <- as.factor(ref.obj$cell_type)
names(cluster) <- colnames(ref.obj)
nUMI <- ref.obj$nCount_RNA
names(nUMI) <- colnames(ref.obj)
nUMI <- colSums(counts)
levels(cluster) <- levels(cluster)
reference <- Reference(counts, cluster, nUMI)


query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")
coords <- GetTissueCoordinates(xenium.obj, which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")



annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
xenium.obj$predicted.celltype <- annotations
ImageDimPlot(xenium.obj, group.by='predicted.celltype', cols = "polychrome", size = 0.75)
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
xenium.obj.sub <- subset(xenium.obj, cells = keep.cells)

xenium.obj.sub <- BuildNicheAssay(object = xenium.obj.sub, group.by = "predicted.celltype",
    niches.k = 5, neighbors.k = 30, fov='fov')

# plot annotations and niches
celltype.plot <- ImageDimPlot(xenium.obj.sub, group.by = "predicted.celltype", size = 1.5, cols = "polychrome",
    dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj.sub, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
    scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot


saveRDS(xenium.obj,'xenium.obj.Rds')
saveRDS(xenium.obj.sub,'xenium.obj.sub.Rds')
saveRDS(ref.obj, 'AorticAneurysm.ref.obj.Rds')
saveRDS(reference, 'AorticAneurysm.reference.Rds')

