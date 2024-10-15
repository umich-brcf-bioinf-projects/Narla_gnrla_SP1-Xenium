## Seurat analysis of Xenium aorta data

## load libraries
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE) # consider using withr option to avoid needing to re-install each time
library(spacexr)
library(Seurat)

library(future)
availableCores()
#plan("multisession", workers = 4) # have BioC already sets
library(tidyverse)
library(svglite)

## increase max global setting
options(future.globals.maxSize = 1000 * 1024^2) # sets max to 1GB

## set up directory for outputs
dir.create("./outputs_Seurat_integrated/", showWarnings = FALSE)
results_folder = "./outputs_Seurat_integrated/"

## load data
raw_all = list()
samples = c('A','B','C','D')


raw_all[['A']] = LoadXenium('inputData/20231116__185301__9703-DM/output-XETG00077__0006740__9703-DM-1_ROI_A__20231116__185322/')

raw_all[['B']] = LoadXenium('inputData/20231116__185301__9703-DM/output-XETG00077__0006740__9301-DM-1_ROI_B__20231116__185322/')

raw_all[['C']] = LoadXenium('inputData/20231116__185301__9703-DM/output-XETG00077__0006740__9703-DM-1_ROI_C__20231116__185322/')

raw_all[['D']] = LoadXenium('inputData/20231116__185301__9703-DM/output-XETG00077__0006740__9703-DM-1_ROI_D__20231116__185322/')

## note - per squidpy documentation - can be useful to calculate QC metrics
## e.g. total counts, n genes by counts, cell_area, nucleus ratio
# question - how best to create those plots if using Seurat?

for (i in samples) {
raw_all[[i]] <- subset(raw_all[[i]], subset = nCount_Xenium > 1)
  Idents(raw_all[[i]]) <- i
  raw_all[[i]]$sample <- i
  raw_all[[i]] <- RenameCells(raw_all[[i]], add.cell.id = i)}

data_merged <- merge(raw_all[['A']], y = c(raw_all[['B']], raw_all[['C']], raw_all[['D']]), project = "xenium_ABCD")
names(data_merged@images)[1] = 'A'
names(data_merged@images)[2] = 'B'
names(data_merged@images)[3] = 'C'
names(data_merged@images)[4] = 'D'

# normalize and run dim reduction
data_merged <- SCTransform(data_merged, assay = "Xenium")
data_merged <- RunPCA(data_merged)

# look at elbow plot to estimate number of PCs to include
ElbowPlot(data_merged, ndims=50)

# set max dims
pcs=35

## Individual UMAP/Clustering (skipped)
# data_merged <- RunUMAP(data_merged, dims = 1:pcs)
# data_merged <- FindNeighbors(data_merged, dims = 1:pcs)
# data_merged <- FindClusters(data_merged, verbose = FALSE, resolution = 0.1)

## Integration
data_integrated <- IntegrateLayers(object = data_merged, method = RPCAIntegration, normalization.method = "SCT", verbose = FALSE)
data_integrated  # check dim reductions

## check integration
dim_checkPlot <- DimPlot(data_integrated, reduction = "integrated.dr") # looks improved but less overlap for section C

## Integrated UMAP/Clustering
data_integrated <- FindNeighbors(data_integrated, reduction = "integrated.dr", dims = 1:pcs)
data_integrated <- FindClusters(data_integrated, resolution = 0.1)
data_integrated <- RunUMAP(data_integrated, dims = 1:pcs, reduction = "integrated.dr")

##### UMAPs with integrated clusters (ImageDimPlot cannot plot 2 samples together. Had to find an alternative way to plot. Weisheng posted this on github: https://github.com/satijalab/seurat/issues/8713)

# add integrated clusters from metadata to individual objects
seurat_clusters_from_int = list()   
for (i in samples) {
seurat_clusters_from_int[[i]] = data_integrated@meta.data %>% rownames_to_column('barcode') %>% select(barcode, seurat_clusters) %>% `colnames<-`(c('barcode','seurat_clusters_from_int')) %>% right_join(raw_all[[i]]@meta.data %>% rownames_to_column('barcode'), by='barcode')
raw_all[[i]]@meta.data$seurat_clusters_from_int = seurat_clusters_from_int[[i]]$seurat_clusters_from_int
Idents(raw_all[[i]]) = raw_all[[i]]@meta.data$seurat_clusters_from_int}

# run workflow on individual samples to ensure have reductions and neighborhoods added
for (i in samples) {
raw_all[[i]] <- SCTransform(raw_all[[i]],assay = "Xenium")
raw_all[[i]] <- RunPCA(raw_all[[i]])
raw_all[[i]] <- RunUMAP(raw_all[[i]], dims = 1:pcs)
raw_all[[i]] <- FindNeighbors(raw_all[[i]], dims = 1:pcs)
raw_all[[i]] <- FindClusters(raw_all[[i]], verbose = FALSE, resolution = 0.1)}

# plot full sections per sample
full_sections_checkAB <- ImageDimPlot(raw_all[['A']],group.by='seurat_clusters_from_int') | ImageDimPlot(raw_all[['B']],group.by='seurat_clusters_from_int')

full_sections_checkCD <- ImageDimPlot(raw_all[['C']],group.by='seurat_clusters_from_int') | ImageDimPlot(raw_all[['D']],group.by='seurat_clusters_from_int')

# new fov by cropping - requires package sf
crop <- Crop(raw_all[['C']][["fov"]], x = c(2500, 3000), y = c(1000, 1600))
raw_all[['C']][["crop"]] <- crop
DefaultBoundary(raw_all[['C']][["crop"]]) <- "segmentation"
FOV_crop_check <- ImageDimPlot(raw_all[['C']], fov='crop', cols = "polychrome", axes=T, molecules=c("MYH11","ACTA2"))


FOV_check <-ImageDimPlot(raw_all[['C']], fov = "fov", molecules = c("MYH11","ACTA2"), nmols = 20000)



# initial coordinates triggered error since not cell centroids - use Xenium viewer to guide selection?
cropped.coords <- Crop(raw_all[['C']][["fov"]], x = c(1500, 2000), y = c(3000, 3500), coords = "plot")
raw_all[['C']][["zoom"]] <- cropped.coords

# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(raw_all[['C']][["zoom"]]) <- "segmentation"
zoom_check <- ImageDimPlot(raw_all[['C']], fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",coord.fixed = FALSE, molecules = c("MYH11","ACTA2"), nmols = 10000,alpha=0.1)
# interesting concentration of genes in one/few clusters


# then after generate integrated clusters, use reference data to annotate cell types
### singleR
aggexp <- AggregateExpression(data_integrated, group.by = "seurat_clusters", return.seurat = T)@assays$SCT@data
write.csv(aggexp, paste0(results_folder, 'xenium.aggexp.csv'))

BiocManager::install(c("SingleR", "celldex"))

library(SingleR)
library(celldex)

hpca.se <- HumanPrimaryCellAtlasData()

pred.hpca.se <- SingleR(test = aggexp, ref = hpca.se, labels = hpca.se$label.main)
View(as.data.frame(pred.hpca.se) %>% select(pruned.labels)) 
# not great predictions for 25 PCs; better with 35 PCs

## using ref data
# count matrix and annotation data downloaded from: https://singlecell.broadinstitute.org/single_cell/study/SCP1909/aortic-cellular-diversity-and-quantitative-genome-wide-association-study-trait-prioritization-through-single-nuclear-rna-sequencing-of-the-aneurysmal-human-aorta#study-download
ref.expression_matrix <- ReadMtx(
  mtx = "inputData/Broad_SCP1265/expression/600f1c04771a5b0d71956e3d/matrix.mtx.gz", features = "inputData/Broad_SCP1265/expression/600f1c04771a5b0d71956e3d/genes.tsv.gz",
  cells = "inputData/Broad_SCP1265/expression/600f1c04771a5b0d71956e3d/barcodes.tsv.gz")
ref.obj <- CreateSeuratObject(counts = ref.expression_matrix)
ref.meta = read.delim('inputData/Broad_SCP1265/metadata/ascending_descending_human_aorta_metadata.txt')[-1,] %>% select(NAME,cell_type__ontology_label)

unique(ref.meta$cell_type__ontology_label)

ref.meta$cluster = sapply(ref.meta$cell_type__ontology_label, function(x){strsplit(gsub('\\.','',x),' ')[[1]][1]})
ref.meta$cell_type = ref.meta$cell_type__ontology_label

head(ref.meta)

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
reference <- Reference(round(counts), cluster, round(nUMI)) # got non integers errors?

# data_integrated <- JoinLayers(data_integrated, assay = "integrated.dr") # issue for layered data access - https://github.com/satijalab/seurat/issues/8304; want to have joined integrated data for predictions but might need to use single slice first
query.counts <- GetAssayData(raw_all[["A"]], assay = "Xenium", slot = "counts")
coords <- GetTissueCoordinates(data_integrated, which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 4)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
data_integrated$predicted.celltype <- annotations
ImageDimPlot(data_integrated, group.by='predicted.celltype', cols = "polychrome", size = 0.75)
# decent predictions but fewer than expected smooth muscle so try increasing PCs and re-running integrated clustering
keep.cells <- Cells(data_integrated)[!is.na(data_integrated$predicted.celltype)]
#data_integrated.sub <- subset(data_integrated, cells = keep.cells) # have error here
## issue with this method is that not every cell is annotated

## add predicted cell types to split out raw data for plotting
## modifying first loop doesn't work since many NAs so make seurat cluster to prediction key instead OR use `pred.hpca.se`predictions instead

predicted_summary <- data_integrated@meta.data %>% select(seurat_clusters, predicted.celltype) %>% 
  group_by(seurat_clusters) %>% table() %>% as_tibble()

predicted_summary <- predicted_summary %>% group_by(seurat_clusters) %>% 
 summarize(max_prediction = predicted.celltype[which.max(n)])

new.cluster.ids <- predicted_summary$max_prediction
names(new.cluster.ids) <- predicted_summary$seurat_clusters
data_integrated <- RenameIdents(data_integrated, new.cluster.ids)
data_integrated@meta.data$clust.pred.celltype <- Idents(data_integrated)


for (i in samples) {
  #print(head(Idents(raw_all[[i]])))
  Idents(raw_all[[i]]) <- raw_all[[i]]@meta.data$seurat_clusters_from_int
  raw_all[[i]] <- RenameIdents(raw_all[[i]], new.cluster.ids)
  raw_all[[i]]@meta.data$clust.pred.celltype <- Idents(raw_all[[i]])
  } # issues with getting C&D identities set to new labels but expect code to work as is, if running from start to end

## manually add some neuron identities to D slide to fix color plotting
# raw_all[["D"]]@meta.data$clust.pred.celltype <- Idents(raw_all[["D"]])
# unique(Idents(raw_all[["A"]]))
# unique(Idents(raw_all[["D"]]))
# unique(raw_all[["D"]]@meta.data$clust.pred.celltype) ## none have neuronal included so unclear how it's being plotted if not one of the cluster identities

# plot full sections per sample
full_sections_checkAB <- ImageDimPlot(raw_all[['A']],group.by='clust.pred.celltype') | ImageDimPlot(raw_all[['B']],group.by='clust.pred.celltype')

full_sections_checkCD <- ImageDimPlot(raw_all[['C']],group.by='clust.pred.celltype') | ImageDimPlot(raw_all[['D']],group.by='clust.pred.celltype')

# write out integrated & annotated plots to file
# output to file
minWidth = 15
minHeight = 10
OutputName <- paste0('SlideOverlays_Annotated_', pcs, "-PCs_", 'A|B') 

png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
  print(full_sections_checkAB)
dev.off()
pdf(file=paste0(results_folder, OutputName,'.pdf'), width = minWidth, height = minHeight) # default = inches
  print(full_sections_checkAB)
dev.off()
## add svg outputs
svglite(file=paste0(results_folder, OutputName,'.svg'),  width = minWidth, height = minHeight) # default = inches
  print(full_sections_checkAB)
dev.off()


OutputName <- paste0('SlideOverlays_Annotated_', pcs, "-PCs_", 'C|D') 

png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
print(full_sections_checkCD)
dev.off()
pdf(file=paste0(results_folder, OutputName,'.pdf'), width = minWidth, height = minHeight) # default = inches
print(full_sections_checkCD)
dev.off()
## add svg outputs
svglite(file=paste0(results_folder, OutputName,'.svg'),  width = minWidth, height = minHeight) # default = inches
print(full_sections_checkCD)
dev.off()


## niche analysis for integrate data
#i="A" # test for single sample
for(i in samples){
  raw_all[[i]] <- BuildNicheAssay(object = raw_all[[i]], group.by = "clust.pred.celltype",
                                    niches.k = 5, neighbors.k = 30, fov='fov')
} # if run this way for multiple sections, the generated niches independently (so numbers don't correspond across sections)

niche.plots <- list()
for(i in samples){
  # plot niches
  niche.plots[[i]] <- ImageDimPlot(raw_all[[i]], group.by = "niches", size = 1.5, dark.background = T) + ggtitle("Niches") + # try named vector to match niche outputs
    scale_fill_manual(values = c("4"="#442288", "5"="#6CA2EA","3"="#B5D33D","2"= "#FED23F","1"= "#EB7D5B"))
}

ImageDimPlot(raw_all[["A"]], group.by = "niches", size = 1.5, dark.background = T) + ggtitle("Niches") + # try named vector to match niche outputs
  # scale_fill_manual(values = c("4"="#442288", "5"="#6CA2EA","3"="#B5D33D","2"= "#FED23F","1"= "#EB7D5B"))

# plot niches
minWidth = 13
minHeight = 10

for(i in samples){
  OutputName <- paste0('Slide_Niche+AnnoClusters_section', i) 
  
  png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
    print(ImageDimPlot(raw_all[[i]],group.by='clust.pred.celltype') + ggtitle("Clusters") | niche.plots[[i]])
  dev.off()
  
  pdf(file=paste0(results_folder, OutputName,'.pdf'), width = minWidth, height = minHeight) # default = inches
    ImageDimPlot(raw_all[[i]],group.by='clust.pred.celltype') + ggtitle("Clusters") | niche.plots[[i]]
  dev.off()
  
  ## add svg outputs
  svglite(file=paste0(results_folder, OutputName,'.svg'),  width = minWidth, height = minHeight) # default = inches
    ImageDimPlot(raw_all[[i]],group.by='clust.pred.celltype') + ggtitle("Clusters") | niche.plots[[i]]
  dev.off()
}

## add plots for genes of interest (do this before running niche, otherwise need to reset assay)
# genes of interest from: https://pubmed.ncbi.nlm.nih.gov/35762613/
DefaultAssay(raw_all[['A']]) <- "SCT"
# molecules vs features
ImageFeaturePlot(raw_all[['A']], features = c("GZMK"), molecules = c("GZMK","IL7R","LIF","GNLY", "VWF", "ACTA2"), axes = TRUE)
ImageFeaturePlot(raw_all[['A']], features = c("GZMK","IL7R","LIF","GNLY"), axes = TRUE, size = 0.75, cols = c("white", "red"))


minWidth = 8
minHeight = 8
for(i in samples){
  DefaultAssay(raw_all[[i]]) <- "SCT"
  
  OutputName <- paste0('Slide_GenesOfInterestByMolecule_section', i) 
  
  png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
  print(ImageFeaturePlot(raw_all[[i]], features = c("GZMK"), molecules = c("GZMK","IL7R","LIF","GNLY", "VWF", "ACTA2"), axes = TRUE))
  dev.off()
}

minWidth = 10
minHeight = 10
for(i in samples){
  
  OutputName <- paste0('Slide_GenesOfInterest_section', i) 
  
  png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
  print(ImageFeaturePlot(raw_all[[i]], features = c("GZMK","IL7R","LIF","GNLY"), axes = TRUE, size = 0.75, cols = c("white", "red")))
  dev.off()
}




# session cleanup
rm(ref.obj, ref.expression_matrix, ref.meta, reference)
rm(hpca.se, Q_mat, query, query.counts, RCTD)
gc()

#############
# save session info & Rdata to reload!
session_initial <- sessionInfo()
dir.create(paste0(results_folder,"Robjs")) # modify folder to have more informative name

## save intermediate files
saveRDS(raw_all,file=paste0(results_folder,'Robjs/','merged_raw_xenium.obj.Rds'))
saveRDS(data_integrated,file=paste0(results_folder,'Robjs/','integrated_xenium.Rds'))

# full session data
save.image(file=paste0(results_folder,"Robjs/integrated_annotated.RData"))
