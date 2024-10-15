# tried to update packages per thread to fix subsettig issue: https://github.com/satijalab/seurat/issues/8000; however this lead to issues localing Seurat package at all

library(Seurat) # using version 5
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyverse)
library(glmGamPoi)
library(DESeq2)
library(scran)
library(data.table)

# install.packages("harmony")
library(harmony)

# pcs optimization function
check_pcs <- function(scObj){
  pct <- scObj@reductions$pca@stdev / sum(scObj@reductions$pca@stdev) * 100
  pct <- scObj@reductions$pca@stdev / sum(scObj@reductions$pca@stdev) * 100
  cum <- cumsum(pct)
  co1 <- which(cum > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > .1), decreasing = T)[1] + 1
  pcs_return <- min(co1, co2)
  return(pcs_return)
}

#setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Raman_dayanidhi.raman@utoledo_SP1_weishwu_8699-DR/outputs')
setwd('/nfs/mm-isilon/bioinfcore/ActiveProjects/Raman_dayanidhi.raman@utoledo_SP1_weishwu_8699-DR/outputs/DMK_alternative_analysis/')

## initialize data
raw_all = list()

raw_all[['normal']] = Load10X_Spatial('../inputs/8699-DR_reprocessed/10x_analysis_8699-DR_reprocessed/Sample_8699-DR-1-GEX_TTTCGCAGGC-GTGCATCGAG/', filename='filtered_feature_bc_matrix.h5')
raw_all[['cancer']] = Load10X_Spatial('../inputs/8699-DR_reprocessed/10x_analysis_8699-DR_reprocessed/Sample_8699-DR-2-GEX_CATTTACCGT-TTTACAGGTG/', filename='filtered_feature_bc_matrix.h5')

samples = names(raw_all)

for (i in samples) {
  Idents(raw_all[[i]]) <- i
  raw_all[[i]]$sample <- i
  raw_all[[i]] <- RenameCells(raw_all[[i]], add.cell.id = i)}


raw_all <- lapply(raw_all, function(x){
  x1 <- PercentageFeatureSet(x, "^MT-", col.name = "percent_mito")
  return(x1)
  })

Assays(raw_all[['cancer']])
DefaultAssay(raw_all[['normal']])
DefaultAssay(raw_all[['cancer']])

# use k-means clustering from 10x processing to exclude normal lymphnode tissue from cancer slide
# use add on strategy (per https://github.com/satijalab/seurat/issues/530); with k=3 and cluster 3 as the lymph node stroma that should be excluded
cancer_10xK3clusters <- read_csv("../inputs/8699-DR_reprocessed/10x_analysis_8699-DR_reprocessed/Sample_8699-DR-2-GEX_CATTTACCGT-TTTACAGGTG/analysis/clustering/gene_expression_kmeans_3_clusters/clusters.csv")
head(raw_all$cancer)
# add cancer pre-fix to allow for merging with cancer object metadata
cancer_10xK3clusters$Barcode <- paste0("cancer_",cancer_10xK3clusters$Barcode)
colnames(cancer_10xK3clusters)[2] <- "k.means.10x_k3"
cancer_10xK3clusters$k.means.10x_k3 <- as.factor(cancer_10xK3clusters$k.means.10x_k3)
meta_alt_cancer <- raw_all$cancer@meta.data
meta_alt_cancer$Barcode <- row.names(meta_alt_cancer)
meta_alt_cancer <- base::merge(meta_alt_cancer, cancer_10xK3clusters, all.x=TRUE, sort=FALSE)
head(meta_alt_cancer)
rownames(meta_alt_cancer) <- meta_alt_cancer$Barcode
raw_all$cancer@meta.data <- meta_alt_cancer
str(raw_all$cancer@meta.data)


# filter
# In the "normal" section, nFeature & nCount are higher around the center gap (where the tissue was folded) and at the southeast edges (which may also had some folding). Used a nFeature filter to remove these spots.
# subset cancer data to exclude likely normal lymph tissue in metastatic cancer sample, based on k-means results from 10x processing
# %mito seems to be associated with anatomy location. No filtering was applied.
data_filt = list()
data_filt[['normal']] = subset(raw_all[['normal']], subset = nFeature_Spatial <= 2000)
idx <- which(raw_all$cancer@meta.data$k.means.10x_k3 != "3")
data_filt[['cancer']] = subset(raw_all[['cancer']], subset = k.means.10x_k3 == "3", invert = TRUE)
# warning for not validating Seurat object

output_Path <- "./DMK_alternative_analysis/"
dir.create(output_Path)
saveRDS(data_filt, paste0(output_Path, 'filtered_data_non-lymph.Rds'))

# QC plots after filtering
pdf(file = paste0(output_Path, 'qc_postFilter.pdf'),width=14,height=7)
print(VlnPlot(data_filt[['normal']], features = c("nFeature_Spatial"), pt.size=0.1) + ggtitle('nFeature:normal') + NoLegend() | VlnPlot(data_filt[['cancer']], features = c("nFeature_Spatial"), pt.size=0.1) + ggtitle('nFeature:cancer') + NoLegend())
print(VlnPlot(data_filt[['normal']], features = c("nCount_Spatial"), pt.size=0.1) + ggtitle('nCount:normal') + NoLegend() | VlnPlot(data_filt[['cancer']], features = c("nCount_Spatial"), pt.size=0.1) + ggtitle('nCount:cancer') + NoLegend())
print(FeatureScatter(data_filt[['normal']], feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + ggtitle('normal') + NoLegend() | FeatureScatter(data_filt[['cancer']], feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + ggtitle('cancer') + NoLegend())
print(VlnPlot(data_filt[['normal']], features = c("percent_mito"), pt.size=0.1) + ggtitle('%mito:normal') + NoLegend() | VlnPlot(data_filt[['cancer']], features = c("percent_mito"), pt.size=0.1) + ggtitle('%mito:cancer') + NoLegend())
print(SpatialFeaturePlot(data_filt[['normal']],'nFeature_Spatial',pt.size.factor=2.5) + ggtitle('nFeature:normal') | SpatialFeaturePlot(data_filt[['cancer']],'nFeature_Spatial',pt.size.factor=2.5) + ggtitle('nFeature:cancer'))
print(SpatialFeaturePlot(data_filt[['normal']],'nCount_Spatial',pt.size.factor=2.5) + ggtitle('nCount:normal') | SpatialFeaturePlot(data_filt[['cancer']],'nCount_Spatial',pt.size.factor=2.5) + ggtitle('nCount:cancer'))
print(SpatialFeaturePlot(data_filt[['normal']],'percent_mito',pt.size.factor=2.5) + ggtitle('%mito:normal') | SpatialFeaturePlot(data_filt[['cancer']],'percent_mito',pt.size.factor=2.5) + ggtitle('%mito:cancer'))
dev.off()

# further filter to exclude high mitochondria or try independent clustering as is?

# normalize data - independently
data_SCTcluster = list()
for (i in samples) {
  data_SCTcluster[[i]] = data_filt[[i]] %>% SCTransform(assay = "Spatial", new.assay.name = "SCT") %>% RunPCA(reduction.name = "pca_SCT")}

# in lieu of SCTransform, can also consider using log transformation
data_logcluster = list()
for (i in samples) {
  data_logcluster[[i]] = data_filt[[i]] %>%  NormalizeData() %>%
    ScaleData(vars.to.regress=c("nCount_Spatial", "nFeature_Spatial", "percent_mito")) %>% FindVariableFeatures(nfeatures = 1000) %>% RunPCA(reduction.name = "pca_logNorm")} # could lower variable features further since median number for normal slide is only 1000 or tailor to each slide

for(i in samples){
  print(ElbowPlot(data_SCTcluster[[i]], reduction = "pca_SCT") + ggtitle(paste0("Variance explained, SCT: ", i)))
  print(ElbowPlot(data_logcluster[[i]], reduction = "pca_logNorm") + ggtitle(paste0("Variance explained, logNorm: ", i)))
}

print(ElbowPlot(data_cluster[[i]]))

## save elbow plots
pdf(file = paste0(output_Path, 'ElbowPlots_individuals.pdf'),width=14,height=7)
for(i in samples){
  print(ElbowPlot(data_SCTcluster[[i]], reduction = "pca_SCT") + ggtitle(paste0("Variance explained, SCT: ", i)))
  print(ElbowPlot(data_logcluster[[i]], reduction = "pca_logNorm") + ggtitle(paste0("Variance explained, logNorm: ", i)))
}
dev.off()


## manually update pcs for each 
pcs_list_SCT <- list()
for(i in samples){
  scData <- data_SCTcluster[[i]]
  pct <- scData@reductions$pca@stdev / sum(scData@reductions$pca@stdev) * 100
  pct <- scData@reductions$pca@stdev / sum(scData@reductions$pca@stdev) * 100
  cum <- cumsum(pct)
  co1 <- which(cum > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > .1), decreasing = T)[1] + 1
  pcs_list_SCT[[i]] <- min(co1, co2)
}

pcs_list_log <- list()
for(i in samples){
  scData <- data_logcluster[[i]]
  pct <- scData@reductions$pca@stdev / sum(scData@reductions$pca@stdev) * 100
  pct <- scData@reductions$pca@stdev / sum(scData@reductions$pca@stdev) * 100
  cum <- cumsum(pct)
  co1 <- which(cum > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > .1), decreasing = T)[1] + 1
  pcs_list_log[[i]] <- min(co1, co2)
}

# check optimization results
pcs_list_SCT
pcs_list_log # much closer number of PCAs, maybe better to use?
res <- 0.3

for (i in samples) {
  pcs <- pcs_list_SCT[[i]]
  data_SCTcluster[[i]] = data_SCTcluster[[i]] %>% RunUMAP(dims = 1:pcs, reduction = "pca_SCT") %>%  FindNeighbors(dims = 1:pcs, reduction = "pca_SCT") %>% FindClusters(verbose = FALSE, resolution = res, )}

for (i in samples) {
  pcs <- pcs_list_log[[i]]
  data_logcluster[[i]] = data_logcluster[[i]] %>% RunUMAP(dims = 1:pcs, reduction = "pca_logNorm") %>%   FindNeighbors(dims = 1:pcs, reduction = "pca_logNorm") %>% FindClusters(verbose = FALSE, resolution = res)}

# check clustering
plot_grid(
  DimPlot(data_SCTcluster[['normal']]), 
  DimPlot(data_SCTcluster[['cancer']]),
  SpatialDimPlot(data_SCTcluster[['normal']], pt.size.factor=2.5),
  SpatialDimPlot(data_SCTcluster[['cancer']], pt.size.factor=2.5),
  nrow=2, align='hv') 
ggsave(paste0('separate_clustering-SCTnormalized.res',res,'.pdf'),width=10,height=8)

plot_grid(
  DimPlot(data_logcluster[['normal']]), 
  DimPlot(data_logcluster[['cancer']]),
  SpatialDimPlot(data_logcluster[['normal']], pt.size.factor=2.5),
  SpatialDimPlot(data_logcluster[['cancer']], pt.size.factor=2.5),
  nrow=2, align='hv') 
ggsave(paste0('separate_clustering-logNormalized.res',res,'.pdf'),width=10,height=8)

saveRDS(data_SCTcluster, paste0(output_Path, 'data_separate_cluster-SCT.Rds'))
saveRDS(data_logcluster, paste0(output_Path, 'data_separate_cluster-logNorm.Rds'))

# # markers provided by PI
# gene_markers = list()
# gene_markers[['Epithelial_normal']] = c('KRT5','KRT6A','KRT7','KRT14','KRT17')
# gene_markers[['Epithelial_cancer']] = c('KRT8','KRT19')
# gene_markers[['B_lymphocyte']] = c('IGHM','CD19')
# gene_markers[['Dentritic']] = c('ITGAX')
# gene_markers[['Fibroblasts']] = c('PDGFRA','CD34','COL1A1','COL1A2','COL5A1','LOXL1','LUM','FBLN1','FBLN2')
# gene_markers[['Adipocytes']] = c('PPARG','ADIPOQ','LEP')
# gene_markers[['T_cytotoxic']] = c('CD8A')
# gene_markers[['T_helper']] = c('CD4')
# gene_markers[['T_regs']] = c('IL2RA','IL7R','FOXP3')
# gene_markers[['T_lymphocyte']] = c('CD3G')
# 
# # marker gene spatial feature plots
# pdf(file = paste0(output_Path,'gene_markers_spatialfeatureplots.pdf'), width=16,height=7)
# gene_markers_spatialfeatureplots = list()
# for (i in names(gene_markers)) {
#   for (j in gene_markers[[i]]) {
#     id = paste(i, j, sep=':')
#     gene_markers_spatialfeatureplots[[id]] = SpatialFeaturePlot(data_cluster[[1]],j,pt.size.factor = 2.5) + ggtitle(paste0(id,':normal')) | SpatialFeaturePlot(data_cluster[[2]],j,pt.size.factor = 2.5) + ggtitle(paste0(id,':cancer'))
#     print(gene_markers_spatialfeatureplots[[id]])}}
# dev.off()
# 
# # marker gene feature plots
# gene_markers_featureplots = list()
# pdf(file = paste0(output_Path, 'gene_markers_featureplots.pdf'), width=16,height=7)
# for (i in names(gene_markers)) {
#   for (j in gene_markers[[i]]) {
#     id = paste(i, j, sep=':')
#     gene_markers_featureplots[[id]] = plot_grid(
#       FeaturePlot(data_cluster[[1]],j,pt.size=0.5) + ggtitle(paste0(id,':normal')), 
#       FeaturePlot(data_cluster[[2]],j,pt.size=0.1) + ggtitle(paste0(id,':cancer')),
#       align='hv', nrow=1)
#     print(gene_markers_featureplots[[id]])}}
# dev.off()
# 
# # marker gene ridge plots
# gene_markers_ridgeplots = list()
# pdf(file = paste0(output_Path, 'gene_markers_ridgeplots.pdf'), width=16,height=7)
# for (i in names(gene_markers)) {
#   for (j in gene_markers[[i]]) {
#     id = paste(i, j, sep=':')
#     gene_markers_ridgeplots[[id]] = plot_grid(
#       RidgePlot(data_cluster[[1]],j) + ggtitle(paste0(id,':normal')), 
#       RidgePlot(data_cluster[[2]],j) + ggtitle(paste0(id,':cancer')),
#       align='hv', nrow=1)
#     print(gene_markers_ridgeplots[[id]])}}
# dev.off() ## can maybe skip until have data merged?

# add labeled dotplot
# markers provided by PI
gene_markers = list()
gene_markers[['Epi_normal']] = c('KRT5','KRT6A','KRT7','KRT14','KRT17')
gene_markers[['Epi_cancer']] = c('KRT8','KRT19')
gene_markers[['B_lymph']] = c('IGHM','CD19')
gene_markers[['Dentritic']] = c('ITGAX')
gene_markers[['Fibroblasts']] = c('PDGFRA','CD34','COL1A1','COL1A2','COL5A1','LOXL1','LUM','FBLN1','FBLN2')
gene_markers[['Adip']] = c('PPARG','ADIPOQ','LEP')
gene_markers[['T_cytotox']] = c('CD8A')
gene_markers[['T_helper']] = c('CD4')
gene_markers[['T_regs']] = c('IL2RA','IL7R','FOXP3')
gene_markers[['T_lympho']] = c('CD3G')

dotplot_list <- list()
for(i in samples){
  dotplot_list[[i]] <- DotPlot(object = data_SCTcluster[[i]], features=gene_markers, cluster.idents=T, dot.min = 0) +
    theme(text=element_text(size=10), axis.text.x = element_text(angle = 90)) + labs(title = paste0("Provided marker expression, SCT - ", i))
  print(dotplot_list[[i]])
}

pdf(file = paste0(output_Path, 'gene_markers_dotplots-SCT.pdf'), width=18,height=7)
for (i in names(dotplot_list)) {
  print(dotplot_list[[i]])}
dev.off()

## repeat for log normed?
dotplot_list <- list()
for(i in samples){
  dotplot_list[[i]] <- DotPlot(object = data_logcluster[[i]], features=gene_markers, cluster.idents=T, dot.min = 0) +
    theme(text=element_text(size=10), axis.text.x = element_text(angle = 90)) + labs(title = paste0("Provided marker expression, logNorm - ", i))
  print(dotplot_list[[i]])
}

pdf(file = paste0(output_Path, 'gene_markers_dotplots-logNorm.pdf'), width=18,height=7)
for (i in names(dotplot_list)) {
  print(dotplot_list[[i]])}
dev.off()


## for next steps - plot initial PCA and then run integration with Harmony and check integration quality as well as if clusters better correspond to expected cell types; logNormalized approach seems to have better scaling matches between slides based on marker plots

# merge/integrate data (e.g.: https://www.10xgenomics.com/analysis-guides/correcting-batch-effects-in-visium-data)
names(data_filt)
#data_merged <- merge(data_filt[[1]], y=data_filt[[2]], project = "Raman")
data_merged <- merge(data_SCTcluster[[1]], y=data_SCTcluster[[2]], project = "Raman", merge.data= FALSE, merge.dr = FALSE)
data_merged$orig.ident
head(colnames(data_merged)); tail(colnames(data_merged))
LayerData(data_merged)
names(data_merged@images) <- c("normal", "cancer")
names(data_merged@meta.data)

# SpatialFeaturePlot(data_merged, "nCount_Spatial")
# SpatialFeaturePlot(data_merged, "nFeature_Spatial")
# SpatialFeaturePlot(data_merged, "percent_mito")

# PCA
experiment.merged <- data_merged %>%  SCTransform(variable.features.n = 1000, assay = "Spatial") %>% RunPCA(reduction.name = "pca")
# re-run PCA since cleared when merged? have counts.1 and counts.2

pcs <- check_pcs(experiment.merged)

# UMAP
experiment.merged <- RunUMAP(experiment.merged,
                             reduction = "pca",
                             dims = 1:pcs,
                             verbose = FALSE)
#rm(data_merged)

Idents(experiment.merged) <- experiment.merged@meta.data$sample

p1 <- DimPlot(experiment.merged,
        reduction = "pca",
        split.by = "sample",
        shuffle = TRUE) +
  scale_color_viridis_d() + NoLegend() #+ ggtitle("PCA")

p2 <- DimPlot(experiment.merged,
        reduction = "umap",
        split.by = "sample",
        shuffle = TRUE) +
  scale_color_viridis_d() #+ ggtitle("UMAP")
# output to file - with merge, normal and cancer seem to have some overlaps but dominated by cancer spots (as expected)

# create combined plots
p1+p2+plot_annotation(
  title = 'Merged data - pre-integration SCT')

ggsave(paste0('UMAP-PCA_pre-integration-SCT.pcs',pcs,'.pdf'),width=13,height=6)

# experiment.merged <- readRDS("/nfs/mm-isilon/bioinfcore/ActiveProjects/Raman_dayanidhi.raman@utoledo_SP1_weishwu_8699-DR/outputs/DMK_alternative_analysis/merged-preIntegration_non-lymph.Rds")
# Idents(experiment.merged) # identities for SCT version are cancer normal not clusters
# 
# p1 <- DimPlot(experiment.merged,
#               reduction = "pca",
#               split.by = "sample",
#               shuffle = TRUE) +
#   scale_color_viridis_d() + NoLegend() #+ ggtitle("PCA")
# 
# p2 <- DimPlot(experiment.merged,
#               reduction = "umap",
#               split.by = "sample",
#               shuffle = TRUE) +
#   scale_color_viridis_d() #+ ggtitle("UMAP")
# # output to file - with merge, normal and cancer seem to have some overlaps but dominated by cancer spots (as expected)
# 
# # create combined plots
# p1+p2+plot_annotation(
#   title = 'Merged data - pre-integration SCT')
# 
# ggsave(paste0('UMAP-PCA_pre-integration-SCT.res',res,'.pdf'),width=13,height=6)


# continue with Harmony integratation
experiment.harmony <- RunHarmony(experiment.merged,
                                 group.by.vars = "sample", # for real batch correction use "slide"
                                 assay.use = "data",
                                 verbose = FALSE,
                                 ncores=2,
                                 theta=1) # default is theta=2; can try lower thetas (1, 0.5) to force less overlaps
## stopping point!! Ended up with two SCT transformations, need to correct
ElbowPlot(experiment.harmony)
# 16 is probably high, 7 is probably too low (based on elbow plot and downstream UMAP)

## output Seurat object to file for integrated results
saveRDS(experiment.harmony, paste0('integrated_harmony_non-lymph-preClustering.Rds'))

check_pcs(experiment.harmony) #14 with log, 16 with SCT

pcs <- 14 # 14 seemed to look the best
experiment.harmony <- RunUMAP(experiment.harmony,
                              reduction = "harmony",
                              dims = 1:pcs,
                              verbose = FALSE)
# also calculateTSNE?
# experiment.harmony <- RunTSNE(experiment.harmony,
#                               reduction="harmony")

p1 <- DimPlot(experiment.harmony,
        reduction = "umap",
        group.by = "sample",
        shuffle = TRUE) +
  scale_color_viridis_d()

p2 <- DimPlot(experiment.harmony,
        reduction = "pca",
        group.by = "sample",
        shuffle = TRUE) +
  scale_color_viridis_d()

p1+p2+plot_annotation(
  title = 'Merged data - harmony integration with SCT input')

ggsave(paste0('UMAP-PCA_SCT_harmony_integrated.pcs',pcs,'.pdf'),width=13,height=6)

# vast improvement in overlaps but not the best UMAP structure

# try clustering
experiment.harmony <- FindNeighbors(experiment.harmony, reduction = "harmony", dims = 1:pcs) 
experiment.harmony <- FindClusters(experiment.harmony,
                                  resolution = seq(0.2, 0.8, 0.1),
                                  verbose = FALSE)
head(experiment.harmony@meta.data)

lapply(grep("snn", colnames(experiment.harmony@meta.data), value = TRUE),
       function(res){
         DimPlot(experiment.harmony,
                 reduction = "umap",
                 group.by = res,
                 shuffle = TRUE) +
           scale_color_viridis_d(option = "turbo")
       }) # with logNorm input, tried pcs 12 and 14 - have not very not very distinct cluster, but res 0.3 looks okay

head(experiment.harmony@meta.data)
res <- 0.3
experiment.harmony@meta.data$seurat_clusters <- experiment.harmony@meta.data$SCT_snn_res.0.3

# "orig.ident" = original identity of cells
DimPlot(experiment.harmony, group.by = "SCT_snn_res.0.3", label = TRUE) + NoLegend() # expect 11 clusters if identities set correctly

minWidth = 5; minHeight=5
ggsave(filename=paste0("UMAP-integrated-SCTharmony_", pcs, "pcs", res, "res", ".png"),width=minWidth, height=minHeight, unit="in")

# "ident" = identity, which are clusters
DimPlot(experiment.harmony, group.by = "SCT_snn_res.0.3", 
        split.by = 'sample')

minWidth = 10.5; minHeight=5
ggsave(filename=paste0("UMAP-integrated-SCTharmony_", pcs, "pcs", res, "res", "-BySample.png"),width=minWidth, height=minHeight, unit="in")


## add plot for clusters overlaid on slides
SpatialDimPlot(experiment.harmony, label=FALSE, label.size = 3, image.alpha = 0.5, group.by = "SCT_snn_res.0.3")

minWidth = 10.5; minHeight=5
ggsave(filename=paste0("ClusterSlideOverlays- integrated_SCT_harmony--", pcs, "pcs", res, "res", ".png"),width=minWidth, height=minHeight, unit="in")

# set idents for selected pcs/res combination
Idents(experiment.harmony) <- experiment.harmony@meta.data$SCT_snn_res.0.3

## output Seurat object to file for integrated results
saveRDS(experiment.harmony, paste0('integrated_SCT_harmony_non-lymph_clustered.Rds'))

# next - create dotplot of marker genes for integrated/batch corrected data
dot.harmony <- DotPlot(object = experiment.harmony, features=gene_markers, cluster.idents=T, dot.min = 0) +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 90)) + labs(title = paste0("Provided marker expression - integrated data"))
print(dot.harmony)


pdf(file = paste0('gene_markers_dotplots-SCT_harmony-integrated.pdf'), width=18,height=7)
print(dot.harmony)
dev.off()

## then run scCATCH or use another tool for label transfer to annotate clusters?

# read in Grey et all marker genes and create signature file
# data from: https://pubmed.ncbi.nlm.nih.gov/35617956/
Grey_majorCellTypes <- read_csv("../../inputs/Grey_etAl_2022_BC-celltypes/MajorCellTypeSignatures_GreyetAl.csv", skip = 2)
head(Grey_majorCellTypes) # columns titles = major type; need to clarify "AV", "HS", and "BA" labels

Grey_minorCellTypes <- read_csv("../../inputs/Grey_etAl_2022_BC-celltypes/CellSubTypeSignatures_GreyetAl.csv", skip = 3)
head(Grey_minorCellTypes)

Grey_majorCellTypes <- as.data.frame(Grey_majorCellTypes)
Grey_minorCellTypes <- as.data.frame(Grey_minorCellTypes)

# re-order tables to match input signature structure file for VISION package (https://yoseflab.github.io/VISION/articles/Signatures.html)
# or install/load package to create signature object (might be easier since will need signed list instead of reworking table to be .gmt format)

colnames(Grey_majorCellTypes)[1:3] <- c("Alveolar", "HormoneSensing", "BasalCells")
# from paper: AV = alveolar cells; HG = Hormone sensing, BA = basal cells
colnames(Grey_majorCellTypes)[5] <- "Vasc.Lymphatic"

# add some of the key minor celltypes
Grey_minorCellTypesSubset <- Grey_minorCellTypes[ ,c(10:17)]
colnames(Grey_minorCellTypesSubset) <- c("Lymph_endo", "Vasc_endo", "Pericyte", "Myeloid", "NK_cell", "Tcell", "Bcell", "PlasmaCell")

## create function/loop to process each signal?
# signatureList <- list()
# sigName <- colnames(Grey_majorCellTypes)
# for(i in sigName){
#   signatureList[[i]] <- na.omit(Grey_majorCellTypes[, i])
#   signatureList[[i]] <- rep(1, length(signatureList[[i]]))
#   names(signatureList[[i]]) <- na.omit(Grey_majorCellTypes[, i])
# }
# head(signatureList) # getting major signature list
#  
# sigName <- colnames(Grey_minorCellTypesSubset)
# for(i in sigName){
#   signatureList[[i]] <- na.omit(Grey_minorCellTypesSubset[, i])
#   signatureList[[i]] <- rep(1, length(signatureList[[i]]))
#   names(signatureList[[i]]) <- na.omit(Grey_minorCellTypesSubset[, i])
# }
# str(signatureList)

# then create individual signatures (e.g. https://yoseflab.github.io/VISION/articles/Signatures.html#creating-a-signature-object-in-r)
# SigObjects <- list()
# for(i in names(signatureList)){
#   print(i) # first 6 = major types; remainder = minor types
#   SigObjects[[i]] <- createGeneSignature(name= paste0("GreyetAl_", i), sigData = signatureList[[i]])
# }
# 
# head(SigObjects)
# experiment.harmony
# 
# require(devtools)
# install_github("YosefLab/VISION")
# library(VISION)
# 
# vision.obj <- Vision(experiment.harmony,
#                      assay = "SCT",
#                      cellsPerPartition = 5,
#                      meta = experiment.harmony@meta.data, 
#                      signatures = SigObjects,
#                      dimRedComponents = pcs,
#                      projection_methods = NULL)
# 
# #viewResults(vision.obj) #need to run with additional options to make accessible outside of VNC
# # Display autocorrelation coefficients, p-values for signatures
# vis <- analyze(vision.obj)
# str(vis)
# 
# # view results within R
# getSignatureAutocorrelation(vis) # returns results but doesn't have existing clusters included
# getSignatureScores(vis) # kind of what I want, but based on cells and not pre-exisiting clusters
# getMetaAutocorrelation(vis)


# depending on results - can maybe try to use same tool as used for mSigDB with Vision signature input
#devtools::install_github("arc85/singleseqgset")
#install.packages("heatmap3")
library(singleseqgset)
library(heatmap3)
library(msigdbr)

# example mSigDB from Raman has similar structure as `signature list` but with genes as list values, so can re-run to match 
signatureList <- list()
sigName <- colnames(Grey_majorCellTypes)
for(i in sigName){
  signatureList[[i]] <- na.omit(Grey_majorCellTypes[, i])
}
head(signatureList) # getting major signature list

signatureList_minor <- list()
sigName <- colnames(Grey_minorCellTypesSubset)
for(i in sigName){
  signatureList_minor[[i]] <- na.omit(Grey_minorCellTypesSubset[, i])
}
str(signatureList)
names(signatureList) <- paste0("Grey-major_", names(signatureList))

names(signatureList_minor) <- paste0("Grey-minor_", names(signatureList_minor))

# logfc.data <- singleseqgset::logFC(cluster.ids = experiment.harmony@meta.data$seurat_clusters,
#                                    expr.mat = experiment.harmony@assays$SCT@data)
# names(logfc.data)

logfc.data <- singleseqgset::logFC(cluster.ids = experiment.harmony@meta.data$seurat_clusters,
                                   expr.mat = experiment.harmony@assays$SCT@data)
names(logfc.data)

gse.grey <- wmw_gsea(expr.mat = experiment.harmony@assays$SCT@data,
                      cluster.cells = logfc.data[[1]],
                      log.fc.cluster = logfc.data[[2]],
                      gene.sets = signatureList)

names(gse.grey)

res.stats <- gse.grey[["GSEA_statistics"]]
res.pvals <- gse.grey[["GSEA_p_values"]]
res.pvals <- apply(res.pvals, 2, p.adjust, method = "fdr") #correct for multiple comparisons
res.stats[order(res.stats[ ,1], decreasing = TRUE)[1:10], ] # top genesets by zscore
res.pvals[order(res.stats[,1],decreasing=TRUE)[1:10],]
rm(gse.grey)

# heatmap3(res.stats, Colv = NA, cexRow = 0.5, cexCol = 1, scale="row", showColDendro = F, showRowDendro = F)
# 
# heatmap3(res.pvals[order(res.stats[,1],decreasing=TRUE)[1:length(names(signatureList))],], Colv = NA, cexRow = 0.5, cexCol = 1, scale="row", showColDendro = F, showRowDendro = F) ## plotting only worked once made plotting window extra large & ensured that 
# # maybe helpful but unclear what scale is (- vs + = smaller pval?)

# ## if looks useful, format as single output table
str(res.stats); str(res.pvals)
res.stats$GreyetAl_Annotations <- row.names(res.stats)
res.pvals <- as.data.frame(res.pvals)
res.pvals$GreyetAl_Annotations <- row.names(res.pvals)
res.merged <- merge(res.stats, res.pvals,
                    by= "GreyetAl_Annotations",
                    suffixes = c("-zscore","-adjpv"),
                    sort = TRUE)
colnames(res.merged)
col_order <- sort(colnames(res.merged)[2:23])
res.merged <- res.merged[,c("GreyetAl_Annotations", col_order)]
# 
# # write table to file
TablePathName <- paste0("./","ClusterAnnotationPredictions/")
dir.create(TablePathName)

fwrite(res.merged, paste0(TablePathName,'GreyetAl-BreastCancerAnnotation_ClusterEnrichments.csv'),
       col.names=T)

res.stats <- as.matrix(res.stats[,c(1:11)])
GSEA_plot <- pheatmap::pheatmap(res.stats[order(res.stats[,1],decreasing=TRUE)[c(1:length(signatureList))],], 
                                scale = "row", cluster_cols = TRUE, angle_col = 0, 
                                fontsize_row = 10, fontsize_col = 12)

# export to file - to same directory as tables
minWidth = 6.5 ; minHeight = 7
png(file=paste(TablePathName, paste0('GreyetAl-Enrichments_', pcs,'PCs_IntegratedData','.png'), sep="_"), width = minWidth*100, height = minHeight*100)
print(GSEA_plot)
dev.off()
pdf(file=paste(TablePathName, paste0('GreyetAl_', pcs,'PCs_IntegratedData','.pdf'), sep="_"), 
    width = minWidth, height = minHeight)
print(GSEA_plot)
dev.off()

## add scCATCH predictions since these predictions don't seem as informative as ideal; next time would probably separate major/minor gene sets and run enrichments instead of pooling together

## ScCatch - see available tissues: https://github.com/ZJUFanLab/scCATCH/wiki
library(scCATCH)

catch_obj = createscCATCH(data = experiment.harmony@assays$SCT@data, cluster = as.character(Idents(experiment.harmony)))

cell_types = c("Breast", "Mammary gland","Mammary epithelium", "Pluripotent stem cell", "Blood", "Peripheral blood", "Serum") # pulled from reference via wiki
new_cellmatch = cellmatch[cellmatch$species == 'Human' & cellmatch$tissue %in% cell_types, ]

catch_obj <- findmarkergene(object = catch_obj,
                            if_use_custom_marker = TRUE,
                            marker = new_cellmatch, 
                            use_method = "2") # finished in < 1 minute for 2 samples
catch_obj = findcelltype(catch_obj)

str(catch_obj)
head(catch_obj@markergene) # marker gene table
head(catch_obj@celltype) # prediction table


# write out scCATCH marker and predictions to file
# use fwrite to have better formatting
fwrite(catch_obj@markergene, paste0(TablePathName,'scCATCH_clusterMarkersGenes.csv'),
       col.names=T)

fwrite(catch_obj@celltype, paste0(TablePathName,'scCATCH_CellTypePredictions.csv'),
       col.names=T)

rm(catch_obj, new_cellmatch)
rm(Grey_majorCellTypes, Grey_minorCellTypes, Grey_minorCellTypesSubset, GSEA_plot, logfc.data, p1, p2, res.pvals, res.stats)
gc()

## add empirical marker genes or hold off?


## Consider feature plot for EPCAM to show cancer/non-cancer
FeaturePlot(experiment.harmony, features = "EPCAM", split.by = "sample")
SpatialFeaturePlot(experiment.harmony, features = "EPCAM")
####
## add plots for genes of interest after integration
# markers provided by PI
gene_markers = list()
gene_markers[['Epi_normal']] = c('KRT5','KRT6A','KRT7','KRT14','KRT17')
gene_markers[['Epi_cancer']] = c('KRT8','KRT19')
gene_markers[['Fibroblasts']] = c('PDGFRA','CD34','COL1A1','COL1A2','COL5A1','LOXL1','LUM','FBLN1','FBLN2')
gene_markers[['CAFibroblasts']] = c('ACTA2','AIFM2') # AIFM2 = FSP1
gene_markers[['Adip']] = c('PPARG','ADIPOQ','LEP')
gene_markers[['B_lymph']] = c('IGHM','CD19')
gene_markers[['Dentritic']] = c('ITGAX')
gene_markers[['T_cytotox']] = c('CD8A')
gene_markers[['T_helper']] = c('CD4')
gene_markers[['T_regs']] = c('IL2RA','IL7R','FOXP3')
gene_markers[['T_lympho']] = c('CD3G')
output_Path <- "./"
# marker gene spatial feature plots
pdf(file = paste0(output_Path,'integrated-gene_markers_spatialfeatureplots.pdf'), width=16,height=7)
gene_markers_spatialfeatureplots = list()
for (i in names(gene_markers)) {
  for (j in gene_markers[[i]]) {
    id = paste(i, j, sep=':')
    gene_markers_spatialfeatureplots[[id]] = SpatialFeaturePlot(experiment.harmony,j,pt.size.factor = 2.5)
    print(gene_markers_spatialfeatureplots[[id]])
    }}
dev.off()

# marker gene feature plots
gene_markers_featureplots = list()
pdf(file = paste0(output_Path, 'integrated-gene_markers_featureplots.pdf'), width=7,height=7)
for (i in names(gene_markers)) {
  for (j in gene_markers[[i]]) {
    id = paste(i, j, sep=':')
    gene_markers_featureplots[[id]] = 
      FeaturePlot(experiment.harmony,j,pt.size=0.5) + ggtitle(paste0(id, " (integrated data)"))
    print(gene_markers_featureplots[[id]])}}
dev.off()

# marker gene ridge plots
gene_markers_ridgeplots = list()
pdf(file = paste0(output_Path, 'integrated-gene_markers_ridgeplots.pdf'), width=16,height=7)
for (i in names(gene_markers)) {
  for (j in gene_markers[[i]]) {
    id = paste(i, j, sep=':')
    gene_markers_ridgeplots[[id]] = 
      RidgePlot(experiment.harmony,j, layer ="data", sort = TRUE, log = TRUE) + ggtitle(paste0(id,' (integrated data)'))
    print(gene_markers_ridgeplots[[id]])}}
dev.off()

# add labeled dotplot
# markers provided by PI
gene_markers = list()
gene_markers[['Epi_normal']] = c('KRT5','KRT6A','KRT7','KRT14','KRT17')
gene_markers[['Epi_cancer']] = c('KRT8','KRT19')
gene_markers[['Fibroblasts']] = c('PDGFRA','CD34','COL1A1','COL1A2','COL5A1','LOXL1','LUM','FBLN1','FBLN2')
gene_markers[['CAFibroblasts']] = c('ACTA2','AIFM2') # AIFM2 = FSP1
gene_markers[['Adip']] = c('PPARG','ADIPOQ','LEP')
gene_markers[['B_lymph']] = c('IGHM','CD19')
gene_markers[['Dentritic']] = c('ITGAX')
gene_markers[['T_cytotox']] = c('CD8A')
gene_markers[['T_helper']] = c('CD4')
gene_markers[['T_regs']] = c('IL2RA','IL7R','FOXP3')
gene_markers[['T_lympho']] = c('CD3G')



pdf(file = paste0(output_Path, 'integrated-gene_markers_dotplots.pdf'), width=18,height=7)
DotPlot(object = experiment.harmony, features=gene_markers, cluster.idents=T, dot.min = 0) +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 90)) + labs(title = paste0("Provided marker expression - integrated data"))
dev.off()

## add tables of empirical markers as well; 
# skipping prep SCT based on errors and recommendation here: https://github.com/satijalab/seurat/issues/7313
#cluster_markers_all <- FindAllMarkers(experiment.harmony, only.pos = TRUE, test.use = "MAST", assay = "SCT")


###
# clean up environment & then save session data
rm(data_cluster, split_seurat, scData, data_merged, data_filt)
gc()

save.image("/nfs/mm-isilon/bioinfcore/ActiveProjects/Raman_dayanidhi.raman@utoledo_SP1_weishwu_8699-DR/outputs/DMK_alternative_analysis/SessionData.RData")
