## Seurat analysis of Xenium aorta data

## load libraries
library(withr)
toolPath <- "/home/damki/R/x86_64-pc-linux-gnu-library/Xenium_patch/"
# with_libpaths(new=toolPath, devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)) # expect to only need single installation; for GL, did not update any dependencies
# with_libpaths(new=toolPath, install.packages("Matrix"))

.libPaths(c(toolPath,"/usr/local/lib/R/*"))
library(Seurat)
library(spacexr)
library(MAST)

library(future)
currentCores <- availableCores()[[1]]
plan("multisession", workers = currentCores) # set if BioC not already set
library(tidyverse)
library(plyr)
library(svglite)
library(matrixStats)

## increase max global setting
options(future.globals.maxSize = 100000 * 1024^2) # sets max to 80GB (needed to increase)

## check/set working directory
setwd("/nfs/mm-isilon/bioinfcore/ActiveProjects/Narla_gnrla_SP1-Xenium_damki_9638-KZ")

## set up directory for outputs
results_folder = "./outputs/"
dir.create(results_folder, showWarnings = FALSE)
dir.create(paste0(results_folder,"Robjs"), showWarnings = FALSE) # modify folder to have more informative name?

## load data
raw_list = list()
samples = c('Rat.PDX', 'Mouse.PDX')

# add wrapper to ensure that data is only loaded once
raw_list[['Rat.PDX']] = LoadXenium('inputs/20240301__220829__9638-KZ/output-XETG00077__0021019__9638-KZ-1_AOI_A__20240301__220842/', fov = "Rat.fov")
gc() # clear cache after loading data

raw_list[['Mouse.PDX']] = LoadXenium('inputs/20240301__220829__9638-KZ/output-XETG00077__0021021__9638-KZ-2_AOI_A__20240301__220842/', fov = "Mouse.fov")
gc() # clear cache after loading data

## note - per squidpy documentation - can be useful to calculate QC metrics, e.g. total counts, n genes by counts, cell_area, nucleus ratio
# question - how best to create those plots if using Seurat?

for (i in samples) {
raw_list[[i]] <- subset(raw_list[[i]], subset = nCount_Xenium > 5) # increase threshold
  Idents(raw_list[[i]]) <- i
  raw_list[[i]]$sample <- i
  raw_list[[i]] <- RenameCells(raw_list[[i]], add.cell.id = i)}


## save intermediate file
saveRDS(raw_list,file=paste0(results_folder,'Robjs/','raw_list_xenium.obj.Rds'))
## loaded intermediate as "raw_list" instead of "raw_all" in original analysis


head(raw_list)
# head(raw_list$Mouse.PDX)
## expect cell Ids to have format of "Mouse.PDX_aaachjfa.."

# Prior to merge, read in 10x kmer results and use to subset PDX data
Rat_Kmean2 <- read_csv("./inputs/20240301__220829__9638-KZ/output-XETG00077__0021019__9638-KZ-1_AOI_A__20240301__220842/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv")
Mouse_Kmean2 <- read_csv("./inputs/20240301__220829__9638-KZ/output-XETG00077__0021021__9638-KZ-2_AOI_A__20240301__220842/analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv")
head(raw_list$Rat.PDX)
# add pre-fix to kmean table allow for merging with existing object metadata
Rat_Kmean2$Barcode <- paste0("Rat.PDX_", Rat_Kmean2$Barcode)
colnames(Rat_Kmean2)[2] <- "Kmeans.2_10x"
Rat_Kmean2$Kmeans.2_10x <- as.factor(Rat_Kmean2$Kmeans.2_10x)

meta_alt_rat <- raw_list$Rat.PDX@meta.data
meta_alt_rat$Barcode <- rownames(meta_alt_rat)
meta_alt_rat <- base::merge(meta_alt_rat, Rat_Kmean2, all.x=TRUE, sort=FALSE)
rownames(meta_alt_rat) <- meta_alt_rat$Barcode
meta_alt_rat$Barcode <- NULL
head(meta_alt_rat)
raw_list$Rat.PDX@meta.data <- meta_alt_rat
str(raw_list$Rat.PDX@meta.data)

Mouse_Kmean2$Barcode <- paste0("Mouse.PDX_", Mouse_Kmean2$Barcode)
colnames(Mouse_Kmean2)[2] <- "Kmeans.2_10x"
Mouse_Kmean2$Kmeans.2_10x <- as.factor(Mouse_Kmean2$Kmeans.2_10x)

meta_alt_mouse <- raw_list$Mouse.PDX@meta.data
meta_alt_mouse$Barcode <- rownames(meta_alt_mouse) 
meta_alt_mouse <- base::merge(meta_alt_mouse, Mouse_Kmean2, all.x = TRUE, sort = FALSE)
rownames(meta_alt_mouse) <- meta_alt_mouse$Barcode
meta_alt_mouse$Barcode <- NULL
head(meta_alt_mouse)
raw_list$Mouse.PDX@meta.data <- meta_alt_mouse
str(raw_list$Mouse.PDX@meta.data)

# save annotated data
saveRDS(raw_list,file=paste0(results_folder,'Robjs/','raw_kmeans-anno_xenium.obj.Rds'))

# check totals
summary(raw_list$Rat.PDX@meta.data$Kmeans.2_10x)
# 1       2 
# 1852167  492756 
summary(raw_list$Mouse.PDX@meta.data$Kmeans.2_10x)
# 1       2 
# 1024741   89557

#### -----


## Create list of raw data that is subset only to the correct kmeans

#### ----- 

## create a plot to summarize total cell numbers and proportions of likely host vs cancer cells
# combine metadata
str()
combined_meta <- rbind(raw_list$Rat.PDX@meta.data, raw_list$Mouse.PDX@meta.data)

# reassign factor values to label (per:https://stackoverflow.com/questions/11810605/replace-contents-of-factor-column-in-r-dataframe)
combined_meta$PutativeOrigin <- combined_meta$Kmeans.2_10x
combined_meta$PutativeOrigin <- revalue(combined_meta$PutativeOrigin, c("1" = "tumor", "2" = "host"))
combined_meta$sample <- as.factor(combined_meta$sample)
str(combined_meta)

# create plot with proportions from each
ggplot(combined_meta, aes(fill=PutativeOrigin, x=sample)) + 
  geom_bar(position="fill", stat="count") +
  labs(title = "Total segmented cells, pre-filtering", 
       subtitle = paste0("Total mouse = ", round(sum(combined_meta$sample == "Mouse.PDX")/1000000, digits =4), "E6, Total rat = ", round(sum(combined_meta$sample == "Rat.PDX")/1000000, digits = 4), "E6"),
       ylab = "proportion of segmented cells per slide") + 
  theme(plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="black"),
        plot.title=element_text(size=14, hjust=0.5, face="bold")) ## add total number of "cells" pre-filtering

outputPath <- "./outputs/QC_plots/"
dir.create(outputPath)
ggsave(paste0(outputPath,'pre-filtering_cell-origin_barplot','.pdf'),width=8,height=10)

#save.image("./outputs/Robjs/SessionData.RData")

## Summarize QC metrics and cell counts for each k-mean cluster 
head(combined_meta);tail(combined_meta)
Count_summary <- combined_meta %>% select(sample, PutativeOrigin) %>% table() %>% as_data_frame() %>% dplyr::rename(CellCount = n)


combined_meta <- combined_meta %>% 
  mutate(CombinedLabel = paste0(sample, "_",PutativeOrigin))
head(combined_meta)

QC_summary <- combined_meta %>% select(sample, PutativeOrigin, nCount_Xenium, nFeature_Xenium) %>% group_by(sample, PutativeOrigin) %>%
  summarise_all(.funs = mean) %>% 
  dplyr::rename(AvgPerCell_nCount = nCount_Xenium, 
                AvgPerCell_nFeature = nFeature_Xenium)

QC_summary; Count_summary

SummaryTable <- merge(Count_summary, QC_summary)

write.table(as.data.frame(SummaryTable), file = paste0(outputPath, "SummaryTable_preFilter.csv"), sep = ",", row.names = FALSE)

#######
## Confirm host vs cancer with human specific gene targets
HumantoRat_scores <- read_csv(file = "./inputs/10x_probeSet_crossSpecies/human_multi_to_rat.csv")
HumantoMouse_scores <- read_csv(file = "./inputs/10x_probeSet_crossSpecies/human_multi_to_mouse.csv")
head(HumantoRat_scores); head(HumantoMouse_scores)

HumantoRat_lowOffTarget <- HumantoRat_scores %>% filter(fraction_probes_binds_ortholog == 0)
HumantoMouse_lowOffTarget <- HumantoMouse_scores %>% filter(fraction_probes_binds_ortholog == 0)

# intersect to find more unique to human genes
ProbeSet_lowOffTarget <- join(HumantoRat_lowOffTarget, HumantoMouse_lowOffTarget, by = c("gene_id", "gene_name"), type= "inner")
str(ProbeSet_lowOffTarget)

# check ACTG2 expression since that was DE in both species 
names(raw_list)
#ImageFeaturePlot(raw_list[["Rat.PDX"]], features = "ACTG2") # doesn't look helpful

RatSlide_expressed <- ProbeSet_lowOffTarget$gene_name[which(ProbeSet_lowOffTarget$gene_name %in% rownames(raw_list[["Rat.PDX"]][["Xenium"]]$counts))]

MouseSlide_expressed <- ProbeSet_lowOffTarget$gene_name[which(ProbeSet_lowOffTarget$gene_name %in% rownames(raw_list[["Mouse.PDX"]][["Xenium"]]$counts))]

lowOffTarget_expressed <- intersect(RatSlide_expressed, MouseSlide_expressed)

# subset cells to make plotting easier (as per https://github.com/satijalab/seurat/issues/2486 & https://satijalab.org/seurat/articles/essential_commands.html#subset-seurat-objects)
str(raw_list$Rat.PDX@meta.data)
Idents(raw_list[["Rat.PDX"]]) <- "Kmeans.2_10x"
Idents(raw_list[["Mouse.PDX"]]) <- "Kmeans.2_10x"

#save.image("./outputs/Robjs/SessionData.RData")

# tried to use "standard" normalization for plotting since SCTransform is very memory intensive but didn't seem to shift data much (still 0-6 scale)
# for(slide in samples){
#   print(slide)
#   raw_list[[slide]] <- NormalizeData(raw_list[[slide]])
#   gc()
# }

# from testing, use a more informative subset of genes (from testing full set); PDX marker source: https://www.nature.com/articles/s42003-021-02562-8#:~:text=In%20the%20PDX%20models%20we,1d%2C%20Table%201).
lowOffTarget_plot <- c("EPCAM", "GPC3", # cancer/PDX markers 
                       "LAMP3", "CXCL2", "SMIM24" # from 10x predictions
) # could also consider adding genes that have high cross-binding bet

for(slides in samples){
  # total cells (at current minimal filter) = ~ 2 million
  use.cells <- subset(raw_list[[slides]], downsample = 1000)
  
  ## Try heatmap since violin plot had issues but with scaled data
  DoHeatmap(use.cells, 
            features = lowOffTarget_plot,
            group.by = "PutativeOrigin",
            slot="counts", 
            angle = 0, hjust = 0.5, vjust = 0.25) + 
    scale_fill_gradientn(colours = c("lightgrey", "red")) +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold")) + 
    guides(color = "none") +
    labs(title = paste0("PDX markers and low cross-reactive targets \n across putative tissue origin (K=2)"),
         caption = paste0("SlideID: ", slides))
  
  ggsave(paste0(outputPath,'heatmap_geneMarkers_putativeTissueOrigin-',slides,'.pdf'),width=5,height=3.5)
  
  rm(use.cells)
  gc()
  
}

###########
## Filter out host tissue for each slide, based on k-means results from 10x processing
data_filt = list()
for(slides in samples){
  data_filt[[slides]] = subset(raw_list[[slides]], subset = PutativeOrigin == "tumor")
}
# warnings for not validating Seurat object/feilds

## save intermediate object after filtering
saveRDS(data_filt,file=paste0(results_folder,'Robjs/','filtered_kmeans-anno_xenium.obj.Rds'))

# cleanup raw data
rm(raw_list,meta_alt_mouse, meta_alt_rat, combined_meta,
   Mouse_Kmean2, HumantoMouse_scores, HumantoRat_scores,
   HumantoMouse_lowOffTarget, HumantoRat_lowOffTarget)
gc()

##save.image("./outputs/Robjs/SessionData.RData")
#######

# set identities back to sample ID
Idents(data_filt$Rat.PDX) <- "sample"
Idents(data_filt$Mouse.PDX) <- "sample"

# generate QC plots for after filtering for sharing
pdf(file = paste0(outputPath, 'postFilter_QC-violinPlots.pdf'),width=14,height=7)
for(slides in samples){
  VlnPlot(data_filt[[slides]], features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0) + labs(caption = paste0(slides,": QC metrics, post-filtering out putative host cells")) + 
    theme(plot.caption = element_text(hjust = 1.5))
}
dev.off()
## stopping point

## Run additional QC filtering ------------
# summarize nFeature and nCount per slide after filtering
tumor_filt_QC <- list()
for(slides in samples){
  tumor_filt_QC[[slides]] <- data_filt[[slides]]@meta.data %>% select(nCount_Xenium, nFeature_Xenium) %>% summary()
}
tumor_filt_QC ## check stats

### nCount - use min = 10, max = 800
## for rat slide, 1st quart = 65; max nCount = 2437 & 3rd quart = 143
## for mouse slide, 1st quart = 66; max nCount = 1142 & 3rd quart = 139
### nFeature - use min = 5
## similar stats between slides, min=1/3 & max=101/138 

data_filt_qc = list()
for(slides in samples){
  data_filt_qc[[slides]] = subset(data_filt[[slides]], subset = nCount_Xenium > 20 & nCount_Xenium < 800 & nFeature_Xenium > 5)
}
data_filt_qc

saveRDS(data_filt_qc,file=paste0(results_folder,'Robjs/','filtered_qc-kmeans-anno_xenium.obj.Rds'))

rm(data_filt, Rat_Kmean2); gc()

# re-load in filtered data if starting from Oct 7th Session Data 
#data_filt_qc <- readRDS("/nfs/mm-isilon/bioinfcore/ActiveProjects/Narla_gnrla_SP1-Xenium_damki_9638-KZ/outputs/Robjs/filtered_qc-kmeans-anno_xenium.obj.Rds")

## merge and process data --------------
data_merged <- merge(data_filt_qc[['Rat.PDX']], y = c(data_filt_qc[['Mouse.PDX']]), project = "Narla_PDX")
names(data_merged@images)[1] = 'Rat.PDX'
names(data_merged@images)[2] = 'Mouse.PDX'
rm(data_filt_qc); gc()

# check object
data_merged[["Xenium"]]$counts.1 # rat so counts.2 = mouse
#names(data_merged@assays$Xenium@layers) # leave as named so functions work w/o specifying layer names

## normalize data -----------
## decide which normalization is most appropriate - SCT vs "standard" lognorm

# # normalize and run dim reduction - SCT
# imaging based normalization paper recommended DEseq2, TMM, & SCTransform: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10491191/; However, having issues getting it to run so if continue to crashes
data_merged <- SCTransform(data_merged, assay = "Xenium",
                           variable.features.n = 250, # had error for second slide so will try to select # of genes to use as well
                           variable.features.rv.th = 1.5,
                           vst.flavor = "v2", 
                           conserve.memory = TRUE)
# once run - save object to file
saveRDS(data_merged,file=paste0(results_folder,'Robjs/','merged_SCT-kmeanFilt_xenium.obj.Rds'))

## to pick up analysis - load SCT normalized data
# data_merged <- readRDS("/nfs/mm-isilon/bioinfcore/ActiveProjects/Narla_gnrla_SP1-Xenium_damki_9638-KZ/outputs/Robjs/merged_SCT-kmeanFilt_xenium.obj.Rds")
#rm(data_filt, Rat_Kmean2); gc()

data_merged <- RunPCA(data_merged, reduction.name = "unintegrated.pca")
data_merged

# # normalize and run dim reduction - logNorm
# data_merged <- NormalizeData(data_merged, assay = "Xenium")
# data_merged <- FindVariableFeatures(data_merged, assay = "Xenium",
#                                    nfeatures = 30) # note - scale number of features to data
# data_merged <- ScaleData(data_merged)
# data_merged <- RunPCA(data_merged, reduction.name = "unintegrated.pca")

#save.image("./outputs/Robjs/SessionData1015.RData")

## check elbow plot?
ElbowPlot(data_merged, ndims = 50, reduction = "unintegrated.pca")

## check recommended # of PCs
# Estimate optimal PCs for clustering with a function -------------------------
optimal_pcs = function(so, reduction) {
  # quantitative check for number of PCs to include
  pct = so@reductions[[reduction]]@stdev / sum(so@reductions[[reduction]]@stdev) * 100
  cum = cumsum(pct)
  co1 = which(cum > 90 & pct < 5)[1]
  co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > .1), decreasing = T)[1] + 1
  pcs = min(co1, co2) 
  
  return(pcs)
}

# Apply function to our data
pcs = optimal_pcs(data_merged, 'unintegrated.pca')
pcs # recommendation is 11 but I could see using less

## check PCA plot?
head(data_merged@meta.data)
DimPlot(data_merged, reduction = "unintegrated.pca", group.by = "sample")

ggsave(file = paste0(outputPath, "qc_pca_plot_unintegrated_sct_sample.png"), width = 7, height = 6, units = "in")

## run unintegrated clustering and then DE comparisons since having issues with integration

data_merged = FindNeighbors(data_merged, dims = 1:pcs, reduction = 'unintegrated.pca')
data_merged = FindClusters(data_merged, resolution = 0.4, cluster.name = 'unintegrated.sct.clusters')

data_merged = RunUMAP(data_merged, dims = 1:pcs, reduction = 'unintegrated.pca', reduction.name = 'umap.unintegrated.pca')

# data_merged = RunTSNE(data_merged, dims = 1:pcs, reduction = 'unintegrated.pca', reduction.name = 'tsne.unintegrated.pca')

# plot UMAP
results_folder <- "./outputs/unintegrated_clustering/"

head(data_merged@meta.data)
pre_integration_umap_plot_clusters = DimPlot(data_merged, group.by = 'seurat_clusters', label = FALSE, reduction = 'umap.unintegrated.pca')
ggsave(filename = paste0(results_folder, 'umap_unintegrated_sct_clusters.png'), plot = pre_integration_umap_plot_clusters, width = 8, height = 6, units = 'in')

pre_integration_umap_plot_orig.ident = DimPlot(data_merged, group.by = 'seurat_clusters', split.by = 'sample', label = FALSE, reduction = 'umap.unintegrated.pca')
ggsave(filename = paste0(results_folder, 'umap_unintegrated_sct_splitByCondition.png'), plot = pre_integration_umap_plot_orig.ident, width = 12, height = 6, units = 'in')

rm(pre_integration_umap_plot_clusters, pre_integration_umap_plot_orig.ident); gc()

## save object to save UMAP reductions
saveRDS(data_merged,file=paste0(results_folder,'Robjs/','unintegrated-clustered_SCT-kmeanFilt_xenium.obj.Rds'))

# to pick up analysis, reload data 
# data_merged <- readRDS("/nfs/mm-isilon/bioinfcore/ActiveProjects/Narla_gnrla_SP1-Xenium_damki_9638-KZ/outputs/Robjs/unintegrated-clustered_SCT-kmeanFilt_xenium.obj.Rds")


## plot full sections per sample with existing cluster identities
## UMAPs with integrated clusters (ImageDimPlot cannot plot 2 samples together. Had to find an alternative way to plot. Weisheng posted this on github: https://github.com/satijalab/seurat/issues/8713) but current solution didn't work well
names(data_merged)
head(data_merged@meta.data)
unique(data_merged@images)

# full_sections_check <- ImageDimPlot(data_merged[['Mouse.PDX']],group.by='unintegrated.sct.clusters') | ImageDimPlot(raw_all[['Rat.PDX']],group.by='unintegrated.sct.clusters')
# 
# write out plots to file
# output to file
# minWidth = 15
# minHeight = 10
# OutputName <- paste0('SlideOverlays_Annotated_', pcs, "-PCs_", 'A|B')
# 
# png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
# print(full_sections_checkAB)
# dev.off()
# pdf(file=paste0(results_folder, OutputName,'.pdf'), width = minWidth, height = minHeight) # default = inches
# print(full_sections_checkAB)
# dev.off()
# ## add svg outputs
# svglite(file=paste0(results_folder, OutputName,'.svg'),  width = minWidth, height = minHeight) # default = inches
# print(full_sections_checkAB)
# dev.off()


### --- 

## Go back to subsetted data
# there are a few clusters that look to be primarily or exclusively present in one species. A cluster count table might be helpful to better see shared vs unique clusters before deciding which ones to combine (skip TSNE since so slow to run). 

head(data_merged@meta.data)
cell_count_per_cluster <- table(data_merged$unintegrated.sct.clusters, data_merged@meta.data$sample)
cell_count_per_cluster <- cell_count_per_cluster %>% as.data.frame() %>% pivot_wider(names_from = c("Var2"), values_from = "Freq")
colnames(cell_count_per_cluster) <- c("ClusterID", "Mouse.PDX_cell_counts", "Rat.PDX_cell_counts")
cell_count_per_cluster$ClusterID <- base::factor(cell_count_per_cluster$ClusterID, 
                                                 levels=c(0:13))
cell_count_per_cluster <- cell_count_per_cluster %>% arrange(ClusterID)

# write to file
write.csv(cell_count_per_cluster , paste0(results_folder, 'unintegrated_clusters_cell_counts_byCondition.csv'))

## identify which clusters are mostly shared 
shared_clusters <- cell_count_per_cluster %>% dplyr::filter(Mouse.PDX_cell_counts > 3000 & Rat.PDX_cell_counts > 3000)
shared_clusters <- as.numeric(shared_clusters$ClusterID)

## After limiting to shared clusters, run DE comparisons on all cells
head(data_merged@meta.data)
data_merged_shared_clusters <- subset(data_merged, subset = seurat_clusters %in% shared_clusters)
gc()

saveRDS(data_merged_shared_clusters,file=paste0('./outputs/Robjs/','unintegrated-sharedClustersSubset_SCT-kmeanFilt_xenium.obj.Rds'))

### ----- 


## Run DE comparisons, ignoring clusters
data_merged_shared_clusters
head(data_merged_shared_clusters@meta.data)
Idents(data_merged_shared_clusters) <- "sample"

## prep object, modifying slot names to allow functions to work per issue:
slot(object = data_merged_shared_clusters@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"SCT"
slot(object = data_merged_shared_clusters@assays$SCT@SCTModel.list[[1]], name="umi.assay")<-"SCT"
SCTResults(object=data_merged_shared_clusters, slot="umi.assay")

# Run prep SCT after changing names in object from Xenium to SCT
data_merged_shared_clusters <- PrepSCTFindMarkers(data_merged_shared_clusters)

# Try DE with MAST
# create thresholds
fcThresh <- 0.10
pvalThresh <- 0.01
DE_Rat_v_Mouse_MAST <- FindMarkers(data_merged_shared_clusters, 
                                   logfc.threshold = fcThresh, 
                                   slot="counts", test.use = "MAST", 
                                   ident.1 = "Rat.PDX", ident.2 = "Mouse.PDX")

saveRDS(DE_Rat_v_Mouse_MAST ,file=paste0('./outputs/Robjs/','DE-MAST.Xenium.Rds'))

# spot check results
DE_Rat_v_Mouse_MAST  %>% filter(p_val_adj < pvalThresh) %>% pull(p_val) %>% length() # 103 DE genes w/ FC 0.25; 158 w/ FC 0.10
OutTable <- tibble::rownames_to_column(DE_Rat_v_Mouse_MAST, "gene")
head(OutTable)
dim(OutTable)

## write out to file and deliver along with UMAP plot highlighting shared clusters
temp_path <- "./outputs/DE_shared_clusters/"
dir.create(temp_path, showWarnings = FALSE)
# Write to file
write.csv(OutTable,
          paste0(temp_path,"Shared-clusters_Rat-v-Mouse_MAST-DE-genes.csv"),
          quote = FALSE,
          row.names = FALSE)

# Write to file
write.table(OutTable, sep = "\t",
            paste0(temp_path,"Shared-clusters_Rat-v-Mouse_MAST-DE-genes",".txt"),
            quote = FALSE, 
            row.names = FALSE)


## Generate full DE results
DE_results.All <- FindMarkers(data_merged_shared_clusters,
                              logfc.threshold = 0, test.use = "MAST", 
                              ident.1 = "Rat.PDX", ident.2 = "Mouse.PDX")

res_tbl <- tibble::rownames_to_column(DE_results.All, "gene")
head(res_tbl); tail(res_tbl) # most genes are DE which could make functional enrichment challenging; GSEA as option with FC rank?
dim(res_tbl)

write.csv(res_tbl,
          paste0(temp_path,"Shared-clusters_Rat-v-Mouse_MAST-DE-allStats.csv"),
          quote = FALSE,
          row.names = FALSE)

# Write to file
write.table(res_tbl, sep = "\t",
            paste0(temp_path,"Shared-clusters_Rat-v-Mouse_MAST-DE-allStats",".txt"),
            quote = FALSE, 
            row.names = FALSE)

#### -----

## Next steps 
# add column to object to label shared and non-shared clusters 
# split up merged data into separate objects to see if plotting can work (vs tyring to read raw data back in, which has crashed sessions)
# plot UMAP and on slide visualizations for original and then shared and non shared clusters

# start with just shared cluster object since would need to re-read in `data_merged` object or load different session data
head(data_merged_shared_clusters)

# merged_split <- list()
# merged_split <- SplitObject(data_merged_shared_clusters, split.by = "sample") # warning for not validating centroid or FOV or Seurat objects
# head(merged_split$Rat.PDX)
# merged_split$Rat.PDX
# 
# full_sections_check <- ImageDimPlot(merged_split[['Mouse.PDX']],group.by='unintegrated.sct.clusters') | ImageDimPlot(merged_split[['Rat.PDX']],group.by='unintegrated.sct.clusters')
# full_sections_check # plot wroked but colors are not not very informative for subset data so clear and load back in full data set
# 
# rm(full_sections_check, merged_split); gc()
# rm(data_merged_shared_clusters); gc()

# load in merged data
data_merged <- readRDS("/nfs/mm-isilon/bioinfcore/ActiveProjects/Narla_gnrla_SP1-Xenium_damki_9638-KZ/outputs/Robjs/unintegrated-clustered_SCT-kmeanFilt_xenium.obj.Rds")

merged_split <- SplitObject(data_merged, split.by = "sample") 
head(merged_split$Rat.PDX)
rm(data_merged); gc()

# on slide plot for full cluster set
full_sections_check <- ImageDimPlot(merged_split[['Mouse.PDX']],group.by='unintegrated.sct.clusters') + labs(title="Mouse.PDX, cancer tissue") | ImageDimPlot(merged_split[['Rat.PDX']],group.by='unintegrated.sct.clusters') + labs(title="Rat.PDX, cancer tissue") ## "|" = patchwork operator
full_sections_check
# also output plots individually

# write out plots to file
# output to file
minWidth = 15
minHeight = 10
OutputName <- paste0('SlideOverlays_allClusters_', pcs, "-PCs_", 'RatvsMouse')

png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
print(full_sections_check)
dev.off()
## add svg outputs
svglite(file=paste0(results_folder, OutputName,'.svg'),  width = minWidth, height = minHeight) # default = inches
print(full_sections_check)
dev.off()
rm(full_sections_check); gc()

## add output for individual slides to improve formatting



## for each object, add column for which clusters are shared or not
print(shared_clusters)
for(i in names(merged_split)){
  print(i)
  merged_split[[i]]$shared_status <- "PerSlide"
  merged_split[[i]]@meta.data[which(merged_split[[i]]@meta.data$seurat_clusters %in% shared_clusters), "shared_status"] <- "Shared"
  merged_split[[i]]$shared_status <- factor(merged_split[[i]]$shared_status, levels = c("Shared", "PerSlide"))
}
#head(merged_split[[i]]@meta.data)

## save later
#saveRDS(merged_split,file=paste0(results_folder,'Robjs/','unintegrated-clustered_SCT-sharedClustersStatus_xenium.obj.Rds'))

full_sections_shared <- ImageDimPlot(merged_split[['Mouse.PDX']],group.by='shared_status') + labs(title="Mouse.PDX, cancer tissue") | ImageDimPlot(merged_split[['Rat.PDX']],group.by='shared_status') + labs(title="Rat.PDX, cancer tissue") ## "|" = patchwork operator
full_sections_shared

# output plots to file and then generate count table to compare
minWidth = 15
minHeight = 10
OutputName <- paste0('SlideOverlays_SharedClusters_', pcs, "-PCs_", 'RatvsMouse')

png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
print(full_sections_shared)
dev.off()
## add svg outputs
svglite(file=paste0(results_folder, OutputName,'.svg'),  width = minWidth, height = minHeight) # default = inches
print(full_sections_shared)
dev.off()
rm(full_sections_shared); gc()

## individual plot for mouse
OutputName <- paste0('SlideOverlays_SharedClusters_', pcs, "-PCs_", 'MouseOnly')
minWidth = 8
minHeight = 4
png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
ImageDimPlot(merged_split[['Mouse.PDX']],group.by='shared_status') + labs(title="Mouse.PDX, cancer tissue") 
dev.off()
## add svg outputs
svglite(file=paste0(results_folder, OutputName,'.svg'),  width = minWidth, height = minHeight) # default = inches
ImageDimPlot(merged_split[['Mouse.PDX']],group.by='shared_status') + labs(title="Mouse.PDX, cancer tissue") 
dev.off()

## repeat for Rat with same dimensions
OutputName <- paste0('SlideOverlays_SharedClusters_', pcs, "-PCs_", 'RatOnly')
png(file=paste0(results_folder, OutputName,'.png'), width = minWidth, height = minHeight, unit = "in", res = 300)
ImageDimPlot(merged_split[['Rat.PDX']],group.by='shared_status') + labs(title="Rat.PDX, cancer tissue")
dev.off()
## add svg outputs
svglite(file=paste0(results_folder, OutputName,'.svg'),  width = minWidth, height = minHeight) # default = inches
ImageDimPlot(merged_split[['Rat.PDX']],group.by='shared_status') + labs(title="Rat.PDX, cancer tissue")
dev.off()
#### -----

# generate summary table of number of shared vs non-shared
shared_counts <- list()
for(i in names(merged_split)){
  shared_counts[[i]] <- merged_split[[i]]$shared_status %>% table()
}
shared_counts

# re-organize list into table and add total number column
shared_table <- as.data.frame(shared_counts)
shared_table$Mouse.PDX.. <- NULL
colnames(shared_table)[1] <- "Status"

# add row of total counts for each slide
shared_table <- shared_table %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))

shared_table <- shared_table %>%
  dplyr::mutate(Total = rowSums(across(where(is.numeric)), na.rm=TRUE))

# output table to file
write.csv(shared_table , paste0(results_folder, 'unintegrated_shared_cell-counts.csv'))


## add step to run PrepSCT to ensure that resulting object is set up for DE comparisons
for(i in names(merged_split)){
  merged_split[[i]] <- PrepSCTFindMarkers(merged_split[[i]])
}

## generate markers for shared vs non-shared clusters within and between slides
head(merged_split$Rat.PDX@meta.data)
Idents(merged_split$Rat.PDX) <- "shared_status"
Idents(merged_split$Mouse.PDX) <- "shared_status"

# start with full set (no FC threshold)
fcThresh <- 0

DE_list <- list()
for(i in names(merged_split)){
  print(i)
  DE_list[[i]] <- FindMarkers(merged_split[[i]], 
                                   logfc.threshold = fcThresh,
                                   slot="counts", test.use = "MAST",
                                   ident.1 = "Shared", ident.2 = "PerSlide")
}

head(DE_list$Rat.PDX)
head(DE_list$Mouse.PDX)

## write out csv and text "full" versions to use for iPG
temp_path <- "./outputs/DE_between_clusters/"
dir.create(temp_path, showWarnings = FALSE)

for(i in names(DE_list)){
  
  res_tbl <- DE_list[[i]]
  res_tbl <- tibble::rownames_to_column(res_tbl, "gene")

  write.csv(res_tbl,
            paste0(temp_path,"Shared-v-Unique-",i,"_MAST-DE-allStats.csv"),
            quote = FALSE,
            row.names = FALSE)

  # Write to file
  write.table(res_tbl, sep = "\t",
              paste0(temp_path,"Shared-v-Unique-",i,"_MAST-DE-allStats",".txt"),
              quote = FALSE,
              row.names = FALSE)
}


## filter to DE results, plus generate summary table
fcThresh <- 0.25
pvalThresh <- 0.01

DE_subset <- list()
for(i in names(DE_list)){
  DE_subset[[i]] <- DE_list[[i]] %>% dplyr::filter(p_val_adj < pvalThresh & abs(avg_log2FC) > fcThresh)
  
  res_tbl <- tibble::rownames_to_column(DE_subset[[i]], "gene")
  write.csv(res_tbl,
            paste0(temp_path,"Shared-v-Unique-",i,"_MAST-DEonly.csv"),
            quote = FALSE,
            row.names = FALSE)
}

## summarize top DE from each comparison as dotplot
head(DE_subset$Rat.PDX)
head(DE_subset$Mouse.PDX)

top_10 <- rownames(head(DE_subset[[i]], n=10))

DotPlot(merged_split[[i]], top_10, group.by = "shared_status") # modify to be long format

## repeat DE directly between non-shared clusters to compare to shared cluster DE?
# this will require re-combining the Seurat object


#### -----
## save dated copy of Rdata, clearing previous versions as needed
version <- Sys.Date(); gc()
save.image(paste0("./outputs/Robjs/SessionData_",version,".RData"))

