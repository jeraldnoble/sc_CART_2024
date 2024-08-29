library(Seurat)
library(BPCells)
library(SeuratDisk)
library(SeuratData)
library(dplyr)
library(cowplot)
library(patchwork)
library(stringr)
library(ggplot2)
library(formattable)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(gridExtra)
library(Azimuth)
library(stringr)
library(Polychrome)
options(scipen=999)
library(celldex)
library(scRNAseq)
library(SingleR)
library(SingleCellExperiment)
library(SeuratWrappers)
library(data.table)
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")

load("/share/lab_freeman/MM_CAR_2022/script/cycle.rda")

load("allSamps_abecShort.Rdata")
### get QC metrics
allSamps$percent.mt = PercentageFeatureSet(object = allSamps, pattern = "^MT-")
allSamps$log10GenesPerUMI = log10(allSamps$nFeature_RNA) / log10(allSamps$nCount_RNA)
allSamps$barcode = rownames(allSamps@meta.data)

load("doublet_labels.Rdata")
allSamps@meta.data = plyr::join(allSamps@meta.data, dubDf, by = "barcode")
rownames(allSamps@meta.data) = allSamps@meta.data$barcode
allSamps@meta.data$doublet = "no"
#allSamps@meta.data$doublet[allSamps@meta.data$doublet_finder == "Doublet" & allSamps@meta.data$scrublet_doublet == "True"] = "yes"
allSamps@meta.data$doublet[allSamps@meta.data$doublet_finder == "Doublet"] = "yes"



resMeta = read.table("/share/lab_freeman/MM_CAR_2022/data/response_batch_meta.tsv", sep = "\t", header = T)
allSamps@meta.data = plyr::join(allSamps@meta.data, resMeta, by = "orig.ident")
rownames(allSamps@meta.data) = allSamps@meta.data$barcode

#allSampsFilt = subset(allSamps, subset = (nCount_RNA >= 1000) & (nFeature_RNA >= 500) & (nFeature_RNA <= 7000) & (log10GenesPerUMI >= 0.8) & (percent.mt < 15) & (doublet == "no"))
allSampsFilt = subset(allSamps, subset = (nCount_RNA >= 1000) & (nFeature_RNA >= 500) & (log10GenesPerUMI >= 0.8) & (percent.mt < 15) & (doublet == "no"))
#allSampsFilt = subset(allSampsFilt, subset = (type == "IP") & (batch == 3))

### remove genes expressed in less than 10 cells
counts = LayerData(object = allSampsFilt, layer = "counts")
nonzero = counts > 0
keep_genes = Matrix::rowSums(nonzero) >= 10
filtered_counts = counts[keep_genes,]
allSampsFilt = CreateSeuratObject(filtered_counts, meta.data = allSampsFilt@meta.data)

### require cd3 expression
#counts = as.data.frame(t(as.matrix(LayerData(object = allSampsFilt, layer = "counts"))))
#allSampsFilt$cd3 = "no"
#allSampsFilt$cd3[counts$CD3D > 0 | counts$CD3E > 0 | counts$CD3G > 0] = "yes"
#allSampsFilt = subset(allSampsFilt, subset = cd3 == "yes")
allSampsFilt = RunAzimuth(allSampsFilt, reference = "pbmcref")

allSampsFilt = NormalizeData(allSampsFilt)
allSampsFilt = FindVariableFeatures(allSampsFilt, selection.method = "vst", nfeatures = 3000)
varGenes = VariableFeatures(allSampsFilt)
tmp = which(!grepl("^TRA.+|^TRB.+|^TRD.+|^TRG.+|^IGL.+|Abecma|Carvykti", varGenes))
VariableFeatures(allSampsFilt) = varGenes[tmp]
allSampsFilt = ScaleData(allSampsFilt)
allSampsFilt = CellCycleScoring(allSampsFilt, g2m.features=g2m_genes, s.features = s_genes)

allSampsFilt[["RNA"]] = split(allSampsFilt[["RNA"]], f = allSampsFilt$orig.ident)
allSampsFilt = NormalizeData(allSampsFilt)
allSampsFilt = FindVariableFeatures(allSampsFilt, selection.method = "vst", nfeatures = 3000)
varGenes = VariableFeatures(allSampsFilt)
tmp = which(!grepl("^TRA.+|^TRB.+|^TRD.+|^TRG.+|^IGL.+|Abecma|Carvykti", varGenes))
VariableFeatures(allSampsFilt) = varGenes[tmp]

allSampsFilt = ScaleData(allSampsFilt, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score", "nCount_RNA"))
#allSampsFilt = ScaleData(allSampsFilt)
allSampsFilt = RunPCA(allSampsFilt)
allSampsFilt = FindNeighbors(allSampsFilt, dims = 1:50, reduction="pca")
allSampsFilt = FindClusters(allSampsFilt, resolution = 1, cluster.name="unint_clusters")
allSampsFilt = RunUMAP(allSampsFilt, dims = 1:50, reduction = "pca", reduction.name= "umap.unintegrated")

#allSampsFilt = IntegrateLayers(allSampsFilt, method = CCAIntegration, new.reduction = "integrated.cca", dims = 1:50, orig.reduction = "pca")
#allSampsFilt = IntegrateLayers(allSampsFilt, method = RPCAIntegration, new.reduction = "integrated.rpca", dims = 1:50, orig.reduction = "pca")
allSampsFilt = IntegrateLayers(allSampsFilt, method = HarmonyIntegration, new.reduction = "integrated.harmony", dims = 1:50, orig.reduction = "pca")

#allSampsFilt = FindNeighbors(allSampsFilt, reduction = "integrated.cca", dims = 1:50)
#allSampsFilt = FindClusters(allSampsFilt, resolution = c(0.5, 0.75, 1))
#allSampsFilt = RunUMAP(allSampsFilt, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")

#allSampsFilt = FindNeighbors(allSampsFilt, reduction = "integrated.rpca", dims = 1:50)
#allSampsFilt = FindClusters(allSampsFilt, resolution = c(0.5), cluster.name = "rpca_clusters.0.5")
#allSampsFilt = FindClusters(allSampsFilt, resolution = c(0.75), cluster.name = "rpca_clusters.0.75")
#allSampsFilt = FindClusters(allSampsFilt, resolution = c(1), cluster.name = "rpca_clusters.1")
#allSampsFilt = RunUMAP(allSampsFilt, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca")

allSampsFilt = FindNeighbors(allSampsFilt, reduction = "integrated.harmony", dims = 1:50)
allSampsFilt = FindClusters(allSampsFilt, resolution = c(0.5), cluster.name = "harmony_clusters.0.5")
allSampsFilt = FindClusters(allSampsFilt, resolution = c(0.75), cluster.name = "harmony_clusters.0.75")
allSampsFilt = FindClusters(allSampsFilt, resolution = c(1), cluster.name = "harmony_clusters.1")
allSampsFilt = RunUMAP(allSampsFilt, reduction = "integrated.harmony", dims = 1:50, reduction.name = "umap.harmony")


allSampsFilt = JoinLayers(allSampsFilt)

marker_list = fread("/share/lab_freeman/MM_CAR_2022/data/marker_list.txt")
marker_list = as.list(marker_list)
for(i in 1:length(marker_list)){marker_list[[i]] = marker_list[[i]][which(marker_list[[i]] != "")]}

other_markers = fread("/share/lab_freeman/MM_CAR_2022/data/other_markers.txt")
other_markers = as.list(other_markers)
for(i in 1:length(other_markers)){other_markers[[i]] = other_markers[[i]][which(other_markers[[i]] != "")]}

lh_deconv = fread("/share/lab_freeman/MM_CAR_2022/data/linghua_deconvolution_markers.txt")
lh_deconv = as.list(lh_deconv)
for(i in 1:length(lh_deconv)){lh_deconv[[i]] = lh_deconv[[i]][which(lh_deconv[[i]] != "")]}

lh_cd4 = fread("/share/lab_freeman/MM_CAR_2022/data/linghua_cd4_markers.txt")
lh_cd4 = as.list(lh_cd4)
for(i in 1:length(lh_cd4)){lh_cd4[[i]] = lh_cd4[[i]][which(lh_cd4[[i]] != "")]}

lh_cd8 = fread("/share/lab_freeman/MM_CAR_2022/data/linghua_cd8_markers.txt")
lh_cd8 = as.list(lh_cd8)
for(i in 1:length(lh_cd8)){lh_cd8[[i]] = lh_cd8[[i]][which(lh_cd8[[i]] != "")]}

anderson_markers = fread("/share/lab_freeman/MM_CAR_2022/data/anderson_markers.txt")
anderson_markers = as.list(anderson_markers)
for(i in 1:length(anderson_markers)){anderson_markers[[i]] = anderson_markers[[i]][which(anderson_markers[[i]] != "")]}

allSampsFilt = AddModuleScore(allSampsFilt, features = marker_list, name = names(marker_list))
allSampsFilt = AddModuleScore(allSampsFilt, features = other_markers, name = names(other_markers))
allSampsFilt = AddModuleScore(allSampsFilt, features = lh_deconv, name = names(lh_deconv))
allSampsFilt = AddModuleScore(allSampsFilt, features = lh_cd4, name = names(lh_cd4))
allSampsFilt = AddModuleScore(allSampsFilt, features = lh_cd8, name = names(lh_cd8))
allSampsFilt = AddModuleScore(allSampsFilt, features = anderson_markers, name = names(anderson_markers))


save(allSampsFilt, file = "v5_sampIntegrated_cca_5kfeat_pc50_Sscale_G2Mscale_rmDub_0321.Rdata")

