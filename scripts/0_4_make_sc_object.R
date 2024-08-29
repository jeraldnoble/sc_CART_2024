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
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")

dirs = list.dirs("/share/lab_freeman/MM_CAR_2022/data/all_data/scrnaseq/cellranger_outs", recursive=T)
tmp = which(grepl("filtered_feature_bc_matrix", dirs))
dirs = dirs[tmp]

### get sample names
sampNames = str_split(dirs, '/', simplify = T)[,9]

### make list of 10x objects
allSamps = sapply(dirs, function(i){
        print(i)
        samp = Read10X(data.dir = i)
        sampName = str_split_1(i, '/')[9]
        print(sampName)
        colnames(samp) = paste0(sampName, "_", colnames(samp))
        samp
})

### collapse list
allSamps.data = do.call("cbind", allSamps)

### make seurat object from list
allSamps = CreateSeuratObject(allSamps.data, project = "abecma", min.cells = 10, min.features = 500, names.field = 1, names.delim = "_")
save(allSamps, file = "allSamps_abecShort.Rdata")

