library(sceasy)
library(reticulate)
loompy <- reticulate::import('loompy')

#sceasy::convertFormat("ana.h5ad", from = "anndata", to = "seurat", main_layer = "scale.data", outFile = "adata.rds")

allSamps = readRDS("adata.rds")

df = read.table("classify_cell_types.tsv", sep = "\t", header = T, row.names =1)

df$barcode = rownames(df)

df$sample = gsub(".*-1-", "", df$barcode)
df$newbar = gsub("*-CAR.*", "", df$barcode)
df$newbar = paste0(df$sample, "_", df$newbar)

meta = allSamps@meta.data
meta$barcode = rownames(meta)
meta$cell_type[meta$leiden == 10] = "Mono"
meta$cell_type = gsub(" T", "", meta$cell_type)
colnames(meta) = c("barcode", "li_cell_type")
write.table(meta, "li_cell_types.tsv", sep = "\t", row.names = F, quote = F)


load("../v5_sampIntegrated_cca_5kfeat_pc50_Sscale_G2Mscale_rmDubStrict_0226.Rdata")
allSampsFilt@meta.data = plyr::join(allSampsFilt@meta.data, meta, by = "barcode")
rownames(allSampsFilt@meta.data) = allSampsFilt@meta.data$barcode
save(allSampsFilt, file = "../v5_sampIntegrated_cca_5kfeat_pc50_Sscale_G2Mscale_rmDubStrict_0226.Rdata")
