library(Seurat)

load("v5_sampIntegrated_cca_5kfeat_pc50_Sscale_G2Mscale_rmDubStrict_0203.Rdata")

save(subset(allSampsFilt, subset = li_cell_type == "CD4"), file = "CD4/v5_sampIntegrated_cca_5kfeat_pc50_Sscale_G2Mscale_rmDubStrict_0203.Rdata")
save(subset(allSampsFilt, subset = li_cell_type == "CD8"), file = "CD8/v5_sampIntegrated_cca_5kfeat_pc50_Sscale_G2Mscale_rmDubStrict_0203.Rdata")
