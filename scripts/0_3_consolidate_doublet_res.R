library(stringr)

dubDf = data.frame()

dirs = list.dirs("/share/lab_freeman/MM_CAR_2022/data/all_data/scrnaseq/cellranger_outs", recursive=T)
tmp = which(grepl("filtered_feature_bc_matrix", dirs))
dirs = dirs[tmp]

sampNames = str_split(dirs, '/', simplify = T)[,9]

for(i in 1:length(dirs))
{
	dir = dirs[i]
	samp = sampNames[i]
	scrub = paste0(dir, "/", samp, "_scrublet_res.tsv")
	scrub = read.table(scrub, header = T, sep = "\t")
	load(paste0(dir, "/", samp, "_doublet_finder.Rdata"))
	df =  merge(seuratObj@meta.data[c("doublet_finder", "barcode")], scrub[c("scrublet_doublet", "barcode")], by = "barcode")
	dubDf = rbind(dubDf, df)
}

save(dubDf, file = "doublet_labels.Rdata")
