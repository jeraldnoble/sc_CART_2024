library(Seurat)
library(DoubletFinder)
library(stringr)

samp = Read10X(data.dir = "PATH")
sampName = str_split_1("PATH", '/')[9]
seuratObj = CreateSeuratObject(samp, project = sampName, min.cells = 10, min.features = 500, names.field = 1, names.delim = "_")

run_doublet_finder = function(seuratObj, expDoublets, sampName)
{

	seuratObj = NormalizeData(seuratObj)
 	seuratObj = FindVariableFeatures(seuratObj, verbose = F)
	seuratObj = ScaleData(seuratObj, verbose = F)
	seuratObj = RunPCA(seuratObj, npcs = 30, verbose = F)

	stdv = seuratObj[['pca']]@stdev
	sum.stdv = sum(seuratObj[['pca']]@stdev)
	percent.stdv = (stdv / sum.stdv) * 100
	cumulative = cumsum(percent.stdv)
	co1 = which(cumulative > 90 & percent.stdv < 5)[1]
	co2 = sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
		      percent.stdv[2:length(percent.stdv)]) > 0.1), 
	     decreasing = T)[1] + 1
	min.pc = min(co1, co2)
	seuratObj = RunUMAP(seuratObj, dims = 1:min.pc, verbose = F)
	seuratObj = FindNeighbors(seuratObj, reduction = "pca", dims = 1:min.pc)
	seuratObj = FindClusters(seuratObj, resolution =c(0.4,0.6,0.8,1.0,1.2))

	sweep.res = paramSweep(seuratObj, PCs = 1:min.pc)
	sweep.stats = summarizeSweep(sweep.res, GT = F)
	bcmvn = find.pK(sweep.stats)
	bcmvn.max = bcmvn[which.max(bcmvn$BCmetric),]
	pk = as.numeric(levels(bcmvn.max$pK))[bcmvn.max$pK]

	annotations = seuratObj@meta.data$seurat_clusters
	homotypic.prop = modelHomotypic(annotations)
	nExp_poi = round(pk * nrow(seuratObj@meta.data))
	nExp_poi.adj = round(nExp_poi * (1 - homotypic.prop))

	### get expected doublets for the sample
	expDoublets = expDoublets$doublets[which(expDoublets$sample == sampName)]
	nExp_man = round(expDoublets * nrow(seuratObj@meta.data))

	seuratObj = doubletFinder(seuratObj, PCs = 1:min.pc, pK = pk, nExp = nExp_poi.adj)

	colnames(seuratObj@meta.data)[grepl("DF.classification", colnames(seuratObj@meta.data))] = 'doublet_finder'
	#eval(call("<-", as.name(sampName), seuratObj))

	### change here for filt/unfilt data
	outFile = paste(sampName, "_doubletFinder.Rdata", sep = "")
	#save(seuratObj, file = outFile)
	#get(sampName) %>% save(file = outFile)
	return(seuratObj)
}


expDoublets = read.table("/share/lab_freeman/MM_CAR_2022/script/expected_doublets_by_sample.txt", sep = "\t", header = T)
seuratObj = run_doublet_finder(seuratObj, expDoublets, sampName)
seuratObj$barcode = paste0(seuratObj$orig.ident, "_", rownames(seuratObj@meta.data))

outFile = paste0(sampName, "_doublet_finder.Rdata")
save(seuratObj, file = outFile)

