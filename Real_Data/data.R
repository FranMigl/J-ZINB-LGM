# see https://nbisweden.github.io/workshop-scRNAseq/labs/seurat/seurat_01_qc.html for quality control 
library(scater)
library(MouseGastrulationData)
library(scuttle)


AtlasSampleMetadata

sce5 <- EmbryoAtlasData(samples = 5)
sce5

sce18 <- EmbryoAtlasData(samples = 18)
sce18


singlets <- which(!(colData(sce5)$doublet | colData(sce5)$stripped))

plot(
  x = reducedDim(sce5, "umap")[singlets, 1],
  y = reducedDim(sce5, "umap")[singlets, 2],
  col = EmbryoCelltypeColours[colData(sce5)$celltype[singlets]],
  pch = 19,
  xaxt = "n", yaxt = "n",
  xlab = "UMAP1", ylab = "UMAP2"
)

sce5=sce5[,singlets]

singlets <- which(!(colData(sce18)$doublet | colData(sce18)$stripped))

plot(
  x = reducedDim(sce18, "umap")[singlets, 1],
  y = reducedDim(sce18, "umap")[singlets, 2],
  col = EmbryoCelltypeColours[colData(sce18)$celltype[singlets]],
  pch = 19,
  xaxt = "n", yaxt = "n",
  xlab = "UMAP1", ylab = "UMAP2"
)


sce18=sce18[,singlets]

is.mito <- as.logical(grepl("^mt-", rownames(sce5)))
summary(is.mito)
wrkrs=15

rowData(sce5)$ave.counts <- calculateAverage(sce5, exprs_values = "counts")
to.keep <- rowData(sce5)$ave.counts > 0
sce5 <- sce5[to.keep,]
summary(to.keep)
dim(sce5)

rowData(sce18)$ave.counts <- calculateAverage(sce18, exprs_values = "counts")
to.keep <- rowData(sce18)$ave.counts > 0
sce18 <- sce18[to.keep,]
summary(to.keep)
dim(sce18)

table(colData(sce5)$celltype)
table(colData(sce18)$celltype)

ind=which(colData(sce18)$celltype=="Primitive Streak")
sce18_p=sce18[,ind]

ind=which(colData(sce5)$celltype=="Primitive Streak")
sce5_p=sce5[,ind]

to.keep <- rowData(sce18_p)$ave.counts > 5
sce18_p <- sce18_p[to.keep,]
summary(to.keep)
dim(sce18_p)


to.keep <- rowData(sce5_p)$ave.counts > 5
sce5_p <- sce5_p[to.keep,]
summary(to.keep)
dim(sce5_p)

a=intersect(rownames(sce5_p),rownames(sce18_p))
ind=match(a,rownames(sce5_p))
sce5_p=sce5_p[ind,]

ind=match(a,rownames(sce18_p))
sce18_p=sce18_p[ind,]

sample5=as.matrix(counts(sce5_p))
sample18=as.matrix(counts(sce18_p))

rownames(sample5)=rowData(sce5_p)$SYMBOL
rownames(sample18)=rowData(sce18_p)$SYMBOL

sample5=t(sample5)
sample18=t(sample18)
ind=grep("mt-",colnames(sample5))
sample5=sample5[,-ind]
sample18=sample18[,-ind]

#togliere i geni Rp e malat1, calcolare il total count senza questi geni



IntegMat=list(sample5,sample18)
saveRDS(IntegMat, file = "MouseData.rds")
saveRDS(sample5_offset, file = "offset_sample5.rds")
saveRDS(sample18_offset, file = "offset_sample18.rds")
