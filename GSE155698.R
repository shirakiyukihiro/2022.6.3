library(GEOquery)

setwd("//DESKTOP-KK97BAV/Raw_data")

#GSEのデータをそのまますべてダウンロードする
dir.create("GSE155698_pdac")
setwd("//DESKTOP-KK97BAV/Raw_data/GSE155698_pdac")

filePaths = getGEOSuppFiles("GSE155698")
filePaths

setwd("//DESKTOP-KK97BAV/Raw_data/GSE155698_pdac")

list.files(list.dirs(),full.names = T)

library(tidyverse)
list.files(list.dirs(),full.names = T)%>%str_subset("matrix.mtx.gz")%>%str_remove_all("matrix.mtx.gz")->list
list.files(list.dirs(),full.names = T)%>%str_subset("matrix.mtx$")%>%str_remove_all("matrix.mtx$")->list1
list.files(list.dirs(),full.names = T)%>%str_subset(".h5$")

library(Seurat)
library(hdf5r)
mtx1 <- lapply(list.files(list.dirs(),full.names = T)%>%str_subset(".h5"),Read10X_h5)
names(mtx1)<-list.files(list.dirs(),full.names = T)%>%str_subset(".h5")%>%str_remove_all("./GSE155698/|/filtered_.+$")

library(Matrix)
library(scater)
mtx2<-lapply(list, function(x){
  barcode.path <- paste0(x,"barcodes.tsv.gz")
  features.path <- paste0(x,"features.tsv.gz")
  matrix.path <- paste0(x,"matrix.mtx.gz")
  mat<-readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V2
  mat
})

names(mtx2)<-list%>%str_remove_all("./GSE155698/|/filtered_.+$")


mtx3<-lapply(list1, function(x){
  barcode.path <- paste0(x,"barcodes.tsv")
  features.path <- paste0(x,"genes.tsv")
  matrix.path <- paste0(x,"matrix.mtx")
  mat<-readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V2
  mat
})

names(mtx3)<-list1%>%str_remove_all("./GSE155698/|/filtered_.+$")

mtx<-c(mtx1,mtx2,mtx3)

saveRDS(mtx,"mtx.rds")







seurat<-lapply(1:43,function(x){CreateSeuratObject(mtx[[x]], project = names(mtx)[x], min.cells = 3, min.features = 200)})
#for (variable in names(seurat)){assign(names(seurat[variable]),seurat[[variable]])}

percent.mt <-lapply(seurat,function(x){PercentageFeatureSet(x, pattern = "^mt-|^MT-")})
percent.mt[[1]]
seurat[[1]][["nFeature_RNA"]]

for (variable in 1:43){seurat[[variable]][["percent.mt"]]<- percent.mt[[variable]]}

VlnPlot(seurat[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat<-lapply(seurat,function(x){subset(x,  subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 20)})

seurat<-lapply(seurat,NormalizeData)
seurat<-lapply(seurat,function(x){FindVariableFeatures(x,nfeatures = 2000)})
seurat<-lapply(seurat,function(x){ScaleData(x,features=rownames(x))})
seurat<-lapply(seurat,function(x){RunPCA(x,features=VariableFeatures(x))})
seurat<-lapply(seurat,function(x){JackStraw(x, num.replicate = 20)})
seurat<-lapply(seurat,function(x){ScoreJackStraw(x, dims = 1:20)})
library(PCAtools)#BiocManager::install("PCAtools")
seurat<-lapply(seurat,function(x){FindNeighbors(x, dims = seq(PCAtools::findElbowPoint(Stdev(x)^2)))})
seurat<-lapply(seurat,function(x){FindClusters(x, resolution = 0.1)})
seurat<-lapply(seurat,function(x){RunUMAP(x, dims = seq(PCAtools::findElbowPoint(Stdev(x)^2)))})


names(seurat)<-names(mtx)

saveRDS(seurat,"seurat.rds")

setwd("//DESKTOP-KK97BAV/Raw_data/GSE155698_pdac")
seurat<-readRDS("seurat.rds")

names(seurat)%>%str_detect("PDAC_TISSUE")


features <- SelectIntegrationFeatures(object.list = seurat[names(seurat)%>%str_detect("PDAC_TISSUE")])

anchors <- FindIntegrationAnchors(object.list = seurat[names(seurat)%>%str_detect("PDAC_TISSUE")], anchor.features = features)
# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)

saveRDS(combined,"combined.rds")

DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:6)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:6)
combined <- FindClusters(combined, resolution = 0.1)
# Visualization
p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(combined, reduction = "umap", split.by = "orig.ident", label=T)

# For performing differential expression after integration, we switch back to the original data
DefaultAssay(combined) <- "RNA"

saveRDS(combined,"combined.rds")

FeaturePlot(combined, features ="COL1A1",label = T)
FeaturePlot(combined, features ="CSPG4",label = T)














###############未実施

combined<-FindSubCluster(combined,"7", graph.name = "integrated_snn", subcluster.name = "sub.cluster",resolution = 0.1, algorithm = 1)
DimPlot(combined,  group.by = "sub.cluster", label = T)
Idents(combined)<-"sub.cluster"


all.markers <- FindAllMarkers(combined, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,"all.markers.csv",quote=F,sep=",",row.names=T,col.names=NA,append=F)
# 本ケースではマーカの中からトップ10を選んで描画する
library(tidyverse)
top30_all.markers <- all.markers %>% group_by(cluster) %>% top_n(n =30, wt=avg_log2FC)




combined <- RenameIdents(combined, `0` = "Cancer_cell", `1` = "Macrophage", `2` = "Endothelium", 
                         `3` = "Perivascular_cell", `4` = "T_cell", `5` = "B_cell", `6` = "??", 
                         `7` = "Endothelium?", `8` = "Cancer_cell?")

combined$cluster2<-Idents(combined)

DimPlot(combined, label = TRUE)

levels(factor(combined$orig.ident))
combined$orig.ident2<-as.factor(str_remove_all(combined$orig.ident,"_sample."))
combined$orig.ident2<-factor(combined$orig.ident2,levels = c("uninj","1dpi","3dpi","7dpi"))

DimPlot(combined, label = TRUE,split.by = "orig.ident2")
FeaturePlot(combined, features ="Islr",label = T,split.by = "orig.ident2")


markers.to.plot <-  c("Col1a1", "Gfap", "Nes", "Cd74", "Pdgfra", "Olig2", "Islr", "Acta2","Stra6",
                      "Ms4a6c","Tek","Cd14","Ptprc","Cd34","Pecam1","Cspg4","Cspg5","Nefl","Rgs5","Mki67")
markers.to.plot <-  c("Acan", "Vcan", "Ncan", "Cspg4", "Cspg5", "Smc3", "Bcan", "Cd44","Sdc2","Hspg2")


DotPlot(combined, features = markers.to.plot, cols = c("black","red","yellow","blue"), dot.scale = 8, split.by = "orig.ident2") +  RotatedAxis()
VlnPlot(combined, features = "Islr", split.by = "orig.ident2", group.by = "cluster",  pt.size = 0, combine = FALSE)

table(combined$orig.ident2,combined$cluster2)
prop.table(table(combined$orig.ident2,combined$cluster2), margin = 1)
write.table(prop.table(table(combined$orig.ident2,combined$cluster2), margin = 1),"table.csv",quote=F,sep=",",row.names=T,col.names=NA,append=F)



DimPlot(combined, label = TRUE)
Idents(combined) <- factor(Idents(combined), levels = c( "Fibroblast","CD14_macrophage","Cd74_macrophage", "Cd83_macrophage",
                                                         "Ccr7_dendritic_macrophage", "Gzmb_macrophage","Ly6d_dendritic_macrophage",
                                                         "Csf1_1_neutrophil", "Csf1_3_neutrophil", "Cxcr2_neutrophil", "endothelium_5", "Grama_NK_T", 
                                                         "Gzmc_T", "Cd8a_T", "Klra1_NK_T","CD4_T", "Mki67_T","Cd79a_B" ,"Krt18_8" ))
markers.to.plot <-  c("Cd3d", "Col1a1", "Cd79a", "Cd74", "Cd8a", "Krt18", "Islr", "Acta2","Pdgfra","Stra6",
                      "Ms4a6c","Csf1","Cd14","Cxcr2","Cd34","Pecam1","Gzma","Gzmb",
                      "Klra1","Krt8","Cd4","Cd83","Ccr7","Mki67","Ly6d")
DotPlot(combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + 
  RotatedAxis()

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
Fibro <- subset(combined, idents = "Fibroblast")
Idents(Fibro) <- "stim"
avg.Fibro <- as.data.frame(log1p(AverageExpression(Fibro, verbose = FALSE)$RNA))
avg.Fibro$gene <- rownames(avg.Fibro)

cd8.T <- subset(combined, idents = "Cd8a_T")
Idents(cd8.T) <- "stim"
avg.cd8.T <- as.data.frame(log1p(AverageExpression(cd8.T, verbose = FALSE)$RNA))
avg.cd8.T$gene <- rownames(avg.cd8.T)

genes.to.label = c("Islr", "Cd8a")
p1 <- ggplot(avg.Fibro, aes(con, AM80)) + geom_point() + ggtitle("Fibroblast")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd8.T, aes(con, AM80)) + geom_point() + ggtitle("Cd8a_T")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2

combined$celltype.stim <- paste(Idents(combined), combined$stim, sep = "_")
combined$celltype <- Idents(combined)

Idents(combined) <- "celltype.stim"
levels(Idents(combined))
fibro.response <- FindMarkers(combined, ident.1 = "Fibroblast_con", ident.2 = "Fibroblast_AM80", verbose = FALSE)
head(fibro.response, n = 15)

ec.response <- FindMarkers(combined, ident.1 = "endothelium_5_con", ident.2 = "endothelium_5_AM80", verbose = FALSE)
head(ec.response, n = 15)

b.response <- FindMarkers(combined, ident.1 = "Cd79a_B_con", ident.2 = "Cd79a_B_AM80", verbose = FALSE)
head(b.response, n = 15)


FeaturePlot(combined, features = c("Islr", "Mxra8", "Stra6"), split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))

plots <- VlnPlot(combined, features = c("Islr", "Mxra8", "Stra6"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
library(patchwork)
wrap_plots(plots = plots, ncol = 1)

saveRDS(combined,"combined.rds")


Idents(combined) <- "celltype"

pdf("seurat.pdf",height =25,width = 40)
DimPlot(combined, reduction = "umap", split.by = "stim", label=T)
FeaturePlot(combined, features = markers.to.plot, min.cutoff = "q9")
DotPlot(combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + 
  RotatedAxis()
FeaturePlot(combined, features = c("Islr", "Mxra8", "Stra6","Acta2"), split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))
plots <- VlnPlot(combined, features = c("Islr", "Mxra8", "Stra6","Acta2"), split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()


Idents(combined) <- "celltype.stim"

levels(Idents(combined))[grep("_T_con",levels(Idents(combined)))] <- rep("T_con")
levels(Idents(combined))[grep("_T_AM80",levels(Idents(combined)))] <- rep("T_AM80")
levels(Idents(combined))[grep("_macrophage_con",levels(Idents(combined)))] <- rep("macrophage_con")
levels(Idents(combined))[grep("_macrophage_AM80",levels(Idents(combined)))] <- rep("macrophage_AM80")
levels(Idents(combined))[grep("_neutrophil_con",levels(Idents(combined)))] <- rep("neutrophil_con")
levels(Idents(combined))[grep("_neutrophil_AM80",levels(Idents(combined)))] <- rep("neutrophil_AM80")

combined$celltype2<-Idents(combined) 
levels(Idents(combined))
T.response <- FindMarkers(combined, ident.1 = "T_con", ident.2 = "T_AM80", verbose = FALSE)
head(T.response, n = 15)

macrophage.response <- FindMarkers(combined, ident.1 = "macrophage_con", ident.2 = "macrophage_AM80", verbose = FALSE)
head(macrophage.response, n = 15)

neutrophil.response <- FindMarkers(combined, ident.1 = "neutrophil_con", ident.2 = "neutrophil_AM80", verbose = FALSE)
head(neutrophil.response, n = 15)

table(Idents(combined))


write.table(fibro.response,"fibro.response.csv",quote=F,sep=",",row.names = T,col.names = NA,append = F)
write.table(T.response,"T.response.csv",quote=F,sep=",",row.names = T,col.names = NA,append = F)
write.table(macrophage.response,"macrophage.response.csv",quote=F,sep=",",row.names = T,col.names = NA,append = F)
write.table(b.response,"b.response.csv",quote=F,sep=",",row.names = T,col.names = NA,append = F)
write.table(ec.response,"ec.response.csv",quote=F,sep=",",row.names = T,col.names = NA,append = F)
write.table(neutrophil.response,"neutrophil.response.csv",quote=F,sep=",",row.names = T,col.names = NA,append = F)
write.table(table(Idents(combined)),"table.csv",quote=F,sep=",",row.names = T,col.names = NA,append = F)

library(EnhancedVolcano)#BiocManager::install('EnhancedVolcano')

EnhancedVolcano(fibro.response,
                lab = rownames(fibro.response),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'con versus AM80',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

pdf("volcano.pdf")
EnhancedVolcano(fibro.response,
                lab = rownames(fibro.response),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Fibroblast con versus AM80',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)
EnhancedVolcano(macrophage.response,
                lab = rownames(macrophage.response),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Macrophage con versus AM80',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)
EnhancedVolcano(T.response,
                lab = rownames(T.response),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'T-cell con versus AM80',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)
dev.off()

ifnb.list <- lapply(X = c(sc_con,sc_am80), FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = c(sc_con,sc_am80), anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT", 
                                  anchor.features = features)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)
p1 <- DimPlot(combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE, 
              repel = TRUE)
p1 + p2


##################################################
sce<-lapply(mtx, function(x){
  SingleCellExperiment(assays = list(counts = x))
})


mt.gene<-lapply(sce, function(x){grep("^mt-",rownames(x))})

sce[[names(mtx)[1]]]
mt.gene[[names(mtx)[1]]]

sce<-lapply(names(sce), function(x){addPerCellQC(sce[[x]], subsets=list(Mito=mt.gene[[x]]))})
names(sce)<-names(mtx)

plotColData(sce$.Mutant1, x="sum", y="detected")
plotColData(sce$.Mutant1, x="sum", y="subsets_Mito_percent")

sce<-lapply(sce, function(x){x$discard <- x$detected <200|x$subsets_Mito_percent > 25
x})

sce_<-lapply(sce, function(x){discard<-x$discard
x[,!discard]})

library(scran)
sce_<-lapply(sce_,function(x){
  x<-computeSumFactors(x, cluster=quickCluster(x, BPPARAM=BiocParallel::MulticoreParam()), min.mean=0.1)
  x<- logNormCounts(x)
  x
})

saveRDS(sce_,"sce_.rds")

genes<-lapply(sce_, rownames)
commongenes<-genes[[1]]
for (variable in genes) {
  commongenes<-intersect(commongenes,variable)
}
sce__<-lapply(sce_, function(x){
  x[match(commongenes,rownames(x)),]
})

library(batchelor)
#値にNAがあると計算できないので注意
set.seed(1000101001)
mnn<-fastMNN(sce__, k=50, BSPARAM=BiocSingular::RandomParam(deferred=TRUE))

gc(reset = TRUE)
gc(reset = TRUE)


library(scater)
library(scran)
snn.gr <- buildSNNGraph(mnn, k=50, use.dimred="corrected")
clusters.mnn <- igraph::cluster_louvain(snn.gr)$membership

set.seed(1100101001)
mnn <- runUMAP(mnn, dimred="corrected")

mnn$batch <- factor(mnn$batch)

colLabels(mnn) <- factor(clusters.mnn)
assay(mnn, "logcounts") <- assay(mnn, "reconstructed")

plotUMAP(mnn, colour_by="batch")

saveRDS(mnn,"mnn.rds")

plotUMAP(mnn,colour_by="Islr",text_by="label")
plotUMAP(mnn,colour_by="Col1a1",text_by="label")

plotUMAP(mnn[,mnn$batch%>%str_detect("WT")],colour_by="Islr",text_by="label")
plotUMAP(mnn[,mnn$batch%>%str_detect("Mutant")],colour_by="Islr",text_by="label")

plotExpression(mnn[,mnn$batch%>%str_detect("WT")],x="label",features = "Islr",colour_by = "label")#11,14
plotExpression(mnn[,mnn$batch%>%str_detect("Mutant")],x="label",features = "Islr",colour_by = "label")#11,14


plotExpression(mnn[,mnn$batch%>%str_detect("WT")],x="label",features = "Col1a1",colour_by = "label")#11,14
plotExpression(mnn[,mnn$batch%>%str_detect("Mutant")],x="label",features = "Col1a1",colour_by = "label")#11,14














