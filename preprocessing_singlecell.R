
##This script is to preprocess the single cell RNA sequencing data from BD Rhapsody seurat objects for paired tissue and blood to study BCG response
##We had data with fastq files converted into downloadable seurat objects demultiplexed using BD Rhapsody pipeline on the seven bridges platform

##Preprocessing of BD Rhas


setwd("/projects/b1042/MeeksLab/Mairah/BD_data")

##Load libraries
library(ggrepel)
library(EnhancedVolcano)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(celldex)
library(dittoSeq)
library(SingleR)
library(dplyr)
library(biomaRt)
library(clusterProfiler)
library(cowplot)
library(ggpubr)
library("glmGamPoi")
library("ggplot2")
library("gridExtra")
library("parallel")
library(dplyr)
library(tidyverse)
library(SeuratWrappers)
library(monocle3)
library(decoupleR)
library(circlize)
library(ComplexHeatmap)
library(Matrix.utils)
library(DESeq2)
library(harmony)
library(clustree)
library(jsonlite)
library(DoubletFinder)
library(RColorBrewer)
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

##Read in the meta data file
metadata<-read.csv("/projects/b1042/MeeksLab/Mairah/BD_data/metadata_tumor.csv",sep=",",stringsAsFactors = FALSE)

###consider only those files where there is no chemotherapy
metadata_nochemo<-metadata[metadata$Chemotherapy=="No",]




##Functions for reading through the seurat objects, normalization and doublet removal

merged<-function(x) {
m<-readRDS(x)
m[["percent.mt"]] <- PercentageFeatureSet(m, pattern = "^MT-")
return(m)
} 


  
  
###First I need to just do a simple merge of the data and see what it looks like

BD_data<-lapply(metadata_nochemo$Files,merged)


###Let me remove the doublets and undetermined counts from the objects which are based on BD Rhapsody algorithm

doublet<-function(x){
	x<- SetIdent(x, value = "Sample_Tag")
	x<-subset(x,idents=c("Undetermined","Multiplet"),invert=TRUE)
}
BD_data[2:length(BD_data)]<-lapply(BD_data[2:length(BD_data)],doublet)

###Ok now remove the mitochondrial reads
##Let me make the volcano plots

plots<-function(x){
plot<- VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+ labs(x = unique(x@meta.data$Patient))
return(plot)}

p <- lapply(BD_data, FUN = plots)


#Remove mitochondrial reads

mito<-c(25,10,20,25,20,20,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25)

sub<-function(x,y){
	m=subset(x, subset = nFeature_RNA>200 & percent.mt < y)
	
}

BD_data_edited<-mapply(sub,BD_data,mito)


##Save the R object
save(BD_data_edited,file="BD_data_edited.RData")

##Add additional metadata

for (i in 1 : nrow(metadata_nochemo)){
Meta=BD_data_edited[[i]]@meta.data
## Add patient
Meta["Patient"] <- metadata_nochemo$Patient[i]
MetaTrim <- subset(Meta, select = c("Patient"))
BD_data_edited[[i]] <- AddMetaData(BD_data_edited[[i]], MetaTrim)

## Add Response
Meta["Response"] <- metadata_nochemo$Response[i]
MetaTrim <- subset(Meta, select = c("Response"))
BD_data_edited[[i]] <- AddMetaData(BD_data_edited[[i]], MetaTrim)

##Add Chemo
Meta["Chemotherapy"] <- metadata_nochemo$Chemotherapy[i]
MetaTrim <- subset(Meta, select = c("Chemotherapy"))
BD_data_edited[[i]] <- AddMetaData(BD_data_edited[[i]], MetaTrim)

##Add Stage
Meta["Stage"] <- metadata_nochemo$Stage[i]
MetaTrim <- subset(Meta, select = c("Stage"))
BD_data_edited[[i]] <- AddMetaData(BD_data_edited[[i]], MetaTrim)
}

###Since there are sample specific tags for paired tissue and blood, we need to add this info in, it was not in the metadata
tag<-function(x){
MetaTrim<-paste0(x@meta.data$Patient, '-',x@meta.data$Sample_Tag)
x <- AddMetaData(x, MetaTrim,col.name="Patient_sampletag")
}
BD_data_edited[2:length(BD_data_edited)]<-lapply(BD_data_edited[2:length(BD_data_edited)],tag)

###Add NA to sample tag for the non sample tagged sample
BD_data_edited[[1]]<- AddMetaData(BD_data_edited[[1]],paste0(BD_data_edited[[1]]@meta.data$Patient, '-',"NA"),col.name="Patient_sampletag")


##Let me subset each seurat object based on sample tag, as for each patient we have tissue and blood info for each of them 

processing<-function(x){
	x<- SetIdent(x, value = "Patient_sampletag")
   obj.list <- SplitObject(x, split.by = "Patient_sampletag")
   
}
BD_data_processed<-lapply(BD_data_edited,processing)
BD_data_processed<-unlist(BD_data_processed)

###I can now do doublet finder 

###Let me check which dataset has less than 50 cells I need to remove
cells<-c()
for (i in 1:length(BD_data_processed)){
	cells[i]= table(BD_data_processed[[i]]@meta.data$Patient)
	names(cells[i])<-names(table(BD_data_processed[[i]]@meta.data$Patient))
}

###Remove the tumor sample with 29 cells

BD_data_processed_1<-BD_data_processed[-9]


###For each of the datasets we can do sctransform separately and this means separating out the blood and the tumor before doing any of the transformation

sample <- function(x) {
x<-unlist(x)
x<-SCTransform(x, vars.to.regress = c("percent.mt"), vst.flavor = "v2")  
x<-RunPCA(x,features = VariableFeatures(x))
stdv <- x[["pca"]]@stdev
sum.stdv <- sum(x[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -  percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = T)[1] + 1
min.pc <- min(co1, co2)
min.pc
# finish pre-processing
x<-RunUMAP(x, dims = 1:min.pc) 
x<-FindNeighbors(x,dims = 1:min.pc)             
x<-FindClusters(x) 
x<-RunUMAP(x,dims = 1:min.pc)
x<-RunTSNE(x,dims = 1:min.pc) 
# pK identification (no ground-truth)
sweep.list <- paramSweep(x, PCs = 1:min.pc, sct=T)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)
# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
nExp_poi <- round(0.08*nrow(x@meta.data))
 # run DoubletFinder
x <- doubletFinder(seu = x,PCs = 1:min.pc,  pK = optimal.pk,nExp = nExp_poi,sct=TRUE)
metadata <- x@meta.data
colnames(metadata)[ncol(metadata)] <- "doublet_finder"
x@meta.data <- metadata   
# subset and save
singlets <- subset(x, doublet_finder == "Singlet")
return(singlets)
  }
  
##Normalise each of the seurat objects 
 BD_data_normalised<-lapply(BD_data_processed_1,sample)

##Merge the seurat objects 
BD<-merge(BD_data_normalised[[1]],y=BD_data_normalised[2:length(BD_data_normalised)],merge.data=TRUE)

###Use singleR to identify the immune cells 
sce<- as.SingleCellExperiment(BD,assay="RNA")
hpca.ref <- celldex::HumanPrimaryCellAtlasData()

annotations.final <- SingleR(test =sce,assay.type.test = 1,ref = hpca.ref,  labels = hpca.ref$label.main)

table(annotations.final$pruned.labels)

BD@meta.data$annotations.final <- annotations.final$pruned.labels
BD<- SetIdent(BD, value = "annotations.final")


###Data is already SCT normalised, now I would like to do clustering using PCA and UMAP
VariableFeatures(BD[["SCT"]]) <- rownames(BD[["SCT"]]@scale.data)


##Run PCA, UMAP and neighbours
BD<-BD %>%
RunPCA()

BD<-BD%>%RunUMAP(dim=1:15) %>%
FindNeighbors(dim=1:15)

BD<- FindClusters(BD, resolution = seq(from = 0.1, to = 1, by = 0.1))

###use clustree to identify the most stable clusters
clustree(BD, node_colour = "sc3_stability")

BD<- FindClusters(BD, resolution = 0.5)


###Now make a way of classifying sample type as tumor and blood with the sample tags


BD$type<-sapply(BD$Patient_sampletag,function(x){
	if(x=="P13-NA"||x=="P2-SampleTag03_hs"||x=="P8-SampleTag05_hs"||x=="P1-P-SampleTag07_hs"||x=="P14-SampleTag11_hs"||x=="P9-SampleTag09_hs"||x=="P10-SampleTag11_hs"||
	x=="P15-SampleTag01_hs"||x=="P4-SampleTag03_hs"||x=="P5-SampleTag05_hs"|| x=="P16-SampleTag09_hs"||x=="P11-SampleTag11_hs"||x=="P6-SampleTag05_hs"||x=="P7-SampleTag07_hs"||x=="P12-SampleTag09_hs"||x=="P17-SampleTag08_hs"||x=="P2-P-SampleTag01_hs"||x=="P4-P-SampleTag03_hs"||x=="P18-SampleTag07_hs"||x=="P19-SampleTag09_hs"){
		y="tissue"
	}else{y="blood"}
	})

BD$Patient_sampletype <- paste0(BD$Patient, '-', BD$type)

##Let us visualised the UMAPs

m<-DimPlot(BD, reduction = "umap",group.by="annotations.final")+theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20), 
    axis.title.x = element_text(size=20, face="bold", colour = "black"),
     axis.text.x = element_text(face="bold", colour = "black",size=20))
     


m1<-DimPlot(BD, reduction = "umap",group.by="Patient_sampletype")
m2<-DimPlot(BD, reduction = "umap",group.by="type")

pdf("Alldata.pdf")
m+m1
dev.off()

##
pdf("blood_tumor.pdf")
m2
dev.off()

