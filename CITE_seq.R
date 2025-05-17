##I had Ab-seq data for two patients along with the single cell RNA seq data which I will be using tyo confirm SingleR annotation


setwd("/projects/b1042/MeeksLab/Mairah/BD_data")
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

###Abseq analyses for T34 and T35 can we recapitulate it
###Let me look through the analyses for the last two samples
T34<-readRDS("/projects/b1042/MeeksLab/Mairah/BD_data/T34-B29_Seurat.rds")
###Separate the Ab-seq from the RNA expression
T35<-readRDS("/projects/b1042/MeeksLab/Mairah/BD_data/T35-B30_Seurat.rds")
###Separate the Ab-seq from the RNA expression

BD_data<-list(T34,T35)


##Merge data


ab_RNA_data<-merge(BD_data[[1]],BD_data[[2]])

ab_RNA_data<- SetIdent(ab_RNA_data, value = "Sample_Tag")
ab_RNA_data<-subset(ab_RNA_data,idents=c("Undetermined","Multiplet"),invert=TRUE) ##Remove multiplets identified by BD pipeline


##Subset Abseq and RNA data
Abseq <- subset(ab_RNA_data,features=rownames(ab_RNA_data)[(grep("pAbO", rownames(ab_RNA_data)))])
RNA <- subset(ab_RNA_data,features=rownames(ab_RNA_data)[(-grep("pAbO", rownames(ab_RNA_data)))])

 
# create new assay
RNA@assays[["AB"]] <- GetAssay(Abseq,assay = "RNA")


###Do immune cell annotations using SingleR with transcriptomic data

DefaultAssay(RNA) <- "RNA"

sce<- as.SingleCellExperiment(RNA,assay="RNA")

hpca.ref <- celldex::HumanPrimaryCellAtlasData()

annotations.final <- SingleR(test =sce,assay.type.test = 1,ref = hpca.ref,  labels = hpca.ref$label.main)

table(annotations.final$pruned.labels)

RNA@meta.data$annotations.final <- annotations.final$pruned.labels

##We can now look at SingleR on the transcriptomic and then group the protein and the mRNA assay according to it..


###Ok now I take the markers for the different immune cells and see if they are actually what Single R says they are and also see whether you get the same immune cells:


Idents(RNA)<-"annotations"
levels(RNA) <- c("B cells", "Monocyte", "Neutrophils","NK cells", "T cells")
DE <- FindAllMarkers(RNA, assay = "AB",random.state=123,max.cells.per.ident=100)

all.markers_sig<-DE[DE$p_val_adj<0.05,]
top10_markers <- as.data.frame(all.markers_sig %>% group_by(cluster) %>% top_n(10, wt= avg_log2FC))
top5_markers <- as.data.frame(all.markers_sig %>% group_by(cluster) %>% top_n(5, wt= avg_log2FC))

pdf("Ab-seq.pdf",height=14,width=12)
DoHeatmap(RNA,features=top5_markers$gene, size = 5,label=TRUE,group.by="annotations",group.colors=c("T cells"="#fd7f6f","Macrophage"="#7eb0d5","Monocyte"="#b2e061","NK cells"="#bd7ebe","Neutrophils"="#fdcce5","B cells"="#beb9db"))+ theme(legend.position="right", text = element_text(size = 14), axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        legend.title = element_text(face="bold")
)
dev.off()
