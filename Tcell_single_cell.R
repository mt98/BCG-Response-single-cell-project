###Looking at T cells within the processed immune cell data to see how they are associated with BCG Response

###We looked at several questions in our analyses
###Separated out CD4 and CD8 T cells for each of the analyses
##These questions were asked separately for tissue and blood to look at BCG response

#1)-Are the T cells different between BCG exposed and BCG naive samples

#2)-Are the T cells different between BCG responders and non-responders 

##We also looked at how are the T cell different between tissue and blood


###Load libraries
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
library(scales)
library(harmony)
library(clustree)
library(jsonlite)
library(DoubletFinder)
library(RColorBrewer)
library(condiments)
library(speckle)
library(slingshot)
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)


###Load the processed object
load("BD_1_FINAL.RData")

###Subset out the T cells from the immune cell popularion
Idents(BD_1)<-"annotations.final"
BD_Tcells<-subset(BD_1,idents="T_cells")


#Let us look at CD4 T cells first to answer these questions:



###Separate CD4 T cells first based on cell type experimental values which are from the BD Rhapsody pipeline
Idents(BD_Tcells)<-"Cell_Type_Experimental"
BD_Tcells_CD4<-subset(BD_Tcells,idents=c("T_CD4_memory","T_CD4_naive"))

##Split T cells based on the sample tag for each patient which is separating issue and blood samples 
BD_Tcells_CD4.list <- SplitObject(BD_Tcells_CD4, split.by="Patient_sampletag")


##Let us define variable features to integrate on it will take all the variable features and see which are common among all datasets to integrate

var.features <- SelectIntegrationFeatures(object.list = BD_Tcells_CD4.list , nfeatures = 3000)

##Let us merge the CD4 T cells for each of the patient samples 
BD_Tcells_CD4.sct <- merge(x = BD_Tcells_CD4.list[[1]], y =BD_Tcells_CD4.list[2:length(BD_Tcells_CD4.list)], merge.data=TRUE)

##Select the variable features for the CD4 T cells 
VariableFeatures(BD_Tcells_CD4.sct ) <- var.features

##Clustering of CD4 T cells after running harmony based on different patients and type: tissue or blood
BD_Tcells_CD4.sct <- BD_Tcells_CD4.sct %>%
  RunPCA(verbose = FALSE) %>%
  RunHarmony(assay.use = "SCT", group.by.vars = c("Patient", "type")) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

clustree(BD_Tcells_CD4.sct, node_colour = "sc3_stability")

save(BD_Tcells_CD4.sct,file="BD_Tcells_CD4.sct.RData") ##save CD4 T cell object

##Select resolution for the CD4 T cell clustering based on clustree algorithm
BD_Tcells_CD4.sct_1<- FindClusters(BD_Tcells_CD4.sct, resolution =  0.5)



##BCG naive is naive group and the responders and non responders after BCG are exposed group, we define these categories

BD_Tcells_CD4.sct_1@meta.data$exposure<-ifelse(Tissue_CD4@meta.data$Response=="BCG_naive","BCG_naive","BCG_exposed")

##Feature Plot to ensure that CD4 T cells have high expression of CD4 as expected, it is a sanity check

FeaturePlot(BD_Tcells_CD4.sct_1,features=c("CD4"))+theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),axis.text.x = elem
ent_text(face="bold", colour = "black",size=20))+guides(x = axis, y = axis) +theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +scale_y_continuous(breaks = NULL)


##Visualise the CD4 T cell seurat clusters in UMAP
DimPlot(BD_Tcells_CD4.sct_1, reduction = "umap",group.by="seurat_clusters") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("CD4")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)


###Visualise the CD4 T cells tissue and blood

DimPlot(BD_Tcells_CD4.sct_1, reduction = "umap",group.by="type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
guides(x = axis, y = axis) +
ggtitle("")+
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)+labs(color = "Type") 


###Use markers to identify clusters

CD<-c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G","CD44")
Treg<-c("FOXP3","IL2RA","IKZF2")
Naive<-c("TCF7","SELL","LEF1","CCR7")
Exhausted<-c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4","LAYN","ENTPD1")
Cytotoxic<-c("GZMA","GZMB","GZMK","GNLY","PRF1","NKG7")
Effector_memory_CD4<-c("CRIP1","ITGB1","S100A11")
Th17<-c("IL1R1","AHR","CSF2","KLRB1","IL17A","MAF","CCR6","NFKBIZ","IL21R","IL21","IL12RB1","RORA","RORC","STAT3")
Th2<-c("CXCR4","CCR8","PTGDR2","HAVCR1","GATA3","STAT6","IL4","IL5","IL13","AREG")
Th1<-c("KLRD1","IFNGR1","CXCR3","CXCR6","CCR1","CCR5","STAT1","STAT4","TBX21","TNF","LTA","IFNG","IL2")
Effector_memory_CD8<-c("CCL4","CD69")



###Let us identify the CD4 T cells using seurat clusters
row_split_T =c(rep("CD",7),rep("Treg",3),rep("Naive_Tcell",4),rep("Exhausted",7),rep("Cytotoxic",6),rep("Effector_memory_CD4",3),rep("Th17",14),rep("Th2",10),rep("Effector_memory_CD8",2),rep("Th1",13))
row_split_T = factor(row_split_T,levels = c("CD","Treg","Naive_Tcell","Exhausted","Cytotoxic","Effector_memory_CD4","Th17","Th2","Effector_memory_CD8","Th1"))

### plot heatmap
features_Tcell<-c(CD,Treg,Naive,Exhausted,Cytotoxic,Effector_memory_CD4,Th17,Th2,Effector_memory_CD8,Th1)


### set colnames order
plot_ord_T <- c("0","1","2","3","4","5","6","7","8","9")

plot_ord_T <- rownames(table(BD_Tcells_CD4.sct_1$Tcell_type))
data.plot_T <- Heat_Dot_data(object=BD_Tcells_CD4.sct_1,features=features_Tcell,group.by="seurat_clusters") ###Function defined and loaded 

exp.mat_T <- data.plot_T %>% dplyr::select(features.plot,id,avg.exp.scaled) %>% spread(id,avg.exp.scaled)
rownames(exp.mat_T) <- exp.mat_T$features.plot
exp.mat_T$features.plot <- NULL
exp.mat_T <- exp.mat_T[,plot_ord_T]
per.mat_T <- data.plot_T %>% dplyr::select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat_T) <- per.mat_T$features.plot
per.mat_T$features.plot <- NULL
per.mat_T <- per.mat_T[,plot_ord_T]/100
min(exp.mat_T);max(exp.mat_T)



# left annotation
annot_T = c("CD","Treg","Naive","Exhausted","Cytotoxic","Effector_memory_CD4","Th17","Th2","Effector_memory_CD8")

pdf("Tcell_cluster_CD4_type.pdf",height=14,width=10)
Heatmap(as.matrix(exp.mat_T),cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,
        column_names_side = "bottom",row_names_side = "left",
        row_split = row_split_T,
        row_gap = unit(3.5, "mm"),
        column_gap = unit(3.5, "mm"), width = ncol(exp.mat_T)*unit(4, "mm"),
    height = nrow(exp.mat_T)*unit(4, "mm"),
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()


###cluster names based on marker genes

clusters = BD_Tcells_CD4.sct_1$seurat_clusters
Tcell_type<-case_when( clusters %in% 0 ~ 'CD4_Th17-like',
               clusters %in% 1 ~ 'CD4_Naive_1',
               clusters %in% 2 ~ 'CD4_Th17-like',
               clusters %in% 3 ~ 'CD4_TNF',
               clusters %in%4 ~ 'Treg',
               clusters %in% 5 ~ 'CD4_TNF' ,
                      clusters %in% 6 ~ 'CD4_TEMRA',
                       clusters %in% 7 ~ 'CD4_STAT4' ,
                       clusters %in% 8 ~ 'CD4_Naive_2',
                        clusters %in% 9 ~ 'CD4_EM'
            )

BD_Tcells_CD4.sct_1$Tcell_type<-Tcell_type

save(BD_Tcells_CD4.sct_1,file="BD_Tcells_CD4.sct_1.RData")



##UMAP plot for the CD4 T cell clusters



DimPlot(BD_Tcells_CD4.sct_1, reduction = "umap",group.by="Tcell_type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("CD4")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)

LabelClusters(plot = t, id = 'Tcell_type',fontface = "bold",size=5)

##UMAP plot for the CD4 T cell clusters by type: tissue or blood


t1<-DimPlot(BD_Tcells_CD4.sct_1, reduction = "umap",group.by="type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("CD4")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)



##UMAP plot for the CD4 T cell clusters
DimPlot(BD_Tcells_CD4.sct_1, reduction = "umap",group.by="Tcell_type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("CD4")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)

LabelClusters(plot = t, id = 'Tcell_type',fontface = "bold",size=5)


t1<-DimPlot(BD_Tcells_CD4.sct_1, reduction = "umap",group.by="type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("CD4")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)



###UMAPs for tissue and blood CD4 T cells separately

Idents(BD_Tcells_CD4.sct_1)<-"type"
Tissue_CD4=subset(BD_Tcells_CD4.sct_1,idents=c("tissue"))


t2<-DimPlot(Tissue_CD4, reduction = "umap",group.by="Tcell_type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("Tissue CD4")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)

LabelClusters(plot = t5, id = 'Tcell_type',fontface = "bold")

Idents(BD_Tcells_CD4.sct_1)<-"type"
Blood_CD4=subset(BD_Tcells_CD4.sct_1,idents=c("blood"))


t3<-DimPlot(Blood_CD4, reduction = "umap",group.by="Tcell_type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("Blood CD4")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)
LabelClusters(plot = t6, id = 'Tcell_type',fontface = "bold")






##We want to make a comparison of tissue and blood CD4 T cell clusters
##Comparison plots for CD4 clusters in tissue and blood
BD_Tcells_CD4.sct_1$type <- factor(x = BD_Tcells_CD4.sct_1$type, levels = c("tissue","blood"))

t4<-ggplot(BD_Tcells_CD4.sct_1@meta.data, aes(x=type, fill=Tcell_type)) + geom_bar(position="fill")+xlab("Type")+ylab("Fraction")+theme(text = element_text(size = 18,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=18), # bold
    axis.title.x = element_text(size=18, face="bold", colour = "black"),
     axis.text.x = element_text(face="bold", colour = "black",size=18), axis.ticks.x=element_blank())+theme(plot.title=element_text(hjust=0.5))+ theme(strip.background = element_blank(),
   strip.text.x = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(fill="CD4")+ scale_x_discrete(labels=c("tissue" = "Tissue", "blood" = "Blood"))


##We are comparing two CD4 clusters naive 1 and naive 2 using find markers function, as one of them is enriched in blood the other in tissue

##Compare naive 1 and  naive 2 CD4 clusters

Idents(BD_Tcells_CD4.sct_1)<-"Tcell_type"
DE_naive1_naive2 <- FindMarkers(BD_Tcells_CD4.sct_1 , assay = "SCT",recorrect_umi=FALSE,max.cells.per.ident=500,random.seed=123,ident.1="CD4_Naive_1",ident.2="CD4_Naive_2")

DE_naive1_naive2_sig_high<-DE_naive1_naive2%>%filter(avg_log2FC>0.58&p_val_adj<0.05)%>%arrange(desc(avg_log2FC))

DE_naive1_naive2_sig_low<-DE_naive1_naive2%>%filter(avg_log2FC<-0.58&p_val_adj<0.05)%>%arrange(avg_log2FC)

###Now we look at BCG responders and non-responders which are categories after BCG exposure by making stacked barplots for visualisation

Idents(Tissue_CD4)<-"Response"
Tissue_CD4_response<-subset(Tissue_CD4, idents=c("BCG_naive"),invert=T) ##Remove BCG naive group

n<-ggplot(Tissue_CD4_response@meta.data, aes(x=Response, fill=Tcell_type)) + geom_bar(position="fill")+xlab("")+ylab("Fraction")+theme(text = element_text(size = 18,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=18), # bold
    axis.title.x = element_text(size=18, face="bold", colour = "black"),
     axis.text.x = element_text(face="bold", colour = "black",size=18), axis.ticks.x=element_blank())+theme(plot.title=element_text(hjust=0.5))+ theme(strip.background = element_blank(),
   strip.text.x = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(fill="CD4")


###Now we look at BCG exposed and BCG naive CD4 T cells by making stacked barplots for visualisation

Tissue_CD4$exposure<- factor(x = Tissue_CD4$exposure, levels = c("BCG_naive","BCG_exposed"))

n1<-ggplot(Tissue_CD4@meta.data, aes(x=exposure, fill=Tcell_type)) + geom_bar(position="fill")+xlab("")+ylab("Fraction")+theme(text = element_text(size = 18,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=18), # bold
    axis.title.x = element_text(size=18, face="bold", colour = "black"),
     axis.text.x = element_text(face="bold", colour = "black",size=18), axis.ticks.x=element_blank())+theme(plot.title=element_text(hjust=0.5))+ theme(strip.background = element_blank(),
   strip.text.x = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(fill="CD4")
Tissue_blood_comparison<-propeller(
  x = BD_Tcells_CD4.sct_1,
  clusters = BD_Tcells_CD4.sct_1@meta.data$Tcell_type,
  sample = BD_Tcells_CD4.sct_1@meta.data$Patient_sampletype,
  group = BD_Tcells_CD4.sct_1@meta.data$type,
  transform = "asin"
)


###Use propeller function to compare cluster proportions of CD4 T cells between tissue and blood 


grp <- rep(c("tissue","blood"), c(20, 20))

props <- getTransformedProps(BD_Tcells_CD4.sct_1@meta.data$Tcell_type,BD_Tcells_CD4.sct_1@meta.data$Patient_sampletype, transform="asin")

df=data.frame(t(props$TransformedProps[,c( "P1-P-tissue","P10-tissue","P11-tissue","P12-tissue","P13-tissue","P14-tissue","P15-tissue","P16-tissue","P17-tissue","P18-tissue", 
 "P19-tissue","P2-P-tissue","P2-tissue","P4-P-tissue","P4-tissue","P5-tissue","P6-tissue","P7-tissue","P8-tissue","P9-tissue", "P1-P-blood", "P10-blood" , "P11-blood" , "P12-blood" ,
  "P14-blood",  "P15-blood" , "P16-blood" , "P17-blood",  "P18-blood" , "P19-blood",  "P2-blood" ,  "P2-P-blood", "P3-blood",  "P4-blood", "P4-P-blood", "P5-blood","P6-blood","P7-blood","P8-blood"   ,"P9-blood" )]),grp)


l5<-ggplot(df, aes(x=clusters, y=Freq)) +
  geom_boxplot(aes(fill=grp)) +geom_point(position=position_dodge(width=0.75),aes(group=grp)) +
labs(y="Cell-type proportion per sample")+labs(x="")+ labs(fill="Type")+theme(plot.title = element_text(hjust = 0.5,face="bold",size = 16))+theme(axis.text.x = element_text(size=12,color="black"),axis.text.y = element_text(size=16,color="black")) +theme(axis.title.y = element_text(size = 16))+ggtitle("")+theme(plot.title = element_text(hjust = 0.5,face="bold",size=16))+theme(strip.background = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(fill="")+geom_point(x=1,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=1.1,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=1.2,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=3,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=3.1,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=3.2,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=3.3,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=4,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=4.1,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=4.2,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)+geom_point(x=4.3,y=0.9,shape = "*", size=7 ,color="black",stroke = 2)


###

###Use propeller function to compare cluster proportions of  tissue CD4 T cells between BCG responders and BCG non-responders

##Proportion analyses
##Tissue

Idents(Tissue_CD4)<-"Patient"

###subset out patient samples with low counts which could bias the proportion analyses
Tissue_CD4_sub<-subset(Tissue_CD4,idents=c("P2","P2-P","P4-P"),invert=T)

##subset further

Idents(Tissue_CD4_sub)<-"Response"
Tissue_CD4_R<-subset(Tissue_CD4_sub,idents=c("BCG_responsive","BCG_unresponsive"))


output_tcell.asin <- propeller(clusters=Tissue_CD4_R@meta.data$Tcell_type, sample=Tissue_CD4_R@meta.data$Patient, group=Tissue_CD4_R@meta.data$Response, transform="asin") 

props <- getTransformedProps(Tissue_CD4_R@meta.data$Tcell_type,Tissue_CD4_R@meta.data$Patient, transform="asin")



###Use propeller function to compare cluster proportions of  tissue CD4 T cells between BCG exposed and BCG naive

Tissue_CD4_sub@meta.data$exposure<-ifelse(Tissue_CD4_sub@meta.data$Response=="BCG_naive","BCG_naive","BCG_exposed")

output_t.asin <- propeller(clusters=Tissue_CD4_sub@meta.data$Tcell_type, sample=Tissue_CD4_sub@meta.data$Patient, group=Tissue_CD4_sub@meta.data$exposure, transform="asin") 


###Proportion analyses blood

###Use propeller function to compare cluster proportions of blood CD4 T cells between BCG responders and BCG non-responders

Idents(Blood_CD4)<-"Patient"

##small nos of cell remove
Blood_CD4_sub<-subset(Blood_CD4,idents=c("P2","P14","P16"),invert=T)

##subset further

Idents(Blood_CD4_sub)<-"Response"
Blood_CD4_R<-subset(Blood_CD4_sub,idents=c("BCG_responsive","BCG_unresponsive"))



output_tcell.asin <- propeller(clusters=Blood_CD4_R@meta.data$Tcell_type, sample=Blood_CD4_R@meta.data$Patient, 
group=Blood_CD4_R@meta.data$Response, transform="asin") 

###Use propeller function to compare cluster proportions of  blood CD4 T cells between BCG exposed and BCG naive

Blood_CD4_sub@meta.data$exposure<-ifelse(Blood_CD4_sub@meta.data$Response=="BCG_naive","BCG_naive","BCG_exposed")

output_t.asin <- propeller(clusters=Blood_CD4_sub@meta.data$Tcell_type, sample=Blood_CD4_sub@meta.data$Patient, group=Blood_CD4_sub@meta.data$exposure, transform="asin") 
props <- getTransformedProps(Blood_CD4_sub@meta.data$Tcell_type,Blood_CD4_sub@meta.data$Patient, transform="asin")

###Differential expression analyses for tissue CD4 T cells between BCG responders and non-responders


###Use Differential expression analyses with MAST using Sex as confounding variable for Responsive vs Unresponsive

Sex<-sapply(Tissue_CD4@meta.data$Patient,function(x){
if(x=="P1-P"|x=="P15"|x=="P19"){y="Female"} else{y="Male"}	
})
Tissue_CD4@meta.data$Sex<-Sex


Idents(Tissue_CD4)<-"Response"
DE_CD4<-FindMarkers(Tissue_CD4, assay="SCT",recorrect_umi=FALSE,ident.1="BCG_responsive",ident.2="BCG_unresponsive",test.use="MAST",latent.vars="Sex")

DE_CD4 <- DE_CD4 %>%
  mutate(Genes = case_when(avg_log2FC>= 0.58 & p_val_adj <= 0.05 ~ "Up in BCG_responsive",
                             avg_log2FC <= -0.58 &  p_val_adj <= 0.05 ~ "Up in BCG_unresponsive",
                               TRUE ~ "NS"))


sig_genes_CD4<-DE_CD4%>% filter(Genes!="NS")


sig_genes_CD4_up<-sig_genes_CD4%>% dplyr::filter(Genes=="Up in BCG_responsive")%>% arrange(desc(avg_log2FC))

sig_genes_CD4_down<-sig_genes_CD4%>% dplyr::filter(Genes=="Up in BCG_unresponsive")%>% arrange(desc(avg_log2FC))

write.csv(sig_genes_CD4,file="CD4_R_UR_tissue.csv",quote=F)


###Look at enriched pathways for tissue CD4 T cells 

go_CD4<-enrichGO(
rownames(sig_genes_CD4_up),
'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")


###what are the pathways
barplot(go_CD4, showCategory=10) 



go_CD4_down<-enrichGO(
rownames(sig_genes_CD4_down),
'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")


###what are the pathways
barplot(go_CD4_down, showCategory=10) 



###Use MAST for differential expression analyses with Sex as confounding exposed vs naive in tissue CD4 T cells

Idents(Tissue_CD4)<-"exposure"
DE_CD4_e <- FindMarkers(Tissue_CD4, assay = "SCT",recorrect_umi=FALSE,ident.1="BCG_exposed",ident.2="BCG_naive",test.use="MAST",latent.vars="Sex")

DE_CD4_e <- DE_CD4_e %>%
  mutate(Genes = case_when(avg_log2FC>= 0.58 & p_val_adj <= 0.05 ~ "Up in BCG exposed",
                             avg_log2FC <= -0.58 &  p_val_adj <= 0.05 ~ "Up in BCG naive",
                               TRUE ~ "NS"))


sig_genes_CD4_e<-DE_CD4_e%>% filter(Genes!="NS")


sig_genes_CD4_e_up<-sig_genes_CD4_e%>% dplyr::filter(Genes=="Up in BCG exposed")%>% arrange(desc(avg_log2FC))

sig_genes_CD4_e_down<-sig_genes_CD4_e%>% dplyr::filter(Genes=="Up in BCG naive")%>% arrange(desc(avg_log2FC))



###Enrichment analyses
go_CD4<-enrichGO(
rownames(sig_genes_CD4_e_up),
'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")


###what are the pathways
barplot(go_CD4, showCategory=10) 



go_CD4_down<-enrichGO(
rownames(sig_genes_CD4_down),
'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")


###what are the pathways
barplot(go_CD4_down, showCategory=10) 

###Repeat same analyses for blood 


Sex<-sapply(Blood_CD4@meta.data$Patient,function(x){
if(x=="P1-P"|x=="P15"|x=="P19"){y="Female"} else{y="Male"}	
})
Blood_CD4@meta.data$Sex<-Sex

Idents(Blood_CD4)<-"Response"
DE_CD4_blood <- FindMarkers(Blood_CD4, assay = "SCT",recorrect_umi=FALSE,ident.1="BCG_responsive",ident.2="BCG_unresponsive",test.use="MAST",latent.vars="Sex")

DE_CD4_blood <- DE_CD4_blood %>%
  mutate(Genes = case_when(avg_log2FC>= 0.58 & p_val_adj <= 0.05 ~ "Up in BCG responsive",
                             avg_log2FC <= -0.58 &  p_val_adj <= 0.05 ~ "Up in BCG unresponsive",
                               TRUE ~ "NS"))


sig_genes_CD4_blood<-DE_CD4_blood%>% filter(Genes!="NS")


sig_genes_CD4_up_blood<-sig_genes_CD4_blood%>% dplyr::filter(Genes=="Up in BCG responsive")%>% arrange(desc(avg_log2FC))

sig_genes_CD4_down_blood<-sig_genes_CD4_blood%>% dplyr::filter(Genes=="Up in BCG unresponsive")%>% arrange(desc(avg_log2FC))


write.csv(sig_genes_CD4_blood,file="sig_genes_CD4_blood.csv")

##Enrichment

go_CD4_blood<-enrichGO(
rownames(sig_genes_CD4_down_blood),
'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")

barplot(go_CD4_blood, showCategory=10) 

##Exposed vs Naive for blood


Idents(Blood_CD4)<-"exposure"
DE_CD4_blood_e <- FindMarkers(Blood_CD4, assay = "SCT",recorrect_umi=FALSE,ident.1="BCG_exposed",ident.2="BCG_naive")

DE_CD4_blood_e <- DE_CD4_blood_e %>%
  mutate(Genes = case_when(avg_log2FC>= 0.58 & p_val_adj <= 0.05 ~ "Up in BCG exposed",
                             avg_log2FC <= -0.58 &  p_val_adj <= 0.05 ~ "Up in BCG naive",
                               TRUE ~ "NS"))


sig_genes_CD4_blood_e<-DE_CD4_blood_e%>% filter(Genes!="NS")


sig_genes_CD4_up_blood_e<-sig_genes_CD4_blood_e%>% dplyr::filter(Genes=="Up in BCG exposed")%>% arrange(desc(avg_log2FC))

sig_genes_CD4_down_blood_e<-sig_genes_CD4_blood_e%>% dplyr::filter(Genes=="Up in BCG naive")%>% arrange(desc(avg_log2FC))

write.csv(sig_genes_CD4_blood_e,file="sig_genes_CD4_blood_exposed.csv")


###Trajectory analyses for the CD4 T cell cluster

  #Making plots more siminmar to ggplot outputs of Seurat
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(16, ...)
    return(pal[cut(cell_vars, breaks = 16)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
#We need color palettes Leiden clusters. These would be the same colors seen in the Seurat plots.

cell_colors_clust <- cell_pal(BD_Tcells_CD4.sct_1$Tcell_type, hue_pal())



set.seed(1)
lineages <- as.SlingshotDataSet(getLineages(
  data           = BD_Tcells_CD4.sct_1_subset@reductions$umap@cell.embeddings,
  clusterLabels  = BD_Tcells_CD4.sct_1_subset$Tcell_type,
  dist.method    = "mnn", 
  start.clus     = "CD4_Naive_2"
)) # define where to START the trajectories


plot(BD_Tcells_CD4.sct_1_subset@reductions$umap@cell.embeddings, col = cell_colors_clust[BD_Tcells_CD4.sct_1_subset$Tcell_type], cex = .5, pch = 16)

lines(lineages, lwd = 1, col = "black", cex = 2)


curves <- as.SlingshotDataSet(getCurves(
  data          = lineages,
  thresh        = 1e-1,
  stretch       = 1e-1,
  allow.breaks  = F,
  approx_points = 100
))
plot(BD_Tcells_CD4.sct_1_subset@reductions$umap@cell.embeddings, col = cell_colors_clust[BD_Tcells_CD4.sct_1_subset$Tcell_type], cex = .5, pch = 16)
lines(curves, lwd = 2, col = "black")


pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)

x <- rowMeans(pseudotime)
x <- x / max(x)
o <- order(x)

plot(BD_Tcells_CD4.sct_1_subset@reductions$umap@cell.embeddings[o, ],
    main = paste0("pseudotime"), pch = 16, cex = 0.4, axes = F, xlab = "", ylab = "",
    col = colorRampPalette(c("grey70", "orange3", "firebrick", "purple4"))(99)[x[o] * 98 + 1]
  )
lines(curves, lwd = 2, col = "black")

col = colorRampPalette(c("grey70", "orange3", "firebrick", "purple4"))(99)[x[o] * 98 + 1]

plot(BD_Tcells_CD4.sct_1_subset@reductions$umap@cell.embeddings[o, ],
    main = paste0("pseudotime"), pch = 16, cex = 0.4, axes = F, xlab = "", ylab = "",
    col = colorRampPalette(c("grey70", "orange3", "firebrick", "purple4"))(99)[x[o] * 98 + 1])
lines(curves, lwd = 2, col = "black")


###Look at density plots within each trajectory to compare conditions
Idents(BD_Tcells_CD4.sct_1_subset)<-"Response"
Responder_cells=colnames(subset(BD_Tcells_CD4.sct_1_subset,idents=c("BCG_responsive")))
Nonresponder_cells=colnames(subset(BD_Tcells_CD4.sct_1_subset,idents=c("BCG_unresponsive")))

ds <- list(Responders = density(pseudotime[Responder_cells,1]),
           NonResponders = density(pseudotime[Nonresponder_cells,1]))
xlim <- range(c(ds$Responders$x, ds$NonResponders$x))
ylim <- range(c(ds$Responders$y, ds$NonResponders$y))

ds1 <- list(Responders = density(pseudotime[Responder_cells,2]),
           NonResponders = density(pseudotime[Nonresponder_cells,2]))
xlim1 <- range(c(ds1$Responders$x, ds1$NonResponders$x))
ylim1 <- range(c(ds1$Responders$y, ds1$NonResponders$y))


ds2 <- list(Responders = density(pseudotime[Responder_cells,3]),
           NonResponders = density(pseudotime[Nonresponder_cells,3]))
xlim2 <- range(c(ds2$Responders$x, ds2$NonResponders$x))
ylim2 <- range(c(ds2$Responders$y, ds2$NonResponders$y))



par(mfrow=c(1,3))
plot(xlim, ylim, col = "white", xlab = "Pseudotime",ylab="Density",xlim=c(-2,14),ylim=c(0,0.4),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main="Lineage 1")
polygon(c(min(ds$Responders$x),ds$Responders$x,max(ds$Responders$x)),
		c(0,ds$Responders$y,0),border="red")
polygon(c(min(ds$NonResponders$x),ds$Responders$x,max(ds$Responders$x)),
		c(0,ds$NonResponders$y,0), border = "black")
legend("topleft", legend=c("BCG responsive", "BCG unresponsive"),
       col=c("red", "black"), lty=1:1,bty="n",cex=1.5)


plot(xlim1, ylim1, col = "white", xlab = "Pseudotime", ylab = "",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main="Lineage 2")
polygon(c(min(ds1$Responders$x),ds1$Responders$x,max(ds1$Responders$x)),
		c(0,ds1$Responders$y,0), border="red")
polygon(c(min(ds1$NonResponders$x),ds1$Responders$x,max(ds1$Responders$x)),
		c(0,ds1$NonResponders$y,0), border="black")
legend("topleft", legend=c("BCG responsive", "BCG unresponsive"),
       col=c("red", "black"), lty=1:1,bty="n")




plot(xlim2, ylim2, col = "white", xlab = "Pseudotime", ylab = "",xlim=c(-2,14),ylim=c(0,0.4),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main="Lineage 3")
polygon(c(min(ds2$Responders$x),ds2$Responders$x,max(ds2$Responders$x)),
		c(0,ds2$Responders$y,0), border="red")
polygon(c(min(ds2$NonResponders$x),ds2$Responders$x,max(ds2$Responders$x)),
		c(0,ds2$NonResponders$y,0), border="black")

###Let us look at distribution differences between pseudotimes between conditions using different tests

###change into counts, pseudotime and cellWeights
counts=as.matrix(BD_Tcells_CD4.sct_1_subset@assays$RNA@counts[,c(Responder_cells,Nonresponder_cells)])
pseudotime_condition<-pseudotime[c(Responder_cells,Nonresponder_cells),]
cell_weight<-cellWeights[c(Responder_cells,Nonresponder_cells),]
conditions=as.factor(c(rep("BCG responsive",length(Responder_cells)), rep("BCG unresponsive",length(Nonresponder_cells))))


prog_res <- progressionTest(pseudotime_condition, conditions = conditions,
                   cellWeights = cell_weight, global = TRUE, lineages = TRUE,method="permutation")


##Look at density plot for lineage 2
condition=c(rep("BCG responsive",length(Responder_cells)), rep("BCG unresponsive",length(Nonresponder_cells)))
density_lineage2<-cbind.data.frame(pseudotime_tradeseq[,1],condition)
colnames(density_lineage2)<-c("Lineage2","conditions")

ggplot(density_lineage2, aes(x = Lineage2, fill = conditions)) +
  geom_density(alpha = .5) +ylab("Density")+xlab("Pseudotime")+
  scale_fill_brewer(type = "qual")+theme(plot.title = element_text(hjust = 0.5,face="bold",size = 20))+ theme(axis.text.x = element_text(size=15,color="black"),axis.text.y = element_text(size=16,color="black")) +theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))+ggtitle("Lineage 2")+theme(strip.background = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#a8ddb5","#43a2ca"))



###Repeat same analyses for CD8 T cells


###Now use the CD8 T cells for same clustering analyses

BD_Tcells_CD8<-subset(BD_Tcells,idents=c("T_CD8_memory","T_CD8_naive"))
BD_Tcells_CD8.list <- SplitObject(BD_Tcells_CD8, split.by="Patient_sampletag")


var.features <- SelectIntegrationFeatures(object.list = BD_Tcells_CD8.list , nfeatures = 3000)

BD_Tcells_CD8.sct <- merge(x = BD_Tcells_CD8.list[[1]], y =BD_Tcells_CD8.list[2:length(BD_Tcells_CD8.list)], merge.data=TRUE)

VariableFeatures(BD_Tcells_CD8.sct) <- var.features
BD_Tcells_CD8.sct<- RunPCA(BD_Tcells_CD8.sct , verbose = FALSE)
BD_Tcells_CD8.sct <- RunHarmony(BD_Tcells_CD8.sct , assay.use="SCT", group.by.vars = c("Patient","type"))
BD_Tcells_CD8.sct <- RunUMAP(BD_Tcells_CD8.sct , reduction = "harmony", dims = 1:30)
BD_Tcells_CD8.sct<- FindNeighbors(BD_Tcells_CD8.sct , reduction = "harmony", dims = 1:30) 
BD_Tcells_CD8.sct <- FindClusters(BD_Tcells_CD8.sct ,resolution = seq(from = 0.1, to = 1, by = 0.1))
clustree(BD_Tcells_CD8.sct, node_colour = "sc3_stability")


save(BD_Tcells_CD8.sct,file="BD_Tcells_CD8.sct.RData")
BD_Tcells_CD8.sct_1<- FindClusters(BD_Tcells_CD8.sct, resolution =  0.4)

##UMAP plots

DimPlot(BD_Tcells_CD8.sct_1, reduction = "umap",group.by="seurat_clusters") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("CD8")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)


DimPlot(BD_Tcells_CD8.sct_1, reduction = "umap",group.by="type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("T cell clusters")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)

###Let me check the clusters and what they would represent

Naive<-c("TCF7","SELL","LEF1","CCR7")
Effector_memory<-c("CCL4","CD69","IFNG")
Exhausted<-c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4","LAYN","ENTPD1")
Cytotoxic<-c("GZMA","GZMB","GZMK","GNLY","PRF1","NKG7")
Transcription_factor<-c("STAT4")

row_split_T =c(rep("Naive",4),rep("Effector\nmemory",3),rep("Exhausted",7),rep("Cytotoxic",6),rep("Transcription\nfactor",1))
row_split_T = factor(row_split_T,levels = c("Naive","Effector\nmemory","Exhausted","Cytotoxic","Transcription\nfactor"))

### plot heatmap
features_Tcell<-c(Naive,Effector_memory,Exhausted,Cytotoxic,Transcription_factor)


### set colnames order

plot_ord_T <-c("CD8_Naive","CD8_EM","CD8_Tex","CD8_Pex_1","CD8_Pex_2","CD8_STAT4","CD8_TEMRA")
data.plot_T <- Heat_Dot_data(object=BD_Tcells_CD8.sct_1_sub,features=features_Tcell,group.by="Tcell_type")

###Look at the heatmap
exp.mat_T <- data.plot_T %>% dplyr::select(features.plot,id,avg.exp.scaled) %>% spread(id,avg.exp.scaled)
rownames(exp.mat_T) <- exp.mat_T$features.plot
exp.mat_T$features.plot <- NULL
exp.mat_T <- exp.mat_T[,plot_ord_T]
per.mat_T <- data.plot_T %>% dplyr::select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat_T) <- per.mat_T$features.plot
per.mat_T$features.plot <- NULL
per.mat_T <- per.mat_T[,plot_ord_T]/100
min(exp.mat_T);max(exp.mat_T)


pdf("Tcell_cluster_CD8_type.pdf",height=14,width=10)
Heatmap(as.matrix(exp.mat_T),cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,
        column_names_side = "bottom",row_names_side = "left",
        row_split = row_split_T,
        row_gap = unit(3.5, "mm"),
        column_gap = unit(3.5, "mm"), width = ncol(exp.mat_T)*unit(7, "mm"),
    height = nrow(exp.mat_T)*unit(7, "mm"),
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()


####Cluster labelling

clusters = BD_Tcells_CD8.sct_1$seurat_clusters
Tcell_type<-case_when( clusters %in% 0 ~ 'CD8_EM',
               clusters %in% 1 ~ 'CD8_Pex_1',
               clusters %in% 2 ~ 'CD8_STAT4',
               clusters %in% 3 ~ 'CD8_Naive',
               clusters %in%4 ~ 'CD8_Pex_2',
               clusters %in% 5 ~ 'CD8_TEMRA' ,
                      clusters %in% 6 ~ 'CD8_Tex',
                       clusters %in% 7 ~ 'CD8_HLA-DR+' ,
                       clusters %in% 8 ~ 'CD8_Naive')

BD_Tcells_CD8.sct_1$Tcell_type<-Tcell_type

##Remove cluster 7 small and labelling not sure
Idents(BD_Tcells_CD8.sct_1)<-"Tcell_type"
BD_Tcells_CD8.sct_1_sub<-subset(BD_Tcells_CD8.sct_1,idents=c("CD8_HLA-DR+"),invert=T)
save(BD_Tcells_CD8.sct_1_sub,file="BD_Tcells_CD8.sct_1_sub.RData") ###Final CD8 clusters


###Visualise tissue and blood CD8 T cells

BD_Tcells_CD8.sct_1_sub$type <- factor(x = BD_Tcells_CD8.sct_1_sub$type, levels = c("tissue","blood"))
t8<-ggplot(BD_Tcells_CD8.sct_1_sub@meta.data, aes(x=type, fill=Tcell_type)) + geom_bar(position="fill")+xlab("Type")+ylab("Fraction")+theme(text = element_text(size = 14,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=14), # bold
    axis.title.x = element_text(size=14, face="bold", colour = "black"),
     axis.text.x = element_text(face="bold", colour = "black",size=10), axis.ticks.x=element_blank())+theme(plot.title=element_text(hjust=0.5))+ theme(strip.background = element_blank(),
   strip.text.x = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

##Final UMAPs
t9<-DimPlot(BD_Tcells_CD8.sct_1_sub, reduction = "umap",group.by="Tcell_type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("CD8")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)

LabelClusters(plot = t9, id = 'Tcell_type',fontface = "bold")

t10<-DimPlot(BD_Tcells_CD8.sct_1_sub, reduction = "umap",group.by="type") +theme(text = element_text(size = 20,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=20),axis.title.x = element_text(size=20, face="bold", colour = "black"),
axis.text.x = element_text(face="bold", colour = "black",size=20))+
ggtitle("T cell clusters")+
guides(x = axis, y = axis) +
theme(axis.line = element_line(arrow = arrow()),
 axis.title = element_text(hjust = 0))+scale_x_continuous(breaks = NULL) +
scale_y_continuous(breaks = NULL)


####DE analyses MAST for tissue

Idents(BD_Tcells_CD8.sct_1_sub)<-"type"
Tissue_CD8=subset(BD_Tcells_CD8.sct_1_sub,idents=c("tissue"))


###CD8 T cells  BCG unresponsive vs BCG responsive


Idents(Tissue_CD8)<-"Response"
DE_CD8 <- FindMarkers(Tissue_CD8, assay ="SCT",recorrect_umi=FALSE,ident.1="BCG_responsive",ident.2="BCG_unresponsive",test.use="MAST",latent.vars="Sex")

DE_CD8 <- DE_CD8 %>%
  mutate(Genes = case_when(avg_log2FC>= 0.58 & p_val_adj <= 0.05 ~ "Up in BCG_responsive",
                             avg_log2FC <= -0.58 &  p_val_adj <= 0.05 ~ "Up in BCG_unresponsive",
                               TRUE ~ "NS"))


sig_genes_CD8<-DE_CD8%>% filter(Genes!="NS")


sig_genes_CD8_up<-sig_genes_CD8%>% dplyr::filter(Genes=="Up in BCG_responsive")%>% arrange(desc(avg_log2FC))

sig_genes_CD8_down<-sig_genes_CD8%>% dplyr::filter(Genes=="Up in BCG_unresponsive")%>% arrange(desc(avg_log2FC))


###Enrichment analyses

go_CD8<-enrichGO(
rownames(sig_genes_CD8_up),
'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")


###what are the pathways
barplot(go_CD8, showCategory=10) 

 ###Enrichment analyses

go_CD8_down<-enrichGO(
rownames(sig_genes_CD8_down),
'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")


###what are the pathways
barplot(go_CD8_down, showCategory=10) 


###differential expression tissue CD8 T cells bcg exposed/naive 

Idents(Tissue_CD8)<-"exposure"
DE_CD8_e <- FindMarkers(Tissue_CD8, assay = "SCT",recorrect_umi=FALSE,ident.1="BCG_exposed",ident.2="BCG_naive",test.use="MAST",latent.vars="Sex")

DE_CD8_e <- DE_CD8_e %>%
  mutate(Genes = case_when(avg_log2FC>= 0.58 & p_val_adj <= 0.05 ~ "Up in BCG_exposed",
                             avg_log2FC <= -0.58 &  p_val_adj <= 0.05 ~ "Up in BCG_naive",
                               TRUE ~ "NS"))


sig_genes_CD8_e<-DE_CD8_e%>% filter(Genes!="NS")


sig_genes_CD8_e_up<-sig_genes_CD8_e%>% dplyr::filter(Genes=="Up in BCG_exposed")%>% arrange(desc(avg_log2FC))

sig_genes_CD8_e_down<-sig_genes_CD8_e%>% dplyr::filter(Genes=="Up in BCG_naive")%>% arrange(desc(avg_log2FC))

###Enrichment analyses
go_CD8_e<-enrichGO(
rownames(sig_genes_CD8_e_up),
'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")


barplot(go_CD8_e, showCategory=10) 



go_CD8_e_down<-enrichGO(
rownames(sig_genes_CD8_e_down),
'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")


barplot(go_CD8_e_down, showCategory=10) 

##Visualise the cluster proportions

##Make stacked barplots
Idents(Tissue_CD8)<-"Response"
Tissue_CD8_response<-subset(Tissue_CD8, idents=c("BCG_naive"),invert=T)

n2<-ggplot(Tissue_CD8_response@meta.data, aes(x=Response, fill=Tcell_type)) + geom_bar(position="fill")+xlab("")+ylab("Fraction")+theme(text = element_text(size = 18,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=18), # bold
    axis.title.x = element_text(size=18, face="bold", colour = "black"),
     axis.text.x = element_text(face="bold", colour = "black",size=18), axis.ticks.x=element_blank())+theme(plot.title=element_text(hjust=0.5))+ theme(strip.background = element_blank(),
   strip.text.x = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(fill="CD8")

Tissue_CD8$exposure<- factor(x = Tissue_CD8$exposure, levels = c("BCG_naive","BCG_exposed"))

n3<-ggplot(Tissue_CD8@meta.data, aes(x=exposure, fill=Tcell_type)) + geom_bar(position="fill")+xlab("")+ylab("Fraction")+theme(text = element_text(size = 18,face="bold"), axis.text.y = element_text(face="bold", colour = "black",size=18), # bold
    axis.title.x = element_text(size=18, face="bold", colour = "black"),
     axis.text.x = element_text(face="bold", colour = "black",size=18), axis.ticks.x=element_blank())+theme(plot.title=element_text(hjust=0.5))+ theme(strip.background = element_blank(),
   strip.text.x = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(fill="CD8")


##comparing proportions using propeller

##Comparing cluster proportions using CD8 T cells between tissue and blood
Tissue_blood_comparison_CD8<-propeller(
  x =BD_Tcells_CD8.sct_1_sub,
  clusters = BD_Tcells_CD8.sct_1_sub@meta.data$Tcell_type,
  sample = BD_Tcells_CD8.sct_1_sub@meta.data$Patient_sampletype,
  group = BD_Tcells_CD8.sct_1_sub@meta.data$type,
  transform = "asin"
)




####Cluster proportion analyses Speckle #tissue
#Between BCG responders and non-responders


Idents(Tissue_CD8)<-"Patient"
Tissue_CD8_sub<-subset(Tissue_CD8,idents=c("P2"),invert=T)

##subset further

Idents(Tissue_CD8_sub)<-"Response"
Tissue_CD8_R<-subset(Tissue_CD8_sub,idents=c("BCG_responsive","BCG_unresponsive"))



output_tcell.asin <- propeller(clusters=Tissue_CD8_R@meta.data$Tcell_type, sample=Tissue_CD8_R@meta.data$Patient, group=Tissue_CD8_R@meta.data$Response, transform="asin") 

props <- getTransformedProps(Tissue_CD8_R@meta.data$Tcell_type,Tissue_CD8_R@meta.data$Patient, transform="asin")


##Between BCG exposed and naive
output_t.asin <- propeller(clusters=Tissue_CD8_sub@meta.data$Tcell_type, sample=Tissue_CD8_sub@meta.data$Patient, group=Tissue_CD8_sub@meta.data$exposure, transform="asin") 



###Do similar analyses for CD8 T cell blood. 

Idents(BD_Tcells_CD8.sct_1_sub)<-"type"
Blood_CD8=subset(BD_Tcells_CD8.sct_1_sub,idents=c("blood"))

Blood_CD8@meta.data$exposure<-ifelse(Blood_CD8@meta.data$Response=="BCG_naive","BCG_naive","BCG_exposed")

###Proportion analyses

Idents(Blood_CD8)<-"Patient"
Blood_CD8_sub<-subset(Tissue_CD8,idents=c("P2","P14","P16"),invert=T)

##subset further

Idents(Blood_CD8_sub)<-"Response"
Blood_CD8_R<-subset(Blood_CD8_sub,idents=c("BCG_responsive","BCG_unresponsive"))

output_tcell.asin <- propeller(clusters=Blood_CD8_R@meta.data$Tcell_type, sample=Blood_CD8_R@meta.data$Patient, group=Blood_CD8_R@meta.data$Response, transform="asin") 
grp.BCG <- rep(c("BCG_responsive","BCG_unresponsive"), c(7, 5))


output_tcell.asin <- propeller(clusters=Blood_CD8_sub@meta.data$Tcell_type, sample=Blood_CD8_sub@meta.data$Patient, group=Blood_CD8_sub@meta.data$exposure, transform="asin") 

props <- getTransformedProps(Blood_CD8_sub@meta.data$Tcell_type,Blood_CD8_sub@meta.data$Patient, transform="asin")

##Blood CD8 differential expression

Idents(Blood_CD8)<-"Response"
DE_CD8_blood <- FindMarkers(Blood_CD8, assay = "SCT",recorrect_umi=FALSE,ident.1="BCG_responsive",ident.2="BCG_unresponsive",test.use="MAST",latent.vars="Sex")

DE_CD8_blood <- DE_CD8_blood %>%
  mutate(Genes = case_when(avg_log2FC>= 0.58 & p_val_adj <= 0.05 ~ "Up in BCG responsive",
                             avg_log2FC <= -0.58 &  p_val_adj <= 0.05 ~ "Up in BCG unresponsive",
                               TRUE ~ "NS"))


sig_genes_CD8_blood<-DE_CD8_blood%>% filter(Genes!="NS")


sig_genes_CD8_up_blood<-sig_genes_CD8_blood%>% dplyr::filter(Genes=="Up in BCG responsive")%>% arrange(desc(avg_log2FC))

sig_genes_CD8_down_blood<-sig_genes_CD8_blood%>% dplyr::filter(Genes=="Up in BCG unresponsive")%>% arrange(desc(avg_log2FC))


Idents(Blood_CD8)<-"exposure"
DE_CD8_blood_e <- FindMarkers(Blood_CD8, assay = "SCT",recorrect_umi=FALSE,ident.1="BCG_exposed",ident.2="BCG_naive")

DE_CD8_blood_e <- DE_CD8_blood_e %>%
  mutate(Genes = case_when(avg_log2FC>= 0.58 & p_val_adj <= 0.05 ~ "Up in BCG exposed",
                             avg_log2FC <= -0.58 &  p_val_adj <= 0.05 ~ "Up in BCG naive",
                               TRUE ~ "NS"))


sig_genes_CD8_blood_e<-DE_CD8_blood_e%>% filter(Genes!="NS")


sig_genes_CD8_up_blood_e<-sig_genes_CD8_blood_e%>% dplyr::filter(Genes=="Up in BCG exposed")%>% arrange(desc(avg_log2FC))

sig_genes_CD8_down_blood_e<-sig_genes_CD8_blood_e%>% dplyr::filter(Genes=="Up in BCG naive")%>% arrange(desc(avg_log2FC))



###Look at trajectory analyses

set.seed(1)
lineages <- as.SlingshotDataSet(getLineages(
  data           = BD_Tcells_CD8.sct_1_sub@reductions$umap@cell.embeddings,
  clusterLabels  = BD_Tcells_CD8.sct_1_sub$Tcell_type,
  dist.method    = "mnn", 
  start.clus     = "CD8_Naive"

)) # define where to START the trajectories


cell_colors_clust <- cell_pal(BD_Tcells_CD8.sct_1_sub$Tcell_type, hue_pal())

plot(BD_Tcells_CD8.sct_1_sub@reductions$umap@cell.embeddings, cex = .5, pch = 16,col=cell_colors_clust)

lines(lineages, lwd = 1, col = "black", cex = 2)


curves <- as.SlingshotDataSet(getCurves(
  data          = lineages,
  thresh        = 1e-1,
  stretch       = 1e-1,
  allow.breaks  = F,
  approx_points = 100
))

plot(BD_Tcells_CD8.sct_1_sub@reductions$umap@cell.embeddings, col = cell_colors_clust, cex = .5, pch = 16)
lines(curves, lwd = 2, col = "black")


pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)

x <- rowMeans(pseudotime)
x <- x / max(x)
o <- order(x)

plot(BD_Tcells_CD8.sct_1_sub@reductions$umap@cell.embeddings[o, ],
    main = paste0("pseudotime"), pch = 16, cex = 0.4, axes = F, xlab = "", ylab = "",
    col = colorRampPalette(c("grey70", "orange3", "firebrick", "purple4"))(99)[x[o] * 98 + 1]
  )
lines(curves, lwd = 2, col = "black")


###Look at density plots within each trajectory to compare conditions
Idents(BD_Tcells_CD8.sct_1_sub)<-"Response"
Responder_cells_CD8=colnames(subset(BD_Tcells_CD8.sct_1_sub,idents=c("BCG_responsive")))
Nonresponder_cells_CD8=colnames(subset(BD_Tcells_CD8.sct_1_sub,idents=c("BCG_unresponsive")))

ds <- list(Responders_CD8 = density(pseudotime[Responder_cells_CD8,1]),
           NonResponders = density(pseudotime[Nonresponder_cells_CD8,1]))
xlim <- range(c(ds$Responders$x, ds$NonResponders$x))
ylim <- range(c(ds$Responders$y, ds$NonResponders$y))

ds1 <- list(Responders = density(pseudotime[Responder_cells_CD8,2]),
           NonResponders = density(pseudotime[Nonresponder_cells_CD8,2]))
xlim1 <- range(c(ds1$Responders$x, ds1$NonResponders$x))
ylim1 <- range(c(ds1$Responders$y, ds1$NonResponders$y))


ds2 <- list(Responders = density(pseudotime[Responder_cells_CD8,3]),
           NonResponders = density(pseudotime[Nonresponder_cells_CD8,3]))
xlim2 <- range(c(ds2$Responders$x, ds2$NonResponders$x))
ylim2 <- range(c(ds2$Responders$y, ds2$NonResponders$y))



par(mfrow=c(1,3))
plot(xlim, ylim, col = "white", xlab = "Pseudotime",ylab="Density",xlim=c(-2,14),ylim=c(0,0.4),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main="Lineage 1")
polygon(c(min(ds$Responders$x),ds$Responders$x,max(ds$Responders$x)),
		c(0,ds$Responders$y,0),border="red")
polygon(c(min(ds$NonResponders$x),ds$Responders$x,max(ds$Responders$x)),
		c(0,ds$NonResponders$y,0), border = "black")
legend("topleft", legend=c("BCG responsive", "BCG unresponsive"),
       col=c("red", "black"), lty=1:1,bty="n",cex=1.5)


plot(xlim1, ylim1, col = "white", xlab = "Pseudotime", ylab = "",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main="Lineage 2")
polygon(c(min(ds1$Responders$x),ds1$Responders$x,max(ds1$Responders$x)),
		c(0,ds1$Responders$y,0), border="red")
polygon(c(min(ds1$NonResponders$x),ds1$Responders$x,max(ds1$Responders$x)),
		c(0,ds1$NonResponders$y,0), border="black")
legend("topleft", legend=c("BCG responsive", "BCG unresponsive"),
       col=c("red", "black"), lty=1:1,bty="n")




plot(xlim2, ylim2, col = "white", xlab = "Pseudotime", ylab = "",xlim=c(-2,14),ylim=c(0,0.4),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main="Lineage 3")
polygon(c(min(ds2$Responders$x),ds2$Responders$x,max(ds2$Responders$x)),
		c(0,ds2$Responders$y,0), border="red")
polygon(c(min(ds2$NonResponders$x),ds2$Responders$x,max(ds2$Responders$x)),
		c(0,ds2$NonResponders$y,0), border="black")


###Within each trajectory can we see what genes vary within each treatment

###change into counts, pseudotime and cellWeights
counts_CD8=as.matrix(BD_Tcells_CD8.sct_1_sub@assays$RNA@counts[,c(Responder_cells_CD8,Nonresponder_cells_CD8)])
pseudotime_CD8<-pseudotime[c(Responder_cells_CD8,Nonresponder_cells_CD8),]
cell_weight_CD8<-cellWeights[c(Responder_cells_CD8,Nonresponder_cells_CD8),]
conditions_CD8=as.factor(c(rep("BCG responsive",length(Responder_cells_CD8)), rep("BCG unresponsive",length(Nonresponder_cells_CD8))))


prog_res <- progressionTest(pseudotime_CD8, conditions = conditions_CD8,method="Permutation",
                   cellWeights = cell_weight_CD8, global = TRUE, lineages = TRUE)


##Look at density plot for lineage 1

density_lineage1<-cbind.data.frame(pseudotime_CD8[,1],conditions_CD8)
colnames(density_lineage1)<-c("Lineage1","conditions")

ggplot(density_lineage1, aes(x = Lineage1, fill = conditions_CD8)) +
  geom_density(alpha = .5) +ylab("Density")+xlab("Pseudotime")+
  scale_fill_brewer(type = "qual")+theme(plot.title = element_text(hjust = 0.5,face="bold",size = 20))+ theme(axis.text.x = element_text(size=15,color="black"),axis.text.y = element_text(size=16,color="black")) +theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))+ggtitle("Lineage 1")+theme(strip.background = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#a8ddb5","#43a2ca")) + scale_x_continuous(breaks=seq(0, 12.5, 2.5))



##Look at density plot for lineage 2

density_lineage2<-cbind.data.frame(pseudotime_CD8[,2],conditions_CD8)
colnames(density_lineage2)<-c("Lineage2","conditions")

ggplot(density_lineage2, aes(x = Lineage2, fill = conditions_CD8)) +
  geom_density(alpha = .5) +ylab("Density")+xlab("Pseudotime")+
  scale_fill_brewer(type = "qual")+theme(plot.title = element_text(hjust = 0.5,face="bold",size = 20))+ theme(axis.text.x = element_text(size=15,color="black"),axis.text.y = element_text(size=16,color="black")) +theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16))+ggtitle("Lineage 2")+theme(strip.background = element_blank())+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual(values=c("#a8ddb5","#43a2ca"))



