#############
# Begin
#############
#load Seurat
library(Seurat)
library(patchwork)
library(future)

#load generic libraries
library(tidyverse)
library(viridis)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(data.table)


#############
#Preferences
#############
seurat_umap_alpha_hex_string = "d0"
#setwd("/home/martin/Documents/20-22/Derek scRNAseq 2022 analysis/Liu Nat Com 2022/Process")
paired_colors <- brewer.pal(name = "Paired", n=12)
beach_colors <- c("#1a748e","#f0df99","#55c4d7","#5f6c24","#992915","#d38e31","#81c7f8","#bbc4af","#393430","#92dccd","#63fba9","#1155d4","#abb2ba","#eccd16","#7b5c52","#063581","#ecdfcf")
beach_colors_alpha <- paste0( beach_colors, seurat_umap_alpha_hex_string )
paired_colors_alpha <- paste0( paired_colors, seurat_umap_alpha_hex_string )



##############################################################################
#
#     RAW DATA INPUT -- download raw data from Liu et al. Nature Communications 2022
#
##############################################################################

###################################
#Download and read in data, including metadata; subset before even creating Seurat object
###################################
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE160741&format=file",method="curl",destfile = "GSM4878385_E18_5.tar.gz")
untar("GSM4878385_E18_5.tar.gz",files="GSM4878385_E18_5_WT_gene_counts.tsv.gz")
file.remove("GSM4878385_E18_5.tar.gz")
complete_bm_matrix <- fread("GSM4878385_E18_5_WT_gene_counts.tsv.gz",sep="\t",header=TRUE,data.table = FALSE)
complete_bm_matrix$gene <- as.factor(complete_bm_matrix$gene)
complete_bm_matrix$cell <- as.factor(complete_bm_matrix$cell)
complete_bm_matrix <- as.matrix(spread(complete_bm_matrix, key="cell", value="count", fill = 0, convert = FALSE))
rownames(complete_bm_matrix) <- complete_bm_matrix[,1]
complete_bm_matrix <- complete_bm_matrix[,-1, drop = FALSE]
#write.table(complete_bm_matrix, file = "/home/martin/Downloads/GSM4878385_matrix.csv", sep="\t", row.names = FALSE)

#############
#Create initial Seurat object
#############
Liu.seurat <- CreateSeuratObject(counts = complete_bm_matrix, project = "Embryonic Bone Marrow")
Liu.seurat[["percent.mt"]] <- PercentageFeatureSet(Liu.seurat, pattern = "^mt-")
VlnPlot(Liu.seurat,features=c("nCount_RNA","nFeature_RNA","percent.mt"))
Liu.seurat <- subset(Liu.seurat, subset = percent.mt <15 & nFeature_RNA > 500 )
save(Liu.seurat,file="Liu.seurat.Rds")


##############################################################################
#
#     END RAW DATA INPUT
#
##############################################################################




#############
# Run the standard workflow for visualization and clustering
#############
load("Liu.seurat.Rds")

Liu.seurat <- NormalizeData(Liu.seurat)
Liu.seurat <- ScaleData(Liu.seurat, verbose = FALSE)
Liu.seurat <- FindVariableFeatures(Liu.seurat, selection.method = "vst", nfeatures = 2000)
Liu.seurat <- RunPCA(Liu.seurat, npcs = 30, verbose = FALSE)
ElbowPlot( Liu.seurat, ndims=30)
Liu.seurat <- RunUMAP(Liu.seurat, reduction = "pca", dims = 1:15)
Liu.seurat <- FindNeighbors(Liu.seurat, reduction = "pca", dims = 1:15)
Liu.seurat <- FindClusters(Liu.seurat, resolution = 0.25)

## Visualization
p1 <- FeaturePlot(Liu.seurat, features=c("Cdh5"))
p2 <- DimPlot(Liu.seurat,  reduction = "umap", cols=beach_colors_alpha, pt.size=0.35, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p2 + p1

#EC cluster is #3
Liu.seurat <- subset( Liu.seurat, idents = 3 )
save(Liu.seurat,file="Liu_endothelial.seurat.Rds")
#ggsave( plot = p1 + p2, height=150, width=300, filename = "He_complete_UMAP.png", units = "mm",dpi = 300 )


#############
# Re-cluster with the standard workflow for visualization and clustering
#############
Liu.seurat <- FindVariableFeatures(Liu.seurat, selection.method = "vst", nfeatures = 2000)
Liu.seurat <- RunPCA(Liu.seurat, npcs = 30, verbose = FALSE)
ElbowPlot( Liu.seurat, ndims=30)
Liu.seurat <- RunUMAP(Liu.seurat, reduction = "pca", dims = 1:15)
Liu.seurat <- FindNeighbors(Liu.seurat, reduction = "pca", dims = 1:15)
Liu.seurat <- FindClusters(Liu.seurat, resolution = 0.25)

## Gross Visualization
p1 <- FeaturePlot(Liu.seurat, features=c("Flt4"))
p3 <- FeaturePlot(Liu.seurat, features=c("Kdr"))
p2 <- DimPlot(Liu.seurat,  reduction = "umap", cols=beach_colors_alpha, pt.size=0.35, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p2 + p1 + p3

#all.markers <- FindAllMarkers(object = Liu.seurat)
#top20 <- all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC ) # * pct.1* pct.1)
#DoHeatmap(object = Liu.seurat, features = top20$gene, label=FALSE) #slim.col.label = TRUE, remove.key = TRUE)
#VlnPlot(Liu.seurat,features=c("Cxcl12","Jag1","Jag2","Dll4","Vegfc","Kitl","Pdgfb","Tgfb2"))

#############
#Cluster naming
#############
Idents(Liu.seurat) <- "seurat_clusters"
Liu.seurat[["old.ident"]] <- Idents(Liu.seurat)
Liu.seurat <- RenameIdents(object = Liu.seurat, 
                            '1' = "SECs", 
                            '0' = "Capillary ECs",
                           '2' = "Arterial BECs"
)
Liu.seurat[["cluster_names"]] <- Idents(Liu.seurat)


#############
#Visualize final object
#############
#p1 <- DimPlot(Liu.seurat, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.65) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(Liu.seurat,  reduction = "umap", cols=beach_colors_alpha, group.by = "cluster_names", pt.size=0.65, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p2
#ggsave( plot = p1 + p2, height=150, width=300, filename = "He_complete_UMAP.pdf", units = "mm",dpi = 300 )


#HM <- DoHeatmap(object = Liu.seurat, features = top20$gene, label=FALSE, group.colors = beach_colors)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(axis.text.y = element_text(size = 15)) + NoLegend()#slim.col.label = TRUE, remove.key = TRUE)
#ggsave( plot = HM, height=380, width=300, filename = "He_complete_clusterHM.pdf", units = "mm",dpi = 300 )


#ggsave( plot = ((p1+p2)/HM) + plot_layout(heights=c(1,3)), height=300, width=300, filename = "He_supplement.pdf", units = "mm",dpi = 300 )
DefaultAssay(Liu.seurat) <- "RNA"
my.markers <- c("Pecam1","Emcn","Stab1","Stab2","Efnb2","Acvrl1","Nrp1","Sox17","Hey1")
DP <- DotPlot(Liu.seurat, features = my.markers, dot.scale = 8 ) + #, cols = c("Blue","Orange")) +
  scale_y_discrete(limits=rev) + RotatedAxis() + theme(axis.title.y = element_blank(),axis.title.x = element_blank())

ggsave( plot = (p2/DP), height=150, width=150, filename = "Liu_supplement_DP.pdf", units = "mm",dpi = 300 )




#############
#Cluster compaction
#############
Liu.compact_seurat <- subset( Liu.seurat, idents = c("Arterial BECs","SECs") )
#Liu.compact_seurat = Liu.seurat
Liu.compact_seurat[["cluster_names"]] <- Idents(Liu.compact_seurat)
Liu.compact_seurat <- RenameIdents(object = Liu.compact_seurat, 
                                   #"Capillary ECs" = "BECs", 
                                   "Arterial BECs" = "BECs", 
                                    "SECs" = "SECs"
)
Liu.compact_seurat[["cluster_names"]] <- Idents(Liu.compact_seurat)

#FeaturePlot(Wang.seurat,features=c("Pecam1","Cdh5","Bcam","Vwf","Sox17","Stab2","Stab1","Flt4"))
#VlnPlot(He_endothelial.compact_seurat, group.by="cluster_names", features=c("Vwf","Lyve1","Kdr","Flt4","Cdh5","Nr2f2","Stab2","Stab1","Bcam","Pecam1","Cd34","Sox17","Emcn"))




#############
#Violin for BEC/SEC/LEC
#############
DefaultAssay(Liu.compact_seurat) <- "RNA"
vln <- VlnPlot(Liu.compact_seurat, cols=c(beach_colors[3],beach_colors[1]), features = c("Cdh5","Kdr","Flt4"), ncol = 3, log=FALSE, combine = FALSE, pt.size=0)
for(i in 1:length(vln)) {
  if ( i == length(vln) ) {
    vln[[i]] <- vln[[i]] + NoAxes() + theme(plot.margin = margin(0, 0, 5, 4, "mm"), plot.title = element_text(size=12,face="plain",vjust = -65 ) ) + stat_compare_means(comparisons=list(levels(Liu.compact_seurat$cluster_names)), method="t.test", label = "p.signif", label.x.npc=0.1, paired=FALSE, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("P<0.0001", "P<0.001", "P<0.01", "P<0.05", "n.s.")), size=2.5 )+ scale_y_continuous(limits=c(-0.1,5))   #comparisons=list(levels(Wang.compact.seurat$cluster_names)), 
  } else {
    vln[[i]] <- vln[[i]] + NoAxes() +NoLegend() + theme(plot.margin = margin(0, 8, 5, 2, "mm"), plot.title = element_text(size=12,face="plain",vjust = -65 ) ) + stat_compare_means(comparisons=list(levels(Liu.compact_seurat$cluster_names)),method="t.test", label = "p.signif", label.x.npc=0.1, paired=FALSE, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("P<0.0001", "P<0.001", "P<0.01", "P<0.05", "n.s.")), size=2.5 )+ scale_y_continuous(limits=c(-0.1,5))   #comparisons=list(levels(Wang.compact.seurat$cluster_names)), 
  }
  vln[[i]]$layers[[1]]$aes_params$size = 0
}
#wrap_plots(vln,nrow=1)
ggsave( plot = wrap_plots(vln,nrow=1) , height=70, width=120, filename = "Liu_Vlns.pdf", units = "mm",dpi = 300,useDingbats = FALSE)



