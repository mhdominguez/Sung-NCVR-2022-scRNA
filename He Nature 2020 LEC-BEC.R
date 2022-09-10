

#############
# Begin
#############
library(Seurat)
library(patchwork)

#load generic libraries
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(data.table)



#############
# Preferences
#############
seurat_umap_alpha_hex_string = "d0"
working_path = "/home/martin/Documents/20-22/Derek scRNAseq 2022 analysis/He Nature 2020/Process"
setwd(working_path)
#write(paste0("TMP = \'",working_path,"\'"), file=file.path('~/.Renviron'))
paired_colors <- brewer.pal(name = "Paired", n=12)
beach_colors <- c("#1a748e","#f0df99","#55c4d7","#5f6c24","#992915","#d38e31","#81c7f8","#bbc4af","#393430","#92dccd","#63fba9","#1155d4","#abb2ba","#eccd16","#7b5c52","#063581","#ecdfcf")
beach_colors_alpha <- paste0( beach_colors, seurat_umap_alpha_hex_string )
paired_colors_alpha <- paste0( paired_colors, seurat_umap_alpha_hex_string )



##############################################################################
#
#     RAW DATA INPUT -- download raw data from He et al. Nature 2020
#
##############################################################################

###################################
#Download and read in data, including metadata; subset before even creating Seurat object
###################################
download.file("https://cells.ucsc.edu/mouse-limb/10x/exprMatrix.tsv.gz",method="curl",destfile = "exprMatrix.tsv.gz")
complete_limb_matrix <- fread("exprMatrix.tsv.gz",sep="\t",header=TRUE,data.table = FALSE,integer64="numeric")
rownames(complete_limb_matrix) <- complete_limb_matrix$gene
complete_limb_matrix = complete_limb_matrix[,!(names(complete_limb_matrix) %in% "gene")]

file.remove("exprMatrix.tsv")
download.file("https://cells.ucsc.edu/mouse-limb/10x/meta.tsv",destfile="meta.tsv")
cell_key = read.csv(file = "meta.tsv", stringsAsFactors = TRUE, header=TRUE, sep="\t" , row.names = 1)
cell_key <- cell_key[ which(cell_key$cell_type=="Endothelial "), ]
cell_key$stage <- paste( "E", cell_key$stage, sep="" )
cells.use <- rownames(cell_key)
endothelial_limb_matrix <- complete_limb_matrix[ , names(complete_limb_matrix) %in% cells.use ]
remove(complete_limb_matrix)
gc()


#############
#Create initial Seurat object
#############
He_endothelial_seurat <- CreateSeuratObject(counts = endothelial_limb_matrix, project = "Embryonic Forelimb Endothelial")
He_endothelial_seurat[["percent.mt"]] <- PercentageFeatureSet(He_endothelial_seurat, pattern = "^mt-")
He_endothelial_seurat <- AddMetaData( object = He_endothelial_seurat, metadata = cell_key   )
VlnPlot(He_endothelial_seurat, features=c("nFeature_RNA","Sox17","nCount_RNA","Cdh5"))
save(He_endothelial_seurat,file="He_endothelial.seurat.Rds")


##############################################################################
#
#     END RAW DATA INPUT
#
##############################################################################




#############
# Run the standard workflow for visualization and clustering
#############
load("He_endothelial.seurat.Rds")

#this data is already normalized
#He_endothelial_seurat <- NormalizeData(He_endothelial_seurat)
He_endothelial_seurat <- subset(He_endothelial_seurat, subset = nCount_RNA > 1000 & nFeature_RNA > 2000 )
He_endothelial_seurat <- ScaleData(He_endothelial_seurat, verbose = FALSE)
He_endothelial_seurat <- FindVariableFeatures(He_endothelial_seurat, selection.method = "vst", nfeatures = 2000)
He_endothelial_seurat <- RunPCA(He_endothelial_seurat, npcs = 30, verbose = FALSE)
ElbowPlot( He_endothelial_seurat, ndims=30)
He_endothelial_seurat <- RunUMAP(He_endothelial_seurat, reduction = "pca", dims = 1:15)
He_endothelial_seurat <- FindNeighbors(He_endothelial_seurat, reduction = "pca", dims = 1:15)
He_endothelial_seurat <- FindClusters(He_endothelial_seurat, resolution = 0.25)

## Visualization
#p1 <- DimPlot(He_endothelial_seurat, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.35) + NoAxes() + labs(title = "Tissue") + theme(plot.title = element_text(hjust = 0.5))
#p2 <- DimPlot(He_endothelial_seurat,  reduction = "umap", cols=beach_colors_alpha, pt.size=0.35, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
#p1 + p2
#ggsave( plot = p1 + p2, height=150, width=300, filename = "He_complete_UMAP.png", units = "mm",dpi = 300 )


#############
#Cluster naming
#############
Idents(He_endothelial_seurat) <- "seurat_clusters"
He_endothelial_seurat[["old.ident"]] <- Idents(He_endothelial_seurat)
He_endothelial_seurat <- RenameIdents(object = He_endothelial_seurat, 
                            '0' = "BECs", 
                            '1' = "Endothelial Progenitors", 
                            '2' = "Mesenchymal Cells", 
                            '3' = "LECs"
)
He_endothelial_seurat[["cluster_names"]] <- Idents(He_endothelial_seurat)
He_endothelial_seurat <- subset( He_endothelial_seurat, idents = c("BECs","Mesenchymal Cells","LECs") ) 

#############
#Visualize final object
#############
p1 <- DimPlot(He_endothelial_seurat, reduction = "umap", cols=paired_colors_alpha, group.by = "stage", pt.size=0.65) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(He_endothelial_seurat,  reduction = "umap", cols=beach_colors_alpha, group.by = "cluster_names", pt.size=0.65, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p1 + p2
#ggsave( plot = p1 + p2, height=150, width=300, filename = "He_complete_UMAP.pdf", units = "mm",dpi = 300 )

#all.markers <- FindAllMarkers(object = He_endothelial_seurat)
#top20 <- all.markers %>% group_by(cluster) %>% top_n(8, avg_log2FC * pct.1* pct.1)
#DoHeatmap(object = Wang.seurat, features = top20$gene, label=FALSE) #slim.col.label = TRUE, remove.key = TRUE)

#HM <- DoHeatmap(object = He_endothelial_seurat, features = top20$gene, label=FALSE, group.colors = beach_colors)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(axis.text.y = element_text(size = 15)) + NoLegend()#slim.col.label = TRUE, remove.key = TRUE)
#ggsave( plot = HM, height=380, width=300, filename = "He_complete_clusterHM.pdf", units = "mm",dpi = 300 )


#ggsave( plot = ((p1+p2)/HM) + plot_layout(heights=c(1,3)), height=300, width=300, filename = "He_supplement.pdf", units = "mm",dpi = 300 )
DefaultAssay(He_endothelial_seurat) <- "RNA"
my.markers <- c("Pecam1","Emcn","Efnb2","Acvrl1","Nrp1","Flt1","Sox17","Hey1","Ephb4","Lyve1","Prox1","Pdpn","Pdgfrb","Gas1")
DP <- DotPlot(He_endothelial_seurat, features = my.markers, dot.scale = 8 ) + #, cols = c("Blue","Orange")) +
  scale_y_discrete(limits=rev) + RotatedAxis() + theme(axis.title.y = element_blank(),axis.title.x = element_blank())

ggsave( plot = ((p1+p2)/DP) + plot_layout(heights=c(2.5,1)), height=180, width=300, filename = "He_supplement_DP.pdf", units = "mm",dpi = 300 )




#############
#Cluster compaction
#############
He_endothelial.compact_seurat <- subset( He_endothelial_seurat, idents = c("BECs","LECs") )
He_endothelial.compact_seurat[["cluster_names"]] <- Idents(He_endothelial.compact_seurat)
He_endothelial.compact_seurat <- RenameIdents(object = He_endothelial.compact_seurat, 
                                    "BECs" = "BECs", 
                                    "LECs" = "LECs"
)
He_endothelial.compact_seurat[["cluster_names"]] <- Idents(He_endothelial.compact_seurat)

#FeaturePlot(Wang.seurat,features=c("Pecam1","Cdh5","Bcam","Vwf","Sox17","Stab2","Stab1","Flt4"))
#VlnPlot(He_endothelial.compact_seurat, group.by="cluster_names", features=c("Vwf","Lyve1","Kdr","Flt4","Cdh5","Nr2f2","Stab2","Stab1","Bcam","Pecam1","Cd34","Sox17","Emcn"))




#############
#Violin for BEC/SEC/LEC
#############
DefaultAssay(He_endothelial.compact_seurat) <- "RNA"
vln <- VlnPlot(He_endothelial.compact_seurat, cols=c(beach_colors[3],beach_colors[2]), features = c("Cdh5","Kdr","Flt4"), ncol = 3, log=FALSE, combine = FALSE, pt.size=0)
for(i in 1:length(vln)) {
  if ( i == length(vln) ) {
    vln[[i]] <- vln[[i]] + NoAxes() + theme(plot.margin = margin(0, 0, 5, 4, "mm"), plot.title = element_text(size=12,face="plain",vjust = -65 ) ) + stat_compare_means(comparisons=list(levels(He_endothelial.compact_seurat$cluster_names)), method="t.test", label = "p.signif", label.x.npc=0.1, paired=FALSE, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("P<0.0001", "P<0.001", "P<0.01", "P<0.05", "n.s.")), size=2.5 )+ scale_y_continuous(limits=c(-0.1,6))   #comparisons=list(levels(Wang.compact.seurat$cluster_names)), 
  } else {
    vln[[i]] <- vln[[i]] + NoAxes() +NoLegend() + theme(plot.margin = margin(0, 8, 5, 2, "mm"), plot.title = element_text(size=12,face="plain",vjust = -65 ) ) + stat_compare_means(comparisons=list(levels(He_endothelial.compact_seurat$cluster_names)),method="t.test", label = "p.signif", label.x.npc=0.1, paired=FALSE, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("P<0.0001", "P<0.001", "P<0.01", "P<0.05", "n.s.")), size=2.5 )+ scale_y_continuous(limits=c(-0.1,7))   #comparisons=list(levels(Wang.compact.seurat$cluster_names)), 
  }
  vln[[i]]$layers[[1]]$aes_params$size = 0
}
ggsave( plot = wrap_plots(vln,nrow=1) , height=70, width=120, filename = "He_Vlns.pdf", units = "mm",dpi = 300,useDingbats = FALSE)



