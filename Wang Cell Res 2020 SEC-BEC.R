#############
# Begin
#############
#load Seurat
library(Seurat)
library(patchwork)

#load generic libraries
library(readxl)
library(tidyverse)
library(viridis)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)


#############
#Preferences
#############
seurat_umap_alpha_hex_string = "d0"
data_path = "/home/martin/Documents/20-22/Derek scRNAseq 2022 analysis/Wang Cell Res 2020/Process"
beach_colors <- c("#1a748e","#f0df99","#55c4d7","#5f6c24","#992915","#d38e31","#81c7f8","#bbc4af","#393430","#92dccd","#63fba9","#1155d4","#abb2ba","#eccd16","#7b5c52","#063581","#ecdfcf")
paired_colors <- brewer.pal(name = "Paired", n=12)
beach_colors_alpha <- paste0( beach_colors, seurat_umap_alpha_hex_string )
paired_colors_alpha <- paste0( paired_colors, seurat_umap_alpha_hex_string )

setwd(data_path)


##############################################################################
#
#     RAW DATA INPUT -- download raw data from Wang et al. Cell Research 2020
#
##############################################################################

download.file("ftp://download.cncb.ac.cn/OMIX/OMIX895/OMIX895-05-04.xlsx",method="curl",destfile="OMIX895-05-04.xlsx")
cell_key <- as.data.frame(read_xlsx("OMIX895-05-04.xlsx",sheet=1,skip=1))
rownames(cell_key) <- cell_key$SampleName
cell_key = cell_key[,which(names(cell_key)!='SampleName')]
#cell_key = read.csv(file = "OMIX895-05-04_modified.tsv", stringsAsFactors = FALSE, header=T, sep="\t" , row.names = 1)
cell_key <- cell_key[ which(cell_key$cluster=='Liver endothelial cell'), ]

download.file("ftp://download.cncb.ac.cn/OMIX/OMIX895/OMIX895-20-01.5.zip",method="curl",destfile="E17.5.zip")
download.file("ftp://download.cncb.ac.cn/OMIX/OMIX895/OMIX895-20-02.0.zip",method="curl",destfile="E11.0.zip")
download.file("ftp://download.cncb.ac.cn/OMIX/OMIX895/OMIX895-20-03.0.zip",method="curl",destfile="E13.0.zip")
download.file("ftp://download.cncb.ac.cn/OMIX/OMIX895/OMIX895-20-06.5.zip",method="curl",destfile="E11.5.zip")
download.file("ftp://download.cncb.ac.cn/OMIX/OMIX895/OMIX895-20-07.0.zip",method="curl",destfile="E16.0.zip")
download.file("ftp://download.cncb.ac.cn/OMIX/OMIX895/OMIX895-20-08.5.zip",method="curl",destfile="E14.5.zip")

sample_list = c("E11.0","E11.5","E13.0","E14.5","E16.0","E17.5")

lapply(X = sample_list, FUN = function(x) {
  unzip( paste0( x, ".zip") )
  dir.create( x )
  file.rename( from=paste0(x,"_barcodes.tsv"), to=file.path(x,"barcodes.tsv"))
  file.rename( from=paste0(x,"_matrix.mtx"), to=file.path(x,"matrix.mtx"))
  file.rename( from=paste0(x,"_genes.tsv"), to=file.path(x,"genes.tsv"))
})


#############
#Read in all the 10X data
#############
load_10x_data_from_folders = function(folders, key=cell_key, path=data_path ){

  sample_name = folders[1]
  folder = file.path(path, sample_name)
  print( sample_name )
  seurat.obj.a <- Read10X(folder)
  colnames(seurat.obj.a) = paste(sample_name, matrix(unlist(strsplit(colnames(seurat.obj.a), "-")),ncol=2,byrow=TRUE,dimnames = NULL)[,1], sep = "_")
  seurat.obj = CreateSeuratObject(counts = seurat.obj.a, min.cells = 0, min.features = 0, project = sample_name)
  cells.use <- rownames(key)
  seurat.obj <- subset(seurat.obj, cells = cells.use)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
  VlnPlot( seurat.obj, features=c("nFeature_RNA","nCount_RNA","percent.mt") )
  # if more the 1 folder load next folder(s)
  for (i in 2:length(folders)){
    sample_name = folders[i]

    folder = file.path(path, sample_name)
    print( sample_name )
    seurat.obj.b <- Read10X(folder)
    colnames(seurat.obj.b) = paste(sample_name, matrix(unlist(strsplit(colnames(seurat.obj.b), "-")),ncol=2,byrow=TRUE,dimnames = NULL)[,1], sep = "_")
    seurat.new = CreateSeuratObject(counts = seurat.obj.b, min.cells = 0, min.features = 0, project = sample_name )
    seurat.new <- subset(seurat.new, cells = cells.use)
    seurat.new[["percent.mt"]] <- PercentageFeatureSet(seurat.new, pattern = "^mt-")
    VlnPlot( seurat.new, features=c("nFeature_RNA","nCount_RNA","percent.mt") )
    seurat.obj <- merge(x = seurat.obj, y = seurat.new)
  }

  return(seurat.obj)
}

split_and_integrate_seurat_object_CCA = function(object=object,min.cells=50,split_by = "orig.ident") {
  Wang.list <- SplitObject(object, split.by = split_by )
  remove(object)

  #Remove sub-objects with too few cells
  for ( i in length(Wang.list):1 ) {
    if ( length(colnames(Wang.list[[i]])) < min.cells ) {
      Wang.list[[i]] <- NULL
    }
  }

  Wang.list <- lapply(X = Wang.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

  features <- SelectIntegrationFeatures(object.list = Wang.list)
  Wang.anchors <- FindIntegrationAnchors(object.list = Wang.list, anchor.features = features )
  Wang.combined <- IntegrateData(anchorset = Wang.anchors)
  DefaultAssay(Wang.combined) <- "integrated"
  remove(Wang.anchors)
  return(Wang.combined)

}


#############
#Load data, split objects by sample, integrate, and run the standard workflow for ECs only
#############
Wang.seurat <- load_10x_data_from_folders( sample_list, key = cell_key )
Wang.seurat <- split_and_integrate_seurat_object_CCA( object = Wang.seurat )

save(Wang.seurat,file="Wang_endothelial.seurat.Rds")

##############################################################################
#
#     END RAW DATA INPUT
#
##############################################################################


#############
#Standard workflow
#############

load("Wang_endothelial.seurat.Rds")

Idents(Wang.seurat) <- "orig.ident"
Wang.seurat <- Wang.compact.seurat <- subset( Wang.seurat, idents = c("E11.5","E13.0","E14.5","E16.0","E17.5") )

DefaultAssay(Wang.seurat) <- "integrated"

Wang.seurat <- ScaleData(Wang.seurat, verbose = FALSE)
Wang.seurat <- FindVariableFeatures(Wang.seurat)

Wang.seurat <- RunPCA(Wang.seurat, npcs = 30, verbose = FALSE)
ElbowPlot( Wang.seurat, ndims=30)
Wang.seurat <- RunUMAP(Wang.seurat, reduction = "pca", dims = 1:15)
Wang.seurat <- FindNeighbors(Wang.seurat, reduction = "pca", dims = 1:15)
Wang.seurat <- FindClusters(Wang.seurat, resolution = 0.25)
#p2 <- DimPlot(Wang.seurat,  reduction = "umap", group.by = "seurat_clusters", pt.size=0.25, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
#p3 <- DimPlot(Wang.seurat, reduction = "umap", group.by = "orig.ident", pt.size=0.25) + NoAxes() + labs(title = "tissue") + theme(plot.title = element_text(hjust = 0.5))
#p2 +p3

#all.markers <- FindAllMarkers(object = Wang.seurat)
#top20 <- all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
#DoHeatmap(object = Wang.seurat, features = top20$gene, label=FALSE) #slim.col.label = TRUE, remove.key = TRUE)



#############
#Cluster naming
#############
Idents(Wang.seurat) <- "seurat_clusters"
Wang.seurat[["old.ident"]] <- Idents(Wang.seurat)
Wang.seurat <- RenameIdents(object = Wang.seurat, 
                            '0' = "SECs", 
                            '1' = "Proliferative SECs", 
                            '3' = "Hematopoetic Cells", 
                            '4' = "Central Vein BECs", 
                            '2' = "Arterial BECs", 
                            '5' = "Mesenchymal Cells"
)
Wang.seurat[["cluster_names"]] <- Idents(Wang.seurat)

#############
#Visualize final object
#############
p1 <- DimPlot(Wang.seurat, reduction = "umap", cols=paired_colors_alpha, group.by = "orig.ident", pt.size=0.65) + NoAxes() + labs(title = "Stage") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(Wang.seurat,  reduction = "umap", cols=beach_colors_alpha, group.by = "cluster_names", pt.size=0.65, label=FALSE) + NoAxes() + labs(title = "Cluster") + theme(plot.title = element_text(hjust = 0.5))
p1 + p2
#ggsave( plot = p1 + p2, height=200, width=500, filename = "Wang_complete_UMAP.png", units = "mm",dpi = 300 )

#all.markers <- FindAllMarkers(object = Wang.seurat)
#top20 <- all.markers %>% group_by(cluster) %>% top_n(8, avg_log2FC*pct.1*pct.1)
#HM <- DoHeatmap(object = Wang.seurat, features = top20$gene, label=FALSE, group.colors = beach_colors)  + scale_fill_gradientn(colors = c("#440154FF", "#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")) + theme(legend.title = element_text(size = 22),legend.text = element_text(size = 20),axis.text.y = element_text(size = 15),legend.key.size = unit(1.5,"line")) + NoLegend() #slim.col.label = TRUE, remove.key = TRUE)
#ggsave( plot = HM, height=500, width=400, filename = "Wang_complete_clusterHM.png", units = "mm",dpi = 300 )

DefaultAssay(Wang.seurat) <- "RNA"
my.markers <- c("Pecam1","Emcn","Rspo3","Gja5","Efnb2","Acvrl1","Flt1","Sox17","Hey1","Stab1","Stab2","Lyve1","Prox1","Pdpn","Cdk1","Mki67","Pdgfrb")
DP <- DotPlot(Wang.seurat, features = my.markers, dot.scale = 8 ) + # cols = c("Blue","Orange")) +
  scale_y_discrete(limits=rev) + RotatedAxis() + theme(axis.title.y = element_blank(),axis.title.x = element_blank())

#ggsave( plot = ((p1+p2)/HM) + plot_layout(heights=c(1,3)), height=360, width=300, filename = "Wang_supplement.pdf", units = "mm",dpi = 300 )
ggsave( plot = ((p1+p2)/DP) + plot_layout(heights=c(2,1)), height=200, width=300, filename = "Wang_supplement_DP.pdf", units = "mm",dpi = 300 )




#############
#Cluster compaction
#############
Wang.compact.seurat <- subset( Wang.seurat, idents = c("SECs","Proliferative SECs","Central Vein BECs","Arterial BECs") )
Wang.compact.seurat[["cluster_names"]] <- Idents(Wang.compact.seurat)
Wang.compact.seurat <- RenameIdents(object = Wang.compact.seurat, 
                                    "SECs" = "SECs", 
                                    "Proliferative SECs" = "SECs", 
                                    "Central Vein BECs" = "BECs", 
                                    "Arterial BECs" = "BECs"
)
Wang.compact.seurat[["cluster_names"]] <- Idents(Wang.compact.seurat)
                            
 
#FeaturePlot(Wang.seurat,features=c(my.markers,"Rspo3","Nts","Lhx6","Fbln2","Dll4"))
#DefaultAssay(Wang.compact.seurat) <- "RNA"
#VlnPlot(Wang.compact.seurat, group.by="cluster_names", features=c("Vwf","Lyve1","Kdr","Flt4","Cdh5","Nr2f2","Stab2","Stab1","Bcam","Pecam1","Cd34","Sox17","Emcn"))
#markers <- FindMarkers( Wang.seurat,ident.1 = "Venous BECs" )


#############
#Violin for BEC/SEC/LEC
#############
DefaultAssay(Wang.compact.seurat) <- "RNA"
levels(Wang.compact.seurat) <- sort(levels(Wang.compact.seurat))
vln <- VlnPlot(Wang.compact.seurat, cols=c(beach_colors[3],beach_colors[1]), features = c("Cdh5","Kdr","Flt4"), ncol = 3, log=FALSE, combine = FALSE, pt.size=0)
for(i in 1:length(vln)) {
  if ( i == length(vln) ) {
    vln[[i]] <- vln[[i]] + NoAxes() + theme(plot.margin = margin(0, 0, 5, 4, "mm"), plot.title = element_text(size=12,face="plain",vjust = -65 ) ) + stat_compare_means(comparisons=list(levels(Wang.compact.seurat$cluster_names)), method="t.test", label = "p.signif", label.x.npc=0.1, paired=FALSE, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("P<0.0001", "P<0.001", "P<0.01", "P<0.05", "n.s.")), size=2.5 )+ scale_y_continuous(limits=c(-0.1,3))   #comparisons=list(levels(Wang.compact.seurat$cluster_names)), 
  } else {
    vln[[i]] <- vln[[i]] + NoAxes() +NoLegend() + theme(plot.margin = margin(0, 8, 5, 2, "mm"), plot.title = element_text(size=12,face="plain",vjust = -65 ) ) + stat_compare_means(comparisons=list(levels(Wang.compact.seurat$cluster_names)),method="t.test", label = "p.signif", label.x.npc=0.1, paired=FALSE, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("P<0.0001", "P<0.001", "P<0.01", "P<0.05", "n.s.")), size=2.5 )+ scale_y_continuous(limits=c(-0.1,3.5))   #comparisons=list(levels(Wang.compact.seurat$cluster_names)), 
  }
  vln[[i]]$layers[[1]]$aes_params$size = 0
}
#wrap_plots(vln,nrow=1)
ggsave( plot = wrap_plots(vln,nrow=1) , height=70, width=120, filename = "Wang_Vlns.pdf", units = "mm",dpi = 300,useDingbats = FALSE)
