##https://satijalab.org/seurat/archive/v3.1/immune_alignment.html
##https://satijalab.org/seurat/articles/multimodal_vignette.html
#https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html (Annotation cluster immune system)
#https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
#https://satijalab.org/seurat/articles/integration_introduction.html

#https://r-charts.com/distribution/violin-plot-group-ggplot2/

#https://www.r-bloggers.com/2013/09/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/

#install.packages('Seurat')
library("Seurat")
library("Matrix")
library("ggplot2")
library("dplyr")
library("viridis")
library("RColorBrewer")
library("data.table")
library("multipanelfigure")

##Setup the Seurat objects
#Day0
data_dir_rnaseq <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-pc-GEX/filtered_feature_bc_matrix/gz/'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
seurat_object_day0 = CreateSeuratObject(counts = data_rnaseq)

data_dir_adt <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/ATD/TD005137_RonaldHoffman/TD005137-pc-ADT_ADT_citeseq/umi_count/gz/'
list.files(data_dir_adt) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_adt <- Read10X(data.dir = data_dir_adt,gene.column = 1,cell.column = 1)
seurat_object_day0[['ADT']] = CreateAssayObject(counts = data_adt)

#Day1 - control
data_dir_rnaseq <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-V-control-day1-GEX/filtered_feature_bc_matrix/gz'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
seurat_object_day1_ctrl = CreateSeuratObject(counts = data_rnaseq)

data_dir_adt <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/ATD/TD005137_RonaldHoffman/TD005137-V-control-day1-ADT_ADT_citeseq/umi_count/gz/'
list.files(data_dir_adt) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_adt <- Read10X(data.dir = data_dir_adt,gene.column = 1,cell.column = 1)
seurat_object_day1_ctrl[['ADT']] = CreateAssayObject(counts = data_adt)

#Day1 - VPA
data_dir_rnaseq <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-V-VPA-day1-GEX/filtered_feature_bc_matrix/gz/'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
seurat_object_day1_vpa = CreateSeuratObject(counts = data_rnaseq)

data_dir_adt <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/ATD/TD005137_RonaldHoffman/TD005137-V-VPA-day1-ADT_ADT_citeseq/umi_count/gz/'
list.files(data_dir_adt) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_adt <- Read10X(data.dir = data_dir_adt,gene.column = 1,cell.column = 1)
seurat_object_day1_vpa[['ADT']] = CreateAssayObject(counts = data_adt)

#Day4 - control
data_dir_rnaseq <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-V-control-day4-GEX/filtered_feature_bc_matrix/gz/'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
dim(data_rnaseq)
seurat_object_day4_ctrl = CreateSeuratObject(counts = data_rnaseq)

data_dir_adt <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/ATD/TD005139_ChristophSchaniel/V-control-day4-ADT_ADT_Ava_093.citeseq/umi_count/gz/'
list.files(data_dir_adt) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_adt <- Read10X(data.dir = data_dir_adt,gene.column = 1,cell.column = 1)
dim(data_adt)
seurat_object_day4_ctrl[['ADT']] = CreateAssayObject(counts = data_adt)

#Day4 - VPA
data_dir_rnaseq <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-V-VPA-day4-GEX/filtered_feature_bc_matrix/gz/'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
seurat_object_day4_vpa = CreateSeuratObject(counts = data_rnaseq)

data_dir_adt <- '/Users/tiphainemartin/Documents/lab_Hoffman/paper2/ATD/TD005139_ChristophSchaniel/V-VPA-day4-ADT_ADT_Ava_093.citeseq/umi_count/gz/'
list.files(data_dir_adt) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_adt <- Read10X(data.dir = data_dir_adt,gene.column = 1,cell.column = 1)
seurat_object_day4_vpa[['ADT']] = CreateAssayObject(counts = data_adt)

#Create a list
seurat_object_day0$treatment <- "day0_none"
seurat_object_day1_vpa$treatment <- "day1_vpa"
seurat_object_day4_vpa$treatment <- "day4_vpa"
seurat_object_day1_ctrl$treatment <- "day1_ctrl"
seurat_object_day4_ctrl$treatment <- "day4_ctrl"

ifnb.list <- list("day0_none"=seurat_object_day0,
                  "day1_vpa"=seurat_object_day1_vpa,
                  "day4_vpa"=seurat_object_day4_vpa,
                  "day1_ctrl"=seurat_object_day1_ctrl,
                  "day4_ctrl"=seurat_object_day4_ctrl)
                  

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

##Perform integration
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, 
                                         anchor.features = features)
# this command creates an 'integrated' data assay (bug)
immune.combined <- IntegrateData(anchorset = immune.anchors)

##Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
Idents(immune.combined) <- immune.combined$seurat_clusters
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "treatment")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

p3 <- DimPlot(immune.combined, reduction = "umap", split.by = "treatment")
p3

##Visualize multiple modalities side-by-side
DefaultAssay(immune.combined) <- "integrated"
# Normalize ADT data,
DefaultAssay(immune.combined) <- "ADT"
immune.combined <- NormalizeData(immune.combined, normalization.method = "CLR",
                                 margin = 2)

# Now, we will visualize CD34, CD90, CD201 levels for RNA and protein level
#By setting the default assay, we can visualize one or the other
####  CD34
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD34-Human-GCAGAAATCTCCCTT", split.by = "treatment",
                  cols = c("#e0e0e0", "#2d004b")) 
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD34-Human-GCAGAAATCTCCCTT",
                  cols = c("#e0e0e0", "#2d004b")) +
  ggtitle("CD34 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "CD34") + ggtitle("CD34 RNA")

# place plots side-by-side
p1 | p2


#### Supplement Fig 1:  CD90
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD90-Human-GCATTGTACGATTCA", split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b")) 
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD90-Human-GCATTGTACGATTCA", 
                  cols = c("#e0e0e0", "#2d004b")) +
  ggtitle("CD90 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "THY1",split.by = "treatment",
                  cols = c("#e0e0e0", "#2d004b")) + ggtitle("CD90 RNA")
p2

# place plots side-by-side
p1 | p2

#### Supplement Fig 1:CD38
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD38-Human-CCTATTCCGATTCCG", split.by = "treatment",
                  cols = c("#e0e0e0", "#2d004b"))
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD38-Human-CCTATTCCGATTCCG", 
                  cols = c("#e0e0e0", "#2d004b")) +
  ggtitle("CD38 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "CD38") + ggtitle("CD38 RNA")

p2_treatment <- FeaturePlot(immune.combined, "CD38", split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b")) + ggtitle("CD38 RNA")
p2_treatment

# place plots side-by-side
p1 | p2

#### Supplement Fig 1: CD201
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD201-Human-GTTTCCTTGACCAAG", split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b"))
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD201-Human-GTTTCCTTGACCAAG", 
                  cols = c("#e0e0e0", "#2d004b")) +
  ggtitle("CD201 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "PROCR") + ggtitle("CD201 RNA")

p2_treatment <- FeaturePlot(immune.combined, "PROCR", split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b"))
p2_treatment

# place plots side-by-side
p1 | p2

#### Supplement Fig 1: CD49F
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD49F-Human-TTCCGAGGATGATCT", split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b"))
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD49F-Human-TTCCGAGGATGATCT", 
                  cols = c("#e0e0e0", "#2d004b")) +
  ggtitle("CD201 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "ITGA6") + ggtitle("CD201 RNA")

p2_treatment <- FeaturePlot(immune.combined, "ITGA6", split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b"))
p2_treatment

# place plots side-by-side
p1 | p2

######## Annotation Cells!
DefaultAssay(immune.combined) <- "ADT"
immune.combined_ADT_data<- GetAssayData(object = immune.combined, slot = "counts")
CD34_ADT <- rep(0,ncol(immune.combined_ADT_data))
CD34_ADT[which(immune.combined_ADT_data[rownames(immune.combined_ADT_data) %in% "CD34-Human-GCAGAAATCTCCCTT",] > 0)] <- 1
CD90_ADT <- rep(0,ncol(immune.combined_ADT_data))
CD90_ADT[which(immune.combined_ADT_data[rownames(immune.combined_ADT_data) %in% "CD90-Human-GCATTGTACGATTCA",] > 0)] <- 1
CD38_ADT <- rep(0,ncol(immune.combined_ADT_data))
CD38_ADT[which(immune.combined_ADT_data[rownames(immune.combined_ADT_data) %in% "CD38-Human-CCTATTCCGATTCCG",] > 0)] <- 1
table(CD38_ADT)
CD201_ADT <- rep(0,ncol(immune.combined_ADT_data))
CD201_ADT[which(immune.combined_ADT_data[rownames(immune.combined_ADT_data) %in% "CD201-Human-GTTTCCTTGACCAAG",] > 0)] <- 1
table(CD201_ADT)


######## Annotation Cells!
DefaultAssay(immune.combined) <- "RNA"
immune.combined_RNA_data<- GetAssayData(object = immune.combined, slot = "counts")
CD34_RNA <- rep(0,ncol(immune.combined_RNA_data))
CD34_RNA[which(immune.combined_RNA_data[rownames(immune.combined_RNA_data) %in% "CD34",] > 0)] <- 1
CD90_RNA <- rep(0,ncol(immune.combined_RNA_data))
CD90_RNA[which(immune.combined_RNA_data[rownames(immune.combined_RNA_data) %in% "THY1",] > 0)] <- 1
CD38_RNA <- rep(0,ncol(immune.combined_RNA_data))
CD38_RNA[which(immune.combined_RNA_data[rownames(immune.combined_RNA_data) %in% "CD38",] > 0)] <- 1
CD201_RNA <- rep(0,ncol(immune.combined_RNA_data))
CD201_RNA[which(immune.combined_RNA_data[rownames(immune.combined_RNA_data) %in% "PROCR",] > 0)] <- 1
table(CD90_RNA)

SIRT1_RNA <- rep(0,ncol(immune.combined_RNA_data))
SIRT1_RNA[which(immune.combined_RNA_data[rownames(immune.combined_RNA_data) %in% "SIRT1",] > 0)] <- 1
table(SIRT1_RNA)
INKA1_RNA <- rep(0,ncol(immune.combined_RNA_data))
INKA1_RNA[which(immune.combined_RNA_data[rownames(immune.combined_RNA_data) %in% "INKA1",] > 0)] <- 1
table(INKA1_RNA)


##CD34_CD90_ hypothese: all cells are CD34+
CD34_CD90_EPCR_RNA <- rep("other",ncol(immune.combined_RNA_data))
CD34_CD90_EPCR_RNA [which( CD90_RNA >0 & CD201_RNA > 0  )] <- "CD34+CD90+EPCR+"
CD34_CD90_EPCR_RNA [which( CD90_RNA == 0 & CD201_RNA == 0 )] <- "CD34+CD90-EPCR-"
table(CD34_CD90_EPCR_RNA)
table(CD34_CD90_EPCR_RNA)*100/ length(CD34_CD90_EPCR_RNA)

CD34_CD90_EPCR_ADT <- rep("other",ncol(immune.combined_ADT_data))
CD34_CD90_EPCR_ADT [which( CD90_ADT > 0 & CD201_ADT > 0 )] <- "CD34+CD90+EPCR+"
CD34_CD90_EPCR_ADT [which( CD90_ADT == 0 & CD201_ADT == 0)] <- "CD34+CD90-EPCR-"
table(CD34_CD90_EPCR_ADT )


#At protein level
DefaultAssay(immune.combined) <- "ADT"
immune.combined$group1 <- CD34_CD90_EPCR_ADT
df_immune <- as.data.frame(cbind(immune.combined$treatment,immune.combined$group1))
colnames(df_immune) <- c("treatment","group1")
table(df_immune$treatment)
#day0_none day1_ctrl  day1_vpa day4_ctrl  day4_vpa 
#11831     11430     14180     13970     12814

#Supplement Fig 1: distribution of groups per treatment at protein level
df_immune %>% dplyr::count(treatment,group1)

((df_immune %>% dplyr::count(treatment,group1))$n/11831*100)[1:3]
((df_immune %>% dplyr::count(treatment,group1))$n/11430*100)[4:6]
((df_immune %>% dplyr::count(treatment,group1))$n/14180*100)[7:9]
((df_immune %>% dplyr::count(treatment,group1))$n/13970*100)[10:12]
((df_immune %>% dplyr::count(treatment,group1))$n/12814*100)[13:15]

#Fig 1: Visualisation of group CD34+CD90+EPCR+/CD34+CD90-EPCR-/other per treatment
# at protein level
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "group1",
              split.by = 'treatment',
              cols = c("darkgreen","purple","#e0e0e0"))
p1


#### Gene expression ####*###**##**#**
DefaultAssay(immune.combined) <- "RNA"
immune.combined$group1 <- CD34_CD90_EPCR_RNA
df_immune <- as.data.frame(cbind(immune.combined$treatment,immune.combined$group1))
colnames(df_immune) <- c("treatment","group1")
table(df_immune$treatment)
#day0_none day1_ctrl  day1_vpa day4_ctrl  day4_vpa 
#11831     11430     14180     13970     12814 

#Supplement Fig 1: distribution of groups per treatment at protein level
df_immune %>% dplyr::count(treatment,group1)

((df_immune %>% dplyr::count(treatment,group1))$n/11831*100)[1:2]
((df_immune %>% dplyr::count(treatment,group1))$n/11430*100)[3:5]
((df_immune %>% dplyr::count(treatment,group1))$n/14180*100)[6:8]
((df_immune %>% dplyr::count(treatment,group1))$n/13970*100)[9:11]
((df_immune %>% dplyr::count(treatment,group1))$n/12814*100)[12:14]

#Fig 1: Visualisation of group CD34+CD90+EPCR+/CD34+CD90-EPCR-/other per treatment
#at gene expression level
DefaultAssay(immune.combined) <- "RNA"
p1 <- DimPlot(immune.combined, reduction = "umap",
              group.by = "group1",split.by = 'treatment',
              cols = c("darkgreen","purple","#e0e0e0"))
p1



########### Heatmap pathway - analysis
#https://blog.bioturing.com/2018/09/24/heatmap-color-scale/
# Fig 3 and supplement Fig 
options(future.globals.maxSize = 4000 * 1024^5)

DefaultAssay(immune.combined) <- "ADT"
immune.combined$group1 <- CD34_CD90_EPCR_ADT

DefaultAssay(immune.combined) <- "RNA"

immune.combined$treatment_group1 <- paste0(immune.combined$group1,"_",immune.combined$treatment)
Idents(immune.combined) <- immune.combined$group1

immune.combined <- ScaleData(object = immune.combined, features = rownames(immune.combined))
Idents(immune.combined) <-  paste0(immune.combined$treatment,"_",immune.combined$group1)

library(dplyr)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library("RColorBrewer")
library("viridis")
library(circlize)
#https://github.com/jokergoo/ComplexHeatmap/issues/148
#https://stackoverflow.com/questions/60543530/getting-error-object-should-be-a-named-list-when-plotting-heatmap-in-r

folder_svg <- "/Users/tiphainemartin/Documents/lab_Hoffman/paper2/results/Day0_Day1_Day4_VPA_control_reviewers/Fig3/"

col_fun = colorRamp2(c(-2, 0, 2), c("#2d004b", "#f7f7f7","#b35806"))

day_colors = c("day0"="#e0ecf4","day1"="#9ebcda","day4"="#8856a7")
treatment_colors = c("none"="#FFFFFF","vpa"="#5ab4ac","ctrl"="#fb9a99")
group_colors = c("CD34+CD90-EPCR-"="#dd1c77","CD34+CD90+EPCR+"="#0570b0")
colz = list(Day = day_colors,
            treatment = treatment_colors,
            Group = group_colors)

####Mitochondrial
markers.to.plot <- c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", 
                     "MT-CO1", "MT-CO2","MT-CO3","MT-CYB","MT-ATP6")
DoHeatmap(immune.combined, features = markers.to.plot) #+ NoLegend()

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",
                             annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)
dev.off()      

#OXPHOS structural genes
markers.to.plot <- c("COX8A", "COX7C", "UQCR11", "NDUFB8", 
                     "ATP5G1", "ATP5G2","NDUFA13","COX7A2","NDUFB2",
                     "ATP5E","UQCRH","ATP5B","NDUFC2",
                     "UQCRB","COX4I1")
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",
                             annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)

#Glycolysis genes (list from Luana)
markers.to.plot <- c("ENO2", "PKM", "HK2", "CITED2", 
                     "SLC16A3", "PDK","PFKM","PGAM","ACLY",
                     "MTCH2","ACAT2")
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)

folder_genes <- "/Users/tiphainemartin/Documents/lab_Hoffman/paper2/data/"

#glycolytic process
#http://www.informatics.jax.org/go/term/GO:0006096
list_genes <- read.table(file=paste0(folder_genes,"GO_term_summary_20220217_001710.txt"),
                         header=TRUE,sep="\t")
markers.to.plot <- toupper(unique(list_genes$Symbol))
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)

#mitochondrial fatty acid beta-oxidation multienzyme complex genes
#http://www.informatics.jax.org/go/term/GO:0016507
markers.to.plot <- c("HADHA", "HADHB")
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)


#FAO genes (luneao)
markers.to.plot <- c("CD36", "CPT1A", "CPT2", "ACAA2", 
                     "ACAT1", "ACAT2","HDAH","SLC25A17")
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)

#glutathione hydrolase activity genes
#http://www.informatics.jax.org/go/term/GO:0036374
markers.to.plot <- c("GGT1", "GGT5","GGT6", "GGT7")
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)

#peroxisome proliferator activated receptor signaling pathway
folder_genes <- "/Users/tiphainemartin/Documents/lab_Hoffman/paper2/data/"
#http://www.informatics.jax.org/go/term/GO:0035357
list_genes <- read.table(file=paste0(folder_genes,"GO_term_summary_20220217_003014.txt"),
                         header=TRUE,sep="\t")
markers.to.plot <- toupper(unique(list_genes$Symbol))
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)

#Antioxidant defense
#http://www.informatics.jax.org/go/term/GO:0016209
list_genes <- read.table(file=paste0(folder_genes,"GO_term_summary_20220217_000048.txt"),
                         header=TRUE,sep="\t")
markers.to.plot <- toupper(unique(list_genes$Symbol))
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)


#TCA cycle
#http://www.informatics.jax.org/go/term/GO:0006099
list_genes <- read.table(file=paste0(folder_genes,"GO_term_summary_20220217_001423.txt"),
                         header=TRUE,sep="\t")
markers.to.plot <- toupper(unique(list_genes$Symbol))
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)


#Cell cycle
#### Need lot of memory
#http://www.informatics.jax.org/go/term/GO:0007049
# list_genes <- read.table(file=paste0(folder_genes,"GO_term_summary_20220217_005621.txt"),
#                          header=TRUE,sep="\t")
# markers.to.plot <- toupper(unique(list_genes$Symbol))
# DoHeatmap(immune.combined, features = markers.to.plot)
# 
# visu <- DoHeatmap(immune.combined, features = markers.to.plot)
# data_tmp <- visu$data
# data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
# data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
# wide = data_short[,1:3] %>% 
#   spread(Cell, Expression)
# #wide
# rownames(wide) <- wide[,1]
# wide <- wide[,-1]
# 
# annotation <- unique(data_short[,c(2,4:6)])
# annotation$general <- paste0(annotation$day,"_",
#                              annotation$treatment,"_",annotation$group)
# 
# ha = HeatmapAnnotation(Day = annotation$day, 
#                        treatment = annotation$treatment,
#                        Group = annotation$group,
#                        col=colz)
# Heatmap(wide, name = "Zscore",
#         column_order = order(annotation$day,annotation$treatment,annotation$group),
#         show_column_names = FALSE, top_annotation = ha,
#         show_column_dend = FALSE, show_row_dend = FALSE,
#         col = col_fun, column_split = annotation$general,
#         column_title = NULL)

#Fig 3 B
markers.to.plot <- c("PKM", "MTCH2","ACLY", "ACAT2","CITED2","SLC16A3",
                     "EN02","HK2","PFKM")
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)

#Fig 3 E
markers.to.plot <- c("ACAT2", "ACAA2","ACAT1", "CPT1A","CPT2",
                     "SLC25A17","CD36","HADHA","HADHB")
DoHeatmap(immune.combined, features = markers.to.plot)

visu <- DoHeatmap(immune.combined, features = markers.to.plot)
data_tmp <- visu$data
data <- data_tmp %>% tidyr::separate(Identity, c("day", "treatment","group"),sep="_")
data_short <- data[-which(data$group %in% c("other","CD34+CD90-EPCR-")),]
wide = data_short[,1:3] %>% 
  spread(Cell, Expression)
#wide
rownames(wide) <- wide[,1]
wide <- wide[,-1]

annotation <- unique(data_short[,c(2,4:6)])
annotation$general <- paste0(annotation$day,"_",annotation$treatment,"_",annotation$group)

ha = HeatmapAnnotation(Day = annotation$day, 
                       treatment = annotation$treatment,
                       Group = annotation$group,
                       col=colz)
Heatmap(wide, name = "Zscore",
        column_order = order(annotation$day,annotation$treatment,annotation$group),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun, column_split = annotation$general,
        column_title = NULL)
