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

folder_results = "/lab_Hoffman/paper2/results/Day0_Day1_Day4_VPA_control_reviewers/Fig2/"
##Setup the Seurat objects
#Day0
data_dir_rnaseq <- '/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-pc-GEX/filtered_feature_bc_matrix/gz/'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
seurat_object_day0 = CreateSeuratObject(counts = data_rnaseq)

data_dir_adt <- '/lab_Hoffman/paper2/ATD/TD005137_RonaldHoffman/TD005137-pc-ADT_ADT_citeseq/umi_count/gz/'
list.files(data_dir_adt) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_adt <- Read10X(data.dir = data_dir_adt,gene.column = 1,cell.column = 1)
seurat_object_day0[['ADT']] = CreateAssayObject(counts = data_adt)

#Day1 - control
data_dir_rnaseq <- '/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-V-control-day1-GEX/filtered_feature_bc_matrix/gz'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
seurat_object_day1_ctrl = CreateSeuratObject(counts = data_rnaseq)

data_dir_adt <- '/lab_Hoffman/paper2/ATD/TD005137_RonaldHoffman/TD005137-V-control-day1-ADT_ADT_citeseq/umi_count/gz/'
list.files(data_dir_adt) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_adt <- Read10X(data.dir = data_dir_adt,gene.column = 1,cell.column = 1)
seurat_object_day1_ctrl[['ADT']] = CreateAssayObject(counts = data_adt)

#Day1 - VPA
data_dir_rnaseq <- '/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-V-VPA-day1-GEX/filtered_feature_bc_matrix/gz/'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
seurat_object_day1_vpa = CreateSeuratObject(couns = data_rnaseq)

data_dir_adt <- '/lab_Hoffman/paper2/ATD/TD005137_RonaldHoffman/TD005137-V-VPA-day1-ADT_ADT_citeseq/umi_count/gz/'
list.files(data_dir_adt) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_adt <- Read10X(data.dir = data_dir_adt,gene.column = 1,cell.column = 1)
seurat_object_day1_vpa[['ADT']] = CreateAssayObject(counts = data_adt)

#Day4 - control
data_dir_rnaseq <- '/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-V-control-day4-GEX/filtered_feature_bc_matrix/gz/'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
dim(data_rnaseq)
seurat_object_day4_ctrl = CreateSeuratObject(counts = data_rnaseq)

data_dir_adt <- '/lab_Hoffman/paper2/ATD/TD005139_ChristophSchaniel/V-control-day4-ADT_ADT_Ava_093.citeseq/umi_count/gz/'
list.files(data_dir_adt) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_adt <- Read10X(data.dir = data_dir_adt,gene.column = 1,cell.column = 1)
dim(data_adt)
seurat_object_day4_ctrl[['ADT']] = CreateAssayObject(counts = data_adt)

#Day4 - VPA
data_dir_rnaseq <- '/lab_Hoffman/paper2/rnaseq/TD005069_RonaldHoffman/TD005069-V-VPA-day4-GEX/filtered_feature_bc_matrix/gz/'
list.files(data_dir_rnaseq) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data_rnaseq <- Read10X(data.dir = data_dir_rnaseq)
colnames(data_rnaseq)<-gsub("-1","",colnames(data_rnaseq))
seurat_object_day4_vpa = CreateSeuratObject(counts = data_rnaseq)

data_dir_adt <- '/lab_Hoffman/paper2/ATD/TD005139_ChristophSchaniel/V-VPA-day4-ADT_ADT_Ava_093.citeseq/umi_count/gz/'
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
immune.combined <- NormalizeData(immune.combined, normalization.method = "CLR", margin = 2)

# Now, we will visualize CD34, CD90, CD201 levels for RNA and protein level
#By setting the default assay, we can visualize one or the other
####  CD34
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD34-Human-GCAGAAATCTCCCTT", split.by = "treatment",
                            cols = c("lightgrey", "#542788")) 
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD34-Human-GCAGAAATCTCCCTT",
                  cols = c("lightgrey", "#542788")) +
  ggtitle("CD34 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "CD34") + ggtitle("CD34 RNA")

# place plots side-by-side
p1 | p2


#### Supplement Fig 1:  CD90
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD90-Human-GCATTGTACGATTCA", split.by = "treatment",
                            cols = c("lightgrey", "#542788")) 
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD90-Human-GCATTGTACGATTCA", 
                  cols = c("lightgrey", "#542788")) +
  ggtitle("CD90 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "THY1",split.by = "treatment",
                  cols = c("lightgrey", "#542788")) + ggtitle("CD90 RNA")
p2

# place plots side-by-side
p1 | p2

#### Supplement Fig 1:CD38
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD38-Human-CCTATTCCGATTCCG", split.by = "treatment",
                            cols = c("lightgrey", "#542788"))
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD38-Human-CCTATTCCGATTCCG", 
                  cols = c("lightgrey", "#542788")) +
  ggtitle("CD38 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "CD38") + ggtitle("CD38 RNA")

p2_treatment <- FeaturePlot(immune.combined, "CD38", split.by = "treatment",
                            cols = c("lightgrey", "#542788")) + ggtitle("CD38 RNA")
p2_treatment

# place plots side-by-side
p1 | p2

#### Supplement Fig 1: CD201
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD201-Human-GTTTCCTTGACCAAG", split.by = "treatment",
                            cols = c("lightgrey", "#542788"))
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD201-Human-GTTTCCTTGACCAAG", 
                  cols = c("lightgrey", "#542788")) +
  ggtitle("CD201 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "PROCR") + ggtitle("CD201 RNA")

p2_treatment <- FeaturePlot(immune.combined, "PROCR", split.by = "treatment",
                            cols = c("lightgrey", "#542788"))
p2_treatment

# place plots side-by-side
p1 | p2

#### Supplement Fig 1: CD49F
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD49F-Human-TTCCGAGGATGATCT", split.by = "treatment",
                            cols = c("lightgrey", "#542788"))
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD49F-Human-TTCCGAGGATGATCT", 
                  cols = c("lightgrey", "#542788")) +
  ggtitle("CD201 protein")
DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(immune.combined, "ITGA6") + ggtitle("CD201 RNA")

p2_treatment <- FeaturePlot(immune.combined, "ITGA6", split.by = "treatment",
                            cols = c("lightgrey", "#542788"))
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

##CD34_CD90_EPCR_CD38 hypothese: all cells are CD34+

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


##CD34_CD90_EPCR_CD38 hypothese: all cells are CD34+ RNA level
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
              cols = c("darkgreen","purple","lightgrey"))
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
              cols = c("darkgreen","purple","lightgrey"))
p1

  


###
### Cell cycle
#https://satijalab.org/seurat/articles/cell_cycle_vignette.html

# A list of cell cycle markers, from Kowalczyk MS*., Tirosh I* et al. (2015). 
#"Single-cell RNA-seq reveals changes in cell cycle and differentiation programs upon aging 
#of hematopoietic stem cells." Genome Res 25(12): 1860-1872. [PDF], is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
#https://www.science.org/doi/10.1126/science.aad0501
#We score single cells based on the scoring strategy described in Tirosh et al. 2016.
# https://genome.cshlp.org/content/25/12/1860.long
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


# Create our Seurat object and complete the initalization steps
immune.combined <- NormalizeData(immune.combined)
immune.combined <- FindVariableFeatures(immune.combined, selection.method = "vst")
immune.combined <- ScaleData(immune.combined, features = rownames(immune.combined))

#Assign Cell-Cycle Scores
immune.combined <- CellCycleScoring(immune.combined, 
                                    s.features = s.genes, g2m.features = g2m.genes,
                                    set.ident = TRUE)
# view cell cycle scores and phase assignments
head(immune.combined[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(immune.combined, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
immune.combined <- RunPCA(immune.combined, features = c(s.genes, g2m.genes))
immune.combined$group1_ADT <- CD34_CD90_EPCR_ADT
DimPlot(immune.combined)
DimPlot(immune.combined,split.by = 'treatment')
DimPlot(immune.combined,group.by = "group1_ADT",
        split.by = 'treatment')
##
# DefaultAssay(immune.combined) <- "integrated"
# immune.combined_integrated <- RunPCA(immune.combined, features = c(s.genes, g2m.genes))
# immune.combined_integrated$group1_ADT <- CD34_CD90_EPCR_ADT
# DimPlot(immune.combined_integrated,reduction = "pca")
# DimPlot(immune.combined_integrated,split.by = 'treatment',reduction = "pca")
# DimPlot(immune.combined_integrated,group.by = "group1_ADT",
#         split.by = 'treatment',reduction = "pca")
# DimPlot(immune.combined_integrated,group.by = "Phase",
#         split.by = 'treatment',reduction = "pca")

data_cellcycle <- cbind(immune.combined$treatment,immune.combined$Phase,
                        immune.combined$group1_ADT)
colnames(data_cellcycle) <- c("treatment","Phase","ADT")
data_cellcycle_df <- as.data.frame(data_cellcycle)
data_cellcycle_df$treatment_ADT <- paste0(immune.combined$treatment,"_",
                                                 immune.combined$group1_ADT)
##Per treatment and cluster and Phase
data_cellcycle_result <- data_cellcycle_df %>% count(treatment, Phase, ADT, sort = TRUE)
data_cellcycle_result_df <- as.data.frame(data_cellcycle_result)

data_cellcycle_result_df$Phase[which(data_cellcycle_result_df$Phase %in% "G1")] <- "1_G1/G0"
data_cellcycle_result_df$Phase[which(data_cellcycle_result_df$Phase %in% "S")] <- "2_S"
data_cellcycle_result_df$Phase[which(data_cellcycle_result_df$Phase %in% "G2M")] <- "3_G2M"

data_cellcycle_result_df$treatment_ADT <- paste0(data_cellcycle_result_df$treatment,"_",
                                                 data_cellcycle_result_df$ADT)

file_cc=paste0(folder_results,"CellCycle_treatment_CD34_CD90_EPCR_proteins_number.csv")
write.table(data_cellcycle_result_df,file=file_cc,sep=",",
            quote=FALSE,row.names = FALSE)

##New visualisation
data_cellcycle_result_df <- as.data.frame(data_cellcycle_result)
data_cellcycle_result_df$Phase[which(data_cellcycle_result_df$Phase %in% "G1")] <- "3_G1/G0"
data_cellcycle_result_df$Phase[which(data_cellcycle_result_df$Phase %in% "S")] <- "2_S"
data_cellcycle_result_df$Phase[which(data_cellcycle_result_df$Phase %in% "G2M")] <- "1_G2M"
data_cellcycle_result_df$treatment_ADT <- paste0(data_cellcycle_result_df$treatment,"_",
                                                 data_cellcycle_result_df$ADT)
data_cellcycle_result_df$Phase_treatment <- paste0(data_cellcycle_result_df$Phase,"_",
                                                   data_cellcycle_result_df$treatment)

#scale_fill_manual(values=c("#F8766D","#619CFF","#00BA38")) +
#https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
# Stacked + number
ggplot(data_cellcycle_result_df, aes(fill=Phase, y=n, x=treatment_ADT)) + 
  geom_bar(position="stack", stat="identity") +
  theme_light() + 
  scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

# Stacked + percent
ggplot(data_cellcycle_result_df, aes(fill=Phase, y=n, x=treatment_ADT)) + 
  geom_bar(position="fill", stat="identity") +
  theme_light() +
  scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

data_cellcycle_result_df$Phase[which(data_cellcycle_result_df$Phase %in% "3_G1/G0")] <- "1_G1/G0"
data_cellcycle_result_df$Phase[which(data_cellcycle_result_df$Phase %in% "1_G2M")] <- "3_G2M"

data_cellcycle_result_df$Phase_treatment <- paste0(data_cellcycle_result_df$Phase,"_",
                                                   data_cellcycle_result_df$treatment)

# Stacked + number
ggplot(data_cellcycle_result_df, aes(fill=ADT, y=n, x=Phase_treatment)) + 
  geom_bar(position="stack", stat="identity") +
  theme_light() + 
  scale_fill_manual(values=c("lightgrey","purple","#b6e2d3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

# Stacked + percent
ggplot(data_cellcycle_result_df, aes(fill=ADT, y=n, x=Phase_treatment)) + 
  geom_bar(position="fill", stat="identity") +
  theme_light() +
  scale_fill_manual(values=c("lightgrey","purple","#b6e2d3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

# Stacked + number
ggplot(data_cellcycle_result_df, aes(fill=ADT, y=n, x=Phase)) + 
  geom_bar(position="stack", stat="identity") +
  theme_light() + 
  scale_fill_manual(values=c("lightgrey","purple","#b6e2d3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

# Stacked + percent
ggplot(data_cellcycle_result_df, aes(fill=ADT, y=n, x=Phase)) + 
  geom_bar(position="fill", stat="identity") +
  theme_light() +
  scale_fill_manual(values=c("lightgrey","purple","#b6e2d3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


####### Violin plot for "SIRT1", "SIRT3", "INKA1"
threshold = 0
immune.combined$PhaseNum <-immune.combined$Phase
immune.combined$PhaseNum[which(immune.combined$PhaseNum %in% "G1")] <- "1_G0/G1"
immune.combined$PhaseNum[which(immune.combined$PhaseNum %in% "S")] <- "2_S"
immune.combined$PhaseNum[which(immune.combined$PhaseNum %in% "G2M")] <- "3_G2M"

immune.combined$treatment_group1_phase <- paste0(immune.combined$treatment,"_",
                                                 immune.combined$group1_ADT,"_",
                                                 immune.combined$PhaseNum)
Idents(immune.combined) <- immune.combined$treatment_group1_phase

markers.to.plot <- c("SIRT1", "SIRT3", "INKA1")

DotPlot(immune.combined, features = markers.to.plot, 
        cols=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"),
        dot.scale = 8,  split.by = "treatment") +
  RotatedAxis()

FeaturePlot(immune.combined, 
            features = markers.to.plot, 
            split.by = "treatment", max.cutoff = 3,
            cols = c("#e0e0e0", "#253494"))

FeaturePlot(immune.combined, features = markers.to.plot, min.cutoff = "q9")

info_dot<- DotPlot(immune.combined, features = markers.to.plot, 
                   cols=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"),
                   dot.scale = 8,  split.by = "treatment")

PoI <- c("SIRT1", "SIRT3", "INKA1")

immune.combined_RNA_data_scale<- GetAssayData(object = immune.combined, slot = "data")

CD34_CD90_EPCR_ADT_Phase <- paste0(CD34_CD90_EPCR_ADT,"_",immune.combined$PhaseNum)
## SIRT1
i<-1
PoI[i]

###Define the range to work on
DefaultAssay(immune.combined) <- "RNA"
data_RNA_gene<- immune.combined_RNA_data_scale[rownames(immune.combined_RNA_data_scale) %in% PoI[i],]
(first_quantile <- quantile(data_RNA_gene, 0.25)) # first quartile
(third_quantile <- quantile(data_RNA_gene, 0.75)) # third quartile
(mean_gene <-  mean(data_RNA_gene))
(max_gene <-  max(data_RNA_gene))


SIRT1_RNA <- rep("",ncol(immune.combined_RNA_data_scale))
SIRT1_RNA[which(data_RNA_gene ==  threshold &
                  CD34_CD90_EPCR_ADT %in% "CD34+CD90+EPCR+" )] <- "SIRT1_0_low"
SIRT1_RNA[which(data_RNA_gene > threshold &
                  CD34_CD90_EPCR_ADT %in% "CD34+CD90+EPCR+" )] <- "SIRT1_2_high"

table(SIRT1_RNA)

## create cluster
CD34_CD90_EPCR_ADT_Phase_SIRT1 <- paste0(CD34_CD90_EPCR_ADT_Phase,"_",SIRT1_RNA)
immune.combined$CD34_CD90_EPCR_ADT_Phase_SIRT1 <- CD34_CD90_EPCR_ADT_Phase_SIRT1
CD34_CD90_EPCR_ADT_Phase_SIRT1_day_tmp <- paste0(CD34_CD90_EPCR_ADT_Phase_SIRT1,
                                                 "_",immune.combined$treatment)
data_SIRT1 <- table(CD34_CD90_EPCR_ADT_Phase_SIRT1_day_tmp)
length(table(CD34_CD90_EPCR_ADT_Phase_SIRT1_day_tmp))

Idents(immune.combined) <- immune.combined$CD34_CD90_EPCR_ADT_Phase_SIRT1

###visualisation SIRT1
data_SIRT1_df_tmp <- as.data.frame(data_SIRT1)
data_SIRT1_df_tmp$features.plot <- "SIRT1"
colnames(data_SIRT1_df_tmp)[1]<-"id"

data_SIRT1_df_tmp2 <- data_SIRT1_df_tmp[-which(data_SIRT1_df_tmp$id %like% c("other")),]
data_SIRT1_df_tmp3 <- data_SIRT1_df_tmp2[-which(data_SIRT1_df_tmp2$id %like% c("CD90-")),]
data_SIRT1_df<- data_SIRT1_df_tmp3[-which(data_SIRT1_df_tmp3$id %like% c("low")),]

dot.PoI.num <-ggplot(data_SIRT1_df, aes(x = id, y = features.plot,
                               size = Freq, fill = Freq)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low="#e0ecf4", high="#762a83",limits=c(0,100))+
  scale_size(limits=c(0,5000)) +
  scale_x_discrete(breaks=data_SIRT1_df$id,
                   labels=data_SIRT1_df$Freq)+
  theme_minimal() + labs(y="", x="") +
  theme(legend.position="none")

info_dot<- DotPlot(immune.combined, features = PoI, 
                   cols=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"),
                   dot.scale = 8,  split.by = "treatment")

data_dot <- info_dot$data
dot_spec_tmp <- data_dot[-which(data_dot$id %like% c("other")),]
dot_spec_tmp2 <- dot_spec_tmp[-which(dot_spec_tmp$id %like% c("CD90-")),]
dot_spec_tmp3 <- dot_spec_tmp2[-which(dot_spec_tmp2$id %like% c("low")),]
dot_spec <- dot_spec_tmp2[which(dot_spec_tmp2$features.plot %in% PoI[i]),]

dot_spec$features.plot <- "SIRT1"
dot.PoI <-ggplot(dot_spec, aes(x = id, y = features.plot,
                               size = pct.exp, fill = pct.exp)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low="#e0ecf4", high="#762a83",limits=c(0,100))+
  scale_size(limits=c(0,100)) +
  scale_x_discrete(breaks=dot_spec$id,
                   labels=paste0(round(dot_spec$pct.exp,2),"%"))+
  theme_minimal() + labs(y="", x="") +
  theme(legend.position="none")


data<-VlnPlot(object = immune.combined, features = PoI[i],
              group.by ='treatment',split.by='CD34_CD90_EPCR_ADT_Phase_SIRT1')

info.avg.PoI_tmp<- data$data %>% group_by(split,ident) %>% dplyr::summarise(n = n(), 
                                                                            average = mean(SIRT1), 
                                                                            maximum = max(SIRT1))
info.avg.PoI_tmp$id <- paste0(info.avg.PoI_tmp$split,"_",info.avg.PoI_tmp$ident)
info.avg.PoI_tmp$features.plot <- "SIRT1"
info.avg.PoI_tmp2 <- info.avg.PoI_tmp[-which(info.avg.PoI_tmp$split %like% c("other_")),]
info.avg.PoI_tmp3 <- info.avg.PoI_tmp2[-which(info.avg.PoI_tmp2$split %like% c("CD90-")),]
info.avg.PoI <- info.avg.PoI_tmp3[-which(info.avg.PoI_tmp3$split %like% c("low")),]

dot.avg.PoI <-ggplot(info.avg.PoI, aes(x = id, y = features.plot, 
                                       size = average, fill = average)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low="#e0ecf4", high="#762a83",limits=c(0,1.25))+
  scale_size(limits=c(0,1.25)) +
  scale_x_discrete(breaks=info.avg.PoI$id,
                   labels=round(info.avg.PoI$average,3))+
  theme_minimal() + labs(y="", x="") +
  theme(legend.position="none")


dim(data$data)


data<-VlnPlot(object = immune.combined, features = PoI[i],
              group.by ='CD34_CD90_EPCR_ADT_Phase_SIRT1',split.by='treatment')
visu_data_tmp <- data$data[-which(data$data$ident %like% "other"),]
visu_data_tmp2 <- visu_data_tmp[-which(visu_data_tmp$ident %like% c("CD90-")),]
visu_data_tmp3 <- visu_data_tmp2[-which(visu_data_tmp2$ident %like% c("low")),]
dim(visu_data_tmp)
#visu_data <- visu_data_tmp[which(visu_data_tmp$SIRT1 >0.01),]
visu_data <- visu_data_tmp3
dim(visu_data)
#scale_fill_brewer() + 
violin.PoI <-ggplot(visu_data, aes(x=ident, y=SIRT1,fill=split)) +
  geom_violin(trim=FALSE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  scale_fill_manual(values = c("#BCE4D8", "#c994c7","#49A4B9","#ce1256", "#2C5985")) +
  labs(y=paste0("Expression of ",PoI[i]), x="",fill ="Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) 
# theme(legend.position = "none")
#  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


figure2 <- multi_panel_figure(columns = 5, rows = 5, panel_label_type = "upper-roman")
figure2 %<>%
  fill_panel(violin.PoI, column = 1:5, row = 1:3) %<>%
  fill_panel(dot.PoI.num, column = 1:4, row = 4) %<>%
  fill_panel(dot.avg.PoI, column = 1:4, row = 5)
figure2


## SIRT3
i<-2
PoI[i]

###Define the range to work on
DefaultAssay(immune.combined) <- "RNA"
data_RNA_gene<- immune.combined_RNA_data_scale[rownames(immune.combined_RNA_data_scale) %in% PoI[i],]
(first_quantile <- quantile(data_RNA_gene, 0.25)) # first quartile
(third_quantile <- quantile(data_RNA_gene, 0.75)) # third quartile
(mean_gene <-  mean(data_RNA_gene))
(max_gene <-  max(data_RNA_gene))

SIRT3_RNA <- rep("",ncol(immune.combined_RNA_data_scale))
SIRT3_RNA[which(data_RNA_gene == threshold &
                  CD34_CD90_EPCR_ADT %in% "CD34+CD90+EPCR+" )] <- "SIRT3_0_low"
SIRT3_RNA[which(data_RNA_gene > threshold &
                  CD34_CD90_EPCR_ADT %in% "CD34+CD90+EPCR+" )] <- "SIRT3_2_high"

table(SIRT3_RNA)

## create cluster
CD34_CD90_EPCR_ADT_Phase_SIRT3 <- paste0(CD34_CD90_EPCR_ADT_Phase,"_",SIRT3_RNA)
immune.combined$CD34_CD90_EPCR_ADT_Phase_SIRT3 <- CD34_CD90_EPCR_ADT_Phase_SIRT3
CD34_CD90_EPCR_ADT_Phase_SIRT3_day_tmp <- paste0(CD34_CD90_EPCR_ADT_Phase_SIRT3,"_",immune.combined$treatment)
data_SIRT3 <- table(CD34_CD90_EPCR_ADT_Phase_SIRT3_day_tmp)
length(table(CD34_CD90_EPCR_ADT_Phase_SIRT3_day_tmp))

Idents(immune.combined) <- immune.combined$CD34_CD90_EPCR_ADT_Phase_SIRT3

###visualisation SIRT3
data_SIRT3_df_tmp <- as.data.frame(data_SIRT3)
data_SIRT3_df_tmp$features.plot <- "SIRT3"
colnames(data_SIRT3_df_tmp)[1]<-"id"

data_SIRT3_df_tmp2 <- data_SIRT3_df_tmp[-which(data_SIRT3_df_tmp$id %like% c("other")),]
data_SIRT3_df_tmp3 <- data_SIRT3_df_tmp2[-which(data_SIRT3_df_tmp2$id %like% c("CD90-")),]
data_SIRT3_df_tmp4 <- data_SIRT3_df_tmp3[-which(data_SIRT3_df_tmp3$id %like% c("low")),]

nodata <-data_SIRT3_df_tmp4[1,] 
nodata$id<-"CD34+CD90+EPCR+_3_G2M_INKA1_2_high_day0_none"
nodata$Freq <- 0
levels(data_SIRT3_df_tmp4$id) <-c(levels(data_SIRT3_df_tmp4$id),
                                  "CD34+CD90+EPCR+_3_G2M_INKA1_2_high_day0_none")
data_SIRT3_df1 <- rbind(data_SIRT3_df_tmp4,nodata)
data_SIRT3_df1$id <- as.character(data_SIRT3_df1$id)
data_SIRT3_df <-data_SIRT3_df1[order(data_SIRT3_df1$id),]
data_SIRT3_df$id <- as.factor(data_SIRT3_df$id)

dot.PoI.num <-ggplot(data_SIRT3_df, aes(x = id, y = features.plot,
                                        size = Freq, fill = Freq)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low="#e0ecf4", high="#762a83",limits=c(0,100))+
  scale_size(limits=c(0,5000)) +
  scale_x_discrete(breaks=data_SIRT3_df$id,
                   labels=data_SIRT3_df$Freq)+
  theme_minimal() + labs(y="", x="") +
  theme(legend.position="none")

info_dot<- DotPlot(immune.combined, features = PoI, 
                   cols=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"),
                   dot.scale = 8,  split.by = "treatment")

data_dot <- info_dot$data
dot_spec_tmp <- data_dot[-which(data_dot$id %like% c("other")),]
dot_spec_tmp2 <- dot_spec_tmp[-which(dot_spec_tmp$id %like% c("CD90-")),]
dot_spec_tmp3 <- dot_spec_tmp2[-which(dot_spec_tmp2$id %like% c("low")),]
dot_spec <- dot_spec_tmp3[which(dot_spec_tmp3$features.plot %in% PoI[i]),]

dot_spec$features.plot <- "SIRT3"
dot.PoI <-ggplot(dot_spec, aes(x = id, y = features.plot,
                               size = pct.exp, fill = pct.exp)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low="#e0ecf4", high="#762a83",limits=c(0,100))+
  scale_size(limits=c(0,100)) +
  scale_x_discrete(breaks=dot_spec$id,
                   labels=paste0(round(dot_spec$pct.exp,2),"%"))+
  theme_minimal() + labs(y="", x="") +
  theme(legend.position="none")


data<-VlnPlot(object = immune.combined, features = PoI[i],
              group.by ='treatment',split.by='CD34_CD90_EPCR_ADT_Phase_SIRT3')

info.avg.PoI_tmp<- data$data %>% group_by(split,ident) %>% dplyr::summarise(n = n(), 
                                                                            average = mean(SIRT3), 
                                                                            maximum = max(SIRT3))
info.avg.PoI_tmp$id <- paste0(info.avg.PoI_tmp$split,"_",info.avg.PoI_tmp$ident)
info.avg.PoI_tmp$features.plot <- "SIRT3"
info.avg.PoI_tmp2 <- info.avg.PoI_tmp[-which(info.avg.PoI_tmp$split %like% c("other_")),]
info.avg.PoI_tmp3 <- info.avg.PoI_tmp2[-which(info.avg.PoI_tmp2$split %like% c("CD90-")),]
info.avg.PoI_tmp4 <- info.avg.PoI_tmp3[-which(info.avg.PoI_tmp3$split %like% c("low")),]

nodata <-info.avg.PoI_tmp4[1,] 
nodata$id<-"CD34+CD90+EPCR+_3_G2M_INKA1_2_high_day0_none"
nodata$split<-"CD34+CD90+EPCR+_3_G2M_INKA1_2_high"
nodata$average <- 0
levels(info.avg.PoI_tmp4$id) <-c(levels(info.avg.PoI_tmp4$id),
                                 "CD34+CD90+EPCR+_3_G2M_INKA1_2_high_day0_none")
info.avg.PoI_df1 <- rbind(info.avg.PoI_tmp4,nodata)
info.avg.PoI_df1$id <- as.character(info.avg.PoI_df1$id)
info.avg.PoI <-info.avg.PoI_df1[order(info.avg.PoI_df1$id),]
info.avg.PoI$id <- as.factor(info.avg.PoI$id)


dot.avg.PoI <-ggplot(info.avg.PoI, aes(x = id, y = features.plot, 
                                       size = average, fill = average)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low="#e0ecf4", high="#762a83",limits=c(0,1.25))+
  scale_size(limits=c(0,1.25)) +
  scale_x_discrete(breaks=info.avg.PoI$id,
                   labels=round(info.avg.PoI$average,3))+
  theme_minimal() + labs(y="", x="") +
  theme(legend.position="none")


dim(data$data)


data<-VlnPlot(object = immune.combined, features = PoI[i],
              group.by ='CD34_CD90_EPCR_ADT_Phase_SIRT3',split.by='treatment')
visu_data_tmp <- data$data[-which(data$data$ident %like% "other"),]
visu_data_tmp2 <- visu_data_tmp[-which(visu_data_tmp$ident %like% c("CD90-")),]
visu_data_tmp3 <- visu_data_tmp2[-which(visu_data_tmp2$ident %like% c("low")),]
dim(visu_data_tmp3)

#create nul value
d4_vpa <-visu_data_tmp3[1,]
d4_vpa$ident<-"CD34+CD90+EPCR+_3_G2M_SIRT3_2_high"
d4_vpa$SIRT3 <-0
levels(visu_data_tmp3$ident) <-c(levels(visu_data_tmp3$ident),"CD34+CD90+EPCR+_3_G2M_SIRT3_2_high")
visu_data <- rbind(visu_data_tmp3,d4_vpa,d4_vpa)

dim(visu_data)
#scale_fill_brewer() + 
violin.PoI <-ggplot(visu_data, aes(x=ident, y=SIRT3,fill=split)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  scale_fill_manual(values = c("#BCE4D8", "#c994c7","#49A4B9","#ce1256", "#2C5985")) +
  labs(y=paste0("Expression of ",PoI[i]), x="",fill ="Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
#  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


figure2 <- multi_panel_figure(columns = 5, rows = 5, panel_label_type = "upper-roman")
figure2 %<>%
  fill_panel(violin.PoI, column = 1:5, row = 1:3) %<>%
  fill_panel(dot.PoI.num, column = 1:4, row = 4) %<>%
  fill_panel(dot.avg.PoI, column = 1:4, row = 5)
figure2


## INKA1
i<-3
PoI[i]

###Define the range to work on
DefaultAssay(immune.combined) <- "RNA"
data_RNA_gene<- immune.combined_RNA_data_scale[rownames(immune.combined_RNA_data_scale) %in% PoI[i],]
(first_quantile <- quantile(data_RNA_gene, 0.25)) # first quartile
(third_quantile <- quantile(data_RNA_gene, 0.75)) # third quartile
(mean_gene <-  mean(data_RNA_gene))
(max_gene <-  max(data_RNA_gene))

INKA1_RNA <- rep("",ncol(immune.combined_RNA_data_scale))
INKA1_RNA[which(data_RNA_gene ==  threshold &
                  CD34_CD90_EPCR_ADT %in% "CD34+CD90+EPCR+" )] <- "INKA1_0_low"
INKA1_RNA[which(data_RNA_gene > threshold &
                  CD34_CD90_EPCR_ADT %in% "CD34+CD90+EPCR+" )] <- "INKA1_2_high"

table(INKA1_RNA)

## create cluster
CD34_CD90_EPCR_ADT_Phase_INKA1 <- paste0(CD34_CD90_EPCR_ADT_Phase,"_",INKA1_RNA)
immune.combined$CD34_CD90_EPCR_ADT_Phase_INKA1 <- CD34_CD90_EPCR_ADT_Phase_INKA1
CD34_CD90_EPCR_ADT_Phase_INKA1_day_tmp <- paste0(CD34_CD90_EPCR_ADT_Phase_INKA1,"_",immune.combined$treatment)
data_INKA1 <- table(CD34_CD90_EPCR_ADT_Phase_INKA1_day_tmp)
length(table(CD34_CD90_EPCR_ADT_Phase_INKA1_day_tmp))

Idents(immune.combined) <- immune.combined$CD34_CD90_EPCR_ADT_Phase_INKA1

###visualisation INKA1
data_INKA1_df_tmp <- as.data.frame(data_INKA1)
data_INKA1_df_tmp$features.plot <- "INKA1"
colnames(data_INKA1_df_tmp)[1]<-"id"

data_INKA1_df_tmp2 <- data_INKA1_df_tmp[-which(data_INKA1_df_tmp$id %like% c("other")),]
data_INKA1_df_tmp3 <- data_INKA1_df_tmp2[-which(data_INKA1_df_tmp2$id %like% c("CD90-")),]
data_INKA1_df_tmp4 <- data_INKA1_df_tmp3[-which(data_INKA1_df_tmp3$id %like% c("low")),]

nodata <-data_INKA1_df_tmp4[1,] 
nodata$id<-"CD34+CD90+EPCR+_3_G2M_INKA1_2_high_day0_none"
nodata$Freq <- 0
levels(data_INKA1_df_tmp4$id) <-c(levels(data_INKA1_df_tmp4$id),
                                  "CD34+CD90+EPCR+_3_G2M_INKA1_2_high_day0_none")
data_INKA1_df1 <- rbind(data_INKA1_df_tmp4,nodata)
data_INKA1_df1$id <- as.character(data_INKA1_df1$id)
data_INKA1_df <-data_INKA1_df1[order(data_INKA1_df1$id),]
data_INKA1_df$id <- as.factor(data_INKA1_df$id)

dot.PoI.num <-ggplot(data_INKA1_df, aes(x = id, y = features.plot,
                                        size = Freq, fill = Freq)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low="#e0ecf4", high="#762a83",limits=c(0,100))+
  scale_size(limits=c(0,5000)) +
  scale_x_discrete(breaks=data_INKA1_df$id,
                   labels=data_INKA1_df$Freq)+
  theme_minimal() + labs(y="", x="") +
  theme(legend.position="none")

info_dot<- DotPlot(immune.combined, features = PoI, 
                   cols=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"),
                   dot.scale = 8,  split.by = "treatment")

data_dot <- info_dot$data
dot_spec_tmp <- data_dot[-which(data_dot$id %like% c("other")),]
dot_spec_tmp2 <- dot_spec_tmp[-which(dot_spec_tmp$id %like% c("CD90-")),]
dot_spec_tmp3 <- dot_spec_tmp2[-which(dot_spec_tmp2$id %like% c("low")),]
dot_spec <- dot_spec_tmp2[which(dot_spec_tmp2$features.plot %in% PoI[i]),]

dot_spec$features.plot <- "INKA1"
dot.PoI <-ggplot(dot_spec, aes(x = id, y = features.plot,
                               size = pct.exp, fill = pct.exp)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low="#e0ecf4", high="#762a83",limits=c(0,100))+
  scale_size(limits=c(0,100)) +
  scale_x_discrete(breaks=dot_spec$id,
                   labels=paste0(round(dot_spec$pct.exp,2),"%"))+
  theme_minimal() + labs(y="", x="") +
  theme(legend.position="none")


data<-VlnPlot(object = immune.combined, features = PoI[i],
              group.by ='treatment',split.by='CD34_CD90_EPCR_ADT_Phase_INKA1')

info.avg.PoI_tmp<- data$data %>% group_by(split,ident) %>% dplyr::summarise(n = n(), 
                                                                            average = mean(INKA1), 
                                                                            maximum = max(INKA1))
info.avg.PoI_tmp$id <- paste0(info.avg.PoI_tmp$split,"_",info.avg.PoI_tmp$ident)
info.avg.PoI_tmp$features.plot <- "INKA1"
info.avg.PoI_tmp2 <- info.avg.PoI_tmp[-which(info.avg.PoI_tmp$split %like% c("other_")),]
info.avg.PoI_tmp3 <- info.avg.PoI_tmp2[-which(info.avg.PoI_tmp2$split %like% c("CD90-")),]
info.avg.PoI_tmp4 <- info.avg.PoI_tmp3[-which(info.avg.PoI_tmp3$split %like% c("low")),]


nodata <-info.avg.PoI_tmp4[1,] 
nodata$id<-"CD34+CD90+EPCR+_3_G2M_INKA1_2_high_day0_none"
nodata$split<-"CD34+CD90+EPCR+_3_G2M_INKA1_2_high"
nodata$average <- 0
levels(info.avg.PoI_tmp4$id) <-c(levels(info.avg.PoI_tmp4$id),
                                  "CD34+CD90+EPCR+_3_G2M_INKA1_2_high_day0_none")
info.avg.PoI_df1 <- rbind(info.avg.PoI_tmp4,nodata)
info.avg.PoI_df1$id <- as.character(info.avg.PoI_df1$id)
info.avg.PoI <-info.avg.PoI_df1[order(info.avg.PoI_df1$id),]
info.avg.PoI$id <- as.factor(info.avg.PoI$id)

dot.avg.PoI <-ggplot(info.avg.PoI, aes(x = id, y = features.plot, 
                                       size = average, fill = average)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low="#e0ecf4", high="#762a83",limits=c(0,1.25))+
  scale_size(limits=c(0,1.25)) +
  scale_x_discrete(breaks=info.avg.PoI$id,
                   labels=round(info.avg.PoI$average,3))+
  theme_minimal() + labs(y="", x="") +
  theme(legend.position="none")


dim(data$data)


data<-VlnPlot(object = immune.combined, features = PoI[i],
              group.by ='CD34_CD90_EPCR_ADT_Phase_INKA1',split.by='treatment')
visu_data_tmp <- data$data[-which(data$data$ident %like% "other"),]
visu_data_tmp2 <- visu_data_tmp[-which(visu_data_tmp$ident %like% c("CD90-")),]
visu_data_tmp3 <- visu_data_tmp2[-which(visu_data_tmp2$ident %like% c("low")),]
dim(visu_data_tmp3)

#for day0/ S phase
dd <- visu_data_tmp3[which(visu_data_tmp3$ident %in% "CD34+CD90+EPCR+_2_S_INKA1_2_high" &
                             visu_data_tmp3$split %in% "day0_none"),]
dd1 <- dd
dd1[1] <- "0.3769456"

#create nul value
d4_vpa <-dd
d4_vpa$ident<-"CD34+CD90+EPCR+_3_G2M_INKA1_2_high"
d4_vpa$INKA1 <-0
levels(visu_data_tmp3$ident) <-c(levels(visu_data_tmp3$ident),"CD34+CD90+EPCR+_3_G2M_INKA1_2_high")
visu_data <- rbind(visu_data_tmp3,dd,d4_vpa,d4_vpa)
dim(visu_data)
#scale_fill_brewer() + 

violin.PoI <-ggplot(visu_data, aes(x=ident, y=INKA1,fill=split)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  scale_fill_manual(values = c("#BCE4D8", "#c994c7","#49A4B9","#ce1256", "#2C5985")) +
  labs(y=paste0("Expression of ",PoI[i]), x="",fill ="Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
#  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


figure2 <- multi_panel_figure(columns = 5, rows = 5, panel_label_type = "upper-roman")
figure2 %<>%
  fill_panel(violin.PoI, column = 1:5, row = 1:3) %<>%
  fill_panel(dot.PoI.num, column = 1:4, row = 4) %<>%
  fill_panel(dot.avg.PoI, column = 1:4, row = 5)
figure2



#### Heatmap GSEA
#https://datacarpentry.org/r-intro-geospatial/01-rstudio-intro/
#https://grunwaldlab.github.io/Population_Genetics_in_R/intro_to_R.html

#http://bioconductor.org/packages/devel/bioc/vignettes/escape/inst/doc/vignette.html
#https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html
#https://bioconductor.org/packages/release/bioc/html/GSEABase.html
#https://www.biostars.org/p/442531/
#https://www.nature.com/articles/s41467-020-15298-6

library(escape)
library(dittoSeq)
library(GSEABase)

DefaultAssay(immune.combined) <- "RNA"

CD34_CD90_EPCR_ADT_treatment_Phase <- paste0(CD34_CD90_EPCR_ADT,"_",
                                             immune.combined$treatment,"_",
                                             immune.combined$PhaseNum)
Idents(immune.combined) <-  CD34_CD90_EPCR_ADT_treatment_Phase

####LT_HSC
markers.to.plot_LT_HSC <- c("Ches1", "SREC", "Blr1", "Procr",
                     "Fzd4", "Igf1r", "Mtap7", "MYO5C", "Pclo", "Sparcl1", "Ocln",
                     "Jcam2", "Jcam3", "Mpdz", "Nbea", "SCOP", "Ptpn21", "Ndr2",
                     "Rras", "Agpt", "Efnb2", "Rbp1", "Aldh2", "Fkbp7", "Smpd1",
                     "Tapbp", "Nnp1", "Htf9c", "Elavl4", "Tcf3", "Pphn", "Hoxa5", "P2rx4",
                     "Slc12a2")
markers.to.plot_LT_HSC <- toupper(markers.to.plot_LT_HSC)

#DoHeatmap(immune.combined, features = markers.to.plot_LT_HSC) #+ NoLegend()

#LT_HST monkey
monkey_gene_LT_HTS <-c("ARHGEF12", "ITSN2", "MEIS1", "SVIL", "STARD9", "ZNF154",
                       "FAM30A", "ZBTB20", "MYEF2", "SCAI", "LRBA", "CAPN2", 
                       "PDZD2", "PRKCH", "ZBTB21", "HEG1", "C7orf49", 
                       "ABLIM1", "ENSG00000189089", "NFAT5", "CCPG1", "C1orf21", 
                       "ZFHX3", "MAML2", "ARMCX4", "MYCT1", "ZNF286A", "LPP", 
                       "ENSG00000281195", "CRIM1", "DSG2", "TMEM136", "SCN9A",
                       "ABCA5", "HOXB3", "SNTB1", "CLU", "GPRASP1", "CHN2",
                       "ENSG00000260244", "FZD6", "NR1D2", "TCEANC2", "PARP11",
                       "DST", "FHL1", "ZKSCAN1", "RBPMS", "NRIP1", "ZMYM1", 
                       "MLLT3", "PLAG1", "LOC100505501", "TFPI", "HMGA2", "SOCS2",
                       "BTB37", "SETBP1", "GUCY1A3", "RAB29", "MAPK11", "ZNF462",
                       "TMEM163", "PLK2", "EVA1C", "BMPR2", "DNAH6", "ALDH1A1",
                       "INPP4B", "TSPAN2", "TOX", "HOPX", "ARHGEF40", "C11orf63",
                       "IDS", "SAMD12", "STYK1", "FAT4", "MIPOL1", "SLC22A17",
                       "MAGI2-AS3", "DLK1", "PRDM16", "GATAD1", "MMRN1", "AR",
                       "AKAP6", "ARMCX2", "MECOM", "COL6A2", "ENSG00000215208",
                       "PARD3B", "PDGFD", "ZNF532", "ADGRG6", "WASF3", "GATA3", 
                       "ENSG00000263345", "HLA-DQA1", "ENSG00000259591", "ZNF165", 
                       "ANK3", "PCDHGB7", "SKIDA1", "GPRASP2", "ENSG00000242795",
                       "LIMCH1", "ZNF483", "HLF", "TCEAL2", "AFDN", "LOC101927577",
                       "ST8SIA1", "PREX2", "MEG3","MME", "CR2", "SEMA3C", "PTGER3",
                       "CD84", "PIK3AP1", "CYP2U1", "ABCA1")
#Cell cycle
#http://www.informatics.jax.org/go/term/GO:0007049
folder_genes <- "/Users/tiphainemartin/Documents/lab_Hoffman/paper2/data/"
list_genes <- read.table(file=paste0(folder_genes,"GO_term_summary_20220217_005621.txt"),
                         header=TRUE,sep="\t")
markers.to.plot_CC <- toupper(unique(list_genes$Symbol))
#DoHeatmap(immune.combined, features = markers.to.plot_CC)

#Kaufmann_CD112lo
list_genes <- read.table(file=paste0(folder_genes,"kauffmann_CD112low.txt"),
                         header=TRUE,sep="\t")
markers.to.plot_cd112lo <- toupper(unique(list_genes$gene))

#Kaufmann_CD112high
list_genes <- read.table(file=paste0(folder_genes,"kauffmann_CD112high.txt"),
                         header=TRUE,sep="\t")
markers.to.plot_cd112high <- toupper(unique(list_genes$gene))

gene.sets_perso <- list(LT_HSC = markers.to.plot_LT_HSC,
                       LT_HSC_Monkey = monkey_gene_LT_HTS,
                       Cell_Cycle = markers.to.plot_CC,
                       Kauffmann_CD112lo = markers.to.plot_cd112lo,
                       Kauffmann_CD112high = markers.to.plot_cd112high)

#Enrichment ssGSEA normalized
ES.seurat_norm <- enrichIt(obj = immune.combined, 
                      gene.sets = gene.sets_perso, 
                      method="ssGSEA",
                      ssGSEA.norm= TRUE,
                      groups = 1000, cores = 2,
                      min.size = 5)

## if working with a Seurat object
immune.combined <- Seurat::AddMetaData(immune.combined, ES.seurat_norm)

#visulization
colors <- colorRampPalette(c("#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF"))

Idents(immune.combined) <- immune.combined$treatment
p1 <- FeaturePlot(object = immune.combined, features = 'Kauffmann_CD112lo',
                  cols=c("#0D0887FF","#F0F921FF"))
p1

p1 <- FeaturePlot(object = immune.combined, features = 'Kauffmann_CD112high',
                  cols=c("#0D0887FF","#F0F921FF"))
p1

library(ComplexHeatmap)
library(circlize)
#https://github.com/jokergoo/ComplexHeatmap/issues/359#issuecomment-531669275
#https://www.biostars.org/p/398775/

df_heatmap_norm<-cbind(ES.seurat_norm)
annotation_norm_short <-paste0(CD34_CD90_EPCR_ADT,"_",
                               immune.combined$PhaseNum)
annotation_norm <- paste0(CD34_CD90_EPCR_ADT,"_",
                     immune.combined$PhaseNum,"_",
                     immune.combined$treatment)

write.table(df_heatmap_norm,file=paste0(folder_results,"df_heatmap_all_norm.csv"),
            sep=",",quote=FALSE,
            row.names=FALSE)
write.table(annotation_norm,file=paste0(folder_results,"annotation_all_norm.csv"),
            sep=",",quote=FALSE,
            row.names=FALSE)


no <- which(annotation_norm %like% "other" | annotation_norm %like% "CD90-")
df_heatmap_final_norm <- df_heatmap_norm[-no,]
annotation_final_norm <- annotation_norm[-no]
write.table(df_heatmap_final_norm,file=paste0(folder_results,"df_heatmap_final_norm.csv"),
            sep=",",quote=FALSE,
            row.names=FALSE)
write.table(annotation_final_norm,file=paste0(folder_results,"annotation_final_norm.csv"),
            sep=",",quote=FALSE,
            row.names=FALSE)

### Visuaisation only CD34+CD90+EPCR+
annotation_norm_short <-paste0(CD34_CD90_EPCR_ADT,"_",
                               immune.combined$PhaseNum)
no_short <- which(annotation_norm_short %like% "other" |
                    annotation_norm_short %like% "CD90-")
df_heatmap_final_norm_short <- df_heatmap_norm[-no_short,]
annotation_final_norm_short <- annotation_norm_short[-no_short]

# Define the number of colors you want
nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)
names(mycolors) <- names(table(annotation_final_norm))

ha = HeatmapAnnotation(Annotation = annotation_final_norm,
                       Annotation_short = annotation_final_norm_short,
                       col = list(Annotation = mycolors),
                       show_annotation_name = FALSE)

#draw(ha)
col_fun = colorRamp2(c(-0.2, 0, 0.25, 0.5, 1), c("#4393c3", "#f7f7f7", 
                                                        "#f4a582","#d6604d",
                                                        "#67001f"))
data <- t(df_heatmap_final_norm)
#dev.new(RStudioGD())
Heatmap(data, name = "ssGSEA",
        column_order = order(annotation_final_norm),
        show_column_names = FALSE, top_annotation = ha,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun)


#CD34+CD90+EPCR+ vs CD34+CD90-EPCR-
cell_type <- CD34_CD90_EPCR_ADT
phase_study <- immune.combined$PhaseNum
treatment <-  immune.combined$treatment

no_short <- which(cell_type %like% "other" )

df_heatmap_final_norm_all_short <- df_heatmap_norm[-no_short,]
cell_type_final <- cell_type[-no_short]
phase_study_final <- phase_study[-no_short]
treatment_final <- treatment[-no_short]
annotation_order<- paste0(cell_type_final,"_",
                          treatment_final,"_",
                          phase_study_final)

ha_all = HeatmapAnnotation("Cell type" = cell_type_final,
                       Phase = phase_study_final,
                       Treatment = treatment_final,
                       show_annotation_name = FALSE)
  
#draw(ha_all)
data_all <- t(df_heatmap_final_norm_all_short)
#dev.new(RStudioGD())
Heatmap(data_all, name = "ssGSEA",
        column_order = order(annotation_order),
        show_column_names = FALSE, top_annotation = ha_all,
        show_column_dend = FALSE, show_row_dend = FALSE,
        col = col_fun)

celltype_phase<- paste0(cell_type_final,"_",
                        phase_study_final)
data_all_annot <- cbind(t(data_all), cell_type_final,
                        phase_study_final, treatment_final,
                        annotation_order,celltype_phase)
df_data_all_annot <- as.data.frame(data_all_annot)
df_data_all_annot$Kauffmann_CD112lo <-as.numeric(as.character(df_data_all_annot$Kauffmann_CD112lo))
df_data_all_annot$Kauffmann_CD112high <-as.numeric(as.character(df_data_all_annot$Kauffmann_CD112high))


#visualisation of Kauffmann_CD112lo
## Per treatment
violin.treatment_cd112low<-ggplot(df_data_all_annot, aes(x=treatment_final, 
                                        y=Kauffmann_CD112lo,fill=treatment_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  scale_fill_manual(values = c("#BCE4D8", "#c994c7","#49A4B9","#ce1256", "#2C5985")) +
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.treatment_cd112low

## Per Cell type
one.way <- aov(Kauffmann_CD112lo ~ cell_type_final, data = df_data_all_annot)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]
violin.cellType_cd112low<-ggplot(df_data_all_annot, aes(x=cell_type_final, 
                                                         y=Kauffmann_CD112lo,
                                                         fill=cell_type_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  scale_fill_manual(values = c("darkgreen", "purple")) +
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Cell type") +
  theme_minimal() +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme(axis.text.x = element_blank())
violin.cellType_cd112low 

## Per phase_study_final
violin.annot_cd112low<-ggplot(df_data_all_annot, aes(x=phase_study_final, 
                                                     y=Kauffmann_CD112lo,
                                                     fill=phase_study_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Phase") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low


## Per Combination annotation
violin.annot_cd112low<-ggplot(df_data_all_annot, aes(x=annotation_order, 
                                                        y=Kauffmann_CD112lo,
                                                        fill=annotation_order)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Annotation") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

## Per Combination annotation
one.way <- aov(Kauffmann_CD112lo ~ celltype_phase, data = df_data_all_annot)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112low<-ggplot(df_data_all_annot, aes(x=celltype_phase, 
                                                     y=Kauffmann_CD112lo,
                                                     fill=celltype_phase)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

#visualisation of Kauffmann_CD112high
## Per treatment
violin.treatment_cd112high<-ggplot(df_data_all_annot, aes(x=treatment_final, 
                                                         y=Kauffmann_CD112high,
                                                         fill=treatment_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  scale_fill_manual(values = c("#BCE4D8", "#c994c7","#49A4B9","#ce1256", "#2C5985")) +
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.treatment_cd112high

## Per Cell type
one.way <- aov(Kauffmann_CD112high ~ cell_type_final, data = df_data_all_annot)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]
violin.cellType_cd112high<-ggplot(df_data_all_annot, aes(x=cell_type_final, 
                                                        y=Kauffmann_CD112high,
                                                        fill=cell_type_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  scale_fill_manual(values = c("darkgreen", "purple")) +
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Cell type") +
  theme_minimal() +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme(axis.text.x = element_blank())
violin.cellType_cd112high 

## Per phase_study_final
violin.annot_cd112high<-ggplot(df_data_all_annot, aes(x=phase_study_final, 
                                                     y=Kauffmann_CD112high,
                                                     fill=phase_study_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Phase") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112high


## Per Combination annotation
violin.annot_cd112high<-ggplot(df_data_all_annot, aes(x=annotation_order, 
                                                     y=Kauffmann_CD112high,
                                                     fill=annotation_order)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Annotation") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112high

## Per Combination annotation
one.way <- aov(Kauffmann_CD112high ~ celltype_phase, data = df_data_all_annot)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112high<-ggplot(df_data_all_annot, aes(x=celltype_phase, 
                                                     y=Kauffmann_CD112high,
                                                     fill=celltype_phase)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112high

#Relationship between Kauffmann CD112 low and high with SIRT1, SIRT3, INKA1
PoI <- c("SIRT1", "SIRT3", "INKA1")
immune.combined_RNA_data_scale<- GetAssayData(object = immune.combined, slot = "data")
threshold =0

#SIRT1
data_RNA_gene<- immune.combined_RNA_data_scale[rownames(immune.combined_RNA_data_scale) %in% PoI[1],]
SIRT1_RNA <- rep("",ncol(immune.combined_RNA_data_scale))
SIRT1_RNA[which(data_RNA_gene ==  threshold )] <- "SIRT1_low"
SIRT1_RNA[which(data_RNA_gene > threshold )] <- "SIRT1_high"
SIRT1_RNA_final <- SIRT1_RNA[-no_short]
SIRT1_RNA_value <- data_RNA_gene[-no_short]

#SIRT3
data_RNA_gene<- immune.combined_RNA_data_scale[rownames(immune.combined_RNA_data_scale) %in% PoI[2],]
SIRT3_RNA <- rep("",ncol(immune.combined_RNA_data_scale))
SIRT3_RNA[which(data_RNA_gene ==  threshold )] <- "SIRT3_low"
SIRT3_RNA[which(data_RNA_gene > threshold )] <- "SIRT3_high"
SIRT3_RNA_final <- SIRT3_RNA[-no_short]
SIRT3_RNA_value  <- data_RNA_gene[-no_short]

#INKA1
data_RNA_gene<- immune.combined_RNA_data_scale[rownames(immune.combined_RNA_data_scale) %in% PoI[3],]
INKA1_RNA <- rep("",ncol(immune.combined_RNA_data_scale))
INKA1_RNA[which(data_RNA_gene ==  threshold )] <- "INKA1_low"
INKA1_RNA[which(data_RNA_gene > threshold )] <- "INKA1_high"
INKA1_RNA_final <- INKA1_RNA[-no_short]
INKA1_RNA_value  <- data_RNA_gene[-no_short]

###Visualisation
df_data_all_annot_v2 <- cbind(df_data_all_annot,INKA1_RNA_final,SIRT3_RNA_final,
                              SIRT1_RNA_final,
                              INKA1_RNA_value,SIRT3_RNA_value,
                              SIRT1_RNA_value)

##SIRT1
#Kauffmann_CD112high
one.way <- aov(Kauffmann_CD112high ~ SIRT1_RNA_final, data = df_data_all_annot_v2)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112high<-ggplot(df_data_all_annot_v2, aes(x=SIRT1_RNA_final, 
                                                     y=Kauffmann_CD112high,
                                                     fill=SIRT1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112high

#Kauffmann_CD112lo
one.way <- aov(Kauffmann_CD112lo ~ SIRT1_RNA_final, data = df_data_all_annot_v2)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112low<-ggplot(df_data_all_annot_v2, aes(x=SIRT1_RNA_final, 
                                                     y=Kauffmann_CD112lo,
                                                     fill=SIRT1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

##SIRT3
#Kauffmann_CD112high
one.way <- aov(Kauffmann_CD112high ~ SIRT3_RNA_final, data = df_data_all_annot_v2)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112high<-ggplot(df_data_all_annot_v2, aes(x=SIRT3_RNA_final, 
                                                      y=Kauffmann_CD112high,
                                                      fill=SIRT3_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112high

#Kauffmann_CD112lo
one.way <- aov(Kauffmann_CD112lo ~ SIRT3_RNA_final, data = df_data_all_annot_v2)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112low<-ggplot(df_data_all_annot_v2, aes(x=SIRT3_RNA_final, 
                                                     y=Kauffmann_CD112lo,
                                                     fill=SIRT3_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

##INKA1
#Kauffmann_CD112high
one.way <- aov(Kauffmann_CD112high ~ INKA1_RNA_final, data = df_data_all_annot_v2)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112high<-ggplot(df_data_all_annot_v2, aes(x=INKA1_RNA_final, 
                                                      y=Kauffmann_CD112high,
                                                      fill=INKA1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112high

#Kauffmann_CD112lo
one.way <- aov(Kauffmann_CD112lo ~ INKA1_RNA_final, data = df_data_all_annot_v2)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112low<-ggplot(df_data_all_annot_v2, aes(x=INKA1_RNA_final, 
                                                     y=Kauffmann_CD112lo,
                                                     fill=INKA1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

#INKA1 and celltype
INKA1_Celltype <- paste0(df_data_all_annot_v2$cell_type_final,
                         df_data_all_annot_v2$INKA1_RNA_final)
df_data_all_annot_v2 <-cbind(df_data_all_annot_v2,INKA1_Celltype)

one.way <- aov(Kauffmann_CD112lo ~ INKA1_Celltype, data = df_data_all_annot_v2)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]
violin.annot_cd112low<-ggplot(df_data_all_annot_v2, aes(x=INKA1_Celltype, 
                                                     y=Kauffmann_CD112lo,
                                                     fill=INKA1_Celltype)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

###Only CD34+CD90+EPCR+
onlypos <- which(cell_type_final %like% "CD90-" )
df_data_all_annot_v3 <- df_data_all_annot_v2[-onlypos,]

##SIRT1
#Kauffmann_CD112high
one.way <- aov(Kauffmann_CD112high ~ SIRT1_RNA_final, data = df_data_all_annot_v3)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112high<-ggplot(df_data_all_annot_v3, aes(x=SIRT1_RNA_final, 
                                                         y=Kauffmann_CD112high,
                                                         fill=SIRT1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112high

#Kauffmann_CD112lo
one.way <- aov(Kauffmann_CD112lo ~ SIRT1_RNA_final, data = df_data_all_annot_v3)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112low<-ggplot(df_data_all_annot_v3, aes(x=SIRT1_RNA_final, 
                                                        y=Kauffmann_CD112lo,
                                                        fill=SIRT1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

##SIRT3
#Kauffmann_CD112high
one.way <- aov(Kauffmann_CD112high ~ SIRT3_RNA_final, data = df_data_all_annot_v3)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112high<-ggplot(df_data_all_annot_v3, aes(x=SIRT3_RNA_final, 
                                                         y=Kauffmann_CD112high,
                                                         fill=SIRT3_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112high

#Kauffmann_CD112lo
one.way <- aov(Kauffmann_CD112lo ~ SIRT3_RNA_final, data = df_data_all_annot_v3)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112low<-ggplot(df_data_all_annot_v3, aes(x=SIRT3_RNA_final, 
                                                        y=Kauffmann_CD112lo,
                                                        fill=SIRT3_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

##INKA1
#Kauffmann_CD112high
one.way <- aov(Kauffmann_CD112high ~ INKA1_RNA_final, data = df_data_all_annot_v3)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112high<-ggplot(df_data_all_annot_v3, aes(x=INKA1_RNA_final, 
                                                         y=Kauffmann_CD112high,
                                                         fill=INKA1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 high"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112high

#Kauffmann_CD112lo
one.way <- aov(Kauffmann_CD112lo ~ INKA1_RNA_final, data = df_data_all_annot_v3)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112low<-ggplot(df_data_all_annot_v3, aes(x=INKA1_RNA_final, 
                                                        y=Kauffmann_CD112lo,
                                                        fill=INKA1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of Kauffmann CD112 low"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

##INKA1 vs SIRT1 
one.way <- aov(SIRT1_RNA_value ~ INKA1_RNA_final, data = df_data_all_annot_v3)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112low<-ggplot(df_data_all_annot_v3, aes(x=INKA1_RNA_final, 
                                                        y=SIRT1_RNA_value,
                                                        fill=INKA1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of SIRT1"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

df_data_all_annot_v3 %>%group_by(INKA1_RNA_final)  %>% summarize(n=n(), avg=mean(SIRT1_RNA_value),sd=sd(SIRT1_RNA_value))

df_data_all_annot_v3 %>%group_by(INKA1_RNA_final,SIRT1_RNA_final)  %>% summarize(n=n(), avg=mean(SIRT1_RNA_value),sd=sd(SIRT1_RNA_value))

##INKA1 vs SIRT3
one.way <- aov(SIRT3_RNA_value ~ INKA1_RNA_final, data = df_data_all_annot_v3)
info_one.way <- summary(one.way)
pvalue <-info_one.way[[1]][["Pr(>F)"]][1]

violin.annot_cd112low<-ggplot(df_data_all_annot_v3, aes(x=INKA1_RNA_final, 
                                                        y=SIRT3_RNA_value,
                                                        fill=INKA1_RNA_final)) +
  geom_violin(trim=TRUE)+ 
  geom_boxplot(width=0.07,position=position_dodge(width=0.9)) + 
  labs(y=paste0("Expression of SIRT3"), x="",fill ="Annotation") +
  ggtitle(paste0("p-value = ", pvalue)) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
violin.annot_cd112low

df_data_all_annot_v3 %>%group_by(INKA1_RNA_final)  %>% summarize(n=n(), avg=mean(SIRT3_RNA_value),sd=sd(SIRT3_RNA_value))

df_data_all_annot_v3 %>%group_by(INKA1_RNA_final,SIRT3_RNA_final)  %>% summarize(n=n(), avg=mean(SIRT3_RNA_value),sd=sd(SIRT3_RNA_value))


### Pseudotime
#http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html
#https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/velocity.html
#https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html
#https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
#https://www.datanovia.com/en/blog/easy-way-to-expand-color-palettes-in-r/

library("slingshot")
library("ggbeeswarm")
library("ggthemes")
library("RColorBrewer")

# Run the standard workflow for visualization and clustering
#Step done in the beginning, but need before running RunVelocity
#immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
#immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
#immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
#immune.combined <- FindClusters(immune.combined, resolution = 0.5)

sce_immune <- as.SingleCellExperiment(immune.combined)
immune.combined_pseudotime <- slingshot(sce_immune, reducedDim = 'PCA')
#No cluster labels provided. Continuing with one cluste

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(7, alpha = 1)
plot(reducedDims(immune.combined_pseudotime)$PCA, 
     col = colors[cut(immune.combined_pseudotime$slingPseudotime_1,breaks=7)],
     pch=16, asp = 1)
lines(SlingshotDataSet(immune.combined_pseudotime), lwd=2)


# Plot Slingshot pseudotime vs cell stage. 
sce_immune$numPhase <- sce_immune$Phase
sce_immune$numPhase[which(sce_immune$Phase %in% "G1")] <- "1_G0/G1"
sce_immune$numPhase[which(sce_immune$Phase %in% "S")] <- "2_S"
sce_immune$numPhase[which(sce_immune$Phase %in% "G2M")] <- "3_G2M"

sce_immune$phase_treatment <- paste0(sce_immune$numPhase,"_",sce_immune$treatment)

nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(as.data.frame(colData(sce_immune)),
       aes(x = immune.combined_pseudotime$slingPseudotime_1, y = phase_treatment, 
                                             colour = phase_treatment)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_fill_manual(values = mycolors)+ theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Phase/treatment timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

# Cluster cells using the Seurat workflow below.
gcdata <- CreateSeuratObject(counts = counts(immune.combined_pseudotime),
                             project = "slingshot")

gcdata <- NormalizeData(gcdata, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
gcdata <- FindVariableFeatures(gcdata, selection.method = "vst", 
                               nfeatures = 2000)
gcdata <- ScaleData(object = gcdata, do.center = T, do.scale = F)

gcdata <- RunPCA(gcdata, features = VariableFeatures(gcdata), npcs = 10, 
                 ndims.print = 1:5, nfeatures.print = 5)

# Cluster the cells using the first ten principal components.
gcdata <- FindNeighbors(gcdata, reduction = "pca", dims = 1:10, k.param = 10)

gcdata <- FindClusters(gcdata, resolution = 0.6, algorithm = 1, random.seed = 100)

# Add clustering information from Seurat to the immune.combined_pseudotime object
immune.combined_pseudotime$slingPseudotime_1 <- NULL  # remove old slingshot pseudotime data
colData(immune.combined_pseudotime)$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character
head(colData(immune.combined_pseudotime))

# Then run Slingshot using these new cluster assignments.
immune.combined_pseudotime <- slingshot(immune.combined_pseudotime, 
                                        clusterLabels = 'Seurat_clusters', 
                                        reducedDim = 'PCA')

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(7, alpha = 1)
mycolors_pseudo <- colors[cut(immune.combined_pseudotime$slingPseudotime_1,breaks=7)]
rD <- reducedDims(immune.combined_pseudotime)$PCA

plot(rD, 
     col = mycolors_pseudo, 
     pch=16, asp = 1,
     xlab="PC 1",
     ylab="PC 2",
     xlim=c(-10,25),
     ylim=c(-25,15))
lines(SlingshotDataSet(immune.combined_pseudotime), lwd=2)
#lines(SlingshotDataSet(immune.combined_pseudotime),type = 'lineages', 
#      lwd=1,col = 'grey')

#Per treatment
treatment <- immune.combined_pseudotime$treatment
table(treatment)

#day0
day0_cells <- which(treatment %in% "day0_none")
plot(rD[day0_cells,], 
     col = mycolors_pseudo[day0_cells], 
     pch=16, asp = 1, 
     xlab="PC 1",
     ylab="PC 2",
     xlim=c(-10,25),
     ylim=c(-25,15))
     title(main="Day 0")
lines(SlingshotDataSet(immune.combined_pseudotime), lwd=2)

#day1_ctrl
day1_ctrl_cells <- which(treatment %in% "day1_ctrl")
plot(rD[day1_ctrl_cells,], 
     col = mycolors_pseudo[day1_ctrl_cells], 
     pch=16, asp = 1, 
     xlab="PC 1",
     ylab="PC 2",
     xlim=c(-10,25),
     ylim=c(-25,15))
title(main="Day 1 ctrl")
lines(SlingshotDataSet(immune.combined_pseudotime), lwd=2)

#day1_vpa
day1_vpa_cells <- which(treatment %in% "day1_vpa")
plot(rD[day1_vpa_cells,], 
     col = mycolors_pseudo[day1_vpa_cells], 
     pch=16, asp = 1, 
     xlab="PC 1",
     ylab="PC 2",
     xlim=c(-10,25),
     ylim=c(-25,15))
title(main="Day 1 vpa")
lines(SlingshotDataSet(immune.combined_pseudotime), lwd=2)

#day4_ctrl
day4_ctrl_cells <- which(treatment %in% "day4_ctrl")
plot(rD[day4_ctrl_cells,], 
     col = mycolors_pseudo[day4_ctrl_cells], 
     pch=16, asp = 1, 
     xlab="PC 1",
     ylab="PC 2",
     xlim=c(-10,25),
     ylim=c(-25,15))
title(main="Day 4 ctrl")
lines(SlingshotDataSet(immune.combined_pseudotime), lwd=2)

#day4_vpa
day4_vpa_cells <- which(treatment %in% "day4_vpa")
plot(rD[day4_vpa_cells,], 
     col = mycolors_pseudo[day4_vpa_cells], 
     pch=16, asp = 1, 
     xlab="PC 1",
     ylab="PC 2",
     xlim=c(-10,25),
     ylim=c(-25,15))
title(main="Day 4 vpa")
lines(SlingshotDataSet(immune.combined_pseudotime), lwd=2)

#color pseudotime cluster defined by Slingshot
# df <- as.data.frame(rD)
# df$mycolors_pseudo <- mycolors_pseudo
# ggplot(df, aes(x=PC_1, y=PC_2,colour = mycolors_pseudo)) +
#   geom_point() +
#   scale_fill_manual(values = mycolors_pseudo)+ theme_classic() +
#   theme(legend.position="none") +
#   xlab("Slingshot pseudotime 1") + ylab("Slingshot pseudotime 2") +
#   facet_grid(. ~ treatment)
# 
# 
# #Color samples based on Phase and treatment
# df$treatment <- immune.combined_pseudotime$treatment
# df$Phase <- immune.combined_pseudotime$Phase
# df$numPhase <- immune.combined_pseudotime$Phase
# df$numPhase[which(immune.combined_pseudotime$Phase %in% "G1")] <- "1_G0/G1"
# df$numPhase[which(immune.combined_pseudotime$Phase %in% "S")] <- "2_S"
# df$numPhase[which(immune.combined_pseudotime$Phase %in% "G2M")] <- "3_G2M"
# 
# df$Phase_treatement <- paste0(df$numPhase,"_",
#                               df$treatment)
# ggplot(df, aes(x=PC_1, y=PC_2,colour = Phase)) +
#   geom_point() + theme_classic() +
#   xlab("Slingshot pseudotime 1") + ylab("Slingshot pseudotime 2") +
# facet_grid(. ~ treatment)


# Plot Slingshot pseudotime vs treatment and cell cycle. 
ggplot(as.data.frame(colData(sce_immune)),
       aes(x = immune.combined_pseudotime$slingPseudotime_1, y = phase_treatment, 
           colour = phase_treatment)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_fill_manual(values = mycolors)+ theme_classic() +
  xlab("Slingshot pseudotime 1") + ylab("Phase/treatment timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime 1")

ggplot(as.data.frame(colData(sce_immune)),
       aes(x = immune.combined_pseudotime$slingPseudotime_2, y = phase_treatment, 
           colour = phase_treatment)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_fill_manual(values = mycolors)+ theme_classic() +
  xlab("Slingshot pseudotime 2") + ylab("Phase/treatment timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime 2")

ggplot(as.data.frame(colData(sce_immune)),
       aes(x = immune.combined_pseudotime$slingPseudotime_3, y = phase_treatment, 
           colour = phase_treatment)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_fill_manual(values = mycolors)+ theme_classic() +
  xlab("Slingshot pseudotime 3") + ylab("Phase/treatment timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime ")


##Plots of gene expression over time.
#BiocManager::install("scater")
library("scater")
#INKA1
plotExpression(immune.combined_pseudotime, "INKA1", x = "slingPseudotime_1", 
               colour_by = "treatment", show_violin = FALSE,
               show_smooth = TRUE)

#SIRT1
plotExpression(immune.combined_pseudotime, "SIRT1", x = "slingPseudotime_1", 
               colour_by = "treatment", show_violin = FALSE,
               show_smooth = TRUE)

#SIRT3
plotExpression(immune.combined_pseudotime, "SIRT3", x = "slingPseudotime_1", 
               colour_by = "treatment", show_violin = FALSE,
               show_smooth = TRUE)

#INKA1 vs SIRT1
SIRT1

#INKA1 vs SIRT3
plotExpression(immune.combined_pseudotime, "INKA1", x = "SIRT3", 
               colour_by = "treatment", show_violin = FALSE,
               show_smooth = TRUE)

#SIRT1 vs SIRT3
plotExpression(immune.combined_pseudotime, "SIRT1", x = "SIRT3", 
               colour_by = "treatment", show_violin = FALSE,
               show_smooth = TRUE)
