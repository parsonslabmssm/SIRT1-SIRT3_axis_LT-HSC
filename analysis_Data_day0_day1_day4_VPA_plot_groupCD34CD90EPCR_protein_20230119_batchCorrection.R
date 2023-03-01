##https://satijalab.org/seurat/archive/v3.1/immune_alignment.html
##https://satijalab.org/seurat/articles/multimodal_vignette.html
#https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html (Annotation cluster immune system)
#https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
#https://satijalab.org/seurat/articles/integration_introduction.html
#https://www.singlecellcourse.org/introduction-to-single-cell-rna-seq.html
#https://ludvigla.github.io/STUtility_web_site/Normalization.html
#https://www.biostars.org/p/492222/
#https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/seurat-pre-process-filtering-confounding-genes.html
#https://learn.gencore.bio.nyu.edu/single-cell-rnaseq/seurat-part-3-data-normalization/
#https://chipster.csc.fi/manual/single-cell-seurat-filter-regress.html
#https://notebook.community/greenelab/GCB535/36_scRNAseq-I/scRNAseq_inclass_1

##Normalisation without anchoir
#https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
#https://satijalab.org/seurat/articles/sctransform_vignette.html 
#https://swaruplab.bio.uci.edu/tutorial/integration/integration_tutorial.html
#https://scalex.readthedocs.io/en/latest/tutorial/Integration_PBMC.html
#https://www.nature.com/articles/s41467-022-33758-z
#https://ucdavis-bioinformatics-training.github.io/2022-July-Single-Cell-RNA-Seq-Analysis/data_analysis/scRNA_Workshop-PART3
#https://rockefelleruniversity.github.io/scRNA-seq/exercises/answers/exercise1_answers.html

#https://r-charts.com/distribution/violin-plot-group-ggplot2/

#https://www.r-bloggers.com/2013/09/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/

##Install update SEURAT
# Install sctransform with V2
# install glmGamPoi
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")
# install sctransform from Github (automatically updates sctransform)
devtools::install_github("satijalab/sctransform", ref = "develop")

set.seed(123456)

#install.packages('Seurat')
library("Seurat")
library("Matrix")
library("ggplot2")
library("dplyr")
library("viridis")
library("RColorBrewer")
library("data.table")
library("multipanelfigure")
library("glmGamPoi")
library("patchwork")

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
seurat_object_day1_vpa = CreateSeuratObject(counts = data_rnaseq)

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

###batch associated with ADT sequencing
seurat_object_day0[["batch.adt"]] <- "TD005137_RonaldHoffman"
seurat_object_day1_vpa[["batch.adt"]] <- "TD005137_RonaldHoffman"
seurat_object_day1_ctrl[["batch.adt"]] <- "TD005137_RonaldHoffman"
seurat_object_day4_vpa[["batch.adt"]] <- "TD005139_ChristophSchaniel"
seurat_object_day4_ctrl[["batch.adt"]] <- "TD005139_ChristophSchaniel"

##number of genes with mitochondrial genes
seurat_object_day0[["percent.mt"]] <- PercentageFeatureSet(seurat_object_day0, pattern = "^MT-")
seurat_object_day1_vpa[["percent.mt"]] <- PercentageFeatureSet(seurat_object_day1_vpa, pattern = "^MT-")
seurat_object_day4_vpa[["percent.mt"]] <- PercentageFeatureSet(seurat_object_day4_vpa, pattern = "^MT-")
seurat_object_day1_ctrl[["percent.mt"]] <- PercentageFeatureSet(seurat_object_day1_ctrl, pattern = "^MT-")
seurat_object_day4_ctrl[["percent.mt"]] <- PercentageFeatureSet(seurat_object_day4_ctrl, pattern = "^MT-")

##number of genes with ribosomal proteins
seurat_object_day0[["percent.rb"]] <- PercentageFeatureSet(seurat_object_day0, pattern = "^RP[SL]")
seurat_object_day1_vpa[["percent.rb"]] <- PercentageFeatureSet(seurat_object_day1_vpa, pattern = "^RP[SL]")
seurat_object_day4_vpa[["percent.rb"]] <- PercentageFeatureSet(seurat_object_day4_vpa, pattern = "^RP[SL]")
seurat_object_day1_ctrl[["percent.rb"]] <- PercentageFeatureSet(seurat_object_day1_ctrl, pattern = "^RP[SL]")
seurat_object_day4_ctrl[["percent.rb"]] <- PercentageFeatureSet(seurat_object_day4_ctrl, pattern = "^RP[SL]")

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

###Visualisation before integration of samples
experiment.test <- immune.combined
DefaultAssay(experiment.test) <- "RNA"
experiment.test <- NormalizeData(object=experiment.test, assay="RNA")
experiment.test <- ScaleData(object=experiment.test, assay="RNA")
experiment.test <- FindVariableFeatures(object=experiment.test, assay="RNA")
experiment.test <- RunPCA(object=experiment.test, assay="RNA")
experiment.test <- RunUMAP(experiment.test, 
                           reduction = "pca", dims = 1:30)
experiment.test <- FindClusters(experiment.test, resolution = 0.5)

DimPlot(object = experiment.test, group.by="treatment", 
        reduction="pca", shuffle=TRUE)
DimPlot(object = experiment.test, group.by="treatment", 
        reduction="umap", shuffle=TRUE)
Idents(experiment.test) <- experiment.test$treatment
DimPlot(experiment.test, reduction = "umap", split.by = "treatment")
DimPlot(experiment.test, reduction = "pca", split.by = "treatment")

##Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- NormalizeData(object=immune.combined)
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- FindVariableFeatures(object=immune.combined)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", 
                                 dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)


DimPlot(object = immune.combined, group.by="treatment", 
        reduction="pca", shuffle=TRUE)
DimPlot(object = immune.combined, group.by="treatment", 
        reduction="umap", shuffle=TRUE)
DimPlot(immune.combined, reduction = "umap", split.by = "treatment")

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

#### Supplement Fig 1:  CD90
#before integration
DefaultAssay(experiment.test) <- "ADT"
experiment.test <- NormalizeData(experiment.test, 
                                 normalization.method = "CLR", 
                                 margin = 2)
p1_before <- FeaturePlot(experiment.test, "CD90-Human-GCATTGTACGATTCA",
                         split.by = "treatment",
                         cols = c("#e0e0e0", "#2d004b")) 

FeaturePlot(experiment.test, "CD90-Human-GCATTGTACGATTCA",
            cols = c("#e0e0e0", "#2d004b"))

#after integration
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD90-Human-GCATTGTACGATTCA",
                            split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b")) 
p1_treatment
FeaturePlot(immune.combined, "CD90-Human-GCATTGTACGATTCA",
            cols = c("#e0e0e0", "#2d004b"))


#### Supplement Fig 1: CD201
DefaultAssay(immune.combined) <- "ADT"
p1_treatment <- FeaturePlot(immune.combined, "CD201-Human-GTTTCCTTGACCAAG", split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b"))
p1_treatment

p1 <- FeaturePlot(immune.combined, "CD201-Human-GTTTCCTTGACCAAG", 
                  cols = c("#e0e0e0", "#2d004b")) +
  ggtitle("CD201 protein")

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


##CD34_CD90_EPCR hypothese: all cells are CD34+ RNA level
##CD34_CD90_ hypothese: all cells are CD34+
CD34_CD90_EPCR_ADT <- rep("other",ncol(immune.combined_ADT_data))
CD34_CD90_EPCR_ADT [which( CD90_ADT > 0 & CD201_ADT > 0 )] <- "CD34+CD90+EPCR+"
CD34_CD90_EPCR_ADT [which( CD90_ADT == 0 & CD201_ADT == 0)] <- "CD34+CD90-EPCR-"
table(CD34_CD90_EPCR_ADT )


#Fig 1: Visualisation of group CD34+CD90+EPCR+/CD34+CD90-EPCR-/other per treatment
# at protein level
#Before integration
experiment.test$group1 <- CD34_CD90_EPCR_ADT
df_immune_before <- as.data.frame(cbind(experiment.test$treatment,
                                        experiment.test$group1))
colnames(df_immune_before) <- c("treatment","group1")
table(df_immune_before$treatment)

#split.by = 'treatment',
p1 <- DimPlot(experiment.test, reduction = "umap", 
              group.by = "group1",
              cols = c("darkgreen","purple","#e0e0e0"))
p1

#After integration
immune.combined$group1 <- CD34_CD90_EPCR_ADT
df_immune <- as.data.frame(cbind(immune.combined$treatment,
                                 immune.combined$group1))
colnames(df_immune) <- c("treatment","group1")
table(df_immune$treatment)

p1 <- DimPlot(immune.combined, reduction = "umap", 
              group.by = "group1",
              cols = c("darkgreen","purple","#e0e0e0"))
p1


#### Integratin single cell datasets without anchor
#https://swaruplab.bio.uci.edu/tutorial/integration/integration_tutorial.html

# merge seurat object
seurat_list <- list("day0_none"=seurat_object_day0,
                    "day1_vpa"=seurat_object_day1_vpa,
                    "day4_vpa"=seurat_object_day4_vpa,
                    "day1_ctrl"=seurat_object_day1_ctrl,
                    "day4_ctrl"=seurat_object_day4_ctrl)

# normalize and identify variable features for each dataset independently
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#Merge data
seurat_obj_merge <- merge(x=seurat_list[[1]], 
                          y=seurat_list[2:length(seurat_list)])

#Check quality metrics and remove low-quality cells
VlnPlot(seurat_obj_merge, group.by="treatment", 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"), 
        ncol = 2, pt.size=0, )

seurat_obj_merge <- subset(seurat_obj_merge, nFeature_RNA > 200 & 
                       nFeature_RNA < 10000)

#Check after correction
VlnPlot(seurat_obj_merge, group.by="treatment", 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"), 
        ncol = 2, pt.size=0, )

#SCTransform
#Not enough memory
#library(sctransform)

#seurat_obj_merge <- SCTransform(
#  seurat_obj_merge,
#  vars.to.regress = c('batch.adt'),
#  verbose=TRUE
#)
#gc()

#saveRDS(seurat_obj_merge, file=paste0(out_data_dir, 'seurat_object_sct.rds'))

##Harmony
library(harmony)
DefaultAssay(seurat_obj_merge) <- "RNA"
seurat_obj_merge <- NormalizeData(object=seurat_obj_merge, assay="RNA")
seurat_obj_merge <- ScaleData(object=seurat_obj_merge, assay="RNA")
seurat_obj_merge <- FindVariableFeatures(object=seurat_obj_merge, assay="RNA")
seurat_obj_merge <- RunPCA(object=seurat_obj_merge, assay="RNA")
seurat_obj_merge <- RunHarmony(seurat_obj_merge, "batch.adt",assay.use = "sct")
seurat_obj_merge <- RunUMAP(seurat_obj_merge, reduction = "harmony",dims = 1:30)
seurat_obj_merge <- FindNeighbors(seurat_obj_merge, reduction='harmony')
seurat_obj_merge <- FindClusters(seurat_obj_merge, resolution = 0.3)

#visualisation
p1 <- DimPlot(seurat_obj_merge, group.by='seurat_clusters', reduction='umap') +
  ggtitle('seurat clusters') 
p1

p1_label <- DimPlot(seurat_obj_merge, group.by='seurat_clusters', reduction='umap',
              label=TRUE) +
  ggtitle('seurat clusters') 
p1_label

p2 <- DimPlot(seurat_obj_merge, group.by='treatment', reduction='umap') +
  ggtitle('Treatment') 
p2

### visualisation of different protein of interest
# Normalize ADT data,
DefaultAssay(seurat_obj_merge) <- "ADT"
seurat_obj_merge <- NormalizeData(seurat_obj_merge,
                                 normalization.method = "CLR", margin = 2)

# Now, we will visualize CD34, CD90, CD201 levels for RNA and protein level
#By setting the default assay, we can visualize one or the other

#### Supplement Fig 1:  CD90
DefaultAssay(seurat_obj_merge) <- "ADT"
p1_treatment <- FeaturePlot(seurat_obj_merge, "CD90-Human-GCATTGTACGATTCA",
                            split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b")) + 
                xlim(c(-10,10)) + ylim(c(-10,10))

FeaturePlot(seurat_obj_merge, "CD90-Human-GCATTGTACGATTCA",
            cols = c("#e0e0e0", "#2d004b"))

DefaultAssay(immune.combined) <- "RNA"
p2 <- FeaturePlot(seurat_obj_merge, "THY1",split.by = "treatment",
                  cols = c("lightgrey", "#542788")) + ggtitle("CD90 RNA") + 
  xlim(c(-10,10)) + ylim(c(-10,10))
p2

#### Supplement Fig 1: CD201
DefaultAssay(seurat_obj_merge) <- "ADT"
p1_treatment <- FeaturePlot(seurat_obj_merge, "CD201-Human-GTTTCCTTGACCAAG", 
                            split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b"))  + 
  xlim(c(-10,10)) + ylim(c(-10,10))
p1_treatment

p1 <- FeaturePlot(seurat_obj_merge, "CD201-Human-GTTTCCTTGACCAAG", 
                  cols = c("#e0e0e0", "#2d004b")) +
  ggtitle("CD201 protein")

DefaultAssay(seurat_obj_merge) <- "RNA"
p2 <- FeaturePlot(seurat_obj_merge, "PROCR", split.by = "treatment",
                  cols = c("#e0e0e0", "#2d004b")) + ggtitle("CD201 RNA")  + 
  xlim(c(-10,10)) + ylim(c(-10,10))
p2


#### Supplement Fig 1: CD49F
DefaultAssay(seurat_obj_merge) <- "ADT"
p1_treatment <- FeaturePlot(seurat_obj_merge, "CD49F-Human-TTCCGAGGATGATCT", split.by = "treatment",
                            cols = c("#e0e0e0", "#2d004b"))+ 
  xlim(c(-10,10)) + ylim(c(-10,10))
p1_treatment

p1 <- FeaturePlot(seurat_obj_merge, "CD49F-Human-TTCCGAGGATGATCT", 
                  cols = c("#e0e0e0", "#2d004b")) +
  ggtitle("CD201 protein")

DefaultAssay(seurat_obj_merge) <- "RNA"
p2 <- FeaturePlot(seurat_obj_merge, "ITGA6") + ggtitle("CD201 RNA")

p2_treatment <- FeaturePlot(seurat_obj_merge, "ITGA6", split.by = "treatment",
                            cols = c("lightgrey", "#542788")) + 
  xlim(c(-10,10)) + ylim(c(-10,10))
p2_treatment

#### Supplement Fig 1: CD49C
DefaultAssay(seurat_obj_merge) <- "RNA"
p2 <- FeaturePlot(seurat_obj_merge, "ITGA3") + ggtitle("CD49c RNA")

p2_treatment <- FeaturePlot(seurat_obj_merge, "ITGA3", split.by = "treatment",
                            cols = c("lightgrey", "#542788")) + 
  xlim(c(-10,10)) + ylim(c(-10,10))
p2_treatment

######## Annotation Cells!
DefaultAssay(seurat_obj_merge) <- "ADT"
seurat_obj_merge_ADT_data<- GetAssayData(object = seurat_obj_merge, slot = "counts")
CD34_ADT <- rep(0,ncol(seurat_obj_merge_ADT_data))
CD34_ADT[which(seurat_obj_merge_ADT_data[rownames(seurat_obj_merge_ADT_data) %in% "CD34-Human-GCAGAAATCTCCCTT",] > 0)] <- 1
CD90_ADT <- rep(0,ncol(seurat_obj_merge_ADT_data))
CD90_ADT[which(seurat_obj_merge_ADT_data[rownames(seurat_obj_merge_ADT_data) %in% "CD90-Human-GCATTGTACGATTCA",] > 0)] <- 1
CD38_ADT <- rep(0,ncol(seurat_obj_merge_ADT_data))
CD38_ADT[which(seurat_obj_merge_ADT_data[rownames(seurat_obj_merge_ADT_data) %in% "CD38-Human-CCTATTCCGATTCCG",] > 0)] <- 1
table(CD38_ADT)
CD201_ADT <- rep(0,ncol(seurat_obj_merge_ADT_data))
CD201_ADT[which(seurat_obj_merge_ADT_data[rownames(seurat_obj_merge_ADT_data) %in% "CD201-Human-GTTTCCTTGACCAAG",] > 0)] <- 1
table(CD201_ADT)


##CD34_CD90_EPCR hypothese: all cells are CD34+ RNA level
##CD34_CD90_ hypothese: all cells are CD34+
CD34_CD90_EPCR_ADT <- rep("other",ncol(seurat_obj_merge_ADT_data))
CD34_CD90_EPCR_ADT [which( CD90_ADT > 0 & CD201_ADT > 0 )] <- "CD34+CD90+EPCR+"
CD34_CD90_EPCR_ADT [which( CD90_ADT == 0 & CD201_ADT == 0)] <- "CD34+CD90-EPCR-"
table(CD34_CD90_EPCR_ADT )


#Fig 1: Visualisation of group CD34+CD90+EPCR+/CD34+CD90-EPCR-/other per treatment
# at protein level
seurat_obj_merge$group1 <- CD34_CD90_EPCR_ADT
df_harmony<- as.data.frame(cbind(seurat_obj_merge$treatment,
                                        seurat_obj_merge$group1))
colnames(df_harmony) <- c("treatment","group1")
table(df_harmony$treatment)

#split.by = 'treatment',
p1 <- DimPlot(seurat_obj_merge, reduction = "umap", 
              group.by = "group1", split.by = 'treatment',
              cols = c("darkgreen","purple","#e0e0e0")) + 
  xlim(c(-10,13)) + ylim(c(-10,10))
p1

p1 <- DimPlot(seurat_obj_merge, reduction = "umap", 
              group.by = "group1", 
              cols = c("darkgreen","purple","#e0e0e0")) + 
  xlim(c(-10,13)) + ylim(c(-10,10))
p1

#### Gene expression ####*###**##**#**
DefaultAssay(seurat_obj_merge) <- "RNA"
seurat_obj_merge_RNA_data<- GetAssayData(object = seurat_obj_merge, slot = "counts")
CD34_RNA <- rep(0,ncol(seurat_obj_merge_RNA_data))
CD34_RNA[which(seurat_obj_merge_RNA_data[rownames(seurat_obj_merge_RNA_data) %in% "CD34",] > 0)] <- 1
CD90_RNA <- rep(0,ncol(seurat_obj_merge_RNA_data))
CD90_RNA[which(seurat_obj_merge_RNA_data[rownames(seurat_obj_merge_RNA_data) %in% "THY1",] > 0)] <- 1
CD38_RNA <- rep(0,ncol(seurat_obj_merge_RNA_data))
CD38_RNA[which(seurat_obj_merge_RNA_data[rownames(seurat_obj_merge_RNA_data) %in% "CD38",] > 0)] <- 1
CD201_RNA <- rep(0,ncol(seurat_obj_merge_RNA_data))
CD201_RNA[which(seurat_obj_merge_RNA_data[rownames(seurat_obj_merge_RNA_data) %in% "PROCR",] > 0)] <- 1
table(CD90_RNA)

SIRT1_RNA <- rep(0,ncol(seurat_obj_merge_RNA_data))
SIRT1_RNA[which(seurat_obj_merge_RNA_data[rownames(seurat_obj_merge_RNA_data) %in% "SIRT1",] > 0)] <- 1
table(SIRT1_RNA)
INKA1_RNA <- rep(0,ncol(seurat_obj_merge_RNA_data))
INKA1_RNA[which(seurat_obj_merge_RNA_data[rownames(seurat_obj_merge_RNA_data) %in% "INKA1",] > 0)] <- 1
table(INKA1_RNA)

##CD34_CD90_ hypothese: all cells are CD34+
CD34_CD90_EPCR_RNA <- rep("other",ncol(seurat_obj_merge_ADT_data))
CD34_CD90_EPCR_RNA [which( CD90_RNA >0 & CD201_RNA > 0  )] <- "CD34+CD90+EPCR+"
CD34_CD90_EPCR_RNA [which( CD90_RNA == 0 & CD201_RNA == 0 )] <- "CD34+CD90-EPCR-"
table(CD34_CD90_EPCR_RNA)
table(CD34_CD90_EPCR_RNA)*100/ length(CD34_CD90_EPCR_RNA)

seurat_obj_merge$group1 <- CD34_CD90_EPCR_RNA
df_seurat_obj_merge <- as.data.frame(cbind(seurat_obj_merge$treatment,
                                 seurat_obj_merge$group1))
colnames(df_seurat_obj_merge) <- c("treatment","group1")
table(df_seurat_obj_merge$treatment)
#day0_none day1_ctrl  day1_vpa day4_ctrl  day4_vpa 
#11831     11430     14180     13970     12814 

#Fig 1: Visualisation of group CD34+CD90+EPCR+/CD34+CD90-EPCR-/other per treatment
#at gene expression level
DefaultAssay(seurat_obj_merge) <- "RNA"
p1 <- DimPlot(seurat_obj_merge, reduction = "umap",
              group.by = "group1",split.by = 'treatment',
              cols = c("darkgreen","purple","#e0e0e0")) + 
  xlim(c(-10,13)) + ylim(c(-10,10))
p1

##Marker genes
#Try to select marker genes for each cluster and generate a heatmap of the top five marker genes for each cluster.
seurat_obj_merge <- SetIdent(seurat_obj_merge,value = "seurat_clusters")
clust.markers <- FindAllMarkers(seurat_obj_merge, 
                                only.pos = TRUE,
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
head(clust.markers)
file_cluster_markers <-paste0("/lab_Hoffman/paper2/results/Day0_Day1_Day4_VPA_control_reviewers/batch_correction/clust.markers_harmony.csv")
write.csv(clust.markers,file=file_cluster_markers)

topG <- clust.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)
#
DoHeatmap(seurat_obj_merge, features = topG$gene) + NoLegend()
DoHeatmap(seurat_obj_merge, features = topG$gene) 

#Please display the top marker gene of each cluster in 
#
topG <- clust.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)
top_gene <- topG$gene

#FeaturePlot 
FeaturePlot(seurat_obj_merge,features = top_gene,pt.size = 0.2)

#RidgePlot
RidgePlot(seurat_obj_merge,features = top_gene)

# Seurat Cluster vs treatment
cluster_var <- 'seurat_clusters'
clusters <- unique(seurat_obj_merge@meta.data[[cluster_var]])
clusters <- clusters[order(clusters)]
df <- data.frame(matrix(ncol = 4, nrow = 0))
total_cells_treatment <- table(seurat_obj_merge$treatment)
cn<-names(total_cells_treatment)
for(i in 1:length(clusters)){
  print(i)
  cur_df_tmp <- seurat_obj_merge@meta.data %>% subset(seurat_clusters == clusters[i]) %>% .$treatment %>% table()
  cn_tmp <- names(cur_df_tmp)
  cur_df_tmp2 <- setNames(data.frame(matrix(ncol = 5, nrow = 1)), cn)
  cur_df_tmp2[1,]<-c(0,0,0,0,0)
  for(j in 1:length(cn)){
    for(k in 1:length(cn_tmp)){
      if(cn[j] %in% cn_tmp[k]){
        cur_df_tmp2[1,j] <- cur_df_tmp[k]
      }
    }
  }
  cur_df <- as.data.frame(t(as.data.frame(cur_df_tmp2[1,]/total_cells_treatment)))
  colnames(cur_df)[1]<-"value"
  cur_df$treatment <- rownames(cur_df)
  cur_df$Freq <- cur_df$value * 1/(sum(cur_df$value))
  cur_df$cluster <- clusters[i]
  
  nd<-paste0(cur_df$treatment,"_",cur_df$cluster)
  df1 <- data.frame(cur_df, row.names=nd)
  
  df3 <- rbind(df,df1)
  df <- df3
  print(i)
}
# Specify columns you want to change
#i <- c(1:ncol(df))   
# Specify own function within apply
#df.num <- apply(df[ , i], 2,            
#                   function(x) as.numeric(as.character(x))) 

p <- ggplot(df, aes(y=Freq, x=cluster, fill=treatment)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  ylab('normalized proportion') +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
p

# Seurat Cluster vs group1
cluster_var <- 'seurat_clusters'
clusters <- unique(seurat_obj_merge@meta.data[[cluster_var]])
clusters <- clusters[order(clusters)]
df <- data.frame(matrix(ncol = 3, nrow = 0))
total_cells_group1 <- table(seurat_obj_merge$group1)
cn<-names(total_cells_group1)
for(i in 1:length(clusters)){
  print(i)
  cur_df_tmp <- seurat_obj_merge@meta.data %>% subset(seurat_clusters == clusters[i]) %>% .$group1 %>% table()
  cn_tmp <- names(cur_df_tmp)
  cur_df_tmp2 <- setNames(data.frame(matrix(ncol = 3, nrow = 1)), cn)
  cur_df_tmp2[1,]<-c(0,0,0)
  for(j in 1:length(cn)){
    for(k in 1:length(cn_tmp)){
      if(cn[j] %in% cn_tmp[k]){
        cur_df_tmp2[1,j] <- cur_df_tmp[k]
      }
    }
  }
  cur_df <- as.data.frame(t(as.data.frame(cur_df_tmp2[1,]/total_cells_group1)))
  colnames(cur_df)[1]<-"value"
  cur_df$group1 <- rownames(cur_df)
  cur_df$Freq <- cur_df$value * 1/(sum(cur_df$value))
  cur_df$cluster <- clusters[i]
  
  nd<-paste0(cur_df$group1,"_",cur_df$cluster)
  df1 <- data.frame(cur_df, row.names=nd)
  
  df3 <- rbind(df,df1)
  df <- df3
  print(i)
}
# Specify columns you want to change
#i <- c(1:ncol(df))   
# Specify own function within apply
#df.num <- apply(df[ , i], 2,            
#                   function(x) as.numeric(as.character(x))) 

p <- ggplot(df, aes(y=Freq, x=cluster, fill=group1)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  ylab('normalized proportion') +
  scale_fill_manual(values= c("darkgreen","purple","#e0e0e0")) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
p
