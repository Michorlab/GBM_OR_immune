## load required libraries
library(tidyverse)
library(viridis)
library(janitor)
library(Seurat)


## read in data (downloaded from here:https://singlecell.broadinstitute.org/single_cell/study/SCP503/gradient-of-developmental-and-injury-reponse-transcriptional-states-define-functional-vulnerabilities-underpinning-glioblastoma-heterogeneity)
pugh.counts <- read.table("../data/pugh/Richards_NatureCancer_GBM_scRNAseq_counts.csv", sep = ",", header = TRUE)
rownames(pugh.counts) <- pugh.counts$X
pugh.counts <- pugh.counts[,-1]

# metadata
pugh.meta <- read.table("../data/pugh/Richards_NatureCancer_GBM_scRNAseq_meta.csv", sep = ",", header = TRUE)
pugh.meta$X <- str_replace_all(pugh.meta$X, "-", ".")
all.equal(pugh.meta$X, colnames(pugh.counts))


## read in Seurat
pugh <- CreateSeuratObject(counts = pugh.counts, min.cells = 10, min.features = 200)
pugh.meta <- pugh.meta[match(colnames(pugh@assays$RNA@counts), pugh.meta$X),]
all.equal(pugh.meta$X, colnames(pugh@assays$RNA@counts))


## add metadata
pugh$nGene <- pugh.meta$nGene
pugh$nUMI <- pugh.meta$nUMI
pugh$orig.ident <- pugh.meta$orig.ident
pugh$SampleID <- pugh.meta$SampleID
pugh$PatientID <- pugh.meta$PatientID
pugh$SampleType <- pugh.meta$SampleType
pugh$Pathology <- pugh.meta$Pathology
pugh$Stage <- pugh.meta$Stage
pugh$percent.mito <- pugh.meta$percent.mito
pugh$res.1.5 <- pugh.meta$res.1.5
pugh$CellType <- pugh.meta$CellType
pugh$UMAP1 <- pugh.meta$UMAP1
pugh$UMAP2 <- pugh.meta$UMAP2
pugh$Brain.GTEx_AUC <- pugh.meta$Brain.GTEx_AUC
pugh$ESTIMATE.Immune_AUC <- pugh.meta$ESTIMATE.Immune_AUC
pugh$BrainCells <- pugh.meta$BrainCells
pugh$ImmuneCells <- pugh.meta$ImmuneCells
pugh$Classification <- pugh.meta$Classification
pugh$BrainCutoff <- pugh.meta$BrainCutoff
pugh$ImmuneCutoff <- pugh.meta$ImmuneCutoff


## preprocess all
pugh[["percent.mt"]] <- PercentageFeatureSet(pugh, pattern = "MT-")
percent.mito.thresh <- 20
pugh <- pugh[,-which(pugh$percent.mt > percent.mito.thresh)]

pugh <- SCTransform(pugh, verbose = TRUE)
pugh <- RunPCA(pugh, features = VariableFeatures(object = pugh, assay = "SCT"), assay = "SCT")
pugh <- FindNeighbors(pugh, dims = 1:30, verbose = FALSE, assay = "SCT", reduction = "pca", graph.name = "SNN.SCT")
pugh <- FindClusters(pugh, verbose = TRUE, graph.name = "SNN.SCT")
pugh <- RunUMAP(pugh, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.SCT", assay = "SCT")

pugh$SampleShort <- dplyr::recode(pugh$SampleID)


## plot
pdf("../plots/pugh/umap-celltype.pdf") #this uses the label in the downloaded annotation
DimPlot(pugh, group.by = "CellType")
dev.off()
pdf("../plots/pugh/umap-clusters.pdf")
DimPlot(pugh, label = TRUE)
dev.off()
pdf("../plots/pugh/umap-sampleID.pdf")
DimPlot(pugh, group.by = "SampleID")
dev.off()
pugh$Sample.short <- sapply(pugh$SampleID, function(x){strsplit(x, "-")[[1]][1]})
pdf("../plots/pugh/umap-sample.short.pdf")
DimPlot(pugh, group.by = "Sample.short")
dev.off()


## assign cell types and macrophages
canonical.markers <- read.csv(file = "../tables/markers.csv")
canonical.markers <- canonical.markers$marker
DefaultAssay(pugh) <- "SCT"
pdf("../plots/pugh/heatmap.markers.pdf", height = 10, width = 20)
DoHeatmap(pugh, features = canonical.markers)
dev.off()

macrophage.markers <- readxl::read_excel(path = "../tables/immune_markers.xlsx")
macrophage.markers <- macrophage.markers %>% filter(organism == "human")
macrophage.M1 <- macrophage.markers %>% filter(type == "M1")
macrophage.M2 <- macrophage.markers %>% filter(type == "M2")
macrophage.markers <- toupper(macrophage.markers$gene)
macrophage.M1 <- toupper(macrophage.M1$gene)
macrophage.M2 <- toupper(macrophage.M2$gene)

DefaultAssay(pugh) <- "SCT"
pdf("../plots/pugh/heatmap.markers.macrophages.pdf", height = 10, width = 20)
DoHeatmap(pugh, features = macrophage.markers)
dev.off()

DefaultAssay(pugh) <- "RNA"
pugh <- NormalizeData(pugh, assay = "RNA")
pugh$M1sign <- colSums(pugh@assays$RNA@data[na.omit(match(macrophage.M1, rownames(pugh@assays$RNA@data))),])
pugh$M2sign <- colSums(pugh@assays$RNA@data[na.omit(match(macrophage.M2, rownames(pugh@assays$RNA@data))),])

pdf("../plots/pugh/ridge.macrophage.scores.pdf", height = 4)
RidgePlot(pugh, features = c("M1sign", "M2sign"))
dev.off()

pugh$celltypes2 <- dplyr::recode(pugh$seurat_clusters, "0" = "immune", "1" = "immune", "2" = "immune", "3" = "non-immune",
              "4" = "immune", "5" = "non-immune", "6" = "non-immune", "7" = "non-immune", "8" = "immune",
              "9" = "non-immune", "10" = "non-immune", "11" = "immune", "12" = "immune", "13" = "immune",
              "14" = "non-immune", "15" = "immune", "16" = "immune", "17" = "non-immune", "18" = "non-immune")
pdf("../plots/pugh/umap-assignment-broad.pdf") # this is the celltype assigned by us
DimPlot(pugh, group.by = "celltypes2")
dev.off()

pugh$immune.detail <- dplyr::recode(pugh$seurat_clusters, "0" = "macrophage", "1" = "macrophage", "2" = "macrophage", "3" = "non-immune",
                                 "4" = "macrophage", "5" = "non-immune", "6" = "non-immune", "7" = "non-immune", "8" = "macrophage",
                                 "9" = "non-immune", "10" = "non-immune", "11" = "macrophage", "12" = "T-cells", "13" = "macrophage",
                                 "14" = "non-immune", "15" = "macrophage", "16" = "macrophage", "17" = "non-immune", "18" = "non-immune")
pdf("../plots/pugh/umap-assignment-immune.pdf")
DimPlot(pugh, group.by = "immune.detail")
dev.off()


pugh$macrophage.detail <- dplyr::recode(pugh$seurat_clusters, "0" = "M1", "1" = "M1", "2" = "M1", "3" = "non-immune",
                                    "4" = "M1", "5" = "non-immune", "6" = "non-immune", "7" = "non-immune", "8" = "M2",
                                    "9" = "non-immune", "10" = "non-immune", "11" = "M2", "12" = "T-cells", "13" = "M2",
                                    "14" = "non-immune", "15" = "M2", "16" = "M2", "17" = "non-immune", "18" = "non-immune")
pdf("../plots/pugh/umap-assignment-macrophage.pdf")
DimPlot(pugh, group.by = "macrophage.detail")
dev.off()


Idents(pugh) <- pugh$macrophage.detail
my_levels <- c("M1", "M2", "T-cells","non-immune")
pugh@active.ident <- factor(x = pugh@active.ident, levels = my_levels)
DefaultAssay(pugh) <- "SCT"
pdf("../plots/pugh/heatmap.markers.percelltype.pdf", height = 10, width = 20)
DoHeatmap(pugh, features = canonical.markers)
dev.off()


## compare among samples
matSC <- cbind("macrophage" = as.character(pugh$macrophage.detail), "sample" = pugh$Sample.short)
matSC <- as.data.frame(matSC)

pdf("../plots/pugh/barplot-macrophages-sample.pdf")
matSC %>% ggplot(aes(x = sample)) + 
  geom_bar(aes(fill = macrophage), stat="count", position = position_stack()) + scale_fill_viridis(discrete = T, option = "A") + theme_bw()
dev.off()

pdf("../plots/pugh/barplot-macrophages-noimmune-sample.pdf")
matSC %>% filter(!macrophage %in% c("non-immune")) %>% ggplot(aes(x = sample)) + 
  geom_bar(aes(fill = macrophage), stat="count", position = position_stack()) + scale_fill_viridis(discrete = T, option = "A") + theme_bw()
dev.off()

pdf("../plots/pugh/barplot-macrophages-noimmune-percentage-sample.pdf")
matSC %>% filter(!macrophage %in% c("non-immune")) %>% ggplot(aes(x = sample, fill = macrophage)) + 
  geom_bar(aes(fill = macrophage), stat="count", position = "fill") + scale_fill_viridis(discrete = T, option = "A") + theme_bw()
dev.off()