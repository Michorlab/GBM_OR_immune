## load required libraries
library(tidyverse)
library(infercnv)
library(viridis)
library(janitor)


## read in data and clean
# counts
pugh.counts <- read.table("../data/pugh/Richards_NatureCancer_GBM_scRNAseq_counts.csv", sep = ",", header = TRUE)
rownames(pugh.counts) <- pugh.counts$X
pugh.counts <- pugh.counts[,-1]

# metadata
pugh.meta <- read.table("../data/pugh/Richards_NatureCancer_GBM_scRNAseq_meta.csv",sep = ",", header = TRUE)
pugh.meta$X <- str_replace_all(colnames(pugh.counts), "-", ".")
all.equal(pugh.meta$X, colnames(pugh.counts))


## separate data per patient
pdf("../plots/pugh/cellassignment.pdf", width = 10)
pugh.meta %>% 
  mutate(CellType = factor(CellType, levels = c("Immune", "NormalBrain", "Tumour"))) %>% 
  group_by(PatientID, CellType) %>% count(CellType) %>% 
  ggplot(aes(x = PatientID, y = n, fill = CellType)) + geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T, option = "E") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
dev.off()


## run inferCNV
# annotation
annotation <- pugh.meta %>% select(X, CellType, PatientID)
annotation$CellType <- apply(annotation,1, function(x){if(x[2] == "Tumour") return(paste(x[2], x[3], sep = "_")) else return(x[2]) })
r <- annotation$X # keep cell names
annotation <- annotation$CellType # retain only annotation as single column
annotation <- as.data.frame(annotation) # transform again to data frame, currently it is a vector
colnames(annotation) <- NULL # get rid of any potential column names
rownames(annotation) <- r # add cell names as rownames
head(annotation)

# file with genomic position
gene <- read.table(file = "../data/pugh/gencode_v19_gene_pos.txt", sep = "\t")
rownames(gene) <- gene$V1
colnames(gene) <- NULL
gene <- gene[,-1]
head(gene)

# create the inferCNV object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = pugh.counts,
                                     annotations_file = annotation,
                                     gene_order_file= gene,
                                     ref_group_names=c("Immune", "NormalBrain"))

# run inferCNV (takes very long because of the large number of reference cells, needs to be done on a cluster)
infercnv_pugh_tpm_i3 = infercnv::run(infercnv_obj,
                                       cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                       out_dir="../data/inferCNV/pughTPMi3/", 
                                       tumor_subcluster_partition_method = "qnorm",
                                       cluster_by_groups = TRUE, 
                                       window_length = 100,
                                       denoise=TRUE,
                                       HMM_type="i3",
                                       HMM = TRUE,
                                       analysis_mode = "subclusters")


## analyze inferCNV results for EGFR and CDK4
genes.i3 <- read.table("../data/infercnv/pughTPMi3/HMM_CNV_predictions.HMMi3.qnorm.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                       header = TRUE)
groups.i3 <- read.table("../data/infercnv/pughTPMi3/17_HMM_predHMMi3.qnorm.hmm_mode-subclusters.cell_groupings",
                        header = TRUE)
genes.i3.short <- genes.i3 %>% filter(gene %in% c("EGFR", "CDK4")) %>% select(cell_group_name, state, gene)

# check which groupings don't have entries for this region at all
setdiff(groups.i3$cell_group_name, genes.i3$cell_group_name)
missing.groups <- grep("Tumour", setdiff(groups.i3$cell_group_name, genes.i3.short$cell_group_name), value = TRUE)
# add these groups for both genes
genes.i3.short <- rbind(genes.i3.short, cbind("cell_group_name" = missing.groups, "state" = 2, "gene" = "EGFR"))
genes.i3.short <- rbind(genes.i3.short, cbind("cell_group_name" = missing.groups, "state" = 2, "gene" = "CDK4"))


# check which groupings only have 1 entry
groups.once <- genes.i3.short %>% group_by(cell_group_name) %>% summarise(p = table(cell_group_name)) %>% filter(p == 1) %>% pull(cell_group_name)
# add the missing gene for these groups
groups.to.add.single <- t(apply(genes.i3.short[match(groups.once, genes.i3.short$cell_group_name),], 1, function(x){
  if (x[3] == "EGFR")
    return(c(x[1], 2, "CDK4"))
  if (x[3] == "CDK4")
    return(c(x[1], 2, "EGFR"))
}))
colnames(groups.to.add.single) <- c("cell_group_name", "state", "gene")
genes.i3.short <- rbind(genes.i3.short, groups.to.add.single)
genes.i3.short$combine <- paste(genes.i3.short$gene, genes.i3.short$state, sep = "_")
genes.i3.short$patient <- sapply(genes.i3.short$cell_group_name, function(x){
  strsplit(x, "_")[[1]][3]
})

# do 2x2 contingency tables
tabyl.i3 <- genes.i3.short %>% arrange(cell_group_name) %>% tabyl(cell_group_name, combine)
tabyl.i3$patient <- sapply(tabyl.i3$cell_group_name, function(x){
  strsplit(x, "_")[[1]][3]
})

# there's no EGFR_1, so no loss of EGFR
tabyl.i3$EC <- apply(tabyl.i3, 1, function(x){
  # if (x["CDK4_1"] == 1 & x["EGFR_1"] == 1)
  #   return("C.loss_E.loss")
  if (x["CDK4_1"] == 1 & x["EGFR_2"] == 1)
    return("C.loss_E.neutral")
  if (x["CDK4_1"] == 1 & x["EGFR_3"] == 1)
    return("C.loss_E.gain")
  # if (x["CDK4_2"] == 1 & x["EGFR_1"] == 1)
  #   return("C.neutral_E.loss")
  if (x["CDK4_2"] == 1 & x["EGFR_2"] == 1)
    return("C.neutral_E.neutral")
  if (x["CDK4_2"] == 1 & x["EGFR_3"] == 1)
    return("C.neutral.E.gain")
  # if (x["CDK4_3"] == 1 & x["EGFR_1"] == 1)
  #   return("C.gain_E.loss")
  if (x["CDK4_3"] == 1 & x["EGFR_2"] == 1)
    return("C.gain_E.neutral")
  if (x["CDK4_3"] == 1 & x["EGFR_3"] == 1)
    return("C.gain_E.gain")
})

# by cell
tabyl.i3.bycell <- apply(tabyl.i3, 1, function(x){
  cells.to.add <- groups.i3$cell[which(groups.i3$cell_group_name == x["cell_group_name"])]
  return(cbind("cell" = cells.to.add, "patient" = rep(x["patient"], length(cells.to.add)), "EC" = rep(x["EC"], length(cells.to.add))))
})
tabyl.i3.bycell <- do.call(rbind, tabyl.i3.bycell)
tabyl.i3.bycell <- as.tibble(tabyl.i3.bycell)
pdf("../plots/pugh/cells.EC.TPMi3.pdf")
tabyl.i3.bycell %>% ggplot(aes(x = patient, fill = EC)) + geom_bar(position="stack") + 
  scale_fill_viridis(discrete = T, option = "E") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
dev.off()

