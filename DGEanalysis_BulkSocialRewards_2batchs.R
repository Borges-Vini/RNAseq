# RNAseq Analysis
# Vini Borges
# Nov 6, 2024

# Packages ----------------------------------------------------------------

install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

install.packages("remotes")
remotes::install_github("rqtl/qtl2")
library(umap)
library(DESeq2)
library(ggplot2)
library(Rsubread)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(reshape2)
library(edgeR)
library(limma)
library(ssviz)
library(biomaRt)
library(gridExtra)
library(tidyverse)
library(pheatmap)
library(grid)
library(gridExtra)
library(apeglm)
library(qtl2)
library(variancePartition)
library(msigdbr)
library(enrichR)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(purrr)

# Loading Samples ---------------------------------------------------------

#Set WD
setwd("/Users/vborges/Library/CloudStorage/OneDrive-MarshallUniversity/Marshall/BulkRNA_Socialrewards/rawfiles_sample_bam/")

#Results Path
path.results <- "/Users/vborges/Library/CloudStorage/OneDrive-MarshallUniversity/Marshall/BulkRNA_Socialrewards/results/"

#Load the sample sheet
sample_info <- read.csv("Social_FASTQ_Sample_Guide.csv", row.names = "Subject")

#Converting to factor
sample_info$Sex <- factor(sample_info$Sex)
sample_info$Age <- gsub("w ", "", sample_info$Age); sample_info$Age <- factor(sample_info$Age)
sample_info$Mouseline <- gsub("Inbred ", "", sample_info$Mouseline); sample_info$Mouseline <- factor(sample_info$Mouseline)

#Set up "Condition" column
sample_info$Condition <- factor(gsub("Social - ", "", sample_info$Experiment))

#Creating object by strain
samples_by_mouseline <- split(sample_info, sample_info$Mouseline)

#Adding colors based on strain
CCcolors <- c(
  AJ = "#F0E442", B6 = "#555555", `129` = "#E69F00", NOD = "#0072B2", NZO = "#56B4E9", 
  CAST = "#009E73", PWK = "#D55E00", WSB = "#CC79A7",
  DBA = "#999999", CC019 = "#CC0000", CC002 = "#6699CC", CC005 = "#66CC99", CC051 = "#FFCC33"
)
sample_info$Colors <- NA
sample_info$Colors <- CCcolors[match(sample_info$Mouseline, names(CCcolors))]

#Checkpoint 1
str(sample_info)
table(sample_info$Mouseline, sample_info$Sex, sample_info$Condition)
write.table(sample_info, paste0(path.results, "samples_info.txt"), sep = "\t", row.names = T, col.names = T, quote = F)
writeLines(capture.output(table(sample_info$Mouseline, sample_info$Sex, sample_info$Condition)), paste0(path.results, "samples_info_tab.txt"))

# Loading BAM files ---------------------------------------------------------------

#Define BAM file paths
bam_files <- list.files(getwd(), pattern = "*.bam$", full.names = TRUE)

#Checking BAMs
readBam(bam_files[[1]], tags = character(0))

#Male Only
#male_samples <- sample_info[sample_info$Sex == "Male",]
#bam_files_male <- bam_files[sub("_aligned.bam$", "", basename(bam_files)) %in% male_samples$Sample.ID]


#sample_info <- sample_info[sample_info$Batch == 1,]

# Counts ------------------------------------------------------------------
#Run featureCounts to get the counts
#fc <- featureCounts(files = bam_files, 
                    # annot.ext = "../Mus_musculus.GRCm38.102.gtf",
                    # isGTFAnnotationFile = TRUE, 
                    # useMetaFeatures = TRUE,
                    # isPairedEnd = TRUE)
#save(fc, file=paste0(path.results,"fc.RData"))
load(file = paste0(path.results,"fc.RData"))

counts <- fc$counts
colnames(counts) <- gsub("^(\\d+).*", "\\1-STRI", colnames(counts))

# Matrix for Limma-Voom ---------------------------------------------------
sample_info <- sample_info[match(colnames(counts), sample_info$Sample_ID), ]
table(sample_info$Mouseline, sample_info$Condition)
sample_info$Mouseline <- factor(sample_info$Mouseline, levels = c("AJ", "B6", "CC002", "CC005", "CC019", "CC051", "DBA", "NOD", "PWK"))
sample_info$Condition <- factor(sample_info$Condition, levels = c("Nonreward", "Reward"))
sample_info$Sex <- factor(sample_info$Sex, levels = c("Female", "Male"))
sample_info <- droplevels(sample_info)

mm <- model.matrix(~0 + Mouseline * Condition + Sex, data = sample_info)

mm <- model.matrix(
  ~ 0 + Batch + Sex + Condition + Mouseline + 
    Sex:Condition + Sex:Mouseline + Condition:Mouseline + 
    Sex:Condition:Mouseline,
  data = sample_info,
  contrasts.arg = list(
    Sex = contrasts(sample_info$Sex, contrasts = FALSE),
    Condition = contrasts(sample_info$Condition, contrasts = FALSE),
    Mouseline = contrasts(sample_info$Mouseline, contrasts = FALSE)
  )
)
head(mm)


dge <- DGEList(counts)
dge <- calcNormFactors(dge)
dim(dge)

#mm <- model.matrix(~0 + Batch + Condition:Sex, data = sample_info)

# Matrix for DESeq2 -------------------------------------------------------

#design <- model.matrix(~ Batch + Sex + Condition + Mouseline + Mouseline:Condition + Sex:Condition + Sex:Mouseline + Sex:Condition:Mouseline, data = sample_info)
 design <- model.matrix(~0 + Batch + Sex:Condition:Mouseline, data = sample_info)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = design,
)

summary(dds)

#Convert DESeq2 object to DGEList 
dge <- DGEList(counts = counts(dds))

# LRT ---------------------------------------------------------------------

lrt_design <- model.matrix(~ Batch + Sex + Condition + Mouseline + Sex:Condition + Sex:Mouseline, data = sample_info)

dds_lrt <- DESeq(dds, test = "LRT",
                 reduced = lrt_design)

res_lrt <- results(dds_lrt)
sig_lrt <- res_lrt[which(res_lrt$padj < 0.05), ]
sig_lrt <- sig_lrt[order(abs(sig_lrt$log2FoldChange), decreasing = TRUE), ]
head(sig_lrt)
sig_lrt <- as.data.frame(sig_lrt)

gene_symbols <- as.data.frame(mapIds(org.Mm.eg.db, 
                                     keys = rownames(sig_lrt), 
                                     column = "SYMBOL", 
                                     keytype = "ENSEMBL", 
                                     multiVals = "first"))
colnames(gene_symbols) <- "Symbol"

sig_lrt$Symbol <- gene_symbols[rownames(sig_lrt), "Symbol"]
sig_lrt <- sig_lrt[, c(7, 1:6)]
#save(sig_lrt, file=paste0(path.results,"Mouseline.Condition_sig_lrt.RData"))
#load(file = paste0(path.results,"Sex.Condition.Mouseline_sig_lrt.RData"))

#rownames(sig_lrt) <- paste(rownames(sig_lrt), gene_symbols[rownames(sig_lrt), "Symbol"], sep = " - ")
#write.csv(as.data.frame(sig_lrt), paste0(path.results, "Mouseline.Condition_LRT_results.csv"))

#PLOT 1
norm_counts <- counts(dds_lrt, normalized = TRUE)
gene_to_plot <- rownames(sig_lrt)[4] 

plot_df <- data.frame(
  expression = as.numeric(norm_counts[gene_to_plot, ]),
  sample = colnames(norm_counts),
  stringsAsFactors = FALSE
)

plot_df <- cbind(plot_df, as.data.frame(colData(dds_lrt)))
plot_df$Group <- with(plot_df, interaction(Sex, Condition, Mouseline, sep = "_"))

meta <- as.data.frame(colData(dds_lrt)[, c("Sex", "Condition", "Mouseline")])
names(meta) <- make.names(names(meta), unique = TRUE)

mouselines <- c("AJ","B6","CC002","CC005","CC019","CC051","DBA","NOD","PWK")
sexes <- c("Female","Male")
conds <- c("Nonreward","Reward")

desired_order <- unlist(lapply(conds, function(con) {
  c(paste("Female","AJ", con, sep = "_"),
    paste("Female","B6", con, sep = "_"),
    paste("Female","CC002", con, sep = "_"),
    paste("Female","CC005", con, sep = "_"),
    paste("Female","CC019", con, sep = "_"),
    paste("Female","CC051", con, sep = "_"),
    paste("Female","DBA", con, sep = "_"),
    paste("Female","NOD", con, sep = "_"),
    paste("Female","PWK", con, sep = "_"),
    paste("Male","AJ", con, sep = "_"),
    paste("Male","B6", con, sep = "_"),
    paste("Male","CC002", con, sep = "_"),
    paste("Male","CC005", con, sep = "_"),
    paste("Male","CC019", con, sep = "_"),
    paste("Male","CC051", con, sep = "_"),
    paste("Male","DBA", con, sep = "_"),
    paste("Male","NOD", con, sep = "_"),
    paste("Male","PWK", con, sep = "_"))
}))
present_groups <- unique(as.character(plot_df$Group))
present_order <- desired_order[desired_order %in% present_groups]

plot_df$Group <- factor(plot_df$Group, levels = present_order)
plot_df$Cluster <- with(plot_df, paste(Sex, Condition, sep = "_"))

is_ggplot <- function(x) inherits(x, "gg")
png(paste0(path.results,"plot1d_",gene_to_plot,".png"), width = 3000, height = 1500, res = 200) 
ggplot(plot_df, aes(x = Mouseline, y = expression, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, aes(color = Condition)) +
  facet_grid(. ~ Cluster, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    title = paste("Expression of", gene_to_plot),
    x = "Mouseline",
    y = "Normalized Counts"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 60)
  )
dev.off()

#PLOT 2
volcano_data <- as.data.frame(sig_lrt)

volcano_data$threshold <- "Not significant"
volcano_data$threshold[volcano_data$log2FoldChange >= 1 & volcano_data$padj < 0.05] <- "Upregulated"
volcano_data$threshold[volcano_data$log2FoldChange <= -1 & volcano_data$padj < 0.05] <- "Downregulated"

colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")

png(paste0(path.results,paste0("plot2d.png")), width = 3000, height = 1500,  res = 200) 
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot (Full model vs. Reduced model)",
    x = "Log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Regulation"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme(legend.position = "top")
dev.off()


#PLOT 3
sig_lrt$baseMean_log <- log10(sig_lrt$baseMean + 1)
png(paste0(path.results,paste0("plot3d.png")), width = 3000, height = 1500,  res = 200) 
ggplot(sig_lrt, aes(x = baseMean_log, y = log2FoldChange)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "MA Plot: Full model vs. Reduced model LRT",
       x = "log10(baseMean)",
       y = "log2 Fold Change")
dev.off()


#PLOT 4 
top_genes <- rownames(head(sig_lrt, 30))

# Extract variance-stabilized expression
vsd <- vst(dds_lrt, blind = FALSE)
mat <- assay(vsd)[top_genes, ]

# Scale by gene
mat_scaled <- t(scale(t(mat)))

# Get annotatiob
annotation_df <- as.data.frame(colData(dds_lrt)[, c("Condition", "Sex", "Mouseline")])
annotation_df <- annotation_df[colnames(mat_scaled), , drop = FALSE]

png(paste0(path.results,paste0("plot4d.png")), width = 5000, height = 2000,  res = 300) 
pheatmap(mat_scaled,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         fontsize_row = 7,
         fontsize_col = 4,
         main = "Top Genes from Full model vs. Reduced model (Mouseline:Condition)")
dev.off()

#PLOT5

top_n <- 50
norm_counts <- counts(dds_lrt, normalized = TRUE)
top_genes <- rownames(head(sig_lrt, top_n))
norm_counts <- norm_counts[top_genes, ]
plot_df <- as.data.frame(norm_counts)
plot_df$gene <- rownames(plot_df)

plot_df <- plot_df %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  left_join(as.data.frame(colData(dds_lrt)) %>%
              rownames_to_column("sample"), by = "sample")

# Calculate Reward - Nonreward per Mouseline and Sex
diff_df <- plot_df %>%
  group_by(gene, Sex, Mouseline, Condition) %>%
  summarise(mean_expr = mean(expression), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = mean_expr) %>%
  mutate(diff = Reward - Nonreward)

# Make matrix for heatmap: rows = genes, columns = Sex:Mouseline
heatmap_mat <- diff_df %>%
  mutate(sex_line = paste(Sex, Mouseline, sep = "_")) %>%
  select(gene, sex_line, diff) %>%
  pivot_wider(names_from = sex_line, values_from = diff) %>%
  column_to_rownames("gene") %>%
  as.matrix()

heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

annotation_col <- data.frame(
  Sex = sub("_.*", "", colnames(heatmap_mat)),
  Mouseline = sub(".*_", "", colnames(heatmap_mat))
)
rownames(annotation_col) <- colnames(heatmap_mat)

# Plot heatmap
png(paste0(path.results, "interaction_heatmap_diff.png"),
    width = 3000, height = 2500, res = 300)

pheatmap(heatmap_mat_scaled,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         fontsize_col = 8,
         main = "Interaction Heatmap: (Reward - Nonreward) per Sex & Mouseline",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# Filtering ---------------------------------------------------------------

#Checkpoint 2
png(paste0(path.results,"boxplot_before_filtering.png"), width = 2000, height = 1500,  res = 300) 
boxplot(log2(cpm(dge) + 1), main = "Before Filtering")
dev.off()

# incorporate automatic way to define min.count

#Filter data on lower count rate
cutoff <- 3
drop <- which(apply(cpm(dge), 1, max) < cutoff)
dge <- dge[-drop,]
dim(dge)

#Filter low counts
keep_low_count <- filterByExpr(dge, design = design, group = sample_info$Mouseline, min.count = 10) 
table(keep_low_count)
dge <- dge[keep_low_count, , keep.lib.sizes = FALSE]
dim(dge)


#Filter median expression
keep_median_expr <- rowMedians(cpm(dge)) >2
table(keep_median_expr)
dge <- dge[keep_median_expr, , keep.lib.sizes = FALSE]
dim(dge)
#[1] 14728   292

v <- voom(dge, design, plot = TRUE)

dds <- dds[-drop, ]
dds <- dds[keep_low_count, ]
dds <- dds[keep_median_expr, ]
dds$Condition <- factor(dds$Condition, levels = c("Nonreward","Reward"))

#Checkpoint 3
png(paste0(path.results,"boxplot_after_filtering.png"), width = 2000, height = 1500,  res = 300) 
boxplot(log2(cpm(dge) + 1), main = "After Filtering")
dev.off()

# Plot: MDS ---------------------------------------------------------------

#1
png(paste0(path.results, "MDS1_by_mouseline_condition.png"), width = 2000, height = 1500, res = 300)
colors <- CCcolors[as.character(sample_info$Mouseline)]  
shapes <- ifelse(sample_info$Condition == "Reward", 16, 15)
mds <- plotMDS(dge, gene.selection = "common", plot = FALSE)

plot(mds$x, mds$y, 
     col = colors, 
     pch = shapes, 
     cex = 1.5, 
     xlim = c( 1.8, 2.8),  # Set X-axis limits
     ylim = c(-1, 1), # Set Y-axis limits
     xlab = "Leading logFC dim 1", 
     ylab = "Leading logFC dim 2", 
     main = "MDS Plot by Mouseline and Condition")

legend("topleft", 
       legend = c("Reward", "Nonreward"), 
       pch = c(16, 15),  # Matches shapes used
       col = "black", 
       title = "Condition",
       cex = 0.8)

legend("bottomright",                              
       legend = names(CCcolors[names(CCcolors) %in% sample_info$Mouseline]),               
       fill = CCcolors[names(CCcolors) %in% sample_info$Mouseline],                        
       title = "Mouseline",                   
       border = NA,                            
       cex = 0.8)

dev.off()

#2
#Multi-MDS per strain

unique_strains <- unique(sample_info$Mouseline)
n_strains <- length(unique_strains)
n_col <- 3
n_row <- ceiling(n_strains / n_col)
mds <- plotMDS(dge, gene.selection = "common", plot = FALSE)

png(paste0(path.results, "AllStrains_MDS_by_condition.png"), width = 7000, height = 2500 * n_row, res = 900)
par(mfrow = c(n_row, n_col), mar = c(4, 4, 4, 2))

for (strain in unique_strains) {
  strain_idx <- which(sample_info$Mouseline == strain)
  x <- mds$x[strain_idx]
  y <- mds$y[strain_idx]
  color <- CCcolors[as.character(strain)]
  shape <- ifelse(sample_info$Condition[strain_idx] == "Reward", 16, 15)
  
  plot(x, y,
       col = color,
       pch = shape,
       cex = 1.2,
       xlim = c(1.8, 2.8),  # Adjust if needed
       ylim = c(-1, 1),
       xlab = "logFC dim 1",
       ylab = "logFC dim 2",
       main = strain)
  
  legend("bottomright",
         legend = c("Reward", "Nonreward"),
         pch = c(16, 15),
         col = "black",
         cex = 0.8,
         bty = "n")
}
dev.off()

#3
png(paste0(path.results, "MDS1_by_sex_conditions.png"), width = 2000, height = 1500, res = 300)
colors <- ifelse(sample_info$Sex == "Female", "red", "blue")  
shapes <- ifelse(sample_info$Condition == "Reward", 16,15)
mds <- plotMDS(dge, gene.selection = "common", plot = FALSE)

plot(mds$x, mds$y, 
     col = colors, 
     pch = shapes, 
     cex = 1.5, 
     xlim = c( 1.8, 2.8),  # Set X-axis limits
     ylim = c(-1, 1), # Set Y-axis limits
     xlab = "Leading logFC dim 1", 
     ylab = "Leading logFC dim 2", 
     main = "MDS Plot by Sex and Condition")

legend("topright",
       legend = c("Male - Reward", "Male - Nonreward", "Female - Reward", "Female - Nonreward"),
       col = c("blue", "blue", "red", "red"),
       pch = c(16, 15, 16, 15),
       pt.cex = 1.1,    # Size of the points
       cex = 0.7,       # Size of the text (this shrinks the legend)
       box.lty = 0,
       inset = 0.02)

dev.off()



# Plot: PCA / Figure 1 ---------------------------------------------------------------

#1
is_ggplot <- function(x) inherits(x, "gg")

group <- interaction(sample_info$Batch)
pca.raw.d <- log2(dge$counts + 0.5)
pca.d <- PCA(t(pca.raw.d), graph = FALSE)

fviz_pca_ind(pca.d,
             geom.ind = "point",         # use points only, no text
             col.ind = group,            # color by Condition
             palette = "jco",            # nicer color palette
             addEllipses = TRUE,         # add ellipses to show groupings
             legend.title = "Batch",
             pointsize = 2,
             repel = TRUE,               # better label spacing (if using labels)
             title = "PCA: Batch") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

#Sex
is_ggplot <- function(x) inherits(x, "gg")

group <- interaction(sample_info$Sex)
pca.raw.d <- log2(dge$counts + 0.5)
batch_corrected <- removeBatchEffect(pca.raw.d, batch = sample_info$Batch)
pca.d.corrected <- PCA(t(batch_corrected), graph = FALSE)
pca_df <- as.data.frame(pca.d.corrected$ind$coord)
pca_df$Sample_ID <- rownames(pca_df)
pca_df$Sex <- sample_info$Sex
pca_df$Color <- ifelse(sample_info$Sex == "Female", "red", "blue")  
eig <- pca.d.corrected$eig

mouseline_palette <- setNames(ifelse(sample_info$Sex == "Female", "red", "blue") , sample_info$Sex)
mouseline_palette <- mouseline_palette[!duplicated(names(mouseline_palette))]

png(paste0(path.results,"PCA_sex.png"), width = 2000, height = 1500,  res = 300) 
ggplot(pca_df, aes(x = Dim.1, y = Dim.2, color = Sex)) +
  geom_point(shape = 19, size = 2) +
  scale_color_manual(values = mouseline_palette) +
  labs(
    title = "PCA: Striatum",
    x = paste0("PC1 (", round(eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(eig[2, 2], 1), "%)"),
    color = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
dev.off()

#Mouseline
is_ggplot <- function(x) inherits(x, "gg")

group <- interaction(sample_info$Mouseline)
pca.raw.d <- log2(dge$counts + 0.5)
batch_corrected <- removeBatchEffect(pca.raw.d, batch = sample_info$Batch)
pca.d.corrected <- PCA(t(batch_corrected), graph = FALSE)
pca_df <- as.data.frame(pca.d.corrected$ind$coord)
pca_df$Sample_ID <- rownames(pca_df)
pca_df$Mouseline <- sample_info$Mouseline
pca_df$Color <- sample_info$Colors
eig <- pca.d.corrected$eig

mouseline_palette <- setNames(sample_info$Colors, sample_info$Mouseline)
mouseline_palette <- mouseline_palette[!duplicated(names(mouseline_palette))]

ggplot(pca_df, aes(x = Dim.1, y = Dim.2, color = Mouseline)) +
  geom_point(shape = 19, size = 2) +
  scale_color_manual(values = mouseline_palette) +
  labs(
    title = "PCA: Striatum",
    x = paste0("PC1 (", round(eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(eig[2, 2], 1), "%)"),
    color = "Mouseline"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")


fviz_pca_ind(pca.d.corrected,
             col.ind = sample_info$Age,  # or Age, if available
             addEllipses = F,
             title = "PCA colored by Age")


set.seed(123)
kclusters <- kmeans(umap_result$layout, centers = 9) 
table(kclusters$cluster, sample_info$Mouseline)

sample_info$Cluster <- factor(kclusters$cluster)
umap_result <- umap(t(batch_corrected))
plot(umap_result$layout, col = kclusters$cluster, pch = 19,
     main = "UMAP with K-means Clusters")


plot(umap_result$layout, col = as.factor(sample_info$Mouseline), pch = 19,
     main = "UMAP colored by Mouseline")
legend("topright", legend = levels(as.factor(sample_info$Mouseline)), 
       col = 1:length(levels(as.factor(sample_info$Mouseline))), pch = 19, cex = 0.6)


#3
pca_coords <- as.data.frame(pca.d.corrected$ind$coord)
pca_coords$SampleID <- rownames(pca_coords)
pca_coords$Condition <- factor(sample_info$Condition, levels = c("Nonreward", "Reward"))

# Calculate distance from center
pca_coords$dist <- sqrt(pca_coords$Dim.1^2 + pca_coords$Dim.2^2)

# Label top 5 most distant points
label_samples <- pca_coords %>%
  arrange(desc(dist)) %>%
  slice(1:5)

# ggplot-based PCA with labeled outliers
ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = Condition)) +
  geom_point(size = 2) +
  geom_text(data = label_samples, aes(label = SampleID), vjust = -1, size = 3, show.legend = FALSE) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA: Reward vs Nonreward",
       x = paste0("PC1 (", round(pca.d$eig[1, 2], 1), "%)"),
       y = paste0("PC2 (", round(pca.d$eig[2, 2], 1), "%)"))

#4
pca_coords$Mouseline <- sample_info$Mouseline

ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = Condition, shape = Mouseline)) +
  geom_point(size = 2) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA: Reward vs Nonreward by Mouseline")


# Normalizing Data --------------------------------------------------------

#Run DESeq normalization
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)
results(dds)


# Voom & Fit --------------------------------------------------------------

#Checkpoint 4
png(paste0(path.results,"Distribution of Gene Counts.png"), width = 2000, height = 1500,  res = 300) 
hist_data <- hist(log2(rowSums(normalized_counts + 1)), main = "Distribution of Gene Counts", xlab = "Log2 Counts", breaks = 100,col = "lightblue", border = "black")
axis(1, at = seq(min(hist_data$breaks), max(hist_data$breaks), length.out = 6), 
     labels = round(seq(min(hist_data$breaks), max(hist_data$breaks), length.out = 6), 2))
dev.off()

#Checkpoint 5
png(paste0(path.results,"voom.png"), width = 2000, height = 1500,  res = 300) 
v <- voom(normalized_counts, design, plot = TRUE)
dev.off()

#Fit the model
fit <- lmFit(v, design)
colnames(fit$coefficients) <- make.names(colnames(fit$coefficients))
colnames(fit$design) <- make.names(colnames(fit$design))
colnames(fit$stdev.unscaled) <- make.names(colnames(fit$stdev.unscaled))
colnames(fit$cov.coefficients) <- make.names(colnames(fit$cov.coefficients))
print(colnames(coef(fit))) 
head(coef(fit))

# Contrast Matrix :Sex:Condition:Mouseline" -------------------------------------------------

# List of strains
strains <- c("B6", "AJ", "CC002", "CC005", "CC019", "CC051", "DBA", "NOD", "PWK")
sexes <- c("Female", "Male")

# 1. Reward vs Non-reward within each sex and strain
within_strain_sex_contrasts <- list()
for (strain in strains) {
  for (sex in sexes) {
    name <- paste0("Reward_vs_Nonreward_", sex, "_", strain)
    coef_reward <- paste0("Sex", sex, ".ConditionReward.Mouseline", strain)
    coef_nonreward <- paste0("Sex", sex, ".ConditionNonreward.Mouseline", strain)
    within_strain_sex_contrasts[[name]] <- paste0(coef_reward, " - ", coef_nonreward)
  }
}

# 2. Sex difference in reward effect within each strain
sex_diff_reward_effects <- list()
for (strain in strains) {
  name <- paste0("SexDiff_RewardEffect_", strain)
  coef_MR <- paste0("SexMale.ConditionReward.Mouseline", strain)
  coef_MN <- paste0("SexMale.ConditionNonreward.Mouseline", strain)
  coef_FR <- paste0("SexFemale.ConditionReward.Mouseline", strain)
  coef_FN <- paste0("SexFemale.ConditionNonreward.Mouseline", strain)
  formula <- paste0("(", coef_MR, " - ", coef_MN, ") - (", coef_FR, " - ", coef_FN, ")")
  sex_diff_reward_effects[[name]] <- formula
}

# Combine all contrasts
all_contrasts <- c(within_strain_sex_contrasts, sex_diff_reward_effects)

# Generate contrast matrix
contrast_matrix <- makeContrasts(
  contrasts = all_contrasts,
  levels = colnames(coef(fit))
)

colnames(contrast_matrix) <- names(all_contrasts)


# Contrast Matrix "Condition:Mouseline" -------------------------------------------------

# List of strains
strains <- c("B6", "AJ", "CC002", "CC005", "CC019", "CC051", "DBA", "NOD", "PWK")

within_strain_contrasts <- sapply(strains, function(s) {
  paste0("ConditionReward.Mouseline", s, " - ConditionNonreward.Mouseline", s)
})
names(within_strain_contrasts) <- paste0("Reward_vs_NonReward_", strains)

pairwise_contrasts <- list()
for (i in 1:(length(strains)-1)) {
  for (j in (i+1):length(strains)) {
    s1 <- strains[i]
    s2 <- strains[j]
    name <- paste0("RewardEffect_", s1, "_vs_", s2)
    formula <- paste0(
      "(ConditionReward.Mouseline", s1, " - ConditionNonreward.Mouseline", s1, ") - ",
      "(ConditionReward.Mouseline", s2, " - ConditionNonreward.Mouseline", s2, ")"
    )
    pairwise_contrasts[[name]] <- formula
  }
}

all_contrasts <- c(as.list(within_strain_contrasts), pairwise_contrasts)

contrast_matrix <- makeContrasts(
  contrasts = all_contrasts,
  levels = colnames(coef(fit))
)

head(contrast_matrix)

colnames(contrast_matrix)[1]  <- "Reward_vs_NonReward_B6"
colnames(contrast_matrix)[2]  <- "Reward_vs_NonReward_AJ"
colnames(contrast_matrix)[3]  <- "Reward_vs_NonReward_CC002"
colnames(contrast_matrix)[4]  <- "Reward_vs_NonReward_CC005"
colnames(contrast_matrix)[5]  <- "Reward_vs_NonReward_CC019"
colnames(contrast_matrix)[6]  <- "Reward_vs_NonReward_CC051"
colnames(contrast_matrix)[7]  <- "Reward_vs_NonReward_DBA"
colnames(contrast_matrix)[8]  <- "Reward_vs_NonReward_NOD"
colnames(contrast_matrix)[9]  <- "Reward_vs_NonReward_PWK"

colnames(contrast_matrix)[10]  <- "B6_vs_AJ"
colnames(contrast_matrix)[11]  <- "B6_vs_CC002"
colnames(contrast_matrix)[12]  <- "B6_vs_CC005"
colnames(contrast_matrix)[13]  <- "B6_vs_CC019"
colnames(contrast_matrix)[14]  <- "B6_vs_CC051"
colnames(contrast_matrix)[15]  <- "B6_vs_DBA"
colnames(contrast_matrix)[16]  <- "B6_vs_NOD"
colnames(contrast_matrix)[17]  <- "B6_vs_PWK"

colnames(contrast_matrix)[18]  <- "AJ_vs_CC002"
colnames(contrast_matrix)[19]  <- "AJ_vs_CC005"
colnames(contrast_matrix)[20]  <- "AJ_vs_CC019"
colnames(contrast_matrix)[21]  <- "AJ_vs_CC051"
colnames(contrast_matrix)[22]  <- "AJ_vs_DBA"
colnames(contrast_matrix)[23]  <- "AJ_vs_NOD"
colnames(contrast_matrix)[24]  <- "AJ_vs_PWK"

colnames(contrast_matrix)[25]  <- "CC002_vs_CC005"
colnames(contrast_matrix)[26]  <- "CC002_vs_CC019"
colnames(contrast_matrix)[27]  <- "CC002_vs_CC051"
colnames(contrast_matrix)[28]  <- "CC002_vs_DBA"
colnames(contrast_matrix)[29]  <- "CC002_vs_NOD"
colnames(contrast_matrix)[30]  <- "CC002_vs_PWK"

colnames(contrast_matrix)[31]  <- "CC005_vs_CC019"
colnames(contrast_matrix)[32]  <- "CC005_vs_CC051"
colnames(contrast_matrix)[33]  <- "CC005_vs_DBA"
colnames(contrast_matrix)[34]  <- "CC005_vs_NOD"
colnames(contrast_matrix)[35]  <- "CC005_vs_PWK"

colnames(contrast_matrix)[36]  <- "CC019_vs_CC051"
colnames(contrast_matrix)[37]  <- "CC019_vs_DBA"
colnames(contrast_matrix)[38]  <- "CC019_vs_NOD"
colnames(contrast_matrix)[39]  <- "CC019_vs_PWK"

colnames(contrast_matrix)[40]  <- "CC051_vs_DBA"
colnames(contrast_matrix)[41]  <- "CC051_vs_NOD"
colnames(contrast_matrix)[42]  <- "CC051_vs_PWK"

colnames(contrast_matrix)[43]  <- "DBA_vs_NOD"
colnames(contrast_matrix)[44]  <- "DBA_vs_PWK"
colnames(contrast_matrix)[45]  <- "NOD_vs_PWK"

# Contrast matrix "Condition:Sex" ------------------------------------------------

contrast_matrix <- makeContrasts(
  Reward_vs_Nonreward_Female = ConditionReward.SexFemale - ConditionNonreward.SexFemale,
  Reward_vs_Nonreward_Male   = ConditionReward.SexMale - ConditionNonreward.SexMale,
  Female_vs_Male_Reward      = ConditionReward.SexFemale - ConditionReward.SexMale,
  Female_vs_Male_Nonreward   = ConditionNonreward.SexFemale - ConditionNonreward.SexMale,
  Condition_by_Sex_Interaction = 
    (ConditionReward.SexFemale - ConditionNonreward.SexFemale) -
    (ConditionReward.SexMale - ConditionNonreward.SexMale),
  levels = colnames(coef(fit))
)

 rownames(contrast_matrix) <- colnames(coef(fit))

# Contrast matrix "Condition" ------------------------------------------------
 
 contrast_matrix <- makeContrasts(
   Reward_vs_Nonreward = ConditionReward - ConditionNonreward,
   levels = colnames(coef(fit))
 )
 
 rownames(contrast_matrix) <- colnames(coef(fit))
 
# Fit2 & Filtering  ----------------------------------------------

#Checkpoint 
all(rownames(contrast_matrix) == colnames(coef(fit)))

#Estimate contrast for each gene and Apply EBayes
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Create a named list of top tables (DEGs per strain)
top_tables_all_strains <- lapply(colnames(contrast_matrix), function(contrast_name) {
  topTable(fit2, coef = contrast_name, adjust = "BH", sort.by = "P", number = Inf)
})
names(top_tables_all_strains) <- colnames(contrast_matrix)

top_tables_all_strains1 <- top_tables_all_strains

# Merge them into one named list
top_tables_all_strains_combined <- c(top_tables_all_strains1, top_tables_all_strains2, top_tables_all_strains3, top_tables_all_strains4)
save(top_tables_all_strains_combined, file=paste0(path.results,"top_tables_all_strains_combined.RData"))
#load(file=paste0(path.results,"top_tables_all_strains_combined.RData"))

#1
top_tables_filt <- lapply(top_tables_all_strains_combined, function(tbl) {
  top_tables_filt <- tbl[tbl$adj.P.Val < 0.1, ]
  top_tables_filt <- top_tables_filt[order(abs(top_tables_filt$logFC), decreasing = TRUE), ]
  return(top_tables_filt)
})

for (name in names(top_tables_filt)) {
  write.csv(top_tables_filt[[name]], file = paste0(path.results, "filt_adjpval0.1/", name, "_DEGs.csv"), row.names = TRUE)
}

#2
top_tables_filt2 <- lapply(top_tables_all_strains_combined, function(tbl) {
  top_tables_filt <- tbl[tbl$P.Value < 0.05, ]
  top_tables_filt <- top_tables_filt[order(abs(top_tables_filt$logFC), decreasing = TRUE), ]
  return(top_tables_filt)
})
for (name in names(top_tables_filt2)) {
  write.csv(top_tables_filt2[[name]], file = paste0(path.results, "filt_pval0.05/", name, "_DEGs.csv"), row.names = TRUE)
}

# Plot: TESTING -----------------------------------------------------------------
hist(top_tables_all_strains[["Reward_vs_NonReward_B6"]]$P.Value, breaks = 1000, main = "P-values - B6", xlab = "P-value", xlim = range(0.01,1.00))

LINE = "B6"
plotMDS(v$E[, sample_info$Mouseline == LINE], labels = sample_info$Condition[sample_info$Mouseline == LINE], col = ifelse(sample_info$Condition == "Reward", "tomato", "steelblue"), main = paste0("MDS - ", LINE))

top_nod <- top_tables_all_strains[["Reward_vs_NonReward_B6"]]
subset_nod <- top_nod[top_nod$P.Value < 0.05, ]
nrow(subset_nod)


# Step 1: Filter significant genes per mouseline (adjusted p-value < 0.01)
sig_genes_list <- lapply(top_tables_all_strains, function(df) {
  df %>%
    filter(adj.P.Val < 0.05) %>%
    mutate(GeneID = rownames(.)) %>%
    select(GeneID, logFC)
})

# Step 2: Get union of all significant genes across mouselines
all_sig_genes <- unique(unlist(lapply(sig_genes_list, function(df) df$GeneID)))

# Step 3: Create matrix of logFC
logfc_matrix <- sapply(sig_genes_list, function(df) {
  vec <- setNames(df$logFC, df$GeneID)
  logfc_values <- vec[all_sig_genes]
  logfc_values[is.na(logfc_values)] <- 0
  return(logfc_values)
})

rownames(logfc_matrix) <- all_sig_genes

# Step 4: Plot heatmap
pheatmap(logfc_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Log Fold Change for Significant Genes Across Mouselines",
         fontsize_row = 2)

# Ensembl Annotation ------------------------------------------------------

#Connect to Ensembl
library(org.Mm.eg.db)

top_tables_with_symbols <- lapply(top_tables_filt, function(df) {
  if (nrow(df) > 0) {
    gene_symbols <- as.data.frame(mapIds(org.Mm.eg.db, 
                                         keys = rownames(df), 
                                         column = "SYMBOL", 
                                         keytype = "ENSEMBL", 
                                         multiVals = "first"))
    colnames(gene_symbols) <- "Symbol"
    
    df$Symbol <- gene_symbols[rownames(df), "Symbol"]
    rownames(df) <- paste(rownames(df), gene_symbols[rownames(df), "Symbol"], sep = " - ")
    df <- df[order(abs(df$logFC), decreasing = TRUE), ]
      }
  return(df)
})
names(top_tables_with_symbols) <- names(top_tables_filt)

#list2env(top_tables_with_symbols, envir = .GlobalEnv)

for (name in names(top_tables_with_symbols)) {
  write.csv(top_tables_with_symbols[[name]], 
            file = paste0(path.results,name, ".csv"), 
            row.names = TRUE)
}


# Plot: DEG counts per strain / Figure 2 ---------------------------------------------

deg_counts <- sapply(top_tables_all_strains_combined, function(tbl) {
  sum(tbl$adj.P.Val < 0.1)
})

deg_df <- data.frame(
  Groups = names(deg_counts),
  Num_DEGs = as.integer(deg_counts)
)

deg_df[order(-deg_df$Num_DEGs), ]

#deg_df$Groups <- gsub("Reward_vs_NonReward_", "", deg_df$Groups)

png(paste0(path.results,"Figure2.png"), width = 4000, height = 1500,  res = 300) 
ggplot(deg_df, aes(x = reorder(Groups, -Num_DEGs), y = Num_DEGs)) +
  geom_bar(stat = "identity", fill = "#2C3E50") +
  labs(title = "Number of DEGs (Adj.P-Value < 0.1) per Group",
       x = "Groups",
       y = "Number of DEGs") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 60))
dev.off()




# Plot: Heatmap -----------------------------------------------------------------

for (name in names(top_tables_filt)) {
  if (nrow(top_tables_filt[[name]]) == 0) next 
  
  strain <- sub("Reward_vs_NonReward_", "", name)
  strain_samples <- sample_info[sample_info$Mouseline == strain, 10]
  strain_expr <- normalized_counts[, strain_samples]
  strain_annot <- sample_info[sample_info$Mouseline == strain, c("Condition", "Mouseline")]
  
  top_genes <- rownames(top_tables_filt[[name]])[1:5]
  top_genes <- top_genes[top_genes %in% rownames(strain_expr)]
  if (length(top_genes) == 0) next
  
  # Convert annotation columns to factors
  strain_annot$Condition <- factor(strain_annot$Condition)
  strain_annot$Mouseline <- factor(strain_annot$Mouseline)
  
  # Set annotation colors
  annot_colors <- list(
    Condition = setNames(c("#009E73", "#D55E00"), levels(strain_annot$Condition)),
    Mouseline = setNames(rainbow(length(unique(strain_annot$Mouseline))), unique(strain_annot$Mouseline))
  )
  
  png(paste0(path.results, name, "_heatmap.png"), width = 6000, height = 6000, res = 600)
  pheatmap(log2(cpm(strain_expr)[top_genes, ] + 1),
           annotation_col = strain_annot,
           annotation_colors = annot_colors,
           fontsize_row = 6, fontsize_col = 6,
           angle_col = 45,
           clustering_method = "ward.D2",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "euclidean",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           main = paste("Top DEGs -", name))
  dev.off()
}


# Plot: Dotplot / Figure 3 -----------------------------------------------------------------

#1
plot_data <- list()

for (name in names(top_tables_filt)) {
  if (nrow(top_tables_filt[[name]]) == 0) next
  
  group <- name
  tbl <- top_tables_filt[[group]]
  
  df_tmp <- data.frame(
    ENSEMBL = rownames(tbl),
    Group = group,
    pvalue = tbl$adj.P.Val,
    abs_log = abs(tbl$logFC)
  )
  plot_data[[name]] <- df_tmp
}

df <- do.call(rbind, plot_data)

df_filtered <- df %>%
  distinct(ENSEMBL, Group) %>%        
  count(ENSEMBL) %>%
  filter(n >= 2) %>%
  inner_join(df, by = "ENSEMBL")

png(paste0(path.results,"Figure3.png"), width = 6000, height = 5000,  res = 200) 
is_ggplot <- function(x) inherits(x, "gg")
ggplot(df_filtered, aes(x = Group, y = ENSEMBL)) +
  geom_point(aes(size = abs_log, color = pvalue)) +
  scale_color_viridis_c(option = "viridis", direction = 1, name = "Adj.P.Value") +
  scale_size_continuous(range = c(2, 8), name = "abs_log") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Striatal genes influenced by Social Rewards\n in Groups (adj.p-val < 0.1; min 2 groups)",
    x = "Groups",
    y = "ENSEMBL"
  )
dev.off()

#2
plot_data <- list()

for (name in names(top_tables_filt2)) {
  if (nrow(top_tables_filt2[[name]]) == 0) next
  
  strain <- sub("Reward_vs_NonReward_", "", name)
  tbl <- top_tables_filt2[[name]]
  
  df_tmp <- data.frame(
    ENSEMBL = rownames(tbl),
    Strain = strain,
    pvalue = tbl$P.Value,
    abs_log = abs(tbl$logFC)
  )
  plot_data[[name]] <- df_tmp
}

df <- do.call(rbind, plot_data)

is_ggplot <- function(x) inherits(x, "gg")
ggplot(df, aes(x = Strain, y = ENSEMBL)) +
  geom_point(aes(size = abs_log, color = pvalue)) +
  scale_color_viridis_c(option = "viridis", direction = 1, name = "p-value") +
  scale_size_continuous(range = c(2, 8), name = "abs_log") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Striatal genes influenced by Social Rewards\n in CC/founder strains (p-val < 0.0001)",
    x = "Strain",
    y = "ENSEMBL"
  )


# Plot: Multi-panel volcano plot / Figure 4 ------------------------------------------------

#Groups
logFC_thresh <- 0.75
pval_thresh <- 0.1

plot_list <- list()
for (name in names(top_tables_filt)) {
  tbl <- top_tables_filt[[name]]
  if (nrow(tbl) == 0) next
  
  groups <- name
  
  df_tmp <- data.frame(
    ENSEMBL = rownames(tbl),
    logFC = tbl$logFC,
    pval = tbl$adj.P.Val,
    log10pval = -log10(tbl$adj.P.Val),
    groups = groups
  )
  
  df_tmp <- df_tmp %>%
    mutate(
      group = case_when(
        logFC >= logFC_thresh & pval < pval_thresh ~ "Up",
        logFC <= -logFC_thresh & pval < pval_thresh ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  plot_list[[name]] <- df_tmp
}

volcano_df <- bind_rows(plot_list)

volcano_df$groups <- factor(volcano_df$groups)

png(paste0(path.results,"Figure4.png"), width = 5000, height = 4000,  res = 300) 
ggplot(volcano_df, aes(x = logFC, y = log10pval)) +
  geom_point(aes(color = group), size = 0.8, alpha = 0.7) +
  facet_wrap(~ groups, scales = "free") +
  geom_hline(yintercept = -log10(pval_thresh), color = "red", linetype = "solid", linewidth = 0.3) +
  geom_vline(xintercept = c(-logFC_thresh, logFC_thresh), color = "red", linetype = "solid", linewidth = 0.3) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
  xlim(-3, 3) +
  ylim(0, 6) +
  labs(
    title = "Strain dependent variation of gene expression as a result of social rewards in the Striatum \n (Adj.P.Value < 0.1)",
    x = "log2FoldChange",
    y = "-log10(P Value)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

#Mouseline
logFC_thresh <- 0.75
pval_thresh <- 0.05

plot_list <- list()
for (name in names(top_tables_all_strains)) {
  tbl <- top_tables_all_strains[[name]]
  if (nrow(tbl) == 0) next
  
  strain <- sub("Reward_vs_NonReward_", "", name)
  
  df_tmp <- data.frame(
    ENSEMBL = rownames(tbl),
    logFC = tbl$logFC,
    pval = tbl$P.Value,
    log10pval = -log10(tbl$P.Value),
    strain = strain
  )
  
  df_tmp <- df_tmp %>%
    mutate(
      group = case_when(
        logFC >= logFC_thresh & pval < pval_thresh ~ "Up",
        logFC <= -logFC_thresh & pval < pval_thresh ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  plot_list[[name]] <- df_tmp
}

volcano_df <- bind_rows(plot_list)

volcano_df$strain <- factor(
  volcano_df$strain,
  levels = c("AJ", "B6", "CC002", "CC005", "CC019", "CC051", "NOD", "PWK", "DBA")
)

ggplot(volcano_df, aes(x = logFC, y = log10pval)) +
  geom_point(aes(color = group), size = 0.8, alpha = 0.7) +
  facet_wrap(~ strain, scales = "free") +
  geom_hline(yintercept = -log10(pval_thresh), color = "red", linetype = "solid", linewidth = 0.3) +
  geom_vline(xintercept = c(-logFC_thresh, logFC_thresh), color = "red", linetype = "solid", linewidth = 0.3) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
  xlim(-3, 3) +
  ylim(0, 6) +
  labs(
    title = "Strain dependent variation of gene expression\nas a result of housing condition in the hippocampus",
    x = "log2FoldChange",
    y = "-log10(P Value)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )


# Plot: Barplot DEG per pairwise interaction / Figure 5 ------------------------------------

#Groups
pval_thresh <- 0.1
pairwise_deg_counts <- list()

pairwise_names <- names(top_tables_filt)

for (name in pairwise_names) {
  tbl <- top_tables_filt[[name]]
  degs <- rownames(tbl[tbl$adj.P.Val < pval_thresh, ])
  pairwise_deg_counts[[name]] <- length(degs)
}

deg_counts_df <- data.frame(
  Interaction = names(pairwise_deg_counts),
  DEG_Count = unlist(pairwise_deg_counts)
)

deg_counts_df <- deg_counts_df[order(-deg_counts_df$DEG_Count), ]

png(paste0(path.results,"Figure5.png"), width = 4000, height = 3000,  res = 300) 
ggplot(deg_counts_df, aes(x = reorder(Interaction, -DEG_Count), y = DEG_Count)) +
  geom_bar(stat = "identity", fill = "gray30") +
  labs(
    title = "Number of genes differentially expressed when comparing the interaction between two \n groups in Striatum (Adj.P-value < 0.1)",
    x = "Sex interaction",
    y = "Count of shared DEGs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

#Mouseline
pval_thresh <- 0.05
pairwise_deg_counts <- list()

pairwise_names <- grep("_vs_", names(top_tables_all_strains), value = TRUE)
pairwise_names <- pairwise_names[!grepl("Reward_vs_NonReward", pairwise_names)] 

for (name in pairwise_names) {
  tbl <- top_tables_all_strains[[name]]
  
  if (is.null(tbl) || nrow(tbl) == 0) next
  degs <- rownames(tbl[tbl$adj.P.Val < pval_thresh, ])
  
  if (length(degs) == 0) next
  
  pairwise_deg_counts[[name]] <- length(degs)
}

deg_counts_df <- data.frame(
  Interaction = names(pairwise_deg_counts),
  DEG_Count = unlist(pairwise_deg_counts)
)

deg_counts_df <- deg_counts_df[order(-deg_counts_df$DEG_Count), ]

ggplot(deg_counts_df, aes(x = reorder(Interaction, -DEG_Count), y = DEG_Count)) +
  geom_bar(stat = "identity", fill = "gray30") +
  labs(
    title = "Number of genes differentially expressed when \n comparing the interaction between two \n strains in Striatum (Adj.P-value < 0.05)",
    x = "Strain X Strain interaction",
    y = "Count of shared DEGs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


# Plot: Expression per gene  / Figure 6 -----------------------------------------------------

df_unique <- df_filtered %>%
  arrange(desc(n), pvalue) %>%
  distinct(ENSEMBL, .keep_all = TRUE) 

# Groups
pval_thresh <- 0.1
table <- top_tables_filt[["Female_vs_Male_Reward"]]
deg <- table[table$adj.P.Val < pval_thresh, ]
ensmusg_ids <- rownames(deg)
length(ensmusg_ids)

gene_id <- ensmusg_ids[1]
gene_symbol <- as.data.frame(mapIds(org.Mm.eg.db, 
                                     keys = gene_id, 
                                     column = "SYMBOL", 
                                     keytype = "ENSEMBL", 
                                     multiVals = "first"))

groups <- unique(sample_info$Mouseline)

summarize_gene_counts <- function(line) {
  sample_subset <- sample_info %>% filter(Mouseline == line)
  dge_sub <- dge[gene_id, sample_subset$Sample_ID]
  gene_counts_vector <- dge_sub$counts[gene_id, ]
  gene_counts_df <- data.frame(
    Sample_ID = names(gene_counts_vector),
    Counts = as.numeric(gene_counts_vector)
  )
  gene_counts_merged <- gene_counts_df %>%
    left_join(sample_subset, by = "Sample_ID")
  
  gene_counts_merged %>%
    group_by(Mouseline, Condition) %>%
    summarise(
      Total_Counts = sum(Counts, na.rm = TRUE),
      Mean_Counts = mean(Counts, na.rm = TRUE),
      SD_Counts = sd(Counts, na.rm = TRUE),
      SE_Counts = sd(Counts, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
}

all_counts_summary <- map_dfr(groups, summarize_gene_counts)
print(all_counts_summary)

png(paste0(path.results,paste0("Figure6_Female_vs_Male_Reward_",gene_id,".png")), width = 3000, height = 2000,  res = 250) 
ggplot(all_counts_summary, aes(x = Condition, y = Mean_Counts)) +
  geom_point(aes(color = Condition), position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = Mean_Counts - SE_Counts, ymax = Mean_Counts + SE_Counts, color = Condition),
                width = 0.2, position = position_dodge(width = 0.3)) +
  geom_line(aes(group = Mouseline), color = "gray70", size = 1) +
  facet_wrap(~ Mouseline, nrow = 1, strip.position = "bottom") +
  scale_color_manual(values = c("Reward" = "#F8766D", "Nonreward" = "#00BFC4"),
                     name = "Social Rewards") +
  labs(
    title = paste0("Expression levels of ", gene_symbol[[1]], " (", gene_id, ") \nby Social Reward Status with a significant\ninteraction across Mouseline (Adj.P-value < 0.1)"),
    x = "Mouseline",
    y = "Mean Counts"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  )
dev.off()

# Plot: Pathway and Enrichment Analysis / Figure 7 -----------------------------------------

deg_df[order(-deg_df$Num_DEGs), ]

pval_thresh <- 0.1
table <- top_tables_filt[["Female_vs_Male_Reward"]]
deg <- table[table$adj.P.Val < pval_thresh, ]
ensmusg_ids <- rownames(deg)

gene_symbols <- as.data.frame(mapIds(org.Mm.eg.db, 
                                         keys = ensmusg_ids, 
                                         column = "SYMBOL", 
                                         keytype = "ENSEMBL", 
                                         multiVals = "first"))
colnames(gene_symbols) <- "Symbol"
ensmusg_ids.df <- as.data.frame(ensmusg_ids)
ensmusg_ids.df$Symbol <- gene_symbols[ensmusg_ids.df$ensmusg_ids, "Symbol"]

  if (nrow(ensmusg_ids.df) >= 1) {
    ego.df <- as.data.frame(enrichGO(
      gene = ensmusg_ids.df$Symbol,
      OrgDb = org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pvalueCutoff = 0.2,
      readable = TRUE
    ))
  } else {
    NULL
  }

top_terms <- ego.df %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)

png(paste0(path.results,"Figure7C.png"), width = 4000, height = 3000,  res = 300) 
ggplot(top_terms, aes(x = Count, y = reorder(Description, Count))) +
  geom_bar(stat = "identity", fill = "#00BFC4") +
  labs(
    title = "GO Enriched for all DEGs found for Sex (Female_vs_Male_Reward) \n interaction in the striatum (Adj.P-value < 0.1)",
    x = "Gene Enriched For Term",
    y = "GO Term"
  ) +
  theme_minimal(base_size = 14)
dev.off()

# Heatmap top 50 genes CC019 (Reward vs Nonreward Female)  / Figure 8 -----------------

tbl_cc019 <- top_tables_all_strains_combined[["Reward_vs_Nonreward_Female_CC019"]]
top_cc019 <- tbl_cc019[order(tbl_cc019$adj.P.Val), ][1:50, ]
genes_top <- rownames(top_cc019)

mat <- normalized_counts[genes_top, ]
meta <- sample_info 

# CC019
mat_cc019 <- mat[, meta$Mouseline == "CC019" & meta$Sex == "Female"]
meta_cc019 <- meta[meta$Mouseline == "CC019" & meta$Sex == "Female", ]

ann <- data.frame(
  Sex = meta_cc019$Sex,
  Condition = meta_cc019$Condition
)
rownames(ann) <- colnames(mat_cc019)

png(paste0(path.results,"Figure8.png"), width = 6000, height = 5000,  res = 200) 
pheatmap(
  mat_cc019,
  annotation_col = ann,
  show_rownames = TRUE,
  show_colnames = FALSE,
  scale = "row",   # z-score por gene
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top 50 DEGs in CC019 (Female Reward vs Nonreward)"
)
dev.off()


# Pathway enrichment DEGs CC019  / Figure 9 -------------------------------------------

deg_cc019 <- rownames(tbl_cc019[tbl_cc019$adj.P.Val < 0.1, ])

gene_symbols <- as.data.frame(mapIds(org.Mm.eg.db, 
                                     keys = deg_cc019, 
                                     column = "SYMBOL", 
                                     keytype = "ENSEMBL", 
                                     multiVals = "first"))
colnames(gene_symbols) <- "Symbol"
deg_cc019.df <- as.data.frame(deg_cc019)
deg_cc019.df$Symbol <- gene_symbols[deg_cc019.df$deg_cc019, "Symbol"]

ego <- enrichGO(
    gene = deg_cc019.df$Symbol,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod= "BH",
    pvalueCutoff = 0.2,
    readable = TRUE
  )
  
png(paste0(path.results,"Figure9.png"), width = 3500, height = 4000,  res = 250) 
dotplot(ego, showCategory = 30) + ggtitle("CC019 Female Reward vs Nonreward")
dev.off()



# Top 10 genes em boxplots / Figure 10 ------------------------------------------------

top10 <- rownames(tbl_cc019[order(tbl_cc019$adj.P.Val), ])[1:10]

for (g in top10) {
  df <- data.frame(
    expr = normalized_counts[g, ],
    meta
  )
  png(filename = file.path(path.results, paste0("Expression_of_", g, ".png")),width = 4000, height = 3000, res = 300)
  p <- ggplot(df, aes(x = interaction(Sex, Condition), y = expr, fill = Condition)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    facet_wrap(~Mouseline) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste("Expression of", g))
  print(p)
  dev.off()
}


# # Top 10 genes em violinplot / Figure 11 --------------------------------

top10 <- rownames(tbl_cc019[order(tbl_cc019$adj.P.Val), ])[1:10]

for (g in top10) {
  df <- data.frame(
    expr = normalized_counts[g, ],
    meta
  )
  
  png(filename = file.path(path.results, paste0("Violin plot of", g, ".png")),width = 4000, height = 3000, res = 300)
  p <- ggplot(df, aes(x = Sex, y = expr, fill = Condition)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    facet_wrap(~Mouseline) +
    theme_bw() +
    ggtitle(paste("Violin plot of", g))
  print(p)
  dev.off()
}



# Plot: hallmark + Cannonical pathways + cell type signature --------------------

# Output directory
dir.create(path.results, showWarnings = FALSE)

# Load MSigDB hallmark and canonical pathways for mouse
msig_hallmark <- msigdbr(species = "Mus musculus", category = "H")
msig_canonical <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")

# Format MSigDB gene sets
hallmark_terms <- msig_hallmark[, c("gs_name", "gene_symbol")]
canonical_terms <- msig_canonical[, c("gs_name", "gene_symbol")]

# Enrichr cell-type signature database
enrichr_db <- c("Mouse_Cell_Type_Signatures_from_PanglaoDB_2021")

# Main enrichment loop
extended_enrichment <- list()

for (strain_name in names(top_tables_with_symbols)[-1]) {
  cat("Running enrichment for:", strain_name, "\n")
  
  df <- top_tables_with_symbols[[strain_name]]
  sig_genes <- df$Symbol[df$adj.P.Val < 0.05 & !is.na(df$Symbol)]
  
  if (length(sig_genes) >= 5) {
    enrichments <- list()
    
    # GO:BP
    enrichments$GO_BP <- tryCatch({
      enrichGO(
        gene = sig_genes,
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pvalueCutoff = 0.2,
        readable = TRUE
      )
    }, error = function(e) NULL)
    
    # Hallmark
    enrichments$Hallmark <- tryCatch({
      enricher(
        gene = sig_genes,
        TERM2GENE = hallmark_terms,
        pvalueCutoff = 0.2
      )
    }, error = function(e) NULL)
    
    # Canonical Pathways
    enrichments$Canonical <- tryCatch({
      enricher(
        gene = sig_genes,
        TERM2GENE = canonical_terms,
        pvalueCutoff = 0.2
      )
    }, error = function(e) NULL)
    
    # Cell Types (via Enrichr)
    enrichments$CellTypes <- tryCatch({
      enrichr(sig_genes, enrichr_db)[[1]]
    }, error = function(e) NULL)
    
    # Save to list
    extended_enrichment[[strain_name]] <- enrichments
    
    # Save CSVs
    if (!is.null(enrichments$GO_BP)) 
      write.csv(enrichments$GO_BP@result, paste0(path.results, strain_name, "_GO_BP.csv"), row.names = FALSE)
    if (!is.null(enrichments$Hallmark)) 
      write.csv(enrichments$Hallmark@result, paste0(path.results, strain_name, "_Hallmark.csv"), row.names = FALSE)
    if (!is.null(enrichments$Canonical)) 
      write.csv(enrichments$Canonical@result, paste0(path.results, strain_name, "_Canonical.csv"), row.names = FALSE)
    if (!is.null(enrichments$CellTypes)) 
      write.csv(enrichments$CellTypes, paste0(path.results, strain_name, "_CellTypes.csv"), row.names = FALSE)
    
    # Save dotplots
    png(paste0(path.results, strain_name, "_GO_BP.png"), width = 6000, height = 6000, res = 600)
    if (!is.null(enrichments$GO_BP) && nrow(enrichments$GO_BP) > 0) {
      print(dotplot(enrichments$GO_BP, showCategory = 20) + ggtitle(paste(strain_name, "GO:BP")))
    }
    dev.off()
    
    png(paste0(path.results, strain_name, "_Hallmark.png"), width = 6000, height = 6000, res = 600)
    if (!is.null(enrichments$Hallmark) && nrow(enrichments$Hallmark) > 0) {
      print(dotplot(enrichments$Hallmark, showCategory = 20) + ggtitle(paste(strain_name, "Hallmark Pathways")))
    }
    dev.off()
    
    png(paste0(path.results, strain_name, "_Canonical.png"), width = 6000, height = 6000, res = 600)
    if (!is.null(enrichments$Canonical) && nrow(enrichments$Canonical) > 0) {
      print(dotplot(enrichments$Canonical, showCategory = 20) + ggtitle(paste(strain_name, "Canonical Pathways")))
    }
    dev.off()
    
  } else {
    cat("Not enough genes for enrichment in:", strain_name, "\n")
    extended_enrichment[[strain_name]] <- NULL
  }
}

###########




