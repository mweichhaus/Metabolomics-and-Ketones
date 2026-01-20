################################################################################
# Paired Differential Expression Analysis of ER+ Breast Cancer vs Normal Tissue
# Dataset: GSE58135 (NCBI GEO)
# Analysis: DESeq2 paired design for ketone body and one-carbon metabolism genes
# 
# Author: Michael Weichhaus
# Date: January 2026
# R Version: 4.5.1
################################################################################

# Install required packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "DESeq2", "org.Hs.eg.db"), update = FALSE, ask = FALSE)

# Load libraries
library(data.table)
library(GEOquery)
library(DESeq2)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
library(reshape2)

################################################################################
# 1. DOWNLOAD AND LOAD DATA
################################################################################

# Download NCBI-generated raw counts file from GEO
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE58135&format=file&file=GSE58135_raw_counts_GRCh38.p13_NCBI.tsv.gz"

temp_file <- tempfile(fileext = ".tsv.gz")
cat("Downloading raw counts file...\n")
download.file(url, temp_file, mode = "wb")

# Read the count matrix
cat("Reading count matrix...\n")
tbl <- as.matrix(fread(temp_file, header = TRUE), rownames = 1)

cat("Count table loaded successfully!\n")
cat("Dimensions:", nrow(tbl), "genes x", ncol(tbl), "samples\n")
cat("First 5 column names:", paste(head(colnames(tbl), 5), collapse = ", "), "\n")

################################################################################
# 2. FETCH METADATA FROM GEO
################################################################################

cat("\nFetching metadata from GEO...\n")
gse_full <- getGEO("GSE58135", getGPL = FALSE, AnnotGPL = FALSE)
pData <- pData(gse_full[[1]])

cat("\nSample information:\n")
cat("Total samples in metadata:", nrow(pData), "\n")
print(table(pData$source_name_ch1))

################################################################################
# 3. IDENTIFY ER+ TUMORS AND ADJACENT NORMALS
################################################################################

# Filter for ER+ breast cancer samples
er_tumor <- pData$source_name_ch1 == "ER+ Breast Cancer Primary Tumor"
er_adj <- pData$source_name_ch1 == "Uninvolved Breast Tissue Adjacent to ER+ Primary Tumor"

cat("\nER+ tumor samples:", sum(er_tumor), "\n")
cat("ER+ adjacent normal samples:", sum(er_adj), "\n")

################################################################################
# 4. MATCH TUMOR-NORMAL PAIRS BY PATIENT ID
################################################################################

# Extract patient IDs from sample titles (e.g., "01-02-A056a" -> "01-02")
pData$patient_id <- sub("-A.*$", "", pData$title)

# Get patient IDs for each group
tumor_patients <- pData$patient_id[er_tumor]
adj_patients <- pData$patient_id[er_adj]

# Find patients with BOTH tumor and adjacent normal
paired_patients <- intersect(tumor_patients, adj_patients)

cat("\nPatients with matched tumor-normal pairs:", length(paired_patients), "\n")
cat("Paired patients:", paste(head(paired_patients, 10), collapse=", "), "\n")

# Select GSMs for paired samples only
er_paired_gsm <- pData$geo_accession[pData$patient_id %in% paired_patients & (er_tumor | er_adj)]

cat("Total paired samples (tumor + normal):", length(er_paired_gsm), "\n")

################################################################################
# 5. FILTER COUNT TABLE TO PAIRED SAMPLES
################################################################################

# Check which samples are in the count matrix
available_samples <- colnames(tbl)[colnames(tbl) %in% er_paired_gsm]
cat("Samples available in count matrix:", length(available_samples), "\n")

# Reload count table and filter to paired samples
cat("\nReloading count matrix...\n")
temp_file <- tempfile(fileext = ".tsv.gz")
download.file(url, temp_file, mode = "wb")
tbl <- as.matrix(fread(temp_file, header = TRUE), rownames = 1)

# Filter to only paired samples
tbl <- tbl[, colnames(tbl) %in% er_paired_gsm]

cat("Filtered count matrix dimensions:", nrow(tbl), "x", ncol(tbl), "\n")

# Filter low-count genes (keep genes with ≥10 counts in ≥2 samples)
cat("\nFiltering low-count genes...\n")
keep <- rowSums(tbl >= 10) >= 2
tbl <- tbl[keep, ]
cat("Genes after filtering:", nrow(tbl), "\n")

################################################################################
# 6. CREATE SAMPLE METADATA (colData)
################################################################################

cat("\nCreating sample metadata...\n")
coldata <- data.frame(
  row.names = colnames(tbl),
  condition = ifelse(
    pData$source_name_ch1[match(colnames(tbl), pData$geo_accession)] == "ER+ Breast Cancer Primary Tumor",
    "tumor", "normal"
  ),
  patient = pData$patient_id[match(colnames(tbl), pData$geo_accession)]
)
coldata$condition <- factor(coldata$condition, levels = c("normal", "tumor"))

cat("Initial sample distribution:\n")
print(table(coldata$condition))
cat("Number of patients:", length(unique(coldata$patient)), "\n")

################################################################################
# 7. APPLY STRICT 1:1 MATCHING
################################################################################

cat("\n========== FILTERING TO STRICT 1:1 PAIRED SAMPLES ==========\n")

# Count tumors and normals per patient
patient_tumor_counts <- aggregate(condition == "tumor" ~ patient, data = coldata, FUN = sum)
patient_normal_counts <- aggregate(condition == "normal" ~ patient, data = coldata, FUN = sum)

colnames(patient_tumor_counts) <- c("patient", "n_tumor")
colnames(patient_normal_counts) <- c("patient", "n_normal")

patient_summary <- merge(patient_tumor_counts, patient_normal_counts, by = "patient")

# Keep only patients with exactly 1 tumor AND 1 normal
valid_patients <- patient_summary$patient[
  patient_summary$n_tumor == 1 & patient_summary$n_normal == 1
]

cat("Patients with exactly 1 tumor and 1 normal:", length(valid_patients), "\n")
cat("Valid patients:", paste(valid_patients, collapse = ", "), "\n")

# Filter coldata and count matrix
coldata <- coldata[coldata$patient %in% valid_patients, ]
tbl <- tbl[, rownames(coldata)]

cat("\n========== FINAL PAIRED DATASET ==========\n")
cat("Total samples:", nrow(coldata), "\n")
print(table(coldata$condition))
cat("Number of patients:", length(unique(coldata$patient)), "\n")

# Verify perfect pairing
cat("\nVerifying perfect 1:1 pairing:\n")
patient_counts <- table(coldata$patient)
print(table(patient_counts))  # Should show all patients have exactly 2 samples

################################################################################
# 8. MAP GENE SYMBOLS TO ENTREZ IDs
################################################################################

cat("\nMapping genes of interest...\n")
genes_of_interest <- c("SHMT1", "SHMT2", "BCAT1", "OXCT1", "OXCT2", 
                       "ACAT1", "ACAT2", "BDH1", "BDH2")
entrez_map <- select(org.Hs.eg.db, keys = genes_of_interest, 
                     columns = "ENTREZID", keytype = "SYMBOL")
entrez_ids <- entrez_map$ENTREZID

cat("Genes found:", sum(entrez_ids %in% rownames(tbl)), "/", 
    length(genes_of_interest), "\n")

################################################################################
# 9. RUN DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS
################################################################################

cat("\nRunning DESeq2 analysis (paired design)...\n")
dds <- DESeqDataSetFromMatrix(countData = tbl, 
                               colData = coldata, 
                               design = ~patient + condition)
dds <- DESeq(dds)

# Get results (tumor vs normal)
res <- results(dds, contrast = c("condition", "tumor", "normal"))

# Extract results for genes of interest
res_genes <- as.data.frame(res[rownames(res) %in% entrez_ids, ])
res_genes$symbol <- entrez_map$SYMBOL[match(rownames(res_genes), entrez_map$ENTREZID)]
res_genes <- res_genes[order(res_genes$symbol), ]

################################################################################
# 10. PRINT RESULTS
################################################################################

cat("\n========== DIFFERENTIAL EXPRESSION RESULTS ==========\n")
cat("log2FoldChange > 0 means higher expression in tumor\n")
cat("log2FoldChange < 0 means higher expression in normal\n\n")
print(res_genes[, c("symbol", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")])

################################################################################
# 11. PREPARE DATA FOR VISUALIZATION
################################################################################

cat("\nGenerating visualizations...\n")

# Get VST-normalized expression
vsd <- vst(dds, blind = FALSE)
expr <- assay(vsd)[rownames(vsd) %in% entrez_ids, ]
rownames(expr) <- entrez_map$SYMBOL[match(rownames(expr), entrez_map$ENTREZID)]

# Prepare data frame for plotting
expr_df <- melt(t(expr))
colnames(expr_df) <- c("sample", "gene", "expression")
expr_df$condition <- coldata$condition[match(expr_df$sample, rownames(coldata))]
expr_df$patient <- coldata$patient[match(expr_df$sample, rownames(coldata))]

# Set gene order for plots
gene_order <- c("BDH1", "BDH2", "OXCT1", "OXCT2", "ACAT1", "ACAT2", 
                "SHMT1", "SHMT2", "BCAT1")
gene_order <- gene_order[gene_order %in% unique(expr_df$gene)]
expr_df$gene <- factor(expr_df$gene, levels = gene_order)

# Create significance labels
sig_df <- data.frame(
  gene = res_genes$symbol,
  padj = res_genes$padj
)

# Tiered significance labels
sig_df$sig_label <- ifelse(sig_df$padj < 0.001, "***",
                            ifelse(sig_df$padj < 0.01, "**",
                                   ifelse(sig_df$padj < 0.05, "*", "")))

sig_df$gene <- factor(sig_df$gene, levels = gene_order)

# Calculate position for significance labels in paired plot
expr_max <- aggregate(expression ~ gene, data = expr_df, FUN = max)
sig_df <- merge(sig_df, expr_max, by = "gene")
sig_df$y_pos <- sig_df$expression * 1.05

################################################################################
# 12. GENERATE VISUALIZATIONS
################################################################################

# Plot 1: Paired line plot with significance
p1 <- ggplot(expr_df, aes(x = condition, y = expression, group = patient, color = patient)) +
  geom_point(size = 1.5) + geom_line(alpha = 0.5) +
  geom_text(data = sig_df, aes(x = 1.5, y = y_pos, label = sig_label), 
            inherit.aes = FALSE, size = 6, fontface = "bold") +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  theme_minimal() + 
  labs(title = "Paired Expression Changes (Tumor vs. Normal)", 
       y = "VST-Normalized Expression") +
  theme(legend.position = "none")
print(p1)

# Plot 2: Boxplot per gene
p2 <- ggplot(expr_df, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  theme_minimal() + 
  labs(title = "Expression Distribution (Tumor vs. Normal)", 
       y = "VST-Normalized Expression") +
  scale_fill_manual(values = c("normal" = "blue", "tumor" = "red"))
print(p2)

# Plot 3: Heatmap
expr_ordered <- expr[gene_order[gene_order %in% rownames(expr)], ]

pheatmap(expr_ordered, 
         annotation_col = coldata[, "condition", drop = FALSE],
         show_colnames = FALSE, 
         main = "Heatmap of Gene Expression (VST-Normalized)",
         clustering_method = "ward.D2", 
         scale = "row",
         cluster_rows = FALSE)

# Plot 4: Fold change barplot with tiered significance
fc_df <- data.frame(
  gene = res_genes$symbol, 
  log2FC = res_genes$log2FoldChange,
  se = res_genes$lfcSE, 
  padj = res_genes$padj
)

fc_df$sig <- ifelse(fc_df$padj < 0.001, "***",
                    ifelse(fc_df$padj < 0.01, "**",
                           ifelse(fc_df$padj < 0.05, "*", "")))

p4 <- ggplot(fc_df, aes(x = reorder(gene, log2FC), y = log2FC, fill = log2FC > 0)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = log2FC - se, ymax = log2FC + se), width = 0.2) +
  geom_text(aes(label = sig), hjust = -0.5, size = 6, fontface = "bold") +
  theme_minimal() + 
  labs(title = "Log2 Fold Change (Tumor vs. Normal)", 
       y = "Log2 Fold Change", 
       x = "Gene") +
  scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red"), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip()
print(p4)

cat("\n========== ANALYSIS COMPLETE ==========\n")

################################################################################
# Optional: Save results and plots
################################################################################

# Uncomment to save results
# write.csv(res_genes, "DESeq2_results_ER_breast_cancer.csv", row.names = FALSE)
# ggsave("Figure1_paired_expression.pdf", p1, width = 10, height = 8)
# ggsave("Figure2_fold_change.pdf", p4, width = 6, height = 6)

################################################################################
# END OF SCRIPT
################################################################################