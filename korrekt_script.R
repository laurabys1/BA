#SCRIPT FOR BACHELOR THESIS

#LOADING OF PROTEOMICS DATA----
setwd("~/Desktop")
df <- read.csv("LASTFILE.csv", stringsAsFactors = FALSE)

#Removing non-significant values (ANOVA)
filtered_df <- df[df$ANOVA.Significant == "+", ]

#To account for inverted sign in Log2 Fold Change Columns, those are multiplied by -1
filtered_df[, c("log2_FC_Ctrl_vs_2h4ME", "log2_FC_Ctrl_vs_20HLPS", "log2_FC_Ctrl_vs_2H4ME20HLPS", "log2_FC_Ctrl_vs_22H4ME20HLPS", "log2_FC_Ctrl_vs_22H4ME")] <- filtered_df[, c("log2_FC_Ctrl_vs_2h4ME", "log2_FC_Ctrl_vs_20HLPS", "log2_FC_Ctrl_vs_2H4ME20HLPS", "log2_FC_Ctrl_vs_22H4ME20HLPS", "log2_FC_Ctrl_vs_22H4ME")] * -1

#LOADING OF TRANSCRIPTOMICS DATA----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # For human genes; replace for other organisms if needed
BiocManager::install("DOSE") 

#Loading everything correctly
library(readr)
RawDat <- read_tsv("/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/Laura_JB_caco2_RNA-seq.tsv", col_names = TRUE)
counts <- RawDat[,-c(1,2)]
gene_info <- RawDat[,c(1:2)]
counts <- as.data.frame(counts)
coldata <-read_tsv("/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/Laura_JB_caco2_RNA-seq_meta.tsv", col_names = TRUE)
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample
coldata$group <- as.factor(coldata$group)
coldata$replicate <- as.factor(coldata$replicate)

#Samples 
files <- c("I223" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I223/quant.sf",
           "I224" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I224/quant.sf",
           "I225" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I225/quant.sf",
           "I226" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I226/quant.sf",
           "I227" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I227/quant.sf",
           "I228" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I228/quant.sf",
           "I229" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I229/quant.sf",
           "I230" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I230/quant.sf",
           "I231" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I231/quant.sf",
           "I232" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I232/quant.sf",
           "I233" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I233/quant.sf",
           "I234" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I234/quant.sf",
           "I235" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I235/quant.sf",
           "I236" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I236/quant.sf",
           "I237" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I237/quant.sf",
           "I238" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I238/quant.sf",
           "I239" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I239/quant.sf",
           "I240" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I240/quant.sf",
           "I241" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I241/quant.sf",
           "I242" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I242/quant.sf",
           "I243" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I243/quant.sf",
           "I244" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I244/quant.sf",
           "I245" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I245/quant.sf",
           "I246" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I246/quant.sf")


#Data import, transcript level quantification using Salmon
txi <- tximport(files,
                type = "salmon",
                txOut = TRUE,
                countsFromAbundance = "lengthScaledTPM")


#Round to integer for DESeq2 analysis
counts_int <- round(txi$counts)

#DESeq2 using groups as experimental design
dds2 <- DESeqDataSetFromMatrix(countData = counts_int,
                               colData = coldata,
                               design = ~ group)

#Keeps only genes with sample counts over >= 10
keep <- apply(counts(dds2), 1, function(x) all(x >= 10))
dds2_filtered <- dds2[keep, ]

#DESeq2, size factor estimation, dispersion estimation and GLM fitting.
dds2 <- DESeq(dds2_filtered)

#Log transformation
rld <- rlog(dds2, blind = FALSE)

#PCA plot using Group factor
pcaPlot <- plotPCA(rld, intgroup = "group")
pcaPlot

#3 bad replicates with extreme outliers are observed, these files need to be removed
tx2gene <- RawDat[,c(1,2)]

files <- c("I223" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I223/quant.sf",
           "I224" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I224/quant.sf",
           "I225" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I225/quant.sf",
           "I226" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I226/quant.sf",
           "I227" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I227/quant.sf",
           "I229" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I229/quant.sf",
           "I230" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I230/quant.sf",
           "I231" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I231/quant.sf",
           "I232" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I232/quant.sf",
           "I233" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I233/quant.sf",
           "I234" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I234/quant.sf",
           "I235" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I235/quant.sf",
           "I237" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I237/quant.sf",
           "I238" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I238/quant.sf",
           "I239" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I239/quant.sf",
           "I240" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I240/quant.sf",
           "I241" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I241/quant.sf",
           "I243" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I243/quant.sf",
           "I244" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I244/quant.sf",
           "I245" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I245/quant.sf",
           "I246" = "/Users/laurabyskov/Desktop/Transcriptomics/Laura - RNA-seq/I246/quant.sf")

#Import
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

#Removal of faulty replicates from the meta data
coldata_mod <- coldata[-c(6,14,20),]
raw_counts_df <- as.data.frame(txi$counts)

#Add a column with gene identifiers
raw_counts_df$gene <- rownames(raw_counts_df)

#Reshape the data from wide to long format
raw_counts_long <- melt(raw_counts_df, id.vars = "gene", 
                        variable.name = "sample", 
                        value.name = "counts")

#Boxplot with a log2 transformation
ggplot(raw_counts_long, aes(x = sample, y = log2(counts))) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "Boxplot of Raw Counts (log2-transformed)",
       x = "Samples",
       y = "log2(raw counts)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#As this boxplot shows a need for normalization, the DESeq2 Pipeline is continued
counts_int <- round(txi$counts)
dds <- DESeqDataSetFromMatrix(countData = counts_int,
                              colData = coldata_mod,
                              design = ~ group)
keep <- apply(counts(dds), 1, function(x) all(x >= 10))
dds_filtered <- dds[keep, ]

#Running DESeq2
dds <- DESeq(dds_filtered)

normalized_counts <- counts(dds,normalized=TRUE)
norm_counts_df <- as.data.frame(normalized_counts)

#Adding Gene Identifiers
norm_counts_df$gene <- rownames(normalized_counts)

#Reshape the data from wide to long format
norm_counts_long <- melt(norm_counts_df, id.vars = "gene", 
                         variable.name = "sample", 
                         value.name = "counts")

#Create the boxplot with a log2 transformation
ggplot(norm_counts_long, aes(x = sample, y = log2(counts))) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Boxplot of Normalized Counts (log2-transformed)",
       x = "Samples",
       y = "log2(Normalized counts)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#The data is now normalized

#PCA plot without the faulty replicates
rld <- rlog(dds, blind = FALSE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
variance_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

#Dataframe containging two principal components and metadata
pca_df <- data.frame(PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2],
                     sample = colnames(rld),
                     group = rld$group,       # assumes 'group' exists in colData(rld)
                     replicate = rld$replicate)  # likewise for 'replicate'

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = replicate), size = 4, box.padding = 0.35, point.padding = 0.5) +
  xlab(paste0("PC1: ", round(variance_explained[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(variance_explained[2], 1), "% variance")) +
  ggtitle("PCA of RNA-seq Samples") +
  theme_bw()


#20H LPS DATA
#VOLCANO PLOTS (20H LPS)----

#Proteomics data already contain Log2 Fold Change data, this needs to be calculated for RNA-seq
#Log2 Fold Change for RNA-seq Data
res_1 <- results(dds, contrast = c("group", "Caco2-4ME-short", "Caco2-ctrl"))
res_2 <- results(dds, contrast = c("group", "Caco2-LPS", "Caco2-ctrl"))
res_3 <- results(dds, contrast = c("group", "Caco2-4ME+LPS", "Caco2-ctrl"))
res_4 <- results(dds, contrast = c("group", "Caco2-4ME-long", "Caco2-ctrl"))
res_5 <- results(dds, contrast = c("group", "Caco2-4ME+LPS+4ME", "Caco2-ctrl"))

#Using res_2 for 20H LPS Stimulation
resOrdered <- res_2[order(res_2$padj), ]
resOrdered <- resOrdered[!is.na(resOrdered$padj), ]
res_df <- as.data.frame(resOrdered)
res_df$gene <- rownames(res_df)

# Volcano plot for RNA_seq
p1 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = threshold), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "firebrick")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = subset(res_df, threshold == "Significant" & (abs(log2FoldChange) > 1.5 | -log10(padj) > 20)),
    aes(label = gene),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 50
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot: 20H LPS RNA-seq",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right"
  )

#And now for proteomics
#Columns of interest
fc_col <- "log2_FC_Ctrl_vs_20HLPS"
pval_col <- "ANOVA.p.value"
gene_col <- "PG.Genes"

#Thresholds 
log2fc_cutoff <- 1
pval_cutoff <- 0.05

#Flag significance
#Everything in the Proteomics is significant in at least one condition
filtered_df$Significance <- ifelse(
  filtered_df[[pval_col]] < pval_cutoff & abs(filtered_df[[fc_col]]) >= log2fc_cutoff,
  "Extremely Significant",
  ifelse(
    filtered_df[[pval_col]] < pval_cutoff,
    "Significant",
    "Not Significant"
  )
)

#Volcano plot for Proteomics
p2 <- ggplot(filtered_df, aes_string(x = fc_col, y = "log10pval", color = "Significance")) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(
    data = filtered_df[filtered_df$Significance == "Extremely Significant", ],
    aes_string(label = gene_col),
    size = 3.5,
    fontface = "bold",
    color = "steelblue",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Extremely Significant" = "steelblue",
      "Significant" = "grey50",
      "Not Significant" = "grey80"
    )
  ) +   # <-- closed parentheses here!
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = paste("Volcano Plot:", "20H LPS Proteomics"),
    x = "log2 Fold Change",
    y = "-log10(ANOVA p-value)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

#Visualizing side by side
library(patchwork)
p1+p2

#BARPLOT (20H LPS)----

library(dplyr)
library(reshape2)

#Top genes in RNA-seq
top_rna <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  arrange(desc(log2FoldChange)) %>%
  select(gene, log2FoldChange) %>%
  rename(RNA_seq = log2FoldChange)

#Top genes in Proteomics
top_prot <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_20HLPS > 1) %>%
  arrange(desc(log2_FC_Ctrl_vs_20HLPS)) %>%
  select(`PG.Genes`, log2_FC_Ctrl_vs_20HLPS) %>%
  rename(
    gene = `PG.Genes`, 
    Proteomics = log2_FC_Ctrl_vs_20HLPS
  )

#Making a merged data frame
merged_df <- full_join(top_rna, top_prot, by = "gene")

#Only keep the genes with a significant p-value < 0.05 and Log2(FC)=>1
merged_df <- merged_df %>% 
  filter(!is.na(RNA_seq) & !is.na(Proteomics))
merged_melted <- melt(merged_df, id.vars = "gene", variable.name = "Dataset", value.name = "log2FC")
merged_melted$Dataset <- recode(merged_melted$Dataset,
                                log2FC_RNA = "RNA_seq",
                                log2FC_Proteomics = "Proteomics")

#Barplot for proteomics and RNA_seq
ggplot(merged_melted, aes(x = reorder(gene, log2FC), y = log2FC, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("RNA_seq" = "firebrick", "Proteomics" = "steelblue")) +
  labs(title = "Most Upregulated Genes in Both Datasets, 20H LPS",
       x = "Gene",
       y = "Log2 Fold Change",
       fill = "Dataset") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#CORRELATION BETWEEN PROTEOMICS AND RNA-SEQ (20H LPS)----
#Finding correlation value
cor_value <- cor(merged_df_filtered$RNA_seq, merged_df_filtered$Proteomics)

#Scatter plot with correlation label
p3 <- ggplot(merged_df_filtered, aes(x = RNA_seq, y = Proteomics, label = gene)) +
  geom_point(color = "darkblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = paste0("RNA vs Proteomics (log2FC > 1) - 20h LPS\nPearson r = ", round(cor_value, 2)),
       x = "RNA log2 Fold Change",
       y = "Proteomics log2 Fold Change") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))


p4 <- ggplot(merged_df_filtered, aes(x = Proteomics, y = RNA_seq, label = gene)) +
  geom_point(color = "darkblue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = paste0("Proteomics vs RNA (log2FC > 1) - 20h LPS\nPearson r = ", round(cor_value, 2)),
       x = "Proteomics log2 Fold Change",
       y = "RNA log2 Fold Change") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

library(patchwork)
p3 + p4

#GO ANALYSIS (20H LPS)----

library(clusterProfiler)
library(org.Hs.eg.db)

#GO ANALYSIS FOR RNA-seq
#Filtering upregulated genes according to threshold
rna_up_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(gene) %>%
  unique()

#Converting SYMBOL to the corresponding ENTREZ ID
rna_gene_df <- bitr(rna_up_genes, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

#GO enrichment for Biological Process (BP)
rna_go <- enrichGO(gene         = rna_gene_df$ENTREZID,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

#Plot
p5 <- barplot(rna_go, showCategory = 10, title = "Top GO Terms - RNA-seq (BP)")

#GO ANALYSIS FOR PROTEOMICS 
#Filtering genes according to threshold
prot_up_genes <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_20HLPS > 1) %>%
  pull(`PG.Genes`) %>%
  unique()

#Converting SYMBOL to the corresponding ENTREZ ID
prot_gene_df <- bitr(prot_up_genes, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

#GO enrichment (Biological Process)
prot_go <- enrichGO(gene         = prot_gene_df$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

#Plot
p6 <- barplot(prot_go, showCategory = 10, title = "Top GO Terms - Proteomics (BP)")

#Combined 
p5+p6

library(dplyr)

#GO PLOTS FOR SPECIFIC TERMS RELATED TO INFLAMMATION
#Conversion to a data frame
rna_go_df <- as.data.frame(rna_go)
prot_go_df <- as.data.frame(prot_go)

#Using inflammation keywords and intresting pathways
inflammation_keywords <- c("inflamm", "cytokine", "interleukin", "immune", "TNF", "NF-kB", "NF-κB", "MAPK", "ERK", "JNK", "p38")
pattern <- paste(inflammation_keywords, collapse = "|")

#Combination into a regex pattern
pattern <- paste(inflammation_keywords, collapse = "|")

#Subsetting GO terms
rna_inflam <- rna_go_df %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))

prot_inflam <- prot_go_df %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))

#Naming datasets 
rna_inflam$Dataset <- "RNA-seq"
prot_inflam$Dataset <- "Proteomics"

#Conbining the datasets
combined_go <- rbind(rna_inflam, prot_inflam)

#Combined plot
ggplot(combined_go, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related GO Terms (RNA-seq vs Proteomics)",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("RNA-seq" = "firebrick", "Proteomics" = "steelblue")) +
  theme_minimal(base_size = 14)

#KEGG ANALYSIS (20H LPS)----

#RNA-seq
rna_up <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(gene)

#Proteomics
prot_up <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_20HLPS > 1) %>%
  pull(`PG.Genes`)

#Convert SYMBOL to ENTREZ
rna_entrez <- bitr(rna_up, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

prot_entrez <- bitr(prot_up, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

#RNA-seq KEGG
rna_kegg <- enrichKEGG(gene = rna_entrez$ENTREZID,
                       organism = "hsa",  # human
                       pvalueCutoff = 0.05)

#Proteomics KEGG
prot_kegg <- enrichKEGG(gene = prot_entrez$ENTREZID,
                        organism = "hsa",
                        pvalueCutoff = 0.05)

#Barplots
p7 <- barplot(rna_kegg, showCategory = 10, title = "Top KEGG Pathways - RNA-seq")
p8 <- barplot(prot_kegg, showCategory = 10, title = "Top KEGG Pathways - Proteomics")

#Combined
p7+p8

#Inflammation Related KEGG
inflammation_kegg_terms <- c("TNF", "cytokine", "inflamm", "MAPK", "NF-kB", "IL-", "JAK", "TLR", "chemokine", "immune")
pattern <- paste(inflammation_kegg_terms, collapse = "|")

#Convert to data frames
rna_kegg_df <- as.data.frame(rna_kegg)
prot_kegg_df <- as.data.frame(prot_kegg)

#Filter for inflammation-related pathways
rna_kegg_inflam <- rna_kegg_df %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))

prot_kegg_inflam <- prot_kegg_df %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))

#Name
rna_kegg_inflam$Dataset <- "RNA-seq"
prot_kegg_inflam$Dataset <- "Proteomics"

#Combine
combined_kegg <- rbind(rna_kegg_inflam, prot_kegg_inflam)

#Combined KEGG
ggplot(combined_kegg, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related KEGG Pathways (RNA-seq vs Proteomics)",
       x = "Pathway",
       y = "Gene Count") +
  scale_fill_manual(values = c("RNA-seq" = "firebrick", "Proteomics" = "steelblue")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))




#PATHVIEW (20H LPS)----

#RNA-seq
rna_map <- bitr(res_df$gene, fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#Expression data
rna_fc <- merge(rna_map, res_df, by.x = "SYMBOL", by.y = "gene") %>%
  dplyr::select(ENTREZID, log2FoldChange)

#Create named vector
rna_vec <- setNames(rna_fc$log2FoldChange, rna_fc$ENTREZID)

#Pathview
pathview(
  gene.data  = rna_vec,
  pathway.id = "hsa04010",
  species    = "hsa",
  out.suffix = "RNAseq_MAPK",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)

#Proteomics
prot_map <- bitr(filtered_df$`PG.Genes`, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#Fold changes
prot_fc <- merge(prot_map, filtered_df, by.x = "SYMBOL", by.y = "PG.Genes") %>%
  dplyr::select(ENTREZID, log2_FC_Ctrl_vs_20HLPS)

#Create named vector
prot_vec <- setNames(prot_fc$log2_FC_Ctrl_vs_20HLPS, prot_fc$ENTREZID)

#Pathwiev
pathview(
  gene.data  = prot_vec,
  pathway.id = "hsa04010",
  species    = "hsa",
  out.suffix = "Proteomics_MAPK",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)

#NF-kB Pathway, RNA-seq
pathview(
  gene.data  = rna_vec,
  pathway.id = "hsa04064",  # NF-κB pathway
  species    = "hsa",
  out.suffix = "RNAseq_NFkB",
  kegg.native = TRUE,
  low  = list(gene = "blue"),
  mid  = list(gene = "white"),
  high = list(gene = "red")
)

#NF-kB Pathway, Proteomics
pathview(
  gene.data  = prot_vec,
  pathway.id = "hsa04064",  # NF-κB pathway
  species    = "hsa",
  out.suffix = "Proteomics_NFkB",
  kegg.native = TRUE,
  low  = list(gene = "blue"),
  mid  = list(gene = "white"),
  high = list(gene = "red")
)




#20H LPS + TREATMENT DATA----

#BARPLOT----

#Proteomics
library(ggplot2)
library(reshape2)

#Using the same genes of interest as before
genes_of_interest <- c("CYP1A1", "FOS", "ID1", "JUN", "MYC", "JUNB", "CDKN1A")

#Subsetting data
selected_df <- filtered_df[filtered_df$`PG.Genes` %in% genes_of_interest,
                           c("PG.Genes",
                             "log2_FC_Ctrl_vs_20HLPS",
                             "log2_FC_Ctrl_vs_2H4ME20HLPS",
                             "log2_FC_Ctrl_vs_22H4ME20HLPS")]

#Melt into long format for plotting
selected_long <- melt(selected_df,
                      id.vars = "PG.Genes",
                      variable.name = "Condition",
                      value.name = "log2FC")

#Condition labels for clearer visualization
selected_long$Condition <- factor(selected_long$Condition,
                                  levels = c("log2_FC_Ctrl_vs_20HLPS",
                                             "log2_FC_Ctrl_vs_2H4ME20HLPS",
                                             "log2_FC_Ctrl_vs_22H4ME20HLPS"),
                                  labels = c("20H LPS",
                                             "2H 4-ME + 20H LPS",
                                             "22H 4-ME + 20H LPS"))

#Barplot Proteomics
ggplot(selected_long, aes(x = PG.Genes, y = log2FC, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  labs(title = "Upregulated Genes Across Conditions: Proteomics",
       x = "Gene",
       y = "log2 Fold Change",
       fill = "Condition") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

#RNA-seq
#Genes of interest
genes_of_interest <- c("CYP1A1", "FOS", "ID1", "JUN", "MYC", "JUNB", "CDKN1A")

#Preparation RNA datasets
#Data frame
res2_df <- as.data.frame(res_2)
res2_df$gene <- rownames(res2_df)

res3_df <- as.data.frame(res_3)
res3_df$gene <- rownames(res3_df)

res5_df <- as.data.frame(res_5)
res5_df$gene <- rownames(res5_df)

#Subset only genes of interest
res2_sel <- res2_df %>% dplyr::filter(gene %in% genes_of_interest) %>% dplyr::select(gene, log2FoldChange)
res3_sel <- res3_df %>% dplyr::filter(gene %in% genes_of_interest) %>% dplyr::select(gene, log2FoldChange)
res5_sel <- res5_df %>% dplyr::filter(gene %in% genes_of_interest) %>% dplyr::select(gene, log2FoldChange)

#Rename columns to match
colnames(res2_sel)[2] <- "20H LPS"
colnames(res3_sel)[2] <- "2H 4-ME + 20H LPS"
colnames(res5_sel)[2] <- "22H 4-ME + 20H LPS"

#Merge
merged_rna <- dplyr::full_join(res2_sel, res3_sel, by = "gene") %>%
  dplyr::full_join(res5_sel, by = "gene")

#Melt to long format
rna_long <- melt(merged_rna, id.vars = "gene", variable.name = "Condition", value.name = "log2FC")

#Plot
ggplot(rna_long, aes(x = gene, y = log2FC, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  labs(title = "Upregulated Genes Across Conditions: RNA-seq",
       x = "Gene",
       y = "log2 Fold Change",
       fill = "Condition") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))




#GO ANALYSIS (20H LPS, 2H 4-ME)----

#RNA-seq (20H LPS + 22H 4-ME)
# Convert to data frame
res5_df <- as.data.frame(res_5)
res5_df$gene <- rownames(res5_df)

# Filter significant upregulated or all DEGs
sig_genes_5 <- res5_df %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

entrez_ids_5 <- bitr(sig_genes_5$gene, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_res_5 <- enrichGO(gene          = entrez_ids_5$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",         # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)

#RNA-seq (20H LPS + 2H 4-ME)
res3_df <- as.data.frame(res_3)
res3_df$gene <- rownames(res3_df)

# Filter significant upregulated or all DEGs
sig_genes_3 <- res3_df %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

entrez_ids_3 <- bitr(sig_genes_3$gene, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_res_3 <- enrichGO(gene          = entrez_ids_3$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",         # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)

p9 <- barplot(go_res_5, showCategory = 10, title = "20H LPS + 22H 4-ME: RNA-seq")
p10 <- barplot(go_res_3, showCategory = 10, title = "20H LPS + 2H 4-ME: RNA-seq")

#Combined plot
p9 + p10

#PROTEOMICS (20H LPS + 22H 4-ME)
#Filtering genes according to threshold
prot_up_genes1 <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_22H4ME20HLPS > 1) %>%
  pull(`PG.Genes`) %>%
  unique()

#Converting SYMBOL to the corresponding ENTREZ ID
prot_gene_df1 <- bitr(prot_up_genes1, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

#GO enrichment (Biological Process)
prot_go1 <- enrichGO(gene         = prot_gene_df1$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

#Plot
p11 <- barplot(prot_go1, showCategory = 10, title = "20H LPS + 22H 4-ME: Proteomics")

#PROTEOMICS (20H LPS + 2H 4-ME)
#Filtering genes according to threshold
prot_up_genes2 <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_2H4ME20HLPS > 1) %>%
  pull(`PG.Genes`) %>%
  unique()

#Converting SYMBOL to the corresponding ENTREZ ID
prot_gene_df2 <- bitr(prot_up_genes2, fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

#GO enrichment (Biological Process)
prot_go2 <- enrichGO(gene         = prot_gene_df2$ENTREZID,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

#Plot
p12 <- barplot(prot_go2, showCategory = 10, title = "20H LPS + 2H 4-ME: Proteomics")

#COMBINED
p11 + p12

#GO ANALYSIS (INFLAMMATION, LPS+TREATMENT)----

#GO PLOTS FOR SPECIFIC TERMS RELATED TO INFLAMMATION
#Conversion to a data frame
rna_go_df_5 <- as.data.frame(go_res_5)
rna_go_df_3 <- as.data.frame(go_res_3)

prot_go_df1 <- as.data.frame(prot_go1)
prot_go_df2 <- as.data.frame(prot_go2)

#Using inflammation keywords and intresting pathways
inflammation_keywords <- c("inflamm", "cytokine", "interleukin", "immune", "TNF", "NF-kB", "NF-κB", "MAPK", "ERK", "JNK", "p38")
pattern <- paste(inflammation_keywords, collapse = "|")

#Combination into a regex pattern
pattern <- paste(inflammation_keywords, collapse = "|")

#Subsetting GO terms
rna_inflam_5 <- rna_go_df_5 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))
rna_inflam_3 <- rna_go_df_3 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))


prot_inflam1 <- prot_go_df1 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))
prot_inflam2 <- prot_go_df2 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))

#Naming datasets 
rna_inflam_5$Dataset <- "20H LPS + 22H 4-ME"
rna_inflam_3$Dataset <- "20H LPS + 2H 4-ME"
prot_inflam1$Dataset <- "20H LPS + 22H 4-ME"
#ERROR CODE IN PROTEOMICS, SINCE WE HAVE NO INFLAMMATION GENES IN 20H LPS + 22H 4-ME
#Have to add this 
if (nrow(prot_inflam1) > 0) {
  prot_inflam1$Dataset <- "20H LPS + 22H 4-ME"
} else {
  message("⚠️ No inflammation-related GO terms found in prot_inflam1.")
}
prot_inflam2$Dataset <- "20H LPS + 2H 4-ME"

#Conbining the datasets
combined_go_RNA.treat <- rbind(rna_inflam_5, rna_inflam_3)
combined_go_Pro.treat <- rbind(prot_inflam1, prot_inflam2)

#Combined plot
ggplot(combined_go_RNA.treat, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related GO Terms: RNA-seq",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF")) +
  theme_minimal(base_size = 14)

ggplot(combined_go_Pro.treat, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related GO Terms: Proteomics",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF")) +
  theme_minimal(base_size = 14)

#One that also have 20H LPS
rna_inflam$Dataset <- "20H LPS"
prot_inflam$Dataset <- "20H LPS"

#Conbining the datasets
combined_go_RNA.treat1 <- rbind(rna_inflam_5, rna_inflam_3, rna_inflam)
combined_go_Pro.treat1 <- rbind(prot_inflam1, prot_inflam2, prot_inflam)

ggplot(combined_go_RNA.treat1, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related GO Terms: RNA-seq",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS" = "#F8766D", "20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF")) +
  theme_minimal(base_size = 14)

ggplot(combined_go_Pro.treat1, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related GO Terms: Proteomics",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS" = "#F8766D", "20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF")) +
  theme_minimal(base_size = 14)

#KEGG (LPS+TREATMENT)----

#KEGG for 20H LPS, 22H 4-ME (RNA-seq)
library(clusterProfiler)
library(org.Hs.eg.db)

# Prepare DESeq2 results
res5_df <- as.data.frame(res_5)
res5_df$gene <- rownames(res5_df)

# Filter significant upregulated genes
sig_res5 <- res5_df %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 1)

# Map to Entrez IDs
res5_entrez <- bitr(sig_res5$gene, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# Run KEGG enrichment
kegg_res5 <- enrichKEGG(gene         = res5_entrez$ENTREZID,
                        organism     = "hsa",
                        pvalueCutoff = 0.05)

# Make readable (convert Entrez back to gene symbols)
kegg_res5 <- setReadable(kegg_res5, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Plot
p13 <- barplot(kegg_res5, showCategory = 15, title = "KEGG - 20H LPS + 22H 4-ME: RNA-seq")

res3_df <- as.data.frame(res_3)
res3_df$gene <- rownames(res3_df)

sig_res3 <- res3_df %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 1)

res3_entrez <- bitr(sig_res3$gene, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

kegg_res3 <- enrichKEGG(gene         = res3_entrez$ENTREZID,
                        organism     = "hsa",
                        pvalueCutoff = 0.05)

kegg_res3 <- setReadable(kegg_res3, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

p14 <- barplot(kegg_res3, showCategory = 15, title = "KEGG - 20H LPS + 2H 4-ME : RNA-seq")

p13 + p14

#KEGG FOR PROTEOMICS

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

#2H4ME + 20H LPS
up_prot_2H <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_2H4ME20HLPS > 1) %>%
  pull(`PG.Genes`) %>%
  unique()

rna_entrez_2H <- bitr(up_prot_2H, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

rna_kegg_2H <- enrichKEGG(gene = rna_entrez_2H$ENTREZID,
                           organism = "hsa",  # human
                           pvalueCutoff = 0.05)

p15 <- barplot(rna_kegg_2H, showCategory = 10, title = "KEGG - 20H LPS + 2H 4-ME: Proteomics")

#22H4ME + 20H LPS
up_prot_22H <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_22H4ME20HLPS > 1) %>%
  pull(`PG.Genes`) %>%
  unique()

rna_entrez_22H <- bitr(up_prot_22H, fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

rna_kegg_22H <- enrichKEGG(gene = rna_entrez_22H$ENTREZID,
                       organism = "hsa",  # human
                       pvalueCutoff = 0.05)

p16 <- barplot(rna_kegg_22H, showCategory = 10, title = "KEGG - 20H LPS + 22H 4-ME: Proteomics")

p15 + p16

#INFLAMMATION KEGG (LPS+TREATMENT)----

#Inflammation Related KEGG
inflammation_kegg_terms <- c("TNF", "cytokine", "inflamm", "MAPK", "NF-kB", "IL-", "JAK", "TLR", "chemokine", "immune")
pattern <- paste(inflammation_kegg_terms, collapse = "|")

#Convert to data frames
rna_kegg_df_3 <- as.data.frame(kegg_res3)
rna_kegg_df_5 <- as.data.frame(kegg_res5)
rna_kegg_22H_df <- as.data.frame(rna_kegg_22H)
rna_kegg_2H_df <- as.data.frame(rna_kegg_2H)

#Filter for inflammation-related pathways
rna_kegg_inflam_3 <- rna_kegg_df_3 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))
rna_kegg_inflam_5 <- rna_kegg_df_5 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))

prot_kegg_inflam_22H <- rna_kegg_22H_df %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))
prot_kegg_inflam_2H <- rna_kegg_2H_df %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))


#Name
rna_kegg_inflam_3$Dataset <- "20H LPS + 2H 4-ME"
rna_kegg_inflam_5$Dataset <- "20H LPS + 22H 4-ME"
prot_kegg_inflam_2H$Dataset <- "20H LPS + 2H 4-ME"
prot_kegg_inflam_22H$Dataset <- "20H LPS + 22H 4-ME"

#Combine
combined_kegg_rna <- rbind(rna_kegg_inflam_3, rna_kegg_inflam_5)
combined_kegg_pro <- rbind(prot_kegg_inflam_2H, prot_kegg_inflam_22H)

#Combined KEGG
ggplot(combined_kegg_rna, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related KEGG Pathways: RNA-seq",
       x = "Pathway",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggplot(combined_kegg_pro, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related KEGG Pathways: Proteomics",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF")) +
  theme_minimal(base_size = 14)

#Naming 20H LPS DATA
rna_kegg_inflam$Dataset <- "20H LPS"
prot_kegg_inflam$Dataset <- "20H LPS"

#Conbining the datasets
combined_go_RNA.treat_kegg <- rbind(rna_kegg_inflam, rna_kegg_inflam_3, rna_kegg_inflam_5)
combined_go_Pro.treat_kegg <- rbind(prot_kegg_inflam, prot_kegg_inflam_2H, prot_kegg_inflam_22H)

ggplot(combined_go_RNA.treat_kegg, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related KEGG Terms: RNA-seq",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS" = "#F8766D", "20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF")) +
  theme_minimal(base_size = 14)

ggplot(combined_go_Pro.treat_kegg, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related KEGG Terms: Proteomics",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS" = "#F8766D", "20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF")) +
  theme_minimal(base_size = 14)

#PATHVIEW (LPS+TREATMENT)----
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(dplyr)

#RNA-seq - res_3 and res_5

# Prepare res_3 (2H 4-ME + 20H LPS)
res3_df <- as.data.frame(res_3)
res3_df$gene <- rownames(res3_df)

rna3_map <- bitr(res3_df$gene, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)

rna3_fc <- merge(rna3_map, res3_df, by.x = "SYMBOL", by.y = "gene") %>%
  dplyr::select(ENTREZID, log2FoldChange)

rna3_vec <- setNames(rna3_fc$log2FoldChange, rna3_fc$ENTREZID)

# Pathview: MAPK for res_3
pathview(
  gene.data  = rna3_vec,
  pathway.id = "hsa04010",  # MAPK signaling
  species    = "hsa",
  out.suffix = "RNAseq_res3_MAPK",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)

# Pathview: NF-kB for res_3
pathview(
  gene.data  = rna3_vec,
  pathway.id = "hsa04064",  # NF-kB signaling
  species    = "hsa",
  out.suffix = "RNAseq_res3_NFkB",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)


# Prepare res_5 (22H 4-ME + 20H LPS)
res5_df <- as.data.frame(res_5)
res5_df$gene <- rownames(res5_df)

rna5_map <- bitr(res5_df$gene, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)

rna5_fc <- merge(rna5_map, res5_df, by.x = "SYMBOL", by.y = "gene") %>%
  dplyr::select(ENTREZID, log2FoldChange)

rna5_vec <- setNames(rna5_fc$log2FoldChange, rna5_fc$ENTREZID)

# Pathview: MAPK for res_5
pathview(
  gene.data  = rna5_vec,
  pathway.id = "hsa04010",
  species    = "hsa",
  out.suffix = "RNAseq_res5_MAPK",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)

# Pathview: NF-kB for res_5
pathview(
  gene.data  = rna5_vec,
  pathway.id = "hsa04064",
  species    = "hsa",
  out.suffix = "RNAseq_res5_NFkB",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)


#Proteomics

# 2H 4-ME + 20H LPS
prot_map_2H <- bitr(filtered_df$`PG.Genes`, fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)

prot_fc_2H <- merge(prot_map_2H, filtered_df, by.x = "SYMBOL", by.y = "PG.Genes") %>%
  dplyr::select(ENTREZID, log2_FC_Ctrl_vs_2H4ME20HLPS)

prot_vec_2H <- setNames(prot_fc_2H$log2_FC_Ctrl_vs_2H4ME20HLPS, prot_fc_2H$ENTREZID)

pathview(
  gene.data  = prot_vec_2H,
  pathway.id = "hsa04010",
  species    = "hsa",
  out.suffix = "Proteomics_2H_MAPK",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)

pathview(
  gene.data  = prot_vec_2H,
  pathway.id = "hsa04064",
  species    = "hsa",
  out.suffix = "Proteomics_2H_NFkB",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)


# 22H 4-ME + 20H LPS
prot_fc_22H <- merge(prot_map_2H, filtered_df, by.x = "SYMBOL", by.y = "PG.Genes") %>%
  dplyr::select(ENTREZID, log2_FC_Ctrl_vs_22H4ME20HLPS)

prot_vec_22H <- setNames(prot_fc_22H$log2_FC_Ctrl_vs_22H4ME20HLPS, prot_fc_22H$ENTREZID)

pathview(
  gene.data  = prot_vec_22H,
  pathway.id = "hsa04010",
  species    = "hsa",
  out.suffix = "Proteomics_22H_MAPK",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)

pathview(
  gene.data  = prot_vec_22H,
  pathway.id = "hsa04064",
  species    = "hsa",
  out.suffix = "Proteomics_22H_NFkB",
  kegg.native = TRUE,
  low = list(gene = "blue"),
  mid = list(gene = "white"),
  high = list(gene = "red")
)


#4-ME TREATMENT DATA----
#BARPLOT (ALL)----

#Proteomics
#Using the same genes of interest as before
genes_of_interest <- c("CYP1A1", "FOS", "ID1", "JUN", "MYC", "JUNB", "CDKN1A")

#Subsetting data
selected_df <- filtered_df[filtered_df$`PG.Genes` %in% genes_of_interest,
                           c("PG.Genes",
                             "log2_FC_Ctrl_vs_20HLPS",
                             "log2_FC_Ctrl_vs_2H4ME20HLPS",
                             "log2_FC_Ctrl_vs_22H4ME20HLPS",
                             "log2_FC_Ctrl_vs_2h4ME",
                             "log2_FC_Ctrl_vs_22H4ME")]

#Melt into long format for plotting
selected_long <- melt(selected_df,
                      id.vars = "PG.Genes",
                      variable.name = "Condition",
                      value.name = "log2FC")

#Condition labels for clearer visualization
selected_long$Condition <- factor(selected_long$Condition,
                                  levels = c("log2_FC_Ctrl_vs_20HLPS",
                                             "log2_FC_Ctrl_vs_2H4ME20HLPS",
                                             "log2_FC_Ctrl_vs_22H4ME20HLPS",
                                             "log2_FC_Ctrl_vs_2h4ME",
                                             "log2_FC_Ctrl_vs_22H4ME"),
                                  labels = c("20H LPS",
                                             "2H 4-ME + 20H LPS",
                                             "22H 4-ME + 20H LPS",
                                             "2H 4-ME",
                                             "22H 4-ME"))

#Barplot Proteomics
ggplot(selected_long, aes(x = PG.Genes, y = log2FC, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c(
    "20H LPS" = "#F8766D",
    "2H 4-ME + 20H LPS" = "#00BA38",
    "22H 4-ME + 20H LPS" = "#619CFF",
    "2H 4-ME" = "yellow2",
    "22H 4-ME" = "orange"
  )) +
  labs(title = "Upregulated Genes Across Conditions: Proteomics",
       x = "Gene",
       y = "log2 Fold Change",
       fill = "Condition") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#RNA-seq
#Genes of interest
genes_of_interest <- c("CYP1A1", "FOS", "ID1", "JUN", "MYC", "JUNB", "CDKN1A")

#Preparation RNA datasets
#Data frame
res2_df <- as.data.frame(res_2)
res2_df$gene <- rownames(res2_df)

res3_df <- as.data.frame(res_3)
res3_df$gene <- rownames(res3_df)

res5_df <- as.data.frame(res_5)
res5_df$gene <- rownames(res5_df)

res1_df <- as.data.frame(res_1)
res1_df$gene <- rownames(res1_df)

res4_df <- as.data.frame(res_4)
res4_df$gene <- rownames(res4_df)

#Subset only genes of interest
res2_sel <- res2_df %>% dplyr::filter(gene %in% genes_of_interest) %>% dplyr::select(gene, log2FoldChange)
res3_sel <- res3_df %>% dplyr::filter(gene %in% genes_of_interest) %>% dplyr::select(gene, log2FoldChange)
res5_sel <- res5_df %>% dplyr::filter(gene %in% genes_of_interest) %>% dplyr::select(gene, log2FoldChange)
res1_sel <- res1_df %>% dplyr::filter(gene %in% genes_of_interest) %>% dplyr::select(gene, log2FoldChange)
res4_sel <- res4_df %>% dplyr::filter(gene %in% genes_of_interest) %>% dplyr::select(gene, log2FoldChange)

#Rename columns to match
colnames(res2_sel)[2] <- "20H LPS"
colnames(res3_sel)[2] <- "2H 4-ME + 20H LPS"
colnames(res5_sel)[2] <- "22H 4-ME + 20H LPS"
colnames(res1_sel)[2] <- "2H 4-ME"
colnames(res4_sel)[2] <- "22H 4-ME"

#Merge
merged_rna <- res2_sel %>%
  full_join(res3_sel, by = "gene") %>%
  full_join(res1_sel, by = "gene") %>%
  full_join(res4_sel, by = "gene") %>%
  full_join(res5_sel, by = "gene")

#Melt to long format
rna_long <- melt(merged_rna, id.vars = "gene", variable.name = "Condition", value.name = "log2FC")

#Plot
ggplot(rna_long, aes(x = gene, y = log2FC, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c(
    "20H LPS" = "#F8766D",
    "2H 4-ME + 20H LPS" = "#00BA38",
    "22H 4-ME + 20H LPS" = "#619CFF",
    "2H 4-ME" = "yellow2",
    "22H 4-ME" = "orange"
  )) +
  labs(title = "Upregulated Genes Across Conditions: RNA-seq",
       x = "Gene",
       y = "log2 Fold Change",
       fill = "Condition") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#GO ANALYSIS (TREATMENT)----

#RNA-seq (2H 4-ME)
# Convert to data frame
res1_df <- as.data.frame(res_1)
res1_df$gene <- rownames(res1_df)

# Filter significant upregulated or all DEGs
sig_genes_1 <- res1_df %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

entrez_ids_1 <- bitr(sig_genes_1$gene, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_res_1 <- enrichGO(gene          = entrez_ids_1$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",         # Biological Process
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)

#RNA-seq (22H 4-ME)
res4_df <- as.data.frame(res_4)
res4_df$gene <- rownames(res4_df)

# Filter significant upregulated or all DEGs
sig_genes_4 <- res4_df %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

entrez_ids_4 <- bitr(sig_genes_4$gene, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_res_4 <- enrichGO(gene          = entrez_ids_4$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",         # Biological Process
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)

p17 <- barplot(go_res_1, showCategory = 10, title = "2H 4-ME: RNA-seq")
p18 <- barplot(go_res_4, showCategory = 10, title = "22H 4-ME: RNA-seq")

#Combined plot
p17 + p18

#PROTEOMICS (2H 4-ME)
#Filtering genes according to threshold
prot_up_genes3 <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_2h4ME > 1) %>%
  pull(`PG.Genes`) %>%
  unique()

#Converting SYMBOL to the corresponding ENTREZ ID
prot_gene_df3 <- bitr(prot_up_genes3, fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

#GO enrichment (Biological Process)
prot_go3 <- enrichGO(gene         = prot_gene_df3$ENTREZID,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

#Plot
p19 <- barplot(prot_go3, showCategory = 10, title = "2H 4-ME: Proteomics")

#PROTEOMICS (22H 4-ME)
#Filtering genes according to threshold
prot_up_genes4 <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_22H4ME > 1) %>%
  pull(`PG.Genes`) %>%
  unique()

#Converting SYMBOL to the corresponding ENTREZ ID
prot_gene_df4 <- bitr(prot_up_genes4, fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

#GO enrichment (Biological Process)
prot_go4 <- enrichGO(gene         = prot_gene_df4$ENTREZID,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

p20 <- barplot(prot_go4, showCategory = 10, title = "22H 4-ME: Proteomics")

p19 + p20

#GO ANALYSIS (INFLAMMATION 4-ME)----

#GO PLOTS FOR SPECIFIC TERMS RELATED TO INFLAMMATION
#Conversion to a data frame
rna_go_df_1 <- as.data.frame(go_res_1)
rna_go_df_4 <- as.data.frame(go_res_4)

prot_go_df3 <- as.data.frame(prot_go3)
prot_go_df4 <- as.data.frame(prot_go4)

#Using inflammation keywords and intresting pathways
inflammation_keywords <- c("inflamm", "cytokine", "interleukin", "immune", "TNF", "NF-kB", "NF-κB", "MAPK", "ERK", "JNK", "p38")
pattern <- paste(inflammation_keywords, collapse = "|")

#Combination into a regex pattern
pattern <- paste(inflammation_keywords, collapse = "|")

#Subsetting GO terms
rna_inflam_1 <- rna_go_df_1 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))
rna_inflam_4 <- rna_go_df_4 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))


prot_inflam3 <- prot_go_df3 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))
prot_inflam4 <- prot_go_df4 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))

#Naming datasets 
rna_inflam_1$Dataset <- "2H 4-ME"
rna_inflam_4$Dataset <- "22H 4-ME"
prot_inflam3$Dataset <- "2H 4-ME"
prot_inflam4$Dataset <- "22H 4-ME"

#Conbining the datasets
combined_go_RNA.treat4ME <- rbind(rna_inflam_1, rna_inflam_4)
combined_go_Pro.treat4ME <- rbind(prot_inflam3, prot_inflam4)

#Combined plot
ggplot(combined_go_RNA.treat4ME, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related GO Terms: RNA-seq",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("2H 4-ME" = "yellow2", "22H 4-ME" = "orange")) +
  theme_minimal(base_size = 14)

ggplot(combined_go_Pro.treat4ME, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related GO Terms: Proteomics",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("2H 4-ME" = "yellow2", "22H 4-ME" = "orange")) +
  theme_minimal(base_size = 14)

#One that has all of them
rna_inflam$Dataset <- "20H LPS"
prot_inflam$Dataset <- "20H LPS"

#Conbining the datasets
combined_go_RNA.treat.all <- rbind(rna_inflam_5, rna_inflam_3, rna_inflam, rna_inflam_1, rna_inflam_4)
combined_go_Pro.treat.all <- rbind(prot_inflam1, prot_inflam2, prot_inflam, prot_inflam3, prot_inflam4)

ggplot(combined_go_RNA.treat.all, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related GO Terms: RNA-seq",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS" = "#F8766D", "20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF", "2H 4-ME" = "yellow2", "22H 4-ME" = "orange")) +
  theme_minimal(base_size = 14)

ggplot(combined_go_Pro.treat.all, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related GO Terms: Proteomics",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS" = "#F8766D", "20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF", "2H 4-ME" = "yellow2", "22H 4-ME" = "orange")) +
  theme_minimal(base_size = 14)

#KEGG ANLYSIS 4-ME----

# Prepare DESeq2 results
res1_df <- as.data.frame(res_1)
res1_df$gene <- rownames(res1_df)

# Filter significant upregulated genes
sig_res1 <- res1_df %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 1)

# Map to Entrez IDs
res1_entrez <- bitr(sig_res1$gene, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# Run KEGG enrichment
kegg_res1 <- enrichKEGG(gene         = res1_entrez$ENTREZID,
                        organism     = "hsa",
                        pvalueCutoff = 0.05)

# Make readable (convert Entrez back to gene symbols)
kegg_res1 <- setReadable(kegg_res1, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Plot
p21 <- barplot(kegg_res1, showCategory = 15, title = "KEGG - 2H 4-ME: RNA-seq")

res4_df <- as.data.frame(res_4)
res4_df$gene <- rownames(res4_df)

sig_res4 <- res4_df %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 1)

res4_entrez <- bitr(sig_res4$gene, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

kegg_res4 <- enrichKEGG(gene         = res4_entrez$ENTREZID,
                        organism     = "hsa",
                        pvalueCutoff = 0.05)

kegg_res4 <- setReadable(kegg_res4, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

p22 <- barplot(kegg_res4, showCategory = 15, title = "22H 4-ME : RNA-seq")

p22 + p21

library(ggplotify)
library(patchwork)

p21 <- as.ggplot(barplot(kegg_res1, showCategory = 15, title = "2H 4-ME: RNA-seq"))
p22 <- as.ggplot(barplot(kegg_res4, showCategory = 15, title = "22H 4-ME: RNA-seq"))

p21 + p22
as.data.frame(kegg_res4)
#NO KEGG enrichment in 22H 4-ME for RNA-seq

p21

#KEGG FOR PROTEOMICS

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

#2H4ME
up_prot_2H4ME <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_2h4ME > 1) %>%
  pull(`PG.Genes`) %>%
  unique()

rna_entrez_2H4ME <- bitr(up_prot_2H4ME, fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

rna_kegg_2H4ME <- enrichKEGG(gene = rna_entrez_2H4ME$ENTREZID,
                          organism = "hsa",  # human
                          pvalueCutoff = 0.05)

p23 <- barplot(rna_kegg_2H4ME, showCategory = 10, title = "KEGG - 2H 4-ME: Proteomics")

#22H4ME
up_prot_22H4ME <- filtered_df %>%
  filter(ANOVA.p.value < 0.05 & log2_FC_Ctrl_vs_22H4ME > 1) %>%
  pull(`PG.Genes`) %>%
  unique()

rna_entrez_22H4ME <- bitr(up_prot_22H4ME, fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)

rna_kegg_22H4ME <- enrichKEGG(gene = rna_entrez_22H4ME$ENTREZID,
                           organism = "hsa",  # human
                           pvalueCutoff = 0.05)

p24 <- barplot(rna_kegg_22H4ME, showCategory = 10, title = "KEGG - 22H 4-ME: Proteomics")

p23 + p24

#KEGG ANALYSIS (INFLAMMATION 4-ME)----

#Inflammation Related KEGG
inflammation_kegg_terms <- c("TNF", "cytokine", "inflamm", "MAPK", "NF-kB", "IL-", "JAK", "TLR", "chemokine", "immune")
pattern <- paste(inflammation_kegg_terms, collapse = "|")

#Convert to data frames
rna_kegg_df_1 <- as.data.frame(kegg_res1)
rna_kegg_df_4 <- as.data.frame(kegg_res4)
rna_kegg_22H4ME_df <- as.data.frame(rna_kegg_22H4ME)
rna_kegg_2H4ME_df <- as.data.frame(rna_kegg_2H4ME)

#Filter for inflammation-related pathways
rna_kegg_inflam_1 <- rna_kegg_df_1 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))
rna_kegg_inflam_4 <- rna_kegg_df_4 %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))

prot_kegg_inflam_22H4ME <- rna_kegg_22H4ME_df %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))
prot_kegg_inflam_2H4ME <- rna_kegg_2H4ME_df %>%
  filter(grepl(pattern, Description, ignore.case = TRUE))


#Name
rna_kegg_inflam_1$Dataset <- "2H 4-ME"
rna_kegg_inflam_4$Dataset <- "22H 4-ME"
prot_kegg_inflam_2H4ME$Dataset <- "2H 4-ME"
prot_kegg_inflam_22H4ME$Dataset <- "22H 4-ME"

#Combine
combined_kegg_rna_4ME <- rbind(rna_kegg_inflam_1, rna_kegg_inflam_4)
combined_kegg_pro_4ME <- rbind(prot_kegg_inflam_2H4ME, prot_kegg_inflam_22H4ME)

#Combined KEGG
ggplot(combined_kegg_rna_4ME, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related KEGG Pathways: RNA-seq",
       x = "Pathway",
       y = "Gene Count") +
  scale_fill_manual(values = c("2H 4-ME" = "yellow2", "22H 4-ME" = "orange")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggplot(combined_kegg_pro_4ME, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related KEGG Pathways: Proteomics",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("2H 4-ME" = "yellow2", "22H 4-ME" = "orange")) +
  theme_minimal(base_size = 14)

#Naming 20H LPS DATA
rna_kegg_inflam$Dataset <- "20H LPS"
prot_kegg_inflam$Dataset <- "20H LPS"

#Conbining the datasets
combined_go_RNA.treat_kegg.all1 <- rbind(rna_kegg_inflam, rna_kegg_inflam_3, rna_kegg_inflam_5, rna_kegg_inflam_1, rna_kegg_inflam_4)
combined_go_Pro.treat_kegg.all1 <- rbind(prot_kegg_inflam, prot_kegg_inflam_2H, prot_kegg_inflam_22H, prot_kegg_inflam_2H4ME, prot_kegg_inflam_22H4ME)

ggplot(combined_go_RNA.treat_kegg.all1, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related KEGG Terms: RNA-seq",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS" = "#F8766D", "20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF", "2H 4-ME" = "yellow2", "22H 4-ME" = "orange")) +
  theme_minimal(base_size = 14)

ggplot(combined_go_Pro.treat_kegg.all1, aes(x = reorder(Description, Count), y = Count, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Inflammation-Related KEGG Terms: Proteomics",
       x = "GO Term",
       y = "Gene Count") +
  scale_fill_manual(values = c("20H LPS" = "#F8766D", "20H LPS + 2H 4-ME" = "#00BA38", "20H LPS + 22H 4-ME" = "#619CFF", "2H 4-ME" = "yellow2", "22H 4-ME" = "orange")) +
  theme_minimal(base_size = 14)

#PATHVIEW (4-ME)----

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(dplyr)


#RNA-seq: res_1 and res_4
# res_1: 2H 4-ME
res1_df <- as.data.frame(res_1)
res1_df$gene <- rownames(res1_df)

rna1_map <- bitr(res1_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
rna1_fc <- merge(rna1_map, res1_df, by.x = "SYMBOL", by.y = "gene") %>%
  dplyr::select(ENTREZID, log2FoldChange)
rna1_vec <- setNames(rna1_fc$log2FoldChange, rna1_fc$ENTREZID)

# Pathview: MAPK and NF-kB
pathview(gene.data = rna1_vec, pathway.id = "hsa04010", species = "hsa", out.suffix = "RNAseq_res1_MAPK",
         kegg.native = TRUE, low = list(gene = "blue"), mid = list(gene = "white"), high = list(gene = "red"))

pathview(gene.data = rna1_vec, pathway.id = "hsa04064", species = "hsa", out.suffix = "RNAseq_res1_NFkB",
         kegg.native = TRUE, low = list(gene = "blue"), mid = list(gene = "white"), high = list(gene = "red"))


# res_4: 22H 4-ME
res4_df <- as.data.frame(res_4)
res4_df$gene <- rownames(res4_df)

rna4_map <- bitr(res4_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
rna4_fc <- merge(rna4_map, res4_df, by.x = "SYMBOL", by.y = "gene") %>%
  dplyr::select(ENTREZID, log2FoldChange)
rna4_vec <- setNames(rna4_fc$log2FoldChange, rna4_fc$ENTREZID)

# Pathview: MAPK and NF-kB
pathview(gene.data = rna4_vec, pathway.id = "hsa04010", species = "hsa", out.suffix = "RNAseq_res4_MAPK",
         kegg.native = TRUE, low = list(gene = "blue"), mid = list(gene = "white"), high = list(gene = "red"))

pathview(gene.data = rna4_vec, pathway.id = "hsa04064", species = "hsa", out.suffix = "RNAseq_res4_NFkB",
         kegg.native = TRUE, low = list(gene = "blue"), mid = list(gene = "white"), high = list(gene = "red"))

#Proteomics: 2h4ME and 22H4ME

# 2H 4-ME
prot_map <- bitr(filtered_df$`PG.Genes`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

prot_fc_2h <- merge(prot_map, filtered_df, by.x = "SYMBOL", by.y = "PG.Genes") %>%
  dplyr::select(ENTREZID, log2_FC_Ctrl_vs_2h4ME)
prot_vec_2h <- setNames(prot_fc_2h$log2_FC_Ctrl_vs_2h4ME, prot_fc_2h$ENTREZID)

pathview(gene.data = prot_vec_2h, pathway.id = "hsa04010", species = "hsa", out.suffix = "Proteomics_2h_MAPK",
         kegg.native = TRUE, low = list(gene = "blue"), mid = list(gene = "white"), high = list(gene = "red"))

pathview(gene.data = prot_vec_2h, pathway.id = "hsa04064", species = "hsa", out.suffix = "Proteomics_2h_NFkB",
         kegg.native = TRUE, low = list(gene = "blue"), mid = list(gene = "white"), high = list(gene = "red"))


# 22H 4-ME
prot_fc_22h <- merge(prot_map, filtered_df, by.x = "SYMBOL", by.y = "PG.Genes") %>%
  dplyr::select(ENTREZID, log2_FC_Ctrl_vs_22H4ME)
prot_vec_22h <- setNames(prot_fc_22h$log2_FC_Ctrl_vs_22H4ME, prot_fc_22h$ENTREZID)

pathview(gene.data = prot_vec_22h, pathway.id = "hsa04010", species = "hsa", out.suffix = "Proteomics_22h_MAPK",
         kegg.native = TRUE, low = list(gene = "blue"), mid = list(gene = "white"), high = list(gene = "red"))

pathview(gene.data = prot_vec_22h, pathway.id = "hsa04064", species = "hsa", out.suffix = "Proteomics_22h_NFkB",
         kegg.native = TRUE, low = list(gene = "blue"), mid = list(gene = "white"), high = list(gene = "red"))

#INTENSITY HD----

setwd("~/Desktop")
int_HD <- read.csv("intensity_HD.csv", stringsAsFactors = FALSE)

#select files and add intensity means 
HD_red_means <- c(
  "Control" = (7262.8569 + 7262.8569 + 9989.7640 + 10529.7312 + 10529.7312 + 17256.1207) / 6,
  "4-ME 22H" = (12023.6910 + 12023.6910 + 10371.0261 + 10952.8985 + 10952.8985 + 17453.5344) / 6,
  "4-ME 2H" = (16502.1418 + 16502.1418 + 13469.1759 + 12173.2681 + 12173.2681 + 15448.2405) / 6,
  "LPS 20H" = (9770.5316 + 9770.5316 + 7561.1165) / 3,
  "22H 4-ME + LPS 20H" = (12088.1662 + 12088.1662 + 14813.4276 + 13666.9273 + 13666.9273 + 12965.7606) / 6,
  "2H 4-ME + LPS 20H" = (11951.3899 + 11951.3899 + 13392.1977 + 12963.5993 + 12963.5993 + 15111.6291) / 6
)

# Convert to data frame
df_red <- data.frame(
  Condition = names(HD_red_means),
  Red_Intensity = as.numeric(HD_red_means)
)

#RED INTENSITY PLOT
ggplot(df_red, aes(x = Condition, y = Red_Intensity)) +
  geom_bar(stat = "identity", fill = "red") +
  theme_minimal() +
  labs(
    title = "Average Red Channel Intensity by Condition, High Desnity",
    x = "Condition",
    y = "Mean Red Intensity"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

HD_green_means <- c(
  "Control" = (5649.6868 + 5649.6868 + 14023.7545 + 12371.8533 + 12371.8533 + 19655.7245) / 6,
  "4-ME 22H" = (16523.2320 + 16523.2320 + 15007.2862 + 15404.1263 + 15404.1263 + 18543.3341) / 6,
  "4-ME 2H" = (20578.1609 + 20578.1609 + 18799.3130 + 16739.4258 + 16739.4258 + 19131.7861) / 6,
  "LPS 20H" = (13340.5572 + 13340.5572 + 11523.9905) / 3,
  "22H 4-ME + LPS 20H" = (17516.3832 + 17516.3832 + 19529.7208 + 15951.2168 + 15951.2168 + 14869.8047) / 6,
  "2H 4-ME + LPS 20H" = (18336.7314 + 18336.7314 + 21186.8029 + 19684.1078 + 19684.1078 + 21370.2398) / 6
)

#GREEN INTENSITY PLOT

#Convert to data frame
df_green <- data.frame(
  Condition = names(HD_green_means),
  Green_Intensity = as.numeric(HD_green_means)
)

ggplot(df_green, aes(x = Condition, y = Green_Intensity)) +
  geom_bar(stat = "identity", fill = "green") +
  theme_minimal() +
  labs(
    title = "Average Green Channel Intensity by Condition, High Density",
    x = "Condition",
    y = "Mean Green Intensity"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine into a data frame
df_combined <- data.frame(
  Condition = names(HD_red_means),
  Red = as.numeric(HD_red_means),
  Green = as.numeric(HD_green_means)
)

# Reshape for plotting
df_long <- pivot_longer(df_combined, cols = c(Red, Green),
                        names_to = "Channel", values_to = "Intensity")


# Combine both
df_long <- rbind(df_red, df_green)


# Plot combined bar chart
ggplot(df_long, aes(x = Condition, y = Intensity, fill = Channel)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Red" = "red", "Green" = "green")) +
  theme_minimal() +
  labs(
    title = "Red and Green Channel Intensities by Condition, High Density",
    x = "Condition",
    y = "Mean Intensity"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




library(ggplot2)
library(dplyr)
library(tidyr)

# Define raw intensity values
HD_data <- list(
  Control = list(
    Red = c(7262.8569, 7262.8569, 9989.7640, 10529.7312, 10529.7312, 17256.1207),
    Green = c(5649.6868, 5649.6868, 14023.7545, 12371.8533, 12371.8533, 19655.7245)
  ),
  "4-ME 22H" = list(
    Red = c(12023.6910, 12023.6910, 10371.0261, 10952.8985, 10952.8985, 17453.5344),
    Green = c(16523.2320, 16523.2320, 15007.2862, 15404.1263, 15404.1263, 18543.3341)
  ),
  "4-ME 2H" = list(
    Red = c(16502.1418, 16502.1418, 13469.1759, 12173.2681, 12173.2681, 15448.2405),
    Green = c(20578.1609, 20578.1609, 18799.3130, 16739.4258, 16739.4258, 19131.7861)
  ),
  "LPS 20H" = list(
    Red = c(9770.5316, 9770.5316, 7561.1165),
    Green = c(13340.5572, 13340.5572, 11523.9905)
  ),
  "22H 4-ME + LPS 20H" = list(
    Red = c(12088.1662, 12088.1662, 14813.4276, 13666.9273, 13666.9273, 12965.7606),
    Green = c(17516.3832, 17516.3832, 19529.7208, 15951.2168, 15951.2168, 14869.8047)
  ),
  "2H 4-ME + LPS 20H" = list(
    Red = c(11951.3899, 11951.3899, 13392.1977, 12963.5993, 12963.5993, 15111.6291),
    Green = c(18336.7314, 18336.7314, 21186.8029, 19684.1078, 19684.1078, 21370.2398)
  )
)

# Summarize into a long-format data frame with means and SE
df_long <- lapply(names(HD_data), function(cond) {
  red_vals <- HD_data[[cond]]$Red
  green_vals <- HD_data[[cond]]$Green
  data.frame(
    Condition = cond,
    Channel = c("Red", "Green"),
    Intensity = c(mean(red_vals), mean(green_vals)),
    SE = c(sd(red_vals)/sqrt(length(red_vals)), sd(green_vals)/sqrt(length(green_vals)))
  )
}) %>% bind_rows()

# Plot combined bar chart with error bars
ggplot(df_long, aes(x = Condition, y = Intensity, fill = Channel)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(
    aes(ymin = Intensity - SE, ymax = Intensity + SE),
    width = 0.3,
    position = position_dodge(0.9)
  ) +
  scale_fill_manual(values = c("Red" = "red", "Green" = "green")) +
  theme_minimal() +
  labs(
    title = "Red and Green Channel Intensities by Condition, High Density",
    x = "Condition",
    y = "Mean Intensity ± SE"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#INTENSITY LD----

setwd("~/Desktop")
int_LD <- read.csv("intensity_LD.csv", stringsAsFactors = FALSE)

LD_red_means <- c(
  "Control" = (4705.8894 + 4705.8894 + 3990.9677 + 5683.0875 + 5683.0875 + 7038.0682) / 6,
  "4-ME 22H" = (6516.7016 + 6516.7016 + 11246.5862) / 3,
  "4-ME 2H" = (3908.6444 + 3908.6444 + 9542.5515 + 3818.3015 + 3818.3015 + 3703.0520) / 6,
  "LPS 20H" = (7307.8814 + 7307.8814 + 998.4332 + 4330.3779 + 4330.3779 + 6521.0375) / 6,
  "22H 4-ME + LPS 20H" = (8868.0156 + 8868.0156 + 10930.3169 + 11824.6885 + 11824.6885 + 11129.8793) / 6,
  "2H 4-ME + LPS 20H" = (8975.9823 + 8975.9823 + 4953.2241) / 3
)

LD_green_means <- c(
  "Control" = (3972.3071 + 3972.3071 + 3528.1796 + 6904.5483 + 6904.5483 + 8040.5768) / 6,
  "4-ME 22H" = (9287.3606 + 9287.3606 + 15705.9303) / 3,
  "4-ME 2H" = (3095.2556 + 3095.2556 + 5034.2518 + 1491.6572 + 1491.6572 + 2195.2078) / 6,
  "LPS 20H" = (7437.9369 + 7437.9369 + 1022.7328 + 5132.9478 + 5132.9478 + 5089.8885) / 6,
  "22H 4-ME + LPS 20H" = (11857.9149 + 11857.9149 + 13031.0167 + 13593.8271 + 13593.8271 + 14167.0661) / 6,
  "2H 4-ME + LPS 20H" = (14088.9713 + 14088.9713 + 8122.5171) / 3
)

# Combine into a data frame
df_combined_LD <- data.frame(
  Condition = names(LD_red_means),
  Red = as.numeric(LD_red_means),
  Green = as.numeric(LD_green_means)
)

# Reshape for plotting
df_long_LD <- pivot_longer(df_combined_LD, cols = c(Red, Green),
                        names_to = "Channel", values_to = "Intensity")

# Plot combined bar chart
ggplot(df_long_LD, aes(x = Condition, y = Intensity, fill = Channel)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Red" = "red", "Green" = "green")) +
  theme_minimal() +
  labs(
    title = "Red and Green Channel Intensities by Condition, Low Density",
    x = "Condition",
    y = "Mean Intensity"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

LD_data <- list(
  Control = list(
    Red = c(4705.8894, 4705.8894, 3990.9677, 5683.0875, 5683.0875, 7038.0682),
    Green = c(3972.3071, 3972.3071, 3528.1796, 6904.5483, 6904.5483, 8040.5768)
  ),
  "4-ME 22H" = list(
    Red = c(6516.7016, 6516.7016, 11246.5862),
    Green = c(9287.3606, 9287.3606, 15705.9303)
  ),
  "4-ME 2H" = list(
    Red = c(3908.6444, 3908.6444, 9542.5515, 3818.3015, 3818.3015, 3703.0520),
    Green = c(3095.2556, 3095.2556, 5034.2518, 1491.6572, 1491.6572, 2195.2078)
  ),
  "LPS 20H" = list(
    Red = c(7307.8814, 7307.8814, 998.4332, 4330.3779, 4330.3779, 6521.0375),
    Green = c(7437.9369, 7437.9369, 1022.7328, 5132.9478, 5132.9478, 5089.8885)
  ),
  "22H 4-ME + LPS 20H" = list(
    Red = c(8868.0156, 8868.0156, 10930.3169, 11824.6885, 11824.6885, 11129.8793),
    Green = c(11857.9149, 11857.9149, 13031.0167, 13593.8271, 13593.8271, 14167.0661)
  ),
  "2H 4-ME + LPS 20H" = list(
    Red = c(8975.9823, 8975.9823, 4953.2241),
    Green = c(14088.9713, 14088.9713, 8122.5171)
  )
)

library(dplyr)

df_long_LD <- lapply(names(LD_data), function(cond) {
  red_vals <- LD_data[[cond]]$Red
  green_vals <- LD_data[[cond]]$Green
  data.frame(
    Condition = cond,
    Channel = c("Red", "Green"),
    Intensity = c(mean(red_vals), mean(green_vals)),
    SE = c(sd(red_vals) / sqrt(length(red_vals)), sd(green_vals) / sqrt(length(green_vals)))
  )
}) %>% bind_rows()

library(ggplot2)

ggplot(df_long_LD, aes(x = Condition, y = Intensity, fill = Channel)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(
    aes(ymin = Intensity - SE, ymax = Intensity + SE),
    width = 0.3,
    position = position_dodge(0.9)
  ) +
  scale_fill_manual(values = c("Red" = "red", "Green" = "green")) +
  theme_minimal() +
  labs(
    title = "Red and Green Channel Intensities by Condition, Low Density",
    x = "Condition",
    y = "Mean Intensity ± SE"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


