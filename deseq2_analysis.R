library(DESeq2)

counts <- read.csv("cleaned_counts.csv", row.names = 1)
metadata <- read.csv("metadata_counts.csv", row.names = 1)

colnames(counts) <- rownames(metadata)
metadata$Conditions <- as.factor(metadata$Conditions)

counts <- round(counts)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Conditions)

dds <- DESeq(dds)

res <- results(dds, contrast = c("Conditions", "Diseased", "Control"))

res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 3)
write.csv(as.data.frame(res), "DESeq2_results.csv")
write.csv(as.data.frame(res_sig), "DESeq2_significant_genes.csv")
