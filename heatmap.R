library(DESeq2)
library(pheatmap)

counts <- read.csv("cleaned_counts.csv", row.names = 1)
metadata <- read.csv("metadata_counts.csv", row.names = 1)
colnames(counts) <- rownames(metadata)
metadata$Conditions <- as.factor(metadata$Conditions)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Conditions)
dds <- DESeq(dds)
res <- results(dds)

topGenes <- rownames(res[order(res$padj),])[1:100]
vsd <- vst(dds)

pheatmap(assay(vsd)[topGenes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE)
