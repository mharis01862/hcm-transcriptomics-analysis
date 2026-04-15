library(ggplot2)
library(ggrepel)
library(dplyr)

res_df <- read.csv("DESeq2_results.csv", row.names = 1)
res_df$Significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 3, 
                             ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"), 
                             "Not Significant")

res_df$Shape <- ifelse(res_df$Significance == "Upregulated", 24,
                       ifelse(res_df$Significance == "Downregulated", 22, 21))

top_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 3) %>%
  arrange(padj) %>%
  slice(18:39)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj),
                   color = Significance, shape = as.factor(Shape))) +
  geom_point(alpha = 0.8, size = 1) +
  theme_minimal() +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "black")) +
  scale_shape_manual(values = c("21" = 21, "22" = 22, "24" = 24)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_genes,
                  aes(label = rownames(top_genes)),
                  size = 2) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")
