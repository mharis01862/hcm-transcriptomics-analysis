# Hypertrophic Cardiomyopathy (HCM) Transcriptomic Analysis
Integrative Transcriptomic and Network-Based Analysis of Hypertrophic Cardiomyopathy for Mechanistic Insights and Drug Repurposing

# Overview
Hypertrophic Cardiomyopathy (HCM) is the most prevalent monogenic cardiac disorder, with no definitive cure and treatments limited to symptomatic management. This project presents a comprehensive RNA-seq–based analysis integrating transcriptomics, network biology, and computational drug repurposing to uncover key molecular mechanisms and identify potential therapeutic candidates.

# Objectives
- Identify differentially expressed genes (DEGs) in HCM myocardial tissue
- Elucidate dysregulated biological pathways
- Construct protein–protein interaction (PPI) networks and identify hub genes
- Explore potential drug repurposing candidates

# Dataset
- Source: NCBI Sequence Read Archive (SRA)
- BioProject ID: PRJNA961804
- Sample Type: Human myocardial tissue

# Methodology
RNA-seq Processing:
- Alignment: Bowtie2
- Quantification: featureCounts
- Reference Genome: GRCh38

# Differential Expression Analysis:
- Tool: DESeq2
- Thresholds: log2FC > 3, padj < 0.05

# Network Analysis:
- STRING database
- Cytoscape (CytoHubba)

# Functional Enrichment:
- ShinyGO
- Enrichr
- KEGG pathways

# Drug Repurposing:
- Computational identification of candidate compounds

# Key Results
- 86 high-confidence DEGs identified

# Dysregulated Pathways:
- Mineral absorption (zinc homeostasis disruption)
- NET formation (immune dysregulation)

# Hub Genes:
FPR1, FPR2, CCR1, SSTR5, EDN2, EDN3, S100A8, S100A9, TLR2

# Drug Candidates:
- Zinc chloride, Zinc acetate
- Rebamipide, Nedocromil

# Biological Significance
Highlights the role of zinc homeostasis and innate immune signaling in HCM and suggests potential therapeutic strategies.

# Reproducibility
Clone the repository, install dependencies, and run scripts.

# Author
Muhammad Haris

# Disclaimer
This study is computational and requires experimental validation.

