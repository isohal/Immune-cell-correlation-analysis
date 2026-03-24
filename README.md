# Immune-cell-correlation-analysis
This page describes the code used to conduct immune cell correlation analysis for the following manuscript:

## Exosome-like biogenesis from the Golgi releases extracellular vesicles lacking conventional tetraspanins that mediate immune evasion in cancer
Ikjot S. Sohal<sup>1,2,3</sup>, Sydney N. Shaw<sup>1</sup>, Andrea P. dos Santos<sup>2,5</sup>, Bennett D. Elzey<sup>2,4,5</sup>, Haley Anne Harper<sup>2</sup>, Carlie McMahan<sup>2</sup>, Lauren N. Meeks<sup>1</sup>, Humna Hasan<sup>1</sup>, Zulaida Soto-Vargas<sup>1</sup>, Noor Abdullah<sup>1</sup>, Subhransu Sekhar Sahoo<sup>6</sup>, Majid Kazemian<sup>2,6,7</sup>, Matthew R. Olson<sup>1,2</sup>, Andrea L. Kasinski<sup>1,2*</sup>

<sup>1</sup>1Department of Biological Sciences, Purdue University, West Lafayette, IN 47907, USA 

<sup>2</sup>Purdue Institute for Cancer Research, Purdue University, West Lafayette, IN 47907, USA 

<sup>3</sup>Department of Biological Sciences, University of North Texas, Denton, TX 76201, USA 

<sup>4</sup>Purdue Institute of Inflammation, Immunology and Infectious Disease, Purdue University, West Lafayette, IN 47907 USA 

<sup>5</sup>Department of Comparative Pathobiology, College of Veterinary Medicine, Purdue University, West Lafayette, IN 47907, USA 

<sup>6</sup>Department of Biochemistry, Purdue University, West Lafayette, IN 47907, USA 

<sup>7</sup>Department of Computer Science, Purdue University, West Lafayette, IN 47907, USA 

*Corresponding author: Andrea L. Kasinski – akasinski@purdue.edu 

The authors declare no potential conflicts of interest

## R script for immune cell correlation heatmap data extraction
The following script describes the steps to extract data from immune cell correlation heatmap from the TCGA (The Cancer Genome Atlas) data:

```
library(TCGAplot)
```

### 1. Define your gene or geneset

    Gene would be "TGOLN2" and geneset would be "top 25 proteins" enriched in TGOLN2+ EVs or CD81+ EVs. The below example uses "TGOLN2" as an example:
```
gene_set <- "TGOLN2"
```

### 2. Get all TCGA cancer types
```
all_cancer_types = rownames(get_cancers())
cancer_types = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
```

### 3. Get immune cell ratios across all cancers
```
immu_ratio_all = get_immu_ratio()
```

### 4. Initialize results storage
```
immune_cells <- colnames(immu_ratio_all)
cor_matrix <- matrix(NA, nrow = length(immune_cells), ncol = length(all_cancer_types), dimnames = list(immune_cells, all_cancer_types))

pval_matrix <- matrix(NA, nrow = length(immune_cells), ncol = length(all_cancer_types), dimnames = list(immune_cells, all_cancer_types))

sample_count <- setNames(rep(NA, length(all_cancer_types)), all_cancer_types)
```

### 5. Loop over cancers
```
for (ct in all_cancer_types) { 
  message("Processing: ", ct) 
  
  # TPM for this cancer 
  tpm <- get_tpm(ct) 
  
  # ---- FILTER TO TUMOR SAMPLES ONLY ----
  tumor_samples <- rownames(tpm)[which(tpm$Group == "Tumor")]
  tpm <- tpm[tumor_samples, , drop = FALSE]
  sample_count[ct] <- length(tumor_samples)
  
  # Skip if too few tumor samples 
  if (ncol(tpm) < 10) next
  
  # Gene-set score 
  common_genes <- intersect(gene_set, colnames(tpm)) 
  if (length(common_genes) == 0) next 
  
  gs_score <- rowMeans(tpm[,common_genes, drop = FALSE])
  
  # Align tumor samples with immune ratios
  common_samples <- intersect(tumor_samples, rownames(immu_ratio_all)) 
  if (length(common_samples) < 10) next
  
  immu_ratio <- immu_ratio_all[common_samples, ]
  
  # Compute correlations for all immune cell types
  for (cell in immune_cells) { 
    test <- cor.test(gs_score, immu_ratio[, cell], method = "pearson")
  
    # Store in matrix  
    cor_matrix[cell, ct] <- test$estimate 
    pval_matrix[cell, ct] <- test$p.value
  }
}
```

### 6. Inspect results & save
```
head(cor_matrix) 
head(pval_matrix)

write.csv(cor_matrix, file = "CD81ev_top25_cor_matrix.csv")
write.csv(pval_matrix, file = "CD81ev_top25_pval_matrix.csv")
write.csv(sample_count, file = "CD81ev_top25_sample_count.csv")
```

### 7. Generate correlation matrix
```
library(dplyr)
library(tidyr)

df <- cor_matrix %>%
  as.data.frame() %>%
  mutate(Immune_Cell = rownames(.)) %>%
  pivot_longer(-Immune_Cell, names_to = "Cancer", values_to = "Correlation") %>%
  left_join(
    pval_matrix %>%
      as.data.frame() %>%
      mutate(Immune_Cell = rownames(.)) %>%
      pivot_longer(-Immune_Cell, names_to = "Cancer", values_to = "Pvalue"),
    by = c("Immune_Cell", "Cancer")
  ) %>%
  mutate(
    SampleSize = sample_count[Cancer],
    Significant = Pvalue < 0.05
  )
```

### 8. Heirarchical clustering of the correlation matrix
```
row_dist  <- dist(cor_matrix)              # distance between immune cells
row_clust <- hclust(row_dist, method = "ward.D2")

col_dist  <- dist(t(cor_matrix))           # distance between cancers
col_clust <- hclust(col_dist, method = "ward.D2")
```

### 9. Generate heatmap
```
library(pheatmap)
library(RColorBrewer)

sig_matrix <- ifelse(pval_matrix < 0.001, "***",
                     ifelse(pval_matrix < 0.01, "**",
                            ifelse(pval_matrix < 0.05, "*", "")))

ph = pheatmap(cor_matrix,
              scale = "column",
              cluster_rows = F,
              clustering_distance_rows = "manhattan",
              cluster_cols = F,
              clustering_distance_cols = "manhattan",
              clustering_method = "ward.D2",
              color = colorRampPalette(rev(c("red","white","blue")))(256), 
              #colorRampPalette(rev(brewer.pal(8,'RdBu')))(256),
              show_rownames = T,
              show_colnames = T,
              kmeans_k = NA,
              cutree_rows = 2,
              cutree_cols = 2,
              silent = F,
              fontsize = 11,
              legend = T,
              treeheight_row = 0,
              annotation_legend = T,
              display_numbers = sig_matrix,
              na_col = "white",
              border_color = "white")
```
