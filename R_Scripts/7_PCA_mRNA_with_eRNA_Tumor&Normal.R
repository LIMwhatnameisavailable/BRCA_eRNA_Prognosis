library(data.table)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggsci)

# 数据读取
eRNA_df <- fread("Data_Source/TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv", data.table = FALSE)
rownames(eRNA_df) <- eRNA_df[,1]; eRNA_df <- eRNA_df[,-1]

mRNA_df <- fread("Data_Source/TCGA-BRCA.star_fpkm.tsv.gz", data.table = FALSE)
rownames(mRNA_df) <- mRNA_df[,1]; mRNA_df <- mRNA_df[,-1]

# 清洗 eRNA 
eRNA_raw_cols <- colnames(eRNA_df)
eRNA_group_vec <- rep("Unknown", length(eRNA_raw_cols))
eRNA_group_vec[grep("_normal", eRNA_raw_cols, ignore.case = TRUE)] <- "Normal"
eRNA_group_vec[grep("_tumor", eRNA_raw_cols, ignore.case = TRUE)] <- "Tumor"
eRNA_group_vec[grep("_tu", eRNA_raw_cols, ignore.case = TRUE)] <- "Tumor"

clean_eRNA_cols <- gsub("_normal|_tumor|_tu", "", eRNA_raw_cols, ignore.case = TRUE)
clean_eRNA_cols <- gsub("\\.", "-", clean_eRNA_cols)
clean_eRNA_cols <- substr(clean_eRNA_cols, 1, 12)
colnames(eRNA_df) <- clean_eRNA_cols

# 清洗 mRNA
clean_mRNA_cols <- substr(colnames(mRNA_df), 1, 12)
colnames(mRNA_df) <- clean_mRNA_cols

# 匹配与提取最终变量
common_patients <- intersect(clean_eRNA_cols, clean_mRNA_cols)
eRNA_idx <- match(common_patients, colnames(eRNA_df))
mRNA_idx <- match(common_patients, colnames(mRNA_df))

final_eRNA <- eRNA_df[, eRNA_idx]
final_mRNA <- mRNA_df[, mRNA_idx]
final_group <- eRNA_group_vec[eRNA_idx]

cat("   - 样本数:", length(final_group), "\n")
cat("   - Normal:", sum(final_group == "Normal"), "\n")
cat("   - Tumor: ", sum(final_group == "Tumor"), "\n")

# 绘图函数 
run_pca_plot <- function(expr_mat, group_vec, title_text) {
  # 数据计算
  expr_mat <- log2(as.matrix(expr_mat) + 1)
  vars <- apply(expr_mat, 1, var, na.rm = TRUE)
  keep_idx <- which(vars > 0 & !is.na(vars))
  expr_mat <- expr_mat[keep_idx, ]
  vars <- vars[keep_idx]
  
  ntop <- min(5000, nrow(expr_mat))
  select <- order(vars, decreasing = TRUE)[1:ntop]
  
  pca_input <- t(expr_mat[select, ])
  # 移除极低方差列
  pca_input <- pca_input[, apply(pca_input, 2, var) > 1e-10]
  
  pca <- prcomp(pca_input, scale. = TRUE)
  pca_data <- as.data.frame(pca$x)
  pca_data$Group <- group_vec
  
  var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, fill = Group)) +
    geom_point(shape = 21, size = 3, color = "white", stroke = 0.2, alpha = 0.8) +
    scale_fill_manual(values = c("Normal" = "#457B9D", "Tumor" = "#E63946")) + 
    labs(
      title = title_text,
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme_classic(base_size = 14) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "top",
      legend.title = element_blank(),
      axis.text = element_text(color = "black"),
      aspect.ratio = 1 
    )
  return(p)
}

# 生成图片
p_eRNA <- run_pca_plot(final_eRNA, final_group, "eRNA (Normal vs Tumor)")
p_mRNA <- run_pca_plot(final_mRNA, final_group, "mRNA (Normal vs Tumor)")

combined_plot <- grid.arrange(p_mRNA, p_eRNA, ncol = 2)

ggsave("Results/Fig_PCA_mRNA&eRNA.svg", combined_plot, width = 12, height = 6, device = "svg")
cat(" 图片已保存为: Fig_PCA_mRNA&eRNA.svg")
