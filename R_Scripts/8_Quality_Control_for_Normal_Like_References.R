library(data.table)
library(dplyr)      
library(ggplot2)    
library(ggpubr)     
library(tibble)     

# Set random seed
set.seed(42) 

if(!dir.exists("Results_QC")) dir.create("Results_QC")

# Global paths and variables
CANCER_CODE <- "BRCA"
PRIMARY_SITE <- "Breast"

TPM_DATA_PATH <- "Data_Source/TcgaTargetGtex_rsem_gene_tpm/TcgaTargetGtex_rsem_gene_tpm.tsv"
METADATA_PATH <- "Data_Source/TcgaTargetGTEX_phenotype.txt/TcgaTargetGTEX_phenotype.txt"

# Data loading and preprocessing
message(">>> 开始加载数据...")

if (file.exists(METADATA_PATH)) {
  phenotype <- fread(METADATA_PATH, data.table = FALSE)
  message("  - 元数据加载完成。")
} else {
  stop("错误: 找不到元数据文件!")
}

# Load expression matrix
if (file.exists(TPM_DATA_PATH)) {
  expression_matrix <- fread(TPM_DATA_PATH, data.table = FALSE)
  rownames(expression_matrix) <- expression_matrix[, 1]
  expression_matrix <- expression_matrix[, -1] 
  message("  - 表达矩阵加载完成。")
} else {
  stop("错误: 找不到表达矩阵文件!")
}

# Data filtering and alignment
message(paste0(">>> 正在筛选 ", CANCER_CODE, " 和 GTEx ", PRIMARY_SITE, " 样本..."))

# Filter metadata
target_samples <- phenotype %>%
  filter(`_primary_site` == PRIMARY_SITE) %>%
  filter(
    (`_sample_type` %in% c("Primary Tumor", "Solid Tissue Normal") & grepl("^TCGA", sample)) |
      (grepl("^GTEX", sample))
  )

if (nrow(target_samples) == 0) stop("未找到符合条件的样本，请检查筛选条件。")

# Align data via intersection
common_samples <- intersect(target_samples$sample, colnames(expression_matrix))
message(paste0("  - 数据对齐完成。共找到 ", length(common_samples), " 个共有样本。"))

# Extract corresponding data and ensure consistent ordering
expr_final <- expression_matrix[, common_samples]
meta_final <- target_samples %>% 
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples)) 

# Filter normal-like samples via Spearman correlation
message(">>> 正在构建 GTEx 标准参考谱并计算相关性...")

gtex_samples <- meta_final$sample[grepl("^GTEX", meta_final$sample)]
gtex_mean_profile <- rowMeans(expr_final[, gtex_samples], na.rm = TRUE)

# Calculate correlation of TCGA adjacent normal samples
adj_samples <- meta_final$sample[meta_final$`_sample_type` == "Solid Tissue Normal"]
adj_mat <- expr_final[, adj_samples]
cor_scores <- cor(adj_mat, gtex_mean_profile, method = "spearman", use = "pairwise.complete.obs")
cor_df <- data.frame(sample = rownames(cor_scores), correlation = as.numeric(cor_scores))

# Define statistical threshold (Mean - 2SD)
mu <- mean(cor_df$correlation, na.rm = TRUE)
sigma <- sd(cor_df$correlation, na.rm = TRUE)
cutoff_value <- mu - 2 * sigma 

# Perform filtering
kept_ids <- cor_df$sample[cor_df$correlation >= cutoff_value]
rejected_ids <- cor_df$sample[cor_df$correlation < cutoff_value]

# Print statistical summary
stats_msg <- paste0(
  "--------------------------------------------------\n",
  "筛选统计报告:\n",
  "  - 癌旁样本总数: ", nrow(cor_df), "\n",
  "  - 平均相关性 (Mean): ", round(mu, 4), "\n",
  "  - 标准差 (SD): ", round(sigma, 4), "\n",
  "  - 设定阈值 (Mean - 2SD): ", round(cutoff_value, 4), "\n",
  "  - 保留样本数 (Selected): ", length(kept_ids), " (", round(length(kept_ids)/nrow(cor_df)*100, 1), "%)\n",
  "  - 剔除样本数 (Rejected): ", length(rejected_ids), " (", round(length(rejected_ids)/nrow(cor_df)*100, 1), "%)\n",
  "--------------------------------------------------"
)
message(stats_msg)

# Update global variables
normal_like_samples <- meta_final %>% filter(sample %in% kept_ids)

# Export results
message(">>> 正在导出筛选结果...")

export_df <- cor_df %>%
  filter(sample %in% kept_ids) %>%
  arrange(desc(correlation)) %>%
  rename(Sample_ID = sample, Spearman_Cor_vs_GTEx = correlation)

write.csv(export_df, 
          file = paste0("Results_QC/Selected_Normal_Samples_", CANCER_CODE, ".csv"), 
          row.names = FALSE, quote = FALSE)

write.table(kept_ids, 
            file = paste0("Results_QC/Selected_Normal_Samples_ID_Only_", CANCER_CODE, ".txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Generate validation plots
message(">>> 正在绘制验证图...")

# Prepare tumor data for comparison
tumor_samples <- meta_final$sample[meta_final$`_sample_type` == "Primary Tumor"]

if(length(tumor_samples) > 500) {
  plot_tumor_ids <- sample(tumor_samples, 500)
} else {
  plot_tumor_ids <- tumor_samples
}

tumor_mat <- expr_final[, plot_tumor_ids]
tumor_cor <- cor(tumor_mat, gtex_mean_profile, method = "spearman", use = "pairwise.complete.obs")

# Merge plotting data
plot_data <- rbind(
  data.frame(Group = "Selected Adjacent", Correlation = cor_df$correlation[cor_df$sample %in% kept_ids]),
  data.frame(Group = "Rejected Adjacent", Correlation = cor_df$correlation[cor_df$sample %in% rejected_ids]),
  data.frame(Group = "TCGA Tumor", Correlation = as.numeric(tumor_cor))
)

plot_data$Group <- factor(plot_data$Group, levels = c("Selected Adjacent", "Rejected Adjacent", "TCGA Tumor"))

# Generate plot
box_colors <- c(
  "Selected Adjacent" = "#E9C46A", 
  "Rejected Adjacent" = "gray30",
  "TCGA Tumor" = "#E63946"         
)

p_valid <- ggplot(plot_data, aes(x = Group, y = Correlation, fill = Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, size = 1, color="#B0B0B0", alpha=0.2) +
  scale_fill_manual(values = box_colors) +
  
  # Statistical testing
  stat_compare_means(comparisons = list(
    c("Selected Adjacent", "Rejected Adjacent"),
    c("Selected Adjacent", "TCGA Tumor") 
  ), method = "wilcox.test", label = "p.signif", size = 6, vjust = 0.5) +
  
  # Plot labels and aesthetics
  labs(
    title = "Validation of Sample Selection Strategy",
    subtitle = paste0("Statistical Cutoff: Mean - 2SD (Rho < ", round(cutoff_value, 3), ")"),
    y = "Spearman Correlation with GTEx Profile",
    x = ""
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30")
  )

# Print and save plots
print(p_valid)
ggsave(paste0("Results_QC/Correlation_Validation_Final_", CANCER_CODE, ".png"), p_valid, width = 5, height = 6, dpi = 300)
ggsave(paste0("Results_QC/Correlation_Validation_Final_", CANCER_CODE, ".pdf"), p_valid, width = 5, height = 6)
