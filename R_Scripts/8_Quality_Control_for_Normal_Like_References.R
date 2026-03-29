library(data.table)
library(dplyr)      
library(ggplot2)    
library(ggpubr)     
library(tibble)     

# 设置随机种子
set.seed(42) 

if(!dir.exists("Results_QC")) dir.create("Results_QC")

# 全局路径与变量设置 
# ------------------------------------------------------------------------------
CANCER_CODE <- "BRCA"
PRIMARY_SITE <- "Breast"

TPM_DATA_PATH <- "Data_Source/TcgaTargetGtex_rsem_gene_tpm/TcgaTargetGtex_rsem_gene_tpm.tsv"
METADATA_PATH <- "Data_Source/TcgaTargetGTEX_phenotype.txt/TcgaTargetGTEX_phenotype.txt"

# ==============================================================================
# 加载数据与预处理
# ==============================================================================
message(">>> 开始加载数据...")

if (file.exists(METADATA_PATH)) {
  phenotype <- fread(METADATA_PATH, data.table = FALSE)
  message("  - 元数据加载完成。")
} else {
  stop("错误: 找不到元数据文件!")
}

# 读取表达矩阵
if (file.exists(TPM_DATA_PATH)) {
  expression_matrix <- fread(TPM_DATA_PATH, data.table = FALSE)
  rownames(expression_matrix) <- expression_matrix[, 1]
  expression_matrix <- expression_matrix[, -1] 
  message("  - 表达矩阵加载完成。")
} else {
  stop("错误: 找不到表达矩阵文件!")
}

# ==============================================================================
# 数据筛选与对齐
# ==============================================================================
message(paste0(">>> 正在筛选 ", CANCER_CODE, " 和 GTEx ", PRIMARY_SITE, " 样本..."))

# 筛选元数据
target_samples <- phenotype %>%
  filter(`_primary_site` == PRIMARY_SITE) %>%
  filter(
    (`_sample_type` %in% c("Primary Tumor", "Solid Tissue Normal") & grepl("^TCGA", sample)) |
      (grepl("^GTEX", sample))
  )

if (nrow(target_samples) == 0) stop("未找到符合条件的样本，请检查筛选条件。")

# 3.2 数据对齐 (取交集)
common_samples <- intersect(target_samples$sample, colnames(expression_matrix))
message(paste0("  - 数据对齐完成。共找到 ", length(common_samples), " 个共有样本。"))

# 提取对应的表达数据和元数据，并确保顺序一致
expr_final <- expression_matrix[, common_samples]
meta_final <- target_samples %>% 
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples)) 

# ==============================================================================
# 基于 Spearman 相关性的 Normal-like 样本筛选 
# ==============================================================================
message(">>> 正在构建 GTEx 标准参考谱并计算相关性...")

gtex_samples <- meta_final$sample[grepl("^GTEX", meta_final$sample)]
gtex_mean_profile <- rowMeans(expr_final[, gtex_samples], na.rm = TRUE)

# 计算 TCGA 癌旁样本的相关性
adj_samples <- meta_final$sample[meta_final$`_sample_type` == "Solid Tissue Normal"]
adj_mat <- expr_final[, adj_samples]
cor_scores <- cor(adj_mat, gtex_mean_profile, method = "spearman", use = "pairwise.complete.obs")
cor_df <- data.frame(sample = rownames(cor_scores), correlation = as.numeric(cor_scores))

# 定义统计学阈值 (Mean - 2SD)
mu <- mean(cor_df$correlation, na.rm = TRUE)
sigma <- sd(cor_df$correlation, na.rm = TRUE)
cutoff_value <- mu - 2 * sigma 

# 执行筛选
kept_ids <- cor_df$sample[cor_df$correlation >= cutoff_value]
rejected_ids <- cor_df$sample[cor_df$correlation < cutoff_value]

# 打印统计信息
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

# 更新全局变量
normal_like_samples <- meta_final %>% filter(sample %in% kept_ids)

# ==============================================================================
# 导出结果
# ==============================================================================
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

# ==============================================================================
# 绘制验证图
# ==============================================================================
message(">>> 正在绘制验证图...")

# 准备 Tumor 对比数据
tumor_samples <- meta_final$sample[meta_final$`_sample_type` == "Primary Tumor"]

if(length(tumor_samples) > 500) {
  plot_tumor_ids <- sample(tumor_samples, 500)
} else {
  plot_tumor_ids <- tumor_samples
}

tumor_mat <- expr_final[, plot_tumor_ids]
tumor_cor <- cor(tumor_mat, gtex_mean_profile, method = "spearman", use = "pairwise.complete.obs")

# 合并绘图数据
plot_data <- rbind(
  data.frame(Group = "Selected Adjacent", Correlation = cor_df$correlation[cor_df$sample %in% kept_ids]),
  data.frame(Group = "Rejected Adjacent", Correlation = cor_df$correlation[cor_df$sample %in% rejected_ids]),
  data.frame(Group = "TCGA Tumor", Correlation = as.numeric(tumor_cor))
)

plot_data$Group <- factor(plot_data$Group, levels = c("Selected Adjacent", "Rejected Adjacent", "TCGA Tumor"))

# 绘图
box_colors <- c(
  "Selected Adjacent" = "#E9C46A", 
  "Rejected Adjacent" = "gray30",
  "TCGA Tumor" = "#E63946"         
)

p_valid <- ggplot(plot_data, aes(x = Group, y = Correlation, fill = Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, size = 1, color="#B0B0B0", alpha=0.2) +
  scale_fill_manual(values = box_colors) +
  
  # 统计检验
  stat_compare_means(comparisons = list(
    c("Selected Adjacent", "Rejected Adjacent"),
    c("Selected Adjacent", "TCGA Tumor") 
  ), method = "wilcox.test", label = "p.signif", size = 6, vjust = 0.5) +
  
  # 标签与修饰
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

# 打印并保存
print(p_valid)
ggsave(paste0("Results_QC/Correlation_Validation_Final_", CANCER_CODE, ".png"), p_valid, width = 5, height = 6, dpi = 300)
ggsave(paste0("Results_QC/Correlation_Validation_Final_", CANCER_CODE, ".pdf"), p_valid, width = 5, height = 6)
