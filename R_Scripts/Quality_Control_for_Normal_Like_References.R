library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)

if(!dir.exists("Results_QC")) dir.create("Results_QC")

gtex_samples <- meta_final$sample[grepl("^GTEX", meta_final$sample)]
# 计算 GTEx 的平均表达量
gtex_mean_profile <- rowMeans(expr_final[, gtex_samples], na.rm = TRUE)

adj_samples <- meta_final$sample[meta_final$`_sample_type` == "Solid Tissue Normal"]
adj_mat <- expr_final[, adj_samples]

# 计算 Spearman 相关性
cor_scores <- cor(adj_mat, gtex_mean_profile, method = "spearman", use = "pairwise.complete.obs")
cor_df <- data.frame(sample = rownames(cor_scores), correlation = as.numeric(cor_scores))

mu <- mean(cor_df$correlation, na.rm = TRUE)
sigma <- sd(cor_df$correlation, na.rm = TRUE)

# 阈值定义平均值下浮 2 倍标准差
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

normal_like_samples <- meta_final %>% filter(sample %in% kept_ids)

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

tumor_samples <- meta_final$sample[meta_final$`_sample_type` == "Primary Tumor"]

if(length(tumor_samples) > 500) {
  plot_tumor_ids <- sample(tumor_samples, 500)
} else {
  plot_tumor_ids <- tumor_samples
}

tumor_mat <- expr_final[, plot_tumor_ids]
tumor_cor <- cor(tumor_mat, gtex_mean_profile, method = "spearman", use = "pairwise.complete.obs")

# 合并所有绘图数据
plot_data <- rbind(
  data.frame(Group = "Selected Adjacent", Correlation = cor_df$correlation[cor_df$sample %in% kept_ids]),
  data.frame(Group = "Rejected Adjacent", Correlation = cor_df$correlation[cor_df$sample %in% rejected_ids]),
  data.frame(Group = "TCGA Tumor", Correlation = as.numeric(tumor_cor))
)

plot_data$Group <- factor(plot_data$Group, levels = c("Selected Adjacent", "Rejected Adjacent", "TCGA Tumor"))

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
