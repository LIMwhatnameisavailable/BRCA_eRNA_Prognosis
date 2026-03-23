# 准备上传到 L1000FWD 网站的基因列表
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 检查前序数据
if (!exists("all_genes_diff")) {
  stop("错误：未找到 'all_genes_diff'。请先加载 04_GSEA_Analysis.R 的结果！")
}

# 提取 Top 150 上调和下调基因
deg_sorted <- all_genes_diff %>% 
  tibble::rownames_to_column("Symbol") %>%
  arrange(desc(t)) 
up_genes_list <- head(deg_sorted$Symbol, 150)
down_genes_list <- tail(deg_sorted$Symbol, 150)

if(!dir.exists("Results/Drug_Analysis")) dir.create("Results/Drug_Analysis", recursive = TRUE)

write.table(up_genes_list, "Results/Drug_Analysis/UP_GENES_FOR_WEB.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down_genes_list, "Results/Drug_Analysis/DOWN_GENES_FOR_WEB.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
# 将网站下载的文件重命名为 'L1000_Result.csv'放入 'Results/Drug_Analysis/' 文件夹


library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyr)
library(grid)
if(!exists("expr_final") | !exists("risk_data") | !exists("pred_ic50")) {
  stop("缺少必要数据！请先加载 expr_final, risk_data, pred_ic50")
}

# 提取公共样本
common_samples <- intersect(intersect(colnames(expr_final), risk_data$sample_id), 
                            substr(rownames(pred_ic50), 1, 12))

cat(">>> 样本对齐完成，公共样本数:", length(common_samples), "\n")

# 对齐数据子集
expr_sub <- expr_final[, common_samples]
risk_sub <- risk_data[match(common_samples, risk_data$sample_id), ]
ic50_sub <- pred_ic50[match(common_samples, substr(rownames(pred_ic50), 1, 12)), ]
# 定义统一配色
my_colors <- c("Low" = "#457B9D", "High" = "#E63946")


#  L1000 预测图
lincs_file <- "Results/Drug_Analysis/L1000_Result.csv"

if(file.exists(lincs_file)){
  raw_lincs <- read.csv(lincs_file)
  lincs_data <- raw_lincs %>%
    dplyr::rename(pert = drug, NCS = similarity_scores) %>% 
    mutate(pert = tolower(pert)) %>% 
    filter(NCS < -0.01) %>% 
    group_by(pert) %>% slice_min(order_by = NCS, n = 1) %>% ungroup() %>%
    arrange(NCS) %>% head(20)
  
  p_a <- ggplot(lincs_data, aes(x = reorder(pert, -NCS), y = NCS)) +
    geom_segment(aes(xend = reorder(pert, -NCS), y = 0, yend = NCS), color = "grey80", size = 0.8) +
    geom_point(aes(size = abs(NCS)), color = "#457B9D", alpha = 1) + 
    scale_size_continuous(range = c(2, 6)) + 
    coord_flip() +
    labs(title = "L1000 Drug Prediction", 
         subtitle = "Broad Spectrum Screening",
         x = NULL, y = "Reverse Score", size="Score") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_line(color="grey95"),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(face="bold", hjust=0.5),
          legend.position = "bottom", legend.box = "horizontal")
} else {
  p_a <- ggplot() + annotate("text", x=1, y=1, label="L1000 Result Missing") + theme_void()
}
ggsave("Results/Drug_Analysis/FigA.pdf", p_a, width = 4, height = 8)
ggsave("Results/Drug_Analysis/FigA.svg", p_a, width = 4, height = 8)



# 耐药药物
target_resist_keywords <- c("Paclitaxel", "Docetaxel", "Cisplatin", "Gemcitabine", 
                            "Palbociclib", "Olaparib", "Doxorubicin", "Irinotecan")
# 在列名中匹配真实名称
found_resist_drugs <- c()
for(k in target_resist_keywords){
  match <- grep(k, colnames(ic50_sub), ignore.case=T, value=T)[1]
  if(!is.na(match)) found_resist_drugs <- c(found_resist_drugs, match)
}

# 补齐8个
if(length(found_resist_drugs) < 8){
  diffs <- colMeans(ic50_sub[risk_sub$risk_group=="High",]) - colMeans(ic50_sub[risk_sub$risk_group=="Low",])
  top_resist <- names(sort(diffs, decreasing = T))
  top_resist <- setdiff(top_resist, found_resist_drugs) 
  found_resist_drugs <- c(found_resist_drugs, head(top_resist, 8 - length(found_resist_drugs)))
}

# 提取数据
plot_dat_b <- ic50_sub[, found_resist_drugs] %>%
  as.data.frame() %>%
  mutate(risk_group = risk_sub$risk_group) %>%
  pivot_longer(cols = -risk_group, names_to = "Drug", values_to = "IC50") %>%
  mutate(Drug = sapply(strsplit(Drug, "_"), `[`, 1)) 

# 绘图
p_b <- ggplot(plot_dat_b, aes(x=risk_group, y=IC50, fill=risk_group)) +
  geom_jitter(width=0.3, size=0.3, color="#B0B0B0", alpha=0.8) + 
  geom_boxplot(outlier.shape = NA, alpha=0.9, width=0.6, size=0.3) +
  stat_compare_means(label = "p.format", label.x = 1.5, size=3, vjust=1) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~Drug, scales = "free", nrow = 2, ncol = 4) + 
  labs(title = "Standard Therapy Resistance", 
       subtitle = "High IC50 in High-Risk Group", 
       x=NULL, y="Predicted IC50") +
  
  theme_classic() + 
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5),
    
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9, color = "black"),
    axis.line = element_line(size = 0.4),
    
    strip.background = element_rect(fill=NA, color=NA), 
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    
    legend.position = "none",
    panel.grid = element_blank()
  )

ggsave("Results/Drug_Analysis/FigB.pdf", p_b, width = 8, height = 8)
ggsave("Results/Drug_Analysis/FigB.svg", p_b, width = 8, height = 8)


# 分子机制 (19个基因相关性棒棒糖图) - 风格融合版
target_genes_c <- c("GATA3", "XBP1", "AR", "FOXA1", "PGR","ESR1", 
                    "BRCA1", "E2F2", "FOSL1", "E2F8", "POU5F1", 
                    "E2F4", "RUNX1", "SNAI2", "CCND1", "TWIST2", 
                    "EGR3", "PBX1", "SMAD3")

# 确定参考药物
if(length(found_resist_drugs) >= 5){
  ref_drug <- found_resist_drugs[5]
} else {
  ref_drug <- found_resist_drugs[1]
}
ref_drug_clean <- strsplit(ref_drug, "_")[[1]][1]

# 2. 计算相关性
cor_list <- list()
for(g in target_genes_c){
  if(g %in% rownames(expr_sub)){
    # 计算 Spearman 相关性
    res <- cor.test(as.numeric(expr_sub[g,]), as.numeric(ic50_sub[, ref_drug]), method="spearman")
    cor_list[[g]] <- data.frame(Gene=g, R=res$estimate, P=res$p.value)
  }
}
# 按相关性 R 从大到小排序
cor_df <- do.call(rbind, cor_list) %>% arrange(R) 
cor_df$Gene <- factor(cor_df$Gene, levels = cor_df$Gene)

p_c <- ggplot(cor_df, aes(x=Gene, y=R)) +
  geom_segment(aes(xend=Gene, yend=0), color="grey80", size=0.8) +
  geom_point(aes(size=abs(R), color=R), alpha=1) +
  scale_color_gradient2(low="#1F77B4", mid="white", high="#D62728", midpoint=0, name="Cor (R)") +
  scale_size_continuous(range = c(3, 8), name = "|R|") +
  scale_y_continuous(limits = c(-0.5, 0.5)) + 
  coord_flip() +
  
  labs(title = "Mechanism of Resistance",
       subtitle = paste0("Correlation with ", ref_drug_clean, " IC50"),
       x=NULL, y="Spearman Correlation") +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face="bold", hjust=0.5, size=14),
    plot.subtitle = element_text(hjust=0.5, color="grey40", size=10),
   
    axis.text.y = element_text(size = 10, color="black"),
    axis.text.x = element_text(size = 10, color="black"),
    
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color="grey95", linetype="dashed"),
   
    panel.border = element_rect(colour = "black", size = 1.2),
    
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )

ggsave("Results/Drug_Analysis/FigC.pdf", p_c, width = 4.5, height = 8)
ggsave("Results/Drug_Analysis/FigC.svg", p_c, width = 4.5, height = 8)



# 免疫小提琴图
my_colors_deep <- c("Low" = "#457B9D", "High" = "#E63946")
my_colors_light <- c("Low" = "#A2D2FF", "High" = "#FFA8A8")

imm_genes <- c("HAVCR2", "CD274", "CD8A", "PDCD1") 

df_imm <- expr_sub[imm_genes, ] %>% t() %>% as.data.frame() %>%
  mutate(risk_group = risk_sub$risk_group) %>%
  pivot_longer(cols = -risk_group, names_to = "ID", values_to = "Value") %>%
  mutate(Type = "Immune")

# 设置因子水平星
df_imm$ID <- factor(df_imm$ID, levels = imm_genes)

# 准备数据：敏感药物
diffs_sens <- colMeans(ic50_sub[risk_sub$risk_group=="High",]) - colMeans(ic50_sub[risk_sub$risk_group=="Low",])
candidate_drugs <- names(sort(diffs_sens))

valid_sens_drugs <- c()
for(drug in candidate_drugs){
  tmp_dat <- data.frame(val = ic50_sub[, drug], group = risk_sub$risk_group)
  res_test <- wilcox.test(val ~ group, data = tmp_dat)
  if(res_test$p.value < 0.05){
    valid_sens_drugs <- c(valid_sens_drugs, drug)
  }
  if(length(valid_sens_drugs) >= 4) break
}

cat("筛选出的敏感药物 (P<0.05):", paste(valid_sens_drugs, collapse=", "), "\n")

df_sens <- ic50_sub[, valid_sens_drugs] %>% as.data.frame() %>%
  mutate(risk_group = risk_sub$risk_group) %>%
  pivot_longer(cols = -risk_group, names_to = "ID", values_to = "Value") %>%
  mutate(Type = "Drug") %>%
  mutate(ID = sapply(strsplit(ID, "_"), `[`, 1)) 

p_imm_violin <- ggplot(df_imm, aes(x=risk_group, y=Value, fill=risk_group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) + 
  geom_boxplot(width = 0.15,size = 0.4, fill = "white", color = "black", outlier.shape = NA, alpha = 1) +
  stat_compare_means(label = "p.format", label.x = 1.5, size=3, vjust=1) +
  scale_fill_manual(values = my_colors_light) +
  facet_wrap(~ID, scales = "free_y", nrow = 1, ncol = 4) + 
  
  labs(title = "Immune Checkpoints (Violin Plot)", 
       subtitle = "Expression Level (TPM/FPKM)", 
       x=NULL, y="Expression") +
  
  theme_classic() + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.line = element_line(size = 0.4),
  )
print(p_imm_violin)
ggsave("Results/Drug_Analysis/FigViolin.svg", p_imm_violin, width = 8, height = 4)

# 敏感药物箱线图
p_sens_box <- ggplot(df_sens, aes(x=risk_group, y=Value, fill=risk_group)) +
  geom_jitter(width=0.25, size=0.3, color="#B0B0B0", alpha=0.6) +
  geom_boxplot(outlier.shape = NA, alpha=0.8, width=0.6, size=0.3) +
  stat_compare_means(label = "p.format", label.x = 1.5, size=3, vjust=1) +
  scale_fill_manual(values = my_colors_deep) +
  facet_wrap(~ID, scales = "free_y", nrow = 1, ncol = 4) + 
  labs(title = "Sensitive Targeted Drugs (Boxplot)", 
       subtitle = "Predicted IC50 (Lower = More Sensitive)", 
       x=NULL, y="Predicted IC50") +
  
  theme_classic() + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.line = element_line(size = 0.4)
  )


p_d_final <- p_imm_violin / p_sens_box + 
  plot_annotation(tag_levels = 'A') 
ggsave("Results/Drug_Analysis/FigD.pdf", p_d_final, width = 8, height = 8)
ggsave("Results/Drug_Analysis/FigD.svg", p_d_final, width = 8, height = 8)


# 生成汇总表格
if(!dir.exists("Results/Drug_Analysis")) dir.create("Results/Drug_Analysis", recursive = TRUE)

# 提取 L1000FWD Top 8 逆转药物
lincs_file <- "Results/Drug_Analysis/L1000_Result.csv"

if(file.exists(lincs_file)){
  raw_lincs <- read.csv(lincs_file)
  
  top_l1000 <- raw_lincs %>%
    dplyr::rename(Drug = drug, Score = similarity_scores) %>%
    mutate(Drug = tolower(Drug)) %>%
    filter(Score < 0) %>% # 只保留负分
    group_by(Drug) %>% 
    slice_min(order_by = Score, n = 1) %>% 
    ungroup() %>%
    arrange(Score) %>% 
    head(8) %>% 
    mutate(
      Candidate_Drug = tools::toTitleCase(Drug), 
      Screening_Approach = "L1000FWD (Signature Reversal)",
      Key_Metric = paste0("Score: ", round(Score, 3)),
      Mechanism_of_Action_Target = "" 
    ) %>%
    dplyr::select(Candidate_Drug, Screening_Approach, Key_Metric, Mechanism_of_Action_Target)
} else {
  stop("❌ 未找到 L1000_Result.csv 文件！")
}


# 提取 GDSC Top 8 高风险敏感药物
common_samples <- intersect(risk_data$sample_id, substr(rownames(pred_ic50), 1, 12))
risk_sub <- risk_data[match(common_samples, risk_data$sample_id), ]
ic50_sub <- pred_ic50[match(common_samples, substr(rownames(pred_ic50), 1, 12)), ]

high_idx <- risk_sub$risk_group == "High"
low_idx <- risk_sub$risk_group == "Low"

# 循环计算所有药物的差异和P值
res_list <- list()
for(d in colnames(ic50_sub)) {
  val_high <- as.numeric(ic50_sub[high_idx, d])
  val_low <- as.numeric(ic50_sub[low_idx, d])
  diff_mean <- mean(val_high, na.rm=TRUE) - mean(val_low, na.rm=TRUE)
  
  if(diff_mean < 0) { 
    p_val <- wilcox.test(val_high, val_low)$p.value
    res_list[[d]] <- data.frame(Drug = d, Diff = diff_mean, Pval = p_val)
  }
}

gdsc_df <- do.call(rbind, res_list)
# 筛选最显著的前 8 个
top_gdsc <- gdsc_df %>%
  arrange(Pval, Diff) %>% 
  head(8) %>%
  mutate(
    Candidate_Drug = sapply(strsplit(Drug, "_"), `[`, 1),
    Screening_Approach = "GDSC (IC50 Prediction)",
    Key_Metric = ifelse(Pval < 0.0001, "p < 0.0001", paste0("p = ", signif(Pval, 2))),
    Mechanism_of_Action_Target = "" # 留空待手动填补
  ) %>%
  dplyr::select(Candidate_Drug, Screening_Approach, Key_Metric, Mechanism_of_Action_Target)

# 合并并导出表格
final_table <- rbind(top_l1000, top_gdsc)

output_path <- "Results/Drug_Analysis/Table1_Potential_Therapeutics.csv"
write.csv(final_table, output_path, row.names = FALSE)
