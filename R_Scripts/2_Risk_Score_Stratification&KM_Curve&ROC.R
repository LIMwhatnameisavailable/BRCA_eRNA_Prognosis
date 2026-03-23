coef.min <- coef(cvfit, s = "lambda.min")
selected_eRNAs <- rownames(coef.min)[which(coef.min != 0)]
print(selected_eRNAs)  #10selected_eRNA

#训练集
colnames(train_data) <- gsub(":", "_", colnames(train_data))
colnames(train_data) <- gsub("/month", "_month", colnames(train_data))
train_data$risk_score <- (-0.12338077) * train_data[,"chr1_155158995"]+
  (0.21888531) * train_data[,"chr3_11236700"]+( -0.24173782) * train_data[,"chr8_22624675"]+
  (-0.85296049) * train_data[,"chr3_138070534"]+( -0.20550119) * train_data[,"chr9_71398719"]+
  (-0.20987822) * train_data[,"chr9_114689796"]+( 0.21679753) * train_data[,"chr10_5531356"]+
  (-0.03772779) * train_data[,"chr9_71398939"]+( 0.35946524) * train_data[,"chr12_13371038"]+
  (0.06413942) * train_data[,"chr10_5528926"]
summary(train_data$OS.time_month)
median_risk <- median(train_data$risk_score,na.rm = TRUE)
print(median_risk)
train_data$risk_group <- ifelse(train_data$risk_score > median_risk, "High", "Low") 

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
class(train_data$OS)
#风险评分分布图 
train_data$OS.time_month <- ifelse(
  train_data$OS.time_month >= 150,
  150,
  train_data$OS.time_month
)
train_data$OS <- ifelse(
  train_data$OS.time_month >= 150,
  0,
  train_data$OS
)
train_data$risk_group <- ifelse(train_data$risk_score > median_risk, "High", "Low")

train_data$status_factor <- train_data$OS
#按风险评分排序样本
train_data <- train_data[order(train_data$risk_score), ] 
train_data$sample_order <- 1:nrow(train_data)  
train_data$id <- rownames(train_data)


median_sample_order <- median(train_data$sample_order)
print(median_sample_order)
p_risk <- ggplot(train_data, aes(x = sample_order, y = risk_score)) +
  geom_point(aes(color = risk_group), size = 2, alpha = 0.8) + 
  geom_hline(yintercept = median_risk, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = median_sample_order, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("High" = "#D62728", "Low" = "#1F77B4")) +
  labs(y = "Risk Score", x = "", color = "Risk Group") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(0.98, 0.02), 
    legend.justification = c(1, 0),  
    legend.box.background = element_rect(color = NA), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 8, unit = "mm")
  )
train_data$status_factor <- factor(
  train_data$status_factor, 
  levels = c(0, 1), 
  labels = c("Alive", "Dead")
)
p_time <- ggplot(train_data, aes(x = sample_order, y = OS.time_month/12)) +  #转换为年
  geom_point(aes(shape = status_factor, color = status_factor), size = 2.5) +
  scale_color_manual(values = c("Alive" = "#1F77B4", "Dead" = "#D62728")) +
  scale_shape_manual(values = c(16, 16)) +
  geom_vline(xintercept = median_sample_order, linetype = "dashed", color = "grey50") +
  labs(y = "Survival Time (Years)",
       x = "",
       color = "OS Status", shape = "OS Status") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.98, 0.98),  
    legend.justification = c(1, 1),   
    legend.box = "horizontal",        
    plot.margin = margin(t = 5, r = 5, b = 5, l = 8, unit = "mm")
  )
library(cowplot)
combined_plot <- plot_grid(
  p_risk, 
  p_time,
  ncol = 1, 
  align = "v", 
  axis = "lr",
  rel_heights = c(0.4, 0.6)
)
#添加统一标题
title <- ggdraw() + 
  draw_label("BRCA Prognostic Risk Stratification(Training Cohort)", 
             fontface = "bold", size = 14)
final_plot <- plot_grid(
  title,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.05, 0.95)
)
svg(filename = "Results/Fig_stratification_training.svg", width = 12, height = 8)
print(final_plot)
dev.off()


library(svglite)
#绘制训练集生存曲线
train_data <- na.omit(train_data)
sur_fit <- survfit(Surv(OS.time_month, OS) ~ risk_group, data = train_data)
table(train_data$risk_group)
p <- ggsurvplot(
  sur_fit,
  title = "BRCA Training Cohort (by eRNA)",  
  font.title = c(14, "bold", "black"),
  data = train_data,
  pval = TRUE,
  pval.method = TRUE,
  pval.method.coord = c(20, 0.25),  
  pval.coord = c(20, 0.18),  
  pval.size = 4.5,   
  pval.method.size = 4.5,  
  xlim = c(0, 150),      
  break.time.by = 30,
  risk.table = TRUE,
  risk.table.title = "Number at Risk",  
  risk.table.y.text = FALSE, 
  risk.table.col = "strata",  
  risk.table.height = 0.18, 
  palette = c("#D62728", "#0072B2"),  
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  legend.title = "Risk Group",
  legend.labs = c("High Risk", "Low Risk"),
  legend = c(0.85, 0.87),     
  font.x = c(12, "bold"),
  font.y = c(12, "bold"),
  font.tickslab = c(11, "plain"),
  ggtheme = theme_classic2() + 
    theme(
      panel.border = element_rect(colour = "black", fill = NA), 
      axis.line = element_line(size = 0.8),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(color = "black"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
)
#调整风险表样式
p$table <- p$table + 
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_rect(
      fill = c("#D62728", "#0072B2")
    ),
    strip.text.y = element_text(
      face = "bold",
      color = "white" 
    )
  )
#合并生存曲线与风险表
svg("Results/Fig_KM_training.svg", width = 12, height = 8)
print(p) 
dev.off()



#测试集检验
colnames(test_data) <- gsub(":", "_", colnames(test_data))
colnames(test_data) <- gsub("/month", "_month", colnames(test_data))

test_data$risk_score <- (-0.12338077) * test_data[,"chr1_155158995"]+
  (0.21888531) * test_data[,"chr3_11236700"]+( -0.24173782) * test_data[,"chr8_22624675"]+
  (-0.85296049) * test_data[,"chr3_138070534"]+( -0.20550119) * test_data[,"chr9_71398719"]+
  (-0.20987822) * test_data[,"chr9_114689796"]+( 0.21679753) * test_data[,"chr10_5531356"]+
  (-0.03772779) * test_data[,"chr9_71398939"]+( 0.35946524) * test_data[,"chr12_13371038"]+
  (0.06413942) * test_data[,"chr10_5528926"]

median_risk_test <- median(test_data$risk_score,na.rm = TRUE) 
print(median_risk_test)
test_data$risk_group <- ifelse(test_data$risk_score > median_risk_test, "High", "Low") 
test_data$status_factor <- test_data$OS
test_data$OS.time_month <- ifelse(
  test_data$OS.time_month >= 120,
  120,
  test_data$OS.time_month
)
test_data$OS <- ifelse(
  test_data$OS.time_month >= 120,
  0,
  test_data$OS
)
summary(test_data$OS.time_month)


#测试集
test_data <- test_data[order(test_data$risk_score), ] 
test_data$sample_order <- 1:nrow(test_data) 
median_sample_order <- median(test_data$sample_order)
print(median_sample_order)
test_data$status_factor <- factor(
  test_data$status_factor, 
  levels = c(0, 1), 
  labels = c("Alive", "Dead")
)
p_risk_test <- ggplot(test_data, aes(x = sample_order, y = risk_score)) +
  geom_point(aes(color = risk_group), size = 2, alpha = 0.8) + 
  geom_hline(yintercept = median_risk, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = median_sample_order, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("High" = "#D62728", "Low" = "#1F77B4")) +
  labs(y = "Risk Score", x = "", color = "Risk Group") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(0.98, 0.02), 
    legend.justification = c(1, 0),  
    legend.box.background = element_rect(color = NA), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 8, unit = "mm")
  )
p_time_test <- ggplot(test_data, aes(x = sample_order, y = OS.time_month/12)) + 
  geom_point(aes(shape = status_factor, color = status_factor), size = 2.5) +
  scale_color_manual(values = c("Alive" = "#1F77B4", "Dead" = "#D62728")) +
  scale_shape_manual(values = c(16, 16)) +
  geom_vline(xintercept = median_sample_order, linetype = "dashed", color = "grey50") +
  labs(y = "Survival Time (Years)",
       x = "",
       color = "OS Status", shape = "OS Status") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.98, 0.98),  
    legend.justification = c(1, 1),   
    legend.box = "horizontal",        
    plot.margin = margin(t = 5, r = 5, b = 5, l = 8, unit = "mm")
  )
combined_plot_test <- plot_grid(
  p_risk_test, 
  p_time_test,
  ncol = 1, 
  align = "v", 
  axis = "lr",
  rel_heights = c(0.4, 0.6)
)
#添加统一标题
title <- ggdraw() + 
  draw_label("BRCA Prognostic Risk Stratification(Testing Cohort)", 
             fontface = "bold", size = 14)
final_plot_test <- plot_grid(
  title,
  combined_plot_test,
  ncol = 1,
  rel_heights = c(0.05, 0.95)
)
svg(filename = "Results/Fig_stratification_testing.svg", width = 12, height = 8)
print(final_plot_test)
dev.off()


#绘制测试集生存曲线
test_data <- na.omit(test_data)
test_data$OS.time_month <- ifelse(
  test_data$OS.time_month >= 120,
  120,
  test_data$OS.time_month
)
test_data$OS <- ifelse(
  test_data$OS.time_month >= 120,
  0,
  test_data$OS
)
sur_fit_test <- survfit(Surv(OS.time_month, OS) ~ risk_group, data = test_data)
table(test_data$risk_group)
p <- ggsurvplot(
  sur_fit_test,
  title = "BRCA Testing Cohort (by eRNA)",  
  font.title = c(14, "bold", "black"),
  data = test_data,
  pval = TRUE,
  pval.method = TRUE,
  pval.method.coord = c(20, 0.25),  
  pval.coord = c(20, 0.18),  
  pval.size = 4.5,   
  pval.method.size = 4.5,  
  risk.table = TRUE,
  risk.table.title = "Number at Risk",  
  risk.table.y.text = FALSE, 
  risk.table.col = "strata",  
  risk.table.height = 0.18, 
  palette = c("#D62728", "#0072B2"),  
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  xlim = c(0,120),
  break.time.by = 30,       
  legend.title = "Risk Group",
  legend.labs = c("High Risk", "Low Risk"),
  legend = c(0.85, 0.87),     
  font.x = c(12, "bold"),
  font.y = c(12, "bold"),
  font.tickslab = c(11, "plain"),
  ggtheme = theme_classic2() + 
    theme(
      panel.border = element_rect(colour = "black", fill = NA), 
      axis.line = element_line(size = 0.8),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(color = "black"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
)
#调整风险表样式
p$table <- p$table + 
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_rect(
      fill = c("#D62728", "#0072B2")
    ),
    strip.text.y = element_text(
      face = "bold",
      color = "white" 
    )
  )
svg("Results/Fig_KM_testing.svg", width = 12, height = 8)
print(p)
dev.off()


#总体检验
colnames(combined_data) <- gsub(":", "_", colnames(combined_data))
colnames(combined_data) <- gsub("/month", "_month", colnames(combined_data))

combined_data$risk_score <- (-0.12338077) * combined_data[,"chr1_155158995"]+
  (0.21888531) * combined_data[,"chr3_11236700"]+( -0.24173782) * combined_data[,"chr8_22624675"]+
  (-0.85296049) * combined_data[,"chr3_138070534"]+( -0.20550119) * combined_data[,"chr9_71398719"]+
  (-0.20987822) * combined_data[,"chr9_114689796"]+( 0.21679753) * combined_data[,"chr10_5531356"]+
  (-0.03772779) * combined_data[,"chr9_71398939"]+( 0.35946524) * combined_data[,"chr12_13371038"]+
  (0.06413942) * combined_data[,"chr10_5528926"]

median_risk_combined <- median(combined_data$risk_score,na.rm = TRUE) 
print(median_risk_combined)
combined_data$risk_group <- ifelse(combined_data$risk_score > median_risk_combined, "High", "Low") 
combined_data$status_factor <- combined_data$OS
combined_data$OS.time_month <- ifelse(
  combined_data$OS.time_month >= 150,
  150,
  combined_data$OS.time_month
)
combined_data$OS <- ifelse(
  combined_data$OS.time_month >= 150,
  0,
  combined_data$OS
)
summary(combined_data$OS.time_month)

#总集
combined_data <- combined_data[order(combined_data$risk_score), ] 
combined_data$sample_order <- 1:nrow(combined_data) 
median_sample_order <- median(combined_data$sample_order)
print(median_sample_order)
combined_data$status_factor <- factor(
  combined_data$status_factor, 
  levels = c(0, 1), 
  labels = c("Alive", "Dead")
)
p_risk_combined <- ggplot(combined_data, aes(x = sample_order, y = risk_score)) +
  geom_point(aes(color = risk_group), size = 2, alpha = 0.8) + 
  geom_hline(yintercept = median_risk, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = median_sample_order, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("High" = "#D62728", "Low" = "#1F77B4")) +
  labs(y = "Risk Score", x = "", color = "Risk Group") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(0.98, 0.02), 
    legend.justification = c(1, 0),  
    legend.box.background = element_rect(color = NA), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 8, unit = "mm")
  )
p_time_combined <- ggplot(combined_data, aes(x = sample_order, y = OS.time_month/12)) + 
  geom_point(aes(shape = status_factor, color = status_factor), size = 2.5) +
  scale_color_manual(values = c("Alive" = "#1F77B4", "Dead" = "#D62728")) +
  scale_shape_manual(values = c(16, 16)) +
  geom_vline(xintercept = median_sample_order, linetype = "dashed", color = "grey50") +
  labs(y = "Survival Time (Years)",
       x = "",
       color = "OS Status", shape = "OS Status") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.98, 0.98),  
    legend.justification = c(1, 1),   
    legend.box = "horizontal",        
    plot.margin = margin(t = 5, r = 5, b = 5, l = 8, unit = "mm")
  )
library(cowplot)
combined_plot_combined <- plot_grid(
  p_risk_combined, 
  p_time_combined,
  ncol = 1, 
  align = "v", 
  axis = "lr",
  rel_heights = c(0.4, 0.6)
)
#添加统一标题
title <- ggdraw() + 
  draw_label("BRCA Prognostic Risk Stratification(Combined Cohort)", 
             fontface = "bold", size = 14)
final_plot_combined <- plot_grid(
  title,
  combined_plot_combined,
  ncol = 1,
  rel_heights = c(0.05, 0.95)
)
svg(filename = "Results/Fig_stratification_combined.svg", width = 12, height = 8)
print(final_plot_combined)
dev.off()

#绘制总体生存曲线
combined_data <- na.omit(combined_data)
sur_fit_combined <- survfit(Surv(OS.time_month, OS) ~ risk_group, data = combined_data)
table(combined_data$risk_group)
p <- ggsurvplot(
  sur_fit_combined,
  title = "BRCA Combined Cohort (by eRNA)",  
  font.title = c(14, "bold", "black"),
  data = combined_data,
  pval = TRUE,
  pval.method = TRUE,
  pval.method.coord = c(20, 0.25),  
  pval.coord = c(20, 0.18),  
  pval.size = 4.5,   
  pval.method.size = 4.5,  
  risk.table = TRUE,
  risk.table.title = "Number at Risk",  
  risk.table.y.text = FALSE, 
  risk.table.col = "strata",  
  risk.table.height = 0.18, 
  palette = c("#D62728", "#0072B2"),  
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  break.time.by = 30,      
  xlim = c(0,150),
  legend.title = "Risk Group",
  legend.labs = c("High Risk", "Low Risk"),
  legend = c(0.85, 0.87),     
  font.x = c(12, "bold"),
  font.y = c(12, "bold"),
  font.tickslab = c(11, "plain"),
  ggtheme = theme_classic2() + 
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      axis.line = element_line(size = 0.8),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(color = "black"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
)

#调整风险表样式
p$table <- p$table + 
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_rect(
      fill = c("#D62728", "#0072B2") ),
    strip.text.y = element_text(
      face = "bold",
      color = "white" 
    )
  )

svg("Results/Fig_KM_combined.svg", width = 12, height = 8)
print(p) 
dev.off()

# 总集热图
selected_genes <- c("chr1_155158995", "chr3_11236700", "chr8_22624675","chr3_138070534","chr9_71398719",
                    "chr9_114689796","chr10_5531356","chr9_71398939","chr12_13371038","chr10_5528926") 
expr_matrix_combined <- t(scale((combined_data[, selected_genes]), center = TRUE, scale = TRUE))
rownames(expr_matrix_combined) <- selected_genes
custom_order <- c(
  "chr10_5531356",
  "chr10_5528926",
  "chr3_11236700",
  "chr12_13371038",
  "chr3_138070534",
  "chr8_22624675",
  "chr1_155158995",
  "chr9_114689796",
  "chr9_71398719",
  "chr9_71398939"
)
gene_index_combined <- match(custom_order, rownames(expr_matrix_combined))
risk_colors <- c("Low" = "#1F77B4", "High" = "#D62728")
top_anno <- HeatmapAnnotation(
  `Risk Group` = combined_data$risk_group,
  col = list(`Risk Group` = risk_colors),
  show_legend = FALSE,  
  annotation_name_side = "right",  
  simple_anno_size = unit(0.5, "cm") 
)
gold_standard_sample_order <- rownames(combined_data)
expr_matrix_synced <- expr_matrix_combined[, gold_standard_sample_order]

p_heatmap_combined <- Heatmap(
  expr_matrix_synced,
  name = "Expression Z-score",
  col = circlize::colorRamp2(c(-3, 0, 3), c("#1F77B4", "white", "#D62728")),
  
  cluster_rows = FALSE,  
  row_order = gene_index_combined, 
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_labels = gsub("_", ":", rownames(expr_matrix_synced)),  
  row_names_gp = gpar(fontsize = 10),
  top_annotation = top_anno,  
  heatmap_legend_param = list(
    title = NULL,  
    legend_direction = "vertical",  
    at = c(-3, 0, 3),  
    labels = c("Low", "Mid", "High"),
    legend_height = unit(3, "cm"),  
    legend_width = unit(0.8, "cm"),  
    grid_border = "black",  
    labels_gp = gpar(fontsize = 10)
  )
)
svg(filename = "Results/Fig_10eRNA_heatmap_combined.svg", width = 12, height = 8)
print(p_heatmap_combined)
dev.off()


# 总集ROC
library(timeROC)
library(ggplot2)
roc_res <- timeROC(
  T = combined_data$OS.time,
  delta = combined_data$OS,
  marker = combined_data$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(12, 36, 60), 
  iid = TRUE
)
dat_roc <- rbind(
  data.frame(FPR = roc_res$FP[,1], TPR = roc_res$TP[,1], Time = "1-Year"),
  data.frame(FPR = roc_res$FP[,2], TPR = roc_res$TP[,2], Time = "3-Year"),
  data.frame(FPR = roc_res$FP[,3], TPR = roc_res$TP[,3], Time = "5-Year")
)
lab_1 <- paste0("1-Year (AUC = ", round(roc_res$AUC[1], 3), ")")
lab_3 <- paste0("3-Year (AUC = ", round(roc_res$AUC[2], 3), ")")
lab_5 <- paste0("5-Year (AUC = ", round(roc_res$AUC[3], 3), ")")

p_roc <- ggplot(dat_roc, aes(x = FPR, y = TPR, color = Time)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", size = 0.8) +
  geom_path(size = 1.2) + 
  scale_color_manual(
    values = c("#E64B35", "#4DBBD5", "#00A087"), 
    labels = c(lab_1, lab_3, lab_5)
  ) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  labs(title = "ROC Curve (Total Cohort)",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = c(0.65, 0.25), 
    legend.background = element_blank(), 
    legend.key = element_blank(),        
    legend.title = element_blank(),     
    legend.text = element_text(size = 11, face = "bold")
  )
svg(filename = "Results/Fig_ROC_combined.svg", width = 12, height = 8)
print(p_roc)
dev.off()
