if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", 
                       "dplyr", "data.table", "ggplot2", "limma", "pheatmap")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(data.table)
library(ggplot2)
library(limma)
library(pheatmap)


#===============================================================================
#                       第一部分：基因表达矩阵预处理                             
#===============================================================================

expr_path <- "Data_Source/TCGA-BRCA.star_fpkm.tsv.gz"
annot_path <- "Data_Source/gencode.v36.annotation.gtf.gene.probemap"
expr_data <- fread(expr_path, data.table = FALSE)
gene_annot <- fread(annot_path, data.table = FALSE)

# 预处理表达矩阵
original_colnames <- colnames(expr_data)
sample_type <- substr(original_colnames, 14, 15)
tumor_cols_idx <- which(sample_type == "01")
stopifnot(length(tumor_cols_idx) > 0)
expr_data_tumor_raw <- expr_data[, c(1, tumor_cols_idx), drop = FALSE]
cat("筛选出", ncol(expr_data_tumor_raw) - 1, "个原发肿瘤样本。\n")

rownames(expr_data_tumor_raw) <- expr_data_tumor_raw[, 1]
expr_data_tumor_raw <- expr_data_tumor_raw[, -1, drop = FALSE]
gene_ids <- rownames(expr_data_tumor_raw)
# 将所有列转换为数值型
expr_data_tumor_raw <- as.data.frame(lapply(expr_data_tumor_raw, as.numeric))
rownames(expr_data_tumor_raw) <- gene_ids

patient_ids <- substr(colnames(expr_data_tumor_raw), 1, 12)
if (any(duplicated(patient_ids))) {
  cat("检测到来自同一病人的重复样本，现在进行合并（取均值）...\n")
  patient_data_list <- split.default(expr_data_tumor_raw, patient_ids)
  merged_data_list <- lapply(patient_data_list, function(patient_df) {
    if (ncol(patient_df) > 1) {
      rowMeans(patient_df, na.rm = TRUE)
    } else {
      patient_df[[1]]
    }
  })
  expr_data_tumor_unique <- do.call(cbind, merged_data_list)
  
  cat("合并完成。\n")
} else {
  cat("未检测到重复样本，直接使用筛选后的数据。\n")
  expr_data_tumor_unique <- as.matrix(expr_data_tumor_raw)
  colnames(expr_data_tumor_unique) <- patient_ids
}
rownames(expr_data_tumor_unique) <- rownames(expr_data_tumor_raw)
expr_data_tumor <- as.data.frame(expr_data_tumor_unique) %>%
  tibble::rownames_to_column("ensembl_id")
stopifnot(!any(duplicated(colnames(expr_data_tumor)[-1])))
cat("最终得到", ncol(expr_data_tumor) - 1, "个唯一病人的肿瘤表达谱。\n")

# 预处理基因注释文件
gene_annot_clean <- gene_annot %>%
  dplyr::select(id, gene) %>%                
  dplyr::rename(ensembl_id = id, gene_name = gene) %>% 
  dplyr::filter(!is.na(gene_name) & gene_name != "")  

# 4合并与去重
cat("表达数据中的前5个ensembl_id:\n")
print(head(expr_data_tumor$ensembl_id, 5))
cat("\n注释文件中的前5个ensembl_id:\n")
print(head(gene_annot_clean$ensembl_id, 5))

# 检查ID重叠数量
overlap_count <- length(intersect(expr_data_tumor$ensembl_id, gene_annot_clean$ensembl_id))
cat("\n在处理前，两个文件共有的ID数量：", overlap_count, "\n")
overlap_count_after <- length(intersect(expr_data_tumor$ensembl_id, gene_annot_clean$ensembl_id))
cat("在处理后，两个文件共有的ID数量：", overlap_count_after, "\n")
stopifnot(overlap_count_after > 0) 
expr_final <- expr_data_tumor %>%
  # 合并注释
  dplyr::inner_join(gene_annot_clean, by = "ensembl_id") %>%
  dplyr::select(-ensembl_id) %>%
  dplyr::mutate(across(-gene_name, as.numeric)) %>%
  dplyr::mutate(row_mean = rowMeans(dplyr::select(., -gene_name), na.rm = TRUE)) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice_max(order_by = row_mean, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(-row_mean) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("gene_name")

cat("\n最终表达矩阵维度：", nrow(expr_final), "个基因 ×", ncol(expr_final), "个样本\n")


#===============================================================================
#                  第二部分：风险模型计算与样本分组                              
#===============================================================================

risk_scores <- 
  (-0.12338077 * combined_data[,"chr1_155158995"]) +
  (0.21888531 * combined_data[,"chr3_11236700"]) +
  (-0.24173782 * combined_data[,"chr8_22624675"]) +
  (-0.85296049 * combined_data[,"chr3_138070534"]) +
  (-0.20550119 * combined_data[,"chr9_71398719"]) +
  (-0.20987822 * combined_data[,"chr9_114689796"]) +
  (0.21679753 * combined_data[,"chr10_5531356"]) +
  (-0.03772779 * combined_data[,"chr9_71398939"]) +
  (0.35946524 * combined_data[,"chr12_13371038"]) +
  (0.06413942 * combined_data[,"chr10_5528926"])

risk_data <- data.frame(
  sample_id = rownames(combined_data),
  risk_score = risk_scores,
  stringsAsFactors = FALSE
)
median_score <- median(risk_data$risk_score, na.rm = TRUE)
risk_data$risk_group <- ifelse(risk_data$risk_score >= median_score, "High", "Low")
cat("风险评分中位数：", median_score, "\n")
cat("高风险组样本数：", sum(risk_data$risk_group == "High"), "\n")
cat("低风险组样本数：", sum(risk_data$risk_group == "Low"), "\n")


#===============================================================================
#                   第三部分：基因排序 (使用 limma)                         
#===============================================================================

# 匹配表达矩阵和风险分组的样本
colnames(expr_final) <- gsub("\\.", "-", colnames(expr_final))
common_samples <- intersect(colnames(expr_final), risk_data$sample_id)
expr_filtered <- expr_final[, common_samples]
risk_filtered <- risk_data[match(common_samples, risk_data$sample_id), ]

# 确保两者的样本顺序完全一致
risk_filtered <- risk_filtered[order(risk_filtered$sample_id),]
expr_filtered <- expr_filtered[, order(colnames(expr_filtered))]
stopifnot(all(colnames(expr_filtered) == risk_filtered$sample_id))

# 使用limma进行差异表达分析
risk_filtered$risk_group <- factor(risk_filtered$risk_group, levels = c("Low", "High"))
design <- model.matrix(~0 + risk_group, data = risk_filtered)
colnames(design) <- c("Low", "High")
contrast.matrix <- makeContrasts(High_vs_Low = High - Low, levels = design)

# 线性模型拟合与贝叶斯平滑
fit <- lmFit(expr_filtered, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取结果并构建GSEA排序列表
all_genes_diff <- topTable(fit2, coef = "High_vs_Low", number = Inf)
# 使用t统计量作为排序指标
gene_ranks <- all_genes_diff$t
names(gene_ranks) <- rownames(all_genes_diff)
gene_ranks <- sort(gene_ranks, decreasing = TRUE) # t>0 在High组高表达


#===============================================================================
#                           第四部分：执行GSEA分析             
#===============================================================================

# 将基因名(Symbol)转换为ENTREZ ID
ids <- bitr(names(gene_ranks), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_ranks_entrez <- gene_ranks[names(gene_ranks) %in% ids$SYMBOL]
names(gene_ranks_entrez) <- ids$ENTREZID[match(names(gene_ranks_entrez), ids$SYMBOL)]

# 执行GSEA分析
options(timeout = 600)
gsea_kegg_result <- gseKEGG(
  geneList = gene_ranks_entrez,
  organism = "hsa", 
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1, 
  verbose = FALSE
)

# 检查是否有显著结果
if (is.null(gsea_kegg_result) || nrow(gsea_kegg_result@result) == 0) {
  stop("GSEA未返回任何结果，请检查基因列表或参数设置。")
}
cat("GSEA分析完成！\n")


#===============================================================================
#                               第五部分：结果可视化              
#===============================================================================

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
library(svglite)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")
library(ComplexHeatmap)
library(circlize)

sig_gsea <- gsea_kegg_result@result %>%
  filter(p.adjust < 0.25) %>%
  arrange(desc(NES))
cat("在p.adjust < 0.25的阈值下，发现", nrow(sig_gsea), "条显著富集的通路。\n")


library(enrichplot)
if (nrow(sig_gsea) > 0) {
  # 气泡图
  gsea_df <- gsea_kegg_result@result
  activated_df <- gsea_df %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    filter(Description != "Olfactory transduction") %>%
    head(20)
  suppressed_df <- gsea_df %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    filter(Description != "Coronavirus disease - COVID-19") %>%
    head(20)
  
  p_dot_activated <- ggplot(activated_df,
                            aes(x = NES, y = reorder(Description, NES), size = setSize, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#D62728", high = "#1F77B4", name = "p.adjust") +
    scale_size(range = c(3, 9), name = "Count") +
    guides(color = guide_colorbar(order = 1),
           size = guide_legend(order = 2)) +
    labs(x = "NES", y = NULL, title = "Activated") +
    theme_bw() +
    scale_x_continuous(limits = c(0, max(activated_df$NES) * 1.1), expand = c(0, 0.05)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(colour = "black", linewidth = 1.2),
      legend.position = c(0.02, 0.99),
      legend.justification = c("left", "top"),
      legend.box.just = "left",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8)
    )

  p_dot_suppressed <- ggplot(suppressed_df,
                             aes(x = NES, y = reorder(Description, NES), size = setSize, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#D62728", high = "#1F77B4", name = "p.adjust") +
    scale_size(range = c(3, 9), name = "Count") +
    guides(color = guide_colorbar(order = 1),
           size = guide_legend(order = 2)) +
    labs(x = "NES", y = NULL, title = "Suppressed") +
    theme_bw() +
    scale_x_continuous(limits = c(min(suppressed_df$NES) * 1.1, 0), expand = c(0, 0.05)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(colour = "black", linewidth = 1.2),
      legend.position = c(0.98, 0.99),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8)
    )
  # 拼接并保存
  p_dotplot_combined <- p_dot_activated + p_dot_suppressed +
    plot_annotation(title = "GSEA Results for High vs. Low Risk Groups",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)))
  print(p_dotplot_combined)
  ggsave("Results/Fig_GSEA_Dotplot.svg", p_dotplot_combined,
         device = svglite::svglite,
         width = 18, height = 9)
  
  
  # 多通路富集图
  high_risk_ids <- filter(sig_gsea, NES > 0) %>%
    arrange(desc(NES)) %>%
    slice_max(order_by = NES, n = 6) %>%
    filter(Description != "Olfactory transduction") %>%
    pull(ID)
  low_risk_ids_to_plot <- sig_gsea %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    slice_min(order_by = NES, n = 6) %>%
    filter(Description != "Coronavirus disease - COVID-19") %>%
    pull(ID)
 
  # 高风险组
  if (length(high_risk_ids) > 0) {
    base_plot_high <- gseaplot2(gsea_kegg_result, geneSetID = high_risk_ids)
    plot_top_high <- base_plot_high[[1]]
    plot_middle_high <- base_plot_high[[2]]
    plot_bottom_high <- base_plot_high[[3]]
    plot_top_high <- plot_top_high +
      theme(
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.key.height = unit(0.8, 'cm'),
        legend.position = c(0.99, 0.9),
        legend.justification = c("right", "top")
      )
    plot_middle_high$layers[[1]]$aes_params$linewidth <- 1.2
    plot_bottom_high <- plot_bottom_high +
      theme(
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12)
      )
    plot_top_high <- ggrastr::rasterise(plot_top_high, layers = "GeomLine", dpi = 300)
    plot_middle_high <- ggrastr::rasterise(plot_middle_high, layers = "GeomSegment", dpi = 300)
    reassembled_plot_high <- plot_top_high / plot_middle_high / plot_bottom_high +
      plot_layout(heights = c(0.55, 0.2, 0.25))
    final_plot_high <- reassembled_plot_high +
      plot_annotation(
        title = "Top Enriched Pathways in High-Risk Group",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
      )
    print(final_plot_high)
    ggsave("Results/Fig_GSEA_Multiplot_HighRisk.svg", final_plot_high,
           device = svglite::svglite,
           width = 12, height = 8)
  }
# 低风险组
if (length(low_risk_ids_to_plot) > 0) {
  base_plot_low <- gseaplot2(gsea_kegg_result, geneSetID = low_risk_ids_to_plot)
  plot_top_low <- base_plot_low[[1]]
  plot_middle_low <- base_plot_low[[2]]
  plot_bottom_low <- base_plot_low[[3]]
  plot_top_low <- plot_top_low +
    theme(
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.key.height = unit(0.8, 'cm'),
      legend.position = c(0.9, 0.9),
      legend.justification = c("right", "top")
    )
  plot_middle_low$layers[[1]]$aes_params$linewidth <- 1.2
  plot_bottom_low <- plot_bottom_low +
    theme(
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12)
    )
  plot_middle_low <- ggrastr::rasterise(plot_middle_low, dpi = 300)
  reassembled_plot_low <- plot_top_low / plot_middle_low / plot_bottom_low +
    plot_layout(heights = c(0.55, 0.2, 0.25))
  final_plot_low <- reassembled_plot_low +
    plot_annotation(
      title = "Top Enriched Pathways in Low-Risk Group",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
    )
  print(final_plot_low)
  ggsave("Results/Fig_GSEA_Multiplot_LowRisk.svg", final_plot_low,
         device = svglite::svglite,
         width = 12, height = 8)
  }
}

# 核心基因热图
# 高风险组
high_risk_top_pathway <- sig_gsea %>% filter(NES > 0) %>% slice_max(NES, n = 1)
if (nrow(high_risk_top_pathway) > 0) {
  pathway_id <- high_risk_top_pathway$ID
  pathway_name <- high_risk_top_pathway$Description
  cat("\n正在为高风险组Top通路生成热图:", pathway_name, "\n")
  core_genes_entrez <- sig_gsea$core_enrichment[sig_gsea$ID == pathway_id] %>% strsplit("/") %>% unlist()
  core_genes_symbol <- bitr(core_genes_entrez, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
  heatmap_data <- expr_filtered[rownames(expr_filtered) %in% core_genes_symbol, ]

  if(nrow(heatmap_data) > 1) {
    annotation_col <- data.frame(RiskGroup = risk_filtered$risk_group)
    rownames(annotation_col) <- risk_filtered$sample_id
    risk_colors <- list(RiskGroup = c("Low" = "#1F77B4", "High" = "#D62728"))
    heatmap_colors <- colorRampPalette(c("#D62728", "white", "#1F77B4"))(100) 
    svglite("Results/Fig_GSEA_Heatmap_HighRisk.svg", width = 12, height = 10)
    pheatmap::pheatmap(
      heatmap_data,
      scale = "row",
      color = heatmap_colors,
      annotation_col = annotation_col,
      annotation_colors = risk_colors,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      fontsize_row = 8,
      breaks = seq(-5, 5, length.out = 101),
      main = pathway_name
    )
    dev.off()
  } else {
    cat("警告: '", pathway_name, "' 通路的热图数据行数不足 (", nrow(heatmap_data), " 行)，已跳过绘图。\n")
  }
}

# 低风险组
low_risk_top_pathway <- sig_gsea %>% filter(NES < 0) %>% slice_min(NES, n = 1)
if (nrow(low_risk_top_pathway) > 0) {
  pathway_id <- low_risk_top_pathway$ID
  pathway_name <- low_risk_top_pathway$Description
  cat("\n正在为低风险组Top通路生成热图:", pathway_name, "\n")
  core_genes_entrez <- sig_gsea$core_enrichment[sig_gsea$ID == pathway_id] %>% strsplit("/") %>% unlist()
  core_genes_symbol <- bitr(core_genes_entrez, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
  heatmap_data <- expr_filtered[rownames(expr_filtered) %in% core_genes_symbol, ]
  gene_fontsize <- if (length(core_genes_symbol) > 100) 6 else 8
  if(nrow(heatmap_data) > 1) {
    annotation_col <- data.frame(RiskGroup = risk_filtered$risk_group)
    rownames(annotation_col) <- risk_filtered$sample_id
    risk_colors <- list(RiskGroup = c("Low" = "#1F77B4", "High" = "#D62728"))
    heatmap_colors <- colorRampPalette(c("#D62728", "white", "#1F77B4"))(100)
    svglite("Results/Fig_GSEA_Heatmap_LowRisk.svg", width = 12, height = 10)
    pheatmap::pheatmap(
      heatmap_data,
      scale = "row",
      color = heatmap_colors,
      annotation_col = annotation_col,
      annotation_colors = risk_colors,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      fontsize_row = gene_fontsize,
      breaks = seq(-5, 5, length.out = 101),
      main = pathway_name
    )
    dev.off()
  } else {
    cat("警告: '", pathway_name, "' 通路的热图数据行数不足 (", nrow(heatmap_data), " 行)，已跳过绘图。\n")
  }
}
