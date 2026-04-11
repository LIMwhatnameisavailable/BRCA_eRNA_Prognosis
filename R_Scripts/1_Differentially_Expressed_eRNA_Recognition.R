library(data.table)
data <- fread("Data_Source/TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv")
data_clinical <- fread("Data_Source/TCGA-BRCA.clinical.tsv")
data_gene <- fread("Data_Source/TCGA-BRCA.methylation27.tsv")
library(survival) 
library(survminer)
library(ggrastr)

# Format column names to TCGA standard and set row names to eRNA names
original_col <- colnames(data)
new_col <- ifelse(
  grepl("_tumor", original_col), gsub("_tumor", "-01A", original_col),
  ifelse(
    grepl("_normal", original_col),gsub("_normal", "-11A", original_col),
    original_col)
)
colnames(data) <- new_col

# Predefine tumor and normal control groups
tumor_col <- colnames(data)[which(substr(colnames(data),14,14) == 0)]
length(tumor_col)
normal_col <- colnames(data)[which(substr(colnames(data),14,14) == 1)]
length(normal_col)
data <- as.data.frame(data)
rownames(data) <- data[, 1]
data <- data[, -1]

# Load survival and clinical data
library(readxl)
survival <- read.delim("Data_Source/TCGA-BRCA.survival.tsv", fileEncoding = "UTF-16")
clinical <- fread("Data_Source/TCGA-BRCA.clinical.tsv")

# Deduplicate patients with both normal and tumor samples (1076 patients)
survival <- unique(survival[,c("X_PATIENT","OS.time","OS")]) 
survival <- as.data.frame(survival)  
rownames(survival) <- survival[, 1]
survival <- survival[, -1] 
clinical <- as.data.frame(clinical) 
rownames(clinical) <- clinical[, 1]
clinical <- clinical[, -1]

# Convert survival days to months and round to two decimal places
survival$OS.time <- round(as.numeric(survival$OS.time)/30,2) 
colnames(survival)[1] <- "OS.time/month"
survival$OS <- as.numeric(survival$OS) 

# Define sample groups and construct design matrix
group <- rep("normal", ncol(data)) 
group[colnames(data) %in% tumor_col] <- "tumor"
table(group)
group <- factor(group)
design <- model.matrix(~0 + group)

# Set design matrix column names to tumor/normal and row names to samples
colnames(design) <- levels(group)
rownames(design) <- colnames(data)

# Split data into training and testing cohorts
library(caret)
set.seed(123)
survival <- as.data.frame(survival)
train_index <- createDataPartition(
  y = survival$'OS.time/month',  
  p = 0.7,         
  list = FALSE
)
class(train_index)
train_cohort <- as.vector(train_index)
train_cohort <- survival[train_index, ]
test_cohort <- survival[-train_index, ]
cat("Number of Training Cohort:", nrow(train_cohort), "\n",
    "Number of Testing Cohort:", nrow(test_cohort))

# Calculate and test mean survival time between cohorts
train_mean <- mean(train_cohort$`OS.time/month`, na.rm = TRUE)
test_mean <- mean(test_cohort$`OS.time/month`, na.rm = TRUE)
t_test <- t.test(train_cohort$`OS.time/month`, test_cohort$`OS.time/month`)
cat("Median OStime of Training Cohort:", round(train_mean, 2), "month\n",
    "Median OStime of Testing Cohort:", round(test_mean, 2), "month\n",
    "p Value:", format.pval(t_test$p.value, digits = 3))

# Filter out lowly expressed eRNAs
library(limma)
keep <- rowSums(data >= 0.3) >= 500  
data.filtered <- data[keep, ]
cat("Original:", nrow(data), "\nFiltered:",nrow(data.filtered))

# Optimize distribution via log2 transformation
data.filtered <- log2(data.filtered + 1)

# Perform differential expression analysis
fit <- lmFit(data.filtered,design)
contrast.matrix <- makeContrasts(tumor_vs_normal = tumor-normal,levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# Adjust variances using empirical Bayes method
fit2 <- eBayes(fit2)
results <- topTable(fit2,coef = "tumor_vs_normal",number = Inf,adjust.method = "BH")

# Filter upregulated and downregulated eRNAs
DEE <- results[results$adj.P.Val < 0.001 & abs(results$logFC) > 0.8,]
dim(DEE)
DEE.up <- results[results$adj.P.Val<0.001 & results$logFC > 0.8,] 
DEE.down <- results[results$adj.P.Val<0.001 & results$logFC < -0.8,]
dim(DEE.up)
dim(DEE.down)

results$diff <- ifelse(
  results$adj.P.Val < 0.001 & abs(results$logFC) > 0.8,
  ifelse(results$logFC > 0.8, "Up", "Down"),
  "NotSig"
)

library(ggplot2)
library(ggview)
library(ggrepel)

# Generate volcano plot
volcano <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), 
                               color = factor(diff,levels = c("Up","Down","NotSig")))) +
  ggrastr::geom_point_rast(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    name = "Expression",
    values = c("Up" = "#E64B35", "Down" = "#457B9D", "NotSig" = "#BDBDBD"),
    labels = c("Up (Tumor High)", "Down (Tumor Low)", "Not Significant")
  ) +
  geom_vline(xintercept = c(-0.8, 0.8), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  labs(
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("Adjusted P-value")),
    title = NULL
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "top", panel.grid.major = element_blank())
print(volcano)
ggsave(
  "Results/Fig_Volcano.svg",  
  plot = volcano,        
  device = "svg",       
  width = 7,           
  height = 6,            
)


library(ComplexHeatmap)
library(circlize) 
library(dplyr)  

rownames(data.filtered) <- gsub(":", "_", rownames(data.filtered))
rownames(DEE.up) <- gsub(":", "_", rownames(DEE.up))
rownames(DEE.down) <- gsub(":", "_", rownames(DEE.down))
top_genes <- rownames(rbind(
  DEE.up[order(DEE.up$logFC, decreasing = TRUE), ][1:50, ],
  DEE.down[order(DEE.down$logFC), ][1:50, ]
))
heatmap_data <- data.filtered[top_genes, ]
heatmap_data_clean <- heatmap_data %>%
  as.matrix() %>%
  {.[!apply(., 1, function(x) any(is.infinite(x) | sd(x)==0)), ]}
annotation_col_df <- data.frame(
  Group = ifelse(colnames(heatmap_data_clean) %in% tumor_col, "T", "N"),
  row.names = colnames(heatmap_data_clean)
)
ann_colors <- list(Group = c(T = "#E63946", N = "#457B9D"))

mat <- heatmap_data_clean
row_means <- apply(mat, 1, mean, na.rm = TRUE)
row_sds <- apply(mat, 1, sd, na.rm = TRUE)
row_sds[row_sds == 0] <- 1 
mat_scaled <- (mat - row_means) / row_sds
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[is.infinite(mat_scaled)] <- 0

top_anno <- HeatmapAnnotation(
  Group = annotation_col_df$Group,
  col = ann_colors,
  show_annotation_name = FALSE 
)

col_fun <- colorRamp2(c(-3, 0, 3), c("royalblue", "white", "red"))

ht_plot <- Heatmap(
  mat_scaled, 
  name = "Z-score",
  col = col_fun,     

  use_raster = TRUE,       
  raster_quality = 5,      
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = top_anno,
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-3, 0, 3),
    labels = c("Low", "Mid", "High")
  )
)

svg("Results/Fig_Heatmap.svg", width = 8, height = 7)
draw(ht_plot) 
dev.off()


# Filter samples with both expression and survival data
head(tumor_col)
head(DEE)
head(data.filtered)
rownames(DEE) <- gsub(":", "_", rownames(DEE))
data.filtered.tumor <- data.filtered[rownames(DEE),tumor_col]
colnames(data.filtered.tumor) <- substr(colnames(data.filtered.tumor),1,12)  # Standardize format
combined_data <- survival
overlap.samples <- intersect(rownames(train_cohort),colnames(data.filtered.tumor))
overlap.samples_combined <- intersect(rownames(combined_data),colnames(data.filtered.tumor))
length(overlap.samples)  
length(overlap.samples_combined) # 1073 samples, 753 in training set
survival.addExp <- cbind(train_cohort[overlap.samples,],
                         t(data.filtered.tumor)[overlap.samples,]) 
train_data <- survival.addExp[,c("OS.time/month", "OS", rownames(DEE))]
overlap_test.samples <- intersect(rownames(test_cohort),colnames(data.filtered.tumor))
length(overlap_test.samples)
survival_test.addExp <- cbind(test_cohort[overlap_test.samples,],
                              t(data.filtered.tumor)[overlap_test.samples,])
test_data <- survival_test.addExp[,c("OS.time/month", "OS", rownames(DEE))]
combined_data <- cbind(survival[overlap.samples_combined,],
                       t(data.filtered.tumor)[overlap.samples_combined,])
combined_data <- combined_data[,c("OS.time/month", "OS", rownames(DEE))]

# Perform univariate Cox regression analysis
library(stringr)
univ_results <- data.frame()
colnames(train_data) <- gsub(":", "_", colnames(train_data))
colnames(train_data) <- gsub("/month", "_month", colnames(train_data))
rownames(DEE) <- gsub(":", "_", rownames(DEE))

for(eRNA in rownames(DEE)){
  formula <- as.formula(paste("Surv(OS.time_month, OS) ~", eRNA, "", sep = ""))
  cox_model <- coxph(formula, data = train_data)
  cox_summary <- summary(cox_model)
  
  univ_results <- rbind(univ_results, data.frame(
    eRNA = eRNA,
    HR = cox_summary$conf.int[1],
    Lower_95CI = cox_summary$conf.int[3],
    Upper_95CI = cox_summary$conf.int[4],
    p.value = cox_summary$coefficients[5]
  ))
}
univ_results <- univ_results %>% arrange(p.value)
cat("单变量Cox回归完成。\n")

# Filter data and format eRNA names
sig_eRNAs_data <- univ_results %>% 
  filter(p.value < 0.01) %>%
  mutate(eRNA_label = str_replace(eRNA, "_", ":")) %>%
  
  # Create grouping column based on HR values
  mutate(EffectGroup = ifelse(HR > 1, 
                              "Risk (P<0.01)", 
                              "Protective (P<0.01)")) %>%
  
  # Sort by HR
  arrange(HR) %>%
  mutate(row_id = row_number()) %>%
  mutate(stripe_group = factor(row_id %% 2)) %>%
  mutate(eRNA_label = factor(eRNA_label, levels = eRNA_label))

cat("筛选出", nrow(sig_eRNAs_data), "个显著eRNA用于绘制森林图。\n")

forest_plot <- ggplot(sig_eRNAs_data, aes(x = HR, y = eRNA_label)) +
  geom_rect(
    aes(ymin = as.numeric(eRNA_label) - 0.5, 
        ymax = as.numeric(eRNA_label) + 0.5,
        xmin = 0.3, xmax = 2.1, 
        fill = stripe_group), # Fill using the created binary grouping
    color = NA # Remove rectangle borders
  ) +
  scale_fill_manual(
    values = c("0" = "white", "1" = "grey95"),
    guide = "none" # Hide legend
  ) +
  geom_errorbarh(aes(xmin = Lower_95CI, xmax = Upper_95CI), 
                 height = 0.25,      
                 linewidth = 0.5,    
                 color = "black") +  
  geom_point(aes(color = EffectGroup), 
             size = 5,   
             shape = 16) + 

  scale_color_manual(
    name = "Effect Type", 
    values = c("Risk (P<0.01)" = "#E63946",  
               "Protective (P<0.01)" = "#0072B2"), 
    labels = c("Protective (P<0.01)", "Risk (P<0.01)")
  ) +
  geom_vline(xintercept = 1, 
             linetype = "dashed",  
             color = "black", 
             linewidth = 0.6) +
  scale_x_log10(
    breaks = c(0.4, 0.6, 1.0, 1.4, 1.8), 
    labels = c("0.4", "0.6", "1.0", "1.4", "1.8"),
    limits = c(0.3, 2.1) 
  ) +
  labs(
    title = NULL,
    x = "Hazard Ratio (HR) with 95% CI",
    y = NULL  
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top" 
  )

print(forest_plot)

ggsave(
  "Results/Fig_Forest_Plot_Cox.svg",
  plot = forest_plot,
  device = "svg",
  width = 8,  
  height = 7, 
)


# Perform LASSO regression
library(glmnet)
library(ggsci) 
library(reshape2) 
library(patchwork) 
sig_eRNAs <- sig_eRNAs_data$eRNA
x <- as.matrix(train_data[, sig_eRNAs])
y <- Surv(train_data$OS.time_month, train_data$OS) 

# Fit LASSO model
fit <- glmnet(x, y, family = "cox", alpha = 1, nlamdba = 300)

# Plot cross-validation error curve
set.seed(123)
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, 
                   nfolds = 10,  
                   type.measure = "deviance") 


coef_path <- as.matrix(fit$beta)
melt_coef <- melt(coef_path, varnames = c("Predictor", "Step"), value.name = "Coefficient")
melt_coef$Lambda <- rep(fit$lambda, each = nrow(coef_path))

# Dynamically assign colors
n_predictors <- nrow(coef_path)     
npg_palette <- pal_npg()(10)       
my_colors <- rep(npg_palette, length.out = n_predictors) 

(p1 <- ggplot(melt_coef, aes(x = log(Lambda), y = Coefficient, 
                             group = Predictor, color = Predictor)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    geom_vline(xintercept = log(cvfit$lambda.min), 
               linetype = "dashed", color = "grey40") +
  
    scale_color_manual(values = my_colors) + 
    
    labs(x = expression(log(lambda)), 
         y = "Coefficients",
         title = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",  # Hide legend
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(color = "black", size = 11)
    ))
cv_data <- data.frame(
  lambda = cvfit$lambda,
  cvm = cvfit$cvm,
  cvsd = cvfit$cvsd
)

(p2 <- ggplot(cv_data, aes(x = log(lambda))) +
    geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), 
                  width = 0.1, color = "grey80") +
    geom_point(aes(y = cvm), size = 2, color = npg_palette[1]) + 
    
    geom_vline(xintercept = log(cvfit$lambda.min), 
               linetype = "dashed", color = "grey40") +
    labs(x = expression(log(lambda)),
         y = "Partial Likelihood Deviance",
         title = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(color = "black", size = 11)
    ))
final_plot <- p1 / p2
plot(final_plot) 

ggsave(
  "Results/Fig_LASSO.svg",
  plot = final_plot,
  device = "svg",
  width = 7,  
  height = 8, 
)

            
# Generate clinical characteristics summary table
library(gtsummary)
library(flextable)
library(forcats) 

# Extract IDs for the 1073 target samples
train_patients <- overlap.samples
test_patients <- overlap_test.samples

# Process and deduplicate clinical data
clinical_df <- clinical
clinical_df$Patient_ID <- substr(rownames(clinical_df), 1, 12)
clinical_df <- clinical_df[!duplicated(clinical_df$Patient_ID), ]

# Filter target samples and assign cohort labels
clinical_filtered <- clinical_df[clinical_df$Patient_ID %in% c(train_patients, test_patients), ]
clinical_filtered$Cohort <- ifelse(clinical_filtered$Patient_ID %in% train_patients, 
                                   "Training Cohort", 
                                   "Testing Cohort")

clinical_ready <- clinical_filtered %>%
  select(
    Cohort,
    Age = age_at_index.demographic,                   
    Stage = ajcc_pathologic_stage.diagnoses,          
    T_classification = ajcc_pathologic_t.diagnoses,   
    N_classification = ajcc_pathologic_n.diagnoses,   
    M_classification = ajcc_pathologic_m.diagnoses,
    Primary_diagnosis = primary_diagnosis.diagnoses
  ) %>%
  mutate(
    Age = as.numeric(Age),
    
    # Clean Stage variable retaining Stage X and setting blanks to NA
    Stage = case_when(
      Stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
      Stage %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
      Stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
      Stage %in% c("Stage IV") ~ "Stage IV",
      Stage == "Stage X" ~ "Stage X",
      TRUE ~ NA_character_ # Convert remaining empty strings to NA
    ),
    
    # Clean T classification retaining TX
    T_classification = case_when(
      T_classification %in% c("T1", "T1a", "T1b", "T1c") ~ "T1",
      T_classification %in% c("T2", "T2a", "T2b") ~ "T2",
      T_classification %in% c("T3", "T3a") ~ "T3",
      T_classification %in% c("T4", "T4b", "T4d") ~ "T4",
      T_classification == "TX" ~ "TX",
      TRUE ~ NA_character_
    ),
    
    # Clean N classification retaining NX
    N_classification = case_when(
      N_classification %in% c("N0", "N0 (i-)", "N0 (i+)", "N0 (mol+)") ~ "N0",
      N_classification %in% c("N1", "N1a", "N1b", "N1c", "N1mi") ~ "N1",
      N_classification %in% c("N2", "N2a") ~ "N2",
      N_classification %in% c("N3", "N3a", "N3b", "N3c") ~ "N3",
      N_classification == "NX" ~ "NX",
      TRUE ~ NA_character_
    ),
    
    # Clean M classification retaining MX
    M_classification = case_when(
      M_classification %in% c("M0", "cM0 (i+)") ~ "M0",
      M_classification %in% c("M1") ~ "M1",
      M_classification == "MX" ~ "MX",
      TRUE ~ NA_character_
    ),
    
    # Convert empty strings in Primary diagnosis to NA
    Primary_diagnosis = na_if(Primary_diagnosis, "")
  ) %>%
  mutate(
    # Convert to factors and drop unused levels
    Stage = fct_drop(factor(Stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV", "Stage X"))),
    T_classification = fct_drop(factor(T_classification, levels = c("T1", "T2", "T3", "T4", "TX"))),
    N_classification = fct_drop(factor(N_classification, levels = c("N0", "N1", "N2", "N3", "NX"))),
    M_classification = fct_drop(factor(M_classification, levels = c("M0", "M1", "MX"))),
  
    Primary_diagnosis = fct_lump_n(Primary_diagnosis, n = 3, other_level = "Other")
  )

# Generate Table 1
theme_gtsummary_journal(journal = "jama")

table1 <- clinical_ready %>%
  tbl_summary(
    by = Cohort, 
    statistic = list(
      all_continuous() ~ "{median} ({min}-{max})", 
      all_categorical() ~ "{n}"                    
    ),
    missing = "no" 
  ) %>%
  add_overall(col_label = "**Combined Cohort**\n(n = {N})") %>% 
  add_p( 
    test = list(
      all_continuous() ~ "wilcox.test", 
      all_categorical() ~ "chisq.test"  
    )
  ) %>% 
  modify_header(label ~ "**Variables**") %>% 
  bold_labels() 

# Print preview
print(table1)

# Export to Word document
table1 %>% 
  as_flex_table() %>% 
  flextable::save_as_docx(path = "Results/Table1_Clinical_Characteristics.docx")
