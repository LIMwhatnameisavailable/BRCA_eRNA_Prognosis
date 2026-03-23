library(dplyr)
library(survival)
dir.create("Results", showWarnings = FALSE, recursive = TRUE)

prepare_cox_data <- function(data_subset, clinical_ref) {
  data_subset$patient_id <- substr(rownames(data_subset), 1, 12)
  merged_data <- inner_join(data_subset, clinical_ref, by = "patient_id")
  merged_data$T_stage <- gsub("[abc]", "", merged_data$T_stage)
  merged_data$N_stage <- gsub("[abc]|\\(mol\\+\\)|mi", "", merged_data$N_stage)
  merged_data$M_stage <- gsub("[abc]", "", merged_data$M_stage)
  
  final_data <- merged_data %>%
    filter(
      T_stage %in% c("T1", "T2", "T3", "T4"),
      N_stage %in% c("N0", "N1", "N2", "N3"),
      M_stage %in% c("M0", "M1")
    ) %>%
    na.omit()

  final_data$T_stage_num <- as.numeric(as.factor(final_data$T_stage))
  final_data$N_stage_num <- as.numeric(as.factor(final_data$N_stage))
  final_data$M_stage_num <- as.numeric(as.factor(final_data$M_stage))
  return(final_data)
}

# 执行 Cox 回归并生成表格
run_cox_analysis <- function(data, cohort_name) {
  covariates <- c(
    "Age" = "age",
    "T Stage" = "T_stage_num",
    "N Stage" = "N_stage_num",
    "M Stage" = "M_stage_num",
    "Risk Group" = "risk_group" 
  )
  data$risk_group <- factor(data$risk_group, levels = c("Low", "High"))
  
  res_df <- data.frame()
  
  # 多因素模型
  multi_form <- as.formula(paste("Surv(OS.time_month, OS) ~", paste(unname(covariates), collapse = "+")))
  multi_fit <- coxph(multi_form, data = data)
  multi_sum <- summary(multi_fit)
  
  for(label in names(covariates)) {
    var <- covariates[[label]]
    
    # 单因素
    uni_fit <- coxph(as.formula(paste("Surv(OS.time_month, OS) ~", var)), data = data)
    uni_sum <- summary(uni_fit)
    uni_hr <- sprintf("%.3f", uni_sum$conf.int[1])
    uni_ci <- paste0("(", sprintf("%.3f", uni_sum$conf.int[3]), "-", sprintf("%.3f", uni_sum$conf.int[4]), ")")
    uni_p <- uni_sum$coefficients[5]
    uni_p_str <- ifelse(uni_p < 0.001, "< 0.001", sprintf("%.3f", uni_p))
    
    # 多因素
    target_row <- grep(var, rownames(multi_sum$coefficients))[1]
    if(!is.na(target_row)) {
      multi_hr_val <- multi_sum$conf.int[target_row, 1]
      multi_ci_low <- multi_sum$conf.int[target_row, 3]
      multi_ci_high <- multi_sum$conf.int[target_row, 4]
      multi_p_val <- multi_sum$coefficients[target_row, 5]
      multi_hr_str <- sprintf("%.3f", multi_hr_val)
      multi_ci_str <- paste0("(", sprintf("%.3f", multi_ci_low), "-", sprintf("%.3f", multi_ci_high), ")")
      multi_p_str <- ifelse(multi_p_val < 0.001, "< 0.001", sprintf("%.3f", multi_p_val))
    } else {
      multi_hr_str <- "NA"; multi_ci_str <- ""; multi_p_str <- "NA"
    }
    
    res_df <- rbind(res_df, data.frame(
      Characteristic = label,
      `Univariate HR (95% CI)` = paste(uni_hr, uni_ci),
      `Univariate P` = uni_p_str,
      `Multivariate HR (95% CI)` = paste(multi_hr_str, multi_ci_str),
      `Multivariate P` = multi_p_str,
      check.names = FALSE
    ))
  }
  return(res_df)
}

clinical_ref <- data_clinical %>%
  dplyr::select(
    patient_id = submitter_id, 
    age = `age_at_earliest_diagnosis_in_years.diagnoses.xena_derived`,
    T_stage = `ajcc_pathologic_t.diagnoses`,
    N_stage = `ajcc_pathologic_n.diagnoses`,
    M_stage = `ajcc_pathologic_m.diagnoses`
  ) %>%
  mutate(patient_id = substr(patient_id, 1, 12)) %>%
  distinct(patient_id, .keep_all = TRUE)

# 处理训练集
cutoff_score <- median(train_data$risk_score)
train_data$risk_group <- ifelse(train_data$risk_score > cutoff_score, "High", "Low")

# 清洗与合并
cox_train_data <- prepare_cox_data(train_data, clinical_ref)
cat("训练集纳入样本数:", nrow(cox_train_data), " (死亡事件:", sum(cox_train_data$OS==1), ")\n")

# 生成表格
table_train <- run_cox_analysis(cox_train_data, "Training Cohort")
print("========= 训练集 Cox 回归表 =========")
print(table_train)
write.csv(table_train, "Results/Table2_Training_Cox.csv", row.names = FALSE)


# 处理测试集
test_cutoff <- median(test_data$risk_score)
test_data$risk_group <- ifelse(test_data$risk_score > test_cutoff, "High", "Low")
cox_test_data <- prepare_cox_data(test_data, clinical_ref)
cat("测试集纳入样本数:", nrow(cox_test_data), " (死亡事件:", sum(cox_test_data$OS==1), ")\n")
if(sum(cox_test_data$OS==1) >= 10) { 
  table_test <- run_cox_analysis(cox_test_data, "Testing Cohort")
  print("========= 测试集 Cox 回归表 =========")
  print(table_test)
  write.csv(table_test, "Results/Table2_Testing_Cox.csv", row.names = FALSE)
} else {
  cat("警告：测试集死亡事件过少，多因素回归结果可能不可靠或报错。\n")
}

