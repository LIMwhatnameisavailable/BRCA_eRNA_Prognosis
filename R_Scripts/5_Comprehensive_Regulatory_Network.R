# ==============================================================================
#                          Step 0: 环境设置
# ==============================================================================
rm(list = ls()); gc() 
options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)
library(tibble)
library(limma)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicRanges)
library(rtracklayer)        
library(BSgenome.Hsapiens.UCSC.hg38)
library(MotifDb)
library(motifbreakR)
library(VariantAnnotation)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(gridExtra)

mRNA_path     <- "Data_Source/TCGA-BRCA.star_fpkm.tsv.gz"
eRNA_path     <- "Data_Source/TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv"
hic_file      <- "Data_Source/loop_info.csv"
clinical_file <- "Data_Source/pcawg_donor_clinical_August2016_v9"
xena_mut_file <- "Data_Source/PCAWG_WGS_mutations.tsv.gz"
xena_cnv_file  <- "Data_Source/TCGA.BRCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
xena_meth_file <- "Data_Source/TCGA.BRCA.sampleMap_HumanMethylation450.gz" 
xena_probe_map <- "Data_Source/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy" 

original_file <- "Data_Source/hg19ToHg38.over.chain.gz"
fixed_file <- "Data_Source/hg19ToHg38.fixed_universal.chain"

lines <- readLines(original_file)
lines_fixed <- gsub(" +", " ", lines) 
lines_fixed <- gsub("^[[:space:]]+|[[:space:]]+$", "", lines_fixed) 
lines_fixed <- gsub("\t", " ", lines_fixed)
writeLines(lines_fixed, fixed_file)

cat("文件已修复并保存为:", fixed_file, "\n")
chain <- import.chain(fixed_file)


# 定义原始 eRNA 列表
target_eRNAs_hg19 <- c(
  "chr1:155158995", "chr3:11236700", "chr8:22624675",
  "chr3:138070534", "chr9:71398719", "chr9:114689796",
  "chr10:5531356", "chr9:71398939", "chr12:13371038",
  "chr10:5528926"  
)

# 候选 TF 列表
candidate_tfs <- c(
  "GATA3", "FOXA1", "ESR1", "PGR", "ERBB2", "BRCA1", "TP53", "MYC",
  "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "SOX2", "POU5F1", "KLF4", "NANOG",
  "E2F1", "E2F2", "E2F3", "E2F4", "E2F6", "E2F7", "E2F8", "MKI67", "CCND1",
  "NFKB1", "RELA", "RELB", "STAT1", "STAT3", "STAT5A", "IRF1", "IRF3", "JUN", "FOS", "JUND", "FOSL1",
  "SP1", "EP300", "CTCF", "YY1", "MAX", "HIF1A", "AR", "SMAD3", "SMAD4", "TEAD1", "TEAD4",
  "ETS1", "ELK1", "RUNX1", "GATA1", "CEBPB", "GRHL2", "OVOL2", "PBX1", "XBP1", "ELF5", "TFAP2A", "TFAP2C", "NR3C1"
)

# ==============================================================================
#                 Step 1: eRNA 坐标转换 (hg19 -> hg38)
# ==============================================================================

# 解析原始坐标
parse_coords <- function(x) {
  parts <- strsplit(x, ":")[[1]]
  return(data.frame(chr=parts[1], start=as.numeric(parts[2]), end=as.numeric(parts[2]), name=x))
}
eRNA_hg19_df <- do.call(rbind, lapply(target_eRNAs_hg19, parse_coords))
eRNA_hg19_gr <- GRanges(seqnames = eRNA_hg19_df$chr, 
                        ranges = IRanges(start = eRNA_hg19_df$start, end = eRNA_hg19_df$end))
names(eRNA_hg19_gr) <- eRNA_hg19_df$name

# 执行转换
eRNA_hg38_list <- liftOver(eRNA_hg19_gr, chain)
eRNA_hg38_gr <- unlist(eRNA_hg38_list)

# 定义核心区域 (hg38)
eRNA_gr <- resize(eRNA_hg38_gr, width = 1000, fix = "center")

cat("转换成功! ", length(eRNA_hg38_gr), "/ 10 个 eRNA 成功映射到 hg38。\n")
valid_eRNAs <- names(eRNA_gr) 

# ==============================================================================
#                        Step 2: 表达量数据处理
# ==============================================================================

mRNA_df <- fread(mRNA_path, data.table = FALSE)
rownames(mRNA_df) <- mRNA_df[,1]; mRNA_df <- mRNA_df[,-1]
eRNA_df <- fread(eRNA_path, data.table = FALSE)
if(ncol(eRNA_df) == length(colnames(eRNA_df)) + 1) {
  rownames(eRNA_df) <- eRNA_df[,1]; eRNA_df <- eRNA_df[,-1]
} else {
  rownames(eRNA_df) <- eRNA_df[,1]; eRNA_df <- eRNA_df[,-1]
}

# 提取目标 eRNA
eRNA_subset <- eRNA_df[rownames(eRNA_df) %in% valid_eRNAs, ]

# 样本对齐
is_tumor_e <- grep("_tumor|_tu", colnames(eRNA_subset), ignore.case = TRUE)
eRNA_tumor <- eRNA_subset[, is_tumor_e]
colnames(eRNA_tumor) <- substr(gsub("\\.", "-", colnames(eRNA_tumor)), 1, 12)

sample_codes <- as.numeric(substr(colnames(mRNA_df), 14, 15))
is_tumor_m <- which(sample_codes < 10)
mRNA_tumor <- mRNA_df[, is_tumor_m]
colnames(mRNA_tumor) <- substr(colnames(mRNA_tumor), 1, 12)

common_samples <- intersect(colnames(eRNA_tumor), colnames(mRNA_tumor))
final_eRNA <- eRNA_tumor[, common_samples]
final_mRNA <- mRNA_tumor[, common_samples]
cat("匹配到的 Tumor 样本数:", length(common_samples), "\n")
rm(mRNA_df, eRNA_df); gc()

# ID 转换 (mRNA Ensembl -> Symbol)
ensembl_ids <- gsub("\\..*", "", rownames(final_mRNA))
gene_map <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
match_idx <- match(ensembl_ids, gene_map$ENSEMBL)
valid_idx <- !is.na(match_idx)
mRNA_matrix <- as.matrix(final_mRNA[valid_idx, ])
gene_symbols <- gene_map$SYMBOL[match_idx[valid_idx]]
mRNA_agg <- avereps(mRNA_matrix, ID = gene_symbols)

clean_eRNA_ids <- gsub("-|_", ":", rownames(final_eRNA))
rownames(final_eRNA) <- clean_eRNA_ids
exp_matrix <- rbind(as.matrix(final_eRNA), mRNA_agg)

valid_tfs <- intersect(candidate_tfs, rownames(exp_matrix))
valid_genes <- setdiff(rownames(exp_matrix), c(clean_eRNA_ids, valid_tfs))


# ==============================================================================
#                        Step 3: 突变机制挖掘
# ==============================================================================

# 初始化结果
motif_breaks <- NULL
mut_overlap_counts <- setNames(rep(0, length(valid_eRNAs)), valid_eRNAs)
if(file.exists(clinical_file) & file.exists(xena_mut_file)) {
  clin_df <- fread(clinical_file, data.table = FALSE)
  if("donor_diagnosis_icd10" %in% colnames(clin_df)) {
    brca_samples <- clin_df %>% filter(grepl("C50", donor_diagnosis_icd10)) %>% pull(icgc_donor_id)
  } else { brca_samples <- clin_df[,1] }
  
  header_preview <- names(fread(xena_mut_file, nrows = 0))
  possible_id_names <- c("Sample", "sample", "icgc_donor_id", "donor_id", "Tumor_Sample_Barcode")
  actual_id_col <- intersect(possible_id_names, header_preview)[1]
  
  mut_df <- fread(xena_mut_file, select = c("chr", "start", "end", "reference", "alt", actual_id_col))
  colnames(mut_df) <- c("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")
  
  mut_df$Chromosome <- as.character(mut_df$Chromosome)
  if(!any(grepl("chr", head(mut_df$Chromosome)))) mut_df$Chromosome <- paste0("chr", mut_df$Chromosome)
  
  mut_df <- mut_df %>% 
    filter(Tumor_Sample_Barcode %in% brca_samples) %>%
    filter(Chromosome %in% paste0("chr", c(1:22, "X", "Y")))
  
  # 坐标转换与重叠
  mut_gr_hg19 <- GRanges(seqnames = mut_df$Chromosome, 
                         ranges = IRanges(start = mut_df$Start_Position, end = mut_df$End_Position),
                         mcols = mut_df[, c("Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")])
  mut_gr_hg38 <- unlist(liftOver(mut_gr_hg19, chain))
  
  overlaps <- findOverlaps(mut_gr_hg38, eRNA_gr)
  
  if(length(overlaps) > 0) {
    candidate_muts <- mut_gr_hg38[queryHits(overlaps)]
    target_names <- names(eRNA_gr)[subjectHits(overlaps)]
    
    counts <- table(target_names)
    mut_overlap_counts[names(counts)] <- as.numeric(counts)
    cat("发现", length(candidate_muts), "个位于 eRNA 区域的 WGS 突变 (hg38)！\n")
    
    if(length(valid_tfs) > 0) {
      cat("   -> 正在准备 Motif PWM 列表...\n")
      library(BSgenome.Hsapiens.UCSC.hg38)
      library(TFBSTools)
      
      motif_list_raw <- query(MotifDb, andStrings=c("sapiens", "jaspar2018"), orStrings=valid_tfs)
      pwm_list_safe <- list() 
      
      if(length(motif_list_raw) > 0) {
        raw_list <- as.list(motif_list_raw)
        for(i in seq_along(raw_list)) {
          mat <- raw_list[[i]]
          name <- names(raw_list)[i]
          if(is.matrix(mat) && nrow(mat) == 4 && !any(is.na(mat))) {
            tryCatch({
              pfm <- PFMatrix(ID = name, name = name, profileMatrix = mat)
              pwm <- toPWM(pfm, pseudocounts = 0.8)
              # 存入普通列表
              pwm_list_safe[[name]] <- pwm
            }, error = function(e) {})
          }
        }
      }
      cat("   -> 准备了", length(pwm_list_safe), "个 PWM 矩阵用于扫描。\n")
      
      # 提取序列
      bsg <- BSgenome.Hsapiens.UCSC.hg38
      scan_window <- 41 
      center_pos <- ceiling(scan_window / 2)
      
      seqlevelsStyle(candidate_muts) <- "UCSC"
      common_seq <- intersect(seqlevels(candidate_muts), seqlevels(bsg))
      candidate_muts <- keepSeqlevels(candidate_muts, common_seq, pruning.mode="coarse")
      
      ranges_expanded <- resize(candidate_muts, width = scan_window, fix = "center")
      ref_seqs <- getSeq(bsg, ranges_expanded)
      
      alt_alleles <- as.character(candidate_muts$mcols.Tumor_Seq_Allele2)
      alt_seqs <- ref_seqs
      for(i in 1:length(alt_seqs)) {
        if(nchar(alt_alleles[i]) == 1) {
          subseq(alt_seqs[[i]], start=center_pos, width=1) <- DNAString(alt_alleles[i])
        }
      }
      
      
      results_list <- list()
      # 外层循环：突变
      for(i in 1:length(candidate_muts)) {
        if(nchar(alt_alleles[i]) != 1) next 
        
        curr_ref_seq <- ref_seqs[[i]]
        curr_alt_seq <- alt_seqs[[i]]
        curr_snp_pos <- paste0(seqnames(candidate_muts)[i], ":", start(candidate_muts)[i])
        curr_erna <- target_names[i]
        
        # 内层循环：Motif
        for(m_name in names(pwm_list_safe)) {
          curr_pwm <- pwm_list_safe[[m_name]]
          
          tryCatch({
            # 扫描 Ref
            hits_ref <- searchSeq(curr_pwm, curr_ref_seq, min.score="80%")
            # 扫描 Alt
            hits_alt <- searchSeq(curr_pwm, curr_alt_seq, min.score="80%")
            
            # 转换为 DF
            df_ref <- as(hits_ref, "data.frame")
            df_alt <- as(hits_alt, "data.frame")
            
            # 过滤位置
            if(nrow(df_ref) > 0) {
              df_ref <- df_ref[df_ref$start <= center_pos & (df_ref$start + df_ref$width - 1) >= center_pos, ]
            }
            if(nrow(df_alt) > 0) {
              df_alt <- df_alt[df_alt$start <= center_pos & (df_alt$start + df_alt$width - 1) >= center_pos, ]
            }
            
            # 比较
            if(nrow(df_ref) > 0) {
              # 取 Ref 最高分
              score_ref <- max(df_ref$score)
              max_score <- curr_pwm@min.score 
              
              if(nrow(df_alt) > 0) {
                score_alt <- max(df_alt$score)
              } else {
                score_alt <- -100
              }
              
              if(score_ref > 5 && (score_ref - score_alt) > 3) {
                res_row <- data.frame(
                  geneSymbol = m_name,
                  seqMatch = as.character(df_ref$siteSeqs[1]), 
                  scoreRef = score_ref,
                  scoreAlt = score_alt,
                  effect = "strong",
                  snpPos = curr_snp_pos,
                  eRNA_ID = curr_erna
                )
                results_list[[length(results_list)+1]] <- res_row
              }
            }
          }, error = function(e) {
          })
        }
      }
      
      if(length(results_list) > 0) {
        motif_breaks <- do.call(rbind, results_list)
        # 去重
        motif_breaks <- motif_breaks[!duplicated(motif_breaks[, c("geneSymbol", "snpPos")]), ]
        cat("识别到", nrow(motif_breaks), "个 Motif 破坏事件。\n")
        print(head(motif_breaks))
      } else {
        cat("分析完成，未发现显著的 Motif 破坏事件。\n")
      }
    }
  } else {
    cat("WGS 突变与 eRNA 无重叠。\n")
  }
  rm(mut_df, clin_df); gc()
  
} else {
  cat("缺失文件。\n")
}

vst_matrix <- exp_matrix
if(!exists("eRNA_gr")) {
  stop(" eRNA_gr 变量丢失，请重新运行 Step 1 代码块。")
} else {
  cat("变量检查通过。\n")
}

# ==============================================================================
#                        Step 4: 驱动机制挖掘
# ==============================================================================

# 确保前置数据存在
if(!exists("vst_matrix") | !exists("eRNA_gr")) {
  stop("错误：缺少 vst_matrix 或 eRNA_gr 变量。请先运行 Step 1-3。")
}

chain_file <- "Data_Source/hg19ToHg38.fixed_universal.chain" 

# 拷贝数变异 (CNV) 关联分析
cnv_hits <- list() # 存储显著结果
if(file.exists(xena_cnv_file)) {
  
  # 读取 CNV 数据
  cnv_df <- fread(xena_cnv_file, data.table = FALSE)
  rownames(cnv_df) <- cnv_df[, 1]
  cnv_mat <- as.matrix(cnv_df[, -1])
  
  # 处理样本名 
  colnames(cnv_mat) <- substr(colnames(cnv_mat), 1, 12)
  cnv_mat <- cnv_mat[, !duplicated(colnames(cnv_mat))]
  
  cat("   -> CNV 数据就绪: ", nrow(cnv_mat), "genes x", ncol(cnv_mat), "samples.\n")
  
  # 寻找 Proxy 基因 
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  all_genes_gr <- genes(txdb)
  
  nearest_idx <- nearest(eRNA_gr, all_genes_gr)
  nearest_gene_ids <- names(all_genes_gr)[nearest_idx]
  
  # ID转换
  gene_map <- select(org.Hs.eg.db, keys=nearest_gene_ids, columns="SYMBOL", keytype="ENTREZID")
  gene_map <- gene_map[!duplicated(gene_map$ENTREZID), ]
  
  # 映射表
  proxy_map <- data.frame(
    eRNA_ID = names(eRNA_gr),
    Proxy_Entrez = nearest_gene_ids,
    stringsAsFactors = FALSE
  )
  proxy_map <- merge(proxy_map, gene_map, by.x="Proxy_Entrez", by.y="ENTREZID", all.x=TRUE)
  proxy_map <- proxy_map[!is.na(proxy_map$SYMBOL), ]
  
  # 关联分析 (Kruskal-Wallis Test)
  expr_t <- t(vst_matrix)
  rownames(expr_t) <- substr(rownames(expr_t), 1, 12)
  
  count_sig_cnv <- 0
  
  for(i in 1:nrow(proxy_map)) {
    curr_erna <- proxy_map$eRNA_ID[i]
    curr_proxy <- proxy_map$SYMBOL[i]
    
    if(!curr_erna %in% colnames(expr_t)) next
    if(!curr_proxy %in% rownames(cnv_mat)) {
      cat("      (跳过) CNV 数据中找不到 Proxy 基因:", curr_proxy, "\n")
      next
    }
    
    # 提取并合并数据
    df_plot <- data.frame(sample = rownames(expr_t), expr = expr_t[, curr_erna])
    cnv_vals <- cnv_mat[curr_proxy, ]
    df_cnv <- data.frame(sample = names(cnv_vals), cnv = as.numeric(cnv_vals))
    merged <- merge(df_plot, df_cnv, by="sample")
    
    # 过滤无效数据 (至少2组不同CNV状态，每组>3样本)
    if(length(unique(merged$cnv)) > 1 && min(table(merged$cnv)) >= 3) {
      
      # 统计检验
      test_res <- kruskal.test(expr ~ as.factor(cnv), data = merged)
      p_val <- test_res$p.value
      
      if(p_val < 0.05) {
        cat("    显著 CNV 驱动: ", curr_erna, "(Proxy:", curr_proxy, ") P =", format.pval(p_val, digits=2), "\n")
        
        # 仅存储统计结果，不画图
        cnv_hits[[curr_erna]] <- list(
          p_value = p_val, 
          proxy = curr_proxy, 
          stats = test_res$statistic
        )
        count_sig_cnv <- count_sig_cnv + 1
      }
    }
  }
  cat("   CNV 分析完成，共发现", count_sig_cnv, "个显著驱动事件。\n")
  
} else {
  cat(" 未找到 CNV 文件，跳过 Part A。\n")
}

# 甲基化分析
library(rtracklayer)

meth_hits <- list() # 存储显著结果

if(file.exists(xena_meth_file) & file.exists(xena_probe_map) & file.exists(chain_file)) {
  cat("\n [甲基化分析] 正在加载探针信息...\n")
  
  probe_info <- fread(xena_probe_map, data.table = FALSE)
  probe_info <- probe_info[!is.na(probe_info$chrom) & probe_info$chrom != "", ]
  
  # 构建 GRanges (hg19)并转换到 hg38
  gr_450k_hg19 <- GRanges(
    seqnames = probe_info$chrom,
    ranges = IRanges(start = probe_info$chromStart, end = probe_info$chromEnd),
    Name = probe_info$id
  )
  
  cat("   -> 正在转换探针坐标 (LiftOver hg19->hg38)...\n")
  chain <- import.chain(chain_file)
  gr_450k_hg38 <- unlist(liftOver(gr_450k_hg19, chain))
  cat("      转换成功！保留了", length(gr_450k_hg38), "个探针。\n")
  
  # 匹配 eRNA 附近的探针
  # 定义区域: eRNA 本体 + 上游 2kb 启动子
  eRNA_promoters <- promoters(eRNA_gr, upstream = 2000, downstream = 200)
  eRNA_regions <- union(eRNA_gr, eRNA_promoters)
  
  overlaps <- findOverlaps(gr_450k_hg38, eRNA_regions)
  
  if(length(overlaps) > 0) {
    target_probes <- gr_450k_hg38[queryHits(overlaps)]$Name
    target_erna_names <- names(eRNA_regions)[subjectHits(overlaps)]
    
    probe_map_df <- data.frame(Probe = target_probes, eRNA = target_erna_names, stringsAsFactors = FALSE)
    cat("   -> 成功匹配到", length(unique(target_probes)), "个相关探针。\n")
    
    # 读取甲基化数据
    cat("   -> 正在读取甲基化矩阵...\n")
    meth_df <- fread(xena_meth_file, data.table = FALSE)
    rownames(meth_df) <- meth_df[,1]
    meth_df <- meth_df[,-1]
    
    # 筛选目标探针
    meth_subset <- meth_df[rownames(meth_df) %in% target_probes, , drop=FALSE]
    rm(meth_df); gc()
    
    colnames(meth_subset) <- substr(colnames(meth_subset), 1, 12)
    meth_subset <- meth_subset[, !duplicated(colnames(meth_subset))]
    
    # 相关性计算
    count_sig_meth <- 0
    
    for(i in seq_len(nrow(probe_map_df))) {
      curr_probe <- probe_map_df$Probe[i]
      curr_erna <- probe_map_df$eRNA[i]
      
      if(!curr_probe %in% rownames(meth_subset)) next
      if(!curr_erna %in% rownames(vst_matrix)) next
      
      common_samps <- intersect(colnames(meth_subset), colnames(vst_matrix))
      if(length(common_samps) < 20) next 
      
      m_vec <- as.numeric(meth_subset[curr_probe, common_samps])
      e_vec <- as.numeric(vst_matrix[curr_erna, common_samps])
      
      valid_idx <- !is.na(m_vec) & !is.na(e_vec)
      m_vec <- m_vec[valid_idx]
      e_vec <- e_vec[valid_idx]
      
      if(length(m_vec) < 20) next
      
      cor_res <- cor.test(m_vec, e_vec, method = "spearman", alternative = "less")
      
      if(cor_res$p.value < 0.05 & cor_res$estimate < -0.3) {
        cat("  显著甲基化沉默: ", curr_erna, "(Probe:", curr_probe, ") Rho =", round(cor_res$estimate, 2), "\n")
        
        # 仅存储结果
        meth_hits[[paste(curr_erna, curr_probe)]] <- list(
          eRNA = curr_erna,
          probe = curr_probe,
          rho = cor_res$estimate,
          p_val = cor_res$p.value
        )
        count_sig_meth <- count_sig_meth + 1
      }
    }
    cat("   甲基化分析完成，共发现", count_sig_meth, "个显著负相关事件。\n")
    
  } else {
    cat("    即使坐标转换后，eRNA 区域依然无甲基化探针覆盖。\n")
  }
} else {
  cat("缺失文件：请检查 chain 文件或数据文件路径。\n")
}

cat("\n Step 4 计算全部完成！\n")


# ==============================================================================
#                         Step 5: Hi-C 验证
# ==============================================================================
hic_validated_pairs <- data.frame(eRNA=character(), Gene=character(), stringsAsFactors=F)

if(file.exists(hic_file)) {
  loops_raw <- read.csv(hic_file, stringsAsFactors = FALSE)
  loops_raw <- loops_raw[loops_raw$ChIP == "H3K27ac", ]
  
  if(nrow(loops_raw) > 0) {
    # 拆分坐标
    loops_bed <- loops_raw %>%
      separate(anchor1, c("chr1", "start1", "end1"), sep = "[:-]", convert = TRUE) %>%
      separate(anchor2, c("chr2", "start2", "end2"), sep = "[:-]", convert = TRUE)
    
    # 构建 GRanges (hg19)
    anchor1_hg19 <- GRanges(seqnames = loops_bed$chr1, ranges = IRanges(start = loops_bed$start1, end = loops_bed$end1))
    anchor2_hg19 <- GRanges(seqnames = loops_bed$chr2, ranges = IRanges(start = loops_bed$start2, end = loops_bed$end2))
    
    # 增加 ID 以便转换后配对
    anchor1_hg19$pair_id <- 1:length(anchor1_hg19)
    anchor2_hg19$pair_id <- 1:length(anchor2_hg19)
    
    # 转换为 hg38
    a1_38 <- unlist(liftOver(anchor1_hg19, chain))
    a2_38 <- unlist(liftOver(anchor2_hg19, chain))
    
    # 准备基因启动子 (hg38)
    genes_gr <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
    promoters_gr <- promoters(genes_gr, upstream=2000, downstream=0)
    
    syms <- AnnotationDbi::select(org.Hs.eg.db, keys=names(promoters_gr), columns="SYMBOL", keytype="ENTREZID")
    syms <- syms[!is.na(syms$SYMBOL),]
    promoters_gr <- promoters_gr[names(promoters_gr) %in% syms$ENTREZID]
    names(promoters_gr) <- syms$SYMBOL[match(names(promoters_gr), syms$ENTREZID)]
    
    # 查找重叠 
    ov_eRNA_a1 <- findOverlaps(eRNA_gr, a1_38)
    ov_eRNA_a2 <- findOverlaps(eRNA_gr, a2_38)
    
    valid_loops_idx <- unique(c(a1_38$pair_id[subjectHits(ov_eRNA_a1)], 
                                a2_38$pair_id[subjectHits(ov_eRNA_a2)]))
    
    if(length(valid_loops_idx) > 0) {
      cand_a1 <- a1_38[a1_38$pair_id %in% valid_loops_idx]
      cand_a2 <- a2_38[a2_38$pair_id %in% valid_loops_idx]
      
      ov_prom_a1 <- findOverlaps(promoters_gr, cand_a1)
      ov_prom_a2 <- findOverlaps(promoters_gr, cand_a2)
      
      res_a1 <- data.frame(Gene = names(promoters_gr)[queryHits(ov_prom_a1)], pair_id = cand_a1$pair_id[subjectHits(ov_prom_a1)])
      res_a2 <- data.frame(Gene = names(promoters_gr)[queryHits(ov_prom_a2)], pair_id = cand_a2$pair_id[subjectHits(ov_prom_a2)])
      
      # 映射 eRNA
      eRNA_map <- rbind(
        data.frame(eRNA=names(eRNA_gr)[queryHits(ov_eRNA_a1)], pair_id=a1_38$pair_id[subjectHits(ov_eRNA_a1)]),
        data.frame(eRNA=names(eRNA_gr)[queryHits(ov_eRNA_a2)], pair_id=a2_38$pair_id[subjectHits(ov_eRNA_a2)])
      )
      
      gene_map_res <- rbind(res_a1, res_a2)
      hic_validated_pairs <- inner_join(eRNA_map, gene_map_res, by="pair_id") %>% dplyr::select(eRNA, Gene) %>% distinct()
      
      cat("坐标统一后，Hi-C 验证找到", nrow(hic_validated_pairs), "对互作。\n")
    }
  }
}


# ==============================================================================
#                             Part 6: 绘图
# ==============================================================================

cnv_hits <- list()      
plot_data_list <- list() 
cnv_colors <- c("-2" = "#457B9D", "-1" = "#7AA6C2", 
                "0"  = "#E0E0E0", 
                "1"  = "#F28E7E", "2"  = "#E64B35")

if(file.exists(xena_cnv_file)) {
  
  if(!exists("cnv_mat")) {
    cnv_df <- fread(xena_cnv_file, data.table = FALSE)
    rownames(cnv_df) <- cnv_df[, 1]
    cnv_mat <- as.matrix(cnv_df[, -1])
    colnames(cnv_mat) <- substr(colnames(cnv_mat), 1, 12) 
    cnv_mat <- cnv_mat[, !duplicated(colnames(cnv_mat))] 
  }
  
  if(!exists("proxy_map")) {
    target_gr <- eRNA_gr[names(eRNA_gr) %in% valid_eRNAs]
    nearest_idx <- nearest(target_gr, genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
    nearest_gene_ids <- names(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))[nearest_idx]
    
    gene_map_ids <- AnnotationDbi::select(org.Hs.eg.db, keys=nearest_gene_ids, columns="SYMBOL", keytype="ENTREZID")
    gene_map_ids <- gene_map_ids[!duplicated(gene_map_ids$ENTREZID), ]
    
    proxy_map <- data.frame(eRNA_ID = names(target_gr), Proxy_Entrez = nearest_gene_ids, stringsAsFactors = F)
    proxy_map <- merge(proxy_map, gene_map_ids, by.x="Proxy_Entrez", by.y="ENTREZID", all.x=TRUE)
    proxy_map <- proxy_map[!is.na(proxy_map$SYMBOL), ]
  }
  
  expr_t <- t(vst_matrix)
  rownames(expr_t) <- substr(rownames(expr_t), 1, 12)
  
  stats_list <- list()
  
  for(i in 1:nrow(proxy_map)) {
    curr_erna <- proxy_map$eRNA_ID[i]
    curr_proxy <- proxy_map$SYMBOL[i]
    
    if(!curr_erna %in% valid_eRNAs) next
    if(!curr_erna %in% colnames(expr_t)) next
    if(!curr_proxy %in% rownames(cnv_mat)) next
    
    # 准备数据
    df_plot <- data.frame(sample = rownames(expr_t), expr = expr_t[, curr_erna])
    cnv_vals <- cnv_mat[curr_proxy, ]
    df_cnv <- data.frame(sample = names(cnv_vals), cnv = as.numeric(cnv_vals))
    merged <- merge(df_plot, df_cnv, by="sample")

    merged$log_expr <- log2(merged$expr + 1)
    
    if(length(unique(merged$cnv)) > 1 && min(table(merged$cnv)) >= 3) {
      test_res <- kruskal.test(expr ~ as.factor(cnv), data = merged)
      p_val <- test_res$p.value
      
      if(p_val < 0.05) {
        stats_list[[length(stats_list)+1]] <- data.frame(eRNA = curr_erna, Proxy = curr_proxy, P_val = p_val, stringsAsFactors = F)
        plot_data_list[[curr_erna]] <- merged
      }
    }
  }
}  
  # 绘图
  if(length(stats_list) > 0) {
    stats_df <- do.call(rbind, stats_list)
    stats_df <- stats_df[order(stats_df$P_val), ] 
    
    final_plots <- list()
    for(g in stats_df$eRNA) {
      curr_data <- plot_data_list[[g]]
      curr_proxy <- stats_df$Proxy[stats_df$eRNA == g]
      curr_p <- stats_df$P_val[stats_df$eRNA == g]
      
      p_label <- ifelse(curr_p < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", curr_p)))
      
      p <- ggplot(curr_data, aes(x=as.factor(cnv), y=log_expr, fill=as.factor(cnv))) +
        geom_jitter(width=0.3, size=0.3, color="#B0B0B0", alpha=0.8) + 
        geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6, size = 0.3) + 
        
        scale_fill_manual(values = cnv_colors) + 
        scale_x_discrete(labels = c("-2"="Del", "-1"="-1", "0"="Neu", "1"="+1", "2"="Amp")) + 
        
        labs(
          title = g,
          subtitle = paste0("Proxy: ", curr_proxy, " | ", p_label),
          x = NULL, 
          y = "Log2 (Expression + 1)"
        ) +
        
        theme_classic() + 
        theme(
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5, color = "black"),
          plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5), 
          axis.text.x = element_text(size = 9, color = "black"),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.title.y = element_text(size = 9, color = "black"),
          axis.line = element_line(size = 0.4),
          legend.position = "none" 
        )
      final_plots[[length(final_plots)+1]] <- p
    }
    
    # 布局计算
    n_plots <- length(final_plots)
    if(n_plots <= 3) { ncol_set <- n_plots } else if (n_plots <= 6) { ncol_set <- 3 } else { ncol_set <- 4 }
    
    grid.arrange(grobs = final_plots, ncol = ncol_set)
    
  } else {
    cat("   未检测到显著结果。\n")
  }

library(svglite)
n_col_save <- 4 
combined_plot <- arrangeGrob(grobs = final_plots, ncol = n_col_save)
ggsave(
  filename = "Results/Fig_CNV_Boxplots.svg", 
  plot = combined_plot, 
  device = "svg", 
  width = 14, 
  height = 7,   
  dpi = 300     
)
cat("文件已保存为: Results/Fig_CNV_Boxplots.svg\n")


# ==============================================================================
#                         Step 7: 网络图
# ==============================================================================
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(readr)

eRNA_mat <- t(exp_matrix[valid_eRNAs, , drop=FALSE])
TF_mat <- t(exp_matrix[valid_tfs, , drop=FALSE])

hic_genes <- if(exists("hic_validated_pairs")) hic_validated_pairs$Gene else c()
final_gene_pool <- intersect(hic_genes, rownames(exp_matrix)) 
if(length(final_gene_pool) == 0) final_gene_pool <- rownames(exp_matrix)[1:10] 

Gene_mat <- t(exp_matrix[final_gene_pool, , drop=FALSE])

# 初始化边列表
edges_list <- list()

# Hi-C (Blue)
if(exists("hic_validated_pairs") && nrow(hic_validated_pairs) > 0) {
  edges_list[["hic"]] <- data.frame(
    from = hic_validated_pairs$eRNA, to = hic_validated_pairs$Gene,
    weight = 1, type = "eRNA-Gene", mechanism = "HiC-Validated", stringsAsFactors = F
  )
}

# Mutation (Red) 
if(exists("motif_breaks") && !is.null(motif_breaks) && nrow(motif_breaks) > 0) {
  edges_list[["mut"]] <- data.frame(
    from = motif_breaks$geneSymbol, to = motif_breaks$eRNA_ID,
    weight = 1, type = "TF-eRNA", mechanism = "Mut/CNV-Driven", stringsAsFactors = F
  )
}

# Co-Exp (Grey) 
cor_tf <- cor(TF_mat, eRNA_mat, method = "spearman")
cor_tf[is.na(cor_tf)] <- 0
edges_tf_bg <- data.frame()
for(e in colnames(cor_tf)) {
  # 取 Top 10 相关 TF
  sigs <- sort(abs(cor_tf[,e]), decreasing = T)
  top_tfs <- names(head(sigs, 10)) 
  if(length(top_tfs) > 0) {
    edges_tf_bg <- rbind(edges_tf_bg, data.frame(
      from = top_tfs, to = e, weight = cor_tf[top_tfs, e], 
      type = "TF-eRNA", mechanism = "Co-Exp", stringsAsFactors = F
    ))
  }
}
edges_list[["tf_bg"]] <- edges_tf_bg

# 合并
all_edges <- do.call(rbind, edges_list)
all_edges <- all_edges %>%
  group_by(from, to) %>%
  arrange(factor(mechanism, levels = c("HiC-Validated", "Mut/CNV-Driven", "Co-Exp"))) %>%
  slice(1) %>%
  ungroup() %>%
  as.data.frame()

if(exists("cnv_hits") && length(cnv_hits) > 0) {
  cnv_drivers <- names(cnv_hits)
  for(dr_erna in cnv_drivers) {
    idx <- which((all_edges$from == dr_erna | all_edges$to == dr_erna) & 
                   all_edges$mechanism == "Co-Exp")
    
    if(length(idx) > 0) {
      sub_weights <- abs(all_edges$weight[idx])
      top_k <- 3
      elite_idx <- idx[order(sub_weights, decreasing = TRUE)][1:min(top_k, length(idx))]
      all_edges$mechanism[elite_idx] <- "Mut/CNV-Driven"
    }
  }
}

cat("网络统计:\n")
print(table(all_edges$mechanism))

all_node_names <- unique(c(all_edges$from, all_edges$to))
nodes_df <- data.frame(name = all_node_names, stringsAsFactors = FALSE) %>%
  mutate(type = case_when(
    name %in% valid_eRNAs ~ "eRNA",
    name %in% valid_tfs ~ "TF",
    TRUE ~ "Gene"
  ))

cat("节点类型统计:\n")
print(table(nodes_df$type))

cytoscape_edges <- all_edges %>%
  dplyr::rename(
    Source = from,       
    Target = to,         
    Interaction = mechanism 
  ) %>%
  dplyr::mutate(
    edge_width = case_when(
      Interaction == "HiC-Validated" ~ 3.0,
      Interaction == "Mut/CNV-Driven" ~ 3.0,
      Interaction == "Co-Exp" ~ 1.0,
      TRUE ~ 1.0
    ),
    edge_color = case_when(
      Interaction == "HiC-Validated" ~ "#377EB8",    
      Interaction == "Mut/CNV-Driven" ~ "#E41A1C",
      Interaction == "Co-Exp" ~ "#999999",          
      TRUE ~ "#000000"
    )
  )

cytoscape_nodes <- nodes_df %>%
  dplyr::rename(ID = name, NodeType = type) %>%
  dplyr::mutate(
    node_shape = case_when(
      NodeType == "TF" ~ "DIAMOND",      
      NodeType == "eRNA" ~ "HEXAGON",    
      NodeType == "Gene" ~ "ELLIPSE",    
      TRUE ~ "ELLIPSE"
    ),
    node_fill_color = case_when(
      NodeType == "TF" ~ "#4DAF4A",     
      NodeType == "eRNA" ~ "#FF7F00",   
      NodeType == "Gene" ~ "#984EA3",   
      TRUE ~ "#CCCCCC"
    ),
    node_size = case_when(
      NodeType == "TF" ~ 50,
      NodeType == "eRNA" ~ 60,
      NodeType == "Gene" ~ 30,
      TRUE ~ 30
    )
  )
nodes_df$type[nodes_df$name == "EGR3"] <- "TF"
cytoscape_nodes <- nodes_df %>%
  dplyr::rename(ID = name, NodeType = type) %>%
  dplyr::mutate(
    node_shape = case_when(
      NodeType == "TF" ~ "DIAMOND",      
      NodeType == "eRNA" ~ "HEXAGON",    
      NodeType == "Gene" ~ "ELLIPSE",    
      TRUE ~ "ELLIPSE"
    ),
    node_fill_color = case_when(
      NodeType == "TF" ~ "#4DAF4A",     
      NodeType == "eRNA" ~ "#FF7F00",   
      NodeType == "Gene" ~ "#984EA3",   
      TRUE ~ "#CCCCCC"
    ),
    node_size = case_when(
      NodeType == "TF" ~ 50,
      NodeType == "eRNA" ~ 60,
      NodeType == "Gene" ~ 30,
      TRUE ~ 30
    )
  )

# 导出
out_dir <- getwd()
write_csv(cytoscape_edges, file.path(out_dir, "Results/Network_Edges_for_Cytoscape.csv"))
write_csv(cytoscape_nodes, file.path(out_dir, "Results/Network_Nodes_for_Cytoscape.csv"))


