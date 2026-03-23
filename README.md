# TCGA-BRCA eRNA Prognostic Signature Project

## 📌 Project Overview
This project identifies and validates a prognostic enhancer RNA (eRNA) signature for Breast Cancer (BRCA) using multi-omics data. The analysis pipeline seamlessly integrates differential expression analysis, LASSO-Cox modeling, GSEA pathway enrichment, multi-omics regulatory network construction (CNV, Mutation, Hi-C), drug sensitivity prediction, and rigorous quality control of normal-like references for further early-discrimination model.

## 📂 Repository Structure
The project is organized to reflect the chronological workflow of the associated manuscript. 
**Note:** `Data_Source` is excluded from the repository due to size limits.

```text
Prognosis_R/
├── Prognosis_R.Rproj         <-- Double-click this to open the project
├── README.md                 <-- Project documentation
├── Data_Source/              <-- Contains all raw input files (Need to be downloaded)
├── R_Scripts/                <-- Sequentially numbered analysis scripts
│   ├── 1_Differentially_Expressed_eRNA_Recognition.R
│   ├── 2_Risk_Score_Stratification&KM_Curve.R
│   ├── 3_Independent_Prognostic_Analysis_Cox.R
│   ├── 4_GSEA_for_Risk_Groups.R
│   ├── 5_Comprehensive_Regulatory_Network_eRNA&TF&Gene.R
│   ├── 6_Candidate_Drug&Immune_Analysis.R
│   ├── 7_PCA_mRNA_with_eRNA_Tumor&Normal.R
│   └── 8_Quality_Control_for_Normal_Like_References.R
└── Results/                  <-- Automatically generated Figures and Tables
    └── Drug_Analysis/        <-- L1000FWD web results and drug prediction plots
```

💾 Data Preparation (Crucial)
⚠️ Action Required: Raw data is NOT included in this repository. To reproduce the analysis, please create a folder named Data_Source in the root directory and download the required datasets.

File Category	Required Filename / Data	Source / Database

eRNA Expr	TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv	TCeA (The Cancer eRNA Atlas)
mRNA Expr	TCGA-BRCA.star_fpkm.tsv.gz	UCSC Xena (TCGA-BRCA)
Clinical	TCGA-BRCA.clinical.tsv & clinical_info.tsv	UCSC Xena
Multi-omics	Mutation, CNV, Methylation, and Hi-C loop data	UCSC Xena / ENCODE
Normal Ref	GTEx breast tissue expression profile	GTEx Portal (For Script 8)
Drug Data	L1000_Result.csv (Place in Results/Drug_Analysis/)	L1000FWD Web Tool (Generated via Script 6)
Note: If you are a member of Southeast University, you can access the pre-compiled data via pan.seu.edu.cn. For inquiries, please contact 213230182@seu.edu.cn.

🚀 How to Run the Pipeline
Prerequisite: Always open the project by double-clicking Prognosis_R.Rproj to ensure the working directory is set correctly. Run the scripts in numerical order:

Part I: Prognostic Model Construction & Validation
Script 1: 1_Differentially_Expressed_eRNA_Recognition.R
Function: Identifies DEEs between tumor and adjacent normal tissues.
Script 2: 2_Risk_Score_Stratification&KM_Curve.R
Function: Applies LASSO-Cox regression to build the multi-eRNA signature, stratifies patients, and evaluates via KM curves and time-dependent ROC.
Script 3: 3_Independent_Prognostic_Analysis_Cox.R
Function: Validates the signature's independence from traditional clinical traits (e.g., age, TNM stage) via multivariate Cox regression.
Part II: Biological Mechanism & Regulatory Network
Script 4: 4_GSEA_for_Risk_Groups.R
Function: Performs Gene Set Enrichment Analysis (GSEA) to uncover pathways driving the high-risk phenotype.
Script 5: 5_Comprehensive_Regulatory_Network_eRNA&TF&Gene.R
Function: Integrates WGS mutations, CNV, and Hi-C validated chromatin loops to map upstream TF hubs and downstream targets.
Part III: Clinical Translation & Therapeutics
Script 6: 6_Candidate_Drug&Immune_Analysis.R
Function: Prepares top DEGs for the L1000FWD platform to identify small-molecule compounds that could reverse the high-risk transcriptional state, and visualizes drug sensitivity predictions.
Part IV: Early-Stage Discrimination & Quality Control
Script 7: 7_PCA_mRNA_with_eRNA_Tumor&Normal.R
Function: Explores global transcriptomic variance between tumor and normal states.
Script 8: 8_Quality_Control_for_Normal_Like_References.R
Function: Calculates Spearman correlation between TCGA adjacent normal tissues and GTEx healthy tissues to rigorously validate the non-tumor baseline.
🛠 Dependencies
R Version: 4.x
Key Packages: survival, survminer, glmnet, timeROC, clusterProfiler, ggplot2, patchwork, ggpubr, dplyr, data.table.
