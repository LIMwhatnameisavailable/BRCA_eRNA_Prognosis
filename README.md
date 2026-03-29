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

## 💾 Data Preparation (Crucial)

**⚠️ Action Required:** Raw data is **NOT** included in this repository. To reproduce the analysis, please create a folder named `Data_Source` in the root directory and download the required datasets.

| File Category | Required Filename (Must be Exact) | Source / Database |
| :--- | :--- | :--- |
| **eRNA Expr** | `TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv` | **TCeA** (The Cancer eRNA Atlas) |
| **mRNA Expr** | `TCGA-BRCA.star_fpkm.tsv.gz` | **UCSC Xena** (TCGA-BRCA) |
| **Clinical** | `TCGA-BRCA.clinical.tsv` | **UCSC Xena** (Phenotype data) |
| **Clinical** | `clinical_info.tsv` | **UCSC Xena** (Survival data) |
| **Mutation** | `PCAWG_WGS_mutations.tsv.gz` | **UCSC Xena** |
| **CNV** | `TCGA.BRCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz` | **UCSC Xena** (Gistic2) |
| **Methylation** | `TCGA.BRCA.sampleMap_HumanMethylation450.gz` | **UCSC Xena** (450k Array) |
| **Annotation** | `probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy` | **UCSC Xena** (Platform Map) |
| **Hi-C** | `loop_info.csv` | **GEO(GSE157381)** |
| **Reference** | `hg19ToHg38.over.chain.gz` | **UCSC Genome Browser** |
| **Reference** | `gencode.v36.annotation.gtf.gene.probemap` | **UCSC Xena** (TCGA TARGET GTEx)|
| **Normal Ref** | `TcgaTargetGtex_rsem_gene_tpm.tsv` | **UCSC Xena** |
| **Drug Data** | `L1000_Result.csv` (Place in `Results/Drug_Analysis/`) | **L1000FWD Web Tool** |

> **Note:** If you are a member of Southeast University, you can access the pre-compiled data via `pan.seu.edu.cn`. For inquiries, please contact 213230182@seu.edu.cn.

## 🚀 How to Run the Pipeline

**Prerequisite:** Always open the project by double-clicking **`Prognosis_R.Rproj`** to ensure the working directory is set correctly. Run the scripts in numerical order:

### Part I: Prognostic Model Construction & Validation
* **Script 1:** `1_Differentially_Expressed_eRNA_Recognition.R`
  * **Function:** Identifies DEEs between tumor and adjacent normal tissues.
* **Script 2:** `2_Risk_Score_Stratification&KM_Curve.R`
  * **Function:** Applies LASSO-Cox regression to build the multi-eRNA signature, stratifies patients, and evaluates via KM curves and time-dependent ROC.
* **Script 3:** `3_Independent_Prognostic_Analysis_Cox.R`
  * **Function:** Validates the signature's independence from traditional clinical traits (e.g., age, TNM stage) via multivariate Cox regression.

### Part II: Biological Mechanism & Regulatory Network
* **Script 4:** `4_GSEA_for_Risk_Groups.R`
  * **Function:** Performs Gene Set Enrichment Analysis (GSEA) to uncover pathways driving the high-risk phenotype.
* **Script 5:** `5_Comprehensive_Regulatory_Network_eRNA&TF&Gene.R`
  * **Function:** Integrates WGS mutations, CNV, and Hi-C validated chromatin loops to map upstream TF hubs and downstream targets.

### Part III: Clinical Translation & Therapeutics
* **Script 6:** `6_Candidate_Drug&Immune_Analysis.R`
  * **Function:** Conducts a comprehensive pharmacogenomic and immunological analysis to identify therapeutic vulnerabilities in high-risk patients. Key executions include:
    * **Signature Reversal:** Extracts DEGs for web-based screening and visualizes compounds capable of reversing the high-risk transcriptional state.
    * **Therapy Resistance:** Evaluates GDSC-predicted IC50 values to demonstrate high-risk patient resistance to standard chemotherapies.
    * **Resistance Mechanisms:** Performs Spearman correlation analysis between drug resistance (IC50) and 19 key genes.
    * **Immune & Targeted Therapy:** Compares immune checkpoint expression across risk groups and identifies targeted drugs.
    * **Therapeutics Summary:** Automatically generates a consolidated table of top candidate drugs.

### Part IV: Early-Stage Discrimination & Quality Control
* **Script 7:** `7_PCA_mRNA_with_eRNA_Tumor&Normal.R`
  * **Function:** Explores global transcriptomic variance between tumor and normal states.
* **Script 8:** `8_Quality_Control_for_Normal_Like_References.R`
  * **Function:** Calculates Spearman correlation between TCGA adjacent normal tissues and GTEx healthy tissues to rigorously validate the non-tumor baseline.

## 🛠 Dependencies
* **R Version:** 4.x
* **Key Packages:** `survival`, `survminer`, `glmnet`, `timeROC`, `clusterProfiler`, `ggplot2`, `patchwork`, `ggpubr`, `dplyr`, `data.table`.
