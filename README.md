# TCGA-BRCA eRNA Prognostic Signature Project

## 📌 Project Overview
This project identifies and validates a prognostic enhancer RNA (eRNA) signature for Breast Cancer (BRCA) using multi-omics data. The analysis pipeline seamlessly integrates differential expression analysis, LASSO-Cox modeling, GSEA pathway enrichment, multi-omics regulatory network construction (CNV, Mutation, Hi-C), drug sensitivity prediction, and rigorous quality control of normal-like references for further early-discrimination model.

## 🌐 Interactive Web Application 
To make our models accessible to clinicians and researchers without programming expertise, we have deployed a user-friendly online tool: eRNACare. You can upload your own eRNA expression profiles to get real-time diagnostic predictions without running any code.

### **Access the Web App here:** https://limseu.shinyapps.io/ernacare/

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

## 💾 Data Preparation & Setup

Some files are already included in this repository (see `Data_Source/` folder).
Large files must be downloaded manually and placed in the same directory.

> **Note:** If you are a member of Southeast University, you can access
> the pre-compiled data via `pan.seu.edu.cn`.
> For inquiries, please contact 213230182@seu.edu.cn.

---

### ✅ Files Already Included in This Repository

| **File** | **Notes** |
| :--- | :--- |
| `TCGA-BRCA.clinical.tsv` | GDC TCGA BRCA phenotype |
| `TCGA-BRCA.survival.tsv` | GDC TCGA BRCA survival data |
| `TCGA.BRCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes` | TCGA BRCA copy number |
| `TcgaTargetGTEX_phenotype.txt.gz` | TCGA TARGET GTEx phenotype |
| `pcawg_donor_clinical_August2016_v9` | PCAWG donor clinical data (original filename) |
| `gencode.v36.annotation.gtf.gene.probemap` | Gencode v36 gene annotation |
| `hg19ToHg38.fixed_universal.chain` | Auto-generated from `hg19ToHg38.over.chain.gz` by the pipeline — do not download separately |
| `hg19ToHg38.over.chain.gz` | UCSC liftOver reference |
| `loop_info.csv` | Hi-C loop data (GEO: GSE157381) |
| `probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy` | Methylation 450k probe annotation |
| Drug Analysis | `L1000_Result.csv` → place in `Results/Drug_Analysis/` | **L1000FWD** | Submit your gene signature at https://maayanlab.cloud/L1000FWD/ and export results |

---

### ⬇️ Files That Must Be Downloaded Manually

Please create a folder named `Data_Source/` in the root directory and
download the following files into it.
**Filenames must be exact.**

| **Category** | **Filename** | **Source** | **How to Get** |
| :--- | :--- | :--- | :--- |
| eRNA Expression | `TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv` | **TCeA** | Navigate to: [TCeA eRNA Quantification](https://bioinformatics.mdanderson.org/Supplements/Super_Enhancer/TCEA_website/parts/3_eRNA_quantification.html) → TCGA-BRCA → eRNAs in ENCODE super-enhancers (~300k) |
| mRNA Expression | `TCGA-BRCA.star_fpkm.tsv.gz` | **UCSC Xena** | Cohort: GDC TCGA Breast Cancer (BRCA) → Gene expression RNAseq |
| Somatic Mutation | `October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS` | **UCSC Xena** | Cohort: PCAWG (donor centric) → Simple somatic mutation (SNVs and indels), non-US donors |
| CNV | `TCGA.BRCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz` | **UCSC Xena** | Cohort: TCGA Breast Cancer (BRCA) → Copy Number (Gistic2) |
| Methylation | `TCGA.BRCA.sampleMap_HumanMethylation450.gz` | **UCSC Xena** | Cohort: TCGA Breast Cancer (BRCA) → DNA Methylation 450k |
| Harmonized Expression | `TcgaTargetGtex_rsem_gene_tpm.tsv` | **UCSC Xena** | Cohort: TCGA TARGET GTEx → Gene expression RNAseq (TOIL) |
| Reference | `hg19ToHg38.over.chain.gz` | **UCSC Genome Browser** | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/ |

> 💡 **UCSC Xena main portal:** https://xenabrowser.net/datapages/
> All Xena files above can be located by searching the filename or
> dataset name directly in the portal.

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
