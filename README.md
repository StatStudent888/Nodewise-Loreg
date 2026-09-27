<a id="readme-top"></a>

<div align="center">
  <h1>Nodewise Loreg Reproducibility Materials</h1>
  <p>
    This file accompanies the reproducibility materials for the paper:<br />
    <strong>Nodewise Loreg: Nodewise L0 Regression for High-dimensional Sparse Precision Matrix Estimation</strong>
  </p>
</div>

---

<details open>
<summary><strong>Table of Contents</strong></summary>

- [Description of the paper](#description-of-the-paper)
- [Repository structure](#repository-structure)
  - [`software/`: method implementation and application example](#software-method-implementation-and-application-example)
  - [`reproduce/`: scripts for reproducing the manuscript results](#reproduce-scripts-for-reproducing-the-manuscript-results)
- [Figure and table reproduction map](#figure-and-table-reproduction-map)
- [Requirements](#requirements)
  - [Main software](#main-software)
  - [R packages](#r-packages)
  - [Python packages](#python-packages)
- [Data access](#data-access)
  - [Raw data and preprocessing code](#raw-data-and-preprocessing-code)
  - [Processed log-transformed data](#processed-log-transformed-data)
- [Approximate runtime](#approximate-runtime)
  - [Quick example](#quick-example)
  - [Main SDAR-based Nodewise Loreg simulations](#main-sdar-based-nodewise-loreg-simulations)
  - [Real-data analysis](#real-data-analysis)
  - [MIO-based and L0BnB-based comparison](#mio-based-and-l0bnb-based-comparison)

</details>

---

## Description of the paper

This paper proposes **Nodewise Loreg**, a nodewise L0-penalized regression method for estimating high-dimensional sparse precision matrices, with SDAR as the default solver. It establishes convergence-rate, support-recovery, and asymptotic-normality results under the stated assumptions. The main simulation studies compare SDAR-based Nodewise Loreg with Nodewise Lasso, GLasso, and CLIME. The supplementary simulations also examine positive definiteness and SDAR convergence, and evaluate alternative MIO-based and L0BnB-based implementations of Nodewise Loreg. The MDA133 breast cancer gene-expression analyses compare SDAR-based Nodewise Loreg, Nodewise Lasso, GLasso, and CLIME using p = 300 and p = 1400 gene probesets.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Repository structure

To make the reproducibility materials easier to navigate, we reorganized the repository into two main folders: `software/` and `reproduce/`.

### `software/`: method implementation and application example

The `software/` folder contains the code implementing the proposed SDAR-based Nodewise Loreg method and a minimal example showing how to apply the method.

- `software/Example/`: contains the main R implementation and a small example.
  - `software/Example/Nodewise_Loreg.R`: main implementation of the SDAR-based Nodewise Loreg method.
  - `software/Example/Test.R`: a small example showing how to call the Nodewise Loreg function.
- `software/use_for_all/`: contains shared helper functions used by the method and the reproduction scripts, including the SDAR C++ implementation, matrix-generation functions, precision-matrix construction functions, FDR/AdaptZ-related functions, and other utility functions.

Users who want to apply Nodewise Loreg to their own data should start from `software/Example/Test.R`.

### `reproduce/`: scripts for reproducing the manuscript results

The `reproduce/` folder contains the scripts and data files used to reproduce the simulation studies, real-data analysis, tables, and figures reported in the manuscript and Supplementary Materials.

- `reproduce/Gaussian-design/`: scripts for reproducing the Gaussian-design simulation studies.
- `reproduce/Sub-Gaussian-design/`: scripts for reproducing the sub-Gaussian-design simulation studies.
- `reproduce/Real_data/`: scripts and data files for reproducing the MDA133 breast cancer gene-expression analysis.
- `reproduce/seed.txt`: random seeds used in the simulation studies.
- `reproduce/realdata_seed.txt`: random seeds used in the real-data analysis.

In short, `software/` contains the reusable implementation of the proposed method, whereas `reproduce/` contains the scripts used to reproduce the numerical results in the paper.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Figure and table reproduction map

For simulations, run the replication-level (`One-replication`) scripts before the corresponding aggregation (`Combine`) scripts. The numerical simulation tables use 100 replications; the Z-score histograms use 400 replications. Select the sample size, dimension, graph structure, and distribution corresponding to the required output.

<details open>
<summary><strong>Detailed figure and table reproduction map</strong></summary>

| Output in paper | Purpose | Scripts to run | Order / notes |
|---|---|---|---|
| Supplementary Tables S.2, S.4, S.13, and S.15 | Matrix norm losses (Tables S.2 and S.13) and support-recovery metrics (Tables S.4 and S.15) under Gaussian settings. | `reproduce/Gaussian-design/CLIME/One-replication/CLIME.R`<br>`reproduce/Gaussian-design/CLIME/Combine/CLIME_combine.R`<br>`reproduce/Gaussian-design/GLasso/One-replication/GLasso_Gaussian.R`<br>`reproduce/Gaussian-design/GLasso/Combine/GLasso_Gaussian_combine.R`<br>`reproduce/Gaussian-design/L0/One-replication/L0_desparsified_unbias_Gaussian_one-two-stage_lfdr.R`<br>`reproduce/Gaussian-design/L0/Combine/l0_desparsified_unbias_Gaussian_one-two-stage_combine.R`<br>`reproduce/Gaussian-design/L1/One-replication/L1_desparsified_Gaussian_one-two-stage_lfdr.R`<br>`reproduce/Gaussian-design/L1/Combine/l1_desparsified_Gaussian_one-two-stage_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 100 replications and n = 200 or 400. Tables S.2/S.4 use p = 200 or 400; Tables S.13/S.15 use p = 1000 or 4000. CLIME is included only for p = 200 or 400; do not run CLIME to reproduce Tables S.13/S.15. |
| Supplementary Tables S.3, S.5, S.14, and S.16 | Matrix norm losses (Tables S.3 and S.14) and support-recovery metrics (Tables S.5 and S.16) under sub-Gaussian settings. | `reproduce/Gaussian-design/CLIME/One-replication/CLIME.R`<br>`reproduce/Gaussian-design/CLIME/Combine/CLIME_combine.R`<br>`reproduce/Sub-Gaussian-design/GLasso/One-replication/GLasso_subGaussian.R`<br>`reproduce/Sub-Gaussian-design/GLasso/Combine/GLasso_subGaussian_combine.R`<br>`reproduce/Sub-Gaussian-design/L0/One-replication/L0_desparsified_unbias_subGaussian_one-two-stage_lfdr.R`<br>`reproduce/Sub-Gaussian-design/L0/Combine/l0_desparsified_unbias_subGaussian_one-two-stage_combine.R`<br>`reproduce/Sub-Gaussian-design/L1/One-replication/L1_desparsified_subGaussian_one-two-stage_lfdr.R`<br>`reproduce/Sub-Gaussian-design/L1/Combine/l1_desparsified_subGaussian_one-two-stage_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 100 replications and n = 200 or 400. Tables S.3/S.5 use p = 200 or 400; Tables S.14/S.16 use p = 1000 or 4000. CLIME is included only for p = 200 or 400. |
| Supplementary Tables S.6--S.8 and S.17--S.19 | Asymptotic normality metrics under Gaussian design. | `reproduce/Gaussian-design/L0/One-replication/L0_desparsified_unbias_Gaussian_normality.R`<br>`reproduce/Gaussian-design/L0/One-replication/L0_true_length_for_otherp.R`<br>`reproduce/Gaussian-design/L0/Combine/l0_Gaussian_normality_combine.R`<br>`reproduce/Gaussian-design/L0/One-replication/L0_true_length_for_p4000.R`<br>`reproduce/Gaussian-design/L0/Combine/l0_true_length_combine_for_p4000.R`<br>`reproduce/Gaussian-design/L1/One-replication/L1_desparsified_Gaussian_normality.R`<br>`reproduce/Gaussian-design/L1/Combine/l1_Gaussian_normality_combine.R` | Run replication-level and required TrueLength calculations before aggregation. Use 100 replications, n = 200 or 400, and 95% confidence intervals. Tables S.6--S.8 use p = 200 or 400; Tables S.17--S.19 use p = 1000 or 4000. |
| Supplementary Tables S.9--S.11 and S.20--S.22 | Asymptotic normality metrics under sub-Gaussian design. | `reproduce/Sub-Gaussian-design/L0/One-replication/L0_desparsified_unbias_subGaussian_normality.R`<br>`reproduce/Gaussian-design/L0/One-replication/L0_true_length_for_otherp.R`<br>`reproduce/Sub-Gaussian-design/L0/Combine/l0_subGaussian_normality_combine.R`<br>`reproduce/Gaussian-design/L0/One-replication/L0_true_length_for_p4000.R`<br>`reproduce/Gaussian-design/L0/Combine/l0_true_length_combine_for_p4000.R`<br>`reproduce/Sub-Gaussian-design/L1/One-replication/L1_desparsified_subGaussian_normality.R`<br>`reproduce/Sub-Gaussian-design/L1/Combine/l1_subGaussian_normality_combine.R` | Run replication-level and required TrueLength calculations before aggregation. Use 100 replications, n = 200 or 400, and 95% confidence intervals. Tables S.9--S.11 use p = 200 or 400; Tables S.20--S.22 use p = 1000 or 4000. |
| Supplementary Tables S.12 and S.23 | Runtime under Gaussian and sub-Gaussian designs. | `reproduce/Gaussian-design/CLIME/One-replication/CLIME.R`<br>`reproduce/Gaussian-design/CLIME/Combine/CLIME_combine.R`<br>`reproduce/Gaussian-design/GLasso/One-replication/GLasso_Gaussian.R`<br>`reproduce/Gaussian-design/GLasso/Combine/GLasso_Gaussian_combine.R`<br>`reproduce/Sub-Gaussian-design/GLasso/One-replication/GLasso_subGaussian.R`<br>`reproduce/Sub-Gaussian-design/GLasso/Combine/GLasso_subGaussian_combine.R`<br>`reproduce/Gaussian-design/L0/One-replication/L0_time.R`<br>`reproduce/Gaussian-design/L0/Combine/l0_time_combine.R`<br>`reproduce/Gaussian-design/L1/One-replication/L1_time.R`<br>`reproduce/Gaussian-design/L1/Combine/l1_time_combine.R` | Run replication-level timing scripts before aggregation. Table S.12 uses p = 200 or 400; Table S.23 uses p = 1000 or 4000 and excludes CLIME. |
| Supplementary Table S.24 | Proportions of positive definite estimates for the two SDAR-based Nodewise Loreg estimators under Gaussian settings. | `reproduce/Gaussian-design/L0/One-replication/L0_desparsified_unbias_Gaussian_one-two-stage_lfdr.R`<br>`reproduce/Gaussian-design/L0/Combine/l0_desparsified_unbias_Gaussian_one-two-stage_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 100 replications; p = 200, 400, 1000, or 4000; n = 200 or 400. |
| Supplementary Table S.25 | Proportions of positive definite estimates for the two SDAR-based Nodewise Loreg estimators under sub-Gaussian settings. | `reproduce/Sub-Gaussian-design/L0/One-replication/L0_desparsified_unbias_subGaussian_one-two-stage_lfdr.R`<br>`reproduce/Sub-Gaussian-design/L0/Combine/l0_desparsified_unbias_subGaussian_one-two-stage_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 100 replications; p = 200, 400, 1000, or 4000; n = 200 or 400. |
| Supplementary Table S.26 | Iteration metrics under Gaussian design. | `reproduce/Gaussian-design/L0/One-replication/L0_Gaussian_NumOfIter.R`<br>`reproduce/Gaussian-design/L0/Combine/l0_Gaussian_NumOfIter_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use p = 200, 400, 1000, or 4000 and n = 200 or 400. |
| Supplementary Table S.27 | Iteration metrics under sub-Gaussian design. | `reproduce/Sub-Gaussian-design/L0/One-replication/L0_subGaussian_NumOfIter.R`<br>`reproduce/Sub-Gaussian-design/L0/Combine/l0_subGaussian_NumOfIter_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use p = 200, 400, 1000, or 4000 and n = 200 or 400. |
| Supplementary Table S.28 | Runtime for Nodewise Loreg based on MIO and L0BnB under Gaussian and sub-Gaussian designs. | `reproduce/Gaussian-design/MIO/One-replication/MIO_Gaussian_time.R`<br>`reproduce/Gaussian-design/MIO/Combine/MIO_Gaussian_time_combine.R`<br>`reproduce/Sub-Gaussian-design/MIO/One-replication/MIO_subGaussian_time.R`<br>`reproduce/Sub-Gaussian-design/MIO/Combine/MIO_subGaussian_time_combine.R`<br>`reproduce/Gaussian-design/BnB/One-replication/BnB_Gaussian.R`<br>`reproduce/Gaussian-design/BnB/Combine/BnB_Gaussian_combine.R`<br>`reproduce/Sub-Gaussian-design/BnB/One-replication/BnB_subGaussian.R`<br>`reproduce/Sub-Gaussian-design/BnB/Combine/BnB_subGaussian_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use (n, p) = (200, 200). |
| Supplementary Tables S.29 and S.31 | Matrix norm losses (Table S.29) and support-recovery metrics (Table S.31) for MIO-based and L0BnB-based Nodewise Loreg under Gaussian settings. | `reproduce/Gaussian-design/MIO/One-replication/MIO_desparsified_unbias_Gaussian_one-two-stage_lfdr.R`<br>`reproduce/Gaussian-design/MIO/Combine/MIO_desparsified_unbias_Gaussian_one-two-stage_combine.R`<br>`reproduce/Gaussian-design/BnB/One-replication/BnB_Gaussian.R`<br>`reproduce/Gaussian-design/BnB/Combine/BnB_Gaussian_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 100 replications and p = 200. MIO uses n = 200 or 400; L0BnB is evaluated only at n = 200. |
| Supplementary Tables S.30 and S.32 | Matrix norm losses (Table S.30) and support-recovery metrics (Table S.32) for MIO-based and L0BnB-based Nodewise Loreg under sub-Gaussian settings. | `reproduce/Sub-Gaussian-design/MIO/One-replication/MIO_desparsified_unbias_subGaussian_one-two-stage_lfdr.R`<br>`reproduce/Sub-Gaussian-design/MIO/Combine/MIO_desparsified_unbias_subGaussian_one-two-stage_combine.R`<br>`reproduce/Sub-Gaussian-design/BnB/One-replication/BnB_subGaussian.R`<br>`reproduce/Sub-Gaussian-design/BnB/Combine/BnB_subGaussian_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 100 replications and p = 200. MIO uses n = 200 or 400; L0BnB is evaluated only at n = 200. |
| Supplementary Tables S.33--S.35 | Asymptotic normality metrics for Nodewise Loreg based on MIO under Gaussian and sub-Gaussian designs. | `reproduce/Gaussian-design/MIO/One-replication/MIO_desparsified_unbias_Gaussian_normality.R`<br>`reproduce/Gaussian-design/MIO/Combine/MIO_Gaussian_normality_combine.R`<br>`reproduce/Sub-Gaussian-design/MIO/One-replication/MIO_desparsified_subGaussian_normality.R`<br>`reproduce/Sub-Gaussian-design/MIO/Combine/MIO_subGaussian_normality_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 100 replications, p = 200, n = 200 or 400, and 95% confidence intervals. |
| Supplementary Table S.36 | Means and standard deviations of pCR classification metrics and sparsity using the 300 most significantly differentially expressed gene probesets. | `reproduce/Real_data/Code_for_p300_BH/l0_desparsified_real_data_lfdr.R`<br>`reproduce/Real_data/Code_for_p300_BH/l0_method3_real_data_lfdr.R`<br>`reproduce/Real_data/Code_for_p300_BH/l0_unbias_real_data_lfdr.R` | Use the p = 300 processed datasets. The table includes SDAR-based Nodewise Loreg, Nodewise Lasso, CLIME, and GLasso. |
| Supplementary Table S.37 | Means and standard deviations of pCR classification metrics and sparsity using all 1400 significantly differentially expressed gene probesets. | `reproduce/Real_data/Code_for_p1400_AdaptZ/l0_desparsified_real_data_p1400.R`<br>`reproduce/Real_data/Code_for_p1400_AdaptZ/l0_method3_real_data_p1400.R`<br>`reproduce/Real_data/Code_for_p1400_AdaptZ/l0_unbias_real_data_p1400.R` | Use the p = 1400 processed datasets. The table includes SDAR-based Nodewise Loreg, Nodewise Lasso, CLIME, and GLasso. |
| Supplementary Figures S.1--S.4 and S.9--S.12 | Z-score histograms under Gaussian design, overlaid with the standard Gaussian density. | `reproduce/Gaussian-design/L0/One-replication/L0_desparsified_unbias_Gaussian_normalplot.R`<br>`reproduce/Gaussian-design/L0/Combine/l0_desparsified_unbias_Gaussian_normalplot_combine.R`<br>`reproduce/Gaussian-design/L1/One-replication/L1_desparsified_Gaussian_normalplot.R`<br>`reproduce/Gaussian-design/L1/Combine/l1_desparsified_Gaussian_normalplot_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 400 replications. Figures S.1--S.4 use p = 200 or 400; Figures S.9--S.12 use p = 1000 or 4000. Use n = 200, 400, or 800. |
| Supplementary Figures S.5--S.8 and S.13--S.16 | Z-score histograms under sub-Gaussian design, overlaid with the standard Gaussian density. | `reproduce/Sub-Gaussian-design/L0/One-replication/L0_desparsified_unbias_subGaussian_normalplot.R`<br>`reproduce/Sub-Gaussian-design/L0/Combine/l0_desparsified_unbias_subGaussian_normalplot_combine.R`<br>`reproduce/Sub-Gaussian-design/L1/One-replication/L1_desparsified_subGaussian_normalplot.R`<br>`reproduce/Sub-Gaussian-design/L1/Combine/l1_desparsified_subGaussian_normalplot_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 400 replications. Figures S.5--S.8 use p = 200 or 400; Figures S.13--S.16 use p = 1000 or 4000. Use n = 200, 400, or 800. |
| Supplementary Figures S.17--S.20 | Z-score histograms for undebiased and debiased MIO-based Nodewise Loreg under Gaussian design, overlaid with the standard Gaussian density. | `reproduce/Gaussian-design/MIO/One-replication/MIO_desparsified_unbias_Gaussian_normalplot.R`<br>`reproduce/Gaussian-design/MIO/Combine/MIO_desparsified_unbias_Gaussian_normalplot_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 400 replications. Use p = 200 and n = 200, 400, or 800. |
| Supplementary Figures S.21--S.24 | Z-score histograms for undebiased and debiased MIO-based Nodewise Loreg under sub-Gaussian design, overlaid with the standard Gaussian density. | `reproduce/Sub-Gaussian-design/MIO/One-replication/MIO_desparsified_unbias_subGaussian_normalplot.R`<br>`reproduce/Sub-Gaussian-design/MIO/Combine/MIO_desparsified_unbias_subGaussian_normalplot_combine.R` | Run `One-replication` scripts before `Combine` scripts. Use 400 replications. Use p = 200 and n = 200, 400, or 800. |
| Supplementary Figure S.25 | Venn diagrams comparing detected edges across SDAR-based Nodewise Loreg, Nodewise Lasso, CLIME, and GLasso, using the 300 most significantly differentially expressed gene probesets. | `reproduce/Real_data/Code_for_p300_BH/Gene_analysis_and_VennDiagram.R` | Use the p = 300 processed data. |
| Supplementary Figure S.26 | Estimated partial-correlation networks from SDAR-based Nodewise Loreg, Nodewise Lasso, CLIME, and GLasso, using the 300 most significantly differentially expressed gene probesets. | `reproduce/Real_data/Code_for_p300_BH/l0_unbias_real_data_network-graph.R` | Use the p = 300 processed data. |
| Supplementary Figure S.27 | Venn diagrams comparing detected edges across SDAR-based Nodewise Loreg, Nodewise Lasso, CLIME, and GLasso, using all 1400 significantly differentially expressed gene probesets. | `reproduce/Real_data/Code_for_p1400_AdaptZ/VennDiagram.R` | Use the p = 1400 processed data. |
| Supplementary Figure S.28 | Estimated partial-correlation networks from SDAR-based Nodewise Loreg, Nodewise Lasso, CLIME, and GLasso, using all 1400 significantly differentially expressed gene probesets. | `reproduce/Real_data/Code_for_p1400_AdaptZ/l0_unbias_real_data_network-graph_p1400.R` | Use the p = 1400 processed data. |

</details>

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Requirements

The code was run under the following software environment.

### Main software

- R version: 4.2.2
- Python version: 3.9.20
- Operating system: Windows with Rtools

### R packages

The following R packages are required:

| Package | Version used |
|---|---|
| Rcpp | 1.0.9 |
| RcppArmadillo | 0.11.4.2.1 |
| huge | 1.3.5 |
| MASS | 7.3-58.1 |
| RcppEigen | 0.3.4.0.2 |
| BiocParallel | 1.32.6 |
| parallel | 4.2.2 |
| glasso | 1.11 |
| flare | 1.7.0.1 |
| glmnet | 4.1-6 |
| reticulate | 1.41.0 |
| caret | 6.0-94 |
| ggraph | 2.1.0 |
| igraph | 1.5.0 |
| extras | 0.6.1 |
| VennDiagram | 1.7.3 |
| limma | 3.54.2 |

### Python packages

The following Python packages are required:

| Package | Version used |
|---|---|
| numpy | 1.23.0 |
| gurobipy | 12.0.1 |

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Data access

The real-data analysis corresponds to Section 5 of the paper, "Analysis of MDA133 breast cancer gene expression," and Supplementary Section S.5, "Additional results and analysis for MDA133 data." The analysis using the 300 most significantly differentially expressed gene probesets is reported in Table S.36 and Figures S.25--S.26. The supplementary analysis using all 1400 significantly differentially expressed gene probesets is reported in Table S.37 and Figures S.27--S.28. Both analyses compare SDAR-based Nodewise Loreg, Nodewise Lasso, CLIME, and GLasso.

The MDA133 dataset was obtained from the MD Anderson Cancer Center public datasets:
https://bioinformatics.mdanderson.org/public-datasets

The dataset contains dChip-normalized model-based expression index (MBEI) values for 22,283 Affymetrix gene probesets from 133 breast cancer patients treated with neoadjuvant chemotherapy. The clinical outcome is treatment response, recorded as pathologic complete response (pCR) or residual disease (RD). In the paper, the RD group corresponds to the non-pCR group.

The real-data files are organized as follows.

### Raw data and preprocessing code

- `reproduce/Real_data/Genetic_data_processing/MDA133.csv`: the original MDA133 gene-expression data. The first column contains gene probeset names, and the remaining columns contain the expression measurements for the 133 patients.

- `reproduce/Real_data/Genetic_data_processing/MDA133_L.csv`: the patient-label file. This file contains the patient IDs and their response labels. In the preprocessing script, patients with label 1 are assigned to the pCR group, and patients with label 0 are assigned to the RD/non-pCR group.

- `reproduce/Real_data/Genetic_data_processing/Genetic_data_processing.R`: the preprocessing script. This script reads `MDA133.csv` and `MDA133_L.csv`, checks missing values, performs the `log2(MBEI + 1)` transformation, separates pCR and RD/non-pCR patients, and screens differentially expressed gene probesets using the `limma` moderated t-test. It generates both the 300-probeset processed datasets used in the main analysis and the 1400-probeset processed datasets used in the supplementary ultra-high-dimensional analysis.

### Processed log-transformed data

The processed data files are stored in `reproduce/Real_data/Genetic_data_processing/Log2/`.

- `PCR_Test_log2_genedata_300.csv`: log-transformed expression data for the pCR patients, restricted to the 300 most significantly differentially expressed gene probesets selected by the `limma` moderated t-test. This file is used in the main real-data analysis with p=300.

- `RD_Test_log2_genedata_300.csv`: log-transformed expression data for the RD/non-pCR patients, restricted to the same 300 selected gene probesets. This file is used together with `PCR_Test_log2_genedata_300.csv` in the main real-data analysis with p=300.

- `PCR_Test_log2_genedata_lfdr1400.csv`: log-transformed expression data for the pCR patients, restricted to the 1400 significantly differentially expressed gene probesets selected with FDR controlled at 0.05 using the AdaptZ procedure. This file is used in the supplementary ultra-high-dimensional real-data analysis with p=1400.

- `RD_Test_log2_genedata_lfdr1400.csv`: log-transformed expression data for the RD/non-pCR patients, restricted to the same 1400 selected gene probesets. This file is used together with `PCR_Test_log2_genedata_lfdr1400.csv` in the supplementary ultra-high-dimensional real-data analysis with p=1400.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Approximate runtime

The following rough runtime guidance is provided to help users plan reproducibility runs.

### Quick example

The example script in `software/Example/Test.R` is intended as a quick demonstration of how to call the main Nodewise Loreg function. On a standard desktop or workstation, this script is expected to finish within a few minutes.

### Main SDAR-based Nodewise Loreg simulations

The expected runtime mainly depends on the dimension \(p\), and whether the user reproduces only the numerical summaries or also the Z-score histogram figures.

For one simulation setting, the approximate runtimes are as follows.

| Dimension | Numerical summaries reproduced | Approximate runtime for one setting |
|---|---|---|
| p=200 | Matrix norm losses and support recovery metrics | 30 minutes--1 hour |
| p=200 | Asymptotic-normality metrics or Z-score histograms | 30 minutes--2 hours |
| p=400 | Matrix norm losses and support recovery metrics | 1--3 hours |
| p=400 | Asymptotic-normality metrics or Z-score histograms | 1--5 hours |
| p=1000 | Matrix norm losses and support recovery metrics | 4--8 hours |
| p=1000 | Asymptotic-normality metrics or Z-score histograms | 4 hours to more than 8 hours |
| p=4000 | Matrix norm losses and support recovery metrics | >8 hours |
| p=4000 | Asymptotic-normality metrics or Z-score histograms | >8 hours |

Here, one setting refers to one fixed combination of dimension, sample size, graph structure, and simulation design, with results aggregated over all replications. The full simulation study reported in the manuscript and Supplementary Materials contains many such settings, including multiple graph structures, dimensions, and sample sizes. Therefore, reproducing the full set of simulation results from scratch can take substantially longer than the runtime for a single setting.

Users who only want to check that the code runs correctly should first run the quick example in `software/Example/Test.R`, or reduce the number of replications in the simulation scripts before running the full experiments.

### Real-data analysis

The real-data analysis scripts in `reproduce/Real_data/Code_for_p300_BH/` reproduce the main MDA133 analysis with p = 300 selected gene probesets. These scripts are expected to take about 2-4 hours.

The scripts in `reproduce/Real_data/Code_for_p1400_AdaptZ/` reproduce the supplementary MDA133 analysis with p = 1400 selected gene probesets. These scripts are expected to take >8 hours.



<p align="right">(<a href="#readme-top">back to top</a>)</p>
