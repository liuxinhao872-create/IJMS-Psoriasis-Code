# IJMS-Psoriasis-Code
Code for the paper: A Diagnostic Model for Psoriasis Based on Apoptosis-Related Genes and Explainable Machine Learning Algorithms.
# Psoriasis Apoptosis Diagnostic Model

This repository contains the official R scripts and analytical workflow for the manuscript: 
"A Diagnostic Model for Psoriasis Based on Apoptosis-Related Genes and Explainable Machine Learning Algorithms".

## Overview
This study integrates apoptosis-associated transcriptomic signatures with explainable machine learning algorithms to construct a molecular classification model for psoriasis. The workflow includes differential expression analysis, protein-protein interaction (PPI) network construction, systematic comparison of eight machine learning algorithms (RF, SVM, GLM, GBM, KNN, NNET, LASSO, and DT), and feature extraction using the DALEX explainability framework. 

The final 5-gene Gradient Boosting Machine (GBM) model (incorporating *CCNB1, KIF11, HDAC1, TPX2, and MELK*) was translated into a clinical nomogram and externally validated.

## Data Availability
The transcriptomic datasets analyzed in this study are publicly available from the NCBI Gene Expression Omnibus (GEO) repository:
* **Training Cohort:** [GSE30999](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30999) and [GSE53552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53552)
* **External Validation Cohort:** [GSE55201](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55201)

Apoptosis-related genes were retrieved from the [GeneCards](https://www.genecards.org/) database.

## System Requirements and Dependencies
* **R version:** >= 4.3.0
* **Key R Packages:** * Data Processing: `limma`
  * Machine Learning: `caret`, `pROC`
  * Explainability: `DALEX`
  * Visualization & Clinical Translation: `rms`, `circlize`, `ggplot2`
* **External Software:** Cytoscape (v3.10.4) with the *cytoHubba* plugin (for PPI network hub-gene identification).

## Repository Structure
The scripts are organized sequentially according to the analytical pipeline described in the manuscript:

* **`01_Data_Preprocessing.R`** Downloads microarray data, performs probe annotation, merges datasets, and applies batch effect correction (`removeBatchEffect`).
* **`02_DEG_and_Intersection.R`** Performs differential expression analysis and intersects DEGs with apoptosis-related genes.
* **`03_Hub_Genes_Analysis.R`** Contains code for Pearson correlation analysis and chromosomal localization visualization (`circlize`) of the top 25 hub genes identified via Cytoscape.
* **`04_Machine_Learning_Comparison.R`** Trains and evaluates the 8 machine learning algorithms using 5-fold cross-validation. Generates ROC curves, standard classification metrics, and residual distribution plots.
* **`05_Explainability_and_Nomogram.R`** Applies the `DALEX` framework to the optimal GBM model to extract permutation-based feature importance. Constructs the diagnostic nomogram, calibration curves, and conducts Decision Curve Analysis (DCA).
* **`06_External_Validation.R`** Evaluates the generalization performance of the final 5-gene model on the independent external validation cohort (GSE55201).

## Usage
1. Clone this repository to your local machine.
2. Ensure all required R packages are installed.
3. Run the scripts sequentially from `01` to `06`. Note: Ensure that the appropriate raw data matrices are downloaded from GEO and placed in the working directory before running `01_Data_Preprocessing.R`.

## Citation
If you use this code or our findings in your research, please cite our paper:
> *(Note: The citation will be updated upon publication)* > Liu X, Fu W, Zhu X, Bo W. A Diagnostic Model for Psoriasis Based on Apoptosis-Related Genes and Explainable Machine Learning Algorithms. *International Journal of Molecular Sciences* (2026).

## Contact
For any questions regarding the code or methodology, please open an issue in this repository or contact the corresponding authors.
