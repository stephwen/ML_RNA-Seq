# Using supervised learning methods for gene selection in RNA-Seq case-control studies

This git repository contains the additional files used in the following research article: ...

### Pre-requisites

* TCGAbiolinks
* dplyr
* DT
* DESeq2
* stringr
* randomForest
* ROCR
* caret
* survival

### Running the code

The source code is currently split into different R source files, corresponding to the following steps:

1. Data download from TCGA, save HTSeq-count files for all samples
2. Run DESeq2 (normalization and differential expression analysis)
3. Run random forests
4. Survival analysis

### Extreme Pseudo-Samples

The EPS method used to extract the other ranking of genes is run in parallel of step 3.
The source code for the EPS method is available at: https://github.com/roohy/Extreme-Pseudo-Sampler/