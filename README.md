# Using supervised learning methods for gene selection in RNA-Seq case-control studies

This git repository contains the additional files used in the following research article: https://www.biorxiv.org/content/early/2018/03/15/282780.1

### Pre-requisites

* TCGAbiolinks
* dplyr
* DT
* DESeq2
* stringr
* ranger
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

Please note that
- Some variables may have to be modified in the R scripts (eg. *baseDir*) according to your system configuration.
- These R scripts rely heavily on external libraries which might be subject to change (especially *TCGAbiolinks*). 
If that were the case, some errors might occur.

### Citation

Using supervised learning methods for gene selection in RNA-Seq case-control studies
Stephane Wenric, Ruhollah Shemirani
bioRxiv 282780; doi: https://doi.org/10.1101/282780 