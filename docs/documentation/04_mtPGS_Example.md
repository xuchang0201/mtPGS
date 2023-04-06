---
layout: page
title: mtPGS Tutorial
description: ~
---
This page provides a tutorial for PGS construction using mtPGS. Before runing the example code, make sure that the mtPGS software is installed and compiled successfully. For instructions for installation, please see the [Installation section](https://xuchang0201.github.io/mtPGS/documentation/02_installation.html).

## mtPGS
The example data for mtPGS tutorial can be downloaded in this [page](https://xuchang0201.github.io/mtPGS/documentation/03_data.html). Here are the details about the required data input illustrated. 
### 1. running mtPGS with four sets of GWAS summary statistics
When four sets of GWAS summary statistics are available (one for overlapped individuals and one for non-overlapped individuals for each trait), and we designate trait 1 as the target trait, the PGS construction for the target trait can be performed using the following command
```
./mtPGS --summstat_int summstat_trait_1_int.assoc.txt summstat_trait_2_int.assoc.txt --summstat_ext summstat_trait_1_ext.assoc.txt summstat_trait_2_ext.assoc.txt --n_s 7000 --n_ext 3000 3000 --block block.txt --target 0 --ref ref_panel_ukbb_500samples --mafMax 0.8 --vg v_g.txt --ve v_e.txt --output trait_1_target_beta trait_2_relevant_beta
```
The essential inputs are:
- Zscore_1: The Zscore matrix of the cis-SNP effect size matrix, each column for one specific gene in eQTL data.
- Zscore_2: The Zscore vector of the cis-SNP effect size vector for one specific trait in GWAS data.
- Sigma1: The LD matrix in eQTL data.
- Sigma2: The LD matrix in GWAS data, both Sigma1 and Sigma2 are often the same from the reference panel.
- R: The estimated correlation matrix of gene expressions.
- n1: The sample size of eQTL data.
- n2: The sample size of GWAS data.
- pindex: A vector with each element represents the number of cis-SNPs for each gene.
- max_iterin: The maximum iteration, which can be determined by users. Default is 1000. 
- epsin: The convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, which can be determined by users. Default is 1e-4. 
- Cores: The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux. Default is 1.
