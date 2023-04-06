---
layout: page
title: mtPGS Tutorial
description: ~
---
This page provides a tutorial for PGS construction using mtPGS. Before runing the example code, make sure that the mtPGS software is installed and compiled successfully. For instructions for installation, please see the [Installation section](https://xuchang0201.github.io/mtPGS/documentation/02_installation.html).

## mtPGS
The example data for mtPGS tutorial can be downloaded in this [page](https://xuchang0201.github.io/mtPGS/documentation/03_data.html). Here are the details about the required data input illustrated. 
### 1. Formats of input data for mtPGS
### 2. Running mtPGS with four sets of GWAS summary statistics
When four sets of GWAS summary statistics are available (one for overlapped individuals and one for non-overlapped individuals for each trait), and we designate trait 1 as the target trait, the PGS construction for the target trait can be performed using the following command
```
./mtPGS --summstat_int summstat_trait_1_int.assoc.txt summstat_trait_2_int.assoc.txt --summstat_ext summstat_trait_1_ext.assoc.txt summstat_trait_2_ext.assoc.txt --n_s 7000 --n_ext 3000 3000 --block block.txt --target 0 --ref ref_panel_ukbb_500samples --mafMax 0.8 --vg v_g.txt --ve v_e.txt --output trait_1_target_beta trait_2_relevant_beta
```
The essential inputs are:
- summstat_int: specify the GWAS summary statistics computed based on overlapped individuals.
- summstat_ext: specify the GWAS summary statistics computed based on non-overlapped individuals.
- n_s: the sample size of overlapped individuals
- n_ext: the sample size of non-overlapped individuals
- block: specify the LD block information.
- target: index for the target trait (0 if the first trait is the target trait, and 1 if the second trait is the target trait)
- ref: specify the prefix of reference panel.
- mafMax: specify the maximium of the allele frequency difference between reference panel and summary data.
- vg: specify the directory of genetic variance components file.
- ve: specify the directory of environmental variance components file.
- output: specify the prefix of output files.