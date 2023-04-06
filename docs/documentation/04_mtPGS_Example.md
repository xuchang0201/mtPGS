---
layout: page
title: mtPGS Tutorial
description: ~
---
This page provides a tutorial for PGS construction using mtPGS. Before runing the example code, make sure that the mtPGS software is installed and compiled successfully. For instructions for installation, please see the [Installation section](https://xuchang0201.github.io/mtPGS/documentation/02_installation.html).

## mtPGS
The example data for mtPGS tutorial can be downloaded in this [page](https://xuchang0201.github.io/mtPGS/documentation/03_data.html). Here are the details about the required data input illustrated. 
### 1. Formats of input data for mtPGS
* GWAS summary statistics: We require the GWAS summary statistics in [GEMMA](https://github.com/genetics-statistics/GEMMA) format, with the following columns: chr, rs, ps, n_mis, n_obs, allele1, allele0, af, beta, se, p_wald. Please be advised that each column should be separated using tab. 
* Reference panel: The reference panel for LD matrix should be in PLINK binary format (bed, bim, and fam).
* LD block information: We use the LD block information from [Berisa et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731402/), and we provide a formatted text file in the [data](https://github.com/xuchang0201/mtPGS/tree/main/data) as input for mtPGS.
* Genetic and environmental variance component matrices: These matrices should be formatted as a plain text, see [example](https://github.com/xuchang0201/mtPGS/blob/main/data/v_g.txt). In this example, we assume that the heritability of the two traits is 0.5, and the genetic covariance is 0.25. We recommend using [GECKO](https://github.com/borangao/GECKO) for the estimation of genetic and environmental variance components.

### 2. Running mtPGS with four sets of GWAS summary statistics
When four sets of GWAS summary statistics are available (one for overlapped individuals and one for non-overlapped individuals for each trait), and we designate trait 1 as the target trait, the PGS construction for the target trait can be performed using the following command
```
workdir=/your/mtPGS/directory #specify the mtPGS directory
${workdir}/mtPGS --summstat_int ${workdir}/data/summstat/summstat_trait_1_int.assoc.txt ${workdir}/data/summstat/summstat_trait_2_int.assoc.txt \
--summstat_ext ${workdir}/data/summstat/summstat_trait_1_ext.assoc.txt ${workdir}/data/summstat/summstat_trait_2_ext.assoc.txt \
--n_s 7000 --n_ext 3000 3000 --block ${workdir}/data/block.txt --target 0 --ref ${workdir}/data/ref --mafMax 0.8 \
--vg v_g.txt --ve v_e.txt --output trait_1_target_beta trait_2_relevant_beta
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

### 3. Running mtPGS with two sets of GWAS summary statistics
When only two sets of summary statistics are available (one for each trait), we perform mtPGS analysis using the following command
```
workdir=/your/mtPGS/directory #specify the mtPGS directory
${workdir}/mtPGS --summstat ${workdir}/data/summstat/summstat_trait_1.assoc.txt ${workdir}/data/summstat/summstat_trait_2.txt \
--n 10000 10000 --block ${workdir}/data/block.txt --target 0 --ref ${workdir}/data/ref --mafMax 0.8 --vg v_g.txt --ve v_e.txt \
--output trait_1_target_beta trait_2_relevant_beta
```
Here, instead of summstat_int, summstat_ext and n_s, n_ext, we use
- summstat: specify the GWAS summary statistics for the two traits
- n: specify sample sizes for the two traits
