---
layout: page
title: mtPGS Tutorial
description: ~
---
This page provides a tutorial for PGS construction using mtPGS. Before runing the example code, make sure that the mtPGS software is installed and compiled successfully. For instructions on installation, please see the [Installation section](https://xuchang0201.github.io/mtPGS/documentation/02_installation.html).

## mtPGS
The example data for mtPGS tutorial can be downloaded in this [page](https://xuchang0201.github.io/mtPGS/documentation/03_data.html). Here are the details about the input data formats and mtPGS commands. 
### 1. Formats of input data for mtPGS
* GWAS summary statistics: We require the GWAS summary statistics in [GEMMA](https://github.com/genetics-statistics/GEMMA) format, with the following columns: chr, rs, ps, n_mis, n_obs, allele1, allele0, af, beta, se, p_wald. Please be advised that each column should be separated using tab. 
* Reference panel: The reference panel for LD matrix should be in PLINK binary format (bed, bim, and fam).
* LD block information: We use the LD block information from [Berisa et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731402/), and we provide a formatted text file in the [data](https://github.com/xuchang0201/mtPGS/tree/main/data) as input for mtPGS.
* Genetic and environmental variance component matrices: These matrices should be formatted as a plain text, see [example](https://github.com/xuchang0201/mtPGS/blob/main/data/v_g.txt). In this example, we assume that the heritability of the two traits is 0.5, and the genetic covariance is 0.25. We recommend using [GECKO](https://github.com/borangao/GECKO) for the estimation of genetic and environmental variance components.

### 2. Running mtPGS with four sets of GWAS summary statistics
When four sets of GWAS summary statistics are available (one for overlapped individuals and one for non-overlapped individuals for each trait), and we designate trait 1 as the target trait, the PGS construction for the target trait can be performed using the following command
```
workdir=/your/mtPGS/directory #specify the mtPGS directory
r2=0.2 #pre-specified LD threshold
pval=1e-6 #pre-specified p value threshold
plink=/usr/cluster/bin/plink-1.9 #the directory of plink 1.9 software
${workdir}/src/mtPGS_int_ext --summstat_int ${workdir}/data/summstat/summstat_trait_1_int.assoc.txt ${workdir}/data/summstat/summstat_trait_2_int.assoc.txt \
--summstat_ext ${workdir}/data/summstat/summstat_trait_1_ext.assoc.txt ${workdir}/data/summstat/summstat_trait_2_ext.assoc.txt \
--n_s 7000 --n_ext 3000 3000 --block ${workdir}/data/EUR_LD_Block.txt --target 0 --ref ${workdir}/data/ref --mafMax 0.8 \
--vg ${workdir}/data/v_g.txt --ve ${workdir}/data/v_e.txt --r2 ${r2} --pval ${pval} --plink ${plink} --c_t ${workdir}/data/c_t_output --output beta_1_target_four_data beta_2_relevant_four_data
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
- r2: specify the LD threshold for C+T procedure
- pval: specify the p value threshold for C+T procedure
- plink: specify the directory of plink 1.9 software
- c_t: specify the prefix for the outputs of C+T procedure
- output: specify the prefix of output files.

### 3. Running mtPGS with two sets of GWAS summary statistics
When only two sets of summary statistics are available (one for each trait), we perform mtPGS analysis using the following command
```
workdir=/your/mtPGS/directory #specify the mtPGS directory
r2=0.2 #pre-specified LD threshold
pval=1e-6 #pre-specified p value threshold
plink=/usr/cluster/bin/plink-1.9 #the directory of plink 1.9 software
${workdir}/src/mtPGS_int_only --summstat ${workdir}/data/summstat/summstat_trait_1.assoc.txt ${workdir}/data/summstat/summstat_trait_2.assoc.txt \
--n 10000 10000 --block ${workdir}/data/EUR_LD_Block.txt --target 0 --ref ${workdir}/data/ref --mafMax 0.8 --vg ${workdir}/data/v_g.txt --ve ${workdir}/data/v_e.txt \
--r2 ${r2} --pval ${pval} --plink ${plink} --c_t ${workdir}/data/c_t_output --output beta_1_target_two_data beta_2_relevant_two_data
```
Here, instead of summstat_int, summstat_ext and n_s, n_ext, we use
- summstat: specify the GWAS summary statistics for the two traits
- n: specify sample sizes for the two traits

### 4. output of mtPGS and PGS computation using PLINK
The example of output file is
```
rs12024068 A 0.0600938 0.118039 1
rs4580536 A -0.0731202 -0.133717 1
rs6856795 A -0.0558476 -0.0954865 1
rs493284 T -0.0583431 -0.0827885 1
rs12354787 T -0.000148759 -0.000365733 0
rs12768213 C -0.000779155 -0.0011521 0
rs11254302 G -0.000137537 -0.000338145 0
rs2942359 C 0.000739596 0.00116475 0
```
The output file has five columns (without header): the first column is the SNP ID, the second column is the effect allele, the third column is the scaled effect sizes, the fourth column is the non-scaled effect sizes (computed using MAF from summary statistics), and the last column is the index of whether this SNP has large or small effect (1 for large-effect SNP and 0 for small-effect SNP). This output file can be directly used to compute PGS using the score function of PLINK. To do this, please use columns 1, 2, and 4 of the mtPGS output and the example code as follows:
```
plink-1.9 --bfile test_genotype_data --score beta.txt 1 2 4 sum --out test_PGS
```
Here, test_genotype_data is the prefix of genotype data that you would like to compute PGS using the estimated effect sizes, beta.txt is the output of mtPGS, and test_PGS is the prefix of the output PGS. 
