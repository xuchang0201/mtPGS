---
layout: page
title: Data Input
description: ~
---
The following are the example data inputs required to construct PGS for a target trait with the aid of one relevant trait. 
1. GWAS summary statistics of the two traits from [GEMMA](https://github.com/genetics-statistics/GEMMA)
When four sets of summary statistics are available (one for overlapped individuals and one for non-overlapped individuals for the two traits)
  * [Summary statistics of overlapped individuals for trait 1](https://github.com/yuanzhongshang/GIFT/blob/main/example/Zx.txt)
  * [Summary statistics of non-overlapped individuals for trait 1](https://github.com/yuanzhongshang/GIFT/blob/main/example/Zy.txt)
  * [Summary statistics of overlapped individuals for trait 2](https://github.com/yuanzhongshang/GIFT/blob/main/example/X.txt)
  * [Summary statistics of non-overlapped individuals for trait 2](https://github.com/yuanzhongshang/GIFT/blob/main/example/Y.txt)
 When two sets of summary statistics are available (one for each trait)
  * [Summary statistics for trait 1](https://github.com/yuanzhongshang/GIFT/blob/main/example/Zx.txt)
  * [Summary statistics for trait 2](https://github.com/yuanzhongshang/GIFT/blob/main/example/X.txt)
  
2. Reference panel for computing LD matrix
  * [Reference panel](https://github.com/yuanzhongshang/GIFT/blob/main/example/Zscore1.txt)

3. LD block information from [Berisa et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731402/)
  * [European](https://github.com/xuchang0201/mtPGS/blob/main/data/EUR_LD_Block.txt)
 
4. Genetic and environmental variance component matrices
  * [Genetic variance component](https://github.com/xuchang0201/mtPGS/blob/main/data/v_g.txt)
  * [Environmental variance component](https://github.com/xuchang0201/mtPGS/blob/main/data/v_e.txt)
