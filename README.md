# Multi-trait assisted Polygenic Scores (mtPGS)

mtPGS is a method that constructs accurate PGS for a target trait of interest through leveraging multiple traits relevant to the target trait. Specifically, mtPGS borrows SNP effect size similarity information between the target trait and its relevant traits to improve the effect size estimation on the target trait, thus achieving accurate PGS. In the process, mtPGS flexibly models the shared genetic architecture between the target and the relevant traits to achieve robust performance, while explicitly accounting for the environmental covariance among them to accommodate different study designs with various sample overlap patterns. In addition, mtPGS uses only summary statistics as input and relies on a deterministic algorithm with several algebraic techniques for scalable computation.


# Getting Started
In order to download CTPR, you should clone this repository via the commands

    git clone https://github.com/wonilchung/CTPR.git
    cd CTPR

Short tutorials describing how to install and run CTPR can be found in the following wiki: https://github.com/wonilchung/CTPR/wiki/CTPR-User-Manual.

# Citations
The CTPR algorithm is described in the following reference:

Wonil Chung, Jun Chen, Constance Turman, Sara Lindstrom, Zhaozhong Zhu, Po-Ru Loh, Peter Kraft and Liming Liang (2019), Efficient cross-trait penalized regression increases prediction accuracy in large cohorts using secondary phenotypes. Nature Communications, 10(1), 569. [Link](https://www.nature.com/articles/s41467-019-08535-0)

# Questions and Requests
If you have any questions on CTPR software, please email to Wonil Chung (wchung@hsph.harvard.edu).
