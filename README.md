# Multi-trait assisted Polygenic Scores (mtPGS)

mtPGS is a method that constructs accurate PGS for a target trait of interest through leveraging multiple traits relevant to the target trait. Specifically, mtPGS borrows SNP effect size similarity information between the target trait and its relevant traits to improve the effect size estimation on the target trait, thus achieving accurate PGS. In the process, mtPGS flexibly models the shared genetic architecture between the target and the relevant traits to achieve robust performance, while explicitly accounting for the environmental covariance among them to accommodate different study designs with various sample overlap patterns. In addition, mtPGS uses only summary statistics as input and relies on a deterministic algorithm with several algebraic techniques for scalable computation.


# Getting Started
For mtPGS installation, please clone this repository via the commands

    git clone https://github.com/xuchang0201/mtPGS.git
    cd mtPGS

# Citations

Chang Xu, Santhi K. Ganesh, and Xiang Zhou (2023). mtPGS: Leverage multiple correlated traits for accurate polygenic score construction.

# Questions 
If you have any questions on mtPGS software, please email to Chang Xu (wchung@hsph.harvard.edu).
