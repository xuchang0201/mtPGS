/*
Multi-trait assisted Polygenic Scores (mtPGS)
*/

#ifndef __MVLMM_MODEL_H__
#define __MVLMM_MODEL_H__

#include <iostream>
#include <vector>
#include <string>
#include <armadillo>

#include "mtPGS_function.hpp"

using namespace std;
using namespace arma;

class mtPGSFIT {
public:
	// Model with both large and small effect SNPs 
	int est(int n_ref, int n_obs, int num_block, mat v_e, mat v_g, vector<int> idv, string bed_str, vector < vector <INFO> > summstat_l,
		vector < vector <INFO> > summstat_s, vector < vector <EFF> > &eff_output_l, vector < vector <EFF> > &eff_output_s);
	// Model with only small effect SNPs
	int est(int n_ref, int n_obs, int num_block, mat v_e, mat v_g, vector<int> idv, string bed_str,
		vector < vector <INFO> > summstat_s, vector < vector <EFF> > &eff_output_s);
	// Effect size estimation within each block for model with both large and small effect SNPs
	int calcBlock(int n_ref, int n_obs, mat v_g, mat v_e, mat U_lambda, mat D_lambda, vector<int> idv, string bed_str, vector < vector <INFO> > info_all_within_l,
		vector < vector <INFO> > info_all_within_s, int num_in_block_l, int num_in_block_s, vector < vector <EFF> > &eff_Block_l, vector < vector <EFF> > &eff_Block_s);
	// Effect size estimation within each block for model with only small effect SNPs
	int calcBlock(int n_ref, int n_obs, mat v_g, mat v_e, mat U_lambda, mat D_lambda, vector<int> idv, string bed_str,
		vector < vector <INFO> > info_all_within_s, int num_in_block_s, vector < vector <EFF> > &eff_Block_s);
	// solve x=Ab
	vec PCGv(mat A, vec b, size_t maxiter, const double tol); 
	// solve x=AB
	mat PCGm(mat A, mat B, size_t maxiter, const double tol);
	// Model fitting for both large and small effect SNPs
	int estBlock(int n_ref, int n_obs, mat v_g, mat v_e, mat U_lambda, mat D_lambda, mat geno_l, mat geno_s,
		mat z_matrix_l, mat z_matrix_s, mat &beta_output_l, mat &beta_output_s);
	// Model fitting for only small effect SNPs
	int estBlock(int n_ref, int n_obs, mat v_g, mat v_e, mat U_lambda, mat D_lambda, mat geno_s, mat z_matrix_s, mat &beta_output_s);
};
#endif
