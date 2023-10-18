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
	// est_w_ext function for general Model with large and small effect SNPs
	int est_w_ext(int n_ref, int n_s, vector <int> n_ext_vec, int num_block, mat v_e, mat v_g, vector<int> idv, string bed_str, vector < vector <INFO> > summstat_l_int,
		vector < vector <INFO> > summstat_s_int, vector < vector <INFO> > summstat_l_ext, vector < vector <INFO> > summstat_s_ext,
		vector < vector <EFF> >& eff_output_l, vector < vector <EFF> >& eff_output_s);
	// calcBlock_w_ext function for general Model with large and small effect SNPs
	int calcBlock_w_ext(int n_ref, int n_s, vector <int> n_ext_vec, mat v_g, mat v_e, mat omega, mat U_lambda, mat D_lambda, mat v_e_large, vector<int> idv, string bed_str,
		vector < vector <INFO> > info_all_within_l_int, vector < vector <INFO> > info_all_within_s_int, vector < vector <INFO> > info_all_within_l_ext,
		vector < vector <INFO> > info_all_within_s_ext, int num_in_block_l, int num_in_block_s, vector < vector <EFF> >& eff_Block_l, vector < vector <EFF> >& eff_Block_s);
	// solve x=Ab
	vec PCGv(mat A, vec b, size_t maxiter, const double tol); 
	// solve x=AB
	mat PCGm(mat A, mat B, size_t maxiter, const double tol);
	// General Model with both large and small effect SNPs, overlapped and non-overlapped individuals
	int estBlock_w_ext(int n_ref, int n_s, vector <int> n_ext_vec, mat v_g, mat v_e, mat omega, mat U_lambda, mat D_lambda,
		mat v_e_large, mat geno_l, mat geno_s, mat z_matrix_l, mat z_matrix_s, mat z_l_external,
		mat z_s_external, mat& beta_output_l, mat& beta_output_s);
};
#endif
