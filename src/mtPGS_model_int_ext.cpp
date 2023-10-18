/*
Multi-trait assisted Polygenic Scores (mtPGS)
*/

#include <iostream>
#include <vector>
#include <armadillo>
#include <string>
#include "omp.h"

#include "mtPGS_function.hpp"
#include "mtPGS_model_int_ext.hpp"

using namespace std;
using namespace arma;

// est_w_ext function for general Model with large and small effect SNPs
int mtPGSFIT::est_w_ext(int n_ref, int n_s, vector <int> n_ext_vec, int num_block, mat v_e, mat v_g, vector<int> idv, string bed_str, vector < vector <INFO> > summstat_l_int,
	                     vector < vector <INFO> > summstat_s_int, vector < vector <INFO> > summstat_l_ext, vector < vector <INFO> > summstat_s_ext, 
	                     vector < vector <EFF> > &eff_output_l, vector < vector <EFF> > &eff_output_s) {

	// number of traits
	int traits_num = v_g.n_rows;

	// get the maximum number of each block
	int count_l = 0, count_s = 0; //index for going through all SNPs
	vec num_each_block_l = zeros<vec>(num_block), num_each_block_s = zeros<vec>(num_block);
	for (int i = 0; i < num_block; i++) {
		for (size_t j = count_l; j < summstat_l_int[0].size(); j++) {
			if (summstat_l_int[0][j].block == i) {
				num_each_block_l(i) += 1;
				count_l++;
			}
			else {
				break;
			}
		}
		for (size_t j = count_s; j < summstat_s_int[0].size(); j++) {
			if (summstat_s_int[0][j].block == i) {
				num_each_block_s(i) += 1;
				count_s++;
			}
			else {
				break;
			}
		}
	}

	double len_max_l = num_each_block_l.max();
	cout << len_max_l << " is the maximum number of large effect SNPs within a block." << endl;

	double len_max_s = num_each_block_s.max();
	cout << len_max_s << " is the maximum number of small effect SNPs within a block." << endl;

	// Specify some matrices
	// Omega matrix
	mat omega = inv(v_e);
	omega *= (double)n_s;
	for (int i = 0; i < (int)omega.n_rows; i++) {
		omega.diag(i) += ((double)n_ext_vec[i] / v_e(i, i));
	}
	omega.print("matrix omega:\n");

	// eigen decomposition
	mat eigen_d_matrix = sqrtmat_sympd(omega) * v_g * sqrtmat_sympd(omega);
	eigen_d_matrix.print("matrix for eigen decomposition:\n");

	vec D_eigenval;
	mat U_lambda;
	eig_sym(D_eigenval, U_lambda, eigen_d_matrix);
	mat D_lambda = diagmat(D_eigenval);
	U_lambda.print("matrix U_lambda:\n");
	D_lambda.print("matrix D_lambda:\n");

	// V_e_large
	int n_snp_total = summstat_l_int[0].size() + summstat_s_int[0].size();
	int n_large_snp = summstat_l_int[0].size();
	mat v_e_large = zeros<mat>(2, 2);
	for (int i = 0; i < (int)v_e_large.n_rows; i++) {
		v_e_large(i, i) = 1.0 - ((double)n_large_snp / (double)n_snp_total) * v_g(i, i);
	}

	double rho_e = v_e(0, 1) / (sqrt(v_e(0, 0)) * sqrt(v_e(1, 1)));
	v_e_large(0, 1) = rho_e * sqrt(v_e_large(0, 0)) * sqrt(v_e_large(1, 1));
	v_e_large(1, 0) = rho_e * sqrt(v_e_large(0, 0)) * sqrt(v_e_large(1, 1));
	v_e_large.print("matrix V_e_large:\n");

	// loop 
	int count_l_int = 0, count_s_int = 0, count_l_ext = 0, count_s_ext = 0; //reset the count of all SNPs

	for (int i = 0; i < num_block; ++i) {

		vector < vector <INFO> > info_all_within_l_int(traits_num), info_all_within_s_int(traits_num);
		vector < vector <INFO> > info_all_within_l_ext(traits_num), info_all_within_s_ext(traits_num);

		// Internal data
		// SNP information of large effect for each trait within the ith block
		for (size_t j = count_l_int; j < summstat_l_int[0].size(); j++) {
			if (summstat_l_int[0][j].block == i) {
				for (int k = 0; k < traits_num; k++) {
					info_all_within_l_int[k].push_back(summstat_l_int[k][j]);
				}
				count_l_int++;
			}
			else {
				break;
			}
		}

		// SNP information of small effect for each trait within the ith block
		for (size_t j = count_s_int; j < summstat_s_int[0].size(); j++) {
			if (summstat_s_int[0][j].block == i) {
				for (int k = 0; k < traits_num; k++) {
					info_all_within_s_int[k].push_back(summstat_s_int[k][j]);
				}
				count_s_int++;
			}
			else {
				break;
			}
		}

		// External data
		// large
		for (size_t j = count_l_ext; j < summstat_l_ext[0].size(); j++) {
			if (summstat_l_ext[0][j].block == i) {
				for (int k = 0; k < traits_num; k++) {
					info_all_within_l_ext[k].push_back(summstat_l_ext[k][j]);
				}
				count_l_ext++;
			}
			else {
				break;
			}
		}

		// small
		for (size_t j = count_s_ext; j < summstat_s_ext[0].size(); j++) {
			if (summstat_s_ext[0][j].block == i) {
				for (int k = 0; k < traits_num; k++) {
					info_all_within_s_ext[k].push_back(summstat_s_ext[k][j]);
				}
				count_s_ext++;
			}
			else {
				break;
			}
		}

		// perform estimation within block
		vector < vector <EFF> > eff_Block_l(traits_num), eff_Block_s(traits_num);
		calcBlock_w_ext(n_ref, n_s, n_ext_vec, v_g, v_e, omega, U_lambda, D_lambda, v_e_large, idv, bed_str, info_all_within_l_int, info_all_within_s_int,
			            info_all_within_l_ext, info_all_within_s_ext, num_each_block_l[i], num_each_block_s[i], eff_Block_l, eff_Block_s);

		// output the estimated effect sizes
		for (int m = 0; m < traits_num; m++) {
			for (int l = 0; l < num_each_block_l[i]; l++) {
				eff_output_l[m].push_back(eff_Block_l[m][l]);
			}
		}

		for (int m = 0; m < traits_num; m++) {
			for (int l = 0; l < num_each_block_s[i]; l++) {
				eff_output_s[m].push_back(eff_Block_s[m][l]);
			}
		}
	}
	return 0;
}

// calcBlock_w_ext function for general Model with large and small effect SNPs
int mtPGSFIT::calcBlock_w_ext(int n_ref, int n_s, vector <int> n_ext_vec, mat v_g, mat v_e, mat omega, mat U_lambda, mat D_lambda, mat v_e_large, vector<int> idv, string bed_str,
	                           vector < vector <INFO> > info_all_within_l_int, vector < vector <INFO> > info_all_within_s_int, vector < vector <INFO> > info_all_within_l_ext, 
	                           vector < vector <INFO> > info_all_within_s_ext, int num_in_block_l, int num_in_block_s, vector < vector <EFF> > &eff_Block_l, vector < vector <EFF> > &eff_Block_s) {
	SNPPROC cSP;
	IO cIO;
	ifstream bed_in(bed_str.c_str(), ios::binary);

	// import the summary statistics of LARGE effect SNPs within the block
	// internal
	vector < vector <INFO> > info_within_block_l_int((int)v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_l; j++) {
			info_within_block_l_int[i].push_back(info_all_within_l_int[i][j]);
		}
	}
	
	// external
	vector < vector <INFO> > info_within_block_l_ext((int)v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_l; j++) {
			info_within_block_l_ext[i].push_back(info_all_within_l_ext[i][j]);
		}
	}


	// z_score matrix for LARGE effect SNPs
	// internal
	mat z_matrix_l_int = zeros<mat>(num_in_block_l, v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int k = 0; k < num_in_block_l; k++)
			z_matrix_l_int(k, i) = info_within_block_l_int[i][k].z;
	}

	// external
	mat z_matrix_l_ext = zeros<mat>(num_in_block_l, v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int k = 0; k < num_in_block_l; k++)
			z_matrix_l_ext(k, i) = info_within_block_l_ext[i][k].z;
	}

	// reference panel genotype matrix of LARGE effect SNPs
	mat geno_l = zeros<mat>(n_ref, num_in_block_l);
	for (int i = 0; i < num_in_block_l; ++i) {
		vec geno = zeros<vec>(n_ref);
		double maf = 0.0;
		cIO.readSNPIm(info_within_block_l_int[0][i].pos, n_ref, idv, bed_in, geno, maf);
		cSP.nomalizeVec(geno);
		geno_l.col(i) = geno;
	}

	// import the summary statistics of SMALL effect SNPs within the block
	// internal
	vector < vector <INFO> > info_within_block_s_int((int)v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_s; j++) {
			info_within_block_s_int[i].push_back(info_all_within_s_int[i][j]);
		}
	}

	// external
	vector < vector <INFO> > info_within_block_s_ext((int)v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_s; j++) {
			info_within_block_s_ext[i].push_back(info_all_within_s_ext[i][j]);
		}
	}

	// z_score matrix for SMALL effect SNPs
	// internal
	mat z_matrix_s_int = zeros<mat>(num_in_block_s, v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int k = 0; k < num_in_block_s; k++)
			z_matrix_s_int(k, i) = info_within_block_s_int[i][k].z;
	}

	// external
	mat z_matrix_s_ext = zeros<mat>(num_in_block_s, v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int k = 0; k < num_in_block_s; k++)
			z_matrix_s_ext(k, i) = info_within_block_s_ext[i][k].z;
	}

	// reference panel genotype matrix of SMALL effect SNPs
	mat geno_s = zeros<mat>(n_ref, num_in_block_s);
	for (int i = 0; i < num_in_block_s; ++i) {
		vec geno = zeros<vec>(n_ref);
		double maf = 0.0;
		cIO.readSNPIm(info_within_block_s_int[0][i].pos, n_ref, idv, bed_in, geno, maf);
		cSP.nomalizeVec(geno);
		geno_s.col(i) = geno;
	}

	//estimation
	mat beta_output_l = zeros<mat>(num_in_block_l, (int)v_g.n_rows);
	mat beta_output_s = zeros<mat>(num_in_block_s, (int)v_g.n_rows);
	estBlock_w_ext(n_ref, n_s, n_ext_vec, v_g, v_e, omega, U_lambda, D_lambda, v_e_large, geno_l, geno_s, z_matrix_l_int, z_matrix_s_int, z_matrix_l_ext, z_matrix_s_ext, beta_output_l, beta_output_s);

	// output effect size for each trait
	// large effect
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_l; j++) {
			EFF eff;
			eff.snp = info_within_block_l_int[i][j].snp;
			eff.a1 = info_within_block_l_int[i][j].a1;
			eff.maf = info_within_block_l_int[i][j].maf;
			eff.beta = beta_output_l(j, i);
			eff_Block_l[i].push_back(eff);
		}
	}

	//small effect
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_s; j++) {
			EFF eff;
			eff.snp = info_within_block_s_int[i][j].snp;
			eff.a1 = info_within_block_s_int[i][j].a1;
			eff.maf = info_within_block_s_int[i][j].maf;
			eff.beta = beta_output_s(j, i);
			eff_Block_s[i].push_back(eff);
		}
	}

	return 0;
}


// solve the equation Ax=b, x=A^{-1}b
vec mtPGSFIT::PCGv(mat A, vec b, size_t maxiter, const double tol){
	vec dA = A.diag();
	// stable checking, using available func to speed up
	for(size_t i = 0; i< dA.n_elem; i++){
		if(dA[i] == 0)
			dA[i] = 1e-4;
	}// end for i
	vec Minv = 1.0/dA;
	// initialize
	vec x = zeros<vec>(b.n_elem);
	vec r = zeros<vec>(b.n_elem);
	vec r1 = zeros<vec>(b.n_elem);
	vec z1 = zeros<vec>(b.n_elem);
	r = b;
	vec z = Minv % r;
	vec p = z;
	size_t iter = 0;
	double sumr2 = norm(r, 2);
	// PCG main loop 
	while( sumr2 > tol && iter < maxiter){
		iter += 1;
		// move direction
		vec Ap = A*p;
		// step size
		double a = dot(r,z)/dot(p,Ap);
		// move
		x = x + a * p;
		r1 = r - a * Ap;
		z1 = Minv % r1;
		double bet = dot(z1, r1)/dot(z, r);
		p = z1 + bet * p;
		z = z1;
		r = r1;
		sumr2 = norm(r, 2);
	}// end while loop
	if (iter >= maxiter){
		cerr << "ERROR: Matrix is Singular!" << endl;
	}
	return(x);
}// end function


mat mtPGSFIT::PCGm(mat A, mat B, size_t maxiter, const double tol){
	
	size_t n_iter = B.n_cols;
	mat x = zeros<mat>(A.n_rows, n_iter);
	for (size_t i = 0; i < n_iter; i++){
		x.col(i) = PCGv(A, B.col(i), maxiter, tol);
	}// end for loop
	return(x);
}// end function


// General Model with both large and small effect SNPs, overlapped and non-overlapped individuals
int mtPGSFIT::estBlock_w_ext(int n_ref, int n_s, vector <int> n_ext_vec, mat v_g, mat v_e, mat omega, mat U_lambda, mat D_lambda,
	                          mat v_e_large, mat geno_l, mat geno_s, mat z_matrix_l, mat z_matrix_s, mat z_l_external, 
	                          mat z_s_external, mat &beta_output_l, mat &beta_output_s) {

	// LD matrix
	mat SIGMA_ref_l = geno_l.t() * geno_l;
	SIGMA_ref_l /= (double)n_ref;

	mat SIGMA_ref_s = geno_s.t() * geno_s;
	SIGMA_ref_s /= (double)n_ref;

	mat SIGMA_ref_sl = geno_s.t() * geno_l;
	SIGMA_ref_sl /= (double)n_ref;

	// large effect SNPs
	// v_e_large inverse and its transpose
	mat v_e_l_inv = inv(v_e_large);
	mat v_e_l_inv_trans = v_e_l_inv.t();

	mat sigma_ll_inv_z_sum = zeros<mat>((int)z_matrix_l.n_rows, (int)v_g.n_rows);

	for (int i = 0; i < (int)v_g.n_rows; i++) {
		int n_ext = n_ext_vec[i];
		vec z_vec_sum_internal = z_matrix_l * v_e_l_inv_trans.col(i);
		z_vec_sum_internal *= sqrt((double)n_s);
		vec z_vec_ext = z_l_external.col(i);
		z_vec_ext *= (sqrt((double)n_ext)/v_e_large(i, i));
		vec z_vec_sum = z_vec_sum_internal + z_vec_ext;
		sigma_ll_inv_z_sum.col(i) = solve(SIGMA_ref_l, z_vec_sum);
	}

	// inverse of the sum of n_s * v_e_l_inv + diag(n_ext) and its transpose
	v_e_l_inv *= (double)n_s;
	for (int i = 0; i < (int)v_e_l_inv.n_rows; i++) {
		v_e_l_inv.diag(i) += ((double)n_ext_vec[i]/v_e_large(i, i));
	}
	mat v_e_inv_times_n_inv = inv(v_e_l_inv);
	mat v_e_inv_times_n_inv_trans = v_e_inv_times_n_inv.t();

	// calculate beta for large effect SNPs
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		beta_output_l.col(i) = sigma_ll_inv_z_sum * v_e_inv_times_n_inv_trans.col(i);
	}

	// small effect SNPs 
	// inverse of v_e and its transpose
	mat v_e_inv = inv(v_e);
	mat v_e_inv_trans = v_e_inv.t();

	// Step 1
	// define block-wise computation matrices for overlap and non-overlap
	mat z_sigma_sl_beta_l_diff_overlap = zeros<mat>((int)z_matrix_s.n_rows, (int)v_g.n_rows);
	mat z_sigmasl_betal_diff_non_over = zeros<mat>((int)z_matrix_s.n_rows, (int)v_g.n_rows);

	for (int i = 0; i < (int)v_g.n_rows; i++) {
		// SIGMA_sl * beta_l
		vec sigma_sl_beta_l = SIGMA_ref_sl * beta_output_l.col(i);

		// overlap
		z_sigma_sl_beta_l_diff_overlap.col(i) = z_matrix_s.col(i) - sqrt((double)n_s) * sigma_sl_beta_l;
		z_sigma_sl_beta_l_diff_overlap.col(i) *= sqrt((double)n_s);

		// non-overlap
		double env_comp = v_e(i, i);
		int n_sample_ext = n_ext_vec[i];
		z_sigmasl_betal_diff_non_over.col(i) = z_s_external.col(i) - sqrt((double)n_sample_ext) * sigma_sl_beta_l;
		z_sigmasl_betal_diff_non_over.col(i) *= (sqrt((double)n_sample_ext)/env_comp);
	}

	// overlap part component (need to times v_e_inverse)
	mat v_e_inv_z_sigma_beta_overlap = zeros<mat>((int)z_matrix_s.n_rows, (int)v_g.n_rows);

	for (int i = 0; i < (int)v_g.n_rows; i++) {
		v_e_inv_z_sigma_beta_overlap.col(i) = z_sigma_sl_beta_l_diff_overlap * v_e_inv_trans.col(i);
	}

	// step 1 result (sum of overlap and non-overlap components)
	mat sum_z_sigma_beta_over_non_over = v_e_inv_z_sigma_beta_overlap + z_sigmasl_betal_diff_non_over;

	// Step 2
	// transpose of vg
	mat vg_trans = v_g.t();

	// t(U_lambda) * sqrt(omega) * v_g
	mat U_lambda_omega_vg = U_lambda.t() * sqrtmat_sympd(omega) * v_g;
	mat U_lambda_omega_vg_trans = U_lambda_omega_vg.t();

	// v_g * sqrt(omega) * U_lambda
	mat v_g_omega_U_lambda = v_g * sqrtmat_sympd(omega) * U_lambda;
	mat v_g_omega_U_lambda_trans = v_g_omega_U_lambda.t();

	// U_lambda*omega*v_g * vector 
	mat U_omega_vg_sum_z_sigma_beta_over_non_over = zeros<mat>((int)z_matrix_s.n_rows, (int)v_g.n_rows);

	for (int i = 0; i < (int)v_g.n_rows; i++) {
		U_omega_vg_sum_z_sigma_beta_over_non_over.col(i) = sum_z_sigma_beta_over_non_over * U_lambda_omega_vg_trans.col(i);
	}

	// SIGMA_ss_D_SIGMA_ss_inv_vector
	mat SIGMA_ss_D_SIGMA_ss_inv_vector = zeros<mat>((int)z_matrix_s.n_rows, (int)v_g.n_rows);

	for (int i = 0; i < (int)v_g.n_rows; i++) {
		SIGMA_ref_s *= (D_lambda(i, i) / 2.0);
		SIGMA_ref_s.diag() += 1.0;
		vec SIGMA_ss_inv_vector = solve(SIGMA_ref_s, U_omega_vg_sum_z_sigma_beta_over_non_over.col(i));
		// reset SIGMA_ss
		SIGMA_ref_s.diag() -= 1.0;
		SIGMA_ref_s /= (D_lambda(i, i) / 2.0);
		// SIGMA_ss * SIGMA_ss_inv_vector
		SIGMA_ss_D_SIGMA_ss_inv_vector.col(i) = SIGMA_ref_s * SIGMA_ss_inv_vector;
	}

	// calculate beta for small effect SNPs
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		vec v_g_sum_z_sigma_beta_over_non_over = sum_z_sigma_beta_over_non_over * vg_trans.col(i);
		v_g_sum_z_sigma_beta_over_non_over /= 2.0;
		vec v_g_omega_U_lambda_SIGMA_ss_D_SIGMA_ss_inv_vector = SIGMA_ss_D_SIGMA_ss_inv_vector * v_g_omega_U_lambda_trans.col(i);
		v_g_omega_U_lambda_SIGMA_ss_D_SIGMA_ss_inv_vector /= 4.0;
		beta_output_s.col(i) = v_g_sum_z_sigma_beta_over_non_over - v_g_omega_U_lambda_SIGMA_ss_D_SIGMA_ss_inv_vector;
	}

	return 0;
}

