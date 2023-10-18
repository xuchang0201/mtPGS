/*
Multi-trait assisted Polygenic Scores (mtPGS)
*/

#include <iostream>
#include <vector>
#include <armadillo>
#include <string>
#include "omp.h"

#include "mtPGS_function.hpp"
#include "mtPGS_model_int_only.hpp"

using namespace std;
using namespace arma;

// Model with both large and small effect SNPs
int mtPGSFIT::est(int n_ref, int n_obs, int num_block, mat v_e, mat v_g, vector<int> idv, string bed_str, vector < vector <INFO> > summstat_l,
	              vector < vector <INFO> > summstat_s, vector < vector <EFF> > &eff_output_l, vector < vector <EFF> > &eff_output_s) {

	// number of traits
	int traits_num = v_g.n_rows;

	// get the maximum number of each block
	int count_l = 0, count_s = 0; //index for going through all SNPs
	vec num_each_block_l = zeros<vec>(num_block), num_each_block_s = zeros<vec>(num_block);
	for (int i = 0; i < num_block; i++) {
		for (size_t j = count_l; j < summstat_l[0].size(); j++) {
			if (summstat_l[0][j].block == i) {
				num_each_block_l(i) += 1;
				count_l++;
			}
			else {
				break;
			}
		}
		for (size_t j = count_s; j < summstat_s[0].size(); j++) {
			if (summstat_s[0][j].block == i) {
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

	// eigen decomposition
	mat eigen_d_matrix = inv(sqrtmat_sympd(v_e)) * v_g * inv(sqrtmat_sympd(v_e));
	eigen_d_matrix.print("matrix for eigen decomposition:\n");

	vec D_eigenval;
	mat U_lambda;
	eig_sym(D_eigenval, U_lambda, eigen_d_matrix);
	mat D_lambda = diagmat(D_eigenval);
	U_lambda.print("matrix U_lambda:\n");
	D_lambda.print("matrix D_lambda:\n");

	// loop 
	count_l = 0, count_s = 0; //reset the count of all SNPs

	for (int i = 0; i < num_block; ++i) {

		vector < vector <INFO> > info_all_within_l(traits_num), info_all_within_s(traits_num);
		vector < vector <EFF> > eff_Block_l(traits_num), eff_Block_s(traits_num);

		// SNP information of large effect for each trait within the ith block
		for (size_t j = count_l; j < summstat_l[0].size(); j++) {
			if (summstat_l[0][j].block == i) {
				for (int k = 0; k < traits_num; k++) {
					info_all_within_l[k].push_back(summstat_l[k][j]);
				}
				count_l++;
			}
			else {
				break;
			}
		}

		// SNP information of small effect for each trait within the ith block
		for (size_t j = count_s; j < summstat_s[0].size(); j++) {
			if (summstat_s[0][j].block == i) {
				for (int k = 0; k < traits_num; k++) {
					info_all_within_s[k].push_back(summstat_s[k][j]);
				}
				count_s++;
			}
			else {
				break;
			}
		}

		calcBlock(n_ref, n_obs, v_g, v_e, U_lambda, D_lambda, idv, bed_str, info_all_within_l, info_all_within_s, 
			num_each_block_l[i], num_each_block_s[i], eff_Block_l, eff_Block_s);

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


// Model with only small effect SNPs 
int mtPGSFIT::est(int n_ref, int n_obs, int num_block, mat v_e, mat v_g, vector<int> idv, string bed_str,
				  vector < vector <INFO> > summstat_s, vector < vector <EFF> > &eff_output_s) {
	
	// number of traits
	int traits_num = v_g.n_rows;

	// get the maximum number of each block
	int count = 0; //index for going through all SNPs
	vec num_each_block = zeros<vec>(num_block); 
	for (int i = 0; i < num_block; i++) {
		for (size_t j = count; j < summstat_s[0].size(); j++) {
			if(summstat_s[0][j].block == i){ 
				num_each_block(i) += 1; 
				count++;
			}else{
				break;
			}
		}
	}
	 
	double len_max = num_each_block.max(); 
	cout << len_max << " is the maximum number of SNPs within a block." << endl;

	// eigen decomposition
	mat eigen_d_matrix = inv(sqrtmat_sympd(v_e)) * v_g * inv(sqrtmat_sympd(v_e));
	eigen_d_matrix.print("matrix for eigen decomposition:\n");

	vec D_eigenval;
	mat U_lambda;
	eig_sym(D_eigenval, U_lambda, eigen_d_matrix);
	mat D_lambda = diagmat(D_eigenval);
	U_lambda.print("matrix U_lambda:\n");
	D_lambda.print("matrix D_lambda:\n");

	// loop 
	count = 0; //reset the count of all SNPs

	for (int i = 0; i < num_block; ++i) {

		vector < vector <INFO> > info_all_within_s(traits_num);
		vector < vector <EFF> > eff_Block_s(traits_num);

		// SNP information for each trait within the ith block
		for (size_t j = count; j < summstat_s[0].size(); j++) {
			if (summstat_s[0][j].block == i) {
				for (int k = 0; k < traits_num; k++) {
					info_all_within_s[k].push_back(summstat_s[k][j]);
				}
				count++;
			}else{
				break;
			}
		}

		calcBlock(n_ref, n_obs, v_g, v_e, U_lambda, D_lambda, idv, bed_str, info_all_within_s, num_each_block[i], eff_Block_s);

		for (int m = 0; m < traits_num; m++) {
			for (int l = 0; l < num_each_block[i]; l++) {
				eff_output_s[m].push_back(eff_Block_s[m][l]);
			}
		}
	}
	return 0;
}


// Effect size estimation within each block for model with both large and small effect SNPs
int mtPGSFIT::calcBlock(int n_ref, int n_obs, mat v_g, mat v_e, mat U_lambda, mat D_lambda, vector<int> idv, string bed_str, vector < vector <INFO> > info_all_within_l,
	                    vector < vector <INFO> > info_all_within_s, int num_in_block_l, int num_in_block_s, vector < vector <EFF> > &eff_Block_l, vector < vector <EFF> > &eff_Block_s) {
	SNPPROC cSP;
	IO cIO;
	ifstream bed_in(bed_str.c_str(), ios::binary);

	// import the summary statistics of large effect SNPs within the block
	vector < vector <INFO> > info_within_block_l((int)v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_l; j++) {
			info_within_block_l[i].push_back(info_all_within_l[i][j]);
		}
	}

	// z_score matrix for large effect SNPs
	mat z_matrix_l = zeros<mat>(num_in_block_l, v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int k = 0; k < num_in_block_l; k++)
			z_matrix_l(k, i) = info_within_block_l[i][k].z;
	}

	// reference panel genotype matrix of large effect SNPs
	mat geno_l = zeros<mat>(n_ref, num_in_block_l);
	for (int i = 0; i < num_in_block_l; ++i) {
		vec geno = zeros<vec>(n_ref);
		double maf = 0.0;
		cIO.readSNPIm(info_within_block_l[0][i].pos, n_ref, idv, bed_in, geno, maf);
		cSP.nomalizeVec(geno);
		geno_l.col(i) = geno;
	}

	// import the summary statistics of small effect SNPs within the block
	vector < vector <INFO> > info_within_block_s((int)v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_s; j++) {
			info_within_block_s[i].push_back(info_all_within_s[i][j]);
		}
	}

	// z_score matrix for small effect SNPs
	mat z_matrix_s = zeros<mat>(num_in_block_s, v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int k = 0; k < num_in_block_s; k++)
			z_matrix_s(k, i) = info_within_block_s[i][k].z;
	}

	// reference panel genotype matrix of small effect SNPs
	mat geno_s = zeros<mat>(n_ref, num_in_block_s);
	for (int i = 0; i < num_in_block_s; ++i) {
		vec geno = zeros<vec>(n_ref);
		double maf = 0.0;
		cIO.readSNPIm(info_within_block_s[0][i].pos, n_ref, idv, bed_in, geno, maf);
		cSP.nomalizeVec(geno);
		geno_s.col(i) = geno;
	}

	//estimation
	mat beta_output_l = zeros<mat>(num_in_block_l, (int)v_g.n_rows);
	mat beta_output_s = zeros<mat>(num_in_block_s, (int)v_g.n_rows);
	estBlock(n_ref, n_obs, v_g, v_e, U_lambda, D_lambda, geno_l, geno_s, z_matrix_l, z_matrix_s, beta_output_l, beta_output_s);

	// output effect size for each trait
	// large effect
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_l; j++) {
			EFF eff;
			eff.snp = info_within_block_l[i][j].snp;
			eff.a1 = info_within_block_l[i][j].a1;
			eff.maf = info_within_block_l[i][j].maf;
			eff.beta = beta_output_l(j, i);
			eff_Block_l[i].push_back(eff);
		}
	}

	//small effect
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_s; j++) {
			EFF eff;
			eff.snp = info_within_block_s[i][j].snp;
			eff.a1 = info_within_block_s[i][j].a1;
			eff.maf = info_within_block_s[i][j].maf;
			eff.beta = beta_output_s(j, i);
			eff_Block_s[i].push_back(eff);
		}
	}

	return 0;
}


// Effect size estimation within each block for model with only small effect SNPs
int mtPGSFIT::calcBlock(int n_ref, int n_obs, mat v_g, mat v_e, mat U_lambda, mat D_lambda, vector<int> idv, string bed_str,
	                    vector < vector <INFO> > info_all_within_s, int num_in_block_s, vector < vector <EFF> > &eff_Block_s) {
	SNPPROC cSP;
	IO cIO;
	ifstream bed_in(bed_str.c_str(), ios::binary);

	// import the summary statistics of small effect SNPs within the block
	vector < vector <INFO> > info_within_block_s((int)v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_s; j++) {
			info_within_block_s[i].push_back(info_all_within_s[i][j]);
		}
	}

	// z_score matrix of small effect SNPs
	mat z_matrix_s = zeros<mat>(num_in_block_s, v_g.n_rows);
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int k = 0; k < num_in_block_s; k++)
			z_matrix_s(k, i) = info_within_block_s[i][k].z;
	}

	// reference panel genotype matrix of small effect SNPs
	mat geno_s = zeros<mat>(n_ref, num_in_block_s);
	for (int i = 0; i < num_in_block_s; ++i) {
		vec geno = zeros<vec>(n_ref);
		double maf = 0.0;
		cIO.readSNPIm(info_within_block_s[0][i].pos, n_ref, idv, bed_in, geno, maf);
		cSP.nomalizeVec(geno);
		geno_s.col(i) = geno;
	}

	//estimation
	mat beta_output_s = zeros<mat>(num_in_block_s, (int)v_g.n_rows);
	estBlock(n_ref, n_obs, v_g, v_e, U_lambda, D_lambda, geno_s, z_matrix_s, beta_output_s);

	// output effect size for each trait
	for (int i = 0; i < (int)v_g.n_rows; i++) {
		for (int j = 0; j < num_in_block_s; j++) {
			EFF eff;
			eff.snp = info_within_block_s[i][j].snp;
			eff.a1 = info_within_block_s[i][j].a1;
			eff.maf = info_within_block_s[i][j].maf;
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


// Model fitting for both large and small effect SNPs
int mtPGSFIT::estBlock(int n_ref, int n_obs, mat v_g, mat v_e, mat U_lambda, mat D_lambda, mat geno_l, mat geno_s,
	                    mat z_matrix_l, mat z_matrix_s, mat &beta_output_l, mat &beta_output_s) {

	// LD matrix
	mat SIGMA_ref_l = geno_l.t() * geno_l;
	SIGMA_ref_l /= (double)n_ref;

	mat SIGMA_ref_s = geno_s.t() * geno_s;
	SIGMA_ref_s /= (double)n_ref;

	mat SIGMA_ref_sl = geno_s.t() * geno_l;
	SIGMA_ref_sl /= (double)n_ref;

	// v_g * v_e inverse
	mat vg_ve_inverse = v_g * inv(v_e);
	mat vg_ve_inverse_trans = vg_ve_inverse.t();

	// t(U_lambda) * inv sqrt v_e and its inverse
	mat U_lambda_v_e = U_lambda.t() * inv(sqrtmat_sympd(v_e));
	mat U_lambda_v_e_trans = U_lambda_v_e.t();

	mat U_lambda_v_e_inv = inv(U_lambda_v_e);
	mat U_lambda_v_e_inv_trans = U_lambda_v_e_inv.t();

	// v_g * inv sqrt v_e * U_lambda
	mat vg_ve_U_lambda = v_g * inv(sqrtmat_sympd(v_e)) * U_lambda;
	mat vg_ve_U_lambda_trans = vg_ve_U_lambda.t();

	// large effect SNPs
	mat sigma_inv_z_matrix_l = zeros<mat>((int)z_matrix_l.n_rows, (int)D_lambda.n_rows);

	for (int i = 0; i < (int)D_lambda.n_rows; i++) {
		SIGMA_ref_s.diag() += 1.0 / (D_lambda(i, i) * (double)n_obs);
		vec z_vector_l = z_matrix_l * U_lambda_v_e_trans.col(i);
		vec z_vector_s = z_matrix_s * U_lambda_v_e_trans.col(i);
		vec z_diff_l_s =  z_vector_l - SIGMA_ref_sl.t() * solve(SIGMA_ref_s, z_vector_s);
		mat sigma_s_inv_sigma_sl = PCGm(SIGMA_ref_s, SIGMA_ref_sl, 1000, 1e-7);
		mat sigma_diff = SIGMA_ref_l - SIGMA_ref_sl.t() * sigma_s_inv_sigma_sl;
		SIGMA_ref_s.diag() -= 1.0 / (D_lambda(i, i) * (double)n_obs);
		sigma_inv_z_matrix_l.col(i) = solve(sigma_diff, z_diff_l_s);
	}

	for (int i = 0; i < (int)D_lambda.n_rows; i++) {
		beta_output_l.col(i) = sigma_inv_z_matrix_l * U_lambda_v_e_inv_trans.col(i);
		beta_output_l.col(i) /= sqrt((double)n_obs);
	}

	// small effect SNPs
	mat sigma_sl_beta_l = zeros<mat>((int)SIGMA_ref_sl.n_rows, (int)D_lambda.n_rows);
	for (int i = 0; i < (int)D_lambda.n_rows; i++) {
		sigma_sl_beta_l.col(i) = SIGMA_ref_sl * beta_output_l.col(i);
	}

	mat sigma_s_sigma_s_inv_diff = zeros<mat>((int)SIGMA_ref_s.n_rows, (int)D_lambda.n_rows);
	for (int i = 0; i < (int)D_lambda.n_rows; i++) {
		vec U_lambda_ve_inv_zs_sigma_sl_beta_l_diff = sqrt((double)n_obs) * z_matrix_s * U_lambda_v_e_trans.col(i) - ((double)n_obs) * sigma_sl_beta_l * U_lambda_v_e_trans.col(i);
		SIGMA_ref_s.diag() += 1.0 / (D_lambda(i, i) * (double)n_obs);
		vec sigma_s_inv_diff = solve(SIGMA_ref_s, U_lambda_ve_inv_zs_sigma_sl_beta_l_diff);
		SIGMA_ref_s.diag() -= 1.0 / (D_lambda(i, i) * (double)n_obs);
		sigma_s_sigma_s_inv_diff.col(i) = SIGMA_ref_s * sigma_s_inv_diff;
	}

	for (int i = 0; i < (int)D_lambda.n_rows; i++) {
		vec vg_ve_inv_zs_sigma_sl_beta_l_diff = sqrt((double)n_obs) * z_matrix_s * vg_ve_inverse_trans.col(i) - ((double)n_obs) * sigma_sl_beta_l * vg_ve_inverse_trans.col(i);
		beta_output_s.col(i) = vg_ve_inv_zs_sigma_sl_beta_l_diff - sigma_s_sigma_s_inv_diff * vg_ve_U_lambda_trans.col(i);
	}

	return 0;
}


// Model fitting for only small effect SNPs
int mtPGSFIT::estBlock(int n_ref, int n_obs, mat v_g, mat v_e, mat U_lambda, mat D_lambda, mat geno_s, mat z_matrix_s, mat &beta_output_s) {

	// LD matrix
	mat SIGMA_ref_s = geno_s.t() * geno_s;
	SIGMA_ref_s /= (double)n_ref;

	// v_g * v_e inverse
	mat vg_ve_inverse = v_g * inv(v_e);
	mat vg_ve_inverse_trans = vg_ve_inverse.t();

	// t(U_lambda) * inv sqrt v_e
	mat U_lambda_v_e = U_lambda.t() * inv(sqrtmat_sympd(v_e));
	mat U_lambda_v_e_trans = U_lambda_v_e.t();

	// v_g * inv sqrt v_e * U_lambda
	mat vg_ve_U_lambda = v_g * inv(sqrtmat_sympd(v_e)) * U_lambda;
	mat vg_ve_U_lambda_trans = vg_ve_U_lambda.t();

	// compute each vectors
	mat sigma_sigma_inv_z_matrix = zeros<mat>((int)z_matrix_s.n_rows, (int)D_lambda.n_rows);

	for (int i = 0; i < (int)D_lambda.n_rows; i++) {
		SIGMA_ref_s.diag() += 1.0 / (D_lambda(i, i) * (double)n_obs);
		vec z_vector = z_matrix_s * U_lambda_v_e_trans.col(i);
		vec sigma_inv_z = solve(SIGMA_ref_s, z_vector);
		SIGMA_ref_s.diag() -= 1.0 / (D_lambda(i, i) * (double)n_obs);
		sigma_sigma_inv_z_matrix.col(i) = SIGMA_ref_s * sigma_inv_z;
	}

	//estimation for each trait
	for (int i = 0; i < (int)D_lambda.n_rows; i++) {
		beta_output_s.col(i) = sqrt((double)n_obs) * (z_matrix_s * vg_ve_inverse_trans.col(i) - sigma_sigma_inv_z_matrix * vg_ve_U_lambda_trans.col(i));
	}

	return 0;
}
