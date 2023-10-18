/*
Multi-trait assisted Polygenic Scores (mtPGS)
*/

#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include <math.h>
#include <armadillo>
#include <algorithm>

#include "mtPGS_function.hpp"
#include "mtPGS_model_int_only.hpp"
#include "mtPGS_int_only.hpp"

using namespace std;
using namespace arma;

mtPGS::mtPGS(void) :
	version("0.0"), date("09/09/2023"), year("2023")
{}

void mtPGS::printHeader(void)
{
	cout << endl;
	cout << "*************************************************************" << endl;
	cout << "  Multi-trait assisted Polygenic Scores (mtPGS)  " << endl;
	cout << "  Chang Xu, Santhi K. Ganesh, and Xiang Zhou (2023)  " << endl;
	cout << "  Visit http://www.xzlab.org/software.html for update  " << endl;
	cout << "*************************************************************" << endl;
	cout << endl;

	return;
}

void mtPGS::printHelp(void) {
	cout << " Software input options for handling two sets of summary statistics:" << endl;
	cout << " --summstat   [filename]   " << " specify the input GWAS summary statistics for each trait." << endl;
	cout << " --ref        [filename]   " << " specify the prefix of reference panel." << endl;
	cout << " --n          [num]        " << " specify the sample size of the summary statistics for each trait." << endl;
	cout << " --mafMax     [num]        " << " specify the maximium of the difference between reference panel and summary data." << endl;
	cout << " --block      [filename]   " << " specify the block information." << endl;
	cout << " --plink      [filename]   " << " specify the directory for the PLINK 1.9 software." << endl;
	cout << " --c_t        [filename]   " << " specify the prefix for the outputs of PLINK C+T procedure." << endl;
	cout << " --output     [filename]   " << " specify the prefix of output files." << endl;
	cout << " --vg         [filename]   " << " specify the directory of genetic variance components file." << endl;
	cout << " --ve         [filename]   " << " specify the directory of environmental variance components file." << endl;
	cout << " --target     [num]        " << " specify the index for target trait." << endl;
	cout << " --r2         [num]        " << " specify the LD threshold for the C+T procedure." << endl;
	cout << " --pval       [num]        " << " specify the p value threshold for the C+T procedure." << endl;
	return;
}

void mtPGS::Assign(int argc, char ** argv, PARAM &cPar) {
	
	string str;
	vector <string> summstat;
	vector <int> sample_size;
	vector <string> output_loc;

	for (int i = 0; i < argc; i++) {

		if (strcmp(argv[i], "--summstat") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][0] != '-') {
				++i;
				summstat.push_back(argv[i]);
			}
			cPar.summstat = summstat;
		}
		else if (strcmp(argv[i], "--target") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.target = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--ref") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.ref = str;
		}
		else if (strcmp(argv[i], "--n") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][0] != '-') {
				++i;
				sample_size.push_back(atoi(argv[i]));
			}
			cPar.n = sample_size;
		}
		else if (strcmp(argv[i], "--mafMax") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.mafMax = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--block") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.block = str;
		}
		else if (strcmp(argv[i], "--plink") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.plink = str;
		}
		else if (strcmp(argv[i], "--c_t") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.c_t = str;
		}
		else if (strcmp(argv[i], "--output") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][0] != '-') {
				++i;
				output_loc.push_back(argv[i]);
			}
			cPar.output = output_loc;
		}
		else if (strcmp(argv[i], "--vg") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.vg = str;
		}
		else if (strcmp(argv[i], "--ve") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.ve = str;
		}
		else if (strcmp(argv[i], "--r2") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.r2 = str;
		}
		else if (strcmp(argv[i], "--pval") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.pval = str;
		}
	}
	return;
}

void mtPGS::BatchRun(PARAM &cPar) {

	SNPPROC cSP;
	IO cIO;
	mtPGSFIT cDBSF;

	cout << "Loading input datasets..." << endl;

	// number of traits included
	int n_traits = cPar.n.size();

	// sample size of summary statistics
	int n_obs = *std::min_element(cPar.n.begin(), cPar.n.end());

	// variance components
	// v_g
	mat v_g;
	v_g.load(cPar.vg);
	v_g.print("Genetic covariance matrix:\n");

	// v_e
	mat v_e;
	v_e.load(cPar.ve);
	v_e.print("Environmental covariance matrix:\n");

	// sample size of reference panel 
	char separate[] = "\t";
	string ref_prefix = cPar.ref;
	string ref_fam_name = ref_prefix + ".fam";
	int n_ref = cIO.getRow(ref_fam_name);
	cout << n_ref << " individuals included in the reference panel." << endl;

	// assign bed file name of reference panel
	string bed_str = ref_prefix + ".bed";

	// number of SNPs in the reference panel
	double mafMax = cPar.mafMax;
	map <string, ALLELE> ref_bim;
	bool constr = true;
	if (abs(mafMax - 1.0) < 1e-10) {
		constr = false;
	}
	cIO.readBim(n_ref, ref_prefix, separate, ref_bim, constr);
	int num_snp_ref = ref_bim.size();
	cout << num_snp_ref << " SNPs included in the reference panel." << endl;

	// input block file
	string block_file = cPar.block;
	vector <BLOCK> block_dat;
	cout << "Loading block information..." << endl;
	cIO.readBlock(block_file, separate, block_dat);
	cout << block_dat.size() << " blocks in the file." << endl;

	// match each summary statistics with reference panel
	vector < vector <POS> > summstat_ref_matched(n_traits);

	for (int i = 0; i < n_traits; i++) {
		string trait_name = cPar.summstat[i];
		vector <SUMM> summ;
		int n_snps = cIO.readSumm(trait_name, separate, summ);
		cout << n_snps << " SNPs included from summary statistics for trait " << i + 1 << "." << endl;
		vector <POS> inter;
		cSP.matchRef(summ, ref_bim, inter, mafMax);
		cout << "After filtering, " << inter.size() << " SNPs of trait " << i + 1 << " are included." << endl;
		summstat_ref_matched[i] = inter;
	}

	// extract overlapping SNPs across the traits
	vector < vector <string> > snps_each_trait_input(n_traits);

	for (int i = 0; i < n_traits; i++) {
		for (size_t j = 0; j < summstat_ref_matched[i].size(); ++j) {
			snps_each_trait_input[i].push_back(summstat_ref_matched[i][j].snp);
		}
	}

	vector < string > overlapping_snps;
	cSP.snps_intersect(snps_each_trait_input, overlapping_snps);
	
	// obtain the corresponding POS information based on overlapping SNPs across traits
	vector < vector < POS > > overlapped_POS_input_block(n_traits);

	for (int i = 0; i < n_traits; i++) {
		map <string, POS> pos_map;
		for (size_t j = 0; j < summstat_ref_matched[i].size(); ++j) {
			if (pos_map.find(summstat_ref_matched[i][j].snp) == pos_map.end()) {
				pos_map[summstat_ref_matched[i][j].snp] = summstat_ref_matched[i][j];
			}
		}

		for (size_t k = 0; k < overlapping_snps.size(); ++k) {
			if (pos_map.find(overlapping_snps[k]) != pos_map.end()) {
				overlapped_POS_input_block[i].push_back(pos_map[overlapping_snps[k]]);
			}
		}
	}

	// match with block file
	vector < vector <INFO> > summstat_all_traits(n_traits);
	int num_block = 0; //initialize the number of block

	for (int i = 0; i < n_traits; i++) {
		vector <INFO> info;
		num_block = cSP.addBlock(overlapped_POS_input_block[i], block_dat, info);
		summstat_all_traits[i] = info;
	}


	// (1/n_snp) * v_g
	int number_snps = summstat_all_traits[0].size();
	mat v_g_per_snp = v_g;
	v_g_per_snp /= (double)number_snps;
	

	// detect large effect SNPs
	cout << "Performing large effect SNPs detection for target trait..." << endl;
	// define the number of input SNPs
	int num_snp_in = summstat_all_traits[0].size();
	int target_index = cPar.target;
	cout << "The data for target trait are " << cPar.summstat[target_index] << "." << endl;
	
	// define a vector for P-values
	vector <double> P_min_vector;

	// extract P-values for target trait
	for (int i = 0; i < num_snp_in; i++) {
		P_min_vector.push_back(summstat_all_traits[target_index][i].P);
	}

	// prepare output for C+T procedure
	string LRTout = cPar.c_t + "_LRT_output.txt";
	ofstream lrt(LRTout.c_str());
	lrt << "SNP" << " " << "P" << endl;
	for (int j = 0; j < num_snp_in; j++) {
		lrt << summstat_all_traits[target_index][j].snp << " " << P_min_vector[j] << endl;
	}
	lrt.close();

	// Perform LD clumping using PLINK
	string plink_str = cPar.plink + " --bfile ";
	string clump_output = cPar.c_t + "_clumping_result";
	plink_str = plink_str + ref_prefix + " --silent --clump " + LRTout + " --clump-r2 " + cPar.r2 + " --clump-p1 " + cPar.pval + " --clump-kb 1000 --out " + clump_output;
	const char *command = plink_str.c_str();
	cout << "Performing LD clumping using command: " << command << endl;
	int ld_command_execute = system(command);

	// Extract large effect SNPs
	string large_list_name = cPar.c_t + "_l_snp_list.txt";
	string large_list = "awk '{print $3}' " + clump_output + ".clumped > " + large_list_name;
	const char *awk_command = large_list.c_str();
	int awk_command_execute = system(awk_command);

	// Read the file of large effect SNPs list
	vector <string> large_snp_list;
	ifstream large_in(large_list_name.c_str());

	if (!large_in) {
		cerr << "ERROR: Cannot open the File : " << large_list_name << endl;
	}
	else {
		string large_str;
		while (getline(large_in, large_str)) {
			if (large_str.size() > 0 && large_str != "SNP")
				large_snp_list.push_back(large_str);
		}
		large_in.close();

		cout << "A total of " << large_snp_list.size() << " large effect SNPs loaded from the PLINK clumping result." << endl;
	}

	// If large effect SNPs detected, proceed with model of both large and small effects
	if (large_snp_list.size() != 0) {
		// load large effect and small effect SNPs
	    // get the list of rsid for SNPs included in the analysis
		vector <string> rsid_list_all;
		for (int i = 0; i < (int)summstat_all_traits[0].size(); i++) {
			rsid_list_all.push_back(summstat_all_traits[0][i].snp);
		}

		// get the index of large effect SNPs
		vector <int> large_snp_index;
		for (int i = 0; i < (int)large_snp_list.size(); i++) {
			auto it = find(rsid_list_all.begin(), rsid_list_all.end(), large_snp_list[i]);
			if (it != rsid_list_all.end()) {
				int index = it - rsid_list_all.begin();
				large_snp_index.push_back(index);
			}
		}

		// sort the index
		sort(large_snp_index.begin(), large_snp_index.end());

		// large
		vector < vector <INFO> > summstat_all_traits_l(n_traits);
		for (int i = 0; i < n_traits; i++) {
			for (int j = 0; j < (int)large_snp_index.size(); j++) {
				summstat_all_traits_l[i].push_back(summstat_all_traits[i][large_snp_index[j]]);
			}
		}

		// small
		vector < vector <INFO> > summstat_all_traits_s(n_traits);
		for (int i = 0; i < n_traits; i++) {
			for (int j = 0; j < (int)summstat_all_traits[0].size(); j++) {
				auto it = find(large_snp_index.begin(), large_snp_index.end(), j);
				if (it == large_snp_index.end()) {
					summstat_all_traits_s[i].push_back(summstat_all_traits[i][j]);
				}
			}
		}

		// fit model with both large and small effect SNPs
		vector < vector <EFF> > effect_output_l(n_traits);
		vector < vector <EFF> > effect_output_s(n_traits);
		vector <int> idv(n_ref);
		for (int i = 0; i < n_ref; i++) idv[i] = 1;
		cout << "Fitting model..." << endl;
		cDBSF.est(n_ref, n_obs, num_block, v_e, v_g_per_snp, idv, bed_str, summstat_all_traits_l, summstat_all_traits_s, effect_output_l, effect_output_s);

		//output results for large and small effect model
		for (int k = 0; k < n_traits; k++) {
			string output_prefix = cPar.output[k];
			string eff_str = output_prefix + ".txt";
			ofstream effFout(eff_str.c_str());

			for (size_t i = 0; i < effect_output_l[k].size(); ++i) {
				double beta_noscl_l = effect_output_l[k][i].beta / sqrt(2 * effect_output_l[k][i].maf * (1 - effect_output_l[k][i].maf));
				if (effect_output_l[k][i].snp != "rs" && isinf(beta_noscl_l) == false)
					effFout << effect_output_l[k][i].snp << " " << effect_output_l[k][i].a1 << " " << effect_output_l[k][i].beta << " " << beta_noscl_l << " " << 1 << endl;
			}

			for (size_t i = 0; i < effect_output_s[k].size(); ++i) {
				double beta_noscl_s = effect_output_s[k][i].beta / sqrt(2 * effect_output_s[k][i].maf * (1 - effect_output_s[k][i].maf));
				if (effect_output_s[k][i].snp != "rs" && isinf(beta_noscl_s) == false)
					effFout << effect_output_s[k][i].snp << " " << effect_output_s[k][i].a1 << " " << effect_output_s[k][i].beta << " " << beta_noscl_s << " " << 0 << endl;
			}
			effFout.close();
		}
	}

	// If no large effect SNP detected, proceed with model of small effect only 
	if (large_snp_list.size() == 0) {
		// fit model with only small effect SNPs
		vector < vector <EFF> > effect_output_s(n_traits);
		vector<int> idv(n_ref);
		for (int i = 0; i < n_ref; i++) idv[i] = 1;
		cout << "Fitting model..." << endl;
		cDBSF.est(n_ref, n_obs, num_block, v_e, v_g_per_snp, idv, bed_str, summstat_all_traits, effect_output_s);

		//output results for small effect model
		for (int k = 0; k < n_traits; k++) {
			string output_prefix = cPar.output[k];
			string eff_str = output_prefix + ".txt";
			ofstream effFout(eff_str.c_str());
			for (size_t i = 0; i < effect_output_s[k].size(); ++i) {
				double beta_noscl_s = effect_output_s[k][i].beta / sqrt(2 * effect_output_s[k][i].maf * (1 - effect_output_s[k][i].maf));
				if (effect_output_s[k][i].snp != "rs" && isinf(beta_noscl_s) == false)
					effFout << effect_output_s[k][i].snp << " " << effect_output_s[k][i].a1 << " " << effect_output_s[k][i].beta << " " << beta_noscl_s << " " << 0 << endl;
			}
			effFout.close();
		}
	}

	return;
}
