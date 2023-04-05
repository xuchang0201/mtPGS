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
#include "mtPGS_model.hpp"
#include "mtPGS.hpp"

using namespace std;
using namespace arma;

mtPGS::mtPGS(void) :
	version("0.0"), date("03/01/2023"), year("2023")
{}

void mtPGS::printHeader(void)
{
	cout << endl;
	cout << "*************************************************************"<< endl;
	cout << "  Multi-trait assisted Polygenic Scores (mtPGS)  " << endl;
	cout << "  Chang Xu, Santhi K. Ganesh, and Xiang Zhou (2023)  " << endl;
	cout << "  Visit http://www.xzlab.org/software.html for update  " << endl;
	cout << "*************************************************************" << endl;
	cout << endl;

	return;
}

void mtPGS::printHelp(void) {
	cout << " Software Input Options:" << endl;
	cout << " --summstat_int      [filename]   " << " specify the GWAS summary statistics from internal dataset with multiple traits measured on the same set of individuals." << endl;
	cout << " --summstat_ext      [filename]   " << " specify the GWAS summary statistics from external datasets (no sample overlap with internal dataset)." << endl;
	cout << " --ref               [filename]   " << " specify the prefix of reference panel." << endl;
	cout << " --n_s               [num]        " << " specify the number of overlapping individuals." << endl;
	cout << " --n_ext             [num]        " << " specify the sample size of each external GWAS summary statistics." << endl;
	cout << " --mafMax            [num]        " << " specify the maximium of the difference between reference panel and summary data." << endl;
	cout << " --block             [filename]   " << " specify the block information." << endl;
	cout << " --output            [filename]   " << " specify the prefix of output files." << endl;
	cout << " --h2                [num]        " << " specify the heritability for each trait." << endl;
	cout << " --gencov            [num]        " << " specify the genetic covariance." << endl;
	cout << " --envcov            [num]        " << " specify the environmental covariance." << endl;
	cout << " --vg                [filename]   " << " specify the directory of genetic variance components file." << endl;
	cout << " --ve                [filename]   " << " specify the directory of environmental variance components file." << endl;
	cout << " --thread            [num]        " << " specify the number of thread." << endl;
	return;
}

void mtPGS::Assign(int argc, char ** argv, PARAM &cPar) {
	
	string str;
	vector <string> summstat_int;
	vector <string> summstat_ext;
	vector <int> n_ext;
	vector <string> output_loc;
	vector <double> h2;
	vector <double> gencov;
	vector <double> envcov;

	for (int i = 0; i < argc; i++) {

		if (strcmp(argv[i], "--summstat_int") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][0] != '-') {
				++i;
				summstat_int.push_back(argv[i]);
			}
			cPar.summstat_int = summstat_int;
		}
		else if (strcmp(argv[i], "--summstat_ext") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][0] != '-') {
				++i;
				summstat_ext.push_back(argv[i]);
			}
			cPar.summstat_ext = summstat_ext;
		}
		else if (strcmp(argv[i], "--ref") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.ref = str;
		}
		else if (strcmp(argv[i], "--n_s") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.n_s = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--n_ext") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][0] != '-') {
				++i;
				n_ext.push_back(atoi(argv[i]));
			}
			cPar.n_ext = n_ext;
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
		else if (strcmp(argv[i], "--output") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][0] != '-') {
				++i;
				output_loc.push_back(argv[i]);
			}
			cPar.output = output_loc;
		}
		else if (strcmp(argv[i], "--h2") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][0] != '-') {
				++i;
				h2.push_back(atof(argv[i]));
			}
			cPar.h2 = h2;
		}
		else if (strcmp(argv[i], "--gencov") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][1] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][1] != '-') {
				++i;
				gencov.push_back(atof(argv[i]));
			}
			cPar.gencov = gencov;
		}
		else if (strcmp(argv[i], "--envcov") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][1] == '-') { continue; }
			while (argv[i + 1] != NULL && argv[i + 1][1] != '-') {
				++i;
				envcov.push_back(atof(argv[i]));
			}
			cPar.envcov = envcov;
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
		else if (strcmp(argv[i], "--thread") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.t = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--target") == 0) {

		    if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
		    ++i;
		    str.clear();
		    str.assign(argv[i]);
		    cPar.target = atoi(str.c_str());
		}
	}
	return;
}

void mtPGS::BatchRun(PARAM &cPar) {

	SNPPROC cSP;
	IO cIO;
	mtPGSFIT cDBSF;

	// number of traits included from internal dataset
	int n_traits = cPar.summstat_int.size();

	// number of traits included from other external datasets
	int n_traits_ext = cPar.summstat_ext.size();

	// number of overlapping individuals
	int n_s = cPar.n_s;
	cout << n_s << " samples were measured for all included phenotypes." << endl;

	// sample sizes of external datasets
	vector <int> n_ext_vec = cPar.n_ext;

	// variance components
	// v_g
	mat v_g;
	v_g.load(cPar.vg);
	v_g.print("matrix g:\n");

	// v_e
	mat v_e;
	v_e.load(cPar.ve);
	v_e.print("matrix e:\n");

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
	cout << "loading block information." << endl;
	cIO.readBlock(block_file, separate, block_dat);
	cout << block_dat.size() << " blocks in the file." << endl;

	// match each summary statistics with reference panel
	// internal data
	vector < vector <POS> > summstat_ref_matched_int(n_traits);

	for (int i = 0; i < n_traits; i++) {
		string trait_name_int = cPar.summstat_int[i];
		vector <SUMM> summ_int;
		int n_snps_int = cIO.readSumm(trait_name_int, separate, summ_int);
		cout << n_snps_int << " SNPs included through readSumm from internal dataset for trait " << i + 1 << "." << endl;
		vector <POS> inter_int;
		cSP.matchRef(summ_int, ref_bim, inter_int, mafMax);
		cout << "After filtering, " << inter_int.size() << " SNPs of trait " << i + 1 << " from internal dataset are included." << endl;
		summstat_ref_matched_int[i] = inter_int;
	}

	// external data
	vector < vector <POS> > summstat_ref_matched_ext(n_traits_ext);

	for (int i = 0; i < n_traits_ext; i++) {
		string trait_name_ext = cPar.summstat_ext[i];
		vector <SUMM> summ_ext;
		int n_snps_ext = cIO.readSumm(trait_name_ext, separate, summ_ext);
		cout << n_snps_ext << " SNPs included through readSumm from other external datasets for trait " << i + 1 << "." << endl;
		vector <POS> inter_ext;
		cSP.matchRef(summ_ext, ref_bim, inter_ext, mafMax);
		cout << "After filtering, " << inter_ext.size() << " SNPs of trait " << i + 1 << " from other external datasets are included." << endl;
		summstat_ref_matched_ext[i] = inter_ext;
	}

	// extract overlapping SNPs across the internal and external datasets
	int n_summ_stat = n_traits + n_traits_ext;
	vector < vector <string> > snps_each_trait_input(n_summ_stat);

	// load internal data
	for (int i = 0; i < n_traits; i++) {
		for (size_t j = 0; j < summstat_ref_matched_int[i].size(); ++j) {
			snps_each_trait_input[i].push_back(summstat_ref_matched_int[i][j].snp);
		}
	}

	// load external data
	for (int k = 0; k < n_traits_ext; k++) {
		int ext_index = k + n_traits;
		for (size_t j = 0; j < summstat_ref_matched_ext[k].size(); ++j) {
			snps_each_trait_input[ext_index].push_back(summstat_ref_matched_ext[k][j].snp);
		}
	}

	// find overlapping SNPs
	vector < string > overlapping_snps;
	cSP.snps_intersect(snps_each_trait_input, overlapping_snps);
	
	// obtain the corresponding POS information based on overlapping SNPs across traits
	// internal
	vector < vector < POS > > overlapped_POS_input_block_int(n_traits);

	for (int i = 0; i < n_traits; i++) {
		map <string, POS> pos_map_int;
		for (size_t j = 0; j < summstat_ref_matched_int[i].size(); ++j) {
			if (pos_map_int.find(summstat_ref_matched_int[i][j].snp) == pos_map_int.end()) {
				pos_map_int[summstat_ref_matched_int[i][j].snp] = summstat_ref_matched_int[i][j];
			}
		}

		for (size_t k = 0; k < overlapping_snps.size(); ++k) {
			if (pos_map_int.find(overlapping_snps[k]) != pos_map_int.end()) {
				overlapped_POS_input_block_int[i].push_back(pos_map_int[overlapping_snps[k]]);
			}
		}
	}

	// external
	vector < vector < POS > > overlapped_POS_input_block_ext(n_traits_ext);

	for (int i = 0; i < n_traits_ext; i++) {
		map <string, POS> pos_map_ext;
		for (size_t j = 0; j < summstat_ref_matched_ext[i].size(); ++j) {
			if (pos_map_ext.find(summstat_ref_matched_ext[i][j].snp) == pos_map_ext.end()) {
				pos_map_ext[summstat_ref_matched_ext[i][j].snp] = summstat_ref_matched_ext[i][j];
			}
		}

		for (size_t k = 0; k < overlapping_snps.size(); ++k) {
			if (pos_map_ext.find(overlapping_snps[k]) != pos_map_ext.end()) {
				overlapped_POS_input_block_ext[i].push_back(pos_map_ext[overlapping_snps[k]]);
			}
		}
	}

	// match with block file
	// internal
	vector < vector <INFO> > summstat_all_traits_int(n_traits);
	int num_block_int = 0; //initialize the number of block

	for (int i = 0; i < n_traits; i++) {
		vector <INFO> info_int;
		num_block_int = cSP.addBlock(overlapped_POS_input_block_int[i], block_dat, info_int);
		cout << num_block_int << " blocks for trait " << i + 1 << " of internal dataset." << endl;
		cout << info_int.size() << " for info_" << i + 1 << " of internal dataset." << endl;
		summstat_all_traits_int[i] = info_int;
	}

	// external
	vector < vector <INFO> > summstat_all_traits_ext(n_traits_ext);
	int num_block_ext = 0; //initialize the number of block

	for (int i = 0; i < n_traits_ext; i++) {
		vector <INFO> info_ext;
		num_block_ext = cSP.addBlock(overlapped_POS_input_block_ext[i], block_dat, info_ext);
		cout << num_block_ext << " blocks for trait " << i + 1 << " of other external datasets." << endl;
		cout << info_ext.size() << " for info_" << i + 1 << " of other external datasets." << endl;
		summstat_all_traits_ext[i] = info_ext;
	}


	// (1/n_snp) * v_g
	int number_snps = overlapping_snps.size();
	mat v_g_per_snp = v_g;
	v_g_per_snp /= (double)number_snps;
	

	// detect large effect SNPs
	cout << "Performing large effect SNPs detection for target trait..." << endl;
	// define the number of input SNPs
	int num_snp_in = summstat_all_traits_int[0].size();
	int target_index = cPar.target;
	cout << "The data for target trait are " << cPar.summstat_int[target_index] << " and " << cPar.summstat_ext[target_index] << "." << endl;

	// prepare METAL input for meta-analysis
	string METAL_output_file_name = cPar.output[0] + "_" + cPar.output[1] + "_meta_analysis_output .tbl";
	string METAL_input = cPar.output[0] + "_" + cPar.output[1] + "_METAL_script.txt";
	ofstream meta(METAL_input.c_str());
	meta << "SEPARATOR" << " " << "TAB" << endl;
	meta << "SCHEME" << " " << "STDERR" << endl;
	meta << "STDERR" << " " << "se" << endl;
	meta << "MARKER" << " " << "rs" << endl;
	meta << "WEIGHT" << " " << "n_obs" << endl;
	meta << "ALLELE" << " " << "allele1" << " " << "allele0" << endl;
	meta << "FREQ" << " " << "af" << endl;
	meta << "EFFECT" << " " << "beta" << endl;
	meta << "STDERR" << " " << "se" << endl;
	meta << "PVAL" << " " << "p_wald" << endl;
	meta << "PROCESS" << " " << cPar.summstat_int[target_index] << endl;
	meta << "PROCESS" << " " << cPar.summstat_ext[target_index] << endl;
	meta << "OUTFILE" << " " << METAL_output_file_name << endl;
	meta << "ANALYZE" << endl;
	meta.close();

	// Perform meta-analysis
	string METAL_command = "metal ";
	METAL_command = METAL_command + METAL_input;
	const char *meta_analysis_command = METAL_command.c_str();
	cout << "Performing meta-analysis using " << meta_analysis_command << endl;
	int METAL_execute = system(meta_analysis_command);

	// Perform LD clumping using PLINK
	string plink_str = "plink-1.9 --bfile ";
	string clump_output = cPar.output[0] + "_" + cPar.output[1] + "_clumping_result";
	string meta_analysis_data = cPar.output[0] + "_" + cPar.output[1] + "_meta_analysis_output1.tbl";
	plink_str = plink_str + ref_prefix + " --silent --clump " + meta_analysis_data + " --clump-snp-field MarkerName --clump-field P-value --clump-r2 0.2 --clump-p1 1e-6 --clump-kb 1000 --out " + clump_output;
	const char *command = plink_str.c_str();
	cout << "Performing LD clumping using " << command << endl;
	int ld_command_execute = system(command);

	// Extract large effect SNPs
	string large_list_name = cPar.output[0] + "_" + cPar.output[1] + "_l_snp_list.txt";
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

		cout << large_snp_list.size() << " large effect SNPs loaded from the PLINK clumping result." << endl;
	}

	// load large effect and small effect SNPs
	// get the list of rsid for SNPs included in the analysis
	vector <string> rsid_list_all;
	for (int i = 0; i < (int)summstat_all_traits_int[0].size(); i++) {
		rsid_list_all.push_back(summstat_all_traits_int[0][i].snp);
	}

	cout << rsid_list_all.size() << " total SNPs included in the rsid list." << endl;

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
	cout << large_snp_index.size() << " total index from the list." << endl;

	// load internal data
	cout << "Loading INFO from internal dataset..." << endl;
	// large  
	vector < vector <INFO> > summstat_all_traits_l_int(n_traits);
	for (int i = 0; i < n_traits; i++) {
		for (int j = 0; j < (int)large_snp_index.size(); j++) {
			summstat_all_traits_l_int[i].push_back(summstat_all_traits_int[i][large_snp_index[j]]);
		}
		cout << summstat_all_traits_l_int[i].size() << " from internal dataset for info_l_" << i + 1 << "." << endl;
	}

	// small 
	vector < vector <INFO> > summstat_all_traits_s_int(n_traits);
	for (int i = 0; i < n_traits; i++) {
		for (int j = 0; j < (int)summstat_all_traits_int[0].size(); j++) {
			auto it = find(large_snp_index.begin(), large_snp_index.end(), j);
			if (it == large_snp_index.end()) {
				summstat_all_traits_s_int[i].push_back(summstat_all_traits_int[i][j]);
			}
		}
		cout << summstat_all_traits_s_int[i].size() << " from internal dataset for info_s_" << i + 1 << "." << endl;
	}

	// load external data
	cout << "Loading INFO from other external datasets..." << endl;
	// large  
	vector < vector <INFO> > summstat_all_traits_l_ext(n_traits);
	for (int i = 0; i < n_traits; i++) {
		for (int j = 0; j < (int)large_snp_index.size(); j++) {
			summstat_all_traits_l_ext[i].push_back(summstat_all_traits_ext[i][large_snp_index[j]]);
		}
		cout << summstat_all_traits_l_ext[i].size() << " from other external datasets for info_l_" << i + 1 << "." << endl;
	}

	// small 
	vector < vector <INFO> > summstat_all_traits_s_ext(n_traits);
	for (int i = 0; i < n_traits; i++) {
		for (int j = 0; j < (int)summstat_all_traits_ext[0].size(); j++) {
			auto it = find(large_snp_index.begin(), large_snp_index.end(), j);
			if (it == large_snp_index.end()) {
				summstat_all_traits_s_ext[i].push_back(summstat_all_traits_ext[i][j]);
			}
		}
		cout << summstat_all_traits_s_ext[i].size() << " from other external datasets for info_s_" << i + 1 << "." << endl;
	}

	// fit model with both large and small effect SNPs
	vector < vector <EFF> > effect_output_l(n_traits);
	vector < vector <EFF> > effect_output_s(n_traits);
	vector <int> idv(n_ref);
	for (int i = 0; i < n_ref; i++) idv[i] = 1;
	cout << "Fitting model..." << endl;
	cDBSF.est_w_ext(n_ref, n_s, n_ext_vec, num_block_int, v_e, v_g_per_snp, idv, bed_str, summstat_all_traits_l_int, summstat_all_traits_s_int,
		            summstat_all_traits_l_ext, summstat_all_traits_s_ext, effect_output_l, effect_output_s);


	// output results for large and small effect model
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

	return;
}
