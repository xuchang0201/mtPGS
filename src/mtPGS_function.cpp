/*
Multi-trait assisted Polygenic Scores (mtPGS)
*/

#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <bitset>
#include <numeric>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <cctype>
#include <boost/math/distributions/students_t.hpp>

#include <armadillo>

#include "mtPGS_function.hpp"

using namespace std;
using namespace arma;
// input block information
int IO::readBlock(string infile, char *separator, vector <BLOCK> &block){

	int block_num = getRow(infile);
	block.resize(block_num);
	string oneblock, element;
	ifstream file_stream(infile.c_str());
	int block_count = 0; 
	while (getline(file_stream, oneblock)) {
		vector<string> block_str;
		stringstream oneblock_stream(oneblock);
		while (getline(oneblock_stream, element, *separator))
			block_str.push_back(element);
		block[block_count].chr = block_str[0]; 
		block[block_count].start = atol(block_str[1].c_str());
		block[block_count].end = atol(block_str[2].c_str()); 
		block_count++;
	}
	file_stream.close();
	return 0;
}

// get row number
int IO::getRow(string infile){
	
	int n_row = 0;
	string row;
	ifstream file_stream(infile.c_str());
	while (getline(file_stream, row))
		n_row++;
	file_stream.close();
	return n_row;
}

// input bim data
int IO::readBim(int n_ref, string ref_str, char *separator, map<string, ALLELE> &bim, bool constr){
	
	string bed_str = ref_str + ".bed", bim_str = ref_str + ".bim";
	ifstream bim_stream(bim_str.c_str());
	int n_snp = getRow(bim_str);
	vec maf = zeros<vec>(n_snp); 
	if (constr == true){
		ifstream bed_stream(bed_str);
		cout << "Calculating MAF of reference panel ..." << endl;
		vector<int> idv(n_ref); 
		for (int i = 0; i < n_ref; i++) idv[i] = 1; 
		for (int i = 0; i < n_snp; i++){
			vec geno = zeros<vec>(n_ref); 
			readSNPIm(i, n_ref, idv, bed_stream, geno, maf(i));
		}
		bed_stream.close();
	} else {
		cout << "[WARNING] Do not consider the difference between reference panel and summary data ..." << endl;
	}
	
	string snp, elem;
	int count = 0;
	while (getline(bim_stream, snp)) {
		vector<string> snp_vec;
		ALLELE allele; 
		stringstream snp_stream(snp);
		while (getline(snp_stream, elem, *separator)) snp_vec.push_back(elem);
		allele.pos = count; 
		allele.a1 = snp_vec[4];
		allele.a2 = snp_vec[5];
		allele.maf = maf(count);
		bim.insert(pair<string, ALLELE>(snp_vec[1], allele));
		count++; 
	}
	bim_stream.close();
	return 0;
}

// calculate P value from GEMMA output
double IO::calP(double beta, double se, int sampleSize){
	double t = beta / se;
	boost::math::students_t tDis(sampleSize);
	double p = 2 * cdf(complement(tDis, fabs(t)));
	return p;
}

// input summary data
int IO::readSumm(string summ_str, char *separator, vector<SUMM > &summ){

	string snp, element, chr_str;
	ifstream summ_stream(summ_str.c_str());

	int num = getRow(summ_str);
	// read the summary data of specific chromosome
	summ.resize(num);
	int summ_count = 0; 
	while (getline(summ_stream, snp)){
		vector<string> snp_summ;
		stringstream summ_stream(snp);
		while (getline(summ_stream, element, *separator)) 
			snp_summ.push_back(element);
		double P = atof(snp_summ[10].c_str()); 
		double se, sei;
		if (isdigit(snp_summ[9].c_str()[0])){
			se = atof(snp_summ[9].c_str());
			if (se - 0.0 > 1e-20){
				summ[summ_count].z = atof(snp_summ[8].c_str()) / se;
				summ[summ_count].P = calP(atof(snp_summ[8].c_str()), se, atoi(snp_summ[4].c_str()));
			} else {
				summ[summ_count].z = 0; 
				summ[summ_count].P = 1;
			}
		} else {
			summ[summ_count].z = 0; 
			summ[summ_count].P = 1;
		}
		summ[summ_count].chr = atoi(snp_summ[0].c_str());
		summ[summ_count].snp = snp_summ[1];
		summ[summ_count].ps = atol(snp_summ[2].c_str());
		summ[summ_count].a1 = snp_summ[5].c_str();
		summ[summ_count].a2 = snp_summ[6].c_str();
		double af = atof(snp_summ[7].c_str());
		summ[summ_count].maf = min(af, 1.0 - af);
		summ_count++; 
	
	}
	summ_stream.close();
	
	return summ_count;
}

// input genotype data
// modify from GEMMA, Xiang Zhou et al.
void IO::readSNPIm(const int pos, int ni_test, const vector<int> &indicator_idv, ifstream &infile, vec &geno, double &maf) {

	// debug_msg("entered");
	size_t ni_total = indicator_idv.size(), n_bit;
	if (ni_total % 4 == 0) {
		n_bit = ni_total / 4;
	}
	else {
		n_bit = ni_total / 4 + 1;
	}

	// n_bit, and 3 is the number of magic numbers.
	infile.seekg(pos * n_bit + 3);

	// Read genotypes.
	char ch[1];
	bitset<8> b;
	
	double geno_mean = 0.0;
	size_t c = 0, c_idv = 0;
	vector<size_t> geno_miss;
	int freq[3];
	freq[0] = freq[1] = freq[2] = 0;
	for (size_t i = 0; i < n_bit; ++i) {
		infile.read(ch, 1);
		b = ch[0];

		// Minor allele homozygous: 2.0; major: 0.0.
		for (size_t j = 0; j < 4; ++j) {
			if ((i == (n_bit - 1)) && c == ni_total) {
				break;
			}
			if (indicator_idv[c] == 0) {
				c++;
				continue;
			}
			c++;

			if (b[2 * j] == 0) {
				if (b[2 * j + 1] == 0) {
					geno(c_idv) = 2.0;
					geno_mean += 2.0;
					freq[2]++;
				}
				else {
					geno(c_idv) = 1.0;
					geno_mean += 1.0;
					freq[1]++;
				}
			}
			else {
				if (b[2 * j + 1] == 1) {
					geno(c_idv) = 0.0;
					geno_mean += 0.0;
					freq[0]++; 
				}
				else {
					geno_miss.push_back(c_idv);
				}
			}

			c_idv++;
		}
	}
	// max imputation
	// int imp_val = distance(freq, max_element(freq, freq + 3));
	// mean imputation
	geno_mean /= (double)(c_idv - geno_miss.size());
	for (size_t i = 0; i < geno_miss.size(); ++i) 
		geno(geno_miss[i]) = geno_mean;
	double af = 0.5 * sum(geno) / geno.n_elem; 
	maf = min(af, 1.0 - af);
	return ;
}

double IO::getWalltime(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void SNPPROC::nomalizeVec(vec &x) {
	
	x -= mean(x); 
	x /= stddev(x);
	return;
}

// match summary result and reference panel
int SNPPROC::matchRef(vector<SUMM> summ, map<string, ALLELE> bim, vector<POS> &inter_snp, double mafMax) {
	
	int dis_count = 0, maf_count = 0; 
	for (size_t i = 0; i < summ.size(); ++i) {
		bool a1_bool = bim[summ[i].snp].a1 == summ[i].a1; 
		bool a2_bool = bim[summ[i].snp].a2 == summ[i].a2; 
		bool maf_bool = fabs(bim[summ[i].snp].maf - summ[i].maf) < mafMax; 
		if (a1_bool == false || a2_bool == false) dis_count++; 
		if (maf_bool == false) maf_count++;
		if (a1_bool == true && a2_bool == true && maf_bool == true && bim.find(summ[i].snp) != bim.end()){
			POS s_summ; 
			s_summ.snp = summ[i].snp; 
			s_summ.ps = summ[i].ps; 
			s_summ.pos = bim[summ[i].snp].pos; 
			s_summ.a1 = summ[i].a1; 
			s_summ.a2 = summ[i].a2;
			s_summ.maf = summ[i].maf; 
			s_summ.z = summ[i].z; 
			s_summ.P = summ[i].P;
			inter_snp.push_back(s_summ);
		}
	}
	cout << "Number of allele discrepency: " << dis_count << endl;
	cout << "Number of maf discrepency:    " << maf_count << endl;
	return 0;
}

int SNPPROC::addBlock(vector<POS> summ, vector <BLOCK> block, vector<INFO> &pos_block){
	
	int num_block = block.size();
	int count = 0; 
	pos_block.resize(summ.size());
	for (int i = 0; i < num_block; i++){
		int start = block[i].start;
		int end = block[i].end;
		// cout << i << " " << begin << " " << end << endl;
		for (size_t j = count; j < summ.size(); j++){
			if (summ[j].ps >= start && summ[j].ps < end) { 
				pos_block[j].snp = summ[j].snp;
				pos_block[j].block = i; 
				pos_block[j].ps = summ[j].ps; // snp position
				pos_block[j].pos = summ[j].pos; // bim file position
				pos_block[j].a1 = summ[j].a1;
				pos_block[j].maf = summ[j].maf;
				pos_block[j].P = summ[j].P;
				pos_block[j].z = summ[j].z;
				count++;
			}else{
				break;
			}
		}
	}
	return num_block; 
}

// extract overlapping SNPs across traits
int SNPPROC::snps_intersect(vector < vector <string> > snps_multi_traits, vector <string> &overlap_snps) {

	for (size_t i = 1; i < snps_multi_traits.size(); ++i) {
		// use the first vector as index, sort the remaining vectors
		sort(snps_multi_traits[i].begin(), snps_multi_traits[i].end());
	}

	for (size_t i = 0; i < snps_multi_traits[0].size(); ++i) {
		bool snp_found = true;
		string snp_i = snps_multi_traits[0][i];

		for (size_t j = 1; j < snps_multi_traits.size(); ++j) {
			if (binary_search(snps_multi_traits[j].begin(), snps_multi_traits[j].end(), snp_i)) {
				snp_found = true;
			}
			else {
				snp_found = false;
				break;
			}
		}

		if (snp_found) {
			overlap_snps.push_back(snp_i);
		}
	}

	return 0;
}

bool sortP(const INFO &snpinfo1, const INFO &snpinfo2){
	return snpinfo1.P < snpinfo2.P;
}

bool sortPOS(const INFO &snpinfo1, const INFO &snpinfo2){
	return snpinfo1.pos < snpinfo2.pos;
}

bool sortPS(const SUMMP &summp1, const SUMMP &summp2)
{
	return summp1.ps < summp2.ps;
}


void printProgBar(int percent) {
	string bar;
	for (int i = 0; i < 50; i++) {
		if (i < (percent / 2)) {
			bar.replace(i, 1, "=");
		}
		else if (i == (percent / 2)) {
			bar.replace(i, 1, ">");
		}
		else {
			bar.replace(i, 1, " ");
		}
	}

	cout << "\r" "[" << bar << "] ";
	cout.width(3);
	cout << percent << "%     " << std::flush;
}



