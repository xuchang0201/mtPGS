/*
Multi-trait assisted Polygenic Scores (mtPGS)
*/

#ifndef __MVLMM_FUNCTION_H__
#define __MVLMM_FUNCTION_H__

#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <fstream>

#include <armadillo>

using namespace std;
using namespace arma;

// ALLELE class
class ALLELE {
public:
	long pos;  // bim file position
	string a1;
	string a2;
	double maf;
};

// ALLELEB class
class ALLELEB {
public:
	long pos; // bim file position
	long ps;  // snp position
	string a1;
	string a2;
	double maf;
};

// BLOCK class
class BLOCK {
public:
	string chr; 
	long start;
	long end;
};

// SUMM class
class SUMM {
public:
	int chr;
	string snp;
	long ps; // snp position
	string a1;
	string a2;
	double maf;
	double z;
	double P;
};

// SUMME class
class SUMMS {
public:
	string snp; 
	string a1; 
	double maf; 
	double z; 
};

// SUMMC class
class SUMMC {
public:
	string snp;
	string a1;
	double maf; 
	double z1; 
	double z2; 
};

// SUMMP class
class SUMMP {
public:
	string snp;
	double z1; 
	double z2; 
	long pos; // bim file position
	long ps;  // snp position
};

// pos class
class POS {
public:
	int chr;
	string snp; 
	long ps; // snp position
	long pos; // bim file position
	string a1;
	string a2;
	double maf;
	double z;
	double P;
};

// block and pos class
class INFO {
public:
	string snp;
	long ps; // snp position
	int pos; // bim file position
	int block;
	string a1;
	double maf;
	double z;
	double P;
};

// effect class (output)
class EFF{
public: 
	string snp;
	string a1;
	double maf;
	double beta;
};

// input and output function class
class IO {
public:
	int getRow(string infile);                                                           // get row number
	int readBlock(string infile, char *separator, vector <BLOCK> &block);                // input block information
	int readBim(int n_ref, string ref_str, char *separator, 
				map<string, ALLELE> &bim, bool constr);                                  // input bim file for mtPGS
	double calP(double beta, double se, int sampleSize);                                 // calculate P value from GEMMA output
	int readSumm(string summ_str, char *separator, vector<SUMM> &summ);                  // input gemma summary data
	void readSNPIm(const int pos, int ni_test, const vector<int> &indicator_idv, 
	               ifstream &infile, vec &geno, double &maf);                            // input genotype data
	double getWalltime();                                                                // get wall time 
};

// SNP process function class
class SNPPROC {
public:
	int addBlock(vector<POS> summ, vector <BLOCK> block, vector<INFO> &pos_block);                        // add bolck information
	void nomalizeVec(vec &x);                                                                             // normalize vector
	int matchRef(vector<SUMM> summ, map<string, ALLELE > bim, vector<POS> &inter_snp, double mafMax);     // match summary statistics to reference panel
	int snps_intersect(vector < vector <string> > snps_multi_traits, vector <string> &overlap_snps);      // extract overlapping SNPs across traits
};

// clear vector
template < class T >
void clearVector( vector< T >& vt ) 
{
    vector< T > vtTemp; 
    vtTemp.swap( vt );
}

bool sortP(const INFO &snpinfo1, const INFO &snpinfo2);         // sort SNP by P
bool sortPOS(const INFO &snpinfo1, const INFO &snpinfo2);       // sort SNP by POS
bool sortPS(const SUMMP &summp1, const SUMMP &summp2);          // sort SNP by POS
void printProgBar(int percent);

#endif 
